// producer_consumer_two_mutex.cpp
//
// Same bidirectional producer/consumer exchange as before, but split
// into two independent mutexes instead of one shared std::mutex:
//   m1 guards the "data" channel (toConsumer, dataReady) -- paired with cvData
//   m2 guards the "result" channel (toProducer, resultReady) -- paired with cvResult
//
// The two brief handoff sections (no wait() involved) safely take
// both locks together with std::lock(lk1, lk2) -- the manual,
// scoped_lock-style deadlock-avoiding acquisition. The two wait()
// sections each take only their own single relevant lock.

#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

std::mutex m1;                    // guards toConsumer + dataReady -- paired with cvData
std::mutex m2;                    // guards toProducer + resultReady -- paired with cvResult
std::condition_variable cvData;
std::condition_variable cvResult;

std::vector<int> toConsumer;
std::vector<int> toProducer;
bool dataReady = false;
bool resultReady = false;

constexpr int BATCH = 5;

void producer() 
{
    int next = 0;
    while (true) 
	{
		std::vector<int> batch(BATCH);
        for (int& v : batch) 
			v = next++;

        {   
			// hand data to the consumer -- no wait() here, safe to take both locks together
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);
            toConsumer = std::move(batch);
            dataReady = true;
        }
        cvData.notify_one();

        std::vector<int> processed;
        {   
			// suspended here until the consumer sends processed data back -- only m2 is involved
            std::unique_lock<std::mutex> lk2(m2);
            cvResult.wait(lk2, [] { return resultReady; });
            processed = std::move(toProducer);
            resultReady = false;              // ready -> false right after reading
        }

        std::cout << "producer received:";
        for (int v : processed) 
			std::cout << ' ' << v;
        std::cout << '\n';
    }
}

void consumer() 
{
    while (true) 
	{
        std::vector<int> batch;
        {   
			// suspended here until data arrives -- only m1 is involved
            std::unique_lock<std::mutex> lk1(m1);
            cvData.wait(lk1, [] { return dataReady; });
            batch = std::move(toConsumer);
            dataReady = false;                // ready -> false right after reading
        }

        std::vector<int> processed;
        processed.reserve(batch.size());
        for (int v : batch)
			processed.push_back(v * v);        // stand-in for the real processing

        {   
			// send processed data back -- no wait() here, safe to take both locks together
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);
            toProducer = std::move(processed);
            resultReady = true;
        }
        cvResult.notify_one();
    }
}

int main() 
{
    std::thread t1(consumer);
    std::thread t2(producer);
    t1.join();   // both loops run forever, by design
    t2.join();
}
