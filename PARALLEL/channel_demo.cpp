// channel_demo.cpp
//
// Reusable single-slot, thread-safe "rendezvous" channel, built on
// std::mutex + std::condition_variable, with symmetric send()/receive():
//   - send()    blocks until any previous value has been consumed
//   - receive() blocks until a value is available, and clears the slot
//               immediately after reading ("ready" -> false right after reading)
//
// Demo below: a producer/consumer pair that hands data back and forth
// forever -- producer sends a raw batch, consumer processes it and sends
// the result back, repeat.

#include <condition_variable>
#include <iostream>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

template <typename T>
class Channel 
{

public:

    void send(T value) 
	{
        std::unique_lock<std::mutex> lk(m_);
        cvSlotFree_.wait(lk, [this] { return !value_.has_value(); });
        value_ = std::move(value);
        lk.unlock();
        cvHasData_.notify_one();
    }

    T receive() 
	{
        std::unique_lock<std::mutex> lk(m_);
        cvHasData_.wait(lk, [this] { return value_.has_value(); });
        T out = std::move(*value_);
        value_.reset();               									// ready -> false right after reading
        lk.unlock();
        cvSlotFree_.notify_one();     								// let a blocked sender proceed
        return out;
    }

private:
    std::mutex m_;
    std::condition_variable cvHasData_;
    std::condition_variable cvSlotFree_;
    std::optional<T> value_;
};

// ---- usage: same two-way producer/consumer as before, now symmetric ----

constexpr int BATCH = 5;

Channel<std::vector<int>> toConsumer;
Channel<std::vector<int>> toProducer;

void producer() 
{
    int next = 0;
    while (true) 
	{
        std::vector<int> batch(BATCH);
        for (int& v : batch) 
			v = next++;

        toConsumer.send(std::move(batch));                  		// suspended if consumer hasn't read the last batch yet
        std::vector<int> processed = toProducer.receive();  	// suspended until consumer sends back

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
        std::vector<int> batch = toConsumer.receive();      	// suspended until data arrives

        std::vector<int> processed;
        processed.reserve(batch.size());
        for (int v : batch)
            processed.push_back(v * v);                      		// stand-in for real processing

        toProducer.send(std::move(processed));                	// suspended if producer hasn't read the last result yet
    }
}

int main() 
{
    std::thread t1(consumer);
    std::thread t2(producer);
    t1.join();
    t2.join();
}
