// deadlock_traced.cpp
//
// Deliberately dangerous version: BOTH wait() sections take {m1,m2}
// together via std::lock(lk1,lk2), same as the working handoff
// sections do. Every lock attempt/acquisition is logged with a
// timestamp so the exact point of the deadlock is visible in the
// output, not just asserted.

#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

std::mutex m1, m2;
std::condition_variable cvData;
std::condition_variable cvResult;
std::mutex coutMutex;

std::vector<int> toConsumer;
std::vector<int> toProducer;
bool dataReady = false;
bool resultReady = false;

using clk = std::chrono::steady_clock;
clk::time_point t_start;

void log(const std::string& who, const std::string& msg) 
{
    std::lock_guard<std::mutex> lk(coutMutex);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(clk::now() - t_start).count();
    std::cout << "[" << ms << "ms] " << who << ": " << msg << std::endl;  // endl flushes -- must not get lost when the process is later killed
}

constexpr int BATCH = 5;

void producer() 
{
    int next = 0;
    while (true) 
	{
        std::vector<int> batch(BATCH);
        for (int& v : batch) 
			v = next++;

        log("producer", "locking {m1,m2} for handoff...");
        {
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);
            log("producer", "acquired {m1,m2}, writing data, releasing both");
            toConsumer = std::move(batch);
            dataReady = true;
        }
        cvData.notify_one();

        std::vector<int> processed;
        log("producer", "locking {m1,m2} before waiting on cvResult...");
        {
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);
            log("producer", "acquired {m1,m2}. Calling cvResult.wait(lk2,...) -- releases ONLY m2. m1 stays LOCKED by me until resultReady.");
            cvResult.wait(lk2, [] { return resultReady; });
            log("producer", "cvResult.wait() returned");
            processed = std::move(toProducer);
            resultReady = false;
        }

        log("producer", "printing result");
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
        log("consumer", "locking {m1,m2} before waiting on cvData...");
        {
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);   // if producer is holding m1 hostage, this is where we get stuck
            log("consumer", "acquired {m1,m2}. Calling cvData.wait(lk1,...)");
            cvData.wait(lk1, [] { return dataReady; });
            log("consumer", "cvData.wait() returned");
            batch = std::move(toConsumer);
            dataReady = false;
        }

        std::vector<int> processed;
        processed.reserve(batch.size());
        for (int v : batch) 
			processed.push_back(v * v);

        log("consumer", "locking {m1,m2} to send result back...");
        {
            std::unique_lock<std::mutex> lk1(m1, std::defer_lock);
            std::unique_lock<std::mutex> lk2(m2, std::defer_lock);
            std::lock(lk1, lk2);
            log("consumer", "acquired {m1,m2}, writing result, releasing both");
            toProducer = std::move(processed);
            resultReady = true;
        }
        cvResult.notify_one();
    }
}

int main() 
{
    t_start = clk::now();
    std::thread t1(consumer);
    std::thread t2(producer);
    t1.join();
    t2.join();
}
