// barrier_demo_cv.cpp
//
// Same two-peer, two-barrier exchange as barrier_demo.cpp, but the
// barrier itself is hand-built from std::condition_variable +
// std::mutex instead of std::barrier -- to show it's possible without
// the C++20 primitive, and roughly what std::barrier does internally.

#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

// A reusable ("cyclic") barrier for a fixed number of threads.
// arrive_and_wait() blocks until 'count' threads have called it, then
// releases all of them together and resets for the next round.
//
// The generation counter is what makes it safely reusable: each
// waiter captures its own round number (gen) and only wakes once
// generation_ has moved past it -- so a thread that races ahead into
// the next round can't corrupt the count for a round still finishing.

class Barrier 
{
public:
    explicit Barrier(unsigned count) : threshold_(count), count_(count), generation_(0) {}

    void arrive_and_wait() 
	{
        std::unique_lock<std::mutex> lk(m_);
        unsigned gen = generation_;
        if (--count_ == 0) 
		{
            // last thread to arrive this round: open the gate, reset for next round
            generation_++;
            count_ = threshold_;
            lk.unlock();
            cv_.notify_all();
        } 
		else 
		{
            // not last: sleep until the round we captured has ended
            cv_.wait(lk, [this, gen] { return gen != generation_; });
        }
    }

private:
    std::mutex m_;
    std::condition_variable cv_;
    const unsigned threshold_;
    unsigned count_;
    unsigned generation_;
};

constexpr int BATCH = 5;

std::vector<int> dataA, dataB;
std::vector<int> resultA, resultB;

Barrier dataBarrier(2);
Barrier resultBarrier(2);

std::mutex coutMutex;

void printResult(const char* who, const std::vector<int>& result) 
{
    std::lock_guard<std::mutex> lk(coutMutex);
    std::cout << who << " received:";
    for (int v : result) 
		std::cout << ' ' << v;
    std::cout << '\n';
}

void threadA() 
{
    int next = 0;
    while (true) 
	{
        dataA.resize(BATCH);
        for (int& v : dataA) 
			v = next++;

        dataBarrier.arrive_and_wait();

        resultA.resize(dataB.size());
        for (size_t i = 0; i < dataB.size(); ++i)
            resultA[i] = dataB[i] * dataB[i];

        resultBarrier.arrive_and_wait();

        printResult("A", resultB);
    }
}

void threadB() 
{
    int next = 100;
    while (true) 
	{
        dataB.resize(BATCH);
        for (int& v : dataB) 
			v = next++;

        dataBarrier.arrive_and_wait();

        resultB.resize(dataA.size());
        for (size_t i = 0; i < dataA.size(); ++i)
            resultB[i] = dataA[i] + 1;

        resultBarrier.arrive_and_wait();

        printResult("B", resultA);
    }
}

int main() {
    std::thread t1(threadA);
    std::thread t2(threadB);
    t1.join();
    t2.join();
}
