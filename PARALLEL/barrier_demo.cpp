// barrier_demo.cpp
//
// Two threads as true peers, not producer/consumer: each produces its
// own data concurrently, then BOTH rendezvous at a barrier before
// either is allowed to read what the other produced. Compute a result
// from it, rendezvous at a second barrier, then read each other's
// result. Repeat forever.
//
// Uses std::barrier (C++20) -- the standard-library primitive for
// exactly "block until everyone has arrived".

#include <barrier>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

constexpr int BATCH = 5;

// Each thread owns and writes only its own vectors. The barriers --
// not a mutex -- are what make it safe for the *other* thread to read
// them: crossing a barrier happens-after everyone's writes from the
// phase before it.
std::vector<int> dataA, dataB;
std::vector<int> resultA, resultB;

std::barrier dataBarrier(2);     // phase 1: both sides done producing data
std::barrier resultBarrier(2);   // phase 2: both sides done producing results

std::mutex coutMutex;            // unrelated to the barriers -- just keeps prints from interleaving

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

        dataBarrier.arrive_and_wait();          // wait for B to finish writing dataB too

        resultA.resize(dataB.size());
        for (size_t i = 0; i < dataB.size(); ++i)
            resultA[i] = dataB[i] * dataB[i];

        resultBarrier.arrive_and_wait();        // wait for B to finish writing resultB too

        printResult("A", resultB);
    }
}

void threadB() 
{
    int next = 100;
    while (true) {
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
