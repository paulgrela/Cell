// barrier_bench.cpp — pure synchronization overhead, no payload work
#include <barrier>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <mutex>
#include <thread>

constexpr long ITERS = 2'000'000;

class CvBarrier 
{
public:
    explicit CvBarrier(unsigned count) : threshold_(count), count_(count), generation_(0) {}
    void arrive_and_wait() 
	{
        std::unique_lock<std::mutex> lk(m_);
        unsigned gen = generation_;
		
        if (--count_ == 0) 
		{
            generation_++;
            count_ = threshold_;
            lk.unlock();
            cv_.notify_all();
        } 
		else 
		{
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

template <typename Barrier>
void bench(const char* name) 
{
    Barrier b(2);
    auto worker = [&b] 
	{
        for (long i = 0; i < ITERS; ++i)
            b.arrive_and_wait();
    };
	
    auto t0 = std::chrono::steady_clock::now();
    std::thread t1(worker), t2(worker);
    t1.join();
    t2.join();
    auto t1_ = std::chrono::steady_clock::now();
    double ns = std::chrono::duration<double, std::nano>(t1_ - t0).count();
    std::printf("%-14s: %8.1f ms total, %7.1f ns / crossing (per thread)\n", name, ns / 1e6, ns / ITERS);
}

int main() 
{
    std::printf("ITERS = %ld per thread, two threads\n\n", ITERS);
	
    for (int run = 0; run < 3; ++run) 
	{
        bench<std::barrier<>>("std::barrier");
        bench<CvBarrier>("cv barrier");
        std::printf("---\n");
    }
}
