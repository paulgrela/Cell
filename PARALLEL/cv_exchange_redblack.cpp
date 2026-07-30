// cv_exchange_redblack.cpp
//
// The fix for cv_exchange_deadlock.cpp. Same two peers, same two
// links, same one-mutex-per-direction design, same condition
// variables -- the ONLY change is that the two sides no longer run
// identical code.
//
// The deadlock came from role symmetry: both peers demanded a reply
// before serving one. So break the symmetry. Colour every cube by
// (gx+gy+gz)%2 and give the two colours opposite orders at each shared
// wall:
//
//   both:  publish my particle list          (non-blocking)
//   red:   wait for my verdict,  THEN judge the neighbour's list
//   black: judge the neighbour's list, THEN wait for my verdict
//
// Now black is always available to serve when red asks, so the cycle
// has an entry point. Red never waits on red, black never waits on
// black, so no cycle can form across the wider neighbour graph either.
//
// This is correct -- but it costs one blocking round-trip per wall,
// six walls per cube, versus three barrier crossings for the whole
// grid. It is the answer to "can option (a) be made to work", not a
// recommendation.
//
// Build: g++ -std=c++20 -O2 -pthread cv_exchange_redblack.cpp -o rb

#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

struct Link {                       // one direction of the A <-> B channel
    std::mutex m;
    std::condition_variable cvRequest;
    std::condition_variable cvReply;
	
    std::vector<int> particles;
    std::vector<char> verdict;
    bool requestReady = false;
    bool replyReady   = false;
};

Link aToB, bToA;

std::mutex coutMutex;
void log(const char* who, const char* msg) {
    std::lock_guard<std::mutex> lk(coutMutex);
    std::cout << who << ": " << msg << std::endl;
}

void publish(Link& out, std::vector<int> mine) {
    {
        std::lock_guard<std::mutex> lk(out.m);
        out.particles = std::move(mine);
        out.requestReady = true;
    }
    out.cvRequest.notify_one();
}

void judgeAndReply(Link& in, const char* who) {
    std::unique_lock<std::mutex> lk(in.m);
    in.cvRequest.wait(lk, [&] { return in.requestReady; });
    in.verdict.assign(in.particles.size(), 1);   // stand-in for the occupancy check
    in.requestReady = false;
    in.replyReady = true;
    lk.unlock();
    in.cvReply.notify_one();
    log(who, "judged the neighbour's particles and replied");
}

std::vector<char> awaitVerdict(Link& out, const char* who) {
    std::unique_lock<std::mutex> lk(out.m);
    out.cvReply.wait(lk, [&] { return out.replyReady; });
    out.replyReady = false;
    log(who, "got my verdict back");
    return std::move(out.verdict);
}

void redPeer() {                    // (gx+gy+gz) even
    for (int round = 0; round < 3; ++round) {
        publish(aToB, {1, 2, 3});
        awaitVerdict(aToB, "A(red)");        // wait first...
        judgeAndReply(bToA, "A(red)");       // ...then serve
    }
}

void blackPeer() {                  // (gx+gy+gz) odd
    for (int round = 0; round < 3; ++round) {
        publish(bToA, {7, 8, 9});
        judgeAndReply(aToB, "B(black)");     // serve first...
        awaitVerdict(bToA, "B(black)");      // ...then wait
    }
}

int main() {
    std::thread t1(redPeer), t2(blackPeer);
    t1.join();
    t2.join();
    std::cout << "finished all rounds -- no deadlock\n";
}
