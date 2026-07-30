// cv_exchange_deadlock.cpp
//
// Option (a) from the design, reduced to its smallest failing case:
// TWO neighbouring cubes, each of which has particles to send to the
// other -- which is the normal situation, not a rare one.
//
// Each thread does exactly what the design says:
//     1. publish my particle list for the neighbour
//     2. notify the neighbour
//     3. wait for the neighbour's confirmation
//     4. delete whatever was accepted
//
// Step 3 is the problem. While a thread sits in wait(), it is not
// running, so it cannot perform step 1..2's counterpart -- judging the
// list the neighbour just sent IT. Both threads reach step 3 and stop.
//
// Build: g++ -std=c++20 -O2 -pthread cv_exchange_deadlock.cpp -o cvdead
// Run:   timeout 3 ./cvdead    ; echo $?     -> 124 (killed, hung)

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

    std::vector<int> particles;     // sender -> receiver
    std::vector<char> verdict;      // receiver -> sender
    bool requestReady = false;
    bool replyReady   = false;
};

Link aToB;   // A sends particles to B on this link
Link bToA;   // B sends particles to A on this link

std::mutex coutMutex;
using clk = std::chrono::steady_clock;
clk::time_point t0;

void log(const char* who, const char* msg) {
    std::lock_guard<std::mutex> lk(coutMutex);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(clk::now() - t0).count();
    std::cout << "[" << ms << "ms] " << who << ": " << msg << std::endl;
}

// send my list on 'out', then block until the neighbour answers
void sendAndAwaitConfirmation(Link& out, const char* who) {
    {
        std::lock_guard<std::mutex> lk(out.m);
        out.particles = {1, 2, 3};
        out.requestReady = true;
    }
    out.cvRequest.notify_one();
    log(who, "sent my particles, now waiting for confirmation");

    std::unique_lock<std::mutex> lk(out.m);
    out.cvReply.wait(lk, [&] { return out.replyReady; });   // <-- both threads stop here
    log(who, "confirmation received");
}

// judge whatever the neighbour sent me on 'in', and answer
void judgeAndReply(Link& in, const char* who) {
    std::unique_lock<std::mutex> lk(in.m);
    in.cvRequest.wait(lk, [&] { return in.requestReady; });
    in.verdict.assign(in.particles.size(), 1);
    in.replyReady = true;
    lk.unlock();
    in.cvReply.notify_one();
    log(who, "judged the neighbour's particles and replied");
}

void threadA() {
    sendAndAwaitConfirmation(aToB, "A");   // blocks until B judges -- but B never gets there
    judgeAndReply(bToA, "A");
}

void threadB() {
    sendAndAwaitConfirmation(bToA, "B");   // blocks until A judges -- but A never gets there
    judgeAndReply(aToB, "B");
}

int main() {
    t0 = clk::now();
    std::thread t1(threadA), t2(threadB);
    t1.join();
    t2.join();
    std::cout << "finished (this line is never reached)\n";
}
