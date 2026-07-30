// cell_diffusion_chessboard.cpp
//
// Lattice particle diffusion over an M x N x L grid of thread-cubes,
// one std::thread per cube, exchanging particles with the 6 full-wall
// face neighbours.
//
// Synchronisation: std::condition_variable + wait/notify ONLY.
// No barriers. Deadlock freedom comes from the 3D chessboard.
//
// ---------------------------------------------------------------- chessboard
//
// colour(cube) = (gx + gy + gz) % 2
//
// Face neighbours differ in exactly one coordinate, so they ALWAYS have
// opposite colours. Every wall therefore joins one white and one black
// cube, and the two sides can be given opposite code:
//
//   half-step A:  white = SENDER      black = JUDGE
//   half-step B:  black = SENDER      white = JUDGE
//
//   sender: publish my 6 lists, wait for 6 verdicts
//   judge : wait for a list, check acceptance, reply
//
// Deadlock freedom, precisely: a sender's publish step blocks on
// nothing, and it happens before the sender's own wait. So when a judge
// waits for a list, the list is already there (or its author is on the
// way to publishing it, blocked on nothing). White never waits on
// white, black never waits on black -- the wait-for graph is bipartite
// and every edge points from a blocked thread to a runnable one.
//
// ---------------------------------------------------------------- two rules
//
// RULE 1: a sender must publish ALL SIX lists before waiting on ANY
// verdict. Publishing blocks on nothing, so doing all six up front is
// free -- and it makes the half-step deadlock-free no matter what order
// the six walls are then serviced in.
//
// The chessboard alone does NOT give you this. Interleaving instead
// (publish wall d, wait on wall d, publish wall d+1, ...) still happens
// to work IF every thread walks its directions in the same global order
// 0..5, because then each waits-for chain is monotone in one grid
// coordinate and terminates at the boundary. That is an accident of the
// loop order, not a property of the algorithm: measured here, giving
// each cube a different rotation of the direction order deadlocks the
// interleaved version within seconds, while this publish-all-first
// version survives the same rotation unharmed. Do not rely on the
// accident -- it will break the first time someone services walls in
// arrival order or sorts them by particle count.
//
// RULE 2: acceptance is judged against a LIVE occupancy array, not a
// snapshot. occ[] still carries every particle waiting to depart (the
// commit has not run), so an arrival can never take the cell of a
// particle that might yet be rejected. And because occ[] is set the
// moment a particle is accepted, two different neighbours cannot both
// be given the same cell -- an edge cell such as (0,0,z) is reachable
// from both the -x and the -y neighbour, and without this it gets
// double-booked.
//
// ---------------------------------------------------------------- round
//
//   diffuse locally, build 6 outgoing lists   (all threads, parallel)
//   white: send-then-judge   |   black: judge-then-send
//   commit: delete the departures that were accepted
//
// There is no global synchronisation point, so cubes may drift ahead of
// each other by a bounded amount. The wall channels supply all the
// ordering that is actually needed.
//
// Build: g++ -std=c++20 -O2 -pthread cell_diffusion_chessboard.cpp -o sim
// Run:   ./sim [M] [N] [L] [cubeEdge] [rounds]
//        ./sim 4 4 4   -> 64 threads     ./sim 4 4 8 -> 128     ./sim 4 8 8 -> 256

#include <array>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

// ---------------------------------------------------------------- config

int GM = 4, GN = 4, GL = 4;   // thread-cube grid
int CS = 16;                  // lattice cells along one cube edge
int ROUNDS = 20;
double FILL = 0.20;

// direction encoding: 0=-x 1=+x 2=-y 3=+y 4=-z 5=+z
inline int opposite(int d) 
{ 
	return d ^ 1; 
}

// ---------------------------------------------------------------- types

struct Particle 
{
    uint32_t id;
    int16_t x, y, z;          // lattice coords, local to the owning cube
};

// One shared wall between a white cube and a black cube.
// Slot 0 carries white -> black traffic, slot 1 carries black -> white.
// Separate slots matter: a single slot reused across the two half-steps
// lets a stale "ready" flag from half-step A satisfy a wait in half-step B.
struct Wall 
{
    std::mutex m;
    std::condition_variable cv;               // one CV, notify_all: only 2 threads here
    std::vector<Particle> proposals[2];
    std::vector<uint8_t>  verdicts[2];
    bool proposalsReady[2] = { false, false };
    bool verdictsReady[2]  = { false, false };
};

struct Cube 
{
    int gx = 0, gy = 0, gz = 0;
    int colour = 0;                            // 0 = white, 1 = black
    std::array<int, 6> wall{};                 // wall index, or -1 at the outer boundary

    std::vector<Particle> particles;
    std::vector<uint8_t>  occ;                 // CS^3 occupancy

    std::array<std::vector<Particle>, 6> outParticles;  // sent to the neighbour
    std::array<std::vector<uint32_t>, 6> outSrc;        // matching index into particles (private)
    std::array<std::vector<uint8_t>, 6>  verdictIn;     // what the neighbour decided
};

std::vector<Cube> cubes;
std::vector<std::unique_ptr<Wall>> walls;      // Wall holds a mutex, so it cannot be moved

inline size_t cellIndex(int x, int y, int z) { return (size_t)((x * CS + y) * CS + z); }
inline int cubeIndex(int gx, int gy, int gz)  { return (gx * GN + gy) * GL + gz; }

// ---------------------------------------------------------------- topology

int OFF_Y = 0, OFF_Z = 0;

int wallIndex(int gx, int gy, int gz, int d) 
{
    switch (d) 
	{
        case 0: return gx > 0 ? ((gx - 1) * GN + gy) * GL + gz : -1;
        case 1: return gx < GM - 1 ? (gx * GN + gy) * GL + gz : -1;
        case 2: return gy > 0 ? OFF_Y + (gx * (GN - 1) + (gy - 1)) * GL + gz : -1;
        case 3: return gy < GN - 1 ? OFF_Y + (gx * (GN - 1) + gy) * GL + gz : -1;
        case 4: return gz > 0 ? OFF_Z + (gx * GN + gy) * (GL - 1) + (gz - 1) : -1;
        case 5: return gz < GL - 1 ? OFF_Z + (gx * GN + gy) * (GL - 1) + gz : -1;
    }
    return -1;
}

// ---------------------------------------------------- local diffusion

// A move crossing more than one face targets an edge/corner neighbour and is
// dropped. A particle proposed for export does NOT free its cell here: we do
// not yet know whether the neighbour will take it. It is freed in commit().
void diffuse(Cube& c, std::mt19937& rng) 
{
    for (int d = 0; d < 6; ++d) 
	{ 
		c.outParticles[d].clear(); 
		c.outSrc[d].clear(); 
	}

    std::uniform_int_distribution<int> step(-1, 1);

    for (uint32_t i = 0; i < c.particles.size(); ++i) 
	{
        Particle& p = c.particles[i];

        int dx = step(rng), dy = step(rng), dz = step(rng);
        if (dx == 0 && dy == 0 && dz == 0) 
			continue;

        int nx = p.x + dx, ny = p.y + dy, nz = p.z + dz;

        int crossings = 0, dir = -1;
        if (nx < 0)        
		{ 
			++crossings; 
			dir = 0; 
		}
        else 
		if (nx >= CS) 
		{ 
			++crossings; 
			dir = 1; 
		}
        
		if (ny < 0)        
		{ 
			++crossings; 
			dir = 2; 
		}
        else 
		if (ny >= CS) 
		{ 
			++crossings; 
			dir = 3; 
		}
		
        if (nz < 0)        
		{ 
			++crossings; 
			dir = 4; 
		}
        else 
		if (nz >= CS) 
		{ 
			++crossings; 
			dir = 5; 
		}

        if (crossings > 1) 
			continue;                     // corner neighbour -> not moved

        if (crossings == 0) 
		{                            	// stays inside this cube
            size_t to = cellIndex(nx, ny, nz);
            if (!c.occ[to]) 
			{                            		// live occ also stops two of my
                c.occ[cellIndex(p.x, p.y, p.z)] = 0;   // own particles taking one cell
                c.occ[to] = 1;
                p.x = (int16_t)nx; p.y = (int16_t)ny; p.z = (int16_t)nz;
            }
			
            continue;
        }

        if (c.wall[dir] < 0) 
			continue;                   // outer boundary of the grid

        // one axis crossed: (v + CS) % CS maps -1 -> CS-1 and CS -> 0
        c.outParticles[dir].push_back(Particle{p.id, (int16_t)((nx + CS) % CS), (int16_t)((ny + CS) % CS), (int16_t)((nz + CS) % CS)});
        c.outSrc[dir].push_back(i);
    }
}

// ---------------------------------------------------- half-step: I am the sender

void sendHalfStep(Cube& c) 
{
    const int slot = c.colour;                 // white sends on slot 0, black on slot 1

    // RULE 1: publish all six first. Publishing blocks on nothing.
    for (int d = 0; d < 6; ++d) 
	{
        if (c.wall[d] < 0)
			continue;
		
        Wall& W = *walls[c.wall[d]];
        {
            std::lock_guard<std::mutex> lk(W.m);
            W.proposals[slot] = c.outParticles[d];
            W.proposalsReady[slot] = true;
        }
		
        W.cv.notify_all();
    }

    // only now block, collecting the six verdicts
    for (int d = 0; d < 6; ++d) 
	{
        if (c.wall[d] < 0) 
		{ 
			c.verdictIn[d].clear(); 
			continue; 
		}
		
        Wall& W = *walls[c.wall[d]];
        std::unique_lock<std::mutex> lk(W.m);
        W.cv.wait(lk, [&]{ return W.verdictsReady[slot]; });
        c.verdictIn[d] = std::move(W.verdicts[slot]);
        W.verdictsReady[slot] = false;         // consume, so the flag is clean next round
    }
}

// ---------------------------------------------------- half-step: I am the judge

void judgeHalfStep(Cube& c) 
{
    const int slot = 1 - c.colour;             // the slot my neighbours send on

    for (int d = 0; d < 6; ++d) 
	{
        if (c.wall[d] < 0) 
			continue;
		
        Wall& W = *walls[c.wall[d]];

        std::unique_lock<std::mutex> lk(W.m);
        W.cv.wait(lk, [&] { return W.proposalsReady[slot]; });
        std::vector<Particle> incoming = std::move(W.proposals[slot]);
        W.proposalsReady[slot] = false;        // consume
        lk.unlock();

        // RULE 2: judged outside the lock, against my own LIVE occ array.
        std::vector<uint8_t> verdict(incoming.size(), 0);
        for (size_t k = 0; k < incoming.size(); ++k) 
		{
            size_t cell = cellIndex(incoming[k].x, incoming[k].y, incoming[k].z);
            if (!c.occ[cell]) 
			{
                c.occ[cell] = 1;               // claim it before another neighbour can
                c.particles.push_back(incoming[k]);
                verdict[k] = 1;
            }
        }

        lk.lock();
        W.verdicts[slot] = std::move(verdict);
        W.verdictsReady[slot] = true;
        lk.unlock();
        W.cv.notify_all();
    }
}

// ---------------------------------------------------- commit

void commit(Cube& c) 
{
    // the judge half-step only appended to c.particles, so the indices
    // recorded in outSrc during diffuse() still address the same particles
    std::vector<uint8_t> departed(c.particles.size(), 0);

    for (int d = 0; d < 6; ++d) 
	{
        if (c.wall[d] < 0) 
			continue;
		
        for (size_t k = 0; k < c.outSrc[d].size(); ++k) 
		{
            if (!c.verdictIn[d][k]) 
				continue;            // rejected -> particle simply stays
			
            uint32_t si = c.outSrc[d][k];
            const Particle& p = c.particles[si];
            c.occ[cellIndex(p.x, p.y, p.z)] = 0;         // now, and only now, free the cell
            departed[si] = 1;
        }
    }

    size_t w = 0;
    for (size_t r = 0; r < c.particles.size(); ++r)
        if (!departed[r]) 
			c.particles[w++] = c.particles[r];
		
    c.particles.resize(w);
}

// ---------------------------------------------------------------- validation

bool validate(size_t expectedTotal) 
{
    size_t total = 0;
    bool ok = true;

    for (size_t ci = 0; ci < cubes.size(); ++ci) 
	{
        Cube& c = cubes[ci];
        total += c.particles.size();

        std::vector<uint8_t> seen(c.occ.size(), 0);
        for (const Particle& p : c.particles) 
		{
            if (p.x < 0 || p.x >= CS || p.y < 0 || p.y >= CS || p.z < 0 || p.z >= CS) 
			{
                std::cout << "  FAIL: cube " << ci << " particle out of bounds\n";
                ok = false; break;
            }
            size_t cell = cellIndex(p.x, p.y, p.z);
            if (seen[cell]) 
			{
                std::cout << "  FAIL: cube " << ci << " has two particles in one cell\n";
                ok = false; break;
            }
            seen[cell] = 1;
            if (!c.occ[cell]) 
			{
                std::cout << "  FAIL: cube " << ci << " particle on a cell occ says is free\n";
                ok = false; break;
            }
        }

        size_t occCount = 0;
        for (uint8_t v : c.occ) occCount += v;
        if (occCount != c.particles.size()) 
		{
            std::cout << "  FAIL: cube " << ci << " occ count " << occCount << " != particle count " << c.particles.size() << "\n";
            ok = false;
        }
    }

    if (total != expectedTotal) 
	{
        std::cout << "  FAIL: total particles " << total << " != expected " << expectedTotal << "\n";
        ok = false;
    }
	
    return ok;
}

// ---------------------------------------------------------------- main

int main(int argc, char** argv) 
{
    if (argc > 1) GM = std::atoi(argv[1]);
    if (argc > 2) GN = std::atoi(argv[2]);
    if (argc > 3) GL = std::atoi(argv[3]);
    if (argc > 4) CS = std::atoi(argv[4]);
    if (argc > 5) ROUNDS = std::atoi(argv[5]);

    const int nCubes = GM * GN * GL;

    const int nWallX = (GM - 1) * GN * GL;
    const int nWallY = GM * (GN - 1) * GL;
    const int nWallZ = GM * GN * (GL - 1);
    OFF_Y = nWallX;
    OFF_Z = nWallX + nWallY;

    walls.reserve(nWallX + nWallY + nWallZ);
    for (int i = 0; i < nWallX + nWallY + nWallZ; ++i)
        walls.push_back(std::make_unique<Wall>());

    cubes.resize(nCubes);
    for (int gx = 0; gx < GM; ++gx)
		for (int gy = 0; gy < GN; ++gy)
			for (int gz = 0; gz < GL; ++gz) 
			{
				Cube& c = cubes[cubeIndex(gx, gy, gz)];
				c.gx = gx; c.gy = gy; c.gz = gz;
				c.colour = (gx + gy + gz) % 2;
				c.occ.assign((size_t)CS * CS * CS, 0);
				
				for (int d = 0; d < 6; ++d) 
					c.wall[d] = wallIndex(gx, gy, gz, d);
			}

    // seed particles
    uint32_t nextId = 0;
    size_t expectedTotal = 0;
	
    for (int ci = 0; ci < nCubes; ++ci) 
	{
        Cube& c = cubes[ci];
        std::mt19937 rng(1234u + ci);
        std::uniform_int_distribution<int> pos(0, CS - 1);
        size_t want = (size_t)(FILL * CS * CS * CS);
		
        while (c.particles.size() < want) 
		{
            int x = pos(rng), y = pos(rng), z = pos(rng);
            size_t cell = cellIndex(x, y, z);
			
            if (c.occ[cell]) 
				continue;
			
            c.occ[cell] = 1;
            c.particles.push_back(Particle{nextId++, (int16_t)x, (int16_t)y, (int16_t)z});
        }
		
        expectedTotal += c.particles.size();
    }

    int nWhite = 0;
    for (const Cube& c : cubes) 
		if (c.colour == 0) 
			++nWhite;

    std::cout << "grid " << GM << "x" << GN << "x" << GL << " = " << nCubes << " threads"<< " (" << nWhite << " white, " << nCubes - nWhite << " black), " << walls.size() << " walls, cube edge " << CS << ", particles " << expectedTotal << ", rounds " << ROUNDS << "\n";

    auto worker = [&](int id) 
	{
        Cube& c = cubes[id];
        std::mt19937 rng(9876u + id);
        for (int r = 0; r < ROUNDS; ++r) 
		{
            diffuse(c, rng);
            if (c.colour == 0) 
			{          // white: send first, then judge
                sendHalfStep(c);
                judgeHalfStep(c);
            } 
			else 
			{                      // black: judge first, then send
                judgeHalfStep(c);
                sendHalfStep(c);
            }
            commit(c);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(nCubes);
	
    for (int i = 0; i < nCubes; ++i) 
		threads.emplace_back(worker, i);
	
    for (auto& t : threads) 
		t.join();

    std::cout << (validate(expectedTotal) ? "OK: particle count conserved, no cell double-occupied\n" : "VALIDATION FAILED\n");
    return 0;
}
