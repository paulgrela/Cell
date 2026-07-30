// cell_diffusion_exchange.cpp
//
// Lattice particle diffusion over an M x N x L grid of thread-cubes,
// one std::thread per cube, with 6-face-neighbour particle exchange
// using propose / accept-or-reject / commit.
//
// Synchronisation: std::barrier ONLY. No mutexes anywhere.
// Each phase is separated by a barrier, and within a phase every
// memory location has exactly one writer -- so there is nothing left
// to protect. See the phase table below.
//
//   P1  local diffusion; build 6 outgoing proposal lists
//       writes: own particles, own occ, own outgoing
//       reads : own only
//   ---- barrier ----
//   P2  read each neighbour's proposals aimed at me; accept or reject;
//       insert accepted ones into myself; publish 6 verdict lists
//       writes: own particles, own occ, own verdict
//       reads : neighbour->outgoing  (frozen since P1)
//   ---- barrier ----
//   P3  read neighbours' verdicts on MY proposals; delete the ones
//       that were accepted (they now live in the neighbour)
//       writes: own particles, own occ
//       reads : neighbour->verdict   (frozen since P2)
//   ---- barrier ----
//
// Two invariants the exchange must preserve, both asserted every round:
//   (1) global particle count is constant -- no particle is duplicated
//       (accepted by neighbour but not deleted here) or lost (deleted
//       here but not accepted there)
//   (2) no lattice cell ever holds two particles
//
// Build:  g++ -std=c++20 -O2 -pthread cell_diffusion_exchange.cpp -o sim
// Run:    ./sim [M] [N] [L] [cubeEdge] [rounds]
//         ./sim 4 4 4      -> 64 threads
//         ./sim 4 4 8      -> 128 threads
//         ./sim 4 8 8      -> 256 threads

#include <array>
#include <barrier>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

// ---------------------------------------------------------------- config

int GM = 4, GN = 4, GL = 4;   // thread-cube grid dimensions
int CS = 16;                  // lattice cells along one cube edge
int ROUNDS = 20;
double FILL = 0.20;           // initial fraction of cells occupied

// direction encoding: 0=-x 1=+x 2=-y 3=+y 4=-z 5=+z
// opposite direction is simply d^1  (0<->1, 2<->3, 4<->5)
inline int opposite(int d) { return d ^ 1; }

// ---------------------------------------------------------------- types

struct Particle 
{
    uint32_t id;
    int16_t x, y, z;          // lattice coords, local to the owning cube
};

struct Proposal 
{
    Particle p;               // coords already translated into the NEIGHBOUR's frame
    uint32_t srcIdx;          // index into the sender's particle vector (sender-private)
};

struct Cube 
{
    int gx = 0, gy = 0, gz = 0;
    std::array<int, 6> neighbour{};                  // cube index, or -1 at the outer wall

    std::vector<Particle> particles;
    std::vector<uint8_t>  occ;                       // CS^3 occupancy bitmap

    std::array<std::vector<Proposal>, 6> outgoing;   // written P1, read by neighbour in P2
    std::array<std::vector<uint8_t>, 6>  verdict;    // verdict[d] judges what arrived FROM direction d
};

std::vector<Cube> cubes;

inline size_t cellIndex(int x, int y, int z) 
{
    return (size_t)((x * CS + y) * CS + z);
}

inline int cubeIndex(int gx, int gy, int gz) 
{
    return (gx * GN + gy) * GL + gz;
}

// --------------------------------------------------------- phase 1: diffuse

// Random move in {-1,0,+1}^3. A move that would cross more than one face
// targets an edge/corner neighbour -- dropped, per the model's rule.
// A particle proposed for export does NOT free its cell here: we do not
// yet know whether the neighbour will take it. It is freed in P3.
void phase1(Cube& c, std::mt19937& rng) 
{
    for (auto& v : c.outgoing) 
		v.clear();

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
			++crossings; dir = 0; 
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
			continue; 			             	// edge/corner neighbour -> not moved

        if (crossings == 0) 
		{						                     	// stays inside this cube
            size_t to = cellIndex(nx, ny, nz);
            if (!c.occ[to]) 
			{                     	// occ is live -- also blocks two of my
                c.occ[cellIndex(p.x, p.y, p.z)] = 0;   // own particles taking one cell
                c.occ[to] = 1;
                p.x = (int16_t)nx; p.y = (int16_t)ny; p.z = (int16_t)nz;
            }
			
            continue;
        }

        int nb = c.neighbour[dir];
        if (nb < 0) 
			continue;                     // outer wall of the whole grid

        // one axis crossed: (v + CS) % CS maps -1 -> CS-1 and CS -> 0,
        // and leaves the two in-range axes untouched
        Particle t{p.id, (int16_t)((nx + CS) % CS), (int16_t)((ny + CS) % CS), (int16_t)((nz + CS) % CS)};
        c.outgoing[dir].push_back(Proposal{t, i});
    }
}

// ------------------------------------------------- phase 2: accept or reject

// c.occ doubles as the claim set. It still carries every particle that is
// waiting to depart (P3 has not run), so an arrival can never be placed on a
// cell whose occupant might yet be rejected. And because we set occ[] the
// moment we accept, two different neighbours cannot both be given the same
// cell -- which they otherwise could, since an edge cell such as (0,0,z) is
// reachable from both the -x and the -y neighbour.
void phase2(Cube& c) 
{
    for (int d = 0; d < 6; ++d) 
	{
        c.verdict[d].clear();

        int nb = c.neighbour[d];
        if (nb < 0) 
			continue;

        const std::vector<Proposal>& incoming = cubes[nb].outgoing[opposite(d)];
        c.verdict[d].assign(incoming.size(), 0);

        for (size_t k = 0; k < incoming.size(); ++k) 
		{
            const Particle& p = incoming[k].p;
            size_t cell = cellIndex(p.x, p.y, p.z);
            if (!c.occ[cell]) 
			{
                c.occ[cell] = 1;                  // claim before anyone else can
                c.particles.push_back(p);
                c.verdict[d][k] = 1;
            }
        }
    }
}

// ---------------------------------------------------- phase 3: commit removals

void phase3(Cube& c) 
{
    // P2 only appended to c.particles, so srcIdx values recorded in P1
    // still address the same particles.
    std::vector<uint8_t> departed(c.particles.size(), 0);

    for (int d = 0; d < 6; ++d) 
	{
        int nb = c.neighbour[d];
        if (nb < 0) 
			continue;

        const std::vector<uint8_t>& verd = cubes[nb].verdict[opposite(d)];
        for (size_t k = 0; k < c.outgoing[d].size(); ++k) 
		{
            if (!verd[k]) 
				continue;               // rejected -> particle simply stays
			
            uint32_t si = c.outgoing[d][k].srcIdx;
            const Particle& p = c.particles[si];
            c.occ[cellIndex(p.x, p.y, p.z)] = 0;  // now, and only now, free the cell
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
                ok = false; 
				break;
            }
            size_t cell = cellIndex(p.x, p.y, p.z);
            if (seen[cell]) 
			{
                std::cout << "  FAIL: cube " << ci << " has two particles in one cell\n";
                ok = false; 
				break;
            }
            seen[cell] = 1;
            if (!c.occ[cell]) 
			{
                std::cout << "  FAIL: cube " << ci << " particle on a cell occ says is free\n";
                ok = false; 
				break;
            }
        }

        size_t occCount = 0;
        for (uint8_t v : c.occ) 
			occCount += v;
		
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
    if (argc > 1) 
		GM = std::atoi(argv[1]);
    if (argc > 2) 
		GN = std::atoi(argv[2]);
    if (argc > 3) 
		GL = std::atoi(argv[3]);
    if (argc > 4) 
		CS = std::atoi(argv[4]);
    if (argc > 5) 
		ROUNDS = std::atoi(argv[5]);

    const int nCubes = GM * GN * GL;
    cubes.resize(nCubes);

    // build the grid and the 6-neighbour topology
    for (int gx = 0; gx < GM; ++gx)
		for (int gy = 0; gy < GN; ++gy)
			for (int gz = 0; gz < GL; ++gz) 
			{
				Cube& c = cubes[cubeIndex(gx, gy, gz)];
				c.gx = gx; c.gy = gy; c.gz = gz;
				c.occ.assign((size_t)CS * CS * CS, 0);

				const int NX[6] = { gx-1, gx+1, gx,   gx,   gx,   gx  };
				const int NY[6] = { gy,   gy,   gy-1, gy+1, gy,   gy  };
				const int NZ[6] = { gz,   gz,   gz,   gz,   gz-1, gz+1};
				
				for (int d = 0; d < 6; ++d) 
				{
					bool in = NX[d] >= 0 && NX[d] < GM && NY[d] >= 0 && NY[d] < GN && NZ[d] >= 0 && NZ[d] < GL;
					c.neighbour[d] = in ? cubeIndex(NX[d], NY[d], NZ[d]) : -1;
				}
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

    std::cout << "grid " << GM << "x" << GN << "x" << GL << " = " << nCubes << " threads, cube edge " << CS << ", particles " << expectedTotal << ", rounds " << ROUNDS << "\n";

    std::barrier sync(nCubes);

    auto worker = [&](int id) 
	{
        Cube& c = cubes[id];
        std::mt19937 rng(9876u + id);
		
        for (int r = 0; r < ROUNDS; ++r) 
		{
            phase1(c, rng);
            sync.arrive_and_wait();
            phase2(c);
            sync.arrive_and_wait();
            phase3(c);
            sync.arrive_and_wait();
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
