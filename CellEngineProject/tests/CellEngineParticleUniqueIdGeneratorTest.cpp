
#include <algorithm>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#include "CellEngineParticleUniqueIdGenerator.h"

namespace
{
    constexpr int SectorsX = 40, SectorsY = 40, SectorsZ = 40;
    constexpr std::uint64_t NumberOfAllSectors = std::uint64_t{ SectorsX } * SectorsY * SectorsZ;
    constexpr int InitialParticlesPerSector = 16;
    constexpr int ReactionsPerSectorPerRound = 4;
    constexpr int Rounds = 20;
    constexpr unsigned NumberOfThreads = 8;

    std::uint64_t GlobalSectorLinearIndex(const int X, const int Y, const int Z)
    {
        return (std::uint64_t(X) * SectorsY + Y) * SectorsZ + Z;
    }

    int Failures = 0;
    void Check(const bool Condition, const char* What)
    {
        std::printf("  [%s] %s\n", Condition ? " OK " : "FAIL", What);
        if (!Condition)
            ++Failures;
    }

    template <class Generator>
    struct SectorData
    {
        Generator IdGenerator;
        std::vector<std::uint64_t> LiveParticles;
    };

    struct RunResult
    {
        std::size_t LiveParticles = 0, IssuedIds = 0, Reactions = 0, ReissuedIds = 0;
        bool LiveUnique = false, NoZeroId = false, CountConserved = false;
    };

    template <class Generator, class InitializeFunction>
    RunResult RunSimulation(InitializeFunction InitializeGenerator)
    {
        std::vector<SectorData<Generator>> Sectors(NumberOfAllSectors);
        for (std::uint64_t SectorIndex = 0; SectorIndex < NumberOfAllSectors; ++SectorIndex)
            InitializeGenerator(Sectors[SectorIndex].IdGenerator, SectorIndex);

        // start-up: each sector creates its particles from its OWN generator — no contains() scan
        std::vector<std::uint64_t> AllIssued;
        for (auto& Sector : Sectors)
            for (int k = 0; k < InitialParticlesPerSector; ++k)
            {
                const std::uint64_t Id = Sector.IdGenerator.GetNewUniqueParticleIndex();
                Sector.LiveParticles.push_back(Id);
                AllIssued.push_back(Id);
            }
        const std::size_t InitialParticles = AllIssued.size();

        std::vector<std::vector<std::uint64_t>> IssuedByThread(NumberOfThreads);
        std::vector<std::size_t> CreatedByThread(NumberOfThreads), DestroyedByThread(NumberOfThreads), ReactionsByThread(NumberOfThreads);

        constexpr int Shapes[5][2] = { { 2, 2 }, { 2, 1 }, { 1, 2 }, { 1, 0 }, { 0, 1 } };   // { reactants, products }

        for (int Round = 0; Round < Rounds; ++Round)
        {
            // REACTIONS, in parallel. Each thread owns a slab of x-planes and is the
            // only one touching those sectors' generators: no atomics, no locks.
            {
                std::vector<std::jthread> Threads;
                for (unsigned t = 0; t < NumberOfThreads; ++t)
                    Threads.emplace_back([&, t]
                    {
                        std::mt19937_64 Rng(1000u * unsigned(Round) + t);
                        const int XBegin = int(SectorsX * t / NumberOfThreads);
                        const int XEnd = int(SectorsX * (t + 1) / NumberOfThreads);

                        for (int X = XBegin; X < XEnd; ++X)
                            for (int Y = 0; Y < SectorsY; ++Y)
                                for (int Z = 0; Z < SectorsZ; ++Z)
                                {
                                    auto& Sector = Sectors[GlobalSectorLinearIndex(X, Y, Z)];
                                    for (int r = 0; r < ReactionsPerSectorPerRound; ++r)
                                    {
                                        const auto& Shape = Shapes[Rng() % 5];
                                        if (Sector.LiveParticles.size() < std::size_t(Shape[0]))
                                            continue;

                                        for (int a = 0; a < Shape[0]; ++a)          // reactants disappear
                                        {
                                            const std::size_t j = Rng() % Sector.LiveParticles.size();
                                            Sector.IdGenerator.ReleaseUniqueParticleIndex(Sector.LiveParticles[j]);
                                            Sector.LiveParticles[j] = Sector.LiveParticles.back();
                                            Sector.LiveParticles.pop_back();
                                            ++DestroyedByThread[t];
                                        }
                                        for (int p = 0; p < Shape[1]; ++p)          // products appear
                                        {
                                            const std::uint64_t Id = Sector.IdGenerator.GetNewUniqueParticleIndex();
                                            Sector.LiveParticles.push_back(Id);
                                            IssuedByThread[t].push_back(Id);
                                            ++CreatedByThread[t];
                                        }
                                        ++ReactionsByThread[t];
                                    }
                                }
                    });
            }   // std::jthread joins here

            // DIFFUSION: ~10% of particles jump to a face neighbour, often one owned
            // by another thread. They keep their ID.
            std::mt19937_64 Rng(777u + unsigned(Round));
            constexpr int DX[6] = { -1, 0, 0, 1, 0, 0 }, DY[6] = { 0, -1, 0, 0, 1, 0 }, DZ[6] = { 0, 0, -1, 0, 0, 1 };
            for (int X = 0; X < SectorsX; ++X)
                for (int Y = 0; Y < SectorsY; ++Y)
                    for (int Z = 0; Z < SectorsZ; ++Z)
                    {
                        auto& Live = Sectors[GlobalSectorLinearIndex(X, Y, Z)].LiveParticles;
                        for (std::size_t i = Live.size(); i-- > 0;)
                        {
                            if (Rng() % 10 != 0)
                                continue;
                            const int d = int(Rng() % 6);
                            const int NX = X + DX[d], NY = Y + DY[d], NZ = Z + DZ[d];
                            if (NX < 0 || NX >= SectorsX || NY < 0 || NY >= SectorsY || NZ < 0 || NZ >= SectorsZ)
                                continue;
                            Sectors[GlobalSectorLinearIndex(NX, NY, NZ)].LiveParticles.push_back(Live[i]);
                            Live[i] = Live.back();
                            Live.pop_back();
                        }
                    }
        }

        RunResult Result;

        std::vector<std::uint64_t> Live;
        for (const auto& Sector : Sectors)
            Live.insert(Live.end(), Sector.LiveParticles.begin(), Sector.LiveParticles.end());
        std::sort(Live.begin(), Live.end());
        Result.LiveParticles = Live.size();
        Result.LiveUnique = std::adjacent_find(Live.begin(), Live.end()) == Live.end();
        Result.NoZeroId = Live.empty() || Live.front() != 0;

        std::size_t Created = 0, Destroyed = 0;
        for (unsigned t = 0; t < NumberOfThreads; ++t)
        {
            Created += CreatedByThread[t];
            Destroyed += DestroyedByThread[t];
            Result.Reactions += ReactionsByThread[t];
            AllIssued.insert(AllIssued.end(), IssuedByThread[t].begin(), IssuedByThread[t].end());
        }
        Result.CountConserved = (Live.size() == InitialParticles + Created - Destroyed);

        Result.IssuedIds = AllIssued.size();
        std::sort(AllIssued.begin(), AllIssued.end());
        const auto DistinctEnd = std::unique(AllIssued.begin(), AllIssued.end());
        Result.ReissuedIds = AllIssued.size() - std::size_t(DistinctEnd - AllIssued.begin());
        return Result;
    }
}

int main()
{
    using Gen = ParticleUniqueIdGenerator;

    std::printf("1. ID layout\n");
    Gen::ConfigureForNumberOfSectors(NumberOfAllSectors);
    Check(Gen::GetSectorBits() == 16 && Gen::GetCounterBits() == 48, "64000 sectors -> 16 sector bits + 48 counter bits");

    Gen Sector0, Sector12345, SectorLast;
    Sector0.Initialize(0);
    Sector12345.Initialize(12345);
    SectorLast.Initialize(NumberOfAllSectors - 1);
    Check(Sector0.GetNewUniqueParticleIndex() == 1, "first ID of sector 0 is 1, so ID 0 is never issued");
    const auto Id = Sector12345.GetNewUniqueParticleIndex();
    Check(Gen::CreatorSectorOf(Id) == 12345 && Gen::LocalCounterOf(Id) == 1, "CreatorSectorOf / LocalCounterOf decode an ID");
    Check(Gen::CreatorSectorOf(SectorLast.GetNewUniqueParticleIndex()) == NumberOfAllSectors - 1, "prefix of the last sector fits");

    Gen Nearly;
    Nearly.Initialize(7);
    Nearly.RestoreState((std::uint64_t{ 1 } << 48) - 2);
    const auto LastId = Nearly.GetNewUniqueParticleIndex();
    bool Threw = false;
    try { (void)Nearly.GetNewUniqueParticleIndex(); }
    catch (const std::overflow_error&) { Threw = true; }
    Check(Gen::CreatorSectorOf(LastId) == 7 && Threw, "last counter value is issued, the next call throws instead of spilling into sector 8");

    Gen::ConfigureForNumberOfSectors(1);
    Check(Gen::GetCounterBits() == 63, "1 sector: counter gets 63 bits, no undefined 64-bit shift");
    Gen::ConfigureForNumberOfSectors(65536);
    const bool SixteenBits = (Gen::GetSectorBits() == 16);
    Gen::ConfigureForNumberOfSectors(65537);
    Check(SixteenBits && Gen::GetSectorBits() == 17, "65536 sectors -> 16 bits, 65537 -> 17 bits");

    Gen::ConfigureForNumberOfSectors(NumberOfAllSectors);   // back to your layout

    std::printf("\n2. Simulation WITHOUT reuse: %u threads, %d rounds, reactions 2->2 2->1 1->2 1->0 0->1, diffusion between sectors\n", NumberOfThreads, Rounds);
    const auto A = RunSimulation<Gen>([](Gen& G, const std::uint64_t SectorIndex) { G.Initialize(SectorIndex); });
    std::printf("  %zu reactions, %zu IDs issued, %zu particles alive at the end\n", A.Reactions, A.IssuedIds, A.LiveParticles);
    Check(A.LiveUnique && A.NoZeroId, "all live IDs unique, none is 0");
    Check(A.ReissuedIds == 0, "no ID was ever issued twice");
    Check(A.CountConserved, "live count == initial + created - destroyed");

    std::printf("\n3. Same simulation WITH reuse (bump pointer + free list, your 100000-wide ranges)\n");
    const auto B = RunSimulation<ParticleUniqueIdGeneratorWithReuse>([](ParticleUniqueIdGeneratorWithReuse& G, const std::uint64_t SectorIndex) { G.Initialize(SectorIndex, 100000); });
    std::printf("  %zu reactions, %zu IDs issued, %zu of them reused, %zu particles alive at the end\n", B.Reactions, B.IssuedIds, B.ReissuedIds, B.LiveParticles);
    Check(B.LiveUnique && B.NoZeroId, "all live IDs unique, none is 0");
    Check(B.ReissuedIds > 0, "released IDs really were handed out again");
    Check(B.CountConserved, "live count == initial + created - destroyed");

    std::printf("\n%s\n", Failures == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED");
    return Failures == 0 ? 0 : 1;
}
