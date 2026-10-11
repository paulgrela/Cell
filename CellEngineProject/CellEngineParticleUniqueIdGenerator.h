
#ifndef CELL_ENGINE_PARTICLE_UNIQUE_ID_GENERATOR_H
#define CELL_ENGINE_PARTICLE_UNIQUE_ID_GENERATOR_H

#include <algorithm>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <vector>

class CellEngineParticleUniqueIdGenerator
{
private:
    std::uint64_t Prefix{ 0 };
    UniqueIdUInt64Type LastIssuedCounter{ 0 };

    static inline unsigned SectorBits{ 16 };
    static inline unsigned CounterBits{ 48 };
    static inline UniqueIdUInt64Type LastCounterValue{ (std::uint64_t{ 1 } << 48) - 1 };

public:
    static void ConfigureForNumberOfSectors(const UnsignedInt NumberOfAllSectors)
    {
        SectorBits = std::max(1u, static_cast<unsigned>(std::bit_width(NumberOfAllSectors - 1)));
        CounterBits = 64u - SectorBits;
        LastCounterValue = (std::uint64_t{ 1 } << CounterBits) - 1;
    }

    void Initialize(const UnsignedInt GlobalSectorLinearIndex)
    {
        Prefix = GlobalSectorLinearIndex << CounterBits;
        LastIssuedCounter = 0;
    }

    [[nodiscard]] UniqueIdUInt64Type GetNewParticleUniqueId()
    {
        if (LastIssuedCounter == LastCounterValue) [[unlikely]]
            throw std::overflow_error("ParticleUniqueIdGenerator: ID space of this sector is exhausted");

        return Prefix | ++LastIssuedCounter;
    }

    static void ReleaseParticleUniqueId(const UniqueIdUInt64Type) noexcept {}

    [[nodiscard]] UniqueIdUInt64Type SaveState() const noexcept
    {
        return LastIssuedCounter;
    }
    void RestoreState(const UniqueIdUInt64Type SavedLastIssuedCounter) noexcept
    {
        LastIssuedCounter = SavedLastIssuedCounter;
    }

    [[nodiscard]] static UniqueIdUInt64Type CreatorSectorOf(const UniqueIdUInt64Type Id) noexcept
    {
        return Id >> CounterBits;
    }
    [[nodiscard]] static UniqueIdUInt64Type LocalCounterOf(const UniqueIdUInt64Type Id) noexcept
    {
        return Id & LastCounterValue;
    }
    [[nodiscard]] static unsigned GetSectorBits() noexcept
    {
        return SectorBits;
    }
    [[nodiscard]] static unsigned GetCounterBits() noexcept
    {
        return CounterBits;
    }
};

class CellEngineParticleUniqueIdGeneratorWithReuse
{
private:
    UniqueIdUInt64Type NextFreshId{ 1 };
    UniqueIdUInt64Type EndOfRange{ 1 };
    std::vector<UniqueIdUInt64Type> FreeList{};

public:
    void Initialize(const UniqueIdUInt64Type GlobalSectorLinearIndex, UniqueIdUInt64Type IdsPerSector)
    {
        NextFreshId = GlobalSectorLinearIndex * IdsPerSector + 1;
        EndOfRange = (GlobalSectorLinearIndex + 1) * IdsPerSector;
        FreeList.clear();
    }

    [[nodiscard]] UniqueIdUInt64Type GetNewParticleUniqueId()
    {
        if (!FreeList.empty())
        {
            const UniqueIdUInt64Type Id = FreeList.back();
            FreeList.pop_back();
            return Id;
        }

        if (NextFreshId == EndOfRange) [[unlikely]]
            throw std::overflow_error("ParticleUniqueIdGeneratorWithReuse: fresh range of this sector is exhausted");

        return NextFreshId++;
    }

    void ReleaseUniqueParticleUniqueId(const UniqueIdUInt64Type Id)
    {
        FreeList.push_back(Id);
    }

    [[nodiscard]] std::size_t NumberOfReleasedIdsWaiting() const noexcept
    {
        return FreeList.size();
    }
};

#endif