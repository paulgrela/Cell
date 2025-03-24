#pragma once

#ifndef CELL_ENGINE_TYPES_H
#define CELL_ENGINE_TYPES_H

#include <stack>
#include <memory>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "vmath.h"

using SignedInt = std::int64_t;
using UnsignedInt = std::uint64_t;

using GeneIdInt = std::uint64_t;
using EntityIdInt = std::uint32_t;
using ChainIdInt = std::uint16_t;
using PositionInt = std::uint16_t;
using UniqueIdInt = std::uint32_t;

using RealType = float;

using MDSRealType = double;

using ElectricChargeType = std::int32_t;

template<class T>
struct vector3
{
    T X, Y, Z;

    bool operator==(const vector3& CTP) const
    {
        return (CTP.X == X && CTP.Y == Y && CTP.Z == Z);
    }
};

using vector3_16 = vector3<PositionInt>;
using vector3_64 = vector3<uint64_t>;
using vector3_Real32 = vector3<RealType>;

struct ThreadPosType
{
    UnsignedInt	ThreadPosX, ThreadPosY, ThreadPosZ;

    bool operator==(const ThreadPosType& CTP) const
    {
        return (CTP.ThreadPosX == ThreadPosX && CTP.ThreadPosY == ThreadPosY && CTP.ThreadPosZ == ThreadPosZ);
    }
};

struct SectorPosType
{
    SignedInt SectorPosX, SectorPosY, SectorPosZ;

    bool operator==(const SectorPosType& CSP) const
    {
        return (CSP.SectorPosX == SectorPosX && CSP.SectorPosY == SectorPosY && CSP.SectorPosZ == SectorPosZ);
    }
};

struct PosType
{
    float PosX, PosY, PosZ;

    bool operator==(const PosType& PT) const
    {
        return (PT.PosX == PosX && PT.PosY == PosY && PT.PosZ == PosZ);
    }
};

using ThreadIdType = std::uint64_t;

template<class Particle>
using ParticlesDetailedContainer = std::unordered_map<UniqueIdInt, Particle>;

template <class Particle>
struct ParticlesContainerInternal
{
public:
    UnsignedInt MPIProcessIndex{ 0 };
    ThreadPosType ThreadPos{ 1, 1, 1 };
public:
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
    ParticlesDetailedContainer<Particle> Particles;
};

template <class Particle>
using ParticlesContainer = std::vector<std::vector<std::vector<ParticlesContainerInternal<Particle>>>>;

template <class SimulationSpace>
using SimulationSpaceForParallelExecutionContainer = std::vector<std::vector<std::vector<std::shared_ptr<SimulationSpace>>>>;

struct SimulationSpaceSectorsRanges
{
public:
    UnsignedInt StartXPos;
    UnsignedInt StartYPos;
    UnsignedInt StartZPos;
    UnsignedInt EndXPos;
    UnsignedInt EndYPos;
    UnsignedInt EndZPos;
public:
    void SetParameters(const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float EndXPosParam, const float EndYPosParam, const float EndZPosParam)
    {
        StartXPos = static_cast<UnsignedInt>(StartXPosParam);
        StartYPos = static_cast<UnsignedInt>(StartYPosParam);
        StartZPos = static_cast<UnsignedInt>(StartZPosParam);
        EndXPos = static_cast<UnsignedInt>(EndXPosParam);
        EndYPos = static_cast<UnsignedInt>(EndYPosParam);
        EndZPos = static_cast<UnsignedInt>(EndZPosParam);
    }
};

struct SimulationSpaceSectorBounds
{
public:
    float StartXPos;
    float StartYPos;
    float StartZPos;
    float SizeX;
    float SizeY;
    float SizeZ;
    float EndXPos;
    float EndYPos;
    float EndZPos;
public:
    void SetParameters(const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam, const float EndXPosParam, const float EndYPosParam, const float EndZPosParam)
    {
        StartXPos = StartXPosParam;
        StartYPos = StartYPosParam;
        StartZPos = StartZPosParam;
        SizeX = SizeXParam;
        SizeY = SizeYParam;
        SizeZ = SizeZParam;
        EndXPos = EndXPosParam;
        EndYPos = EndYPosParam;
        EndZPos = EndZPosParam;
    }
public:
    void AddToStartParameters(const float AddToStartXPosParam, const float AddToStartYPosParam, const float AddToStartZPosParam)
    {
        StartXPos += AddToStartXPosParam;
        StartYPos += AddToStartYPosParam;
        StartZPos += AddToStartZPosParam;
    }
public:
    SimulationSpaceSectorBounds SetParametersForParallelExecutionSectors(const ThreadPosType& CurrentThreadPos, const UnsignedInt SizeOfXInOneThreadInSimulationSpace, const UnsignedInt SizeOfYInOneThreadInSimulationSpace, const UnsignedInt SizeOfZInOneThreadInSimulationSpace)
    {
        StartXPos = static_cast<float>((CurrentThreadPos.ThreadPosX - 1) * SizeOfXInOneThreadInSimulationSpace);
        EndXPos = static_cast<float>((CurrentThreadPos.ThreadPosX - 1) * SizeOfXInOneThreadInSimulationSpace + SizeOfXInOneThreadInSimulationSpace);
        StartYPos = static_cast<float>((CurrentThreadPos.ThreadPosY - 1) * SizeOfYInOneThreadInSimulationSpace);
        EndYPos = static_cast<float>((CurrentThreadPos.ThreadPosY - 1) * SizeOfYInOneThreadInSimulationSpace + SizeOfYInOneThreadInSimulationSpace);
        StartZPos = static_cast<float>((CurrentThreadPos.ThreadPosZ - 1) * SizeOfZInOneThreadInSimulationSpace);
        EndZPos = static_cast<float>((CurrentThreadPos.ThreadPosZ - 1) * SizeOfZInOneThreadInSimulationSpace + SizeOfZInOneThreadInSimulationSpace);
        SizeX = EndXPos - StartXPos;
        SizeY = EndYPos - StartYPos;
        SizeZ = EndZPos - StartZPos;

        return *this;
    }
public:
    SimulationSpaceSectorBounds SetParametersForChosenSector(const UnsignedInt SectorXIndex, const UnsignedInt SectorYIndex, const UnsignedInt SectorZIndex, const RealType ShiftCenterX, const RealType ShiftCenterY, const RealType ShiftCenterZ, const RealType SizeOfParticlesSectorX, const RealType SizeOfParticlesSectorY, const RealType SizeOfParticlesSectorZ)
    {
        StartXPos = - ShiftCenterX + static_cast<RealType>(SectorXIndex) * SizeOfParticlesSectorX;
        StartYPos = - ShiftCenterY + static_cast<RealType>(SectorYIndex) * SizeOfParticlesSectorY;
        StartZPos = - ShiftCenterZ + static_cast<RealType>(SectorZIndex) * SizeOfParticlesSectorZ;
        SizeX = SizeOfParticlesSectorX;
        SizeY = SizeOfParticlesSectorY;
        SizeZ = SizeOfParticlesSectorZ;
        EndXPos = StartXPos + SizeX;
        EndYPos = StartYPos + SizeY;
        EndZPos = StartZPos + SizeZ;

        return *this;
    }
};

struct __attribute__ ((packed)) MPIParticleSenderStruct
{
    UniqueIdInt ParticleIndex{ 0 };
    EntityIdInt ParticleKindId{ 0 };
    int SenderProcessIndex{ 0 };
    int ReceiverProcessIndex{ 0 };
    vector3_16 SectorPos{ 0, 0, 0 };
    vector3_Real32 NewPosition{ 0, 0, 0 };
};

enum class TypesOfLookingForParticlesInProximity : UnsignedInt
{
    FromChosenParticleAsCenter = 1,
    InChosenSectorOfSimulationSpace = 2
};

template <typename MutexT>
class conditional_lock_guard
{
private:
    bool condition;
    MutexT* const mtx;
public:
    explicit conditional_lock_guard(const bool condition_param, MutexT* mtx_param) : condition(condition_param), mtx{ mtx_param }
    {
        if (condition == true)
        {
            mtx->lock();
        }
    }
public:
    ~conditional_lock_guard() noexcept
    {
        if (condition == true)
        {
            mtx->unlock();
        }
    }
public:
    conditional_lock_guard(const conditional_lock_guard&) = delete;
    conditional_lock_guard& operator=(const conditional_lock_guard&) = delete;
public:
    [[nodiscard]] bool owns_lock() const noexcept
    {
        return condition;
    }
public:
    explicit operator bool() const noexcept
    {
        return owns_lock();
    }
public:
    MutexT* mutex() const noexcept
    {
        return mtx;
    }
};

#endif