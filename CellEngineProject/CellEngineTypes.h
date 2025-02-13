#pragma once

#ifndef CELL_ENGINE_TYPES_H
#define CELL_ENGINE_TYPES_H

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

using ElectricChargeType = std::int32_t;

template<class T>
struct vector3
{
    T X, Y, Z;
};

using vector3_16 = vector3<PositionInt>;
using vector3_64 = vector3<uint64_t>;
using vector3_Real32 = vector3<RealType>;

struct CurrentThreadPosType
{
    UnsignedInt	ThreadPosX, ThreadPosY, ThreadPosZ;

    bool operator==(const CurrentThreadPosType& CTP) const
    {
        return (CTP.ThreadPosX == ThreadPosX && CTP.ThreadPosY == ThreadPosY && CTP.ThreadPosZ == ThreadPosZ);
    }
};

struct CurrentSectorPosType
{
    UnsignedInt	SectorPosX, SectorPosY, SectorPosZ;

    bool operator==(const CurrentSectorPosType& CSP) const
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

template <class Particle>
using ParticlesContainerInternal = std::unordered_map<UniqueIdInt, Particle>;

template <class Particle>
using ParticlesContainer = std::vector<std::vector<std::vector<ParticlesContainerInternal<Particle>>>>;

template <class SimulationSpace>
using SimulationSpaceForParallelExecutionContainer = std::vector<std::vector<std::vector<std::shared_ptr<SimulationSpace>>>>;

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
    SimulationSpaceSectorBounds SetParametersForParallelExecutionSectors(const CurrentThreadPosType& CurrentThreadPos, const UnsignedInt SizeOfXInOneThreadInSimulationSpace, const UnsignedInt SizeOfYInOneThreadInSimulationSpace, const UnsignedInt SizeOfZInOneThreadInSimulationSpace)
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