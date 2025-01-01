#pragma once

#ifndef CELL_ENGINE_TYPES_H
#define CELL_ENGINE_TYPES_H

#include <cstdint>
#include "vmath.h"

using SignedInt = std::int64_t;
using UnsignedInt = std::uint64_t;

using GeneIdInt = std::uint64_t;
using EntityIdInt = std::uint32_t;
using ChainIdInt = std::uint16_t;
using PositionInt = std::uint16_t;
using UniqueIdInt = std::uint32_t;

using ElectricChargeType = std::int32_t;

template<class T>
struct vector3
{
    T X, Y, Z;
};

using vector3_16 = vector3<PositionInt>;
using vector3_64 = vector3<uint64_t>;

struct CurrentThreadPosType
{
    UnsignedInt	ThreadPosX, ThreadPosY, ThreadPosZ;

    bool operator==(const CurrentThreadPosType& CTP) const
    {
        return (CTP.ThreadPosX == ThreadPosX && CTP.ThreadPosY == ThreadPosY && CTP.ThreadPosZ == ThreadPosZ);
    }
};

struct SimulationSpaceSectorBounds
{
    UnsignedInt StartXPos;
    UnsignedInt StartYPos;
    UnsignedInt StartZPos;
    UnsignedInt SizeX;
    UnsignedInt SizeY;
    UnsignedInt SizeZ;
    UnsignedInt EndXPos;
    UnsignedInt EndYPos;
    UnsignedInt EndZPos;
};

enum class TypesOfLookingForParticlesInProximity : UnsignedInt
{
    FromChosenParticleAsCenter = 1,
    InChosenSectorOfSimulationSpace = 2
};

#endif