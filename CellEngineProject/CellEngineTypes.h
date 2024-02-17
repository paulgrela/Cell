#pragma once

#ifndef CELL_ENGINE_TYPES_H
#define CELL_ENGINE_TYPES_H

#include <cstdint>
#include "vmath.h"

using SignedInt = std::int64_t;
using UnsignedInt = std::uint64_t;

using EntityIdInt = std::uint32_t;
using ChainIdInt = std::uint16_t;
using UniqueIdInt = std::uint32_t;

using ElectricChargeType = std::int16_t;

template<class T>
struct vector3
{
    T X, Y, Z;
};

using vector3_16 = vector3<uint16_t>;
using vector3_64 = vector3<uint64_t>;

enum class TypesOfLookingForParticlesInProximity : UnsignedInt
{
    FromChosenParticleAsCenter = 1,
    InChosenVoxelSpace = 2
};

#endif