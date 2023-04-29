#pragma once

#ifndef CELL_ENGINE_TYPES_H
#define CELL_ENGINE_TYPES_H

#include <cstdint>
#include "vmath.h"

using SignedInt = std::int64_t;
using UnsignedInt = std::uint64_t;

using EntityIdInt = std::uint16_t;
using ChainIdInt = std::uint16_t;
using UniqueIdInt = std::uint32_t;

template<class T>
struct vector3
{
    T X, Y, Z;
};

using vector3_16 = vector3<uint16_t>;
using vector3_64 = vector3<uint64_t>;

inline vector3_16 GetVector3FormVMathVec3(const vmath::vec3& Color)
{
    return {static_cast<uint16_t>(Color.X() * 100.00f), static_cast<uint16_t>(Color.Y() * 100.00f), static_cast<uint16_t>(Color.Z() * 100.00f) };
}

inline vmath::vec3 GetVMathVec3FromVector3(vector3_16 Color)
{
    return { float(Color.X) / 100.00f, float(Color.Y) / 100.00f, float(Color.Z) / 100.00f };
}

#endif