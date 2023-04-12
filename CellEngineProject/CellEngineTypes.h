#pragma once

#ifndef _CELL_ENGINE_TYPES_H_
#define _CELL_ENGINE_TYPES_H_

#include <cstdint>
#include "vmath.h"

using Int = std::int64_t;
using UnsignedInt = std::uint64_t;

struct vector3
{
    uint16_t X, Y, Z;
};

inline vector3 GetVector3FormVMathVec3(vmath::vec3 Color)
{
    return {static_cast<uint16_t>(Color.X() * 100.00f), static_cast<uint16_t>(Color.Y() * 100.00f), static_cast<uint16_t>(Color.Z() * 100.00f) };
}

inline vmath::vec3 GetVMathVec3FromVector3(vector3 Color)
{
    return { float(Color.X) / 100.00f, float(Color.Y) / 100.00f, float(Color.Z) / 100.00f };
}

#endif