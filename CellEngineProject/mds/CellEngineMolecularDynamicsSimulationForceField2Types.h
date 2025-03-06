
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_2_TYPES_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_2_TYPES_H

#include "../CellEngineTypes.h"

struct Vec3
{
    MDSRealType x, y, z;
    Vec3() : x(0), y(0), z(0)
    {
    }
    Vec3(MDSRealType x, MDSRealType y, MDSRealType z) : x(x), y(y), z(z)
    {
    }
    Vec3 operator+(const Vec3& other) const
    {
        return { x + other.x, y + other.y, z + other.z };
    }
    Vec3 operator-(const Vec3& other) const
    {
        return{ x - other.x, y - other.y, z - other.z };
    }
    Vec3 operator*(MDSRealType s) const
    {
        return { x * s, y * s, z * s };
    }
    Vec3 operator/(MDSRealType s) const
    {
        return { x / s, y / s, z / s };
    }
    Vec3& operator+=(const Vec3& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;

        return *this;
    }
    Vec3& operator-=(const Vec3& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;

        return *this;
    }
    Vec3& operator*=(const Vec3& other)
    {
        x *= other.x;
        y *= other.y;
        z *= other.z;

        return *this;
    }
    Vec3& operator*=(const MDSRealType s)
    {
        x *= s;
        y *= s;
        z *= s;

        return *this;
    }
    [[nodiscard]] MDSRealType dot(const Vec3& other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }
    [[nodiscard]] Vec3 cross(const Vec3& other) const
    {
        return { y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x };
    }
    [[nodiscard]] MDSRealType length() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }
    [[nodiscard]] Vec3 normalized() const
    {
        MDSRealType l = length();

        return { x / l, y / l, z / l };
    }
};

enum class AtomTypes : int
{
    Oxygen = 0,
    Nitrogen = 1,
    Phosphorous = 2,
    Carbon = 3
};

class Atom
{
public:
    AtomTypes AtomType;
    Vec3 Position, Velocity, Force;
    MDSRealType Charge, Mass, SigmaValue, EpsilonValue;
public:
    Atom(const AtomTypes AtomTypeParam, const Vec3& pos, const MDSRealType charge) : AtomType(AtomTypeParam), Position(pos), Charge(charge)
    {
        switch(AtomType)
        {
            case AtomTypes::Oxygen: Mass = 16.0; SigmaValue = 3.0; EpsilonValue = 0.2; break;
            case AtomTypes::Nitrogen: Mass = 14.0; SigmaValue = 3.2; EpsilonValue = 0.1; break;
            case AtomTypes::Phosphorous: Mass = 31.0; SigmaValue = 3.5; EpsilonValue = 0.15; break;
            case AtomTypes::Carbon: Mass = 12.0; SigmaValue = 3.4; EpsilonValue = 0.05; break;
            default: Mass = 12.0; SigmaValue = 3.0; EpsilonValue = 0.1;
        }
    }
};

struct Bond
{
public:
    UnsignedInt i, j;
    MDSRealType r0; // equilibrium distance
    MDSRealType k; // spring constant
public:
    Bond(const UnsignedInt i, const UnsignedInt j, const MDSRealType r0, const MDSRealType k) : i(i), j(j), r0(r0), k(k)
    {
    }
};

struct Angle
{
public:
    UnsignedInt i, j, k;
    MDSRealType theta0; // equilibrium angle in radians
    MDSRealType k_theta; // force constant
public:
    Angle(const UnsignedInt i, const UnsignedInt j, const UnsignedInt k, const MDSRealType theta0, const MDSRealType k_theta) : i(i), j(j), k(k), theta0(theta0), k_theta(k_theta)
    {
    }
};

struct Dihedral
{
public:
    UnsignedInt i, j, k, l;
    MDSRealType k_phi; // force constant
    UnsignedInt n; // multiplicity
    MDSRealType phi0; // phase angle in radians
public:
    Dihedral(const UnsignedInt i, const UnsignedInt j, const UnsignedInt k, const UnsignedInt l, const MDSRealType k_phi, const UnsignedInt n, const MDSRealType phi0) : i(i), j(j), k(k), l(l), k_phi(k_phi), n(n), phi0(phi0)
    {
    }
};

#endif
