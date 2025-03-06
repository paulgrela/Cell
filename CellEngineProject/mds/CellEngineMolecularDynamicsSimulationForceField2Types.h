
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_2_TYPES_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_2_TYPES_H

#include "../CellEngineTypes.h"

struct Vec3
{
public:
    MDSRealType x, y, z;
public:
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

class AtomMDS
{
public:
    AtomTypes AtomType;
    Vec3 Position, Velocity, Force;
    MDSRealType Charge, Mass, SigmaValue, EpsilonValue;
public:
    AtomMDS(const AtomTypes AtomTypeParam, const Vec3& pos, const MDSRealType charge) : AtomType(AtomTypeParam), Position(pos), Charge(charge)
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

struct BondMDS
{
public:
    UnsignedInt Atom1Index, Atom2Index;
    MDSRealType r0_EquilibrumDistance;
    MDSRealType k_SpringConstant;
public:
    BondMDS(const UnsignedInt Atom1IndexParam, const UnsignedInt Atom2IndexParam, const MDSRealType r0_EquilibrumDistanceParam, const MDSRealType k_SpringConstantParam) : Atom1Index(Atom1IndexParam), Atom2Index(Atom2IndexParam), r0_EquilibrumDistance(r0_EquilibrumDistanceParam), k_SpringConstant(k_SpringConstantParam)
    {
    }
};

struct AngleMDS
{
public:
    UnsignedInt Atom1Index, Atom2Index, Atom3Index;
    MDSRealType theta0_EquilibriumAngleInRadians;
    MDSRealType k_theta_ForceConstant;
public:
    AngleMDS(const UnsignedInt Atom1IndexParam, const UnsignedInt Atom2IndexParam, const UnsignedInt Atom3IndexParam, const MDSRealType theta0, const MDSRealType k_theta) : Atom1Index(Atom1IndexParam), Atom2Index(Atom2IndexParam), Atom3Index(Atom3IndexParam), theta0_EquilibriumAngleInRadians(theta0), k_theta_ForceConstant(k_theta)
    {
    }
};

struct DihedralMDS
{
public:
    UnsignedInt Atom1Index, Atom2Index, Atom3Index, Atom4Index;
    MDSRealType k_phi_ForceConstant;
    UnsignedInt N_Multiplicity;
    MDSRealType phi0_PhaseAngleInRadians;
public:
    DihedralMDS(const UnsignedInt Atom1IndexParam, const UnsignedInt Atom2IndexParam, const UnsignedInt Atom3IndexParam,  const UnsignedInt Atom4IndexParam, const MDSRealType k_phi, const UnsignedInt N_MultiplicityParam, const MDSRealType phi0) : Atom1Index(Atom1IndexParam), Atom2Index(Atom2IndexParam), Atom3Index(Atom3IndexParam), Atom4Index(Atom4IndexParam), k_phi_ForceConstant(k_phi), N_Multiplicity(N_MultiplicityParam), phi0_PhaseAngleInRadians(phi0)
    {
    }
};

#endif
