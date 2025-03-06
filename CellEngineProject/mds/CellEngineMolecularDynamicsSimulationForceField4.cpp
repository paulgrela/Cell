
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include "../CellEngineTypes.h"
#include "CellEngineMolecularDynamicsSimulationForceField4.h"
#include "CellEngineMolecularDynamicsSimulationForceField4Types.h"
#include "CellEngineMolecularDynamicsSimulationForceField4Constants.h"

using namespace std;

void ComputeLennardJonesForce(const AtomMDS& a, const AtomMDS& b, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = a.PositionX - b.PositionX, dy = a.PositionY - b.PositionY, dz = a.PositionZ - b.PositionZ;
    const MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    const MDSRealType r = sqrt(r2);
    const MDSRealType r6 = pow(sigma / r, 6);
    const MDSRealType r12 = r6 * r6;
    const MDSRealType force = 24 * epsilon * (2 * r12 - r6) / r;

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeCoulombForce(const AtomMDS& a, const AtomMDS& b, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = a.PositionX - b.PositionX, dy = a.PositionY - b.PositionY, dz = a.PositionZ - b.PositionZ;
    const MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    const MDSRealType r = sqrt(r2);
    const MDSRealType force = (a.Charge * b.Charge) / (4 * pi * epsilon0 * r * r);

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeBondedHookesLawForce(const AtomMDS& a, const AtomMDS& b, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = a.PositionX - b.PositionX, dy = a.PositionY - b.PositionY, dz = a.PositionZ - b.PositionZ;
    const MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
    const MDSRealType force = -k_bond * (r - r0);

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeAngleBendingForce(const AtomMDS& a, const AtomMDS& b, const AtomMDS& c, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType abx = a.PositionX - b.PositionX, aby = a.PositionY - b.PositionY, abz = a.PositionZ - b.PositionZ;
    const MDSRealType cbx = c.PositionX - b.PositionX, cby = c.PositionY - b.PositionY, cbz = c.PositionZ - b.PositionZ;
    const MDSRealType dot = abx * cbx + aby * cby + abz * cbz;
    const MDSRealType mag_ab = sqrt(abx * abx + aby * aby + abz * abz);
    const MDSRealType mag_cb = sqrt(cbx * cbx + cby * cby + cbz * cbz);
    const MDSRealType theta = acos(dot / (mag_ab * mag_cb));
    const MDSRealType force = -k_angle * (theta - theta0);

    const MDSRealType dtheta_dr1 = (cos(theta) * abx - cbx) / (mag_ab * mag_cb);
    const MDSRealType dtheta_dr3 = (cos(theta) * cbx - abx) / (mag_ab * mag_cb);

    ForceX += force * dtheta_dr1;
    ForceY += force * dtheta_dr1;
    ForceZ += force * dtheta_dr1;
    ForceX += force * dtheta_dr3;
    ForceY += force * dtheta_dr3;
    ForceZ += force * dtheta_dr3;
}

void CalculateBondStretchingForces(const std::vector<BondMDS>& Bonds, std::vector<AtomMDS>& Atoms)
{
    for (const auto& Bond : Bonds)
    {
        AtomMDS& a1 = Atoms[Bond.Atom1Index];
        AtomMDS& a2 = Atoms[Bond.Atom2Index];

        const MDSRealType dx = a2.PositionX - a1.PositionX;
        const MDSRealType dy = a2.PositionY - a1.PositionY;
        const MDSRealType dz = a2.PositionZ - a1.PositionZ;
        const MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
        const MDSRealType dr = r - Bond.r0;
        const MDSRealType force_magnitude = -Bond.k_bond * dr;
        const MDSRealType ForceX = force_magnitude * dx / r;
        const MDSRealType ForceY = force_magnitude * dy / r;
        const MDSRealType ForceZ = force_magnitude * dz / r;

        a1.ForceX -= ForceX;
        a1.ForceY -= ForceY;
        a1.ForceZ -= ForceZ;
        a2.ForceX += ForceX;
        a2.ForceY += ForceY;
        a2.ForceZ += ForceZ;
    }
}

void ComputeDihedralTorsionForceBetween4Atoms_Version1(const AtomMDS& a, const AtomMDS& b, const AtomMDS& c, const AtomMDS& d, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    // Vectors for the dihedral
    const MDSRealType b1x = b.PositionX - a.PositionX, b1y = b.PositionY - a.PositionY, b1z = b.PositionZ - a.PositionZ;
    const MDSRealType b2x = c.PositionX - b.PositionX, b2y = c.PositionY - b.PositionY, b2z = c.PositionZ - b.PositionZ;
    const MDSRealType b3x = d.PositionX - c.PositionX, b3y = d.PositionY - c.PositionY, b3z = d.PositionZ - c.PositionZ;

    // Cross products
    MDSRealType n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    MDSRealType n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize cross products
    const MDSRealType n1_mag = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    const MDSRealType n2_mag = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_mag; n1y /= n1_mag; n1z /= n1_mag;
    n2x /= n2_mag; n2y /= n2_mag; n2z /= n2_mag;

    // Compute dihedral angle (phi)
    const MDSRealType cos_phi = n1x * n2x + n1y * n2y + n1z * n2z;
    const MDSRealType sin_phi = (b2x * (n1y * n2z - n1z * n2y) + b2y * (n1z * n2x - n1x * n2z) + b2z * (n1x * n2y - n1y * n2x)) / (b2x * b2x + b2y * b2y + b2z * b2z);
    const MDSRealType phi = atan2(sin_phi, cos_phi);

    // Compute force magnitude
    const MDSRealType force_magnitude = -k_dihedral * n * sin(n * phi);

    // Compute forces on atoms
    //MDSRealType dphi_dr1x = ...; // Derivative of phi with respect to atom1's position
    //MDSRealType dphi_dr1y = ...;
    //MDSRealType dphi_dr1z = ...;
    //MDSRealType dphi_dr4x = ...; // Derivative of phi with respect to atom4's position
    //MDSRealType dphi_dr4y = ...;
    //MDSRealType dphi_dr4z = ...;

    //ForceX += force_magnitude * dphi_dr1x;
    //ForceY += force_magnitude * dphi_dr1y;
    //ForceZ += force_magnitude * dphi_dr1z;
    //ForceX += force_magnitude * dphi_dr4x;
    //ForceY += force_magnitude * dphi_dr4y;
    //ForceZ += force_magnitude * dphi_dr4z;
}

MDSRealType ComputeDihedralTorsionForceBetween4Atoms_Version2(const AtomMDS& a1, const AtomMDS& a2, const AtomMDS& a3, const AtomMDS& a4)
{
    // Vectors between atoms
    const MDSRealType b1x = a2.PositionX - a1.PositionX, b1y = a2.PositionY - a1.PositionY, b1z = a2.PositionZ - a1.PositionZ;
    const MDSRealType b2x = a3.PositionX - a2.PositionX, b2y = a3.PositionY - a2.PositionY, b2z = a3.PositionZ - a2.PositionZ;
    const MDSRealType b3x = a4.PositionX - a3.PositionX, b3y = a4.PositionY - a3.PositionY, b3z = a4.PositionZ - a3.PositionZ;

    // Cross products
    MDSRealType n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    MDSRealType n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize vectors
    const MDSRealType n1_norm = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    const MDSRealType n2_norm = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_norm; n1y /= n1_norm; n1z /= n1_norm;
    n2x /= n2_norm; n2y /= n2_norm; n2z /= n2_norm;

    // Dihedral angle
    const MDSRealType dot = n1x * n2x + n1y * n2y + n1z * n2z;
    const MDSRealType cross_x = n1y * n2z - n1z * n2y;
    const MDSRealType cross_y = n1z * n2x - n1x * n2z;
    const MDSRealType cross_z = n1x * n2y - n1y * n2x;
    const MDSRealType cross_norm = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    const MDSRealType phi = atan2(cross_norm, dot);

    return phi;
}

void ComputePositionsAndVelocitiesByVelocityVerletAlgorithms(AtomMDS& AtomObject)
{
    AtomObject.VelocityX += AtomObject.ForceX / AtomObject.Mass * dt;
    AtomObject.VelocityY += AtomObject.ForceY / AtomObject.Mass * dt;
    AtomObject.VelocityZ += AtomObject.ForceZ / AtomObject.Mass * dt;
    AtomObject.PositionX += AtomObject.VelocityX * dt;
    AtomObject.PositionY += AtomObject.VelocityY * dt;
    AtomObject.PositionZ += AtomObject.VelocityZ * dt;
}

void Simulate(vector<AtomMDS>& Atoms, const vector<BondMDS>& Bonds, const vector<AngleMDS>& Angles, const vector<DihedralMDS>& Dihedrals, const int Steps)
{
    for (int Step = 0; Step < Steps; ++Step)
    {
        for (AtomMDS& AtomObject : Atoms)
            AtomObject.ForceX = AtomObject.ForceY = AtomObject.ForceZ = 0.0;

        for (UnsignedInt Atom1Index = 0; Atom1Index < Atoms.size(); ++Atom1Index)
            for (UnsignedInt Atom2Index = Atom1Index + 1; Atom2Index < Atoms.size(); ++Atom2Index)
            {
                ComputeLennardJonesForce(Atoms[Atom1Index], Atoms[Atom2Index], Atoms[Atom1Index].ForceX, Atoms[Atom1Index].ForceY, Atoms[Atom1Index].ForceZ);
                ComputeCoulombForce(Atoms[Atom1Index], Atoms[Atom2Index], Atoms[Atom1Index].ForceX, Atoms[Atom1Index].ForceY, Atoms[Atom1Index].ForceZ);
            }

        for (const BondMDS& BondObject : Bonds)
            ComputeBondedHookesLawForce(Atoms[BondObject.Atom1Index], Atoms[BondObject.Atom2Index], Atoms[BondObject.Atom1Index].ForceX, Atoms[BondObject.Atom1Index].ForceY, Atoms[BondObject.Atom1Index].ForceZ);

        for (const AngleMDS& AngleObject : Angles)
            ComputeAngleBendingForce(Atoms[AngleObject.Atom1Index], Atoms[AngleObject.Atom2Index], Atoms[AngleObject.Atom3Index], Atoms[AngleObject.Atom1Index].ForceX, Atoms[AngleObject.Atom1Index].ForceY, Atoms[AngleObject.Atom1Index].ForceZ);

        for (const DihedralMDS& DihedralObject : Dihedrals)
            ComputeDihedralTorsionForceBetween4Atoms_Version1(Atoms[DihedralObject.Atom1Index], Atoms[DihedralObject.Atom2Index], Atoms[DihedralObject.Atom3Index], Atoms[DihedralObject.Atom4Index], Atoms[DihedralObject.Atom1Index].ForceX, Atoms[DihedralObject.Atom1Index].ForceY, Atoms[DihedralObject.Atom1Index].ForceZ);

        for (AtomMDS& AtomObject : Atoms)
            ComputePositionsAndVelocitiesByVelocityVerletAlgorithms(AtomObject);
    }
}

int ComputeMolecularDynamicsSimulationForceField4()
{
    std::mt19937 RandomNumberGenerator(std::chrono::system_clock::now().time_since_epoch().count());

    std::uniform_real_distribution<MDSRealType> PositionsDistributionInNM(-1e-9, 1e-9);
    std::uniform_real_distribution<MDSRealType> VelocityDistributionInMPerS(-1e3, 1e3);
    std::uniform_real_distribution<MDSRealType> ChargesInElementaryChargeUnitsDistribution(-1.0, 1.0);

    vector<AtomMDS> Atoms(1000);
    for (AtomMDS& atom : Atoms)
    {
        atom.PositionX = PositionsDistributionInNM(RandomNumberGenerator);
        atom.PositionY = PositionsDistributionInNM(RandomNumberGenerator);
        atom.PositionZ = PositionsDistributionInNM(RandomNumberGenerator);
        atom.VelocityX = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.VelocityY = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.VelocityZ = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.Charge = ChargesInElementaryChargeUnitsDistribution(RandomNumberGenerator) * e_charge;
        atom.Mass = 1.67e-27; // Approximate mass of a proton (kg)
    }

    constexpr vector<BondMDS> Bonds;
    constexpr vector<AngleMDS> Angles;
    constexpr vector<DihedralMDS> Dihedrals;

    Simulate(Atoms, Bonds, Angles, Dihedrals, 1000);
    return 0;
}

