
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField4.h"
#include "CellEngineMolecularDynamicsSimulationForceField4Constants.h"
#include "CellEngineMolecularDynamicsSimulationForceFieldCommonTypes.h"

using namespace std;

void ComputeLennardJonesForce(const AtomMDS& Atom1, const AtomMDS& Atom2, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = Atom1.PositionX - Atom2.PositionX, dy = Atom1.PositionY - Atom2.PositionY, dz = Atom1.PositionZ - Atom2.PositionZ;
    const MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    const MDSRealType r = sqrt(r2);
    const MDSRealType r6 = pow(sigma_LennardJonesPotentialDistanceParameterAngstromPerM / r, 6);
    const MDSRealType r12 = r6 * r6;
    const MDSRealType force = 24 * epsilon_LennardJonesPotentialDepthJouls * (2 * r12 - r6) / r;

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeCoulombForce(const AtomMDS& Atom1, const AtomMDS& Atom2, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = Atom1.PositionX - Atom2.PositionX, dy = Atom1.PositionY - Atom2.PositionY, dz = Atom1.PositionZ - Atom2.PositionZ;
    const MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    const MDSRealType r = sqrt(r2);
    const MDSRealType force = (Atom1.Charge * Atom2.Charge) / (4 * pi * epsilon0_VacuumPermittivityFPerM * r * r);

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeBondedHookesLawForce(const AtomMDS& Atom1, const AtomMDS& Atom2, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType dx = Atom1.PositionX - Atom2.PositionX, dy = Atom1.PositionY - Atom2.PositionY, dz = Atom1.PositionZ - Atom2.PositionZ;
    const MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
    const MDSRealType force = -k_bond_HookesLawBondedInteractionSpringConstantNPerM * (r - r0_EquilibriumBondLengthInM);

    ForceX += force * dx / r;
    ForceY += force * dy / r;
    ForceZ += force * dz / r;
}

void ComputeAngleBendingForce(const AtomMDS& Atom1, const AtomMDS& Atom2, const AtomMDS& Atom3, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    const MDSRealType abx = Atom1.PositionX - Atom2.PositionX, aby = Atom1.PositionY - Atom2.PositionY, abz = Atom1.PositionZ - Atom2.PositionZ;
    const MDSRealType cbx = Atom3.PositionX - Atom2.PositionX, cby = Atom3.PositionY - Atom2.PositionY, cbz = Atom3.PositionZ - Atom2.PositionZ;
    const MDSRealType dot = abx * cbx + aby * cby + abz * cbz;
    const MDSRealType mag_ab = sqrt(abx * abx + aby * aby + abz * abz);
    const MDSRealType mag_cb = sqrt(cbx * cbx + cby * cby + cbz * cbz);
    const MDSRealType theta = acos(dot / (mag_ab * mag_cb));
    const MDSRealType force = -k_angle_AngleBendingSpringConstant * (theta - theta0_EquilibriumAngle90Degree);

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
        AtomMDS& Atom1 = Atoms[Bond.Atom1Index];
        AtomMDS& Atom2 = Atoms[Bond.Atom2Index];

        const MDSRealType dx = Atom2.PositionX - Atom1.PositionX;
        const MDSRealType dy = Atom2.PositionY - Atom1.PositionY;
        const MDSRealType dz = Atom2.PositionZ - Atom1.PositionZ;
        const MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
        const MDSRealType dr = r - Bond.r0;
        const MDSRealType force_magnitude = -Bond.k_bond * dr;
        const MDSRealType ForceX = force_magnitude * dx / r;
        const MDSRealType ForceY = force_magnitude * dy / r;
        const MDSRealType ForceZ = force_magnitude * dz / r;

        Atom1.ForceX -= ForceX;
        Atom1.ForceY -= ForceY;
        Atom1.ForceZ -= ForceZ;
        Atom2.ForceX += ForceX;
        Atom2.ForceY += ForceY;
        Atom2.ForceZ += ForceZ;
    }
}

void ComputeDihedralTorsionForceBetween4Atoms_Version1(const AtomMDS& Atom1, const AtomMDS& Atom2, const AtomMDS& Atom3, const AtomMDS& Atom4, MDSRealType& ForceX, MDSRealType& ForceY, MDSRealType& ForceZ)
{
    // Vectors for the dihedral
    const MDSRealType b1x = Atom2.PositionX - Atom1.PositionX, b1y = Atom2.PositionY - Atom1.PositionY, b1z = Atom2.PositionZ - Atom1.PositionZ;
    const MDSRealType b2x = Atom3.PositionX - Atom2.PositionX, b2y = Atom3.PositionY - Atom2.PositionY, b2z = Atom3.PositionZ - Atom2.PositionZ;
    const MDSRealType b3x = Atom4.PositionX - Atom3.PositionX, b3y = Atom4.PositionY - Atom3.PositionY, b3z = Atom4.PositionZ - Atom3.PositionZ;

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
    const MDSRealType force_magnitude = -k_dihedral_DiherdralTorsionSpringConstantInJoules * PeriodicityOfDihedralTorsionPotential * sin(PeriodicityOfDihedralTorsionPotential * phi);

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

MDSRealType ComputeDihedralTorsionForceBetween4Atoms_Version2(const AtomMDS& Atom1, const AtomMDS& Atom2, const AtomMDS& Atom3, const AtomMDS& Atom4)
{
    // Vectors between atoms
    const MDSRealType b1x = Atom2.PositionX - Atom1.PositionX, b1y = Atom2.PositionY - Atom1.PositionY, b1z = Atom2.PositionZ - Atom1.PositionZ;
    const MDSRealType b2x = Atom3.PositionX - Atom2.PositionX, b2y = Atom3.PositionY - Atom2.PositionY, b2z = Atom3.PositionZ - Atom2.PositionZ;
    const MDSRealType b3x = Atom4.PositionX - Atom3.PositionX, b3y = Atom4.PositionY - Atom3.PositionY, b3z = Atom4.PositionZ - Atom3.PositionZ;

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
    AtomObject.VelocityX += AtomObject.ForceX / AtomObject.Mass * dt_TimeStepForSimulation1Fs;
    AtomObject.VelocityY += AtomObject.ForceY / AtomObject.Mass * dt_TimeStepForSimulation1Fs;
    AtomObject.VelocityZ += AtomObject.ForceZ / AtomObject.Mass * dt_TimeStepForSimulation1Fs;
    AtomObject.PositionX += AtomObject.VelocityX * dt_TimeStepForSimulation1Fs;
    AtomObject.PositionY += AtomObject.VelocityY * dt_TimeStepForSimulation1Fs;
    AtomObject.PositionZ += AtomObject.VelocityZ * dt_TimeStepForSimulation1Fs;
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
    for (AtomMDS& AtomObject : Atoms)
    {
        AtomObject.PositionX = PositionsDistributionInNM(RandomNumberGenerator);
        AtomObject.PositionY = PositionsDistributionInNM(RandomNumberGenerator);
        AtomObject.PositionZ = PositionsDistributionInNM(RandomNumberGenerator);
        AtomObject.VelocityX = VelocityDistributionInMPerS(RandomNumberGenerator);
        AtomObject.VelocityY = VelocityDistributionInMPerS(RandomNumberGenerator);
        AtomObject.VelocityZ = VelocityDistributionInMPerS(RandomNumberGenerator);
        AtomObject.Charge = ChargesInElementaryChargeUnitsDistribution(RandomNumberGenerator) * ElementaryCharge;
        AtomObject.Mass = 1.67e-27; // Approximate mass of a proton (kg)
    }

    constexpr vector<BondMDS> Bonds;
    constexpr vector<AngleMDS> Angles;
    constexpr vector<DihedralMDS> Dihedrals;

    Simulate(Atoms, Bonds, Angles, Dihedrals, 1000);
    return 0;
}

