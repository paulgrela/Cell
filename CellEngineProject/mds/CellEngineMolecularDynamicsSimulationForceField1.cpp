
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField1.h"
#include "CellEngineMolecularDynamicsSimulationForceField1Constants.h"
#include "CellEngineMolecularDynamicsSimulationForceFieldCommonTypes.h"

using namespace std;

MDSRealType distance_squared(const AtomMDS &p1, const AtomMDS &p2)
{
    const MDSRealType dx = p1.PositionX - p2.PositionX;
    const MDSRealType dy = p1.PositionY - p2.PositionY;
    const MDSRealType dz = p1.PositionZ - p2.PositionZ;

    return dx * dx + dy * dy + dz * dz;
}

void InitializeAtoms1(std::vector<AtomMDS>& Atoms, const UnsignedInt NumberOfAtoms)
{
    std::random_device RandomDevice;
    std::mt19937 RandomGenerator(RandomDevice());
    std::uniform_real_distribution<> RandomPositionsInNMDistribution(-1.0e-9, 1.0e-9);
    std::uniform_real_distribution<> RandomVelocitiesDistribution(-1.0, 1.0);
    std::uniform_int_distribution<> RandomChargeDistribution(-1, 1);
    std::uniform_int_distribution<> RandomAtomTypeDistribution(0, 3); // Randomly choose AtomObject type (O, N, P, C)

    for (UnsignedInt AtomIndex = 0; AtomIndex < NumberOfAtoms; ++AtomIndex)
    {
        AtomMDS AtomObject{};

        AtomObject.PositionX = RandomPositionsInNMDistribution(RandomGenerator);
        AtomObject.PositionY = RandomPositionsInNMDistribution(RandomGenerator);
        AtomObject.PositionZ = RandomPositionsInNMDistribution(RandomGenerator);
        AtomObject.VelocityX = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB_BoltzmanConstantJPerK * TemperatureInKelvins / mass_C_kg);
        AtomObject.VelocityY = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB_BoltzmanConstantJPerK * TemperatureInKelvins / mass_C_kg);
        AtomObject.VelocityZ = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB_BoltzmanConstantJPerK * TemperatureInKelvins / mass_C_kg);
        AtomObject.Charge = RandomChargeDistribution(RandomGenerator) * 1.602176634e-19;
        AtomObject.Mass = (RandomAtomTypeDistribution(RandomGenerator) == 0) ? mass_O_kg : ((RandomAtomTypeDistribution(RandomGenerator) == 1) ? mass_N_kg : ((RandomAtomTypeDistribution(RandomGenerator) == 2) ? mass_P_kg : mass_C_kg));
        AtomObject.ForceX = AtomObject.ForceY = AtomObject.ForceZ = 0.0;

        Atoms.emplace_back(AtomObject);
    }
}

void InitializeAtoms2(std::vector<AtomMDS>& Atoms, const UnsignedInt NumberOfAtoms)
{
    Atoms.resize(NumberOfAtoms);

    std::mt19937 ReandomGenerator(std::random_device{}());
    std::uniform_real_distribution<> Distribution(-5.0, 5.0);

    for (auto &AtomObject : Atoms)
    {
        AtomObject.PositionX = Distribution(ReandomGenerator);
        AtomObject.PositionY = Distribution(ReandomGenerator);
        AtomObject.PositionZ = Distribution(ReandomGenerator);
        AtomObject.VelocityX = Distribution(ReandomGenerator);
        AtomObject.VelocityY = Distribution(ReandomGenerator);
        AtomObject.VelocityZ = Distribution(ReandomGenerator);
        AtomObject.Mass = 1.0;
        AtomObject.ForceX = AtomObject.ForceY = AtomObject.ForceZ = 0.0;
    }
}

void InitializeBondsBetweenAtoms(const vector<AtomMDS>& Atoms, vector<BondMDS>& Bonds)
{
    std::random_device RandomDevice;
    std::mt19937 RandomGenerator(RandomDevice());
    std::uniform_int_distribution<UnsignedInt> AtomIndexDistribution(0, Atoms.size() - 1);

    const UnsignedInt NumberOfBonds = Atoms.size() / 4; // Roughly 25% of atoms are bonded

    for (UnsignedInt BindIndex = 0; BindIndex < NumberOfBonds; ++BindIndex)
    {
        UnsignedInt Atom1Index = AtomIndexDistribution(RandomGenerator);
        UnsignedInt Atom2Index = AtomIndexDistribution(RandomGenerator);

        if (Atom1Index != Atom2Index)
            Bonds.push_back({ Atom1Index, Atom2Index });
    }
}

void ComputeLennardJonesPotentialAndForces(AtomMDS &Atom1, AtomMDS &Atom2)
{
    const MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
    const MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
    const MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;
    const MDSRealType r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < r_cutoff_LennardJonesCutoffDistance * r_cutoff_LennardJonesCutoffDistance)
    {
        MDSRealType r6 = (sigma_LennardJonesDistanceParameter * sigma_LennardJonesDistanceParameter) / r2;
        r6 = r6 * r6 * r6;  // (sigma/r)^6
        const MDSRealType r12 = r6 * r6; // (sigma/r)^12

        const MDSRealType ForceScalar = 48 * epsilon_LennardJonesPotentialDepth * (r12 - 0.5 * r6) / r2;

        Atom1.ForceX += ForceScalar * dx;
        Atom1.ForceY += ForceScalar * dy;
        Atom1.ForceZ += ForceScalar * dz;

        Atom2.ForceX -= ForceScalar * dx;
        Atom2.ForceY -= ForceScalar * dy;
        Atom2.ForceZ -= ForceScalar * dz;
    }
}

void ComputeLennardJonesPotentialAndForces(AtomMDS &Atom1, AtomMDS &Atom2, MDSRealType r2)
{
    if (r2 < r_cutoff_LennardJonesCutoffDistance * r_cutoff_LennardJonesCutoffDistance)
    {
        MDSRealType r6 = (sigma_LennardJonesDistanceParameter * sigma_LennardJonesDistanceParameter) / r2;
        r6 = r6 * r6 * r6; // (sigma/r)^6
        const MDSRealType r12 = r6 * r6; // (sigma/r)^12

        const MDSRealType ForceScalar = 48 * epsilon_LennardJonesPotentialDepth * (r12 - 0.5 * r6) / r2;
        const MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
        const MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
        const MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;

        Atom1.ForceX += ForceScalar * dx;
        Atom1.ForceY += ForceScalar * dy;
        Atom1.ForceZ += ForceScalar * dz;
        Atom2.ForceX -= ForceScalar * dx;
        Atom2.ForceY -= ForceScalar * dy;
        Atom2.ForceZ -= ForceScalar * dz;
    }
}

void ComputeCoulombElectrostaticForces(AtomMDS &Atom1, AtomMDS &Atom2, MDSRealType r2)
{
    if (r2 < r_cutoff_LennardJonesCutoffDistance * r_cutoff_LennardJonesCutoffDistance)
    {
        MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
        MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
        MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;
        //MDSRealType r2 = dx * dx + dy * dy + dz * dz;
        MDSRealType r = std::sqrt(r2);

        MDSRealType ForceScalar = (k_e_CoulombConstant * Atom1.Charge * Atom2.Charge) / (r2 * r);

        Atom1.ForceX += ForceScalar * dx;
        Atom1.ForceY += ForceScalar * dy;
        Atom1.ForceZ += ForceScalar * dz;

        Atom2.ForceX -= ForceScalar * dx;
        Atom2.ForceY -= ForceScalar * dy;
        Atom2.ForceZ -= ForceScalar * dz;
    }
}

void ComputeBondStretchingHookesLaw(AtomMDS &Atom1, AtomMDS &Atom2)
{
    MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
    MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
    MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;
    MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);

    MDSRealType ForceScalar = -k_bond_HookesLawBondStretchingForceConstant * (r - r0_bond_HookesLawBondStretchingEquilibriumBondLength) / r;
    Atom1.ForceX += ForceScalar * dx;
    Atom1.ForceY += ForceScalar * dy;
    Atom1.ForceZ += ForceScalar * dz;

    Atom2.ForceX -= ForceScalar * dx;
    Atom2.ForceY -= ForceScalar * dy;
    Atom2.ForceZ -= ForceScalar * dz;
}

void ComputeAngleBendingForcesFor3Atoms(AtomMDS &Atom1, AtomMDS &Atom2, AtomMDS &Atom3, const MDSRealType k_theta, const MDSRealType theta0)
{
    // Vector p1 -> p2
    MDSRealType dx1 = Atom1.PositionX - Atom2.PositionX;
    MDSRealType dy1 = Atom1.PositionY - Atom2.PositionY;
    MDSRealType dz1 = Atom1.PositionZ - Atom2.PositionZ;
    MDSRealType r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

    // Vector p3 -> p2
    MDSRealType dx2 = Atom3.PositionX - Atom2.PositionX;
    MDSRealType dy2 = Atom3.PositionY - Atom2.PositionY;
    MDSRealType dz2 = Atom3.PositionZ - Atom2.PositionZ;
    MDSRealType r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    // Compute the cosine of the angle
    MDSRealType cos_theta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (r1 * r2);
    MDSRealType theta = std::acos(cos_theta);  // Angle in radians

    // Compute the force based on angle deviation from equilibrium
    MDSRealType angle_force = -k_theta * (theta - theta0);

    // Compute forces to apply to each ParticleMDS (force projection along each axis)
    MDSRealType fx1 = angle_force * (dx1 / r1);
    MDSRealType fy1 = angle_force * (dy1 / r1);
    MDSRealType fz1 = angle_force * (dz1 / r1);

    MDSRealType fx3 = angle_force * (dx2 / r2);
    MDSRealType fy3 = angle_force * (dy2 / r2);
    MDSRealType fz3 = angle_force * (dz2 / r2);

    // Apply forces to the atoms
    Atom1.ForceX += fx1;
    Atom1.ForceY += fy1;
    Atom1.ForceZ += fz1;

    Atom3.ForceX += fx3;
    Atom3.ForceY += fy3;
    Atom3.ForceZ += fz3;

    // Apply equal and opposite force to the central atom (p2)
    Atom2.ForceX -= (fx1 + fx3);
    Atom2.ForceY -= (fy1 + fy3);
    Atom2.ForceZ -= (fz1 + fz3);
}

void ComputeDihedralTorsionForcesFor4Atoms(AtomMDS &Atom1, AtomMDS &Atom2, AtomMDS &Atom3, AtomMDS &Atom4, const MDSRealType Vn, UnsignedInt N, MDSRealType Gamma)
{
    // Compute bond vectors (p2 -> p1, p2 -> p3, p3 -> p4)
    MDSRealType dx21 = Atom1.PositionX - Atom2.PositionX;
    MDSRealType dy21 = Atom1.PositionY - Atom2.PositionY;
    MDSRealType dz21 = Atom1.PositionZ - Atom2.PositionZ;

    MDSRealType dx23 = Atom3.PositionX - Atom2.PositionX;
    MDSRealType dy23 = Atom3.PositionY - Atom2.PositionY;
    MDSRealType dz23 = Atom3.PositionZ - Atom2.PositionZ;

    MDSRealType dx34 = Atom4.PositionX - Atom3.PositionX;
    MDSRealType dy34 = Atom4.PositionY - Atom3.PositionY;
    MDSRealType dz34 = Atom4.PositionZ - Atom3.PositionZ;

    // Calculate the normal vectors to the planes formed by p1-p2-p3 and p2-p3-p4
    MDSRealType nx1 = dy21 * dz23 - dz21 * dy23;
    MDSRealType ny1 = dz21 * dx23 - dx21 * dz23;
    MDSRealType nz1 = dx21 * dy23 - dy21 * dx23;

    MDSRealType nx2 = dy23 * dz34 - dz23 * dy34;
    MDSRealType ny2 = dz23 * dx34 - dx23 * dz34;
    MDSRealType nz2 = dx23 * dy34 - dy23 * dx34;

    // Calculate the magnitude of these normal vectors
    MDSRealType n1_mag = sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
    MDSRealType n2_mag = sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);

    // Compute the vector between p2 and p3
    MDSRealType dx = Atom3.PositionX - Atom2.PositionX;
    MDSRealType dy = Atom3.PositionY - Atom2.PositionY;
    MDSRealType dz = Atom3.PositionZ - Atom2.PositionZ;
    MDSRealType r23 = sqrt(dx * dx + dy * dy + dz * dz);

    // Calculate the dihedral angle (phi)
    MDSRealType cos_phi = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (n1_mag * n2_mag);
    MDSRealType sin_phi = r23 * (nx1 * dx34 + ny1 * dy34 + nz1 * dz34) / (n1_mag * n2_mag);
    MDSRealType phi = atan2(sin_phi, cos_phi);  // Dihedral angle in radians

    // Torsion potential force: V(phi) = 0.5 * Vn * (1 + cos(n * phi - gamma))
    MDSRealType TorsionEnergy = 0.5 * Vn * (1 + std::cos(N * phi - Gamma));
    MDSRealType TorsionForce = -0.5 * Vn * N * sin(N * phi - Gamma);

    // Apply forces on p1, p2, p3, and p4
    // Normalize the normals
    nx1 /= n1_mag;
    ny1 /= n1_mag;
    nz1 /= n1_mag;

    nx2 /= n2_mag;
    ny2 /= n2_mag;
    nz2 /= n2_mag;

    // Calculate forces acting on each atom
    // Force on atom p1 (affects normal vector 1)
    MDSRealType fx1 = TorsionForce * nx1;
    MDSRealType fy1 = TorsionForce * ny1;
    MDSRealType fz1 = TorsionForce * nz1;

    // Force on atom p4 (affects normal vector 2)
    MDSRealType fx4 = -TorsionForce * nx2;
    MDSRealType fy4 = -TorsionForce * ny2;
    MDSRealType fz4 = -TorsionForce * nz2;

    // Force on atom p2 and p3 (distribute between them equally)
    MDSRealType fx23 = TorsionForce * (nx1 - nx2);
    MDSRealType fy23 = TorsionForce * (ny1 - ny2);
    MDSRealType fz23 = TorsionForce * (nz1 - nz2);

    // Apply forces to the atoms
    Atom1.ForceX += fx1;
    Atom1.ForceY += fy1;
    Atom1.ForceZ += fz1;

    Atom2.ForceX += 0.5 * fx23;
    Atom2.ForceY += 0.5 * fy23;
    Atom2.ForceZ += 0.5 * fz23;

    Atom3.ForceX += 0.5 * fx23;
    Atom3.ForceY += 0.5 * fy23;
    Atom3.ForceZ += 0.5 * fz23;

    Atom4.ForceX += fx4;
    Atom4.ForceY += fy4;
    Atom4.ForceZ += fz4;
}

void UpdatePositionsAndVelocities(std::vector<AtomMDS>& Atoms, const std::vector<std::vector<MDSRealType>>& Forces)
{
    for (UnsignedInt AtomIndex = 0; AtomIndex < Atoms.size(); ++AtomIndex)
    {
        Atoms[AtomIndex].PositionX += Atoms[AtomIndex].VelocityX * dt_TimeStepFirstInSeconds + 0.5 * Forces[AtomIndex][0] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds * dt_TimeStepFirstInSeconds;
        Atoms[AtomIndex].PositionY += Atoms[AtomIndex].VelocityY * dt_TimeStepFirstInSeconds + 0.5 * Forces[AtomIndex][1] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds * dt_TimeStepFirstInSeconds;
        Atoms[AtomIndex].PositionZ += Atoms[AtomIndex].VelocityZ * dt_TimeStepFirstInSeconds + 0.5 * Forces[AtomIndex][2] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds * dt_TimeStepFirstInSeconds;

        Atoms[AtomIndex].VelocityX += 0.5 * Forces[AtomIndex][0] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds;
        Atoms[AtomIndex].VelocityY += 0.5 * Forces[AtomIndex][1] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds;
        Atoms[AtomIndex].VelocityZ += 0.5 * Forces[AtomIndex][2] / Atoms[AtomIndex].Mass * dt_TimeStepFirstInSeconds;
    }
}

void PerformVelocityVerletIntegration(std::vector<AtomMDS> &Atoms, const std::vector<BondMDS>& Bonds)
{
    for (auto &AtomObject : Atoms)
    {
        AtomObject.VelocityX += 0.5 * AtomObject.ForceX / AtomObject.Mass * dt_TimeStepFirstInSeconds;
        AtomObject.VelocityY += 0.5 * AtomObject.ForceY / AtomObject.Mass * dt_TimeStepFirstInSeconds;
        AtomObject.VelocityZ += 0.5 * AtomObject.ForceZ / AtomObject.Mass * dt_TimeStepFirstInSeconds;

        AtomObject.PositionX += AtomObject.VelocityX * dt_TimeStepFirstInSeconds;
        AtomObject.PositionY += AtomObject.VelocityY * dt_TimeStepFirstInSeconds;
        AtomObject.PositionZ += AtomObject.VelocityZ * dt_TimeStepFirstInSeconds;
    }

    for (UnsignedInt AtomIndex1 = 0; AtomIndex1 < Atoms.size(); ++AtomIndex1)
    {
        for (UnsignedInt AtomIndex2 = AtomIndex1 + 1; AtomIndex2 < Atoms.size(); ++AtomIndex2)
        {
            MDSRealType r2 = distance_squared(Atoms[AtomIndex1], Atoms[AtomIndex2]);
            ComputeLennardJonesPotentialAndForces(Atoms[AtomIndex1], Atoms[AtomIndex2], r2);
            //ComputeLennardJonesPotentialAndForces((particles[i], particles[j]);
            ComputeCoulombElectrostaticForces(Atoms[AtomIndex1], Atoms[AtomIndex2], r2);
        }
    }

    // bont streching produce here error of counting
    // for (const auto &Bond : Bonds)
    //     ComputeBondStretchingHookesLaw(Atoms[Bond.Atom1Index], Atoms[Bond.Atom2Index]);

    for (auto &ParticleObject : Atoms)
    {
        ParticleObject.VelocityX += 0.5 * ParticleObject.ForceX / ParticleObject.Mass * dt_TimeStepFirstInSeconds;
        ParticleObject.VelocityY += 0.5 * ParticleObject.ForceY / ParticleObject.Mass * dt_TimeStepFirstInSeconds;
        ParticleObject.VelocityZ += 0.5 * ParticleObject.ForceZ / ParticleObject.Mass * dt_TimeStepFirstInSeconds;
    }
}

void ApplyBerendsenThermostatToControlTemperature(std::vector<AtomMDS> &Atoms, MDSRealType CurrentTemperature)
{
    MDSRealType Lambda = sqrt(1 + dt_TimeStepFirstInSeconds / (TargetTempratureForTheThermostat - CurrentTemperature) * (TargetTempratureForTheThermostat / CurrentTemperature - 1));

    for (auto &AtomObject : Atoms)
    {
        AtomObject.VelocityX *= Lambda;
        AtomObject.VelocityY *= Lambda;
        AtomObject.VelocityZ *= Lambda;
    }
}

MDSRealType ComputeCurrentKineticTemperature(const std::vector<AtomMDS> &Atoms)
{
    MDSRealType KineticEnergy = 0.0;

    for (const auto &AtomObject : Atoms)
        KineticEnergy += 0.5 * AtomObject.Mass * (AtomObject.VelocityX * AtomObject.VelocityX + AtomObject.VelocityY * AtomObject.VelocityY + AtomObject.VelocityZ * AtomObject.VelocityZ);

    return (2.0 * KineticEnergy) / (3.0 * static_cast<MDSRealType>(Atoms.size()));
}

UnsignedInt ComputeMolecularDynamicsSimulationForceField1()
{
    constexpr UnsignedInt NumberOfAtoms = 100;

    vector<AtomMDS> Particles;
    vector<BondMDS> Bonds;

    InitializeAtoms1(Particles, NumberOfAtoms);
    InitializeBondsBetweenAtoms(Particles, Bonds);

    for (UnsignedInt SimulationStep = 0; SimulationStep < 10000; ++SimulationStep)
    {
        PerformVelocityVerletIntegration(Particles, Bonds);

        const MDSRealType CurrentTemperature = ComputeCurrentKineticTemperature(Particles);

        ApplyBerendsenThermostatToControlTemperature(Particles, CurrentTemperature);

        if (SimulationStep % 100 == 0)
        {
           cout << "Step " << SimulationStep << endl;

           for (const auto& AtomObject : Particles)
               cout << "Position: (" << AtomObject.PositionX << ", " << AtomObject.PositionY << ", " << AtomObject.PositionZ << ")" << " Velocity: (" << AtomObject.VelocityX << ", " << AtomObject.VelocityY << ", " << AtomObject.VelocityZ << ")" << endl;
        }
    }

    cout << "Simulation complete." << endl;

    return 0;
}
