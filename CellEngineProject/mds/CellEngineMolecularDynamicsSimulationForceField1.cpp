
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField1.h"

#include "CellEngineMolecularDynamicsSimulationForceField1Constants.h"

using namespace std;

struct AtomMDS
{
    MDSRealType PositionX, PositionY, PositionZ;
    MDSRealType VelocityX, VelocityY, VelocityZ;
    MDSRealType ForceX, ForceY, ForceZ;
    MDSRealType Charge;
    MDSRealType Mass;
};

struct Bond
{
    UnsignedInt Atom1Index;
    UnsignedInt Atom2Index;
};

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
    std::uniform_int_distribution<> RandomAtomTypeDistribution(0, 3); // Randomly choose atom type (O, N, P, C)

    for (UnsignedInt AtomIndex = 0; AtomIndex < NumberOfAtoms; ++AtomIndex)
    {
        AtomMDS Atom{};

        Atom.PositionX = RandomPositionsInNMDistribution(RandomGenerator);
        Atom.PositionY = RandomPositionsInNMDistribution(RandomGenerator);
        Atom.PositionZ = RandomPositionsInNMDistribution(RandomGenerator);
        Atom.VelocityX = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB * T / mass_C);
        Atom.VelocityY = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB * T / mass_C);
        Atom.VelocityZ = RandomVelocitiesDistribution(RandomGenerator) * sqrt(kB * T / mass_C);
        Atom.Charge = RandomChargeDistribution(RandomGenerator) * 1.602176634e-19;
        Atom.Mass = (RandomAtomTypeDistribution(RandomGenerator) == 0) ? mass_O : ((RandomAtomTypeDistribution(RandomGenerator) == 1) ? mass_N : ((RandomAtomTypeDistribution(RandomGenerator) == 2) ? mass_P : mass_C));
        Atom.ForceX = Atom.ForceY = Atom.ForceZ = 0.0;

        Atoms.emplace_back(Atom);
    }
}

void InitializeAtoms2(std::vector<AtomMDS>& Atoms, const UnsignedInt NumberOfAtoms)
{
    Atoms.resize(NumberOfAtoms);

    std::mt19937 ReandomGenerator(std::random_device{}());
    std::uniform_real_distribution<> Distribution(-5.0, 5.0);

    for (auto &Atom : Atoms)
    {
        Atom.PositionX = Distribution(ReandomGenerator);
        Atom.PositionY = Distribution(ReandomGenerator);
        Atom.PositionZ = Distribution(ReandomGenerator);
        Atom.VelocityX = Distribution(ReandomGenerator);
        Atom.VelocityY = Distribution(ReandomGenerator);
        Atom.VelocityZ = Distribution(ReandomGenerator);
        Atom.Mass = 1.0;
        Atom.ForceX = Atom.ForceY = Atom.ForceZ = 0.0;
    }
}

void InitializeBondsBetweenAtoms(const vector<AtomMDS>& Atoms, vector<Bond>& Bonds)
{
    std::random_device RandomDevice;
    std::mt19937 RandomGenerator(RandomDevice());
    std::uniform_int_distribution<> AtomIndexDistribution(0, Atoms.size() - 1);

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
    MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
    MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
    MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;
    MDSRealType r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6;  // (sigma/r)^6
        const MDSRealType r12 = r6 * r6; // (sigma/r)^12

        const MDSRealType ForceScalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;

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
    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6; // (sigma/r)^6
        const MDSRealType r12 = r6 * r6; // (sigma/r)^12

        const MDSRealType ForceScalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;
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
    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
        MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
        MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;
        //MDSRealType r2 = dx * dx + dy * dy + dz * dz;
        MDSRealType r = std::sqrt(r2);

        MDSRealType ForceScalar = (k_e * Atom1.Charge * Atom2.Charge) / (r2 * r);

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

    MDSRealType ForceScalar = -k_bond * (r - r0_bond) / r;
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
        Atoms[AtomIndex].PositionX += Atoms[AtomIndex].VelocityX * dt + 0.5 * Forces[AtomIndex][0] / Atoms[AtomIndex].Mass * dt * dt;
        Atoms[AtomIndex].PositionY += Atoms[AtomIndex].VelocityY * dt + 0.5 * Forces[AtomIndex][1] / Atoms[AtomIndex].Mass * dt * dt;
        Atoms[AtomIndex].PositionZ += Atoms[AtomIndex].VelocityZ * dt + 0.5 * Forces[AtomIndex][2] / Atoms[AtomIndex].Mass * dt * dt;

        Atoms[AtomIndex].VelocityX += 0.5 * Forces[AtomIndex][0] / Atoms[AtomIndex].Mass * dt;
        Atoms[AtomIndex].VelocityY += 0.5 * Forces[AtomIndex][1] / Atoms[AtomIndex].Mass * dt;
        Atoms[AtomIndex].VelocityZ += 0.5 * Forces[AtomIndex][2] / Atoms[AtomIndex].Mass * dt;
    }
}

void PerformVelocityVerletIntegration(std::vector<AtomMDS> &Atoms, const std::vector<Bond>& Bonds)
{
    // bont streching produce here error of counting
    // for (const auto &Bond : Bonds)
    //     ComputeBondStretchingHookesLaw(Particles[Bond.Atom1Index], Particles[Bond.Atom2Index]);

    for (auto &Atom : Atoms)
    {
        Atom.VelocityX += 0.5 * Atom.ForceX / Atom.Mass * dt;
        Atom.VelocityY += 0.5 * Atom.ForceY / Atom.Mass * dt;
        Atom.VelocityZ += 0.5 * Atom.ForceZ / Atom.Mass * dt;

        Atom.PositionX += Atom.VelocityX * dt;
        Atom.PositionY += Atom.VelocityY * dt;
        Atom.PositionZ += Atom.VelocityZ * dt;
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

    // for (const auto &Bond : Bonds)
    //     ComputeBondStretchingHookesLaw(Particles[Bond.Atom1Index], Particles[Bond.Atom2Index]);

    for (auto &ParticleObject : Atoms)
    {
        ParticleObject.VelocityX += 0.5 * ParticleObject.ForceX / ParticleObject.Mass * dt;
        ParticleObject.VelocityY += 0.5 * ParticleObject.ForceY / ParticleObject.Mass * dt;
        ParticleObject.VelocityZ += 0.5 * ParticleObject.ForceZ / ParticleObject.Mass * dt;
    }
}

void ApplyBerendsenThermostatToControlTemperature(std::vector<AtomMDS> &Atoms, MDSRealType CurrentTemperature)
{
    MDSRealType Lambda = sqrt(1 + dt / (target_temp - CurrentTemperature) * (target_temp / CurrentTemperature - 1));

    for (auto &Atom : Atoms)
    {
        Atom.VelocityX *= Lambda;
        Atom.VelocityY *= Lambda;
        Atom.VelocityZ *= Lambda;
    }
}

MDSRealType ComputeCurrentKineticTemperature(const std::vector<AtomMDS> &Atoms)
{
    MDSRealType KineticEnergy = 0.0;

    for (const auto &p : Atoms)
        KineticEnergy += 0.5 * p.Mass * (p.VelocityX * p.VelocityX + p.VelocityY * p.VelocityY + p.VelocityZ * p.VelocityZ);

    return (2.0 * KineticEnergy) / (3.0 * Atoms.size());
}

UnsignedInt ComputeMolecularDynamicsSimulationForceField1()
{
    constexpr UnsignedInt NumberOfAtoms = 100;

    vector<AtomMDS> Particles;
    vector<Bond> Bonds;

    InitializeAtoms1(Particles, NumberOfAtoms);
    InitializeBondsBetweenAtoms(Particles, Bonds);

    for (UnsignedInt SimulationStep = 0; SimulationStep < 10000; ++SimulationStep)
    {
        PerformVelocityVerletIntegration(Particles, Bonds);

        MDSRealType CurrentTemperature = ComputeCurrentKineticTemperature(Particles);

        ApplyBerendsenThermostatToControlTemperature(Particles, CurrentTemperature);

        if (SimulationStep % 100 == 0)
        {
           cout << "Step " << SimulationStep << endl;

           for (const auto& atom : Particles)
               cout << "Position: (" << atom.PositionX << ", " << atom.PositionY << ", " << atom.PositionZ << ")" << " Velocity: (" << atom.VelocityX << ", " << atom.VelocityY << ", " << atom.VelocityZ << ")" << endl;
        }
    }

    cout << "Simulation complete." << endl;

    return 0;
}
