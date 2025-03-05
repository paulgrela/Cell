
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField1.h"

using namespace std;

// Constants for the physical simulation
const MDSRealType kB = 1.380649e-23; // Boltzmann constant (J/K)
const MDSRealType T = 303.15; // Temperature in Kelvin (30°C)
const MDSRealType dt = 1e-15; // Time step in seconds (1 fs)
const MDSRealType e0 = 8.854187817e-12; // Permittivity of free space
const MDSRealType charge_constant = 8.9875517923e9; // Coulomb constant (1/(4*pi*e0))
const MDSRealType mass_O = 2.656e-26; // Mass of an oxygen atom (kg)
const MDSRealType mass_N = 2.326e-26; // Mass of a nitrogen atom (kg)
const MDSRealType mass_P = 5.14e-26; // Mass of a phosphorus atom (kg)
const MDSRealType mass_C = 1.994e-26; // Mass of a carbon atom (kg)
const MDSRealType bond_strength_k = 100.0; // Bond strength constant (N/m)
const MDSRealType bond_equilibrium_distance = 1.0e-10; // Equilibrium bond distance (1 Å)

const MDSRealType r_cutoff = 2.5e-10; // Lennard-Jones cutoff distance
const MDSRealType epsilon = 1.0; // Lennard-Jones potential depth
const MDSRealType sigma = 1.0e-10; // Lennard-Jones distance parameter


// Constants for Lennard-Jones potential
//const MDSRealType epsilon = 1.0; // Depth of the potential well
//const MDSRealType sigma = 1.0;   // Distance at which potential is zero
//const MDSRealType r_cutoff = 2.5 * sigma; // Cutoff distance for Lennard-Jones interactions

// Constants for Hooke's Law (bond stretching)
const MDSRealType k_bond = 100.0; // Bond force constant
const MDSRealType r0_bond = 1.0;  // Equilibrium bond length

// Constants for angle bending
const MDSRealType k_theta = 50.0;  // Force constant for angle bending
const MDSRealType theta0 = M_PI / 2; // Equilibrium angle (90 degrees)

// Constants for dihedral torsion
const MDSRealType V_n = 1.0;
const MDSRealType n_period = 3.0;
//const MDSRealType gamma = 0.0;

// Coulomb constant
const MDSRealType k_e = 8.99e9;  // Coulomb's constant

// Time step and temperature constants
// const MDSRealType dt = 0.001;  // Time step for the simulation
const MDSRealType target_temp = 300.0;  // Target temperature for the thermostat

// Particle structure
struct Particle
{
    MDSRealType x, y, z;      // Position
    MDSRealType vx, vy, vz;   // Velocity
    MDSRealType fx, fy, fz;   // Force
    MDSRealType charge;       // Charge (for Coulomb forces)
    MDSRealType mass;         // Mass
};

struct Bond
{
    int atom1_index;
    int atom2_index;
};

MDSRealType distance_squared(const Particle &p1, const Particle &p2) {
    MDSRealType dx = p1.x - p2.x;
    MDSRealType dy = p1.y - p2.y;
    MDSRealType dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}

void initialize_atoms(std::vector<Particle>& atoms)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(-1.0e-9, 1.0e-9); // Random positions in nm
    std::uniform_real_distribution<> vel_dist(-1.0, 1.0); // Random velocities
    std::uniform_int_distribution<> charge_dist(-1, 1); // Random charges (-1, 0, 1)
    std::uniform_int_distribution<> atom_type_dist(0, 3); // Randomly choose atom type (O, N, P, C)

    for (int i = 0; i < 100; ++i)
    {
        Particle atom;
        atom.x = pos_dist(gen);
        atom.y = pos_dist(gen);
        atom.z = pos_dist(gen);
        atom.vx = vel_dist(gen) * sqrt(kB * T / mass_C);
        atom.vy = vel_dist(gen) * sqrt(kB * T / mass_C);
        atom.vz = vel_dist(gen) * sqrt(kB * T / mass_C);
        atom.charge = charge_dist(gen) * 1.602176634e-19; // Elementary charge
        atom.mass = (atom_type_dist(gen) == 0) ? mass_O : ((atom_type_dist(gen) == 1) ? mass_N : ((atom_type_dist(gen) == 2) ? mass_P : mass_C));
        atom.fx = atom.fy = atom.fz = 0.0; // Initialize forces to zero
        atoms.push_back(atom);
    }
}

// Initialize random bonds between atoms
void initialize_bonds(const std::vector<Particle>& atoms, std::vector<Bond>& bonds)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> atom_index_dist(0, atoms.size() - 1);
    int num_bonds = atoms.size() / 4; // Roughly 25% of atoms are bonded

    for (int i = 0; i < num_bonds; ++i)
    {
        int atom1_index = atom_index_dist(gen);
        int atom2_index = atom_index_dist(gen);
        if (atom1_index != atom2_index)
        {
            bonds.push_back({atom1_index, atom2_index});
        }
    }
}

// Function to compute Lennard-Jones potential and forces
void computeLennardJones(Particle &p1, Particle &p2)
{
    MDSRealType dx = p1.x - p2.x;
    MDSRealType dy = p1.y - p2.y;
    MDSRealType dz = p1.z - p2.z;
    MDSRealType r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6;  // (sigma/r)^6
        MDSRealType r12 = r6 * r6; // (sigma/r)^12

        MDSRealType force_scalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;
        p1.fx += force_scalar * dx;
        p1.fy += force_scalar * dy;
        p1.fz += force_scalar * dz;

        p2.fx -= force_scalar * dx;
        p2.fy -= force_scalar * dy;
        p2.fz -= force_scalar * dz;
    }
}

void computeLennardJones(Particle &p1, Particle &p2, MDSRealType r2)
{
    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6; // (sigma/r)^6
        MDSRealType r12 = r6 * r6; // (sigma/r)^12

        MDSRealType force_scalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;
        MDSRealType dx = p1.x - p2.x;
        MDSRealType dy = p1.y - p2.y;
        MDSRealType dz = p1.z - p2.z;
        //MDSRealType dx = p2.x - p1.x, dy = p2.y - p1.y, dz = p2.z - p1.z;

        p1.fx += force_scalar * dx;
        p1.fy += force_scalar * dy;
        p1.fz += force_scalar * dz;
        p2.fx -= force_scalar * dx;
        p2.fy -= force_scalar * dy;
        p2.fz -= force_scalar * dz;
    }
}

// Function to compute Hooke's law (bond stretching)
void computeBondStretching(Particle &p1, Particle &p2)
{
    MDSRealType dx = p1.x - p2.x;
    MDSRealType dy = p1.y - p2.y;
    MDSRealType dz = p1.z - p2.z;
    MDSRealType r = std::sqrt(dx * dx + dy * dy + dz * dz);

    MDSRealType force_scalar = -k_bond * (r - r0_bond) / r;
    p1.fx += force_scalar * dx;
    p1.fy += force_scalar * dy;
    p1.fz += force_scalar * dz;

    p2.fx -= force_scalar * dx;
    p2.fy -= force_scalar * dy;
    p2.fz -= force_scalar * dz;
}

// Function to compute Coulomb forces (electrostatic)
void computeCoulomb(Particle &p1, Particle &p2, MDSRealType r2)
{
    if (r2 < r_cutoff * r_cutoff)
    {
        MDSRealType dx = p1.x - p2.x;
        MDSRealType dy = p1.y - p2.y;
        MDSRealType dz = p1.z - p2.z;
        //MDSRealType r2 = dx * dx + dy * dy + dz * dz;
        MDSRealType r = std::sqrt(r2);

        MDSRealType force_scalar = (k_e * p1.charge * p2.charge) / (r2 * r);
        p1.fx += force_scalar * dx;
        p1.fy += force_scalar * dy;
        p1.fz += force_scalar * dz;

        p2.fx -= force_scalar * dx;
        p2.fy -= force_scalar * dy;
        p2.fz -= force_scalar * dz;
    }
}

// Function to compute angle bending forces (3 atoms)
void computeAngleBending(Particle &p1, Particle &p2, Particle &p3)
{
    MDSRealType dx1 = p1.x - p2.x;
    MDSRealType dy1 = p1.y - p2.y;
    MDSRealType dz1 = p1.z - p2.z;
    MDSRealType r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

    MDSRealType dx2 = p3.x - p2.x;
    MDSRealType dy2 = p3.y - p2.y;
    MDSRealType dz2 = p3.z - p2.z;
    MDSRealType r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    MDSRealType cos_theta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (r1 * r2);
    MDSRealType theta = std::acos(cos_theta);
    MDSRealType angle_force = -k_theta * (theta - theta0);

    // TODO: Compute force vectors and apply them to p1, p2, and p3
}

// Function to compute angle bending forces (3 atoms)
void computeAngleBending(Particle &p1, Particle &p2, Particle &p3, MDSRealType k_theta, MDSRealType theta0)
{
    // Vector p1 -> p2
    MDSRealType dx1 = p1.x - p2.x;
    MDSRealType dy1 = p1.y - p2.y;
    MDSRealType dz1 = p1.z - p2.z;
    MDSRealType r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

    // Vector p3 -> p2
    MDSRealType dx2 = p3.x - p2.x;
    MDSRealType dy2 = p3.y - p2.y;
    MDSRealType dz2 = p3.z - p2.z;
    MDSRealType r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    // Compute the cosine of the angle
    MDSRealType cos_theta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (r1 * r2);
    MDSRealType theta = std::acos(cos_theta);  // Angle in radians

    // Compute the force based on angle deviation from equilibrium
    MDSRealType angle_force = -k_theta * (theta - theta0);

    // Compute forces to apply to each particle (force projection along each axis)
    MDSRealType fx1 = angle_force * (dx1 / r1);
    MDSRealType fy1 = angle_force * (dy1 / r1);
    MDSRealType fz1 = angle_force * (dz1 / r1);

    MDSRealType fx3 = angle_force * (dx2 / r2);
    MDSRealType fy3 = angle_force * (dy2 / r2);
    MDSRealType fz3 = angle_force * (dz2 / r2);

    // Apply forces to the atoms
    p1.fx += fx1;
    p1.fy += fy1;
    p1.fz += fz1;

    p3.fx += fx3;
    p3.fy += fy3;
    p3.fz += fz3;

    // Apply equal and opposite force to the central atom (p2)
    p2.fx -= (fx1 + fx3);
    p2.fy -= (fy1 + fy3);
    p2.fz -= (fz1 + fz3);
}

// Function to compute dihedral torsion forces (4 atoms)
void computeDihedralTorsion(Particle &p1, Particle &p2, Particle &p3, Particle &p4, MDSRealType Vn, int n, MDSRealType gamma)
{
    // Compute bond vectors (p2 -> p1, p2 -> p3, p3 -> p4)
    MDSRealType dx21 = p1.x - p2.x;
    MDSRealType dy21 = p1.y - p2.y;
    MDSRealType dz21 = p1.z - p2.z;

    MDSRealType dx23 = p3.x - p2.x;
    MDSRealType dy23 = p3.y - p2.y;
    MDSRealType dz23 = p3.z - p2.z;

    MDSRealType dx34 = p4.x - p3.x;
    MDSRealType dy34 = p4.y - p3.y;
    MDSRealType dz34 = p4.z - p3.z;

    // Calculate the normal vectors to the planes formed by p1-p2-p3 and p2-p3-p4
    MDSRealType nx1 = dy21 * dz23 - dz21 * dy23;
    MDSRealType ny1 = dz21 * dx23 - dx21 * dz23;
    MDSRealType nz1 = dx21 * dy23 - dy21 * dx23;

    MDSRealType nx2 = dy23 * dz34 - dz23 * dy34;
    MDSRealType ny2 = dz23 * dx34 - dx23 * dz34;
    MDSRealType nz2 = dx23 * dy34 - dy23 * dx34;

    // Calculate the magnitude of these normal vectors
    MDSRealType n1_mag = std::sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
    MDSRealType n2_mag = std::sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);

    // Compute the vector between p2 and p3
    MDSRealType dx = p3.x - p2.x;
    MDSRealType dy = p3.y - p2.y;
    MDSRealType dz = p3.z - p2.z;
    MDSRealType r23 = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Calculate the dihedral angle (phi)
    MDSRealType cos_phi = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (n1_mag * n2_mag);
    MDSRealType sin_phi = r23 * (nx1 * dx34 + ny1 * dy34 + nz1 * dz34) / (n1_mag * n2_mag);
    MDSRealType phi = std::atan2(sin_phi, cos_phi);  // Dihedral angle in radians

    // Torsion potential force: V(phi) = 0.5 * Vn * (1 + cos(n * phi - gamma))
    MDSRealType torsion_energy = 0.5 * Vn * (1 + std::cos(n * phi - gamma));
    MDSRealType torsion_force = -0.5 * Vn * n * std::sin(n * phi - gamma);

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
    MDSRealType fx1 = torsion_force * nx1;
    MDSRealType fy1 = torsion_force * ny1;
    MDSRealType fz1 = torsion_force * nz1;

    // Force on atom p4 (affects normal vector 2)
    MDSRealType fx4 = -torsion_force * nx2;
    MDSRealType fy4 = -torsion_force * ny2;
    MDSRealType fz4 = -torsion_force * nz2;

    // Force on atom p2 and p3 (distribute between them equally)
    MDSRealType fx23 = torsion_force * (nx1 - nx2);
    MDSRealType fy23 = torsion_force * (ny1 - ny2);
    MDSRealType fz23 = torsion_force * (nz1 - nz2);

    // Apply forces to the atoms
    p1.fx += fx1;
    p1.fy += fy1;
    p1.fz += fz1;

    p2.fx += 0.5 * fx23;
    p2.fy += 0.5 * fy23;
    p2.fz += 0.5 * fz23;

    p3.fx += 0.5 * fx23;
    p3.fy += 0.5 * fy23;
    p3.fz += 0.5 * fz23;

    p4.fx += fx4;
    p4.fy += fy4;
    p4.fz += fz4;
}

// Velocity Verlet integration
void PerformVelocityVerletIntegration(std::vector<Particle> &particles, const std::vector<Bond>& bonds)
{
    // bont streching produce here error of counting
    /*
    for (const auto &bond : bonds)
    {
        computeBondStretching(particles[bond.atom1_index], particles[bond.atom2_index]);
    }
    */

    for (auto &p : particles)
    {
        // Update velocities (half step)
        p.vx += 0.5 * p.fx / p.mass * dt;
        p.vy += 0.5 * p.fy / p.mass * dt;
        p.vz += 0.5 * p.fz / p.mass * dt;

        // Update positions
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }

    // Compute forces: Lennard-Jones, Coulomb, bonds, angles, dihedrals
    for (size_t i = 0; i < particles.size(); ++i)
    {
        for (size_t j = i + 1; j < particles.size(); ++j)
        {
            MDSRealType r2 = distance_squared(particles[i], particles[j]);
            computeLennardJones(particles[i], particles[j], r2);
            //computeLennardJones(particles[i], particles[j]);
            computeCoulomb(particles[i], particles[j], r2);
        }
    }

    /*
    for (const auto &bond : bonds)
    {
        computeBondStretching(particles[bond.atom1_index], particles[bond.atom2_index]);
    }
    */

    // Update velocities (second half step)
    for (auto &p : particles)
    {
        p.vx += 0.5 * p.fx / p.mass * dt;
        p.vy += 0.5 * p.fy / p.mass * dt;
        p.vz += 0.5 * p.fz / p.mass * dt;
    }
}

// Function to control temperature using Berendsen thermostat
void ApplyBerendsenThermostatToControlTemperature(std::vector<Particle> &particles, MDSRealType current_temp)
{
    MDSRealType lambda = std::sqrt(1 + dt / (target_temp - current_temp) * (target_temp / current_temp - 1));
    for (auto &p : particles)
    {
        p.vx *= lambda;
        p.vy *= lambda;
        p.vz *= lambda;
    }
}

// Function to compute the kinetic temperature of the system
MDSRealType ComputeCurrentTemperature(const std::vector<Particle> &particles)
{
    MDSRealType kinetic_energy = 0.0;
    for (const auto &p : particles) {
        kinetic_energy += 0.5 * p.mass * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    }
    return (2.0 * kinetic_energy) / (3.0 * particles.size());
}

void update_positions_and_velocities(std::vector<Particle>& atoms, const std::vector<std::vector<MDSRealType>>& forces)
{
    for (size_t i = 0; i < atoms.size(); ++i)
    {
        // Update positions
        atoms[i].x += atoms[i].vx * dt + 0.5 * forces[i][0] / atoms[i].mass * dt * dt;
        atoms[i].y += atoms[i].vy * dt + 0.5 * forces[i][1] / atoms[i].mass * dt * dt;
        atoms[i].z += atoms[i].vz * dt + 0.5 * forces[i][2] / atoms[i].mass * dt * dt;

        // Update velocities
        atoms[i].vx += 0.5 * forces[i][0] / atoms[i].mass * dt;
        atoms[i].vy += 0.5 * forces[i][1] / atoms[i].mass * dt;
        atoms[i].vz += 0.5 * forces[i][2] / atoms[i].mass * dt;
    }
}

int main1()
{
    // Number of atoms and bonds
    const int num_atoms = 100;

    // Initialize atoms and bonds
    std::vector<Particle> particles;
    std::vector<Bond> bonds;
    initialize_atoms(particles);
    initialize_bonds(particles, bonds);
    /*
    std::vector<Particle> particles(100); // Create 100 particles randomly
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dist(-5.0, 5.0);

    // Initialize particles
    for (auto &p : particles) {
        p.x = dist(gen);
        p.y = dist(gen);
        p.z = dist(gen);
        p.vx = dist(gen);
        p.vy = dist(gen);
        p.vz = dist(gen);
        p.mass = 1.0;  // All particles have mass = 1 for simplicity
        p.fx = p.fy = p.fz = 0.0;  // Reset forces
    }
    */


    for (UnsignedInt SimulationStep = 0; SimulationStep < 10000; ++SimulationStep)
    {
        PerformVelocityVerletIntegration(particles, bonds);

        MDSRealType CurrentTemperature = ComputeCurrentTemperature(particles);

        ApplyBerendsenThermostatToControlTemperature(particles, CurrentTemperature);

        if (SimulationStep % 100 == 0)
        {
           cout << "Step " << SimulationStep << std::endl;
           for (const auto& atom : particles)
           {
               cout << "Position: (" << atom.x << ", " << atom.y << ", " << atom.z << ")" << " Velocity: (" << atom.vx << ", " << atom.vy << ", " << atom.vz << ")" << endl;
           }
        }
    }

    cout << "Simulation complete." << endl;

    return 0;
}
