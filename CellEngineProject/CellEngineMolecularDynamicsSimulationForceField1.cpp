
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "CellEngineMolecularDynamicsSimulationForceField1.h"

// Constants for the physical simulation
const double kB = 1.380649e-23; // Boltzmann constant (J/K)
const double T = 303.15; // Temperature in Kelvin (30°C)
const double dt = 1e-15; // Time step in seconds (1 fs)
const double e0 = 8.854187817e-12; // Permittivity of free space
const double charge_constant = 8.9875517923e9; // Coulomb constant (1/(4*pi*e0))
const double mass_O = 2.656e-26; // Mass of an oxygen atom (kg)
const double mass_N = 2.326e-26; // Mass of a nitrogen atom (kg)
const double mass_P = 5.14e-26; // Mass of a phosphorus atom (kg)
const double mass_C = 1.994e-26; // Mass of a carbon atom (kg)
const double bond_strength_k = 100.0; // Bond strength constant (N/m)
const double bond_equilibrium_distance = 1.0e-10; // Equilibrium bond distance (1 Å)

const double r_cutoff = 2.5e-10; // Lennard-Jones cutoff distance
const double epsilon = 1.0; // Lennard-Jones potential depth
const double sigma = 1.0e-10; // Lennard-Jones distance parameter


// Constants for Lennard-Jones potential
//const double epsilon = 1.0; // Depth of the potential well
//const double sigma = 1.0;   // Distance at which potential is zero
//const double r_cutoff = 2.5 * sigma; // Cutoff distance for Lennard-Jones interactions

// Constants for Hooke's Law (bond stretching)
const double k_bond = 100.0; // Bond force constant
const double r0_bond = 1.0;  // Equilibrium bond length

// Constants for angle bending
const double k_theta = 50.0;  // Force constant for angle bending
const double theta0 = M_PI / 2; // Equilibrium angle (90 degrees)

// Constants for dihedral torsion
const double V_n = 1.0;
const double n_period = 3.0;
//const double gamma = 0.0;

// Coulomb constant
const double k_e = 8.99e9;  // Coulomb's constant

// Time step and temperature constants
// const double dt = 0.001;  // Time step for the simulation
const double target_temp = 300.0;  // Target temperature for the thermostat

// Particle structure
struct Particle
{
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double fx, fy, fz;   // Force
    double charge;       // Charge (for Coulomb forces)
    double mass;         // Mass
};

struct Bond
{
    int atom1_index;
    int atom2_index;
};

double distance_squared(const Particle &p1, const Particle &p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
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
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    double r2 = dx * dx + dy * dy + dz * dz;

    if (r2 < r_cutoff * r_cutoff)
    {
        double r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6;  // (sigma/r)^6
        double r12 = r6 * r6; // (sigma/r)^12

        double force_scalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;
        p1.fx += force_scalar * dx;
        p1.fy += force_scalar * dy;
        p1.fz += force_scalar * dz;

        p2.fx -= force_scalar * dx;
        p2.fy -= force_scalar * dy;
        p2.fz -= force_scalar * dz;
    }
}

void computeLennardJones(Particle &p1, Particle &p2, double r2)
{
    if (r2 < r_cutoff * r_cutoff)
    {
        double r6 = (sigma * sigma) / r2;
        r6 = r6 * r6 * r6; // (sigma/r)^6
        double r12 = r6 * r6; // (sigma/r)^12

        double force_scalar = 48 * epsilon * (r12 - 0.5 * r6) / r2;
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        //double dx = p2.x - p1.x, dy = p2.y - p1.y, dz = p2.z - p1.z;

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
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    double r = std::sqrt(dx * dx + dy * dy + dz * dz);

    double force_scalar = -k_bond * (r - r0_bond) / r;
    p1.fx += force_scalar * dx;
    p1.fy += force_scalar * dy;
    p1.fz += force_scalar * dz;

    p2.fx -= force_scalar * dx;
    p2.fy -= force_scalar * dy;
    p2.fz -= force_scalar * dz;
}

// Function to compute Coulomb forces (electrostatic)
void computeCoulomb(Particle &p1, Particle &p2, double r2)
{
    if (r2 < r_cutoff * r_cutoff)
    {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        //double r2 = dx * dx + dy * dy + dz * dz;
        double r = std::sqrt(r2);

        double force_scalar = (k_e * p1.charge * p2.charge) / (r2 * r);
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
    double dx1 = p1.x - p2.x;
    double dy1 = p1.y - p2.y;
    double dz1 = p1.z - p2.z;
    double r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

    double dx2 = p3.x - p2.x;
    double dy2 = p3.y - p2.y;
    double dz2 = p3.z - p2.z;
    double r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    double cos_theta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (r1 * r2);
    double theta = std::acos(cos_theta);
    double angle_force = -k_theta * (theta - theta0);

    // TODO: Compute force vectors and apply them to p1, p2, and p3
}

// Function to compute angle bending forces (3 atoms)
void computeAngleBending(Particle &p1, Particle &p2, Particle &p3, double k_theta, double theta0)
{
    // Vector p1 -> p2
    double dx1 = p1.x - p2.x;
    double dy1 = p1.y - p2.y;
    double dz1 = p1.z - p2.z;
    double r1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

    // Vector p3 -> p2
    double dx2 = p3.x - p2.x;
    double dy2 = p3.y - p2.y;
    double dz2 = p3.z - p2.z;
    double r2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

    // Compute the cosine of the angle
    double cos_theta = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (r1 * r2);
    double theta = std::acos(cos_theta);  // Angle in radians

    // Compute the force based on angle deviation from equilibrium
    double angle_force = -k_theta * (theta - theta0);

    // Compute forces to apply to each particle (force projection along each axis)
    double fx1 = angle_force * (dx1 / r1);
    double fy1 = angle_force * (dy1 / r1);
    double fz1 = angle_force * (dz1 / r1);

    double fx3 = angle_force * (dx2 / r2);
    double fy3 = angle_force * (dy2 / r2);
    double fz3 = angle_force * (dz2 / r2);

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
void computeDihedralTorsion(Particle &p1, Particle &p2, Particle &p3, Particle &p4, double Vn, int n, double gamma) {
    // Compute bond vectors (p2 -> p1, p2 -> p3, p3 -> p4)
    double dx21 = p1.x - p2.x;
    double dy21 = p1.y - p2.y;
    double dz21 = p1.z - p2.z;

    double dx23 = p3.x - p2.x;
    double dy23 = p3.y - p2.y;
    double dz23 = p3.z - p2.z;

    double dx34 = p4.x - p3.x;
    double dy34 = p4.y - p3.y;
    double dz34 = p4.z - p3.z;

    // Calculate the normal vectors to the planes formed by p1-p2-p3 and p2-p3-p4
    double nx1 = dy21 * dz23 - dz21 * dy23;
    double ny1 = dz21 * dx23 - dx21 * dz23;
    double nz1 = dx21 * dy23 - dy21 * dx23;

    double nx2 = dy23 * dz34 - dz23 * dy34;
    double ny2 = dz23 * dx34 - dx23 * dz34;
    double nz2 = dx23 * dy34 - dy23 * dx34;

    // Calculate the magnitude of these normal vectors
    double n1_mag = std::sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
    double n2_mag = std::sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);

    // Compute the vector between p2 and p3
    double dx = p3.x - p2.x;
    double dy = p3.y - p2.y;
    double dz = p3.z - p2.z;
    double r23 = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Calculate the dihedral angle (phi)
    double cos_phi = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (n1_mag * n2_mag);
    double sin_phi = r23 * (nx1 * dx34 + ny1 * dy34 + nz1 * dz34) / (n1_mag * n2_mag);
    double phi = std::atan2(sin_phi, cos_phi);  // Dihedral angle in radians

    // Torsion potential force: V(phi) = 0.5 * Vn * (1 + cos(n * phi - gamma))
    double torsion_energy = 0.5 * Vn * (1 + std::cos(n * phi - gamma));
    double torsion_force = -0.5 * Vn * n * std::sin(n * phi - gamma);

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
    double fx1 = torsion_force * nx1;
    double fy1 = torsion_force * ny1;
    double fz1 = torsion_force * nz1;

    // Force on atom p4 (affects normal vector 2)
    double fx4 = -torsion_force * nx2;
    double fy4 = -torsion_force * ny2;
    double fz4 = -torsion_force * nz2;

    // Force on atom p2 and p3 (distribute between them equally)
    double fx23 = torsion_force * (nx1 - nx2);
    double fy23 = torsion_force * (ny1 - ny2);
    double fz23 = torsion_force * (nz1 - nz2);

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
void velocityVerlet(std::vector<Particle> &particles, const std::vector<Bond>& bonds)
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
            double r2 = distance_squared(particles[i], particles[j]);
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
void berendsenThermostat(std::vector<Particle> &particles, double current_temp)
{
    double lambda = std::sqrt(1 + dt / (target_temp - current_temp) * (target_temp / current_temp - 1));
    for (auto &p : particles)
    {
        p.vx *= lambda;
        p.vy *= lambda;
        p.vz *= lambda;
    }
}

// Function to compute the kinetic temperature of the system
double computeTemperature(const std::vector<Particle> &particles)
{
    double kinetic_energy = 0.0;
    for (const auto &p : particles) {
        kinetic_energy += 0.5 * p.mass * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    }
    return (2.0 * kinetic_energy) / (3.0 * particles.size());
}

void update_positions_and_velocities(std::vector<Particle>& atoms, const std::vector<std::vector<double>>& forces)
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


    // Simulation loop
    for (int step = 0; step < 10000; ++step)
    {
        // Perform velocity Verlet integration
        velocityVerlet(particles, bonds);

        // Compute the current temperature
        double current_temp = computeTemperature(particles);

        // Apply thermostat to control temperature
        berendsenThermostat(particles, current_temp);

        // Periodically print atom positions (for debugging)
        if (step % 100 == 0)
        {
           std::cout << "Step " << step << std::endl;
           for (const auto& atom : particles)
           {
               std::cout << "Position: (" << atom.x << ", " << atom.y << ", " << atom.z << ")" << " Velocity: (" << atom.vx << ", " << atom.vy << ", " << atom.vz << ")" << std::endl;
           }
        }

        // Output data or track particle states for visualization (optional)
    }

    std::cout << "Simulation complete.\n";

    return 0;
}
