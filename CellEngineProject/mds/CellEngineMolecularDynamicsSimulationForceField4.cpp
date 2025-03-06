

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

#include "CellEngineMolecularDynamicsSimulationForceField4.h"

// Constants
const double kB = 1.380649e-23; // Boltzmann constant (J/K)
const double e_charge = 1.602176634e-19; // Elementary charge (C)
const double epsilon0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
const double pi = M_PI;

// Parameters for Lennard-Jones potential
const double sigma = 3.4e-10; // Ångström to meters
const double epsilon = 1.65e-21; // Joules

// Parameters for Hooke's law (bonded interactions)
const double k_bond = 500.0; // Bond spring constant (N/m)
const double r0 = 1.5e-10; // Equilibrium bond length (m)

// Parameters for angle bending
const double k_angle = 100.0; // Angle spring constant (J/rad^2)
const double theta0 = pi / 2.0; // Equilibrium angle (radians)

// Parameters for dihedral torsion
const double k_dihedral = 10.0; // Dihedral spring constant (J)
const double n = 3; // Periodicity of the dihedral potential

// Time step for simulation
const double dt = 1e-15; // Time step (1 fs)

// Structure to represent an atom
struct Atom
{
    double x, y, z; // Position (m)
    double vx, vy, vz; // Velocity (m/s)
    double fx, fy, fz; // Force (N)
    double Charge; // Electric charge (C)
    double Mass; // Mass (kg)
};

// Structure to represent a bond between two atoms
struct Bond
{
    int atom1, atom2; // Indices of bonded atoms
    double k_bond, r0;
};

// Structure to represent an angle between three atoms
struct Angle
{
    int atom1, atom2, atom3; // Indices of atoms forming the angle
};

// Structure to represent a dihedral torsion between four atoms
struct Dihedral
{
    int atom1, atom2, atom3, atom4; // Indices of atoms forming the dihedral
};

// Function to compute the Lennard-Jones force
void lennard_jones_force(const Atom& a, const Atom& b, double& fx, double& fy, double& fz)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    double r = sqrt(r2);
    double r6 = pow(sigma / r, 6);
    double r12 = r6 * r6;
    double force = 24 * epsilon * (2 * r12 - r6) / r;
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the Coulomb force
void coulomb_force(const Atom& a, const Atom& b, double& fx, double& fy, double& fz)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    double r2 = dx * dx + dy * dy + dz * dz;
    double r = sqrt(r2);
    double force = (a.Charge * b.Charge) / (4 * pi * epsilon0 * r * r);
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the Hooke's law force
void hooke_force(const Atom& a, const Atom& b, double& fx, double& fy, double& fz)
{
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    double force = -k_bond * (r - r0);
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the angle bending force
void angle_bending_force(const Atom& a, const Atom& b, const Atom& c, double& fx, double& fy, double& fz)
{
    double abx = a.x - b.x, aby = a.y - b.y, abz = a.z - b.z;
    double cbx = c.x - b.x, cby = c.y - b.y, cbz = c.z - b.z;
    double dot = abx * cbx + aby * cby + abz * cbz;
    double mag_ab = sqrt(abx * abx + aby * aby + abz * abz);
    double mag_cb = sqrt(cbx * cbx + cby * cby + cbz * cbz);
    double theta = acos(dot / (mag_ab * mag_cb));
    double force = -k_angle * (theta - theta0);

    // Compute forces on atoms
    double dtheta_dr1 = (cos(theta) * abx - cbx) / (mag_ab * mag_cb);
    double dtheta_dr3 = (cos(theta) * cbx - abx) / (mag_ab * mag_cb);
    fx += force * dtheta_dr1;
    fy += force * dtheta_dr1;
    fz += force * dtheta_dr1;
    fx += force * dtheta_dr3;
    fy += force * dtheta_dr3;
    fz += force * dtheta_dr3;
}
















// Function to compute the dihedral torsion force
void dihedral_torsion_force(const Atom& a, const Atom& b, const Atom& c, const Atom& d, double& fx, double& fy, double& fz)
{
    // Vectors for the dihedral
    double b1x = b.x - a.x, b1y = b.y - a.y, b1z = b.z - a.z;
    double b2x = c.x - b.x, b2y = c.y - b.y, b2z = c.z - b.z;
    double b3x = d.x - c.x, b3y = d.y - c.y, b3z = d.z - c.z;

    // Cross products
    double n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    double n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize cross products
    double n1_mag = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    double n2_mag = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_mag; n1y /= n1_mag; n1z /= n1_mag;
    n2x /= n2_mag; n2y /= n2_mag; n2z /= n2_mag;

    // Compute dihedral angle (phi)
    double cos_phi = n1x * n2x + n1y * n2y + n1z * n2z;
    double sin_phi = (b2x * (n1y * n2z - n1z * n2y) + b2y * (n1z * n2x - n1x * n2z) + b2z * (n1x * n2y - n1y * n2x)) / (b2x * b2x + b2y * b2y + b2z * b2z);
    double phi = atan2(sin_phi, cos_phi);

    // Compute force magnitude
    double force_magnitude = -k_dihedral * n * sin(n * phi);

    // Compute forces on atoms
    //double dphi_dr1x = ...; // Derivative of phi with respect to atom1's position
    //double dphi_dr1y = ...;
    //double dphi_dr1z = ...;
    //double dphi_dr4x = ...; // Derivative of phi with respect to atom4's position
    //double dphi_dr4y = ...;
    //double dphi_dr4z = ...;

    //fx += force_magnitude * dphi_dr1x;
    //fy += force_magnitude * dphi_dr1y;
    //fz += force_magnitude * dphi_dr1z;
    //fx += force_magnitude * dphi_dr4x;
    //fy += force_magnitude * dphi_dr4y;
    //fz += force_magnitude * dphi_dr4z;
}


































/*
// Function to calculate Lennard-Jones potential
double lj_potential(double r)
{
    double sr6 = pow(sigma / r, 6);
    double sr12 = sr6 * sr6;
    return 4 * epsilon * (sr12 - sr6);
}

// Function to calculate Lennard-Jones force
double lj_force(double r)
{
    double sr6 = pow(sigma / r, 6);
    double sr12 = sr6 * sr6;
    return 24 * epsilon * (2 * sr12 - sr6) / r;
}

// Function to calculate Coulomb potential
double coulomb_potential(double q1, double q2, double r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0 * r);
}

// Function to calculate Coulomb force
double coulomb_force(double q1, double q2, double r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0 * r * r);
}

// Function to calculate angle bending potential
double angle_potential(double theta, double theta0)
{
    return 0.5 * k_angle * pow(theta - theta0, 2);
}

// Function to calculate angle bending force
double angle_force(double theta, double theta0)
{
    return k_angle * (theta - theta0);
}

// Function to calculate dihedral torsion potential
double dihedral_potential(double phi)
{
    return k_dihedral * (1 + cos(3 * phi - M_PI));
}

// Function to calculate dihedral torsion force
double dihedral_force(double phi)
{
    return -3 * k_dihedral * sin(3 * phi - M_PI);
}

// Function to calculate Hooke's law potential
double bond_potential(double r, double r0)
{
    return 0.5 * k_bond * pow(r - r0, 2);
}

// Function to calculate Hooke's law force
double bond_force(double r, double r0)
{
    return k_bond * (r - r0);
}

// Function to calculate the distance between two atoms
double distance(const Atom& a1, const Atom& a2) {
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}
// Function to calculate the angle between three atoms
double calculate_angle(const Atom& a1, const Atom& a2, const Atom& a3)
{
    double r12 = distance(a1, a2);
    double r23 = distance(a2, a3);
    double r13 = distance(a1, a3);
    return acos((r12 * r12 + r23 * r23 - r13 * r13) / (2 * r12 * r23));
}
*/










// Function to calculate the dihedral angle between four atoms
double calculate_dihedral(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4)
{
    // Vectors between atoms
    double b1x = a2.x - a1.x, b1y = a2.y - a1.y, b1z = a2.z - a1.z;
    double b2x = a3.x - a2.x, b2y = a3.y - a2.y, b2z = a3.z - a2.z;
    double b3x = a4.x - a3.x, b3y = a4.y - a3.y, b3z = a4.z - a3.z;

    // Cross products
    double n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    double n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize vectors
    double n1_norm = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    double n2_norm = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_norm; n1y /= n1_norm; n1z /= n1_norm;
    n2x /= n2_norm; n2y /= n2_norm; n2z /= n2_norm;

    // Dihedral angle
    double dot = n1x * n2x + n1y * n2y + n1z * n2z;
    double cross_x = n1y * n2z - n1z * n2y;
    double cross_y = n1z * n2x - n1x * n2z;
    double cross_z = n1x * n2y - n1y * n2x;
    double cross_norm = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    double phi = atan2(cross_norm, dot);

    return phi;
}

// Calculate bond stretching forces
void calculate_bond_forces(const std::vector<Bond>& bonds, std::vector<Atom>& atoms)
{
    for (const auto& bond : bonds)
    {
        Atom& a1 = atoms[bond.atom1];
        Atom& a2 = atoms[bond.atom2];
        double dx = a2.x - a1.x;
        double dy = a2.y - a1.y;
        double dz = a2.z - a1.z;
        double r = sqrt(dx * dx + dy * dy + dz * dz);
        double dr = r - bond.r0;
        double force_magnitude = -bond.k_bond * dr;
        double fx = force_magnitude * dx / r;
        double fy = force_magnitude * dy / r;
        double fz = force_magnitude * dz / r;
        a1.fx -= fx;
        a1.fy -= fy;
        a1.fz -= fz;
        a2.fx += fx;
        a2.fy += fy;
        a2.fz += fz;
    }
}


























// Main simulation function
void simulate(std::vector<Atom>& atoms, const std::vector<Bond>& bonds, const std::vector<Angle>& angles, const std::vector<Dihedral>& dihedrals, int steps)
{
    for (int step = 0; step < steps; ++step)
    {
        // Reset forces
        for (Atom& atom : atoms)
        {
            atom.fx = atom.fy = atom.fz = 0.0;
        }

        // Compute non-bonded forces (Lennard-Jones and Coulomb)
        for (size_t i = 0; i < atoms.size(); ++i)
        {
            for (size_t j = i + 1; j < atoms.size(); ++j)
            {
                lennard_jones_force(atoms[i], atoms[j], atoms[i].fx, atoms[i].fy, atoms[i].fz);
                coulomb_force(atoms[i], atoms[j], atoms[i].fx, atoms[i].fy, atoms[i].fz);
            }
        }

        // Compute bonded forces (Hooke's law)
        for (const Bond& bond : bonds)
        {
            hooke_force(atoms[bond.atom1], atoms[bond.atom2], atoms[bond.atom1].fx, atoms[bond.atom1].fy, atoms[bond.atom1].fz);
        }

        // Compute angle bending forces
        for (const Angle& angle : angles)
        {
            angle_bending_force(atoms[angle.atom1], atoms[angle.atom2], atoms[angle.atom3], atoms[angle.atom1].fx, atoms[angle.atom1].fy, atoms[angle.atom1].fz);
        }

        // Compute dihedral torsion forces
        for (const Dihedral& dihedral : dihedrals)
        {
            dihedral_torsion_force(atoms[dihedral.atom1], atoms[dihedral.atom2], atoms[dihedral.atom3], atoms[dihedral.atom4], atoms[dihedral.atom1].fx, atoms[dihedral.atom1].fy, atoms[dihedral.atom1].fz);
        }



        // Update positions and velocities (Velocity Verlet algorithm)
        for (Atom& atom : atoms)
        {
            atom.vx += atom.fx / atom.Mass * dt;
            atom.vy += atom.fy / atom.Mass * dt;
            atom.vz += atom.fz / atom.Mass * dt;
            atom.x += atom.vx * dt;
            atom.y += atom.vy * dt;
            atom.z += atom.vz * dt;
        }
    }
}

int main4()
{
    // Random number generator
    std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> pos_dist(-1e-9, 1e-9); // Positions in nm
    std::uniform_real_distribution<double> vel_dist(-1e3, 1e3); // Velocities in m/s
    std::uniform_real_distribution<double> charge_dist(-1.0, 1.0); // Charges in elementary charge units

    // Create 1000 atoms
    std::vector<Atom> atoms(1000);
    for (Atom& atom : atoms)
    {
        atom.x = pos_dist(rng);
        atom.y = pos_dist(rng);
        atom.z = pos_dist(rng);
        atom.vx = vel_dist(rng);
        atom.vy = vel_dist(rng);
        atom.vz = vel_dist(rng);
        atom.Charge = charge_dist(rng) * e_charge;
        atom.Mass = 1.67e-27; // Approximate mass of a proton (kg)
    }

    // Create bonds, angles, and dihedrals (example)
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;

    // Run simulation for 1000 steps
    simulate(atoms, bonds, angles, dihedrals, 1000);

    return 0;
}

