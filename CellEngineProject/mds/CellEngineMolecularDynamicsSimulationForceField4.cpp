
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

// Function to compute the Lennard-Jones force
void lennard_jones_force(const Atom& a, const Atom& b, MDSRealType& fx, MDSRealType& fy, MDSRealType& fz)
{
    MDSRealType dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    MDSRealType r = sqrt(r2);
    MDSRealType r6 = pow(sigma / r, 6);
    MDSRealType r12 = r6 * r6;
    MDSRealType force = 24 * epsilon * (2 * r12 - r6) / r;
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the Coulomb force
void coulomb_force(const Atom& a, const Atom& b, MDSRealType& fx, MDSRealType& fy, MDSRealType& fz)
{
    MDSRealType dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    MDSRealType r2 = dx * dx + dy * dy + dz * dz;
    MDSRealType r = sqrt(r2);
    MDSRealType force = (a.Charge * b.Charge) / (4 * pi * epsilon0 * r * r);
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the Hooke's law force
void hooke_force(const Atom& a, const Atom& b, MDSRealType& fx, MDSRealType& fy, MDSRealType& fz)
{
    MDSRealType dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
    MDSRealType force = -k_bond * (r - r0);
    fx += force * dx / r;
    fy += force * dy / r;
    fz += force * dz / r;
}

// Function to compute the angle bending force
void angle_bending_force(const Atom& a, const Atom& b, const Atom& c, MDSRealType& fx, MDSRealType& fy, MDSRealType& fz)
{
    MDSRealType abx = a.x - b.x, aby = a.y - b.y, abz = a.z - b.z;
    MDSRealType cbx = c.x - b.x, cby = c.y - b.y, cbz = c.z - b.z;
    MDSRealType dot = abx * cbx + aby * cby + abz * cbz;
    MDSRealType mag_ab = sqrt(abx * abx + aby * aby + abz * abz);
    MDSRealType mag_cb = sqrt(cbx * cbx + cby * cby + cbz * cbz);
    MDSRealType theta = acos(dot / (mag_ab * mag_cb));
    MDSRealType force = -k_angle * (theta - theta0);

    // Compute forces on atoms
    MDSRealType dtheta_dr1 = (cos(theta) * abx - cbx) / (mag_ab * mag_cb);
    MDSRealType dtheta_dr3 = (cos(theta) * cbx - abx) / (mag_ab * mag_cb);
    fx += force * dtheta_dr1;
    fy += force * dtheta_dr1;
    fz += force * dtheta_dr1;
    fx += force * dtheta_dr3;
    fy += force * dtheta_dr3;
    fz += force * dtheta_dr3;
}
















// Function to compute the dihedral torsion force
void dihedral_torsion_force(const Atom& a, const Atom& b, const Atom& c, const Atom& d, MDSRealType& fx, MDSRealType& fy, MDSRealType& fz)
{
    // Vectors for the dihedral
    MDSRealType b1x = b.x - a.x, b1y = b.y - a.y, b1z = b.z - a.z;
    MDSRealType b2x = c.x - b.x, b2y = c.y - b.y, b2z = c.z - b.z;
    MDSRealType b3x = d.x - c.x, b3y = d.y - c.y, b3z = d.z - c.z;

    // Cross products
    MDSRealType n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    MDSRealType n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize cross products
    MDSRealType n1_mag = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    MDSRealType n2_mag = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_mag; n1y /= n1_mag; n1z /= n1_mag;
    n2x /= n2_mag; n2y /= n2_mag; n2z /= n2_mag;

    // Compute dihedral angle (phi)
    MDSRealType cos_phi = n1x * n2x + n1y * n2y + n1z * n2z;
    MDSRealType sin_phi = (b2x * (n1y * n2z - n1z * n2y) + b2y * (n1z * n2x - n1x * n2z) + b2z * (n1x * n2y - n1y * n2x)) / (b2x * b2x + b2y * b2y + b2z * b2z);
    MDSRealType phi = atan2(sin_phi, cos_phi);

    // Compute force magnitude
    MDSRealType force_magnitude = -k_dihedral * n * sin(n * phi);

    // Compute forces on atoms
    //MDSRealType dphi_dr1x = ...; // Derivative of phi with respect to atom1's position
    //MDSRealType dphi_dr1y = ...;
    //MDSRealType dphi_dr1z = ...;
    //MDSRealType dphi_dr4x = ...; // Derivative of phi with respect to atom4's position
    //MDSRealType dphi_dr4y = ...;
    //MDSRealType dphi_dr4z = ...;

    //fx += force_magnitude * dphi_dr1x;
    //fy += force_magnitude * dphi_dr1y;
    //fz += force_magnitude * dphi_dr1z;
    //fx += force_magnitude * dphi_dr4x;
    //fy += force_magnitude * dphi_dr4y;
    //fz += force_magnitude * dphi_dr4z;
}


































/*
// Function to calculate Lennard-Jones potential
MDSRealType lj_potential(MDSRealType r)
{
    MDSRealType sr6 = pow(sigma / r, 6);
    MDSRealType sr12 = sr6 * sr6;
    return 4 * epsilon * (sr12 - sr6);
}

// Function to calculate Lennard-Jones force
MDSRealType lj_force(MDSRealType r)
{
    MDSRealType sr6 = pow(sigma / r, 6);
    MDSRealType sr12 = sr6 * sr6;
    return 24 * epsilon * (2 * sr12 - sr6) / r;
}

// Function to calculate Coulomb potential
MDSRealType coulomb_potential(MDSRealType q1, MDSRealType q2, MDSRealType r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0 * r);
}

// Function to calculate Coulomb force
MDSRealType coulomb_force(MDSRealType q1, MDSRealType q2, MDSRealType r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0 * r * r);
}

// Function to calculate angle bending potential
MDSRealType angle_potential(MDSRealType theta, MDSRealType theta0)
{
    return 0.5 * k_angle * pow(theta - theta0, 2);
}

// Function to calculate angle bending force
MDSRealType angle_force(MDSRealType theta, MDSRealType theta0)
{
    return k_angle * (theta - theta0);
}

// Function to calculate dihedral torsion potential
MDSRealType dihedral_potential(MDSRealType phi)
{
    return k_dihedral * (1 + cos(3 * phi - M_PI));
}

// Function to calculate dihedral torsion force
MDSRealType dihedral_force(MDSRealType phi)
{
    return -3 * k_dihedral * sin(3 * phi - M_PI);
}

// Function to calculate Hooke's law potential
MDSRealType bond_potential(MDSRealType r, MDSRealType r0)
{
    return 0.5 * k_bond * pow(r - r0, 2);
}

// Function to calculate Hooke's law force
MDSRealType bond_force(MDSRealType r, MDSRealType r0)
{
    return k_bond * (r - r0);
}

// Function to calculate the distance between two atoms
MDSRealType distance(const Atom& a1, const Atom& a2) {
    MDSRealType dx = a1.x - a2.x;
    MDSRealType dy = a1.y - a2.y;
    MDSRealType dz = a1.z - a2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}
// Function to calculate the angle between three atoms
MDSRealType calculate_angle(const Atom& a1, const Atom& a2, const Atom& a3)
{
    MDSRealType r12 = distance(a1, a2);
    MDSRealType r23 = distance(a2, a3);
    MDSRealType r13 = distance(a1, a3);
    return acos((r12 * r12 + r23 * r23 - r13 * r13) / (2 * r12 * r23));
}
*/










// Function to calculate the dihedral angle between four atoms
MDSRealType calculate_dihedral(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4)
{
    // Vectors between atoms
    MDSRealType b1x = a2.x - a1.x, b1y = a2.y - a1.y, b1z = a2.z - a1.z;
    MDSRealType b2x = a3.x - a2.x, b2y = a3.y - a2.y, b2z = a3.z - a2.z;
    MDSRealType b3x = a4.x - a3.x, b3y = a4.y - a3.y, b3z = a4.z - a3.z;

    // Cross products
    MDSRealType n1x = b1y * b2z - b1z * b2y, n1y = b1z * b2x - b1x * b2z, n1z = b1x * b2y - b1y * b2x;
    MDSRealType n2x = b2y * b3z - b2z * b3y, n2y = b2z * b3x - b2x * b3z, n2z = b2x * b3y - b2y * b3x;

    // Normalize vectors
    MDSRealType n1_norm = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    MDSRealType n2_norm = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    n1x /= n1_norm; n1y /= n1_norm; n1z /= n1_norm;
    n2x /= n2_norm; n2y /= n2_norm; n2z /= n2_norm;

    // Dihedral angle
    MDSRealType dot = n1x * n2x + n1y * n2y + n1z * n2z;
    MDSRealType cross_x = n1y * n2z - n1z * n2y;
    MDSRealType cross_y = n1z * n2x - n1x * n2z;
    MDSRealType cross_z = n1x * n2y - n1y * n2x;
    MDSRealType cross_norm = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    MDSRealType phi = atan2(cross_norm, dot);

    return phi;
}

// Calculate bond stretching forces
void calculate_bond_forces(const std::vector<Bond>& bonds, std::vector<Atom>& atoms)
{
    for (const auto& bond : bonds)
    {
        Atom& a1 = atoms[bond.Atom1];
        Atom& a2 = atoms[bond.Atom2];
        MDSRealType dx = a2.x - a1.x;
        MDSRealType dy = a2.y - a1.y;
        MDSRealType dz = a2.z - a1.z;
        MDSRealType r = sqrt(dx * dx + dy * dy + dz * dz);
        MDSRealType dr = r - bond.r0;
        MDSRealType force_magnitude = -bond.k_bond * dr;
        MDSRealType fx = force_magnitude * dx / r;
        MDSRealType fy = force_magnitude * dy / r;
        MDSRealType fz = force_magnitude * dz / r;
        a1.fx -= fx;
        a1.fy -= fy;
        a1.fz -= fz;
        a2.fx += fx;
        a2.fy += fy;
        a2.fz += fz;
    }
}


























// Main simulation function
void Simulate(std::vector<Atom>& Atoms, const std::vector<Bond>& Bonds, const std::vector<Angle>& Angles, const std::vector<Dihedral>& Dihedrals, const int Steps)
{
    for (int step = 0; step < Steps; ++step)
    {
        // Reset forces
        for (Atom& atom : Atoms)
            atom.fx = atom.fy = atom.fz = 0.0;

        // Compute non-bonded forces (Lennard-Jones and Coulomb)
        for (size_t i = 0; i < Atoms.size(); ++i)
            for (size_t j = i + 1; j < Atoms.size(); ++j)
            {
                lennard_jones_force(Atoms[i], Atoms[j], Atoms[i].fx, Atoms[i].fy, Atoms[i].fz);
                coulomb_force(Atoms[i], Atoms[j], Atoms[i].fx, Atoms[i].fy, Atoms[i].fz);
            }

        // Compute bonded forces (Hooke's law)
        for (const Bond& bond : Bonds)
            hooke_force(Atoms[bond.Atom1], Atoms[bond.Atom2], Atoms[bond.Atom1].fx, Atoms[bond.Atom1].fy, Atoms[bond.Atom1].fz);

        // Compute angle bending forces
        for (const Angle& angle : Angles)
            angle_bending_force(Atoms[angle.Atom1], Atoms[angle.Atom2], Atoms[angle.Atom3], Atoms[angle.Atom1].fx, Atoms[angle.Atom1].fy, Atoms[angle.Atom1].fz);

        // Compute dihedral torsion forces
        for (const Dihedral& dihedral : Dihedrals)
            dihedral_torsion_force(Atoms[dihedral.Atom1], Atoms[dihedral.Atom2], Atoms[dihedral.Atom3], Atoms[dihedral.Atom4], Atoms[dihedral.Atom1].fx, Atoms[dihedral.Atom1].fy, Atoms[dihedral.Atom1].fz);

        // Update positions and velocities (Velocity Verlet algorithm)
        for (Atom& atom : Atoms)
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
    std::mt19937 RandomNumberGenerator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<MDSRealType> PositionsDistributionInNM(-1e-9, 1e-9); // Positions in nm
    std::uniform_real_distribution<MDSRealType> VelocityDistributionInMPerS(-1e3, 1e3); // Velocities in m/s
    std::uniform_real_distribution<MDSRealType> ChargesDistribution(-1.0, 1.0); // Charges in elementary charge units

    // Create 1000 atoms
    vector<Atom> Atoms(1000);
    for (Atom& atom : Atoms)
    {
        atom.x = PositionsDistributionInNM(RandomNumberGenerator);
        atom.y = PositionsDistributionInNM(RandomNumberGenerator);
        atom.z = PositionsDistributionInNM(RandomNumberGenerator);
        atom.vx = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.vy = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.vz = VelocityDistributionInMPerS(RandomNumberGenerator);
        atom.Charge = ChargesDistribution(RandomNumberGenerator) * e_charge;
        atom.Mass = 1.67e-27; // Approximate mass of a proton (kg)
    }

    // Create bonds, angles, and dihedrals (example)
    vector<Bond> Bonds;
    vector<Angle> Angles;
    vector<Dihedral> Dihedrals;

    // Run simulation for 1000 steps
    Simulate(Atoms, Bonds, Angles, Dihedrals, 1000);

    return 0;
}

