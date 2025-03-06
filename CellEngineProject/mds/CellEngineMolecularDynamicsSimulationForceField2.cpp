
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <algorithm>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField2.h"
#include "CellEngineMolecularDynamicsSimulationForceField2Types.h"

class Simulation
{
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::unordered_set<UnsignedInt> bonded_pairs;
    MDSRealType box_size, temperature, time_step;
    const MDSRealType k_coulomb = 332.0636, k_boltzmann = 0.0019872041;

    void computeBondForces()
    {
        for (const auto& b : bonds)
        {
            Atom &a = atoms[b.i], &b_atom = atoms[b.j];
            Vec3 delta = b_atom.Position - a.Position;
            MDSRealType r = delta.length();
            MDSRealType force_mag = -b.k * (r - b.r0) / r;
            Vec3 force = delta * force_mag;
            a.Force += force;
            b_atom.Force -= force;
        }
    }

    void computeAngleForces()
    {
        for (const auto& a : angles)
        {
            Atom &i = atoms[a.i], &j = atoms[a.j], &k = atoms[a.k];
            Vec3 r_ij = i.Position - j.Position, r_kj = k.Position - j.Position;
            MDSRealType len_ij = r_ij.length(), len_kj = r_kj.length();
            MDSRealType cos_theta = r_ij.dot(r_kj) / (len_ij * len_kj);
            MDSRealType theta = std::acos(cos_theta);
            MDSRealType dtheta = theta - a.theta0;
            MDSRealType dE_dtheta = a.k_theta * dtheta;

            Vec3 dcos_di = (r_kj * len_ij - r_ij * (r_ij.dot(r_kj)/len_ij)) / (len_ij*len_ij*len_kj);
            Vec3 dcos_dk = (r_ij * len_kj - r_kj * (r_ij.dot(r_kj)/len_kj)) / (len_kj*len_kj*len_ij);
            Vec3 dcos_dj = dcos_di + dcos_dk;
            dcos_dj *= -1.0f;

            Vec3 force_i = dcos_di * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_j = dcos_dj * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_k = dcos_dk * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));

            i.Force += force_i;
            j.Force += force_j;
            k.Force += force_k;
        }
    }

    void computeDihedralForces()
    {
        for (const auto& d : dihedrals)
        {
            Atom &i = atoms[d.i], &j = atoms[d.j], &k = atoms[d.k], &l = atoms[d.l];
            Vec3 r_ij = i.Position - j.Position, r_kj = k.Position - j.Position, r_kl = l.Position - k.Position;
            Vec3 n1 = r_ij.cross(r_kj), n2 = r_kj.cross(r_kl);
            MDSRealType n1_len = n1.length(), n2_len = n2.length();
            if (n1_len == 0 || n2_len == 0) continue;

            MDSRealType cos_phi = n1.dot(n2) / (n1_len * n2_len);
            MDSRealType phi = std::acos(cos_phi);
            MDSRealType sign = n1.cross(n2).dot(r_kj) < 0 ? -1 : 1;
            phi *= sign;

            MDSRealType dV_dphi = d.k_phi * d.n * std::sin(d.n * phi - d.phi0);

            Vec3 f1 = (n1 * (d.n * d.k_phi / n1_len)).cross(r_kj);
            Vec3 f4 = (n2 * (d.n * d.k_phi / n2_len)).cross(r_kj);
            Vec3 f2 = (f1.cross(r_ij) - f4.cross(r_kl)) / r_kj.dot(r_kj);
            Vec3 f3 = (f2 + f4) * -1.0;

            i.Force += f1 * dV_dphi;
            j.Force += (f1 * -1.0 + f2) * dV_dphi;
            k.Force += (f3 + f4) * -1.0 * dV_dphi;
            l.Force += f4 * dV_dphi;
        }
    }

    void computeNonBondedForces()
    {
        size_t n = atoms.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = i+1; j < n; ++j)
            {
                if (bonded_pairs.count(i*n + j))
                   continue;

                Atom &a = atoms[i], &b = atoms[j];
                Vec3 delta = a.Position - b.Position;
                MDSRealType r = delta.length();
                if (r == 0)
                   continue;

                // Coulomb interaction
                MDSRealType coulomb = k_coulomb * a.Charge * b.Charge / (r*r*r);
                Vec3 f_coulomb = delta * coulomb;
                a.Force += f_coulomb;
                b.Force -= f_coulomb;

                // Lennard-Jones interaction
                MDSRealType sigma_avg = (a.SigmaValue + b.SigmaValue)/2;
                MDSRealType epsilon_avg = std::sqrt(a.EpsilonValue * b.EpsilonValue);
                MDSRealType sr6 = std::pow(sigma_avg/r, 6);
                MDSRealType lj = 24 * epsilon_avg * (2*sr6*sr6 - sr6) / (r*r);
                Vec3 f_lj = delta * lj;
                a.Force += f_lj;
                b.Force -= f_lj;
            }
        }
    }

public:
    Simulation(MDSRealType box, MDSRealType temp, MDSRealType dt) : box_size(box), temperature(temp), time_step(dt)
    {
    }

    void initializeAtoms1(UnsignedInt n)
    {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<> pos(0, box_size), charge(-1, 1);
        std::uniform_int_distribution<> AtomType(0, 3);

        for (UnsignedInt i = 0; i < n; ++i)
        {
            atoms.emplace_back(static_cast<AtomTypes>(AtomType(gen)), Vec3(pos(gen), pos(gen), pos(gen)), charge(gen));
            MDSRealType scale = std::sqrt(3 * k_boltzmann * temperature / atoms.back().Mass);
            std::normal_distribution<> vel(0, scale);
            atoms.back().Velocity = Vec3(vel(gen), vel(gen), vel(gen));
        }

        // Create molecular topology (example: linear chain)
        for (UnsignedInt i = 0; i < n/2; ++i)
        {
            bonds.emplace_back(i, i+1, 1.5, 100);
            bonded_pairs.insert(i*n + i+1);
            if (i < n/2 - 2)
            {
                angles.emplace_back(i, i+1, i+2, 1.5708, 100);
                dihedrals.emplace_back(i, i+1, i+2, i+3, 5, 2, 0);
            }
        }
    }

    void initializeAtoms2(int num_atoms)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> pos_dist(0.0, box_size);
        std::uniform_int_distribution<> type_dist(0, 3); // 0=O, 1=N, 2=P, 3=C
        std::uniform_real_distribution<> charge_dist(-1.0, 1.0);

        for (UnsignedInt i = 0; i < num_atoms; ++i)
        {
            Vec3 pos(pos_dist(gen), pos_dist(gen), pos_dist(gen));
            UnsignedInt type = type_dist(gen);
            MDSRealType charge = charge_dist(gen);
            atoms.emplace_back(static_cast<AtomTypes>(type), pos, charge);
        }

        // Assign random velocities based on Maxwell-Boltzmann distribution
        std::normal_distribution<> vel_dist(0.0, std::sqrt(k_boltzmann * temperature));
        for (auto& atom : atoms)
        {
            atom.Velocity.x = vel_dist(gen);
            atom.Velocity.y = vel_dist(gen);
            atom.Velocity.z = vel_dist(gen);
            // Subtract center-of-mass velocity to avoid drift
            // (omitted for brevity)
        }

        // Generate random bonds (example: 500 bonds)
        std::uniform_int_distribution<> atom_dist(0, num_atoms - 1);

        for (UnsignedInt b = 0; b < 500; ++b)
        {
            UnsignedInt i = atom_dist(gen);
            UnsignedInt j = atom_dist(gen);
            if (i != j)
            {
                // Ensure i < j to avoid duplicates
                if (i > j)
                   std::swap(i, j);
                // Check if bond already exists
                UnsignedInt key = i * num_atoms + j;
                if (bonded_pairs.find(key) == bonded_pairs.end())
                {
                    bonded_pairs.insert(key);
                    // Assign random r0 and k (example values)
                    MDSRealType r0 = 1.5; // Angstrom
                    MDSRealType k = 100.0; // kcal/(mol·Å²)
                    //bonds.emplace_back(i, j, r0, k);
                }
            }
        }

        // Generate angles and dihedrals (example)
        // For angles, find triplets i-j-k where i-j and j-k are bonded
        // This is a simplified approach; realistically, need to find connected triplets
        // For demonstration, randomly select triplets

        for (UnsignedInt a = 0; a < 300; ++a)
        {
            UnsignedInt i = atom_dist(gen);
            UnsignedInt j = atom_dist(gen);
            UnsignedInt k = atom_dist(gen);
            if (i != j && j != k && i != k) {
                angles.emplace_back(i, j, k,
                1.5708 // 90 degree
                , 100.0 // k_theta
                );
            }
        }

        // Similarly for dihedrals
        for (UnsignedInt d = 0; d < 200; ++d)
        {
            UnsignedInt i = atom_dist(gen);
            UnsignedInt j = atom_dist(gen);
            UnsignedInt k = atom_dist(gen);
            UnsignedInt l = atom_dist(gen);
            if (i != j && j != k && k != l && i != k && i != l && j != l) {
                dihedrals.emplace_back(i, j, k, l,
                5.0 // k_phi
                ,
                2 // n
                ,
                0.0 // phi0
                );
            }
        }

    }

    void computeForces()
    {
        for (auto& a : atoms)
            a.Force = Vec3();

        computeBondForces();
        computeAngleForces();
        computeDihedralForces();
        computeNonBondedForces();
    }

    void integrate()
    {
        // Velocity half-step
        for (auto& a : atoms)
        {
            a.Velocity += a.Force * (time_step / (2 * a.Mass));
        }

        // Position update with periodic boundary conditions
        for (auto& a : atoms)
        {
            a.Position += a.Velocity * time_step;
            a.Position.x = fmod(a.Position.x + box_size, box_size);
            a.Position.y = fmod(a.Position.y + box_size, box_size);
            a.Position.z = fmod(a.Position.z + box_size, box_size);

            if (a.Position.x < 0) a.Position.x += box_size;
            if (a.Position.y < 0) a.Position.y += box_size;
            if (a.Position.z < 0) a.Position.z += box_size;
        }

        // Compute new forces
        computeForces();

        // Velocity second half-step
        for (auto& a : atoms)
        {
            a.Velocity += a.Force * (time_step / (2 * a.Mass));
        }
    }

    void applyBerendsenThermostat(MDSRealType target_temp, MDSRealType tau)
    {
        MDSRealType current_temp = 0;
        for (const auto& a : atoms)
        {
            current_temp += a.Mass * a.Velocity.dot(a.Velocity);
        }
        current_temp /= (3 * atoms.size() * k_boltzmann);

        MDSRealType lambda = std::sqrt(1 + time_step/tau * (target_temp/current_temp - 1));
        for (auto& a : atoms)
        {
            a.Velocity = a.Velocity * lambda;
        }
    }
};

UnsignedInt main2()
{
    const MDSRealType box_size = 100.0;    // Angstroms
    const MDSRealType temperature = 303.15; // 30°C in Kelvin
    const MDSRealType time_step = 0.001;    // picoseconds
    const UnsignedInt num_steps = 1000;
    const UnsignedInt output_interval = 100;

    Simulation sim(box_size, temperature, time_step);
    sim.initializeAtoms1(1000);
    sim.computeForces();

    for (UnsignedInt step = 0; step < num_steps; ++step)
    {
        sim.integrate();
        sim.applyBerendsenThermostat(temperature, 0.1);

        if (step % output_interval == 0)
        {
            std::cout << "Step " << step << " completed" << std::endl;
        }
    }

    return 0;
}

