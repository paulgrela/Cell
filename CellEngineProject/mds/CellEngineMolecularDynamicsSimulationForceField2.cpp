
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
private:
    std::vector<Atom> Atoms;
    std::vector<Bond> Bonds;
    std::vector<Angle> Angles;
    std::vector<Dihedral> Dihedrals;
    std::unordered_set<UnsignedInt> BondedPairs;
private:
    MDSRealType BoxSize;
    MDSRealType Temperature;
    MDSRealType TimeStep;
private:
    const MDSRealType K_Coulomb = 332.0636;
    const MDSRealType K_Boltzmann = 0.0019872041;
public:
    void ComputeBondForces()
    {
        for (const auto& Bond : Bonds)
        {
            Atom& Atom1 = Atoms[Bond.Atom1Index];
            Atom& Atom2 = Atoms[Bond.Atom2Index];
            Vec3 delta = Atom2.Position - Atom1.Position;
            MDSRealType r = delta.length();
            MDSRealType force_mag = -Bond.k_SpringConstant * (r - Bond.r0_EquilibrumDistance) / r;
            const Vec3 Force = delta * force_mag;
            Atom1.Force += Force;
            Atom2.Force -= Force;
        }
    }
public:
    void ComputeAngleForces()
    {
        for (const auto& Angle : Angles)
        {
            Atom& Atom1 = Atoms[Angle.Atom1Index];
            Atom& Atom2 = Atoms[Angle.Atom2Index];
            Atom& Atom3 = Atoms[Angle.Atom3Index];

            Vec3 r_ij = Atom1.Position - Atom2.Position;
            Vec3 r_kj = Atom3.Position - Atom2.Position;

            MDSRealType len_ij = r_ij.length(), len_kj = r_kj.length();
            MDSRealType cos_theta = r_ij.dot(r_kj) / (len_ij * len_kj);
            MDSRealType theta = std::acos(cos_theta);
            MDSRealType dtheta = theta - Angle.theta0_EquilibriumAngleInRadians;
            MDSRealType dE_dtheta = Angle.k_theta_ForceConstant * dtheta;

            Vec3 dcos_di = (r_kj * len_ij - r_ij * (r_ij.dot(r_kj)/len_ij)) / (len_ij*len_ij*len_kj);
            Vec3 dcos_dk = (r_ij * len_kj - r_kj * (r_ij.dot(r_kj)/len_kj)) / (len_kj*len_kj*len_ij);
            Vec3 dcos_dj = dcos_di + dcos_dk;
            dcos_dj *= -1.0f;

            Vec3 force_i = dcos_di * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_j = dcos_dj * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_k = dcos_dk * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));

            Atom1.Force += force_i;
            Atom2.Force += force_j;
            Atom3.Force += force_k;
        }
    }
public:
    void ComputeDihedralForces()
    {
        for (const auto& Dihedral : Dihedrals)
        {
            Atom& Atom1 = Atoms[Dihedral.Atom1Index];
            Atom& Atom2 = Atoms[Dihedral.Atom2Index];
            Atom& Atom3 = Atoms[Dihedral.Atom3Index];
            Atom& Atom4 = Atoms[Dihedral.Atom4Index];

            Vec3 r_ij = Atom1.Position - Atom2.Position;
            Vec3 r_kj = Atom3.Position - Atom2.Position;
            Vec3 r_kl = Atom4.Position - Atom3.Position;

            Vec3 n1 = r_ij.cross(r_kj);
            Vec3 n2 = r_kj.cross(r_kl);

            const MDSRealType n1_length = n1.length();
            const MDSRealType n2_length = n2.length();
            if (n1_length == 0 || n2_length == 0)
                continue;

            MDSRealType cos_phi = n1.dot(n2) / (n1_length * n2_length);
            MDSRealType phi = std::acos(cos_phi);
            MDSRealType sign = n1.cross(n2).dot(r_kj) < 0 ? -1 : 1;
            phi *= sign;

            MDSRealType dV_dphi = Dihedral.k_phi_ForceConstant * static_cast<MDSRealType>(Dihedral.N_Multiplicity) * std::sin(static_cast<MDSRealType>(Dihedral.N_Multiplicity) * phi - Dihedral.phi0_PhaseAngleInRadians);

            Vec3 f1 = (n1 * (static_cast<MDSRealType>(Dihedral.N_Multiplicity) * Dihedral.k_phi_ForceConstant / n1_length)).cross(r_kj);
            Vec3 f4 = (n2 * (static_cast<MDSRealType>(Dihedral.N_Multiplicity) * Dihedral.k_phi_ForceConstant / n2_length)).cross(r_kj);
            Vec3 f2 = (f1.cross(r_ij) - f4.cross(r_kl)) / r_kj.dot(r_kj);
            Vec3 f3 = (f2 + f4) * -1.0;

            Atom1.Force += f1 * dV_dphi;
            Atom2.Force += (f1 * -1.0 + f2) * dV_dphi;
            Atom3.Force += (f3 + f4) * -1.0 * dV_dphi;
            Atom4.Force += f4 * dV_dphi;
        }
    }
public:
    void ComputeCoulombForces(Atom& Atom1, Atom& Atom2, const Vec3& delta, const MDSRealType r) const
    {
        const MDSRealType coulomb = K_Coulomb * Atom1.Charge * Atom2.Charge / (r * r * r);
        const Vec3 f_coulomb = delta * coulomb;

        Atom1.Force += f_coulomb;
        Atom2.Force -= f_coulomb;
    }
public:
    static void ComputeLenardJonesInteractionForces(Atom& Atom1, Atom& Atom2, const Vec3& delta, const MDSRealType r)
    {
        const MDSRealType sigma_avg = (Atom1.SigmaValue + Atom2.SigmaValue) / 2;
        const MDSRealType epsilon_avg = sqrt(Atom1.EpsilonValue * Atom2.EpsilonValue);
        const MDSRealType sr6 = pow(sigma_avg / r, 6);
        const MDSRealType lj = 24 * epsilon_avg * (2*sr6*sr6 - sr6) / (r*r);
        const Vec3 f_lj = delta * lj;

        Atom1.Force += f_lj;
        Atom2.Force -= f_lj;
    }
public:
    void ComputeNonBondedForces()
    {
        const UnsignedInt NumberOfAtoms = Atoms.size();
        for (UnsignedInt Atom1Index = 0; Atom1Index < NumberOfAtoms; ++Atom1Index)
            for (UnsignedInt Atom2Index = Atom1Index + 1; Atom2Index < NumberOfAtoms; ++Atom2Index)
            {
                if (BondedPairs.count(Atom1Index * NumberOfAtoms + Atom2Index))
                   continue;

                Atom& Atom1 = Atoms[Atom1Index];
                Atom& Atom2 = Atoms[Atom2Index];
                Vec3 Delta = Atom1.Position - Atom2.Position;
                MDSRealType R_DeltaLength = Delta.length();

                if (R_DeltaLength == 0)
                   continue;

                // Coulomb interaction
                ComputeCoulombForces(Atom1, Atom2, Delta, R_DeltaLength);
                // MDSRealType coulomb = K_Coulomb * a.Charge * b.Charge / (r*r*r);
                // Vec3 f_coulomb = delta * coulomb;
                // a.Force += f_coulomb;
                // b.Force -= f_coulomb;

                // Lennard-Jones interaction
                ComputeLenardJonesInteractionForces(Atom1, Atom2, Delta, R_DeltaLength);
                // MDSRealType sigma_avg = (a.SigmaValue + b.SigmaValue)/2;
                // MDSRealType epsilon_avg = std::sqrt(a.EpsilonValue * b.EpsilonValue);
                // MDSRealType sr6 = std::pow(sigma_avg/r, 6);
                // MDSRealType lj = 24 * epsilon_avg * (2*sr6*sr6 - sr6) / (r*r);
                // Vec3 f_lj = delta * lj;
                // a.Force += f_lj;
                // b.Force -= f_lj;
            }
    }
public:
    Simulation(MDSRealType box, MDSRealType temp, MDSRealType dt) : BoxSize(box), Temperature(temp), TimeStep(dt)
    {
    }
public:
    void initializeAtoms1(UnsignedInt n)
    {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<> pos(0, BoxSize), charge(-1, 1);
        std::uniform_int_distribution<> AtomType(0, 3);

        for (UnsignedInt i = 0; i < n; ++i)
        {
            Atoms.emplace_back(static_cast<AtomTypes>(AtomType(gen)), Vec3(pos(gen), pos(gen), pos(gen)), charge(gen));
            MDSRealType scale = std::sqrt(3 * K_Boltzmann * Temperature / Atoms.back().Mass);
            std::normal_distribution<> vel(0, scale);
            Atoms.back().Velocity = Vec3(vel(gen), vel(gen), vel(gen));
        }

        // Create molecular topology (example: linear chain)
        for (UnsignedInt i = 0; i < n/2; ++i)
        {
            Bonds.emplace_back(i, i+1, 1.5, 100);
            BondedPairs.insert(i*n + i+1);
            if (i < n/2 - 2)
            {
                Angles.emplace_back(i, i+1, i+2, 1.5708, 100);
                Dihedrals.emplace_back(i, i+1, i+2, i+3, 5, 2, 0);
            }
        }
    }

    void initializeAtoms2(int num_atoms)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> pos_dist(0.0, BoxSize);
        std::uniform_int_distribution<> type_dist(0, 3);
        std::uniform_real_distribution<> charge_dist(-1.0, 1.0);

        for (UnsignedInt i = 0; i < num_atoms; ++i)
        {
            Vec3 pos(pos_dist(gen), pos_dist(gen), pos_dist(gen));
            UnsignedInt type = type_dist(gen);
            MDSRealType charge = charge_dist(gen);
            Atoms.emplace_back(static_cast<AtomTypes>(type), pos, charge);
        }

        // Assign random velocities based on Maxwell-Boltzmann distribution
        std::normal_distribution<> vel_dist(0.0, std::sqrt(K_Boltzmann * Temperature));
        for (auto& atom : Atoms)
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
                if (BondedPairs.find(key) == BondedPairs.end())
                {
                    BondedPairs.insert(key);
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
                Angles.emplace_back(i, j, k,
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
                Dihedrals.emplace_back(i, j, k, l,
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
        for (auto& a : Atoms)
            a.Force = Vec3();

        ComputeBondForces();
        ComputeAngleForces();
        ComputeDihedralForces();
        ComputeNonBondedForces();
    }

    void integrate()
    {
        // Velocity half-step
        for (auto& a : Atoms)
        {
            a.Velocity += a.Force * (TimeStep / (2 * a.Mass));
        }

        // Position update with periodic boundary conditions
        for (auto& a : Atoms)
        {
            a.Position += a.Velocity * TimeStep;
            a.Position.x = fmod(a.Position.x + BoxSize, BoxSize);
            a.Position.y = fmod(a.Position.y + BoxSize, BoxSize);
            a.Position.z = fmod(a.Position.z + BoxSize, BoxSize);

            if (a.Position.x < 0) a.Position.x += BoxSize;
            if (a.Position.y < 0) a.Position.y += BoxSize;
            if (a.Position.z < 0) a.Position.z += BoxSize;
        }

        // Compute new forces
        computeForces();

        // Velocity second half-step
        for (auto& a : Atoms)
        {
            a.Velocity += a.Force * (TimeStep / (2 * a.Mass));
        }
    }

    void applyBerendsenThermostat(MDSRealType target_temp, MDSRealType tau)
    {
        MDSRealType current_temp = 0;
        for (const auto& a : Atoms)
        {
            current_temp += a.Mass * a.Velocity.dot(a.Velocity);
        }
        current_temp /= (3 * Atoms.size() * K_Boltzmann);

        MDSRealType lambda = std::sqrt(1 + TimeStep/tau * (target_temp/current_temp - 1));
        for (auto& a : Atoms)
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

