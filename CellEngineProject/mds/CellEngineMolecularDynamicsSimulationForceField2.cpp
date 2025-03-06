
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <algorithm>

#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceField2.h"
#include "CellEngineMolecularDynamicsSimulationForceField2Types.h"

using namespace std;

class Simulation
{
private:
    std::vector<AtomMDS> Atoms;
    std::vector<BondMDS> Bonds;
    std::vector<AngleMDS> Angles;
    std::vector<DihedralMDS> Dihedrals;
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
            AtomMDS& Atom1 = Atoms[Bond.Atom1Index];
            AtomMDS& Atom2 = Atoms[Bond.Atom2Index];
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
            AtomMDS& Atom1 = Atoms[Angle.Atom1Index];
            AtomMDS& Atom2 = Atoms[Angle.Atom2Index];
            AtomMDS& Atom3 = Atoms[Angle.Atom3Index];

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
            AtomMDS& Atom1 = Atoms[Dihedral.Atom1Index];
            AtomMDS& Atom2 = Atoms[Dihedral.Atom2Index];
            AtomMDS& Atom3 = Atoms[Dihedral.Atom3Index];
            AtomMDS& Atom4 = Atoms[Dihedral.Atom4Index];

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
    void ComputeCoulombForces(AtomMDS& Atom1, AtomMDS& Atom2, const Vec3& delta, const MDSRealType r) const
    {
        const MDSRealType coulomb = K_Coulomb * Atom1.Charge * Atom2.Charge / (r * r * r);
        const Vec3 f_coulomb = delta * coulomb;

        Atom1.Force += f_coulomb;
        Atom2.Force -= f_coulomb;
    }
public:
    static void ComputeLenardJonesInteractionForces(AtomMDS& Atom1, AtomMDS& Atom2, const Vec3& delta, const MDSRealType r)
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

                AtomMDS& Atom1 = Atoms[Atom1Index];
                AtomMDS& Atom2 = Atoms[Atom2Index];
                Vec3 Delta = Atom1.Position - Atom2.Position;
                MDSRealType R_DeltaLength = Delta.length();

                if (R_DeltaLength == 0)
                    continue;

                ComputeCoulombForces(Atom1, Atom2, Delta, R_DeltaLength);

                ComputeLenardJonesInteractionForces(Atom1, Atom2, Delta, R_DeltaLength);
            }
    }
public:
    void ComputeAllForces()
    {
        for (auto& AtomMDS : Atoms)
            AtomMDS.Force = Vec3();

        ComputeBondForces();
        ComputeAngleForces();
        ComputeDihedralForces();
        ComputeNonBondedForces();
    }
public:
    void UpdatePositionWithPeriodicBoundaryConditions()
    {
        for (auto& AtomObject : Atoms)
        {
            AtomObject.Position += AtomObject.Velocity * TimeStep;

            AtomObject.Position.x = fmod(AtomObject.Position.x + BoxSize, BoxSize);
            AtomObject.Position.y = fmod(AtomObject.Position.y + BoxSize, BoxSize);
            AtomObject.Position.z = fmod(AtomObject.Position.z + BoxSize, BoxSize);

            if (AtomObject.Position.x < 0)
                AtomObject.Position.x += BoxSize;
            if (AtomObject.Position.y < 0)
                AtomObject.Position.y += BoxSize;
            if (AtomObject.Position.z < 0)
                AtomObject.Position.z += BoxSize;
        }
    }
public:
    void IntegrateAllComputations()
    {
        for (auto& AtomObject : Atoms)
            AtomObject.Velocity += AtomObject.Force * (TimeStep / (2 * AtomObject.Mass));

        UpdatePositionWithPeriodicBoundaryConditions();

        ComputeAllForces();

        for (auto& AtomObject : Atoms)
            AtomObject.Velocity += AtomObject.Force * (TimeStep / (2 * AtomObject.Mass));
    }
public:
    void ApplyBerendsenThermostat(const MDSRealType TargetTemperature, const MDSRealType Tau)
    {
        MDSRealType CurrentTemperature = 0;

        for (const auto& AtomObject : Atoms)
            CurrentTemperature += AtomObject.Mass * AtomObject.Velocity.dot(AtomObject.Velocity);

        CurrentTemperature /= (3 * Atoms.size() * K_Boltzmann);

        MDSRealType Lambda = std::sqrt(1 + TimeStep / Tau * (TargetTemperature / CurrentTemperature - 1));
        for (auto& AtomObject : Atoms)
            AtomObject.Velocity = AtomObject.Velocity * Lambda;
    }
public:
    Simulation(const MDSRealType BoxSizeParam, const MDSRealType TemperaturParam, const MDSRealType TimeStepParam) : BoxSize(BoxSizeParam), Temperature(TemperaturParam), TimeStep(TimeStepParam)
    {
    }
public:
    void InitializeAtoms1(const UnsignedInt N)
    {
        std::mt19937 RandomGenerator(std::random_device{}());
        std::uniform_real_distribution<> PositionDistribution(0, BoxSize);
        std::uniform_real_distribution<> ChargeDistribution(-1, 1);
        std::uniform_int_distribution<> AtomTypeDistribution(0, 3);

        for (UnsignedInt Index = 0; Index < N; ++Index)
        {
            Atoms.emplace_back(static_cast<AtomTypes>(AtomTypeDistribution(RandomGenerator)), Vec3(PositionDistribution(RandomGenerator), PositionDistribution(RandomGenerator), PositionDistribution(RandomGenerator)), ChargeDistribution(RandomGenerator));
            MDSRealType Scale = sqrt(3 * K_Boltzmann * Temperature / Atoms.back().Mass);
            std::normal_distribution<> VelocityDistribution(0, Scale);
            Atoms.back().Velocity = Vec3(VelocityDistribution(RandomGenerator), VelocityDistribution(RandomGenerator), VelocityDistribution(RandomGenerator));
        }

        // Create molecular topology (example: linear chain)
        for (UnsignedInt Index = 0; Index < N / 2; ++Index)
        {
            Bonds.emplace_back(Index, Index + 1, 1.5, 100);
            BondedPairs.insert(Index * N + Index + 1);
            if (Index < N / 2 - 2)
            {
                Angles.emplace_back(Index, Index + 1, Index + 2, 1.5708, 100);
                Dihedrals.emplace_back(Index, Index + 1, Index + 2, Index + 3, 5, 2, 0);
            }
        }
    }
private:
    void GenerateRandomAtoms(mt19937& RandomGenerator, const UnsignedInt NumberOfAtoms)
    {
        std::uniform_real_distribution<> PositionsDistribution(0.0, BoxSize);
        std::uniform_int_distribution<> AtomTypeDistribution(0, 3);
        std::uniform_real_distribution<> ChargeDistribution(-1.0, 1.0);

        for (UnsignedInt Index = 0; Index < NumberOfAtoms; ++Index)
        {
            Vec3 pos(PositionsDistribution(RandomGenerator), PositionsDistribution(RandomGenerator), PositionsDistribution(RandomGenerator));
            UnsignedInt type = AtomTypeDistribution(RandomGenerator);
            MDSRealType charge = ChargeDistribution(RandomGenerator);
            Atoms.emplace_back(static_cast<AtomTypes>(type), pos, charge);
        }
    }
    void GenerateAssignRandomVelocitiesBasedOnMaxwellBoltzmannDistribution(mt19937& RandomGenerator)
    {
        std::normal_distribution<> VelocityDistribution(0.0, sqrt(K_Boltzmann * Temperature));
        for (auto& AtomObject : Atoms)
        {
            AtomObject.Velocity.x = VelocityDistribution(RandomGenerator);
            AtomObject.Velocity.y = VelocityDistribution(RandomGenerator);
            AtomObject.Velocity.z = VelocityDistribution(RandomGenerator);
            // Subtract center-of-mass velocity to avoid drift
            // (omitted for brevity)
        }
    }
    void GenerateRandomBonds(uniform_int_distribution<>& AtomsIndexesDistribution, mt19937& RandomGenerator, const UnsignedInt NumberOfAtoms)
    {
        for (UnsignedInt Bond = 0; Bond < 500; ++Bond)
        {
            UnsignedInt Atom1Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom2Index = AtomsIndexesDistribution(RandomGenerator);
            if (Atom1Index != Atom2Index)
            {
                // Ensure i < j to avoid duplicates
                if (Atom1Index > Atom2Index)
                    swap(Atom1Index, Atom2Index);
                // Check if bond already exists
                UnsignedInt key = Atom1Index * NumberOfAtoms + Atom2Index;
                if (BondedPairs.find(key) == BondedPairs.end())
                {
                    BondedPairs.insert(key);
                    // Assign random r0 and k (example values)
                    MDSRealType r0 = 1.5; // Angstrom
                    MDSRealType k = 100.0; // kcal/(mol·Å²)
                    Bonds.emplace_back(Atom1Index, Atom2Index, r0, k);
                }
            }
        }
    }
    void GenerateRandomAngles(uniform_int_distribution<>& AtomsIndexesDistribution, mt19937& RandomGenerator)
    {
        for (UnsignedInt Angle = 0; Angle < 300; ++Angle)
        {
            UnsignedInt Atom1Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom2Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom3Index = AtomsIndexesDistribution(RandomGenerator);
            if (Atom1Index != Atom2Index && Atom2Index != Atom3Index && Atom1Index != Atom3Index)
                Angles.emplace_back(Atom1Index, Atom2Index, Atom3Index, 1.5708, 100.0);
        }
    }
    void GenerateRandomDihedrals(uniform_int_distribution<>& AtomsIndexesDistribution, mt19937& RandomGenerator)
    {
        for (UnsignedInt Dihedral = 0; Dihedral < 200; ++Dihedral)
        {
            UnsignedInt Atom1Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom2Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom3Index = AtomsIndexesDistribution(RandomGenerator);
            UnsignedInt Atom4Index = AtomsIndexesDistribution(RandomGenerator);
            if (Atom1Index != Atom2Index && Atom2Index != Atom3Index && Atom3Index != Atom4Index && Atom1Index != Atom3Index && Atom1Index != Atom4Index && Atom2Index != Atom4Index)
                Dihedrals.emplace_back(Atom1Index, Atom2Index, Atom3Index, Atom4Index, 5.0, 2, 0.0);
        }
    }
public:
    void InitializeAtoms2(const UnsignedInt NumberOfAtoms)
    {
        std::random_device RandomDevice;
        std::mt19937 RandomGenerator(RandomDevice());

        GenerateRandomAtoms(RandomGenerator, NumberOfAtoms);
        GenerateAssignRandomVelocitiesBasedOnMaxwellBoltzmannDistribution(RandomGenerator);

        std::uniform_int_distribution<> AtomsIndexesDistribution(0, static_cast<int>(NumberOfAtoms) - 1);
        GenerateRandomBonds(AtomsIndexesDistribution, RandomGenerator, NumberOfAtoms);
        // Generate angles and dihedrals (example)
        // For angles, find triplets i-j-k where i-j and j-k are bonded
        // This is a simplified approach; realistically, need to find connected triplets
        // For demonstration, randomly select triplets
        GenerateRandomAngles(AtomsIndexesDistribution, RandomGenerator);
        GenerateRandomDihedrals(AtomsIndexesDistribution, RandomGenerator);
    }
};

UnsignedInt ComputeMolecularDynamicsSimulationForceField2()
{
    constexpr MDSRealType BoxSizeInAngstroms = 100.0;
    constexpr MDSRealType TemperatureInKelvins = 303.15;
    constexpr MDSRealType TimeStep = 0.001;    // picoseconds
    constexpr UnsignedInt NumberOfSteps = 1000;
    constexpr UnsignedInt OutputInterval = 100;

    Simulation SimulationObject(BoxSizeInAngstroms, TemperatureInKelvins, TimeStep);
    SimulationObject.InitializeAtoms1(1000);
    SimulationObject.ComputeAllForces();

    for (UnsignedInt Step = 0; Step < NumberOfSteps; ++Step)
    {
        SimulationObject.IntegrateAllComputations();
        SimulationObject.ApplyBerendsenThermostat(TemperatureInKelvins, 0.1);

        if (Step % OutputInterval == 0)
            cout << "Step " << Step << " completed" << endl;
    }

    return 0;
}
