
#include "../CellEngineTypes.h"

#include "CellEngineMolecularDynamicsSimulationForceFieldCommonTypes.h"
#include "CellEngineMolecularDynamicsSimulationForceField4Constants.h"
#include "CellEngineMolecularDynamicsSimulationBasicComputations.h"

MDSRealType CountLeonardJonesPotential(const MDSRealType r)
{
    const MDSRealType sr6 = pow(sigma_LennardJonesPotentialDistanceParameterAngstromPerM / r, 6);
    const MDSRealType sr12 = sr6 * sr6;

    return 4 * epsilon_LennardJonesPotentialDepthJouls * (sr12 - sr6);
}

MDSRealType CountLeonardJoneslForce(const MDSRealType r)
{
    const MDSRealType sr6 = pow(sigma_LennardJonesPotentialDistanceParameterAngstromPerM / r, 6);
    const MDSRealType sr12 = sr6 * sr6;

    return 24 * epsilon_LennardJonesPotentialDepthJouls * (2 * sr12 - sr6) / r;
}

MDSRealType CountCoulombPotential(const MDSRealType q1, const MDSRealType q2, const MDSRealType r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0_VacuumPermittivityFPerM * r);
}

MDSRealType CountCoulombForce(const MDSRealType q1, const MDSRealType q2, const MDSRealType r)
{
    return (q1 * q2) / (4 * M_PI * epsilon0_VacuumPermittivityFPerM * r * r);
}

MDSRealType CountAngleBendingPotential(const MDSRealType theta, const MDSRealType theta0)
{
    return 0.5 * k_angle_AngleBendingSpringConstant * pow(theta - theta0, 2);
}

MDSRealType CountAngleBendingForce(const MDSRealType theta, const MDSRealType theta0)
{
    return k_angle_AngleBendingSpringConstant * (theta - theta0);
}

MDSRealType CountDihedralTorsionPotential(const MDSRealType phi)
{
    return k_dihedral_DiherdralTorsionSpringConstantInJoules * (1 + cos(3 * phi - M_PI));
}

MDSRealType CountDihedralTorsionForce(const MDSRealType phi)
{
    return -3 * k_dihedral_DiherdralTorsionSpringConstantInJoules * sin(3 * phi - M_PI);
}

MDSRealType CountBondHookesLawPotential(const MDSRealType r, const MDSRealType r0)
{
    return 0.5 * k_bond_HookesLawBondedInteractionSpringConstantNPerM * pow(r - r0, 2);
}

MDSRealType CountHookesLawBondForce(const MDSRealType r, const MDSRealType r0)
{
    return k_bond_HookesLawBondedInteractionSpringConstantNPerM * (r - r0);
}

MDSRealType CountDistanceBetweenAtoms(const AtomMDS& Atom1, const AtomMDS& Atom2)
{
    const MDSRealType dx = Atom1.PositionX - Atom2.PositionX;
    const MDSRealType dy = Atom1.PositionY - Atom2.PositionY;
    const MDSRealType dz = Atom1.PositionZ - Atom2.PositionZ;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

MDSRealType CalculateAngleBetweenThreeAtoms(const AtomMDS& Atom1, const AtomMDS& Atom2, const AtomMDS& Atom3)
{
    const MDSRealType r12 = CountDistanceBetweenAtoms(Atom1, Atom2);
    const MDSRealType r23 = CountDistanceBetweenAtoms(Atom2, Atom3);
    const MDSRealType r13 = CountDistanceBetweenAtoms(Atom1, Atom3);

    return acos((r12 * r12 + r23 * r23 - r13 * r13) / (2 * r12 * r23));
}
