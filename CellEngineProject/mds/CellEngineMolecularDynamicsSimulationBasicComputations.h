
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_BASIC_FORCES_COMPUTATIONS_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_BASIC_FORCES_COMPUTATIONS_H

#include "../CellEngineTypes.h"

MDSRealType CountLeonardJonesPotential(MDSRealType r);
MDSRealType CountLeonardJoneslForce(MDSRealType r);
MDSRealType CountCoulombPotential(MDSRealType q1, MDSRealType q2, MDSRealType r);
MDSRealType CountCoulombForce(MDSRealType q1, MDSRealType q2, MDSRealType r);
MDSRealType CountAngleBendingPotential(MDSRealType theta, MDSRealType theta0);
MDSRealType CountAngleBendingForce(MDSRealType theta, MDSRealType theta0);
MDSRealType CountDihedralTorsionPotential(MDSRealType phi);
MDSRealType CountDihedralTorsionForce(MDSRealType phi);
MDSRealType CountBondHookesLawPotential(MDSRealType r, MDSRealType r0);
MDSRealType CountHookesLawBondForce(MDSRealType r, MDSRealType r0);
MDSRealType CountDistanceBetweenAtoms(const AtomMDS& Atom1, const AtomMDS& Atom2);
MDSRealType CalculateAngleBetweenThreeAtoms(const AtomMDS& Atom1, const AtomMDS& Atom2, const AtomMDS& Atom3);

#endif
