
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_CONSTANTS_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_CONSTANTS_H

#include "../CellEngineTypes.h"

constexpr MDSRealType kB_BoltzmanConstantJPerK = 1.380649e-23;
constexpr MDSRealType ElementaryCharge = 1.602176634e-19;
constexpr MDSRealType epsilon0_VacuumPermittivityFPerM = 8.8541878128e-12;
constexpr MDSRealType pi = M_PI;

constexpr MDSRealType sigma_LennardJonesPotentialDistanceParameterAngstromPerM = 3.4e-10;
constexpr MDSRealType epsilon_LennardJonesPotentialDepthJouls = 1.65e-21;

constexpr MDSRealType k_bond_HookesLawBondedInteractionSpringConstantNPerM = 500.0;
constexpr MDSRealType r0_EquilibriumBondLengthInM = 1.5e-10;

constexpr MDSRealType k_angle_AngleBendingSpringConstant = 100.0;
constexpr MDSRealType theta0_EquilibriumAngle90Degree = pi / 2.0;

constexpr MDSRealType k_dihedral_DiherdralTorsionSpringConstantInJoules = 10.0;
constexpr MDSRealType PeriodicityOfDihedralTorsionPotential = 3;

constexpr MDSRealType dt_TimeStepForSimulation1Fs = 1e-15;

#endif
