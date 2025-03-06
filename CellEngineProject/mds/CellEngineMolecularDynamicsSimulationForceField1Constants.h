
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_1_CONSTANTS_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_1_CONSTANTS_H

#include "../CellEngineTypes.h"

constexpr MDSRealType kB_BoltzmanConstantJPerK = 1.380649e-23;
constexpr MDSRealType TemperatureInKelvins = 303.15;
constexpr MDSRealType TargetTempratureForTheThermostat = 300.0;
constexpr MDSRealType dt_TimeStepFirstInSeconds = 1e-15;
constexpr MDSRealType dt_TimeStepSecond = 0.001;
constexpr MDSRealType e0 = 8.854187817e-12;
constexpr MDSRealType ChargeCoulombConstant_1Per4MulPIMule0 = 8.9875517923e9;
constexpr MDSRealType k_e_CoulombConstant = 8.99e9;
constexpr MDSRealType mass_O_kg = 2.656e-26;
constexpr MDSRealType mass_N_kg = 2.326e-26;
constexpr MDSRealType mass_P_kg = 5.14e-26;
constexpr MDSRealType mass_C_kg = 1.994e-26;
constexpr MDSRealType BondStrengthConstantK_NPerM = 100.0;
constexpr MDSRealType BondEquilibriumDistance = 1.0e-10;

constexpr MDSRealType epsilon_LennardJonesPotentialDepth = 1.0;
constexpr MDSRealType sigma_LennardJonesDistanceParameter = 1.0e-10;
constexpr MDSRealType r_cutoff_LennardJonesCutoffDistance = 2.5e-10;

constexpr MDSRealType epsilon_LennardJonesPotentialDepthWell = 1.0;
constexpr MDSRealType sigma_LennardJonesDistanceWherePotenitialIsZeroParameter = 1.0;
constexpr MDSRealType r_cutoff_LennardJonesCutoffDistance1 = 2.5 * sigma_LennardJonesDistanceWherePotenitialIsZeroParameter;

constexpr MDSRealType k_bond_HookesLawBondStretchingForceConstant = 100.0;
constexpr MDSRealType r0_bond_HookesLawBondStretchingEquilibriumBondLength = 1.0;

constexpr MDSRealType k_theta_AngleBendingForceConstant = 50.0;
constexpr MDSRealType theta0_EquilibriumAngle90Degree = M_PI / 2;

constexpr MDSRealType V_n_ForDihedralTorsionConstant = 1.0;
constexpr MDSRealType n_period_ForDihedralTorsionConstant = 3.0;
constexpr MDSRealType Gamma_ForDihedralTorsionConstant = 0.0;

#endif
