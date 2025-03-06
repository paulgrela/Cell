
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_CONSTANTS_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_CONSTANTS_H

#include "../CellEngineTypes.h"

// Constants
const MDSRealType kB = 1.380649e-23; // Boltzmann constant (J/K)
const MDSRealType e_charge = 1.602176634e-19; // Elementary charge (C)
const MDSRealType epsilon0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
const MDSRealType pi = M_PI;

// Parameters for Lennard-Jones potential
const MDSRealType sigma = 3.4e-10; // Ångström to meters
const MDSRealType epsilon = 1.65e-21; // Joules

// Parameters for Hooke's law (bonded interactions)
const MDSRealType k_bond = 500.0; // Bond spring constant (N/m)
const MDSRealType r0 = 1.5e-10; // Equilibrium bond length (m)

// Parameters for angle bending
const MDSRealType k_angle = 100.0; // Angle spring constant (J/rad^2)
const MDSRealType theta0 = pi / 2.0; // Equilibrium angle (radians)

// Parameters for dihedral torsion
const MDSRealType k_dihedral = 10.0; // Dihedral spring constant (J)
const MDSRealType n = 3; // Periodicity of the dihedral potential

// Time step for simulation
const MDSRealType dt = 1e-15; // Time step (1 fs)

#endif
