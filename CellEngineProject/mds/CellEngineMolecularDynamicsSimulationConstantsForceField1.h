
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_CONSTANTS_FORCE_FIELD_1_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_CONSTANTS_FORCE_FIELD_1_H

#include "../CellEngineTypes.h"

// Constants for the physical simulation
const MDSRealType kB = 1.380649e-23; // Boltzmann constant (J/K)
const MDSRealType T = 303.15; // Temperature in Kelvin (30°C)
const MDSRealType dt = 1e-15; // Time step in seconds (1 fs)
const MDSRealType e0 = 8.854187817e-12; // Permittivity of free space
const MDSRealType charge_constant = 8.9875517923e9; // Coulomb constant (1/(4*pi*e0))
const MDSRealType mass_O = 2.656e-26; // Mass of an oxygen atom (kg)
const MDSRealType mass_N = 2.326e-26; // Mass of a nitrogen atom (kg)
const MDSRealType mass_P = 5.14e-26; // Mass of a phosphorus atom (kg)
const MDSRealType mass_C = 1.994e-26; // Mass of a carbon atom (kg)
const MDSRealType bond_strength_k = 100.0; // Bond strength constant (N/m)
const MDSRealType bond_equilibrium_distance = 1.0e-10; // Equilibrium bond distance (1 Å)

const MDSRealType r_cutoff = 2.5e-10; // Lennard-Jones cutoff distance
const MDSRealType epsilon = 1.0; // Lennard-Jones potential depth
const MDSRealType sigma = 1.0e-10; // Lennard-Jones distance parameter


// Constants for Lennard-Jones potential
//const MDSRealType epsilon = 1.0; // Depth of the potential well
//const MDSRealType sigma = 1.0;   // Distance at which potential is zero
//const MDSRealType r_cutoff = 2.5 * sigma; // Cutoff distance for Lennard-Jones interactions

// Constants for Hooke's Law (bond stretching)
const MDSRealType k_bond = 100.0; // Bond force constant
const MDSRealType r0_bond = 1.0;  // Equilibrium bond length

// Constants for angle bending
const MDSRealType k_theta = 50.0;  // Force constant for angle bending
const MDSRealType theta0 = M_PI / 2; // Equilibrium angle (90 degrees)

// Constants for dihedral torsion
const MDSRealType V_n = 1.0;
const MDSRealType n_period = 3.0;
//const MDSRealType gamma = 0.0;

// Coulomb constant
const MDSRealType k_e = 8.99e9;  // Coulomb's constant

// Time step and temperature constants
// const MDSRealType dt = 0.001;  // Time step for the simulation
const MDSRealType target_temp = 300.0;  // Target temperature for the thermostat

#endif
