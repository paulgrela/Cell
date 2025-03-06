
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_TYPES_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_TYPES_H

// Structure to represent an atom
struct Atom
{
    MDSRealType x, y, z; // Position (m)
    MDSRealType vx, vy, vz; // Velocity (m/s)
    MDSRealType fx, fy, fz; // Force (N)
    MDSRealType Charge; // Electric charge (C)
    MDSRealType Mass; // Mass (kg)
};

// Structure to represent a bond between two atoms
struct Bond
{
    int Atom1, Atom2; // Indices of bonded atoms
    MDSRealType k_bond, r0;
};

// Structure to represent an angle between three atoms
struct Angle
{
    int Atom1, Atom2, Atom3; // Indices of atoms forming the angle
};

// Structure to represent a dihedral torsion between four atoms
struct Dihedral
{
    int Atom1, Atom2, Atom3, Atom4; // Indices of atoms forming the dihedral
};

#endif
