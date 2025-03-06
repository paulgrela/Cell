
#ifndef CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_TYPES_H
#define CELL_ENGINE_MOLECULAR_DYNAMICS_SIMULATION_FORCE_FIELD_4_TYPES_H

struct AtomMDS
{
    MDSRealType PositionX, PositionY, PositionZ;
    MDSRealType VelocityX, VelocityY, VelocityZ;
    MDSRealType ForceX, ForceY, ForceZ;
    MDSRealType Charge;
    MDSRealType Mass;
};

struct BondMDS
{
    UnsignedInt Atom1Index;
    UnsignedInt Atom2Index;
    MDSRealType k_bond, r0;
};

struct AngleMDS
{
    int Atom1Index, Atom2Index, Atom3Index;
};

struct DihedralMDS
{
    int Atom1Index, Atom2Index, Atom3Index, Atom4Index;
};

#endif
