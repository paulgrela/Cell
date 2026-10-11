
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_SHAPES_GENERATOR_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_SHAPES_GENERATOR_H

#include "CellEngineParticlesFullAtomOperations.h"
#include "CellEngineBasicParticlesOperations.h"

class CellEngineParticlesFullAtomShapesGenerator : virtual public CellEngineParticlesFullAtomOperations
{
protected:
    virtual Particle& GetParticleFromIndexForGenerator(UniqueIdUInt64Type ParticleIndex) = 0;
    virtual void ClearFullAtomSpaceAndParticles() = 0;
protected:
    typedef bool (*CheckFreeSpaceForSelectedSpaceType)(const ParticlesContainer<Particle>& ParticlesParam, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, UniqueIdUInt64Type);
    typedef void (*SetValueToAtomsForSelectedSpaceType)(ListOfAtomsType&, UniqueIdUInt32Type, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType);
protected:
    bool GenerateParticleAtomsWhenSelectedSpaceIsFree(const ParticlesContainer<Particle>& ParticlesParam, UnsignedInt LocalNewParticleIndex, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType SizeOfParticleX, RealType SizeOfParticleY, RealType SizeOfParticleZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToAtomsForSelectedSpaceType SetValueToFullAtomsForSelectedSpace);
public:
    static bool CheckFreeSpaceInCuboidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType SizeOfParticleX, RealType SizeOfParticleY, RealType SizeOfParticleZ, UniqueIdUInt64Type ValueToCheck);
    static bool CheckFreeSpaceForEllipsoidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType RadiusXParam, RealType RadiusYParam, RealType RadiusZParam, UniqueIdUInt64Type ValueToCheck);
public:
    static void SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, UniqueIdUInt32Type ParticleKindId, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType StepXParam, RealType StepYParam, RealType StepZParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam);
    static void SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, UniqueIdUInt32Type ParticleKindId, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType RadiusXParam, RealType RadiusYParam, RealType RadiusZParam);
protected:
    explicit CellEngineParticlesFullAtomShapesGenerator(ParticlesContainer<Particle>& ParticlesParam)
    {
    }
    ~CellEngineParticlesFullAtomShapesGenerator() = default;
};

#endif
