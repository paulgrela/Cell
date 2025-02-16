
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_SHAPES_GENERATOR_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_SHAPES_GENERATOR_H

#include "CellEngineParticlesFullAtomOperations.h"
#include "CellEngineBasicParticlesOperations.h"

class CellEngineParticlesFullAtomShapesGenerator : virtual public CellEngineParticlesFullAtomOperations
{
protected:
    virtual Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) = 0;
    virtual void ClearFullAtomSpaceAndParticles() = 0;
protected:
    typedef bool (CellEngineParticlesFullAtomShapesGenerator::*CheckFreeSpaceForSelectedSpaceType)(const ParticlesContainer<Particle>& ParticlesParam, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, UniqueIdInt);
    typedef void (CellEngineParticlesFullAtomShapesGenerator::*SetValueToAtomsForSelectedSpaceType)(ListOfAtomsType&, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType, RealType);
protected:
    bool GenerateParticleAtomsWhenSelectedSpaceIsFree(const ParticlesContainer<Particle>& ParticlesParam, UnsignedInt LocalNewParticleIndex, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType SizeOfParticleX, RealType SizeOfParticleY, RealType SizeOfParticleZ, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToAtomsForSelectedSpaceType SetValueToFullAtomsForSelectedSpace);
public:
    bool CheckFreeSpaceInCuboidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType SizeOfParticleX, RealType SizeOfParticleY, RealType SizeOfParticleZ, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForEllipsoidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType RadiusXParam, RealType RadiusYParam, RealType RadiusZParam, UniqueIdInt ValueToCheck);
public:
    void SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, RealType StartXPosParam, RealType StartYPosParam, RealType StartZPosParam, RealType StepXParam, RealType StepYParam, RealType StepZParam, RealType SizeXParam, RealType SizeYParam, RealType SizeZParam);
    void SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, RealType PosXStart, RealType PosYStart, RealType PosZStart, RealType StepX, RealType StepY, RealType StepZ, RealType RadiusXParam, RealType RadiusYParam, RealType RadiusZParam);
protected:
    explicit CellEngineParticlesFullAtomShapesGenerator(ParticlesContainer<Particle>& ParticlesParam)
    {
    }
    ~CellEngineParticlesFullAtomShapesGenerator() = default;
};

#endif
