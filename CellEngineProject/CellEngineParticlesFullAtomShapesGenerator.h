
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
    typedef bool (CellEngineParticlesFullAtomShapesGenerator::*CheckFreeSpaceForSelectedSpaceType)(float, float, float, float, float, float, float, float, float, UniqueIdInt );
    typedef void (CellEngineParticlesFullAtomShapesGenerator::*SetValueToAtomsForSelectedSpaceType)(ListOfAtomsType*, UniqueIdInt, float, float, float, float, float, float, float, float, float);
protected:
    bool GenerateParticleAtomsWhenSelectedSpaceIsFree(float LocalNewParticleIndex, float PosXStart, float PosYStart, float PosZStart, float SizeOfParticleX, float SizeOfParticleY, float SizeOfParticleZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToAtomsForSelectedSpaceType SetValueToFullAtomsForSelectedSpace);
public:
    bool CheckFreeSpaceInCuboidSelectedSpace(float PosXStart, float PosYStart, float PosZStart, float StepX, float StepY, float StepZ, float SizeOfParticleX, float SizeOfParticleY, float SizeOfParticleZ, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForEllipsoidSelectedSpace(float PosXStart, float PosYStart, float PosZStart, float StepX, float StepY, float StepZ, float RadiusXParam, float RadiusYParam, float RadiusZParam, UniqueIdInt ValueToCheck);
public:
    void SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, UniqueIdInt FullAtomValue, float StartXPosParam, float StartYPosParam, float StartZPosParam, float StepXParam, float StepYParam, float StepZParam, float SizeXParam, float SizeYParam, float SizeZParam);
    void SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, UniqueIdInt FullAtomValue, float PosXStart, float PosYStart, float PosZStart, float StepX, float StepY, float StepZ, float RadiusXParam, float RadiusYParam, float RadiusZParam);
protected:
    explicit CellEngineParticlesFullAtomShapesGenerator(ParticlesContainer<Particle>& ParticlesParam)
    {
    }
    ~CellEngineParticlesFullAtomShapesGenerator() = default;
};

#endif
