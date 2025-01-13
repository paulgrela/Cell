
#ifndef CELL_ENGINE_PARTICLES_VOXELS_SHAPES_GENERATOR_H
#define CELL_ENGINE_PARTICLES_VOXELS_SHAPES_GENERATOR_H

#include "CellEngineBasicParticlesOperations.h"

class CellEngineParticlesVoxelsShapesGenerator : virtual public CellEngineParticlesVoxelsOperations
{
protected:
    virtual Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) = 0;
    virtual void ClearVoxelSpaceAndParticles() = 0;
protected:
    typedef bool (CellEngineParticlesVoxelsShapesGenerator::*CheckFreeSpaceForSelectedSpaceType)(UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UniqueIdInt );
    typedef void (CellEngineParticlesVoxelsShapesGenerator::*SetValueToVoxelsForSelectedSpaceType)(ListOfElements*, UniqueIdInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt);
protected:
    bool GenerateParticleVoxelsWhenSelectedSpaceIsFree(UnsignedInt LocalNewParticleIndex, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToVoxelsForSelectedSpaceType SetValueToVoxelsForSelectedSpace);
protected:
    void SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(ListOfElements* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ);
public:
    bool CheckFreeSpaceInCuboidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForSphereSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForEllipsoidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
public:
    void SetValueToVoxelsForCuboidSelectedSpace(ListOfElements* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void SetValueToVoxelsForSphereSelectedSpace(ListOfElements* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
    void SetValueToVoxelsForEllipsoidSelectedSpace(ListOfElements* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
protected:
    explicit CellEngineParticlesVoxelsShapesGenerator(ParticlesContainer<Particle>& ParticlesParam)
    {
    }
    ~CellEngineParticlesVoxelsShapesGenerator() = default;
};

#endif
