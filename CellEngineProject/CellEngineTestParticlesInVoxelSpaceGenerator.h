
#ifndef CELL_ENGINE_TEST_PARTICLES_GENERATOR_H
#define CELL_ENGINE_TEST_PARTICLES_GENERATOR_H

#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"

class CellEngineTestParticlesInVoxelSpaceGenerator : public CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator
{
public:
    void GeneratePlanedEllipsoidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GeneratePlanedCuboidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    explicit CellEngineTestParticlesInVoxelSpaceGenerator(ParticlesContainer<Particle>& ParticlesParam) : CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator(ParticlesParam)
    {
    }
};

#endif
