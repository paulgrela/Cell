
#ifndef CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_VOXEL_SPACE_GENERATOR_H
#define CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_VOXEL_SPACE_GENERATOR_H

#include <unordered_set>

#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"
#include "CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator.h"

class CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator : public CellEngineRealRandomParticlesInVoxelSpaceGenerator, public CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator
{
public:
    void AddPairedNucleotidesToExistingNucleotidesFromGenome(UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ);
public:
    void GenerateRandomDNAInWholeCell1or2RandomTurn(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5, bool Paired);
    void GenerateRandomDNAInWholeCell1or2Vertical(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5, bool Paired);
protected:
    static void UpdateRandomPositions(UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, UnsignedInt Size);
    static bool TestFormerForbiddenPositions(const MainSetType<std::string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt Size);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> EraseLastRandomDNAParticle(std::vector<UniqueIdInt>& Genome);
protected:
    void GenerateParticleElementsWhenSelectedSpaceIsFreeBasic(UnsignedInt ParticleIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleStepX, UnsignedInt ParticleStepY, UnsignedInt ParticleStepZ) override;
protected:
    explicit CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator(ParticlesContainer<Particle>& ParticlesParam) : CellEngineRealRandomParticlesInVoxelSpaceGenerator(ParticlesParam), CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator()
    {
    }
    ~CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator() override = default;
};

#endif
