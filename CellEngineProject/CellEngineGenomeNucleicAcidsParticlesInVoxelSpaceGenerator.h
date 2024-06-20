
#ifndef CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_VOXEL_SPACE_GENERATOR_H
#define CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_VOXEL_SPACE_GENERATOR_H

#include <unordered_set>

#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"

class CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator : public CellEngineRealRandomParticlesInVoxelSpaceGenerator
{
protected:
    std::vector<std::string> GenomesLines;
    std::vector<std::vector<UniqueIdInt>> Genomes;
protected:
    void GenerateOneStrand(EntityIdInt EntityId, std::string_view Sequence, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleStepX, UnsignedInt ParticleStepY, UnsignedInt ParticleStepZ);
    Particle* GenerateNucleotideParticle(Particle* ParticlePrev, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeThread, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, bool AddToGenome, std::vector<UniqueIdInt>& Genome, vector3_16 UniqueColorParam, bool LinkWithPreviousNucleotide);
    std::tuple<Particle*, Particle*> GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, EntityIdInt EntityId, ChainIdInt ChainId, UnsignedInt GenomeIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt AddSizeX, UnsignedInt AddSizeY, UnsignedInt AddSizeZ, vector3_16 UniqueColorParam, bool Linked, bool LinkWithPreviousNucleotide);
public:
    void EraseAllDNAParticles();
public:
    void GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5);
protected:
    static void UpdateRandomPositions(UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, UnsignedInt Size);
    static bool TestFormerForbiddenPositions(std::unordered_set<std::string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt Size);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> EraseLastRandomDNAParticle(std::vector<UniqueIdInt>& Genome);
protected:
    void GetMinMaxCoordinatesForDNA();
public:
    void SaveGenomeDataToFile(UnsignedInt ParticleSize);
    void ReadGenomeDataFromFile(bool Paired);
    void ReadGenomeSequenceFromFile();
    void TestGeneratedGenomeCorrectness(UnsignedInt ParticleSize);
protected:
    explicit CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineRealRandomParticlesInVoxelSpaceGenerator(ParticlesParam)
    {
    }
};

#endif
