
#ifndef CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_SPACE_GENERATOR_H
#define CELL_ENGINE_GENOME_NUCLEIC_ACIDS_PARTICLES_IN_SPACE_GENERATOR_H

#include <unordered_set>

#include "CellEngineRealRandomParticlesInFullAtomSpaceGenerator.h"

class CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator : virtual public CellEngineBasicParticlesOperations
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
    void TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms() const;
    static void FindInterGenesSequences();
public:
    virtual void GenerateParticleVoxelsWhenSelectedSpaceIsFreeBasic(UnsignedInt ParticleIndex, UnsignedInt StartPosX, UnsignedInt StartPosY, UnsignedInt StartPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleStepX, UnsignedInt ParticleStepY, UnsignedInt ParticleStepZ)
    {
    }
public:
    void GetMinMaxCoordinatesForDNA();
public:
    void SaveGenomeDataToFile(UnsignedInt ParticleSize);
    void ReadGenomeDataFromFile(bool PairedBool);
    void ReadGenomeSequenceFromFile(bool Paired);
    void ReadGenomeSequenceFromFastaFile();
    void TestGeneratedGenomeCorrectness(UnsignedInt ParticleSize);
protected:
    ~CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator() override = default;
};

#endif
