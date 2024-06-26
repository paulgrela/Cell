
#ifndef CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H
#define CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H

#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"

class CellEngineRealRandomParticlesInVoxelSpaceGenerator : public CellEngineParticlesVoxelsShapesGenerator, virtual public CellEngineRandomDeviceEngine, virtual public CellEngineBasicParticlesOperations
{
private:
    const UnsignedInt MaxNumberOfTriesToInsertNewParticle = 1000;
    const UnsignedInt MaxNumberOfTriesToFindARandomPositionInCell = 1000;
private:
    UnsignedInt TotalNumberOfAllParticles = 0;
public:
    void GenerateAllRealRandomParticles();
protected:
    std::tuple<UnsignedInt, UnsignedInt> GetNumberOfParticlesKind(ParticlesTypes ParticleTypeParam, bool AddToTotalNumberOfAllParticles);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetRandomPositionInsideSphere(UnsignedInt Radius, UnsignedInt RadiusSize);
public:
    void InsertNewRandomParticlesForType(ParticlesTypes ParticleTypeParam, UnsignedInt Radius, UnsignedInt RadiusSize);
    void TryToGenerateRandomParticlesForType(const std::pair<const EntityIdInt, ParticleKind>& ParticleKindObject, UnsignedInt Radius, UnsignedInt RadiusSize, UnsignedInt& NumberOfErrors, GeneIdInt GeneNumId, const std::string& GeneSequence, UnsignedInt GeneSequenceLength);
public:
    void PrintNumberOfParticlesForAllMainTypesOfParticles();
    static void PrintInformationAboutRibosomesProteins();
protected:
    explicit CellEngineRealRandomParticlesInVoxelSpaceGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineParticlesVoxelsShapesGenerator(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
