
#ifndef CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H
#define CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H

#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"

class CellEngineRealRandomParticlesInVoxelSpaceGenerator : public CellEngineParticlesVoxelsShapesGenerator, virtual public CellEngineRandomDeviceEngine, virtual public CellEngineBasicParticlesOperations
{
private:
    const UnsignedInt MaxNumberOfTriesToInsertNewParticle = 100;
    const UnsignedInt MaxNumberOfTriesToFindARandomPositionInCell = 1000;
private:
    UnsignedInt TotalNumberOfAllParticles = 0;
public:
    void GenerateAllRealRandomParticles();
public:
    void GenerateRealRandomOtherParticles();
    void GenerateRealRandomBasicParticles();
    void GenerateRealRandomLipidParticles();
    void GenerateRealRandom_tRNAParticles();
    void GenerateRealRandom_mRNAParticles();
    void GenerateRealRandom_rRNAParticles();
    void GenerateRealRandomMembraneParticles();
    void GenerateRealRandomRibosomesParticles();
    void GenerateRealRandomPolymeraseParticles();
    void GenerateRealRandomRNAPolymeraseParticles();
protected:
    std::tuple<UnsignedInt, UnsignedInt> GetNumberOfParticlesKind(ParticlesTypes ParticleTypeParam);
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetRandomPositionInsideSphere();
protected:
    void PrintNumberOfParticlesForAllMainTypesOfParticles();
protected:
    explicit CellEngineRealRandomParticlesInVoxelSpaceGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineParticlesVoxelsShapesGenerator(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
