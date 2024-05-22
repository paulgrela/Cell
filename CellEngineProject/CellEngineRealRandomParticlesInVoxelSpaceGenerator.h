
#ifndef CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H
#define CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H

#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"

class CellEngineRealRandomParticlesInVoxelSpaceGenerator : public CellEngineParticlesVoxelsShapesGenerator, virtual public CellEngineRandomDeviceEngine, virtual public CellEngineBasicParticlesOperations
{
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
    explicit CellEngineRealRandomParticlesInVoxelSpaceGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineParticlesVoxelsShapesGenerator(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

#endif
