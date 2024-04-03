
#ifndef CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H
#define CELL_ENGINE_REAL_RANDOM_PARTICLES_GENERATOR_H

#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesVoxelsShapesGenerator.h"

class CellEngineRealRandomParticlesGenerator : public CellEngineParticlesVoxelsShapesGenerator,  virtual public CellEngineRandomDeviceEngine
{
public:
    void GenerateAllRealRandomParticles();
    void GenerateRealRandomMembraneParticles();
    void GenerateRealRandomRibosomesParticles();
protected:
    explicit CellEngineRealRandomParticlesGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineParticlesVoxelsShapesGenerator(ParticlesParam)
    {
    }
};

#endif
