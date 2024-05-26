
#ifndef CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H
#define CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H

#include "CellEngineBasicParticlesOperations.h"

class CellEngineSimulationSpaceStatistics : virtual public CellEngineBasicParticlesOperations
{
protected:
    std::vector<std::vector<Particle>> ParticlesSnapshots;
    std::vector<std::unordered_map<UniqueIdInt, Particle>> ParticlesSnapshotsMap;
protected:
    void SaveParticles()
    {
        ParticlesSnapshotsMap.emplace_back(Particles);
    }
};

#endif
