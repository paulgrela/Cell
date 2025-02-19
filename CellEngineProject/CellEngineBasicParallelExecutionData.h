
#ifndef CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H
#define CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

class CellEngineBasicParallelExecutionData
{
    friend class CellEngineSimulationParallelExecutionManager;
protected:
    ThreadIdType CurrentThreadIndex{ 0 };
    CurrentThreadPosType CurrentThreadPos{ 1, 1, 1 };
protected:
    SectorPosType CurrentSectorPos{ 0, 0, 0 };
    SimulationSpaceSectorBounds ActualSimulationSpaceSectorBoundsObject{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
protected:
    std::mutex MainExchangeParticlesMutexObject;
protected:
    ParticlesContainerInternal<Particle> ParticlesForThreads;
protected:
    UnsignedInt ErrorCounter = 0;
    UnsignedInt AddedParticlesInReactions = 0;
protected:
    ParticlesContainerInternal<Particle> FormerParticlesIndexes;
    ParticlesContainerInternal<UniqueIdInt> CancelledParticlesIndexes;
};

#endif
