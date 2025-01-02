
#ifndef CELLENGINEBASICPARALLELEXECUTIONDATA_H
#define CELLENGINEBASICPARALLELEXECUTIONDATA_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

class CellEngineBasicParallelExecutionData
{
    friend class CellEngineSimulationParallelExecutionManager;
protected:
    ThreadIdType CurrentThreadIndex{ 0 };
    CurrentThreadPosType CurrentThreadPos{ 1, 1, 1 };
protected:
    SimulationSpaceSectorBounds ActualSimulationSpaceSectorBoundsObject{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
protected:
    std::mutex MainExchangeParticlesMutexObject;
protected:
    ParticlesContainer<Particle> ParticlesForThreads;
protected:
    UnsignedInt ErrorCounter = 0;
    UnsignedInt AddedParticlesInReactions = 0;
};

#endif
