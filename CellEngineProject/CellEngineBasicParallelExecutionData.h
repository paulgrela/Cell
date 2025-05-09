
#ifndef CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H
#define CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H

#include <queue>
#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

class CellEngineBasicParallelExecutionData
{
    friend class CellEngineSimulationParallelExecutionManager;
public:
    UnsignedInt GetMPIProcessIndex() const
    {
        return MPIProcessIndex;
    }
protected:
    UnsignedInt MPIProcessIndex{ 0 };
    UnsignedInt ProcessGroupNumber;
    SignedInt NeighbourProcessesIndexes[NumberOfAllNeighbours];
    UnsignedInt NumberOfActiveNeighbours;
    std::vector<MPIParticleSenderStruct> VectorOfParticlesToSendToNeighbourProcesses[NumberOfAllNeighbours];
    SimulationSpaceSectorsRanges CurrentMPIProcessSimulationSpaceSectorsRanges;
protected:
    ThreadIdType CurrentThreadIndex{ 0 };
    ThreadPosType CurrentThreadPos{ 1, 1, 1 };
protected:
    SectorPosType CurrentSectorPos{ 0, 0, 0 };
    SimulationSpaceSectorBounds ActualSimulationSpaceSectorBoundsObject{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
public:
    std::mutex MainExchangeParticlesMutexObject;
protected:
    ParticlesDetailedContainer<Particle> ParticlesForThreads;
protected:
    UnsignedInt ErrorCounter = 0;
    UnsignedInt NumberOfExecutedReactions = 0;
    UnsignedInt NumberOfCancelledReactions = 0;
    UnsignedInt NumberOfCancelledAReactions = 0;
    UnsignedInt NumberOfCancelledBReactions = 0;
    UnsignedInt AddedParticlesInReactions = 0;
    UnsignedInt RemovedParticlesInReactions = 0;
    UnsignedInt RestoredParticlesInCancelledReactions = 0;
protected:
    ParticlesDetailedContainer<Particle> FormerParticlesIndexes;
    ParticlesDetailedContainer<UniqueIdInt> CancelledParticlesIndexes;
};

#endif
