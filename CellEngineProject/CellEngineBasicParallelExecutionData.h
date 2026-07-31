
#ifndef CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H
#define CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H

#include <queue>
#include <barrier>
#include <condition_variable>

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
protected:
    UnsignedInt ProcessGroupNumber;
    UnsignedInt NumberOfActiveNeighbors;
public:
    SignedInt NeighborProcessesIndexes[NumberOfAllNeighbors];
protected:
    SimulationSpaceSectorsRanges CurrentMPIProcessSimulationSpaceSectorsRanges;
public:
    ThreadIdType CurrentThreadIndex{ 0 };
protected:
    ThreadPosType CurrentThreadPos{ 1, 1, 1 };
    ThreadPosType NeighborThreadsIndexes[NumberOfAllNeighbors];
public:
    std::vector<ParticleSenderStruct> VectorOfParticlesToSendToNeighborProcessesOrThreads[NumberOfAllNeighbors];
protected:
    std::vector<ParticleSenderStruct> ReceivedParticlesToInsertFromAllNeighborProcessesOrThreads[NumberOfAllNeighbors];
protected:
    std::unique_ptr<std::barrier<>> TwoThreadsWallSychronizationBarriers;
protected:
    std::vector<UniqueIdInt> ConfirmationOfParticlesToRemoveToSent[NumberOfAllNeighbors];
protected:
    std::condition_variable ProposalConditionalVariable;
    std::condition_variable VerdictConditionalVariable;
    bool ProposalsReady = false;
    bool VerdictsReady = false;
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
