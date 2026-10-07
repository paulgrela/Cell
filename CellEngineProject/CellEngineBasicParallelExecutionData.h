
#ifndef CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H
#define CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H

#include <queue>
#include <condition_variable>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

struct __attribute__ ((packed)) ParticleSenderStructMPIMultiProcess
{
    UniqueIdUnsignedInt ParticleIndex{ 0 };
    EntityIdInt ParticleKindId{ 0 };
    int SenderProcessIndex{ 0 };
    int ReceiverProcessIndex{ 0 };
    ThreadPosType SenderThreadPos{ .ThreadPosX = 0, .ThreadPosY = 0, .ThreadPosZ = 0 };
    ThreadPosType ReceiverThreadPos{ .ThreadPosX = 0, .ThreadPosY = 0, .ThreadPosZ = 0 };
    vector3_16 SenderSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_16 ReceiverSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 NewPosition{ .X = 0, .Y = 0, .Z = 0 };
};

struct ParticleSenderStructMultiThreaded
{
    UniqueIdUnsignedInt ParticleIndex{ 0 };
    EntityIdInt ParticleKindId{ 0 };
    int SenderProcessIndex{ 0 };
    int ReceiverProcessIndex{ 0 };
    ThreadPosType SenderThreadPos{ .ThreadPosX = 0, .ThreadPosY = 0, .ThreadPosZ = 0 };
    ThreadPosType ReceiverThreadPos{ .ThreadPosX = 0, .ThreadPosY = 0, .ThreadPosZ = 0 };
    vector3_16 SenderSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_16 ReceiverSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 NewPosition{ .X = 0, .Y = 0, .Z = 0 };

    Particle ParticleObject;
};

struct ConfirmationOfParticlesToRemoveToSentStruct
{
    UniqueIdUnsignedInt ParticleIndex{ 0 };
    vector3_16 SenderSectorPos{ .X = 0, .Y = 0, .Z = 0 };
};

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
public:
    ThreadPosType CurrentThreadPos{ .ThreadPosX = 1, .ThreadPosY = 1, .ThreadPosZ = 1 };
public:
    ThreadPosType NeighborThreadsIndexes[NumberOfAllNeighbors];
public:
    std::vector<ParticleToBeMovedFromOneSectorToAnotherSector> ListOfParticlesToChangeSectors;
public:
    std::vector<ParticleSenderStructMultiThreaded> VectorOfParticlesToSendToNeighborThreads[NumberOfAllNeighbors];
    std::vector<ParticleSenderStructMPIMultiProcess> VectorOfParticlesToSendToNeighborProcesses[NumberOfAllNeighbors];
protected:
    SignedInt TwoThreadsWallSychronizationBarriersIndexes[NumberOfAllNeighbors]{ -1, -1, -1, -1, -1, -1 };
protected:
    std::vector<ConfirmationOfParticlesToRemoveToSentStruct> ConfirmationOfParticlesToRemoveToSent[NumberOfAllNeighbors];
public:
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
    ParticlesDetailedContainer<UniqueIdUnsignedInt> CancelledParticlesIndexes;
};

#endif
