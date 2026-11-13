
#ifndef CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H
#define CELL_ENGINE_BASIC_PARALLEL_EXECUTION_DATA_H

#include <mpi.h>
#include <queue>
#include <barrier>
#include <condition_variable>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

struct ParticleSenderStructForParallelMultiTaskingExecution
{
    UniqueIdUInt64Type ParticleUniqueId{ 0 };
    UniqueIdUInt32Type ParticleKindId{ 0 };
    int SenderParallelTaskIndex{ 0 };
    int ReceiverParallelTaskIndex{ 0 };
    ParallelTaskPosType SenderParallelTaskPos{ .PosX = 0, .PosY = 0, .PosZ = 0 };
    ParallelTaskPosType ReceiverParallelTaskPos{ .PosX = 0, .PosY = 0, .PosZ = 0 };
    vector3_16 SenderSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_16 ReceiverSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 NewParticleCenter{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 Vector{ .X = 0, .Y = 0, .Z = 0 };

    Particle ParticleObject;
};

struct ConfirmationOfParticlesToRemoveToSentStruct
{
    UniqueIdUInt64Type ParticleUniqueId{ 0 };
    vector3_16 SenderSectorPos{ .X = 0, .Y = 0, .Z = 0 };
};

struct ParticleSenderStructForMPIExchange
{
    UniqueIdUInt64Type ParticleUniqueId{ 0 };
    UniqueIdUInt32Type ParticleKindId{ 0 };
    int SenderParallelTaskIndex{ 0 };
    int ReceiverParallelTaskIndex{ 0 };
    vector3_16 ReceiverSectorPos{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 NewParticleCenter{ .X = 0, .Y = 0, .Z = 0 };
    vector3_Real32 Vector{ .X = 0, .Y = 0, .Z = 0 };
};

struct MPIExchangeBuffersStruct
{
    std::vector<ParticleSenderStructForMPIExchange> ParticlesToSend[NumberOfAllNeighbors];
    std::vector<ParticleSenderStructForMPIExchange> ReceivedParticles[NumberOfAllNeighbors];
    std::vector<std::uint8_t> InsertedFlagsToSend[NumberOfAllNeighbors];
    std::vector<std::uint8_t> InsertedFlagsReceived[NumberOfAllNeighbors];
    int NumberOfParticlesSent[NumberOfAllNeighbors]{};
};

struct MPIExchangeRequestsStruct
{
    MPI_Request ReceiveParticlesRequests[NumberOfAllNeighbors];
    MPI_Request ReceiveInsertedFlagsRequests[NumberOfAllNeighbors];
    MPI_Request SendParticlesRequests[NumberOfAllNeighbors];
    MPI_Request SendInsertedFlagsRequests[NumberOfAllNeighbors];

    MPIExchangeRequestsStruct()
    {
        SetAllRequestsToNull();
    }

    void SetAllRequestsToNull()
    {
        for (UnsignedInt NeighborIndex = 0; NeighborIndex < NumberOfAllNeighbors; NeighborIndex++)
        {
            ReceiveParticlesRequests[NeighborIndex] = MPI_REQUEST_NULL;
            ReceiveInsertedFlagsRequests[NeighborIndex] = MPI_REQUEST_NULL;
            SendParticlesRequests[NeighborIndex] = MPI_REQUEST_NULL;
            SendInsertedFlagsRequests[NeighborIndex] = MPI_REQUEST_NULL;
        }
    }
};

struct WallChannel
{
    std::mutex MutexNumber1;
    std::mutex MutexNumber2;
    std::condition_variable ProposalConditionalVariable;
    std::condition_variable VerdictConditionalVariable;

    std::condition_variable WallConditionalVariable;

    std::vector<ParticleSenderStructForParallelMultiTaskingExecution> Proposals[2];
    std::vector<ConfirmationOfParticlesToRemoveToSentStruct> Verdicts[2];
    bool ProposalsReady[2]{ false, false };
    bool VerdictsReady[2]{ false, false };
};

inline UnsignedInt FirstDiffusionStepStampOfCurrentRun = 1;

class CellEngineBasicParallelTasksExecutionData
{
    friend class CellEngineSimulationParallelTasksExecutionManager;
protected:
        MPIExchangeBuffersStruct MPIExchangeBuffers;
        MPIExchangeRequestsStruct MPIExchangeRequests;

public:
    UnsignedInt MPIProcessIndex{ 0 };
public:
    ThreadIdType CurrentThreadIndex{ 0 };
    ParallelTaskPosType CurrentThreadPos{ .PosX = 1, .PosY = 1, .PosZ = 1 };
protected:
    UnsignedInt ParallelTaskGroupNumber;
    UnsignedInt NumberOfActiveNeighbors;
public:
    SignedInt NeighborParallelTasksIndexes[NumberOfAllNeighbors];
protected:
    SimulationSpaceSectorsRanges SimulationSpaceSectorsRangesForCurrentParallelTask;
public:
    #ifdef NEIGHBOR_PARALLEL_TASKS_POSITIONS
    ParallelTaskPosType NeighborParallelTasksPositions[NumberOfAllNeighbors];
    #endif
public:
    SectorPosType CurrentSectorPos{ .PosX = 0, .PosY = 0, .PosZ = 0 };
    SimulationSpaceSectorBounds ActualSimulationSpaceSectorBoundsObject{ .StartXPos = 0, .StartYPos = 0, .StartZPos = 0, .SizeX = 0, .SizeY = 0, .SizeZ = 0, .EndXPos = 0, .EndYPos = 0, .EndZPos = 0 };
public:
    std::vector<ParticleToBeMovedFromOneSectorToAnotherSector> ListOfParticlesToChangeSectors;
public:
    std::vector<ParticleSenderStructForParallelMultiTaskingExecution> VectorOfParticlesToSendToNeighborParallelTasks[NumberOfAllNeighbors];
protected:
    std::vector<ConfirmationOfParticlesToRemoveToSentStruct> ConfirmationOfParticlesToRemoveToSent[NumberOfAllNeighbors];
protected:
    std::vector<std::unique_ptr<WallChannel>> TwoThreadsWallChannels;
    std::vector<std::unique_ptr<std::barrier<>>> TwoThreadsWallSychronizationBarriers;
    SignedInt TwoThreadsWallSychronizationBarriersIndexes[NumberOfAllNeighbors]{ -1, -1, -1, -1, -1, -1 };
public:
    UnsignedInt CurrentDiffusionStepStamp = 0;
protected:
    std::vector<int> CPUsForThreads;
protected:
    std::deque<std::binary_semaphore> StartGatesForThreads;
public:
    std::mutex MainExchangeParticlesMutexObject;
protected:
    ParticlesDetailedContainer<Particle> ParticlesForThreads;
public:
    UnsignedInt NumberOfDNAParticles;
    UnsignedInt NumberOfAllParticlesInMPISimulation;
public:
    UnsignedInt NumberOfParticlesMovedInDiffusion = 0;
    UnsignedInt NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask = 0;
    UnsignedInt NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask = 0;
    UnsignedInt NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected = 0;
public:
    UnsignedInt NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask = 0;
    UnsignedInt NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask = 0;
    UnsignedInt NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask = 0;
    UnsignedInt NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask = 0;
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
    ParticlesDetailedContainer<UniqueIdUInt64Type> CancelledParticlesIndexes;
};

#endif
