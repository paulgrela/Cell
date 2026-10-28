
#include <mpi.h>
#include <chrono>

#include "DateTimeUtils.h"

#include "CellEngineDataFile.h"
#include "CellEngineImGuiMenu.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineSimulationParallelExecutionManager.h"

#ifdef USE_OPENGL
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"
#endif

constexpr bool PrintAdditionalInformation = false;

using namespace std;

namespace
{
    struct WallChannel
    {
        std::mutex MutexNumber1;
        std::mutex MutexNumber2;
        std::condition_variable ProposalConditionalVariable;
        std::condition_variable VerdictConditionalVariable;

        std::condition_variable WallConditionalVariable;

        std::vector<ParticleSenderStructMultiThreaded> Proposals[2];
        std::vector<UniqueIdUInt64Type> Verdicts[2];
        bool ProposalsReady[2]{ false, false };
        bool VerdictsReady[2]{ false, false };
    };
}

namespace
{
    static std::vector<std::unique_ptr<WallChannel>> TwoThreadsWallChannels;
}

namespace
{
    static std::vector<std::unique_ptr<std::barrier<>>> TwoThreadsWallSychronizationBarriers;
}

CellEngineSimulationParallelExecutionManager::CellEngineSimulationParallelExecutionManager() : SimulationSpaceDataForThreads(CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer)
{
}

template void CellEngineSimulationParallelExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineVoxelSimulationSpace>(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);
template void CellEngineSimulationParallelExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineFullAtomSimulationSpace>(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);

template <class SimulationSpaceType>
void CellEngineSimulationParallelExecutionManager::CreateSimulationSpaceForParallelExecution(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles)
{
    try
    {
        UnsignedInt ThreadIndexPos = 1;

        CellEngineSimulationSpaceForThreadsObjectsPointer.clear();
        CellEngineSimulationSpaceForThreadsObjectsPointer.resize(CellEngineConfigDataObject.NumberOfXThreadsInSimulation);
        UnsignedInt ThreadXPos = 1;
        for (auto& ThreadLocalParticlesInProximityXPos : CellEngineSimulationSpaceForThreadsObjectsPointer)
        {
            ThreadLocalParticlesInProximityXPos.resize(CellEngineConfigDataObject.NumberOfYThreadsInSimulation);
            UnsignedInt ThreadYPos = 1;
            for (auto& ThreadLocalParticlesInProximityYPos : ThreadLocalParticlesInProximityXPos)
            {
                ThreadLocalParticlesInProximityYPos.resize(CellEngineConfigDataObject.NumberOfZThreadsInSimulation);
                UnsignedInt ThreadZPos = 1;
                for (auto& ThreadLocalParticlesInProximityZPos : ThreadLocalParticlesInProximityYPos)
                {
                    LoggersManagerObject.Log(STREAM("THREAD INDEXES = " << ThreadIndexPos << " (" << ThreadXPos << ", " << ThreadYPos << ", " << ThreadZPos << ")"));

                    ThreadLocalParticlesInProximityZPos = make_shared<SimulationSpaceType>(Particles, false, ThreadIndexPos, ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(ThreadXPos), .ThreadPosY = static_cast<SignedInt>(ThreadYPos), .ThreadPosZ = static_cast<SignedInt>(ThreadZPos) });
                    ThreadIndexPos++;
                    ThreadZPos++;
                }
                ThreadYPos++;
            }
            ThreadXPos++;
        }
    }
    CATCH("creating simulation space for parallel execution")
}

static UnsignedInt GetProcessGroupNumberVer2(const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    if (ThreadZIndex % 2 == 0)
        return ThreadXIndex % 2 != ThreadYIndex % 2;
    else
        return ThreadXIndex % 2 == ThreadYIndex % 2;
}

SignedInt CellEngineSimulationParallelExecutionManager::GetProcessPrevNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (ThreadXIndex >= 0 && ThreadYIndex >= 0 && ThreadZIndex >= 0)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex - 1;
    else
        return -1;
}

SignedInt CellEngineSimulationParallelExecutionManager::GetProcessNextNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (static_cast<SignedInt>(ThreadXIndex) < CellEngineConfigDataObject.NumberOfXThreadsInSimulation && static_cast<SignedInt>(ThreadYIndex) < CellEngineConfigDataObject.NumberOfYThreadsInSimulation && static_cast<SignedInt>(ThreadZIndex) < CellEngineConfigDataObject.NumberOfZThreadsInSimulation)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex - 1;
    else
        return -1;
}

void CellEngineSimulationParallelExecutionManager::CreateDataForEveryMPIProcessForParallelExecution()
{
    try
    {
        SignedInt MPIProcessIndex = 0;

        for (UnsignedInt MPIProcessXIndex = 1; MPIProcessXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; MPIProcessXIndex++)
            for (UnsignedInt MPIProcessYIndex = 1; MPIProcessYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; MPIProcessYIndex++)
                for (UnsignedInt MPIProcessZIndex = 1; MPIProcessZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; MPIProcessZIndex++)
                {
                    if (MPIProcessIndex == MPIProcessDataObject.CurrentMPIProcessIndex)
                    {
                        CurrentMPIProcessSimulationSpaceSectorsRanges.SetParameters((MPIProcessXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, (MPIProcessYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, (MPIProcessZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation, MPIProcessXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, MPIProcessYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, MPIProcessZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation);

                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Bounds = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " (" << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << ") (" << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << ")"));

                        ProcessGroupNumber = GetProcessGroupNumberVer2(MPIProcessXIndex - 1, MPIProcessYIndex - 1, MPIProcessZIndex - 1);

                        LoggersManagerObject.Log(STREAM("ProcessGroupNumber = " << ProcessGroupNumber << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ")"));

                        NeighborProcessesIndexes[0] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 2, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighborProcessesIndexes[1] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 2, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighborProcessesIndexes[2] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 2);

                        NeighborProcessesIndexes[3] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighborProcessesIndexes[4] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighborProcessesIndexes[5] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex));

                        ListOfParticlesToChangeSectors.reserve(65536);

                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                            if (NeighborProcessesIndexes[NeighborProcessIndex] != -1)
                            {
                                VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].reserve(16384);
                                VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].reserve(16384);
                            }

                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Neighbors = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " " << NeighborProcessesIndexes[0] << " " << NeighborProcessesIndexes[1] << " " << NeighborProcessesIndexes[2] << " " << NeighborProcessesIndexes[3] << " " << NeighborProcessesIndexes[4] << " " << NeighborProcessesIndexes[5]));
                    }

                    MPIProcessIndex++;
                }

        NumberOfActiveNeighbors = count_if(NeighborProcessesIndexes, NeighborProcessesIndexes + NumberOfAllNeighbors, [](const auto& Element){ return Element != -1; });

        LoggersManagerObject.Log(STREAM("NumberOfActiveNeighbors = " << NumberOfActiveNeighbors));
    }
    CATCH("creating data for every mpi process for parallel execution")
}

void CellEngineSimulationParallelExecutionManager::CreateDataForEveryThreadForParallelExecution() const
{
    try
    {
        SignedInt MPIProcessIndex = 0;

        for (UnsignedInt MPIProcessXIndex = 1; MPIProcessXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; MPIProcessXIndex++)
            for (UnsignedInt MPIProcessYIndex = 1; MPIProcessYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; MPIProcessYIndex++)
                for (UnsignedInt MPIProcessZIndex = 1; MPIProcessZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; MPIProcessZIndex++)
                {
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->CurrentMPIProcessSimulationSpaceSectorsRanges.SetParameters((MPIProcessXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, (MPIProcessYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, (MPIProcessZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation, MPIProcessXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, MPIProcessYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, MPIProcessZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation);

                    if constexpr (PrintAdditionalInformation == true)
                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Bounds = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " (" << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << ") (" << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << ")"));

                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->MPIProcessIndex = MPIProcessIndex - 1;
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->ProcessGroupNumber = GetProcessGroupNumberVer2(MPIProcessXIndex - 1, MPIProcessYIndex - 1, MPIProcessZIndex - 1);

                    if constexpr (PrintAdditionalInformation == true)
                        LoggersManagerObject.Log(STREAM("ProcessGroupNumber = " << ProcessGroupNumber << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ")"));

                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[0] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 2, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[1] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 2, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[2] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 2);

                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[3] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[4] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex) - 1);
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes[5] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex));

                    #ifdef NEIGHBOR_THREAD_INDEXES
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[0] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex) - 2, .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex) - 1, .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) - 1 };
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[1] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex) - 1, .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex) - 2, .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) - 1 };
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[2] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex) - 1, .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex) - 1, .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) - 2 };

                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[3] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex), .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex) - 1, .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) - 1 };
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[4] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex) - 1, .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex), .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) - 1 };
                    SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes[5] = ThreadPosType{ .ThreadPosX = static_cast<SignedInt>(MPIProcessXIndex) - 1, .ThreadPosY = static_cast<SignedInt>(MPIProcessYIndex) - 1, .ThreadPosZ = static_cast<SignedInt>(MPIProcessZIndex) };
                    #endif

                    if constexpr (PrintAdditionalInformation == true)
                    {
                        const auto NPI = SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborProcessesIndexes;
                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Neighbors = " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << NPI[0] << " " << NPI[1] << " " << NPI[2] << " " << NPI[3] << " " << NPI[4] << " " << NPI[5]));
                        #ifdef NEIGHBOR_THREAD_INDEXES
                        const auto NTI = SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->NeighborThreadsIndexes;
                        LoggersManagerObject.Log(STREAM("Thread Neighbors Indexes = " << SimulationSpaceDataForThreads[MPIProcessXIndex - 1][MPIProcessYIndex - 1][MPIProcessZIndex - 1]->CurrentThreadIndex << " (" << NTI[0].ThreadPosX << " " << NTI[0].ThreadPosY << " " << NTI[0].ThreadPosZ << ") (" << NTI[1].ThreadPosX << " " << NTI[1].ThreadPosY << " " << NTI[1].ThreadPosZ << ") (" << NTI[2].ThreadPosX << " " << NTI[2].ThreadPosY << " " << NTI[2].ThreadPosZ << ") (" << NTI[3].ThreadPosX << " " << NTI[3].ThreadPosY << " " << NTI[3].ThreadPosZ << ") (" << NTI[4].ThreadPosX << " " << NTI[4].ThreadPosY << " " << NTI[4].ThreadPosZ << ") (" << NTI[5].ThreadPosX << " " << NTI[5].ThreadPosY << " " << NTI[5].ThreadPosZ << ")"));
                        #endif
                    }

                    MPIProcessIndex++;
                }
    }
    CATCH("creating data for every mpi process for parallel execution")
}

constexpr SignedInt NeighborOffsetX[NumberOfAllNeighbors] = { -1,  0,  0, +1,  0,  0 };
constexpr SignedInt NeighborOffsetY[NumberOfAllNeighbors] = {  0, -1,  0,  0, +1,  0 };
constexpr SignedInt NeighborOffsetZ[NumberOfAllNeighbors] = {  0,  0, -1,  0,  0, +1 };

SignedInt CellEngineSimulationParallelExecutionManager::GetProcessNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (ThreadXIndex >= 0 && ThreadXIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfXThreadsInSimulation) && ThreadYIndex >= 0 && ThreadYIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation) && ThreadZIndex >= 0 && ThreadZIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation))
        return static_cast<SignedInt>(SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex);
    else
        return -1;
}

void CellEngineSimulationParallelExecutionManager::CreateDataForEveryThreadForParallelExecutionADD() const
{
    try
    {
        cout << "STEP 0" << endl;

        SignedInt MPIProcessIndex = 0;

        for (UnsignedInt MPIProcessXIndex = 0; MPIProcessXIndex < CellEngineConfigDataObject.NumberOfXThreadsInSimulation; MPIProcessXIndex++)
            for (UnsignedInt MPIProcessYIndex = 0; MPIProcessYIndex < CellEngineConfigDataObject.NumberOfYThreadsInSimulation; MPIProcessYIndex++)
                for (UnsignedInt MPIProcessZIndex = 0; MPIProcessZIndex < CellEngineConfigDataObject.NumberOfZThreadsInSimulation; MPIProcessZIndex++)
                {
                    if constexpr (PrintAdditionalInformation == true)
                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Bounds = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " (" << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << ") (" << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << "," << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << ")"));

                    SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex]->MPIProcessIndex = MPIProcessIndex;

                    if constexpr (PrintAdditionalInformation == true)
                        LoggersManagerObject.Log(STREAM("ProcessGroupNumber = " << ProcessGroupNumber << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ")"));

                    for (UnsignedInt DirectionIndex = 0; DirectionIndex < NumberOfAllNeighbors; DirectionIndex++)
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[DirectionIndex] = -1;

                    MPIProcessIndex++;
                }

        cout << "STEP 1" << endl;

        for (SignedInt ThreadXIndex = 0; ThreadXIndex < CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (SignedInt ThreadYIndex = 0; ThreadYIndex < CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (SignedInt ThreadZIndex = 0; ThreadZIndex < CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    const auto& ThreadData = SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex];
                    ThreadData->NumberOfActiveNeighbors = 0;

                    for (UnsignedInt DirectionIndex = 0; DirectionIndex < NumberOfAllNeighbors; DirectionIndex++)
                    {
                        const SignedInt NeighborXIndex = ThreadXIndex + NeighborOffsetX[DirectionIndex];
                        const SignedInt NeighborYIndex = ThreadYIndex + NeighborOffsetY[DirectionIndex];
                        const SignedInt NeighborZIndex = ThreadZIndex + NeighborOffsetZ[DirectionIndex];

                        ThreadData->NeighborProcessesIndexes[DirectionIndex] = GetProcessNeighbor(NeighborXIndex, NeighborYIndex, NeighborZIndex);
                        #ifdef NEIGHBOR_THREAD_INDEXES
                        ThreadData->NeighborThreadsIndexes[DirectionIndex] = ThreadPosType{ .ThreadPosX = NeighborXIndex, .ThreadPosY = NeighborYIndex, .ThreadPosZ = NeighborZIndex };
                        #endif

                        if (ThreadData->NeighborProcessesIndexes[DirectionIndex] != -1)
                            ThreadData->NumberOfActiveNeighbors++;
                    }
                }

        cout << "STEP 2" << endl;

        TwoThreadsWallSychronizationBarriers.clear();
        UnsignedInt TwoThreadsWallSychronizationBarriersNumber = 0;

        for (UnsignedInt MPIProcessXIndex = 0; MPIProcessXIndex < CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ++MPIProcessXIndex)
            for (UnsignedInt MPIProcessYIndex = 0; MPIProcessYIndex < CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ++MPIProcessYIndex)
                for (UnsignedInt MPIProcessZIndex = 0; MPIProcessZIndex < CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ++MPIProcessZIndex)
                {
                    if (MPIProcessXIndex + 1 < CellEngineConfigDataObject.NumberOfXThreadsInSimulation)
                    {
                        TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                        TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[3] = TwoThreadsWallSychronizationBarriersNumber;
                        SimulationSpaceDataForThreads[MPIProcessXIndex + 1][MPIProcessYIndex][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[0] = TwoThreadsWallSychronizationBarriersNumber;
                        cout << "THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ") with CUBE (" << MPIProcessXIndex + 1 << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ")" << endl;
                        ++TwoThreadsWallSychronizationBarriersNumber;
                    }
                    if (MPIProcessYIndex + 1 < CellEngineConfigDataObject.NumberOfYThreadsInSimulation)
                    {
                        TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                        TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[4] = TwoThreadsWallSychronizationBarriersNumber;
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex + 1][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[1] = TwoThreadsWallSychronizationBarriersNumber;
                        cout << "THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ") with CUBE (" << MPIProcessXIndex << "," << MPIProcessYIndex + 1 << "," << MPIProcessZIndex << ")" << endl;
                        ++TwoThreadsWallSychronizationBarriersNumber;
                    }
                    if (MPIProcessZIndex + 1 < CellEngineConfigDataObject.NumberOfZThreadsInSimulation)
                    {
                        TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                        TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex]->TwoThreadsWallSychronizationBarriersIndexes[5] = TwoThreadsWallSychronizationBarriersNumber;
                        SimulationSpaceDataForThreads[MPIProcessXIndex][MPIProcessYIndex][MPIProcessZIndex + 1]->TwoThreadsWallSychronizationBarriersIndexes[2] = TwoThreadsWallSychronizationBarriersNumber;
                        cout << "THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ") with CUBE (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex + 1 << ")" << endl;
                        ++TwoThreadsWallSychronizationBarriersNumber;
                    }
                }

        cout << "WALLS NUMBER 2 = " << TwoThreadsWallSychronizationBarriersNumber << endl;
    }
    CATCH("creating data for every mpi process for parallel execution")
}

#define FOR_EACH_THREAD_IN_XYZ \
        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++) \
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++) \
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)

#ifdef SHORTER_CODE
void CellEngineSimulationParallelExecutionManager::JoinReactionsStatisticsFromThreads(vector<map<UnsignedInt, ChemicalReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true)
        {
            if (CellEngineConfigDataObject.OpenGLGraphicsSwitchedOff == false && MPIProcessDataObject.CurrentMPIProcessIndex == CellEngineConfigDataObject.MainMPIProcessNumber)
            {
                int ValueToSend = 3;
                MPI_Bcast(&ValueToSend, 1, MPI_INT, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);
            }

            for (const auto& [ReactionStatisticsDataFirst, ReactionStatisticsDataSecond] : SavedReactionsMap.back())
            {
                LoggersManagerObject.Log(STREAM("ReactionStatisticsData S = " << ReactionStatisticsDataFirst << " " << ReactionStatisticsDataSecond.Counter << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

                CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].emplace_back(ReactionStatisticsDataFirst, ReactionStatisticsDataSecond.Counter);
            }

            if (CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].empty() == true)
                CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].emplace_back(10001, 1);

            const int SavedReactionsMapForMPILength = CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].size();
            int SavedReactionsMapForMPILengths[MPIProcessDataObject.NumberOfMPIProcesses];
            MPI_Gather(&SavedReactionsMapForMPILength, 1, MPI_INT, SavedReactionsMapForMPILengths, 1, MPI_INT, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);

            LoggersManagerObject.Log(STREAM("SavedReactionsMapForMPILength = " << SavedReactionsMapForMPILength << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

            int MaximumOfAllSavedReactionsMapForMPILengths;
            if (MPIProcessDataObject.CurrentMPIProcessIndex == CellEngineConfigDataObject.MainMPIProcessNumber)
                MaximumOfAllSavedReactionsMapForMPILengths = *max_element(SavedReactionsMapForMPILengths, SavedReactionsMapForMPILengths + MPIProcessDataObject.NumberOfMPIProcesses);

            MPI_Bcast(&MaximumOfAllSavedReactionsMapForMPILengths, 1, MPI_INT, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);

            LoggersManagerObject.Log(STREAM("MaximumOfAllSavedReactionsMapForMPILengths = " << MaximumOfAllSavedReactionsMapForMPILengths << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

            const unique_ptr<ChemicalReactionStatistics[]> ReactionStatisticsVectorGatheringFromAllMPIProcessPointer(new ChemicalReactionStatistics[MaximumOfAllSavedReactionsMapForMPILengths * MPIProcessDataObject.NumberOfMPIProcesses + 1]);

            const int NumberOfBytesToGatherFromEveryProcess = static_cast<int>(MaximumOfAllSavedReactionsMapForMPILengths * sizeof(ChemicalReactionStatistics));
            MPI_Gather(CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].data(), NumberOfBytesToGatherFromEveryProcess, MPI_BYTE, ReactionStatisticsVectorGatheringFromAllMPIProcessPointer.get(), NumberOfBytesToGatherFromEveryProcess, MPI_BYTE, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);

            if (MPIProcessDataObject.CurrentMPIProcessIndex == CellEngineConfigDataObject.MainMPIProcessNumber)
            {
                SavedReactionsMap[SimulationStepNumber - 1].clear();
                for (UniqueIdUInt64Type LocalMPIProcessIndex = 0; LocalMPIProcessIndex < MPIProcessDataObject.NumberOfMPIProcesses; LocalMPIProcessIndex++)
                    for (UnsignedInt ReactionStatisticsDataIndex = MaximumOfAllSavedReactionsMapForMPILengths * LocalMPIProcessIndex; ReactionStatisticsDataIndex < MaximumOfAllSavedReactionsMapForMPILengths * LocalMPIProcessIndex + SavedReactionsMapForMPILengths[LocalMPIProcessIndex]; ReactionStatisticsDataIndex++)
                    {
                        const auto& [ReactionId, Counter] = ReactionStatisticsVectorGatheringFromAllMPIProcessPointer.get()[ReactionStatisticsDataIndex];
                        LoggersManagerObject.Log(STREAM("ReactionStatisticsData 1 = " << ReactionStatisticsDataIndex << " " << ReactionId << " " << Counter << " LocalMPIProcessIndex = " << LocalMPIProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                        if (ReactionId != 0 && Counter != 0)
                        {
                            LoggersManagerObject.Log(STREAM("ReactionStatisticsData 2 = " << ReactionId << " " << ReactionStatisticsDataIndex << " " << Counter << " LocalMPIProcessIndex = " << LocalMPIProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                            SavedReactionsMap[SimulationStepNumber - 1][ReactionId] = { ReactionId, SavedReactionsMap[SimulationStepNumber - 1][ReactionId].Counter += Counter };
                        }
                    }
            }
        }
        else
            FOR_EACH_THREAD_IN_XYZ
                for (const auto& [ReactionStatisticsDataFirst, ReactionStatisticsDataSecond] : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap.back())
                    SavedReactionsMap[SimulationStepNumber - 1][ReactionStatisticsDataFirst].Counter += ReactionStatisticsDataSecond.Counter;
    }
    CATCH("joining statistics from threads")
}
#else
void CellEngineSimulationParallelExecutionManager::JoinStatisticsFromThreads(vector<map<UnsignedInt, ChemicalReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ReactionData : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap[SimulationStepNumber - 1])
                if (SavedReactionsMap[SimulationStepNumber - 1].contains(ReactionData.second.ChemicalReactionId))
                    SavedReactionsMap[SimulationStepNumber - 1].find(ReactionData.second.ChemicalReactionId)->second.Counter += ReactionData.second.Counter;
                else
                    SavedReactionsMap[SimulationStepNumber - 1].insert(ReactionData);
    }
    CATCH("joining statistics from threads extended")
}
#endif

static void LogCenterOfParticleWithThreadIndex(const Particle& ParticleObject, const ThreadIdType ThreadXIndex, const ThreadIdType ThreadYIndex, const ThreadIdType ThreadZIndex)
{
    LoggersManagerObject.Log(STREAM("Center: " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << endl));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl));
    LoggersManagerObject.Log(STREAM(endl));
}

void CellEngineSimulationParallelExecutionManager::FirstSendParticlesForThreads(const bool PrintCenterOfParticleWithThreadIndex, const bool PrintTime)
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("THREAD = " << CurrentThreadIndex << " X = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInX << " Y = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInY << " Z = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInZ));

        UnsignedInt GoodParticlesCounter = 0;
        UnsignedInt BadParticlesCounter = 0;

        SaveFormerParticlesAsVectorElements();

        const auto start_time = chrono::high_resolution_clock::now();

        FOR_EACH_THREAD_IN_XYZ
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.clear();

        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
            if (ParticleObject.second.Center.X < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && ParticleObject.second.Center.Y < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && ParticleObject.second.Center.Z < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
            {
                const UnsignedInt ThreadXIndex = floor(ParticleObject.second.Center.X / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
                const UnsignedInt ThreadYIndex = floor(ParticleObject.second.Center.Y / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
                const UnsignedInt ThreadZIndex = floor(ParticleObject.second.Center.Z / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if (PrintCenterOfParticleWithThreadIndex == true)
                    LogCenterOfParticleWithThreadIndex(ParticleObject.second, ThreadXIndex, ThreadYIndex, ThreadZIndex);

                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert(ParticleObject);

                GoodParticlesCounter++;
            }
            else
                BadParticlesCounter++;

        FOR_EACH_THREAD_IN_XYZ
            if (PrintCenterOfParticleWithThreadIndex == true)
                LoggersManagerObject.Log(STREAM("THREAD[" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << "] SIZE = " << SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.size()));

        const auto stop_time = chrono::high_resolution_clock::now();

        if (PrintTime == true)
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "First sending particles to threads has taken time = ","Execution in threads")));

        LoggersManagerObject.Log(STREAM("All particles counter = " << GoodParticlesCounter << " Bad particles counter = " << BadParticlesCounter));
    }
    CATCH("first sending particles for threads")
}

static void PrintThreadIndexes(const Particle& ParticleObject, const ThreadIdType CurrentThreadIndex, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const UnsignedInt ThreadXIndexNew, const UnsignedInt ThreadYIndexNew, const UnsignedInt ThreadZIndexNew)
{
    LoggersManagerObject.Log(STREAM("Particle Center: " << ParticleObject.ParticleUniqueId << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticlesKindsManagerObject.GetParticleKind(ParticleObject.ParticleKindId).ParticleKindSpecialDataSector[0].ParticleType)));
    LoggersManagerObject.Log(STREAM("THREAD POS NEW = " << ThreadXIndexNew << ", " << ThreadYIndexNew << ", " << ThreadZIndexNew));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex));
    LoggersManagerObject.Log(STREAM(endl));
}

static void PrintDataAboutParticlesExchange(const bool PrintInfo, const chrono::time_point<chrono::system_clock, chrono::system_clock::duration>& start_time, const chrono::time_point<chrono::system_clock, chrono::system_clock::duration>& stop_time, const UnsignedInt ExchangedParticleCounter, const bool StateOfSimulationSpaceDivisionForThreads)
{
    if (PrintInfo == true)
    {
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Exchanging particles between threads has taken time = ","Execution in threads")));
        LoggersManagerObject.Log(STREAM("Number of exchanged particles: " << ExchangedParticleCounter));
        LoggersManagerObject.Log(STREAM("StateOfSimulationSpaceDivisionForThreads: " << StateOfSimulationSpaceDivisionForThreads));
    }
}

vector<vector<vector<ParticlesContainerInternal<Particle>>>> CellEngineSimulationParallelExecutionManager::GatherParticlesToExchangeBetweenThreads(const UnsignedInt TypeOfGet, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, UnsignedInt& ExchangedParticleCounter, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<ParticlesContainerInternal<Particle>>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<ParticlesContainerInternal<Particle>>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

    try
    {
        auto ParticleIter = SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.begin();

        while (ParticleIter != SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.end())
        {
            const UnsignedInt SizeOfXInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2));
            const UnsignedInt SizeOfYInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2));
            const UnsignedInt SizeOfZInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2));

            if (StateOfSimulationSpaceDivisionForThreads == true)
                if (ParticleIter->second.Center.X < SizeOfXInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Y < SizeOfYInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Z < SizeOfZInOneThreadInSimulationSpaceDiv2)
                {
                    if (PrintInfo == true)
                        LogCenterOfParticleWithThreadIndex(ParticleIter->second, ThreadXIndex, ThreadYIndex, ThreadZIndex);
                    ++ParticleIter;
                    continue;
                }

            const UnsignedInt ThreadXIndexNew = floor((ParticleIter->second.Center.X - SizeOfXInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
            const UnsignedInt ThreadYIndexNew = floor((ParticleIter->second.Center.Y - SizeOfYInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
            const UnsignedInt ThreadZIndexNew = floor((ParticleIter->second.Center.Z - SizeOfZInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

            if (PrintInfo == true)
                PrintThreadIndexes(ParticleIter->second, CurrentThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex, ThreadXIndexNew, ThreadYIndexNew, ThreadZIndexNew);

            if ((ThreadXIndex - 1 != ThreadXIndexNew) || (ThreadYIndex - 1 != ThreadYIndexNew) || (ThreadZIndex - 1 != ThreadZIndexNew))
            {
                switch (TypeOfGet)
                {
                    case 0 : SimulationSpaceDataForThreads[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew]->ParticlesForThreads.insert(SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.extract(ParticleIter++)); break;
                    case 1 : ParticlesToExchange[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew].Particles.insert(SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.extract(ParticleIter++)); break;
                    default : break;
                }

                ExchangedParticleCounter++;
            }
            else
                ++ParticleIter;
        }
    }
    CATCH("gathering particles to exchange between threads")

    return ParticlesToExchange;
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreads(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        UnsignedInt ExchangedParticleCounter = 0;

        FOR_EACH_THREAD_IN_XYZ
            GatherParticlesToExchangeBetweenThreads(0, ThreadXIndex, ThreadYIndex, ThreadZIndex, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads, PrintInfo);

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintDataAboutParticlesExchange(PrintInfo, start_time, stop_time, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads);
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsParallelInsert(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCounter = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange;

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject };

            ParticlesToExchange = GatherParticlesToExchangeBetweenThreads(1, CurrentThreadPos.ThreadPosX, CurrentThreadPos.ThreadPosY, CurrentThreadPos.ThreadPosZ, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads, PrintInfo);
        }

        FOR_EACH_THREAD_IN_XYZ
            if (ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles.empty() == false)
            {
                lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->MainExchangeParticlesMutexObject };

                for (const auto& ParticleToExchange : ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles)
                    SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.insert(ParticleToExchange);
            }

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintDataAboutParticlesExchange(PrintInfo, start_time, stop_time, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads);
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsParallelExtract(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCounter = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange;

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject };

            ParticlesToExchange = GatherParticlesToExchangeBetweenThreads(1, CurrentThreadPos.ThreadPosX, CurrentThreadPos.ThreadPosY, CurrentThreadPos.ThreadPosZ, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads, PrintInfo);
        }

        FOR_EACH_THREAD_IN_XYZ
            if (ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles.empty() == false)
            {
                lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->MainExchangeParticlesMutexObject };

                auto ParticleIterLocalInside = ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles.begin();
                while (ParticleIterLocalInside != ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles.end())
                    SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.insert(ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].Particles.extract(ParticleIterLocalInside++));
            }

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintDataAboutParticlesExchange(PrintInfo, start_time, stop_time, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads);
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::GatherParticlesFromThreads()
{
    try
    {
        GetParticles().clear();
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ParticleObject : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads)
                GetParticles().insert(ParticleObject);
    }
    CATCH("gathering cancelled from threads")
}

void CellEngineSimulationParallelExecutionManager::GatherCancelledParticlesIndexesFromThreads()
{
    try
    {
        CancelledParticlesIndexes.clear();

        FOR_EACH_THREAD_IN_XYZ
        {
            for (const auto& CancelledParticleUniqueId : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->CancelledParticlesIndexes)
                CancelledParticlesIndexes.insert(CancelledParticleUniqueId);

            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->CancelledParticlesIndexes.clear();
        }

        LoggersManagerObject.Log(STREAM("NUMBER OF CANCELLED PARTICLES IN CANCELLED REACTIONS = " << CancelledParticlesIndexes.size()));
    }
    CATCH("gathering particles from threads")
}

void CellEngineSimulationParallelExecutionManager::GatherParticlesFromThreadsToParticlesInMainThread()
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        sleep(1);
        {
            GatherParticlesFromThreads();

            GatherCancelledParticlesIndexesFromThreads();
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Gathering particles from threads to main thread has taken time = ","Execution in threads")));
    }
    CATCH("gathering particles from threads to particles in main thread")
}

void CellEngineSimulationParallelExecutionManager::SychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        switch (CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInThreads)
        {
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByOneGlobalBarrier : SyncPoint->arrive_and_wait(); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByLocalBarriersOnWallsBetweenThreads : SynchronizeWithNeighborByLocalBarrier(CurrentThreadLocalSimulationSpaceData); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithTwoMutexes : SyncPoint->arrive_and_wait(); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithOneMutex : SyncPoint->arrive_and_wait(); break;
            default : break;
        };
    }
    CATCH("synchronize simulation executed in threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsAndSychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam)
{
    try
    {
        switch (CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInThreads)
        {
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByOneGlobalBarrier : ExchangeParticlesBetweenThreadsOneGlobalBarrier(SyncPoint, CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByLocalBarriersOnWallsBetweenThreads : ExchangeParticlesBetweenThreadsLocalBarrier(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithTwoMutexes : ExchangeParticlesBetweenThreadsConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithOneMutex : ExchangeParticlesBetweenThreadsConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData); break;
            default : break;
        }
    }
    CATCH("exchange particles and synchronize simulation executed in threads")
}

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                SimulationSpaceSectorBounds SimulationSpaceSectorBoundsObject = SimulationSpaceSectorBoundsObject.SetParametersForParallelExecutionSectors({ static_cast<SignedInt>(ThreadXIndex), static_cast<SignedInt>(ThreadYIndex), static_cast<SignedInt>(ThreadZIndex) }, CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if (StateOfSimulationSpaceDivisionForThreads == true)
                    SimulationSpaceSectorBoundsObject.AddToStartParameters(CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2);

                for (UnsignedInt PosX = SimulationSpaceSectorBoundsObject.StartXPos; PosX < SimulationSpaceSectorBoundsObject.EndXPos; PosX += CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace)
                    for (UnsignedInt PosY = SimulationSpaceSectorBoundsObject.StartYPos; PosY < SimulationSpaceSectorBoundsObject.EndYPos; PosY += CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace)
                        for (UnsignedInt PosZ = SimulationSpaceSectorBoundsObject.StartZPos; PosZ < SimulationSpaceSectorBoundsObject.EndZPos; PosZ += CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace)
                        {
                            DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("XStart = " << SimulationSpaceSectorBoundsObject.StartXPos << " YStart = " << SimulationSpaceSectorBoundsObject.StartYPos << " ZStart = " << SimulationSpaceSectorBoundsObject.StartZPos << " XEnd = " << SimulationSpaceSectorBoundsObject.EndXPos << " YEnd = " << SimulationSpaceSectorBoundsObject.EndYPos << " ZEnd = " << SimulationSpaceSectorBoundsObject.EndZPos << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInThreads(CurrentThreadLocalSimulationSpaceData, true, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObject.StartXPos, SimulationSpaceSectorBoundsObject.StartYPos, SimulationSpaceSectorBoundsObject.StartZPos, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, true);
                        }
            }
        else
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                for (UnsignedInt ParticleSectorXIndex = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex < ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex < ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex < ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("XStart = " << (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " XEnd = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, true);
                        }

                SychronizeSimulationExecutedInThreads(SyncPoint, CurrentThreadLocalSimulationSpaceData);

                for (UnsignedInt ParticleSectorXIndex = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex < ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex < ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex < ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("XStart = " << (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " YStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace << " ZStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace << " XEnd = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " YEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace << " ZEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("RANGE -> START = (" << (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << ", " <<  (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," <<  (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") END = (" <<  ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << "," <<  ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," <<  ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInThreads(CurrentThreadLocalSimulationSpaceData, true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, (ThreadZIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ, ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ);
                        }

                SychronizeSimulationExecutedInThreads(SyncPoint, CurrentThreadLocalSimulationSpaceData);
            }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

static inline UnsignedInt StepToChangeSimulationSpaceDivisionForThreads(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads)
{
    return ((StepOutside % CellEngineConfigDataObject.StepToChangeSpaceDivisionForThreads == 0) ? !StateOfSimulationSpaceDivisionForThreads : StateOfSimulationSpaceDivisionForThreads);
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
        {
            CurrentThreadLocalSimulationSpaceData->GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(CurrentThreadLocalSimulationSpaceData, NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam, *StateOfSimulationSpaceDivisionForThreads, SyncPoint);

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            {
                SyncPoint->arrive_and_wait();

                if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert || CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                {
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert)
                        CurrentThreadLocalSimulationSpaceData->ExchangeParticlesBetweenThreadsParallelInsert(StepOutside, *StateOfSimulationSpaceDivisionForThreads, false);
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                        CurrentThreadLocalSimulationSpaceData->ExchangeParticlesBetweenThreadsParallelExtract(StepOutside, *StateOfSimulationSpaceDivisionForThreads, false);

                    SyncPoint->arrive_and_wait();

                    if (CurrentThreadIndexParam == 1)
                        *StateOfSimulationSpaceDivisionForThreads = StepToChangeSimulationSpaceDivisionForThreads(StepOutside, *StateOfSimulationSpaceDivisionForThreads);
                }
                else
                if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::InMainThread && CurrentThreadIndexParam == 1)
                {
                    ExchangeParticlesBetweenThreads(StepOutside, *StateOfSimulationSpaceDivisionForThreads, false);

                    *StateOfSimulationSpaceDivisionForThreads = StepToChangeSimulationSpaceDivisionForThreads(StepOutside, *StateOfSimulationSpaceDivisionForThreads);
                }

                SyncPoint->arrive_and_wait();
            }
            else
                ExchangeParticlesBetweenThreadsAndSychronizeSimulationExecutedInThreads(SyncPoint, CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
        }
    }
    CATCH("generating n steps of simulation for whole cell space in one thread")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInThreads(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = true;

        SetZeroForAllParallelExecutionVariables();

        bool StateOfSimulationSpaceDivisionForThreads = false;

        std::barrier SyncPoint(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation);

        vector<vector<vector<thread*>>> Threads(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<thread*>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<thread*>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        const UnsignedInt NumberOfAllParticlesBeforeSimulation = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL PARTICLES BEFORE SIMULATION = " << NumberOfAllParticlesBeforeSimulation));
        LoggersManagerObject.Log(STREAM("START SIMULATION IN THREADS"));
        LoggersManagerObject.Log(STREAM("MaxParticleUniqueId = " << MaxParticleUniqueId));

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        ThreadIdType ThreadIndex = 1;
        FOR_EACH_THREAD_IN_XYZ
        {
            Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1] = new thread(&CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread, this, &SyncPoint, &StateOfSimulationSpaceDivisionForThreads, NumberOfStepsOutside, NumberOfStepsInside, ThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex, SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]);
            ThreadIndex++;
        }

        FOR_EACH_THREAD_IN_XYZ
        {
            Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->join();
            delete Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1];
        }

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        const string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));

        const UnsignedInt NumberOfAllParticlesAfterSimulation = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        LoggersManagerObject.Log(STREAM("END SIMULATION IN THREADS"));
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL PARTICLES AFTER SIMULATION = " << NumberOfAllParticlesAfterSimulation));
        LoggersManagerObject.Log(STREAM("DIFFERENCE OF ALL PARTICLES BEFORE AND AFTER SIMULATION = " << NumberOfAllParticlesAfterSimulation - NumberOfAllParticlesBeforeSimulation));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsByBarrierPhase1(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        LoggersManagerObject.Log(STREAM("NumberOfActiveNeighbors = " << NumberOfActiveNeighbors));

        int Counter = 0;
        for (const auto& ReceivedParticlesToInsert : CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads)
            if (ReceivedParticlesToInsert.empty() == false)
                Counter++;

        LoggersManagerObject.Log(STREAM("NUMBER OF NEIGHBORS THAT SENT PARTCILES TO ACTUAL THREAD = " << Counter << " ACTUAL THREAD = " << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex));
        #endif

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (CurrentThreadLocalSimulationSpaceData->NeighborProcessesIndexes[NeighborProcessIndex] != -1)
            {
                #if defined(NEIGHBOR_THREAD_INDEXES) && defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION)
                const auto [ThreadPosX, ThreadPosY, ThreadPosZ] = CurrentThreadLocalSimulationSpaceData->NeighborThreadsIndexes[NeighborProcessIndex];
                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("LocalNeighborThreadsIndexes = " << ThreadPosX << " , " << ThreadPosY << " , " << ThreadPosZ << " CurrentThreadIndex = [" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << "]"));)
                #endif

                if (CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].empty() == false)
                {
                    const auto& ReceivedParticlesToInsert = CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex];

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("WANT SENDING CONFIRMATION TO NEIGHBOR = " << ReceivedParticlesToInsert[0].SenderProcessIndex << " " << ReceivedParticlesToInsert[0].ReceiverProcessIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex - 1 << ")"));)

                    #if defined(NEIGHBOR_THREAD_INDEXES) && defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION)
                    LoggersManagerObject.LogUnconditional(STREAM("WANT SENDING CONFIRMATION TO NEIGHBOR = (" << ThreadPosX << " " << ThreadPosY << " " << ThreadPosZ << ") SIZE = " << ReceivedParticlesToInsert.size()));
                    #endif

                    for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                        if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("PARTCLE_KIND_ID TO CHECK = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)

                            CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .SectorPosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .SectorPosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .SectorPosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                            if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                            {
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderProcessIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverProcessIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosZ << ")"));)

                                auto& ParticleObject = CurrentThreadLocalSimulationSpaceData->GetParticleFromUniqueId(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);
                                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                                CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                                bool FoundNeighbor = false;
                                UnsignedInt LocalNeighborProcessIndex = 0;
                                for (LocalNeighborProcessIndex = 0; LocalNeighborProcessIndex < NumberOfAllNeighbors; LocalNeighborProcessIndex++)
                                    if (SimulationSpaceDataForThreads[ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosX - 1][ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosY - 1][ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosZ - 1]->NeighborProcessesIndexes[LocalNeighborProcessIndex] == CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex - 1)
                                    {
                                        FoundNeighbor = true;
                                        break;
                                    }

                                if (FoundNeighbor == true)
                                {
                                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("EXPORT CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ") (" << ReceivedParticleUniqueIdToInsert.Vector.X << "," << ReceivedParticleUniqueIdToInsert.Vector.Y << "," << ReceivedParticleUniqueIdToInsert.Vector.Z << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") PSHIFT = (" << ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X << ", " << ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y << ", " << ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z << ")" << " (" << ReceivedParticleUniqueIdToInsert.SenderProcessIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverProcessIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosZ << ")"));)

                                    SimulationSpaceDataForThreads[ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosX - 1][ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosY - 1][ReceivedParticleUniqueIdToInsert.SenderThreadPos.ThreadPosZ - 1]->ConfirmationOfParticlesToRemoveToSent[LocalNeighborProcessIndex].emplace_back(ConfirmationOfParticlesToRemoveToSentStruct{ .ParticleUniqueId = ReceivedParticleUniqueIdToInsert.ParticleUniqueId, .SenderSectorPos = { .X = ReceivedParticleUniqueIdToInsert.SenderSectorPos.X, .Y = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Y, .Z = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Z }});

                                    CurrentThreadLocalSimulationSpaceData->NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion++;
                                    CurrentThreadLocalSimulationSpaceData->NumberOfParticlesMovedInDiffusion++;
                                }
                                #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                                else
                                {
                                    LoggersManagerObject.LogUnconditional(STREAM("ERROR SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderProcessIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverProcessIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosZ << ")"));
                                }
                                #endif
                            }
                            else
                            {
                                CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread++;

                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderProcessIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverProcessIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosZ << ")"));)
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                            }
                        }

                    CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].clear();
                }
            }
    }
    CATCH("exchange particles between threads ver 2")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsByBarrierPhase2(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (CurrentThreadLocalSimulationSpaceData->NeighborProcessesIndexes[NeighborProcessIndex] != -1)
            {
                for (const auto& [ParticleUniqueId, SenderSectorPos] : CurrentThreadLocalSimulationSpaceData->ConfirmationOfParticlesToRemoveToSent[NeighborProcessIndex])
                {
                    CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .SectorPosX = static_cast<SignedInt>(SenderSectorPos.X), .SectorPosY = static_cast<SignedInt>(SenderSectorPos.Y), .SectorPosZ = static_cast<SignedInt>(SenderSectorPos.Z) });

                    CurrentThreadLocalSimulationSpaceData->RemoveParticle(ParticleUniqueId, true, false);

                    CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread++;

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO REMOVE THIS PARTCILE = " << ParticleUniqueId << " (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")"));)
                }

                CurrentThreadLocalSimulationSpaceData->ConfirmationOfParticlesToRemoveToSent[NeighborProcessIndex].clear();
            }
    }
    CATCH("exchange particles threads processes group 3")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsOneGlobalBarrier(barrier<>* SyncPoint, const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        ExchangeParticlesBetweenThreadsByBarrierPhase1(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);

        SyncPoint->arrive_and_wait();

        ExchangeParticlesBetweenThreadsByBarrierPhase2(CurrentThreadLocalSimulationSpaceData);

        SyncPoint->arrive_and_wait();
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelExecutionManager::SynchronizeWithNeighborByLocalBarrier(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        std::array<std::optional<std::barrier<>::arrival_token>, 6> Tokens;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex > 0)
                Tokens[NeighborProcessIndex] = TwoThreadsWallSychronizationBarriers[WallIndex]->arrive();
        }

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (Tokens[NeighborProcessIndex].has_value() == true)
            {
                const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
                TwoThreadsWallSychronizationBarriers[WallIndex]->wait(std::move(*Tokens[NeighborProcessIndex]));
            }
    }
    CATCH("synchronizing with neighbor by local barrier")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsLocalBarrier(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        ExchangeParticlesBetweenThreadsByBarrierPhase1(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);

        SynchronizeWithNeighborByLocalBarrier(CurrentThreadLocalSimulationSpaceData);

        ExchangeParticlesBetweenThreadsByBarrierPhase2(CurrentThreadLocalSimulationSpaceData);

        SynchronizeWithNeighborByLocalBarrier(CurrentThreadLocalSimulationSpaceData);
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        const UnsignedInt SlotIndex = CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];

            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1, std::defer_lock);
                unique_lock<mutex> LockGuardScopedLock2(Wall.MutexNumber2, std::defer_lock);
                lock(LockGuardScopedLock1, LockGuardScopedLock2);

                Wall.Proposals[SlotIndex] = std::move(CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex]);
                Wall.ProposalsReady[SlotIndex] = true;
            }

            CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].clear();

            Wall.ProposalConditionalVariable.notify_all();
        }

        vector<UniqueIdUInt64Type> ReceivedConfirmationOfParticlesToRemove;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];

            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];
            unique_lock<mutex> LockGuardScopedLock2(Wall.MutexNumber2);

            Wall.VerdictConditionalVariable.wait(LockGuardScopedLock2, [&] { return Wall.VerdictsReady[SlotIndex]; });

            for (const auto& ConfirmationOfParticlesToRemoveToSentObject : Wall.Verdicts[SlotIndex])
                ReceivedConfirmationOfParticlesToRemove.emplace_back(ConfirmationOfParticlesToRemoveToSentObject);

            Wall.VerdictsReady[SlotIndex] = false;
        }

        for (const auto& ParticleToRemoveConfirmedIndex : ReceivedConfirmationOfParticlesToRemove)
            RemoveParticle(ParticleToRemoveConfirmedIndex, true, false);
    }
    CATCH("exchange particles threads processes")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        const UnsignedInt SlotIndex = 1 - CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            std::vector<ParticleSenderStructMultiThreaded> ReceivedParticlesToInsert;
            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1);
                Wall.ProposalConditionalVariable.wait(LockGuardScopedLock1, [&] { return Wall.ProposalsReady[SlotIndex]; });

                ReceivedParticlesToInsert = std::move(Wall.Proposals[SlotIndex]);
                Wall.Proposals[SlotIndex].clear();
                Wall.ProposalsReady[SlotIndex] = false;
            }

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("SENDING CONFIRMATION TO Neighbor = " << ReceivedParticlesToInsert[0].SenderProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

            vector<UniqueIdUInt64Type> ConfirmationOfParticlesToRemoveToSent;
            for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                    if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                        ConfirmationOfParticlesToRemoveToSent.emplace_back(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);

            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1, std::defer_lock);
                unique_lock<mutex> LockGuardScopedLock2(Wall.MutexNumber2, std::defer_lock);

                lock(LockGuardScopedLock1, LockGuardScopedLock2);
                Wall.Verdicts[SlotIndex] = std::move(ConfirmationOfParticlesToRemoveToSent);
                Wall.VerdictsReady[SlotIndex] = true;
            }

            Wall.VerdictConditionalVariable.notify_all();
        }
    }
    CATCH("exchange particles between threads ver 2")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        if (CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber == 0)
        {
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
        }
        else
        {
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
        }
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        const UnsignedInt SlotIndex = CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];

            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            {
                lock_guard<mutex> LockGuardScopedLock1(Wall.MutexNumber1);

                Wall.Proposals[SlotIndex] = std::move(CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex]);
                Wall.ProposalsReady[SlotIndex] = true;
            }

            CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].clear();

            Wall.WallConditionalVariable.notify_all();
        }

        vector<UniqueIdUInt64Type> ReceivedConfirmationOfParticlesToRemove;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1);
            Wall.WallConditionalVariable.wait(LockGuardScopedLock1, [&] { return Wall.VerdictsReady[SlotIndex]; });

            for (const auto& ConfirmationOfParticlesToRemoveToSentObject : Wall.Verdicts[SlotIndex])
                ReceivedConfirmationOfParticlesToRemove.emplace_back(ConfirmationOfParticlesToRemoveToSentObject);

            Wall.Verdicts[SlotIndex].clear();
            Wall.VerdictsReady[SlotIndex] = false;
        }

        for (const auto& ParticleToRemoveConfirmedIndex : ReceivedConfirmationOfParticlesToRemove)
            CurrentThreadLocalSimulationSpaceData->RemoveParticle(ParticleToRemoveConfirmedIndex, true, false);
    }
    CATCH("exchange particles between threads group 1")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        const UnsignedInt SlotIndex = 1 - CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            std::vector<ParticleSenderStructMultiThreaded> ReceivedParticlesToInsert;
            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1);
                Wall.WallConditionalVariable.wait(LockGuardScopedLock1, [&] { return Wall.ProposalsReady[SlotIndex]; });

                ReceivedParticlesToInsert = std::move(Wall.Proposals[SlotIndex]);
                Wall.Proposals[SlotIndex].clear();
                Wall.ProposalsReady[SlotIndex] = false;
            }

            std::vector<UniqueIdUInt64Type> ConfirmationOfParticlesToRemove;
            for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                    if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                        ConfirmationOfParticlesToRemove.emplace_back(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);

            {
                lock_guard<mutex> LockGuardScopedLock1(Wall.MutexNumber1);

                Wall.Verdicts[SlotIndex] = std::move(ConfirmationOfParticlesToRemove);
                Wall.VerdictsReady[SlotIndex] = true;
            }

            Wall.WallConditionalVariable.notify_all();
        }
    }
    CATCH("exchange particles between threads group 2 ver 2")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        if (CurrentThreadLocalSimulationSpaceData->ProcessGroupNumber == 0)
        {
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
        }
        else
        {
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
        }
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                for (UnsignedInt ParticleSectorXIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos; ParticleSectorXIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos; ParticleSectorYIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos; ParticleSectorZIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("XStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << " YStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << " ZStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << " XEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << " YEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << " ZEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, true);
                        }

                MPI_Barrier(MPI_COMM_WORLD);

                for (UnsignedInt ParticleSectorXIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos; ParticleSectorXIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos; ParticleSectorYIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos; ParticleSectorZIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("XStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << " YStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << " ZStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << " XEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << " YEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << " ZEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInMPIProcesses(true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ, CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ);
                        }

                MPI_Barrier(MPI_COMM_WORLD);
            }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam)
{
    try
    {
        for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
        {
            LoggersManagerObject.Log(STREAM("NUMBER OF STEPS OUTSIDE = " << NumberOfStepsInside << " " << StepOutside));

            GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);

            const auto start_time_exchange = chrono::high_resolution_clock::now();

            ExchangeParticlesBetweenMPIProcessesVer2();

            const auto stop_time_exchange = chrono::high_resolution_clock::now();

            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time_exchange, stop_time_exchange, "Only exchanging of particles has taken time","Execution in mpi")));
        }
    }
    CATCH("generating n steps of simulation for whole cell space in one thread")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = true;

        SetZeroForAllParallelExecutionVariables();

        LoggersManagerObject.Log(STREAM("START MPI SIMULATION"));
        LoggersManagerObject.Log(STREAM("MaxParticleUniqueId = " << MaxParticleUniqueId));

        CellEngineUseful::SwitchOffLogs();

        if (CellEngineConfigDataObject.OpenGLGraphicsSwitchedOff == false && MPIProcessDataObject.CurrentMPIProcessIndex == CellEngineConfigDataObject.MainMPIProcessNumber)
        {
            int ValueToSend = 1;
            MPI_Bcast(&ValueToSend, 1, MPI_INT, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);
        }

        MPI_Bcast(&NumberOfStepsOutside, 1, MPI_UINT64_T, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);
        MPI_Bcast(&NumberOfStepsInside, 1, MPI_UINT64_T, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);

        CountAllParticlesInMPI();

        LoggersManagerObject.LogUnconditional(STREAM("NUMBER OF ALL PARTICLES IN MPI SIMULATION = " << NumberOfAllParticlesInMPISimulation << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << "START"));

        const auto start_time = chrono::high_resolution_clock::now();

        GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsOutside, NumberOfStepsInside, MPIProcessDataObject.CurrentMPIProcessIndex, 0, 0, 0);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        const string ResultText = "Execution in MPI processes for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in mpi")));

        LoggersManagerObject.Log(STREAM("END MPI SIMULATION"));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in mpi process")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup1()
{
    CellEngineUseful::SwitchOnLogs();

    try
    {
        MainMapType<UniqueIdUInt64Type, SectorPosType> AllParticlesSentToAnotherProcess;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (NeighborProcessesIndexes[NeighborProcessIndex] != -1)
            {
                if (VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].empty() == true)
                {
                    VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].emplace_back(ParticleSenderStructMPIMultiProcess{ 0, 0, static_cast<int>(MPIProcessDataObject.CurrentMPIProcessIndex), static_cast<int>(NeighborProcessesIndexes[NeighborProcessIndex]), 0, 0, 0, 0, 0, 0 });

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Index to EMPTY send = " << NeighborProcessesIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex][0].ReceiverProcessIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)
                }
                else
                {
                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Index to send = " << NeighborProcessesIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex][0].ReceiverProcessIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)
                }

                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Length Message SEND = " << VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].size() << " " << VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].size() * sizeof(ParticleSenderStructMPIMultiProcess) << " " << NeighborProcessesIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex][0].ReceiverProcessIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                char BufferToSend[MaxMPIMessageSize];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].size();

                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI PACKED SIZE = " << NumberOfPackedStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);

                for (const auto& ParticleToSendElement : VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex])
                {
                    MPI_Pack(&ParticleToSendElement.ParticleUniqueId, 1, MPI_UINT64_T, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ParticleKindId, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SenderProcessIndex, 1, MPI_INT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverProcessIndex, 1, MPI_INT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverSectorPos.X, 1, MPI_UNSIGNED_SHORT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverSectorPos.Y, 1, MPI_UNSIGNED_SHORT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverSectorPos.Z, 1, MPI_UNSIGNED_SHORT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewParticleCenter.X, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewParticleCenter.Y, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewParticleCenter.Z, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.Vector.X, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.Vector.Y, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.Vector.Z, 1, MPI_FLOAT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);

                    AllParticlesSentToAnotherProcess.insert(make_pair(ParticleToSendElement.ParticleUniqueId, SectorPosType(ParticleToSendElement.SenderSectorPos.X, ParticleToSendElement.SenderSectorPos.Y, ParticleToSendElement.SenderSectorPos.Z)));
                }

                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, static_cast<int>(NeighborProcessesIndexes[NeighborProcessIndex]), 0, MPI_COMM_WORLD);

                VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].clear();
            }

        int NumberOfReceivedMessages = 0;
        while (NumberOfReceivedMessages < NumberOfActiveNeighbors)
        {
            MPI_Status MPIMessageStatus;

            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int NumberOfBytesReceived;
            MPI_Get_count(&MPIMessageStatus, MPI_UNSIGNED, &NumberOfBytesReceived);

            vector<UniqueIdUInt64Type> ReceivedConfirmationOfParticlesToRemove;

            char ReceivedParticlesToInsert1[MaxMPIMessageSize];
            MPI_Recv(&ReceivedParticlesToInsert1, MaxMPIMessageSize, MPI_PACKED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int PositionInBuffer = 0;
            unsigned int NumberOfUnpackedParticleStructures;
            MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &NumberOfUnpackedParticleStructures, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI UNPACKED SIZE TO REMOVE = " << NumberOfUnpackedParticleStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

            for (UnsignedInt UnpackedParticleStructureIndex = 0; UnpackedParticleStructureIndex < NumberOfUnpackedParticleStructures; UnpackedParticleStructureIndex++)
            {
                UniqueIdUInt64Type MPIParticleSenderStructElementLocalObject;
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject, 1, MPI_UINT64_T, MPI_COMM_WORLD);
                ReceivedConfirmationOfParticlesToRemove.emplace_back(MPIParticleSenderStructElementLocalObject);
            }

            if (ReceivedConfirmationOfParticlesToRemove[0] != 0)
                for (const auto& ParticleToRemoveConfirmedIndex : ReceivedConfirmationOfParticlesToRemove)
                {
                    SetCurrentSectorPos(AllParticlesSentToAnotherProcess.find(ParticleToRemoveConfirmedIndex)->second);

                    RemoveParticle(ParticleToRemoveConfirmedIndex, true, false);

                    NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread++;
                }

            NumberOfReceivedMessages++;
        }
    }
    CATCH("exchange particles between mpi processes")

    CellEngineUseful::SwitchOffLogs();
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup2()
{
    CellEngineUseful::SwitchOnLogs();

    try
    {
        int NumberOfReceivedMessages = 0;

        vector<ParticleSenderStructMPIMultiProcess> ReceivedParticlesToInsertFromAllNeigbhours[NumberOfAllNeighbors];

        while (NumberOfReceivedMessages < NumberOfActiveNeighbors)
        {
            MPI_Status MPIMessageStatus;

            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int NumberOfBytesReceived;
            MPI_Get_count(&MPIMessageStatus, MPI_CHAR, &NumberOfBytesReceived);

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Length Message RECEIVE = " << NumberOfBytesReceived << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << NumberOfBytesReceived / sizeof(ParticleSenderStructMPIMultiProcess)));)

            char ReceivedParticlesToInsert1[MaxMPIMessageSize];

            #ifdef SIMULATION_DETAILED_DE
            BUG_LOG_DIFFUSION
            if (NumberOfBytesReceived > MaxMPIMessageSize)
                LoggersManagerObject.Log(STREAM("ERROR NumberOfBytesReceived = " << NumberOfBytesReceived << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
            #endif

            MPI_Recv(&ReceivedParticlesToInsert1, MaxMPIMessageSize, MPI_PACKED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int PositionInBuffer = 0;
            unsigned int NumberOfUnpackedParticleStructures;
            MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &NumberOfUnpackedParticleStructures, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI UNPACKED SIZE = " << NumberOfUnpackedParticleStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

            for (UnsignedInt UnpackedParticleStructureIndex = 0; UnpackedParticleStructureIndex < NumberOfUnpackedParticleStructures; UnpackedParticleStructureIndex++)
            {
                ParticleSenderStructMPIMultiProcess MPIParticleSenderStructElementLocalObject;
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleUniqueId, 1, MPI_UINT64_T, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleKindId, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SenderProcessIndex, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverProcessIndex, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverSectorPos.X, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverSectorPos.Y, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverSectorPos.Z, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewParticleCenter.X, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewParticleCenter.Y, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewParticleCenter.Z, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.Vector.X, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.Vector.Y, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.Vector.Z, 1, MPI_FLOAT, MPI_COMM_WORLD);

                bool Found = false;
                UnsignedInt NeighborProcessIndex = 0;
                for (NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                    if (NeighborProcessesIndexes[NeighborProcessIndex] == MPIParticleSenderStructElementLocalObject.SenderProcessIndex)
                    {
                        ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex].emplace_back(MPIParticleSenderStructElementLocalObject);
                        Found = true;
                        break;
                    }

                #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                if (Found == false)
                    LoggersManagerObject.Log(STREAM("GET PROCESS INDEX FROM = " << NeighborProcessIndex << " " << NeighborProcessesIndexes[NeighborProcessIndex] << " " << MPIParticleSenderStructElementLocalObject.SenderProcessIndex << " TO = " << MPIParticleSenderStructElementLocalObject.ReceiverProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex << endl));
                #endif
            }

            NumberOfReceivedMessages++;
        }

        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        int Counter = 0;
        for (const auto& ReceivedParticlesToInsert : ReceivedParticlesToInsertFromAllNeigbhours)
            if (ReceivedParticlesToInsert.empty() == false)
                Counter++;
        LoggersManagerObject.Log(STREAM("FROM NeighborS = " << Counter << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
        #endif

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (NeighborProcessesIndexes[NeighborProcessIndex] != -1)
            {
                vector<UniqueIdUInt64Type> ConfirmationOfParticlesToRemoveToSent;

                if (ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex].empty() == false)
                {
                    const auto ReceivedParticlesToInsert = ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex];

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("SENDING CONFIRMATION TO Neighbor = " << ReceivedParticlesToInsert[0].SenderProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                    for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                        if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                        {
                            SetCurrentSectorPos(SectorPosType{ .SectorPosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .SectorPosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .SectorPosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                            UniqueIdUInt64Type UniqueIdOfNewCreatedParticle;
                            if (CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMPIMultiProcessDiffusion(ReceivedParticleUniqueIdToInsert, UniqueIdOfNewCreatedParticle) == true)
                            {
                                auto& ParticleObject = GetParticleFromUniqueId(UniqueIdOfNewCreatedParticle);
                                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                                CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                                ConfirmationOfParticlesToRemoveToSent.emplace_back(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);

                                NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion++;
                                NumberOfParticlesMovedInDiffusion++;
                            }
                        }
                }
                else
                {
                    NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread++;

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderProcessIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverProcessIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.ThreadPosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverThreadPos.ThreadPosZ << ")"));)
                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                }

                if (ConfirmationOfParticlesToRemoveToSent.empty() == true)
                    ConfirmationOfParticlesToRemoveToSent.emplace_back(0);

                char BufferToSend[MaxMPIMessageSize];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = ConfirmationOfParticlesToRemoveToSent.size();

                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                for (const auto& ParticleToSendElement : ConfirmationOfParticlesToRemoveToSent)
                    MPI_Pack(&ParticleToSendElement, 1, MPI_UINT64_T, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, NeighborProcessesIndexes[NeighborProcessIndex], 0, MPI_COMM_WORLD);
            }
    }
    CATCH("exchange particles between mpi processes ver 2")

    CellEngineUseful::SwitchOffLogs();
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesVer2()
{
    try
    {
        if (ProcessGroupNumber == 0)
        {
            ExchangeParticlesBetweenMPIProcessesGroup1();

            MPI_Barrier(MPI_COMM_WORLD);

            ExchangeParticlesBetweenMPIProcessesGroup2();

            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            ExchangeParticlesBetweenMPIProcessesGroup2();

            MPI_Barrier(MPI_COMM_WORLD);

            ExchangeParticlesBetweenMPIProcessesGroup1();

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    CATCH("exchange particles between mpi processes")
}

void CellEngineSimulationParallelExecutionManager::SetZeroForAllParallelExecutionVariables()
{
    try
    {
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL PARTICLES BEFORE SIMULATION = " << NumberOfAllParticlesInMPISimulation));

        NumberOfDNAParticles = 0;

        ErrorCounter = 0;
        NumberOfExecutedReactions = 0;
        NumberOfCancelledReactions = 0;
        NumberOfCancelledAReactions = 0;
        NumberOfCancelledBReactions = 0;
        AddedParticlesInReactions = 0;
        RemovedParticlesInReactions = 0;
        RestoredParticlesInCancelledReactions = 0;

        NumberOfParticlesMovedInDiffusion = 0;
        NumberOfParticlesToBeMovedInDiffusionRejectedLocally = 0;
        NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread = 0;
        NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion = 0;
        NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread  = 0;
        NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfDNAParticles = 0;

            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfExecutedReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledAReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledBReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RemovedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RestoredParticlesInCancelledReactions = 0;

            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesMovedInDiffusion = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionRejectedLocally = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread = 0;
        }
    }
    CATCH("setting zero for all parallel execution variables")
}

void CellEngineSimulationParallelExecutionManager::PrintAllExecutionVariables(const UnsignedInt PhaseId) const
{
    try
    {
        LoggersManagerObject.Log(STREAM("DIFFUSION STATISTICS:"));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesMovedInDiffusion = " + to_string(NumberOfParticlesMovedInDiffusion) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionRejectedLocally = " + to_string(NumberOfParticlesToBeMovedInDiffusionRejectedLocally) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread = " + to_string(NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread = " + to_string(NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion = " + to_string(NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread = " + to_string(NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));

        LoggersManagerObject.Log(STREAM("CHEMICAL REACTIONS STATISTICS:"));
        LoggersManagerObject.Log(STREAM("RestoredParticlesInCancelledReactions = " + to_string(RestoredParticlesInCancelledReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("RemovedParticlesInReactions = " + to_string(RemovedParticlesInReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("AddedParticlesInReactions = " + to_string(AddedParticlesInReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledAReactions = " + to_string(NumberOfCancelledAReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledBReactions = " + to_string(NumberOfCancelledBReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledReactions = " + to_string(NumberOfCancelledReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfExecutedReactions = " + to_string(NumberOfExecutedReactions) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("ErrorCounter = " + to_string(ErrorCounter) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));

        LoggersManagerObject.Log(STREAM("GENERAL STATISTICS:"));
        LoggersManagerObject.Log(STREAM("NumberOfDNAParticles = " << NumberOfDNAParticles << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfAllParticlesInMPISimulation = " << NumberOfAllParticlesInMPISimulation << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
    }
    CATCH("")
}

void CellEngineSimulationParallelExecutionManager::CountAllParticlesInMPI()
{
    try
    {
        NumberOfAllParticlesInMPISimulation = 0;

        for (UnsignedInt ParticleSectorXIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos; ParticleSectorXIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos; ParticleSectorXIndex++)
            for (UnsignedInt ParticleSectorYIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos; ParticleSectorYIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos; ParticleSectorYIndex++)
                for (UnsignedInt ParticleSectorZIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos; ParticleSectorZIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos; ParticleSectorZIndex++)
                    NumberOfAllParticlesInMPISimulation += CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(SectorPosType(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex));

        MPI_Allreduce(MPI_IN_PLACE, &NumberOfAllParticlesInMPISimulation, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    }
    CATCH("counting all particles in MPI")
}

void CellEngineSimulationParallelExecutionManager::GatherAllParallelExecutionVariables()
{
    try
    {
        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true)
        {
            PrintAllExecutionVariables(1);

            MPI_Allreduce(MPI_IN_PLACE, &NumberOfDNAParticles, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

            CountAllParticlesInMPI();

            MPI_Allreduce(MPI_IN_PLACE, &ErrorCounter, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfExecutedReactions, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfCancelledReactions, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfCancelledAReactions , 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfCancelledBReactions , 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &AddedParticlesInReactions, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &RemovedParticlesInReactions, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &RestoredParticlesInCancelledReactions, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesMovedInDiffusion, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionRejectedLocally, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        }
        else
        if (CellEngineConfigDataObject.MultiThreaded == true)
            FOR_EACH_THREAD_IN_XYZ
            {
                NumberOfDNAParticles += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfDNAParticles;

                ErrorCounter += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter;
                NumberOfExecutedReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfExecutedReactions;
                NumberOfCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledReactions;
                NumberOfCancelledAReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledAReactions;
                NumberOfCancelledBReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledBReactions;
                AddedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions;
                RemovedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RemovedParticlesInReactions;
                RestoredParticlesInCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RestoredParticlesInCancelledReactions;

                NumberOfParticlesMovedInDiffusion += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesMovedInDiffusion;
                NumberOfParticlesToBeMovedInDiffusionRejectedLocally += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionRejectedLocally;
                NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherThread;
                NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread  += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherThread;
                NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesSentToAnotherThreadAndMovedInDiffusion;
                NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherThread;
            }

        PrintAllExecutionVariables(2);
    }
    CATCH("gathering all parallel execution variables")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, bool PrintTime)
{
    try
    {
        FirstSendParticlesForThreads(false, true);
        GenerateNStepsOfSimulationForWholeCellSpaceInThreads(NumberOfStepsOutside, NumberOfStepsInside);
        GatherParticlesFromThreadsToParticlesInMainThread();
    }
    CATCH("generate n steps of simulation with sending particles to threads and gathering particles to main threads for whole cell space")
}

void CellEngineSimulationParallelExecutionManager::SaveFormerParticlesAsVectorElements()
{
    try
    {
        FormerParticlesIndexes = GetParticles();
    }
    CATCH("saving former particles as vector elements")
}

void CellEngineSimulationParallelExecutionManager::CheckParticlesCenters(const bool PrintAllParticles)
{
    try
    {
        UnsignedInt NumberOfZeroSizedParticles = 0;
        UnsignedInt NumberOfBadCenterParticles = 0;

        FOR_EACH_PARTICLE_IN_SECTORS_XYZ_CONST
        {
            if (PrintAllParticles ==  true)
                LoggersManagerObject.Log(STREAM("ParticleUniqueId = " << ParticleObject.second.ParticleUniqueId << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));

            if (ParticleObject.second.Center.X > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Y > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Z > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                LoggersManagerObject.Log(STREAM("Wrong Particle Center -> ParticleUniqueId = " << ParticleObject.first << " " << ParticleObject.second.ParticleUniqueId << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z  << " ParticleKindId = " << ParticleObject.second.ParticleKindId));

            if (ParticleObject.second.Center.X == 0 || ParticleObject.second.Center.Y == 0 || ParticleObject.second.Center.Z == 0)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ZERO -> ParticleUniqueId = " << ParticleObject.first << " " << ParticleObject.second.ParticleUniqueId << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z << " ParticleKindId = " << ParticleObject.second.ParticleKindId));
                NumberOfBadCenterParticles++;
            }
            if (ParticleObject.second.Center.X == 1 || ParticleObject.second.Center.Y == 1 || ParticleObject.second.Center.Z == 1)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ONE -> ParticleUniqueId = " << ParticleObject.first << " " << ParticleObject.second.ParticleUniqueId << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z << " ParticleKindId = " << ParticleObject.second.ParticleKindId));
                NumberOfBadCenterParticles++;
            }
            if (ParticleObject.second.ListOfVoxels.empty() == true)
                NumberOfZeroSizedParticles++;
        }

        LoggersManagerObject.Log(STREAM("Number of erased particles with list of voxels sized zero = " << erase_if(GetParticles(), [](auto& P){ return P.second.ListOfVoxels.empty() == true; })));
        LoggersManagerObject.Log(STREAM("Number of erased particles with bad center = " << erase_if(GetParticles(), [](auto& P){ return P.second.Center.X == 0 || P.second.Center.Y == 0 || P.second.Center.Z == 0; })));

        LoggersManagerObject.Log(STREAM("All Particles Centers Checked. Number Of Bad Center Particles = " << NumberOfBadCenterParticles));
        LoggersManagerObject.Log(STREAM("All Particles Size Checked. Number Of Zero Sized Particles = " << NumberOfZeroSizedParticles));

        LoggersManagerObject.LogError(STREAM("Number of erased particles with bad center = " << erase_if(GetParticles(), [](auto& P){ return P.second.Center.X == 0 || P.second.Center.Y == 0 || P.second.Center.Z == 0; })));
        LoggersManagerObject.LogError(STREAM("Number of erased particles with list of voxels sized zero = " << erase_if(GetParticles(), [](auto& P){ return P.second.ListOfVoxels.empty() == true; })));

        LoggersManagerObject.Log(STREAM(""));
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}
