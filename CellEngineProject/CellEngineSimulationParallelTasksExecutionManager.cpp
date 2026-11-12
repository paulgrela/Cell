
#include <mpi.h>
#include <chrono>

#include "DateTimeUtils.h"

#include "CellEngineDataFile.h"
#include "CellEngineImGuiMenu.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineThreadPinningToRealCPUCore.h"
#include "CellEngineSimulationParallelTasksExecutionManager.h"

#ifdef USE_OPENGL
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"
#endif

constexpr bool PrintAdditionalInformation = true;

using namespace std;

CellEngineSimulationParallelTasksExecutionManager::CellEngineSimulationParallelTasksExecutionManager() : SimulationSpaceDataForThreads(CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer)
{
}

template void CellEngineSimulationParallelTasksExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineVoxelSimulationSpace>(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);
template void CellEngineSimulationParallelTasksExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineFullAtomSimulationSpace>(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles);

template <class SimulationSpaceType>
void CellEngineSimulationParallelTasksExecutionManager::CreateSimulationSpaceForParallelExecution(SimulationSpaceForParallelExecutionContainer<CellEngineSimulationSpace>& CellEngineSimulationSpaceForThreadsObjectsPointer, ParticlesContainer<Particle>& Particles)
{
    try
    {
        UnsignedInt ThreadIndexPos = 0;

        CellEngineSimulationSpaceForThreadsObjectsPointer.clear();
        CellEngineSimulationSpaceForThreadsObjectsPointer.resize(CellEngineConfigDataObject.NumberOfXThreadsInSimulation);
        UnsignedInt ThreadXPos = 0;
        for (auto& ThreadLocalParticlesInProximityXPos : CellEngineSimulationSpaceForThreadsObjectsPointer)
        {
            ThreadLocalParticlesInProximityXPos.resize(CellEngineConfigDataObject.NumberOfYThreadsInSimulation);
            UnsignedInt ThreadYPos = 0;
            for (auto& ThreadLocalParticlesInProximityYPos : ThreadLocalParticlesInProximityXPos)
            {
                ThreadLocalParticlesInProximityYPos.resize(CellEngineConfigDataObject.NumberOfZThreadsInSimulation);
                UnsignedInt ThreadZPos = 0;
                for (auto& ThreadLocalParticlesInProximityZPos : ThreadLocalParticlesInProximityYPos)
                {
                    LoggersManagerObject.Log(STREAM("THREAD INDEXES = " << ThreadIndexPos << " (" << ThreadXPos << ", " << ThreadYPos << ", " << ThreadZPos << ")"));

                    ThreadLocalParticlesInProximityZPos = make_shared<SimulationSpaceType>(Particles, false, ThreadIndexPos, ParallelTaskPosType{ .PosX = static_cast<SignedInt>(ThreadXPos), .PosY = static_cast<SignedInt>(ThreadYPos), .PosZ = static_cast<SignedInt>(ThreadZPos) });
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

void CellEngineSimulationParallelTasksExecutionManager::SetProperParallelTasksPositionsAndProperMPIProcessIndexesForEveryParticlesSector(ParticlesContainer<Particle>& ParticlesSectors)
{
    try
    {
        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            const ParallelTaskPosType ParallelTaskPos = { .PosX = static_cast<SignedInt>(ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation), .PosY = static_cast<SignedInt>(ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation), .PosZ = static_cast<SignedInt>(ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation) };

            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].ParallelTaskPos = ParallelTaskPos;

            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ParallelTaskPos.PosX][ParallelTaskPos.PosY][ParallelTaskPos.PosZ]->MPIProcessIndex;

            #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
            if constexpr (PrintAdditionalInformation == true)
                LoggersManagerObject.Log(STREAM("ThreadPos = " << ParallelTaskPos.PosX << "'" << ParallelTaskPos.PosY << "'" << ParallelTaskPos.PosZ << " ProcessIndex = " << ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex));
            #endif
        }
    }
    CATCH("setting proper thread index for every particles sector")
}

static UnsignedInt GetParallelTaskGroupNumber(const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    if (ThreadZIndex % 2 == 0)
        return ThreadXIndex % 2 != ThreadYIndex % 2;
    else
        return ThreadXIndex % 2 == ThreadYIndex % 2;
}

SignedInt CellEngineSimulationParallelTasksExecutionManager::GetProcessPrevNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (ThreadXIndex >= 0 && ThreadYIndex >= 0 && ThreadZIndex >= 0)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex;
    else
        return -1;
}

SignedInt CellEngineSimulationParallelTasksExecutionManager::GetProcessNextNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (static_cast<SignedInt>(ThreadXIndex) < CellEngineConfigDataObject.NumberOfXThreadsInSimulation && static_cast<SignedInt>(ThreadYIndex) < CellEngineConfigDataObject.NumberOfYThreadsInSimulation && static_cast<SignedInt>(ThreadZIndex) < CellEngineConfigDataObject.NumberOfZThreadsInSimulation)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex;
    else
        return -1;
}

void CellEngineSimulationParallelTasksExecutionManager::CreateDataForEveryMPIProcessForParallelExecution()
{
    try
    {
        SignedInt MPIProcessIndex = 0;

        for (UnsignedInt MPIProcessXIndex = 0; MPIProcessXIndex < CellEngineConfigDataObject.NumberOfXThreadsInSimulation; MPIProcessXIndex++)
            for (UnsignedInt MPIProcessYIndex = 0; MPIProcessYIndex < CellEngineConfigDataObject.NumberOfYThreadsInSimulation; MPIProcessYIndex++)
                for (UnsignedInt MPIProcessZIndex = 0; MPIProcessZIndex < CellEngineConfigDataObject.NumberOfZThreadsInSimulation; MPIProcessZIndex++)
                {
                    if (MPIProcessIndex == MPIProcessDataObject.CurrentMPIProcessIndex)
                    {
                        SimulationSpaceSectorsRangesForCurrentParallelTask.SetParameters(MPIProcessXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, MPIProcessYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, MPIProcessZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation, (MPIProcessXIndex + 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, (MPIProcessYIndex + 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, (MPIProcessZIndex + 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation);

                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Bounds = " << MPIProcessDataObject.CurrentMPIProcessIndex << " (" << SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos << ") (" << SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos << ")"));

                        ParallelTaskGroupNumber = GetParallelTaskGroupNumber(MPIProcessXIndex, MPIProcessYIndex, MPIProcessZIndex);

                        LoggersManagerObject.Log(STREAM("ParallelTaskGroupNumber = " << ParallelTaskGroupNumber << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " (" << MPIProcessXIndex << "," << MPIProcessYIndex << "," << MPIProcessZIndex << ")"));

                        NeighborParallelTasksIndexes[0] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex));
                        NeighborParallelTasksIndexes[1] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex));
                        NeighborParallelTasksIndexes[2] = GetProcessPrevNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex) - 1);

                        NeighborParallelTasksIndexes[3] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex) + 1, static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex));
                        NeighborParallelTasksIndexes[4] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex) + 1, static_cast<SignedInt>(MPIProcessZIndex));
                        NeighborParallelTasksIndexes[5] = GetProcessNextNeighbor(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex) + 1);

                        ListOfParticlesToChangeSectors.reserve(65536);

                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                            if (NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
                            {
                                VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].reserve(16384);
                                VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].reserve(16384);
                            }

                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Neighbors = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " " << NeighborParallelTasksIndexes[0] << " " << NeighborParallelTasksIndexes[1] << " " << NeighborParallelTasksIndexes[2] << " " << NeighborParallelTasksIndexes[3] << " " << NeighborParallelTasksIndexes[4] << " " << NeighborParallelTasksIndexes[5]));
                    }

                    MPIProcessIndex++;
                }

        NumberOfActiveNeighbors = count_if(NeighborParallelTasksIndexes, NeighborParallelTasksIndexes + NumberOfAllNeighbors, [](const auto& Element){ return Element != -1; });

        LoggersManagerObject.Log(STREAM("NumberOfActiveNeighbors = " << NumberOfActiveNeighbors));
    }
    CATCH("creating data for every mpi process for parallel execution")
}

constexpr SignedInt NeighborOffsetX[NumberOfAllNeighbors] = { -1,  0,  0, +1,  0,  0 };
constexpr SignedInt NeighborOffsetY[NumberOfAllNeighbors] = {  0, -1,  0,  0, +1,  0 };
constexpr SignedInt NeighborOffsetZ[NumberOfAllNeighbors] = {  0,  0, -1,  0,  0, +1 };

SignedInt CellEngineSimulationParallelTasksExecutionManager::GetProcessNeighbor(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (ThreadXIndex >= 0 && ThreadXIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfXThreadsInSimulation) && ThreadYIndex >= 0 && ThreadYIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation) && ThreadZIndex >= 0 && ThreadZIndex < static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation))
        return static_cast<SignedInt>(SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->CurrentThreadIndex);
    else
        return -1;
}

void CellEngineSimulationParallelTasksExecutionManager::CreateDataForEveryThreadForParallelExecution()
{
    try
    {
        LoggersManagerObject.Log(STREAM("STEP 1"));

        SignedInt ThreadIndex = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            const auto& ThreadData = SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex];

            ThreadData->SimulationSpaceSectorsRangesForCurrentParallelTask.SetParameters(ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation, (ThreadXIndex + 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation, (ThreadYIndex + 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation, (ThreadZIndex + 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation);

            ThreadData->ParallelTaskGroupNumber = GetParallelTaskGroupNumber(ThreadXIndex, ThreadYIndex, ThreadZIndex);

            for (UnsignedInt DirectionIndex = 0; DirectionIndex < NumberOfAllNeighbors; DirectionIndex++)
                ThreadData->TwoThreadsWallSychronizationBarriersIndexes[DirectionIndex] = -1;

            if constexpr (PrintAdditionalInformation == true)
            {
                LoggersManagerObject.Log(STREAM("Current Thread Index = " << ThreadData->CurrentThreadIndex << " Thread Bounds = (" << SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos << ") (" << SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos << ")"));
                LoggersManagerObject.Log(STREAM("ParallelTaskGroupNumber = " << ThreadData->ParallelTaskGroupNumber << " Thread Position = (" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << ")"));
            }

            ThreadIndex++;
        }

        LoggersManagerObject.Log(STREAM("STEP 2"));

        FOR_EACH_THREAD_IN_XYZ
        {
            const auto& ThreadData = SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex];

            ThreadData->NumberOfActiveNeighbors = 0;

            for (UnsignedInt DirectionIndex = 0; DirectionIndex < NumberOfAllNeighbors; DirectionIndex++)
            {
                const SignedInt NeighborXIndex = ThreadXIndex + NeighborOffsetX[DirectionIndex];
                const SignedInt NeighborYIndex = ThreadYIndex + NeighborOffsetY[DirectionIndex];
                const SignedInt NeighborZIndex = ThreadZIndex + NeighborOffsetZ[DirectionIndex];

                ThreadData->NeighborParallelTasksIndexes[DirectionIndex] = GetProcessNeighbor(NeighborXIndex, NeighborYIndex, NeighborZIndex);
                #ifdef NEIGHBOR_PARALLEL_TASKS_POSITIONS
                ThreadData->NeighborParallelTasksPositions[DirectionIndex] = ParallelTaskPosType{ .PosX = NeighborXIndex, .PosY = NeighborYIndex, .PosZ = NeighborZIndex };
                #endif

                if (ThreadData->NeighborParallelTasksIndexes[DirectionIndex] != -1)
                    ThreadData->NumberOfActiveNeighbors++;
            }

            if constexpr (PrintAdditionalInformation == true)
            {
                const auto NPI = ThreadData->NeighborParallelTasksIndexes;
                LoggersManagerObject.Log(STREAM("Thread Neighbors = " << ThreadData->CurrentThreadIndex << " " << NPI[0] << " " << NPI[1] << " " << NPI[2] << " " << NPI[3] << " " << NPI[4] << " " << NPI[5]));
                #ifdef NEIGHBOR_PARALLEL_TASKS_POSITIONS
                const auto NTI = ThreadData->NeighborParallelTasksPositions;
                LoggersManagerObject.Log(STREAM("Thread Neighbors Indexes = " << ThreadData->CurrentThreadIndex << " (" << NTI[0].PosX << " " << NTI[0].PosY << " " << NTI[0].PosZ << ") (" << NTI[1].PosX << " " << NTI[1].PosY << " " << NTI[1].PosZ << ") (" << NTI[2].PosX << " " << NTI[2].PosY << " " << NTI[2].PosZ << ") (" << NTI[3].PosX << " " << NTI[3].PosY << " " << NTI[3].PosZ << ") (" << NTI[4].PosX << " " << NTI[4].PosY << " " << NTI[4].PosZ << ") (" << NTI[5].PosX << " " << NTI[5].PosY << " " << NTI[5].PosZ << ")"));
                #endif
            }

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ListOfParticlesToChangeSectors.reserve(65536);

            for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                if (ThreadData->NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
                {
                    ThreadData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].reserve(16384);
                    ThreadData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].reserve(16384);
                }
        }

        LoggersManagerObject.Log(STREAM("STEP 3"));

        TwoThreadsWallSychronizationBarriers.clear();
        UnsignedInt TwoThreadsWallSychronizationBarriersNumber = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            if (ThreadXIndex + 1 < CellEngineConfigDataObject.NumberOfXThreadsInSimulation)
            {
                TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));

                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->TwoThreadsWallSychronizationBarriersIndexes[3] = TwoThreadsWallSychronizationBarriersNumber;
                SimulationSpaceDataForThreads[ThreadXIndex + 1][ThreadYIndex][ThreadZIndex]->TwoThreadsWallSychronizationBarriersIndexes[0] = TwoThreadsWallSychronizationBarriersNumber;

                LoggersManagerObject.Log(STREAM("THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << ") with CUBE (" << ThreadXIndex + 1 << "," << ThreadYIndex << "," << ThreadZIndex << ")"));
                ++TwoThreadsWallSychronizationBarriersNumber;
            }
            if (ThreadYIndex + 1 < CellEngineConfigDataObject.NumberOfYThreadsInSimulation)
            {
                TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));

                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->TwoThreadsWallSychronizationBarriersIndexes[4] = TwoThreadsWallSychronizationBarriersNumber;
                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex + 1][ThreadZIndex]->TwoThreadsWallSychronizationBarriersIndexes[1] = TwoThreadsWallSychronizationBarriersNumber;

                LoggersManagerObject.Log(STREAM("THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << ") with CUBE (" << ThreadXIndex << "," << ThreadYIndex + 1 << "," << ThreadZIndex << ")"));
                ++TwoThreadsWallSychronizationBarriersNumber;
            }
            if (ThreadZIndex + 1 < CellEngineConfigDataObject.NumberOfZThreadsInSimulation)
            {
                TwoThreadsWallChannels.emplace_back(make_unique<WallChannel>());
                TwoThreadsWallSychronizationBarriers.emplace_back(make_unique<barrier<>>(2));

                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->TwoThreadsWallSychronizationBarriersIndexes[5] = TwoThreadsWallSychronizationBarriersNumber;
                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex + 1]->TwoThreadsWallSychronizationBarriersIndexes[2] = TwoThreadsWallSychronizationBarriersNumber;

                LoggersManagerObject.Log(STREAM("THE WALL " << TwoThreadsWallSychronizationBarriersNumber << " CONNECTS CUBE (" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << ") with CUBE (" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex + 1 << ")"));
                ++TwoThreadsWallSychronizationBarriersNumber;
            }
        }

        LoggersManagerObject.Log(STREAM("WALLS NUMBER 2 = " << TwoThreadsWallSychronizationBarriersNumber));
    }
    CATCH("creating data for every mpi process for parallel execution")
}

#ifdef SHORTER_CODE
void CellEngineSimulationParallelTasksExecutionManager::JoinReactionsStatisticsFromThreads(vector<map<UnsignedInt, ChemicalReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
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
                for (const auto& [ReactionStatisticsDataFirst, ReactionStatisticsDataSecond] : SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->SavedReactionsMap.back())
                    SavedReactionsMap[SimulationStepNumber - 1][ReactionStatisticsDataFirst].Counter += ReactionStatisticsDataSecond.Counter;
    }
    CATCH("joining statistics from threads")
}
#else
void CellEngineSimulationParallelTasksExecutionManager::JoinStatisticsFromThreads(vector<map<UnsignedInt, ChemicalReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ReactionData : SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->SavedReactionsMap[SimulationStepNumber - 1])
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

void CellEngineSimulationParallelTasksExecutionManager::FirstSendParticlesForThreads(const bool PrintCenterOfParticleWithThreadIndex, const bool PrintTime)
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("THREAD = " << CurrentThreadIndex << " X = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInX << " Y = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInY << " Z = " << CellEngineConfigDataObject.NumberOfParticlesSectorsInZ));

        UnsignedInt GoodParticlesCounter = 0;
        UnsignedInt BadParticlesCounter = 0;

        SaveFormerParticlesAsVectorElements();

        const auto start_time = chrono::high_resolution_clock::now();

        FOR_EACH_THREAD_IN_XYZ
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.clear();

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
                LoggersManagerObject.Log(STREAM("THREAD[" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << "] SIZE = " << SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.size()));

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

vector<vector<vector<ParticlesContainerInternal<Particle>>>> CellEngineSimulationParallelTasksExecutionManager::GatherParticlesToExchangeBetweenThreads(const UnsignedInt TypeOfGet, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, UnsignedInt& ExchangedParticleCounter, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<ParticlesContainerInternal<Particle>>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<ParticlesContainerInternal<Particle>>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

    try
    {
        auto ParticleIter = SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.begin();

        while (ParticleIter != SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.end())
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

            if ((ThreadXIndex != ThreadXIndexNew) || (ThreadYIndex != ThreadYIndexNew) || (ThreadZIndex != ThreadZIndexNew))
            {
                switch (TypeOfGet)
                {
                    case 0 : SimulationSpaceDataForThreads[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew]->ParticlesForThreads.insert(SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.extract(ParticleIter++)); break;
                    case 1 : ParticlesToExchange[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew].Particles.insert(SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.extract(ParticleIter++)); break;
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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreads(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsParallelInsert(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCounter = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange;

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.PosX][CurrentThreadPos.PosY][CurrentThreadPos.PosZ]->MainExchangeParticlesMutexObject };

            ParticlesToExchange = GatherParticlesToExchangeBetweenThreads(1, CurrentThreadPos.PosX, CurrentThreadPos.PosY, CurrentThreadPos.PosZ, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads, PrintInfo);
        }

        FOR_EACH_THREAD_IN_XYZ
            if (ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles.empty() == false)
            {
                lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MainExchangeParticlesMutexObject };

                for (const auto& ParticleToExchange : ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles)
                    SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert(ParticleToExchange);
            }

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintDataAboutParticlesExchange(PrintInfo, start_time, stop_time, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads);
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsParallelExtract(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCounter = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<ParticlesContainerInternal<Particle>>>> ParticlesToExchange;

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.PosX][CurrentThreadPos.PosY][CurrentThreadPos.PosZ]->MainExchangeParticlesMutexObject };

            ParticlesToExchange = GatherParticlesToExchangeBetweenThreads(1, CurrentThreadPos.PosX, CurrentThreadPos.PosY, CurrentThreadPos.PosZ, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads, PrintInfo);
        }

        FOR_EACH_THREAD_IN_XYZ
            if (ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles.empty() == false)
            {
                lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MainExchangeParticlesMutexObject };

                auto ParticleIterLocalInside = ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles.begin();
                while (ParticleIterLocalInside != ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles.end())
                    SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert(ParticlesToExchange[ThreadXIndex][ThreadYIndex][ThreadZIndex].Particles.extract(ParticleIterLocalInside++));
            }

        const auto stop_time = chrono::high_resolution_clock::now();

        PrintDataAboutParticlesExchange(PrintInfo, start_time, stop_time, ExchangedParticleCounter, StateOfSimulationSpaceDivisionForThreads);
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GatherParticlesFromThreads()
{
    try
    {
        GetParticles().clear();
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ParticleObject : SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads)
                GetParticles().insert(ParticleObject);
    }
    CATCH("gathering cancelled from threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GatherCancelledParticlesIndexesFromThreads()
{
    try
    {
        CancelledParticlesIndexes.clear();

        FOR_EACH_THREAD_IN_XYZ
        {
            for (const auto& CancelledParticleUniqueId : SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->CancelledParticlesIndexes)
                CancelledParticlesIndexes.insert(CancelledParticleUniqueId);

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->CancelledParticlesIndexes.clear();
        }

        LoggersManagerObject.Log(STREAM("NUMBER OF CANCELLED PARTICLES IN CANCELLED REACTIONS = " << CancelledParticlesIndexes.size()));
    }
    CATCH("gathering particles from threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GatherParticlesFromThreadsToParticlesInMainThread()
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

void CellEngineSimulationParallelTasksExecutionManager::SychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        switch (CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInThreads)
        {
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByOneGlobalBarrier : SyncPoint->arrive_and_wait(); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByLocalBarriersOnWallsBetweenThreads : SyncPoint->arrive_and_wait(); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithTwoMutexes : SyncPoint->arrive_and_wait(); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithOneMutex : SyncPoint->arrive_and_wait(); break;
            default : break;
        };
    }
    CATCH("synchronize simulation executed in threads")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsAndSychronizeSimulationExecutedInThreads(barrier<>* SyncPoint, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        switch (CellEngineConfigDataObject.TypeOfSychronizationfSimulationExecutedInThreads)
        {
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByOneGlobalBarrier : ExchangeParticlesBetweenThreadsOneGlobalBarrier(SyncPoint, CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByLocalBarriersOnWallsBetweenThreads : ExchangeParticlesBetweenThreadsLocalBarrier(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithTwoMutexes : ExchangeParticlesBetweenThreadsConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            case CellEngineConfigData::TypesOfSychronizationfSimulationExecutedInThreads::SynchronizeSimulationExecutedInThreadsByConditionalVariableWithOneMutex : ExchangeParticlesBetweenThreadsConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam); break;
            default : break;
        }
    }
    CATCH("exchange particles and synchronize simulation executed in threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                SimulationSpaceSectorBounds SimulationSpaceSectorBoundsObject;

                SimulationSpaceSectorBoundsObject.SetParametersForParallelExecutionSectors({ .PosX = static_cast<SignedInt>(ThreadXIndex), .PosY = static_cast<SignedInt>(ThreadYIndex), .PosZ = static_cast<SignedInt>(ThreadZIndex) }, CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

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
                CurrentThreadLocalSimulationSpaceData->CurrentDiffusionStepStamp = FirstDiffusionStepStampOfCurrentRun + (StepOutside - 1) * NumberOfStepsInside + (Step2 - 1);

                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                for (UnsignedInt ParticleSectorXIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos; ParticleSectorXIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos; ParticleSectorYIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos; ParticleSectorZIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("XStart = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YStart = " << ThreadYIndex * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZStart = " << (ThreadZIndex) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " XEnd = " << (ThreadXIndex + 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YEnd = " << (ThreadZIndex + 1)  * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZEnd = " << (ThreadZIndex + 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, true);
                        }

                SychronizeSimulationExecutedInThreads(SyncPoint, CurrentThreadLocalSimulationSpaceData);

                for (UnsignedInt ParticleSectorXIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos; ParticleSectorXIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos; ParticleSectorYIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos; ParticleSectorZIndex < CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("RANGE1 -> START = (" << ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << "," <<  ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," <<  (ThreadZIndex) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") END = (" <<  (ThreadXIndex + 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << "," << (ThreadYIndex + 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," << (ThreadZIndex + 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("RANGE2 -> START = (" << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") END = (" << SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY << "," << SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ << ") PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("RANGE3 -> START = (" << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos << "," << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos << "," << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos << ") END = (" << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos << "," << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos << "," << CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos << ") PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("RANGE4 -> XStart = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " YStart = " << ThreadYIndex * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace << " ZStart = " << (ThreadZIndex) * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace << " XEnd = " << (ThreadXIndex + 1) * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " YEnd = " << (ThreadYIndex + 1) * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace << " ZEnd = " << (ThreadZIndex + 1) * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInThreads(CurrentThreadLocalSimulationSpaceData, true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, CurrentThreadLocalSimulationSpaceData->SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ);
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

void CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam, const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, latch* AllThreadsReady)
{
    try
    {
        CellEngineThreadPinning::BindCurrentThreadToCPU(CPUsForThreads[(CurrentThreadIndexParam - 1) % CPUsForThreads.size()]);

        AllThreadsReady->count_down();
        StartGatesForThreads[CurrentThreadIndexParam - 1].acquire();

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

void CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInThreads(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        CPUsForThreads = CellEngineThreadPinning::GetCPUsForThreads();

        CellEngineConfigDataObject.MultiThreaded = true;

        SetZeroForAllParallelExecutionVariables();

        bool StateOfSimulationSpaceDivisionForThreads = false;

        std::barrier SyncPoint(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation);

        vector<vector<vector<jthread*>>> Threads(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<jthread*>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<jthread*>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        const UnsignedInt NumberOfAllParticlesBeforeSimulation = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL PARTICLES BEFORE SIMULATION = " << NumberOfAllParticlesBeforeSimulation));
        LoggersManagerObject.Log(STREAM("START SIMULATION IN THREADS"));
        LoggersManagerObject.Log(STREAM("MaxParticleUniqueId = " << MaxParticleUniqueId));

        CellEngineUseful::SwitchOffLogs();

        std::latch AllThreadsReady{ static_cast<ptrdiff_t>(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation) };

        StartGatesForThreads.clear();
        FOR_EACH_THREAD_IN_XYZ
            StartGatesForThreads.emplace_back(0);

        ThreadIdType ThreadIndex = 1;
        FOR_EACH_THREAD_IN_XYZ
        {
            Threads[ThreadXIndex][ThreadYIndex][ThreadZIndex] = new jthread(&CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread, this, &SyncPoint, &StateOfSimulationSpaceDivisionForThreads, NumberOfStepsOutside, NumberOfStepsInside, ThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex, SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex], &AllThreadsReady);
            ThreadIndex++;
        }

        AllThreadsReady.wait();

        const auto start_time = chrono::high_resolution_clock::now();
        const auto StartTime = std::chrono::steady_clock::now();

        ThreadIndex = 0;
        FOR_EACH_THREAD_IN_XYZ
        {
            StartGatesForThreads[ThreadIndex].release();
            ThreadIndex++;
        }

        FOR_EACH_THREAD_IN_XYZ
            Threads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->join();

        const auto StopTime = std::chrono::steady_clock::now();
        const auto stop_time = chrono::high_resolution_clock::now();

        FirstDiffusionStepStampOfCurrentRun += NumberOfStepsOutside * NumberOfStepsInside;

        FOR_EACH_THREAD_IN_XYZ
            delete Threads[ThreadXIndex][ThreadYIndex][ThreadZIndex];

        CellEngineUseful::SwitchOnLogs();

        const string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(StartTime, StopTime, ResultText.c_str(),"Execution in threads")));

        const UnsignedInt NumberOfAllParticlesAfterSimulation = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        LoggersManagerObject.Log(STREAM("END SIMULATION IN THREADS"));
        LoggersManagerObject.Log(STREAM("NUMBER OF ALL PARTICLES AFTER SIMULATION = " << NumberOfAllParticlesAfterSimulation));
        LoggersManagerObject.Log(STREAM("DIFFERENCE OF ALL PARTICLES BEFORE AND AFTER SIMULATION = " << NumberOfAllParticlesAfterSimulation - NumberOfAllParticlesBeforeSimulation));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsByBarrierPhase1(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        LoggersManagerObject.Log(STREAM("NumberOfActiveNeighbors = " << NumberOfActiveNeighbors));

        int Counter = 0;
        for (const auto& ReceivedParticlesToInsert : CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks)
            if (ReceivedParticlesToInsert.empty() == false)
                Counter++;

        LoggersManagerObject.Log(STREAM("NUMBER OF NEIGHBORS THAT SENT PARTCILES TO ACTUAL THREAD = " << Counter << " ACTUAL THREAD = " << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex));
        #endif

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (CurrentThreadLocalSimulationSpaceData->NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
            {
                #if defined(NEIGHBOR_PARALLEL_TASKS_POSITIONS) && defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION)
                const auto [ThreadPosX, ThreadPosY, ThreadPosZ] = CurrentThreadLocalSimulationSpaceData->NeighborParallelTasksPositions[NeighborProcessIndex];
                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("LocalNeighborThreadsIndexes = " << ThreadPosX << " , " << ThreadPosY << " , " << ThreadPosZ << " CurrentThreadIndex = [" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << "]"));)
                #endif

                if (CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].empty() == false)
                {
                    const auto& ReceivedParticlesToInsert = CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex];

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("WANT SENDING CONFIRMATION TO NEIGHBOR = " << ReceivedParticlesToInsert[0].SenderParallelTaskIndex << " " << ReceivedParticlesToInsert[0].ReceiverParallelTaskIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << ")"));)

                    #if defined(NEIGHBOR_PARALLEL_TASKS_POSITIONS) && defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION)
                    LoggersManagerObject.LogUnconditional(STREAM("WANT SENDING CONFIRMATION TO NEIGHBOR = (" << ThreadPosX << " " << ThreadPosY << " " << ThreadPosZ << ") SIZE = " << ReceivedParticlesToInsert.size()));
                    #endif

                    for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                        if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("PARTCLE_KIND_ID TO CHECK = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)

                            CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .PosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .PosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                            if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                            {
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                                auto& ParticleObject = CurrentThreadLocalSimulationSpaceData->GetParticleFromUniqueId(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);
                                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                                CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                                bool FoundNeighbor = false;
                                UnsignedInt LocalNeighborProcessIndex = 0;
                                for (LocalNeighborProcessIndex = 0; LocalNeighborProcessIndex < NumberOfAllNeighbors; LocalNeighborProcessIndex++)
                                    if (SimulationSpaceDataForThreads[ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosX][ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosY][ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosZ]->NeighborParallelTasksIndexes[LocalNeighborProcessIndex] == CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex)
                                    {
                                        FoundNeighbor = true;
                                        break;
                                    }

                                if (FoundNeighbor == true)
                                {
                                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("EXPORT CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ") (" << ReceivedParticleUniqueIdToInsert.Vector.X << "," << ReceivedParticleUniqueIdToInsert.Vector.Y << "," << ReceivedParticleUniqueIdToInsert.Vector.Z << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") PSHIFT = (" << ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X << ", " << ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y << ", " << ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z << ")" << " (" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                                    SimulationSpaceDataForThreads[ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosX][ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosY][ReceivedParticleUniqueIdToInsert.SenderParallelTaskPos.PosZ]->ConfirmationOfParticlesToRemoveToSent[LocalNeighborProcessIndex].emplace_back(ConfirmationOfParticlesToRemoveToSentStruct{ .ParticleUniqueId = ReceivedParticleUniqueIdToInsert.ParticleUniqueId, .SenderSectorPos = { .X = ReceivedParticleUniqueIdToInsert.SenderSectorPos.X, .Y = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Y, .Z = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Z }});

                                    CurrentThreadLocalSimulationSpaceData->NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask++;
                                }
                                #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                                else
                                {
                                    LoggersManagerObject.LogUnconditional(STREAM("ERROR SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));
                                }
                                #endif
                            }
                            else
                            {
                                CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask++;

                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                            }
                        }

                    CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].clear();
                }
            }
    }
    CATCH("exchange particles between threads ver 2")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsByBarrierPhase2(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData)
{
    try
    {
        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (CurrentThreadLocalSimulationSpaceData->NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
            {
                for (const auto& [ParticleUniqueId, SenderSectorPos] : CurrentThreadLocalSimulationSpaceData->ConfirmationOfParticlesToRemoveToSent[NeighborProcessIndex])
                {
                    CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(SenderSectorPos.X), .PosY = static_cast<SignedInt>(SenderSectorPos.Y), .PosZ = static_cast<SignedInt>(SenderSectorPos.Z) });

                    CurrentThreadLocalSimulationSpaceData->RemoveParticle(ParticleUniqueId, true, false);

                    CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask++;

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO REMOVE THIS PARTCILE = " << ParticleUniqueId << " (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")"));)
                }

                CurrentThreadLocalSimulationSpaceData->ConfirmationOfParticlesToRemoveToSent[NeighborProcessIndex].clear();
            }
    }
    CATCH("exchange particles threads processes group 3")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsOneGlobalBarrier(barrier<>* SyncPoint, const std::shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
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

void CellEngineSimulationParallelTasksExecutionManager::SynchronizeWithNeighborByLocalBarrier(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData) const
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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsLocalBarrier(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData) const
{
    try
    {
        const UnsignedInt SlotIndex = CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber;

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

                Wall.Proposals[SlotIndex] = std::move(CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex]);
                Wall.ProposalsReady[SlotIndex] = true;
            }

            CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].clear();

            Wall.ProposalConditionalVariable.notify_all();
        }

        vector<ConfirmationOfParticlesToRemoveToSentStruct> ReceivedConfirmationOfParticlesToRemove;

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

            Wall.Verdicts[SlotIndex].clear();
            Wall.VerdictsReady[SlotIndex] = false;
        }

        for (const auto& [ParticleUniqueId, SenderSectorPos] : ReceivedConfirmationOfParticlesToRemove)
        {
            CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(SenderSectorPos.X), .PosY = static_cast<SignedInt>(SenderSectorPos.Y), .PosZ = static_cast<SignedInt>(SenderSectorPos.Z) });

            CurrentThreadLocalSimulationSpaceData->RemoveParticle(ParticleUniqueId, true, false);

            CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask++;

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO REMOVE THIS PARTCILE = " << ParticleUniqueId << " (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")"));)
        }
    }
    CATCH("exchange particles threads processes")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        const UnsignedInt SlotIndex = 1 - CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            std::vector<ParticleSenderStructForParallelMultiTaskingExecution> ReceivedParticlesToInsert;
            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1);
                Wall.ProposalConditionalVariable.wait(LockGuardScopedLock1, [&] { return Wall.ProposalsReady[SlotIndex]; });

                ReceivedParticlesToInsert = std::move(Wall.Proposals[SlotIndex]);
                Wall.Proposals[SlotIndex].clear();
                Wall.ProposalsReady[SlotIndex] = false;
            }

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("SENDING CONFIRMATION TO Neighbor = " << ReceivedParticlesToInsert[0].SenderParallelTaskIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

            std::vector<ConfirmationOfParticlesToRemoveToSentStruct> ConfirmationOfParticlesToRemoveToSent;

            for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                {
                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("PARTCLE_KIND_ID TO CHECK = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)

                    CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .PosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .PosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                    if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                    {
                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                        auto& ParticleObject = CurrentThreadLocalSimulationSpaceData->GetParticleFromUniqueId(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);
                        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                        CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                        ConfirmationOfParticlesToRemoveToSent.emplace_back(ConfirmationOfParticlesToRemoveToSentStruct{ .ParticleUniqueId = ReceivedParticleUniqueIdToInsert.ParticleUniqueId, .SenderSectorPos = { .X = ReceivedParticleUniqueIdToInsert.SenderSectorPos.X, .Y = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Y, .Z = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Z }});

                       DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("EXPORT CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ") (" << ReceivedParticleUniqueIdToInsert.Vector.X << "," << ReceivedParticleUniqueIdToInsert.Vector.Y << "," << ReceivedParticleUniqueIdToInsert.Vector.Z << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") PSHIFT = (" << ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X << ", " << ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y << ", " << ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z << ")" << " (" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                        CurrentThreadLocalSimulationSpaceData->NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask++;
                    }
                    else
                    {
                        CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask++;

                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)
                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                    }
                }
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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsConditionalVariableTwoMutexes(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        if (CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber == 0)
        {
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
        }
        else
        {
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableTwoMutexes(CurrentThreadLocalSimulationSpaceData);
        }
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData) const
{
    try
    {
        const UnsignedInt SlotIndex = CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];

            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            {
                lock_guard<mutex> LockGuardScopedLock1(Wall.MutexNumber1);

                Wall.Proposals[SlotIndex] = std::move(CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex]);
                Wall.ProposalsReady[SlotIndex] = true;
            }

            CurrentThreadLocalSimulationSpaceData->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].clear();

            Wall.WallConditionalVariable.notify_all();
        }

        vector<ConfirmationOfParticlesToRemoveToSentStruct> ReceivedConfirmationOfParticlesToRemove;

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

        for (const auto& [ParticleUniqueId, SenderSectorPos] : ReceivedConfirmationOfParticlesToRemove)
        {
            CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(SenderSectorPos.X), .PosY = static_cast<SignedInt>(SenderSectorPos.Y), .PosZ = static_cast<SignedInt>(SenderSectorPos.Z) });

            CurrentThreadLocalSimulationSpaceData->RemoveParticle(ParticleUniqueId, true, false);

            CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask++;

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO REMOVE THIS PARTCILE = " << ParticleUniqueId << " (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")"));)
        }
    }
    CATCH("exchange particles between threads group 1")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        const UnsignedInt SlotIndex = 1 - CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
        {
            const SignedInt WallIndex = CurrentThreadLocalSimulationSpaceData->TwoThreadsWallSychronizationBarriersIndexes[NeighborProcessIndex];
            if (WallIndex < 0)
                continue;

            WallChannel& Wall = *TwoThreadsWallChannels[WallIndex];

            std::vector<ParticleSenderStructForParallelMultiTaskingExecution> ReceivedParticlesToInsert;
            {
                unique_lock<mutex> LockGuardScopedLock1(Wall.MutexNumber1);
                Wall.WallConditionalVariable.wait(LockGuardScopedLock1, [&] { return Wall.ProposalsReady[SlotIndex]; });

                ReceivedParticlesToInsert = std::move(Wall.Proposals[SlotIndex]);
                Wall.Proposals[SlotIndex].clear();
                Wall.ProposalsReady[SlotIndex] = false;
            }

            std::vector<ConfirmationOfParticlesToRemoveToSentStruct> ConfirmationOfParticlesToRemove;

            for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                {
                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("PARTCLE_KIND_ID TO CHECK = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)

                    CurrentThreadLocalSimulationSpaceData->SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .PosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .PosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                    if (CurrentThreadLocalSimulationSpaceData->CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMultiThreadedDiffusion(ReceivedParticleUniqueIdToInsert, ReceivedParticleUniqueIdToInsert.ParticleObject) == true)
                    {
                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                        auto& ParticleObject = CurrentThreadLocalSimulationSpaceData->GetParticleFromUniqueId(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);
                        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                        CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                        ConfirmationOfParticlesToRemove.emplace_back(ConfirmationOfParticlesToRemoveToSentStruct{ .ParticleUniqueId = ReceivedParticleUniqueIdToInsert.ParticleUniqueId, .SenderSectorPos = { .X = ReceivedParticleUniqueIdToInsert.SenderSectorPos.X, .Y = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Y, .Z = ReceivedParticleUniqueIdToInsert.SenderSectorPos.Z }});

                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("EXPORT CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ") (" << ReceivedParticleUniqueIdToInsert.Vector.X << "," << ReceivedParticleUniqueIdToInsert.Vector.Y << "," << ReceivedParticleUniqueIdToInsert.Vector.Z << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") PSHIFT = (" << ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X << ", " << ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y << ", " << ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z << ")" << " (" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadIndexParam << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)

                        CurrentThreadLocalSimulationSpaceData->NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask++;
                    }
                    else
                    {
                        CurrentThreadLocalSimulationSpaceData->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask++;

                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosX << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosY << " " << CurrentThreadLocalSimulationSpaceData->CurrentThreadPos.PosZ << ")["  << CurrentThreadLocalSimulationSpaceData->CurrentThreadIndex << " " << ThreadXIndexParam << " " << ThreadYIndexParam << " " << ThreadZIndexParam << "] NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)
                        DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                    }
                }

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

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenThreadsConditionalVariableOneMutex(const shared_ptr<CellEngineSimulationSpace>& CurrentThreadLocalSimulationSpaceData, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        if (CurrentThreadLocalSimulationSpaceData->ParallelTaskGroupNumber == 0)
        {
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
        }
        else
        {
            ExchangeParticlesBetweenThreadsGroup2ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData, CurrentThreadIndexParam, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
            ExchangeParticlesBetweenThreadsGroup1ConditionalVariableOneMutex(CurrentThreadLocalSimulationSpaceData);
        }
    }
    CATCH("exchange particles between threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                CurrentDiffusionStepStamp = FirstDiffusionStepStampOfCurrentRun + (StepOutside - 1) * NumberOfStepsInside + (Step2 - 1);

                DEBUG_LOG_DIFFUSION_OR_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));)

                for (UnsignedInt ParticleSectorXIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos; ParticleSectorXIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos; ParticleSectorYIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos; ParticleSectorZIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("XStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos << " YStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos << " ZStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos << " XEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos << " YEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos << " ZEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType, true);
                        }

                MPI_Barrier(MPI_COMM_WORLD);

                for (UnsignedInt ParticleSectorXIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos; ParticleSectorXIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos; ParticleSectorYIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos; ParticleSectorZIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos; ParticleSectorZIndex++)
                        {
                            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("XStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos << " YStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos << " ZStart = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos << " XEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos << " YEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos << " ZEnd = " <<  SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));)

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpaceForExecutionInMPIProcesses(true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ, SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos * CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterX, SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos * CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterY, SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos * CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace - CellEngineConfigDataObject.ShiftCenterZ);
                        }

                MPI_Barrier(MPI_COMM_WORLD);
            }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam)
{
    try
    {
        for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
        {
            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("NUMBER OF STEPS OUTSIDE = " << NumberOfStepsInside << " " << StepOutside));)

            GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);

            const auto start_time_exchange = chrono::high_resolution_clock::now();

            ExchangeParticlesBetweenMPIProcesses();

            const auto stop_time_exchange = chrono::high_resolution_clock::now();

            ParallelTaskGroupNumber = 1 - ParallelTaskGroupNumber;

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time_exchange, stop_time_exchange, "Only exchanging of particles has taken time","Execution in mpi")));)
        }
    }
    CATCH("generating n steps of simulation for whole cell space in one thread")
}

void CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(UnsignedInt NumberOfStepsOutside, UnsignedInt NumberOfStepsInside)
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

        //this_thread::sleep_for(chrono::milliseconds(300));

        MPI_Barrier(MPI_COMM_WORLD);

        const double StartTime = MPI_Wtime();

        const auto start_time = chrono::high_resolution_clock::now();

        GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsOutside, NumberOfStepsInside, MPIProcessDataObject.CurrentMPIProcessIndex, 0, 0, 0);

        FirstDiffusionStepStampOfCurrentRun += NumberOfStepsOutside * NumberOfStepsInside;

        MPI_Barrier(MPI_COMM_WORLD);

        const double MyTime = MPI_Wtime() - StartTime;

        const auto stop_time = chrono::high_resolution_clock::now();

        double MaxTime = 0.0;
        double MinTime = 0.0;
        MPI_Reduce(&MyTime, &MaxTime, 1, MPI_DOUBLE, MPI_MAX, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);
        MPI_Reduce(&MyTime, &MinTime, 1, MPI_DOUBLE, MPI_MIN, CellEngineConfigDataObject.MainMPIProcessNumber, MPI_COMM_WORLD);

        CellEngineUseful::SwitchOnLogs();

        if (MPIProcessDataObject.CurrentMPIProcessIndex == CellEngineConfigDataObject.MainMPIProcessNumber)
            LoggersManagerObject.Log(STREAM("Execution in MPI processes for steps outside -> parallel time = " << MaxTime << " s (all processes measured between " << MinTime << " s and " << MaxTime << " s)"));

        const string ResultText = "Execution in MPI processes for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in mpi")));

        LoggersManagerObject.Log(STREAM("END MPI SIMULATION"));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in mpi process")
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup1()
{
    CellEngineUseful::SwitchOnLogs();

    try
    {
        MainMapType<UniqueIdUInt64Type, SectorPosType> AllParticlesSentToAnotherProcess;

        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
            if (NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
            {
                if (VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].empty() == true)
                {
                    VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].emplace_back(ParticleSenderStructForParallelMultiTaskingExecution{ 0, 0, static_cast<int>(MPIProcessDataObject.CurrentMPIProcessIndex), static_cast<int>(NeighborParallelTasksIndexes[NeighborProcessIndex]), 0, 0, 0, 0, 0, 0 });

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Index to EMPTY send = " << NeighborParallelTasksIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex][0].ReceiverParallelTaskIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)
                }
                else
                {
                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Index to send = " << NeighborParallelTasksIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex][0].ReceiverParallelTaskIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)
                }

                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Length Message SEND = " << VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].size() << " " << VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].size() * sizeof(ParticleSenderStructForParallelMultiTaskingExecution) << " " << NeighborParallelTasksIndexes[NeighborProcessIndex] << " " << VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex][0].ReceiverParallelTaskIndex << " " << NeighborProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                char BufferToSend[MaxMPIMessageSize];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].size();

                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI PACKED SIZE = " << NumberOfPackedStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);

                for (const auto& ParticleToSendElement : VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex])
                {
                    MPI_Pack(&ParticleToSendElement.ParticleUniqueId, 1, MPI_UINT64_T, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ParticleKindId, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SenderParallelTaskIndex, 1, MPI_INT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverParallelTaskIndex, 1, MPI_INT, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
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

                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, static_cast<int>(NeighborParallelTasksIndexes[NeighborProcessIndex]), 0, MPI_COMM_WORLD);

                VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].clear();
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

                    NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask++;
                }

            NumberOfReceivedMessages++;
        }
    }
    CATCH("exchange particles between mpi processes")

    CellEngineUseful::SwitchOffLogs();
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup2()
{
    CellEngineUseful::SwitchOnLogs();

    try
    {
        int NumberOfReceivedMessages = 0;

        vector<ParticleSenderStructForParallelMultiTaskingExecution> ReceivedParticlesToInsertFromAllNeigbhours[NumberOfAllNeighbors];

        while (NumberOfReceivedMessages < NumberOfActiveNeighbors)
        {
            MPI_Status MPIMessageStatus;

            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int NumberOfBytesReceived;
            MPI_Get_count(&MPIMessageStatus, MPI_CHAR, &NumberOfBytesReceived);

            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("MPI Process Length Message RECEIVE = " << NumberOfBytesReceived << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << NumberOfBytesReceived / sizeof(ParticleSenderStructForParallelMultiTaskingExecution)));)

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
                ParticleSenderStructForParallelMultiTaskingExecution MPIParticleSenderStructElementLocalObject;
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleUniqueId, 1, MPI_UINT64_T, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleKindId, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SenderParallelTaskIndex, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(ReceivedParticlesToInsert1, MaxMPIMessageSize, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverParallelTaskIndex, 1, MPI_INT, MPI_COMM_WORLD);
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
                    if (NeighborParallelTasksIndexes[NeighborProcessIndex] == MPIParticleSenderStructElementLocalObject.SenderParallelTaskIndex)
                    {
                        ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex].emplace_back(MPIParticleSenderStructElementLocalObject);
                        Found = true;
                        break;
                    }

                #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                if (Found == false)
                    LoggersManagerObject.Log(STREAM("GET PROCESS INDEX FROM = " << NeighborProcessIndex << " " << NeighborParallelTasksIndexes[NeighborProcessIndex] << " " << MPIParticleSenderStructElementLocalObject.SenderParallelTaskIndex << " TO = " << MPIParticleSenderStructElementLocalObject.ReceiverParallelTaskIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex << endl));
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
            if (NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
            {
                vector<UniqueIdUInt64Type> ConfirmationOfParticlesToRemoveToSent;

                if (ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex].empty() == false)
                {
                    const auto ReceivedParticlesToInsert = ReceivedParticlesToInsertFromAllNeigbhours[NeighborProcessIndex];

                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("SENDING CONFIRMATION TO Neighbor = " << ReceivedParticlesToInsert[0].SenderParallelTaskIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));)

                    for (const auto& ReceivedParticleUniqueIdToInsert : ReceivedParticlesToInsert)
                        if (ReceivedParticleUniqueIdToInsert.ParticleUniqueId != 0)
                        {
                            SetCurrentSectorPos(SectorPosType{ .PosX = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.X), .PosY = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Y), .PosZ = static_cast<SignedInt>(ReceivedParticleUniqueIdToInsert.ReceiverSectorPos.Z) });

                            UniqueIdUInt64Type UniqueIdOfNewCreatedParticle;
                            if (CheckPossibilityOfInsertingParticleToCurrentSectorAndInsertIfPossibleInMPIMultiProcessDiffusion(ReceivedParticleUniqueIdToInsert, UniqueIdOfNewCreatedParticle) == true)
                            {
                                auto& ParticleObject = GetParticleFromUniqueId(UniqueIdOfNewCreatedParticle);
                                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + ReceivedParticleUniqueIdToInsert.Vector.X, ParticleObject.Center.Y + ReceivedParticleUniqueIdToInsert.Vector.Y, ParticleObject.Center.Z + ReceivedParticleUniqueIdToInsert.Vector.Z);
                                CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, ReceivedParticleUniqueIdToInsert.Vector.X, ReceivedParticleUniqueIdToInsert.Vector.Y, ReceivedParticleUniqueIdToInsert.Vector.Z);

                                ConfirmationOfParticlesToRemoveToSent.emplace_back(ReceivedParticleUniqueIdToInsert.ParticleUniqueId);

                                NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask++;
                            }
                            else
                            {
                                NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask++;

                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("NOT SENDING CONFIRMATION TO NEIGHBOR = (PARTICLE_INDEX = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId << ")(" << ReceivedParticleUniqueIdToInsert.SenderParallelTaskIndex << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskIndex << ") (" << MPIProcessIndex << " ) NEIGHBOR = (" << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosX << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosY << " " << ReceivedParticleUniqueIdToInsert.ReceiverParallelTaskPos.PosZ << ")"));)
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("PARTCLE_INDEX TO NOT SEND CONFIRMATION = " << ReceivedParticleUniqueIdToInsert.ParticleUniqueId));)
                            }
                        }
                }

                if (ConfirmationOfParticlesToRemoveToSent.empty() == true)
                    ConfirmationOfParticlesToRemoveToSent.emplace_back(0);

                char BufferToSend[MaxMPIMessageSize];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = ConfirmationOfParticlesToRemoveToSent.size();

                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                for (const auto& ParticleToSendElement : ConfirmationOfParticlesToRemoveToSent)
                    MPI_Pack(&ParticleToSendElement, 1, MPI_UINT64_T, BufferToSend, MaxMPIMessageSize, &PositionInBuffer, MPI_COMM_WORLD);
                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, NeighborParallelTasksIndexes[NeighborProcessIndex], 0, MPI_COMM_WORLD);
            }
    }
    CATCH("exchange particles between mpi processes ver 2")

    CellEngineUseful::SwitchOffLogs();
}

void CellEngineSimulationParallelTasksExecutionManager::ExchangeParticlesBetweenMPIProcesses()
{
    try
    {
        if (ParallelTaskGroupNumber == 0)
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

void CellEngineSimulationParallelTasksExecutionManager::SetZeroForAllParallelExecutionVariables()
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
        NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask = 0;
        NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask = 0;
        NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected = 0;

        NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask = 0;
        NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask = 0;
        NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask  = 0;
        NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfDNAParticles = 0;

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ErrorCounter = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfExecutedReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledAReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledBReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->AddedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->RemovedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->RestoredParticlesInCancelledReactions = 0;

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusion = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected = 0;

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask = 0;
            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask = 0;
        }
    }
    CATCH("setting zero for all parallel execution variables")
}

void CellEngineSimulationParallelTasksExecutionManager::PrintAllExecutionVariables(const UnsignedInt PhaseId)
{
    try
    {
        LoggersManagerObject.Log(STREAM("DIFFUSION STATISTICS:"));

        NumberOfParticlesMovedInDiffusion = NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask + NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask + NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask;

        LoggersManagerObject.Log(STREAM("NumberOfParticlesMovedInDiffusion = " + to_string(NumberOfParticlesMovedInDiffusion) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask = " + to_string(NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask = " + to_string(NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected = " + to_string(NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));

        LoggersManagerObject.Log(STREAM("NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask = " + to_string(NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask = " + to_string(NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask = " + to_string(NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));
        LoggersManagerObject.Log(STREAM("NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask = " + to_string(NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask) << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << PhaseId));

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

void CellEngineSimulationParallelTasksExecutionManager::CountAllParticlesInMPI()
{
    try
    {
        NumberOfAllParticlesInMPISimulation = 0;

        for (UnsignedInt ParticleSectorXIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartXPos; ParticleSectorXIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndXPos; ParticleSectorXIndex++)
            for (UnsignedInt ParticleSectorYIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartYPos; ParticleSectorYIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndYPos; ParticleSectorYIndex++)
                for (UnsignedInt ParticleSectorZIndex = SimulationSpaceSectorsRangesForCurrentParallelTask.StartZPos; ParticleSectorZIndex < SimulationSpaceSectorsRangesForCurrentParallelTask.EndZPos; ParticleSectorZIndex++)
                    NumberOfAllParticlesInMPISimulation += CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(SectorPosType(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex));

        MPI_Allreduce(MPI_IN_PLACE, &NumberOfAllParticlesInMPISimulation, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    }
    CATCH("counting all particles in MPI")
}

void CellEngineSimulationParallelTasksExecutionManager::GatherAllParallelExecutionVariables()
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

            NumberOfParticlesMovedInDiffusion = NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask + NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask + NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask;

            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesMovedInDiffusion, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        }
        else
        if (CellEngineConfigDataObject.MultiThreaded == true)
            FOR_EACH_THREAD_IN_XYZ
            {
                NumberOfDNAParticles += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfDNAParticles;

                ErrorCounter += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ErrorCounter;
                NumberOfExecutedReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfExecutedReactions;
                NumberOfCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledReactions;
                NumberOfCancelledAReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledAReactions;
                NumberOfCancelledBReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfCancelledBReactions;
                AddedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->AddedParticlesInReactions;
                RemovedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->RemovedParticlesInReactions;
                RestoredParticlesInCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->RestoredParticlesInCancelledReactions;

                NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionWithinInOneSectorWithinOneParallelTask;
                NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTask;
                NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionSentFromOneSectorToAnotherSectorWithinOneParallelTaskRejected;

                NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesMovedInDiffusionSentToAnotherParallelTask;
                NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask;
                NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask  += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionRejectedByAnotherParallelTask;
                NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask += SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->NumberOfParticlesToBeMovedInDiffusionConfirmedAfterSentInAnotherParallelTask;
            }

        PrintAllExecutionVariables(2);
    }
    CATCH("gathering all parallel execution variables")
}

void CellEngineSimulationParallelTasksExecutionManager::GenerateNStepsOfSimulationWithSendingParticlesToThreadsAndGatheringParticlesToMainThreadForWholeCellSpace(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, bool PrintTime)
{
    try
    {
        FirstSendParticlesForThreads(false, true);
        GenerateNStepsOfSimulationForWholeCellSpaceInThreads(NumberOfStepsOutside, NumberOfStepsInside);
        GatherParticlesFromThreadsToParticlesInMainThread();
    }
    CATCH("generate n steps of simulation with sending particles to threads and gathering particles to main threads for whole cell space")
}

void CellEngineSimulationParallelTasksExecutionManager::SaveFormerParticlesAsVectorElements()
{
    try
    {
        FormerParticlesIndexes = GetParticles();
    }
    CATCH("saving former particles as vector elements")
}

void CellEngineSimulationParallelTasksExecutionManager::CheckParticlesCenters(const bool PrintAllParticles)
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
