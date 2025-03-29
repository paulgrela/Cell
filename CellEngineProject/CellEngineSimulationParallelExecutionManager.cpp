
#include <mpi.h>
#include <chrono>

#include "DateTimeUtils.h"

#include "CellEngineDataFile.h"
#include "CellEngineImGuiMenu.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineSimulationParallelExecutionManager.h"

#include <mpi.h>

#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

using namespace std;

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

                    ThreadLocalParticlesInProximityZPos = std::make_shared<SimulationSpaceType>(Particles, false, ThreadIndexPos, ThreadPosType{ ThreadXPos, ThreadYPos, ThreadZPos });
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

UnsignedInt GetProcessGroupNumberVer2(const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    if (ThreadZIndex % 2 == 0)
        return ThreadXIndex % 2 != ThreadYIndex % 2;
    else
        return ThreadXIndex % 2 == ThreadYIndex % 2;
}

SignedInt CellEngineSimulationParallelExecutionManager::GetProcessPrevNeighbour(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (ThreadXIndex >= 0 && ThreadYIndex >= 0 && ThreadZIndex >= 0)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex - 1;
    else
        return -1;
}

SignedInt CellEngineSimulationParallelExecutionManager::GetProcessNextNeighbour(const SignedInt ThreadXIndex, const SignedInt ThreadYIndex, const SignedInt ThreadZIndex) const
{
    if (static_cast<SignedInt>(ThreadXIndex) < CellEngineConfigDataObject.NumberOfXThreadsInSimulation && static_cast<SignedInt>(ThreadYIndex) < CellEngineConfigDataObject.NumberOfYThreadsInSimulation && static_cast<SignedInt>(ThreadZIndex) < CellEngineConfigDataObject.NumberOfZThreadsInSimulation)
        return SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->MPIProcessIndex - 1;
    else
        return -1;
}

void CellEngineSimulationParallelExecutionManager::CreateDataEveryMPIProcessForParallelExecution()
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

                        NeighbourProcessesIndexes[0] = GetProcessPrevNeighbour(static_cast<SignedInt>(MPIProcessXIndex) - 2, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighbourProcessesIndexes[1] = GetProcessPrevNeighbour(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 2, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighbourProcessesIndexes[2] = GetProcessPrevNeighbour(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 2);

                        NeighbourProcessesIndexes[3] = GetProcessNextNeighbour(static_cast<SignedInt>(MPIProcessXIndex), static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighbourProcessesIndexes[4] = GetProcessNextNeighbour(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex), static_cast<SignedInt>(MPIProcessZIndex) - 1);
                        NeighbourProcessesIndexes[5] = GetProcessNextNeighbour(static_cast<SignedInt>(MPIProcessXIndex) - 1, static_cast<SignedInt>(MPIProcessYIndex) - 1, static_cast<SignedInt>(MPIProcessZIndex));

                        LoggersManagerObject.Log(STREAM("MPIProcessIndex Neighbours = " <<  MPIProcessDataObject.CurrentMPIProcessIndex << " " << NeighbourProcessesIndexes[0] << " " << NeighbourProcessesIndexes[1] << " " << NeighbourProcessesIndexes[2] << " " << NeighbourProcessesIndexes[3] << " " << NeighbourProcessesIndexes[4] << " " << NeighbourProcessesIndexes[5]));
                    }

                    MPIProcessIndex++;
                }

        NumberOfActiveNeighbours = count_if(NeighbourProcessesIndexes, NeighbourProcessesIndexes + NumberOfAllNeighbours, [](const auto& Element){ return Element != -1; });

        LoggersManagerObject.Log(STREAM("NumberOfActiveNeighbours = " << NumberOfActiveNeighbours));
    }
    CATCH("creating data for every mpi process for parallel execution")
}

#define FOR_EACH_THREAD_IN_XYZ \
        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++) \
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++) \
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)

#ifdef SHORTER_CODE
void CellEngineSimulationParallelExecutionManager::JoinReactionsStatisticsFromThreads(vector<map<UnsignedInt, ReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true)
        {
            if (CellEngineConfigDataObject.OpenGLGraphicsSwitchedOff == false && MPIProcessDataObject.CurrentMPIProcessIndex == 0)
            {
                int ValueToSend = 3;
                MPI_Bcast(&ValueToSend, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }

            for (const auto& ReactionStatisticsData : SavedReactionsMap.back())
            {
                LoggersManagerObject.LogUnconditional(STREAM("ReactionStatisticsData S = " << ReactionStatisticsData.first << " " << ReactionStatisticsData.second.Counter << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

                CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].emplace_back(ReactionStatisticsData.first, ReactionStatisticsData.second.Counter);
            }

            const int SavedReactionsMapForMPILength = CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].size();
            int SavedReactionsMapForMPILengths[MPIProcessDataObject.NumberOfMPIProcesses];
            MPI_Gather(&SavedReactionsMapForMPILength, 1, MPI_INT, SavedReactionsMapForMPILengths, 1, MPI_INT, 0, MPI_COMM_WORLD);

            LoggersManagerObject.LogUnconditional(STREAM("SavedReactionsMapForMPILength = " << SavedReactionsMapForMPILength << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

            int MaximumOfAllSavedReactionsMapForMPILengths;
            if (MPIProcessDataObject.CurrentMPIProcessIndex == 0)
                MaximumOfAllSavedReactionsMapForMPILengths = *max_element(SavedReactionsMapForMPILengths, SavedReactionsMapForMPILengths + MPIProcessDataObject.NumberOfMPIProcesses);

            MPI_Bcast(&MaximumOfAllSavedReactionsMapForMPILengths, 1, MPI_INT, 0, MPI_COMM_WORLD);

            LoggersManagerObject.LogUnconditional(STREAM("MaximumOfAllSavedReactionsMapForMPILengths = " << MaximumOfAllSavedReactionsMapForMPILengths << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

            unique_ptr<ReactionStatistics[]> ReactionStatisticsVectorGatheringFromAllMPIProcessPointer(new ReactionStatistics[MaximumOfAllSavedReactionsMapForMPILengths * MPIProcessDataObject.NumberOfMPIProcesses + 1]);

            const int NumberOfBytesToGatherFromEveryProcess = static_cast<int>(MaximumOfAllSavedReactionsMapForMPILengths * sizeof(ReactionStatistics));
            MPI_Gather(CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->SavedReactionsMapForMPI[SimulationStepNumber - 1].data(), NumberOfBytesToGatherFromEveryProcess, MPI_BYTE, ReactionStatisticsVectorGatheringFromAllMPIProcessPointer.get(), NumberOfBytesToGatherFromEveryProcess, MPI_BYTE, 0, MPI_COMM_WORLD);

            if (MPIProcessDataObject.CurrentMPIProcessIndex == 0)
            {
                SavedReactionsMap[SimulationStepNumber - 1].clear();
                for (UniqueIdInt LocalMPIProcessIndex = 0; LocalMPIProcessIndex < MPIProcessDataObject.NumberOfMPIProcesses; LocalMPIProcessIndex++)
                    for (UnsignedInt ReactionStatisticsDataIndex = MaximumOfAllSavedReactionsMapForMPILengths * LocalMPIProcessIndex; ReactionStatisticsDataIndex < MaximumOfAllSavedReactionsMapForMPILengths * LocalMPIProcessIndex + SavedReactionsMapForMPILengths[LocalMPIProcessIndex]; ReactionStatisticsDataIndex++)
                    {
                        const auto& ReactionStatisticsData = ReactionStatisticsVectorGatheringFromAllMPIProcessPointer.get()[ReactionStatisticsDataIndex];
                        LoggersManagerObject.LogUnconditional(STREAM("ReactionStatisticsData 1 = " << ReactionStatisticsDataIndex << " " << ReactionStatisticsData.ReactionId << " " << ReactionStatisticsData.Counter << " LocalMPIProcessIndex = " << LocalMPIProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                        if (ReactionStatisticsData.ReactionId != 0 && ReactionStatisticsData.Counter != 0)
                        {
                            LoggersManagerObject.LogUnconditional(STREAM("ReactionStatisticsData 2 = " << ReactionStatisticsData.ReactionId << " " << ReactionStatisticsDataIndex << " " << ReactionStatisticsData.Counter << " LocalMPIProcessIndex = " << LocalMPIProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                            SavedReactionsMap[SimulationStepNumber - 1][ReactionStatisticsData.ReactionId] = { ReactionStatisticsData.ReactionId, SavedReactionsMap[SimulationStepNumber - 1][ReactionStatisticsData.ReactionId].Counter += ReactionStatisticsData.Counter };
                        }
                    }
            }
        }
        else
            FOR_EACH_THREAD_IN_XYZ
                for (const auto& ReactionStatisticsData : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap.back())
                    SavedReactionsMap[SimulationStepNumber - 1][ReactionStatisticsData.first].Counter += ReactionStatisticsData.second.Counter;
    }
    CATCH("joining statistics from threads")
}
#else
void CellEngineSimulationParallelExecutionManager::JoinStatisticsFromThreads(vector<map<UnsignedInt, ReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ReactionData : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap[SimulationStepNumber - 1])
                if (SavedReactionsMap[SimulationStepNumber - 1].contains(ReactionData.second.ReactionId))
                    SavedReactionsMap[SimulationStepNumber - 1].find(ReactionData.second.ReactionId)->second.Counter += ReactionData.second.Counter;
                else
                    SavedReactionsMap[SimulationStepNumber - 1].insert(ReactionData);
    }
    CATCH("joining statistics from threads extended")
}
#endif

void LogCenterOfParticleWithThreadIndex(const Particle& ParticleObject, const ThreadIdType ThreadXIndex, const ThreadIdType ThreadYIndex, const ThreadIdType ThreadZIndex)
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
                UnsignedInt ThreadXIndex = floor(ParticleObject.second.Center.X / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
                UnsignedInt ThreadYIndex = floor(ParticleObject.second.Center.Y / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
                UnsignedInt ThreadZIndex = floor(ParticleObject.second.Center.Z / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if (PrintCenterOfParticleWithThreadIndex == true)
                    LogCenterOfParticleWithThreadIndex(ParticleObject.second, ThreadXIndex, ThreadYIndex, ThreadZIndex);

                SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert(ParticleObject);

                GoodParticlesCounter++;
            }
            else
                BadParticlesCounter++;

        FOR_EACH_THREAD_IN_XYZ
        {
            InitiateFreeParticleIndexes(SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads, false);

            if (PrintCenterOfParticleWithThreadIndex == true)
                LoggersManagerObject.Log(STREAM("THREAD[" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << "] SIZE = " << SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.size()));
        }

        const auto stop_time = chrono::high_resolution_clock::now();

        if (PrintTime == true)
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "First sending particles to threads has taken time = ","Execution in threads")));

        LoggersManagerObject.Log(STREAM("All particles counter = " << GoodParticlesCounter << " Bad particles counter = " << BadParticlesCounter));
    }
    CATCH("first sending particles for threads")
}

void PrintThreadIndexes(const Particle& ParticleObject, const ThreadIdType CurrentThreadIndex, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const UnsignedInt ThreadXIndexNew, const UnsignedInt ThreadYIndexNew, const UnsignedInt ThreadZIndexNew)
{
    LoggersManagerObject.Log(STREAM("Particle Center: " << ParticleObject.Index << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId).ParticleKindSpecialDataSector[0].ParticleType)));
    LoggersManagerObject.Log(STREAM("THREAD POS NEW = " << ThreadXIndexNew << ", " << ThreadYIndexNew << ", " << ThreadZIndexNew));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex));
    LoggersManagerObject.Log(STREAM(endl));
}

void PrintDataAboutParticlesExchange(const bool PrintInfo, const chrono::time_point<chrono::system_clock, chrono::system_clock::duration>& start_time, const chrono::time_point<chrono::system_clock, chrono::system_clock::duration>& stop_time, const UnsignedInt ExchangedParticleCounter, const bool StateOfSimulationSpaceDivisionForThreads)
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

            UnsignedInt ThreadXIndexNew = floor((ParticleIter->second.Center.X - SizeOfXInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
            UnsignedInt ThreadYIndexNew = floor((ParticleIter->second.Center.Y - SizeOfYInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
            UnsignedInt ThreadZIndexNew = floor((ParticleIter->second.Center.Z - SizeOfZInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

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
            for (const auto& CancelledParticleIndex : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->CancelledParticlesIndexes)
                CancelledParticlesIndexes.insert(CancelledParticleIndex);

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

        auto RememberStateOfUseMutexBetweenMainScreenThreadAndMenuThreads = CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads;
        CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads = true;
        sleep(1);
        {
            lock_guard<recursive_mutex> LockGuard{ CellEngineOpenGLVisualiserOfVoxelSimulationSpace::RenderMenuAndVoxelSimulationSpaceMutexObject };

            GatherParticlesFromThreads();

            InitiateFreeParticleIndexes(GetParticles(), false);

            GatherCancelledParticlesIndexesFromThreads();
        }
        CellEngineConfigDataObject.UseMutexBetweenMainScreenThreadAndMenuThreads = RememberStateOfUseMutexBetweenMainScreenThreadAndMenuThreads;

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Gathering particles from threads to main thread has taken time = ","Execution in threads")));
    }
    CATCH("gathering particles from threads to particles in main thread")
}

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads, barrier<>* SyncPoint)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));

                SimulationSpaceSectorBounds SimulationSpaceSectorBoundsObject = SimulationSpaceSectorBoundsObject.SetParametersForParallelExecutionSectors({ ThreadXIndex, ThreadYIndex, ThreadZIndex }, CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if (StateOfSimulationSpaceDivisionForThreads == true)
                    SimulationSpaceSectorBoundsObject.AddToStartParameters(CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2, CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2, CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2);

                for (UnsignedInt PosX = SimulationSpaceSectorBoundsObject.StartXPos; PosX < SimulationSpaceSectorBoundsObject.EndXPos; PosX += CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace)
                    for (UnsignedInt PosY = SimulationSpaceSectorBoundsObject.StartYPos; PosY < SimulationSpaceSectorBoundsObject.EndYPos; PosY += CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace)
                        for (UnsignedInt PosZ = SimulationSpaceSectorBoundsObject.StartZPos; PosZ < SimulationSpaceSectorBoundsObject.EndZPos; PosZ += CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << SimulationSpaceSectorBoundsObject.StartXPos << " YStart = " << SimulationSpaceSectorBoundsObject.StartYPos << " ZStart = " << SimulationSpaceSectorBoundsObject.StartZPos << " XEnd = " << SimulationSpaceSectorBoundsObject.EndXPos << " YEnd = " << SimulationSpaceSectorBoundsObject.EndYPos << " ZEnd = " << SimulationSpaceSectorBoundsObject.EndZPos << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ));

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpace(true, PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
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
                LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));

                for (UnsignedInt ParticleSectorXIndex = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex < ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex < ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex < ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " XEnd = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpace(true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
                        }

                SyncPoint->arrive_and_wait();

                for (UnsignedInt ParticleSectorXIndex = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex < ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex < ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex < ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " XEnd = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, true);
                        }
            }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

inline UnsignedInt StepToChangeSimulationSpaceDivisionForThreads(const UnsignedInt StepOutside, bool StateOfSimulationSpaceDivisionForThreads)
{
    return ((StepOutside % CellEngineConfigDataObject.StepToChangeSpaceDivisionForThreads == 0) ? !StateOfSimulationSpaceDivisionForThreads : StateOfSimulationSpaceDivisionForThreads);
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread(barrier<>* SyncPoint, bool* StateOfSimulationSpaceDivisionForThreads, const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam) const
{
    try
    {
        for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
        {
            SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam, *StateOfSimulationSpaceDivisionForThreads, SyncPoint);

            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
            {
                SyncPoint->arrive_and_wait();

                if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert || CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                {
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert)
                        SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->ExchangeParticlesBetweenThreadsParallelInsert(StepOutside, *StateOfSimulationSpaceDivisionForThreads, false);
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                        SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->ExchangeParticlesBetweenThreadsParallelExtract(StepOutside, *StateOfSimulationSpaceDivisionForThreads, false);

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
                SyncPoint->arrive_and_wait();
        }
    }
    CATCH("generating n steps of simulation for whole cell space in one thread")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInThreads(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = true;

        LoggersManagerObject.Log(STREAM("MaxParticleIndex = " << MaxParticleIndex));

        SetZeroForAllParallelExecutionVariables();

        bool StateOfSimulationSpaceDivisionForThreads = false;

        std::barrier SyncPoint(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation);

        vector<vector<vector<thread*>>> Threads(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<thread*>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<thread*>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        LoggersManagerObject.Log(STREAM("START THREADS"));

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        ThreadIdType ThreadIndex = 1;
        FOR_EACH_THREAD_IN_XYZ
        {
            Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1] = new thread(&CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInOneThread, this, &SyncPoint, &StateOfSimulationSpaceDivisionForThreads, NumberOfStepsOutside, NumberOfStepsInside, ThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex);
            ThreadIndex++;
        }

        FOR_EACH_THREAD_IN_XYZ
        {
            Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->join();
            delete Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1];
        }

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));

        LoggersManagerObject.Log(STREAM("END THREADS"));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex)
{
    try
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
            {
                LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));

                for (UnsignedInt ParticleSectorXIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos; ParticleSectorXIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos; ParticleSectorYIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos; ParticleSectorZIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos; ParticleSectorZIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << " YStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << " ZStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << " XEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << " YEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << " ZEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, false);
                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                                GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, true);
                        }

                MPI_Barrier(MPI_COMM_WORLD);

                for (UnsignedInt ParticleSectorXIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos; ParticleSectorXIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos; ParticleSectorYIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos; ParticleSectorZIndex < CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos; ParticleSectorZIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartXPos << " YStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartYPos << " ZStart = " << CurrentMPIProcessSimulationSpaceSectorsRanges.StartZPos << " XEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndXPos << " YEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndYPos << " ZEnd = " << CurrentMPIProcessSimulationSpaceSectorsRanges.EndZPos << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));

                            if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                                GenerateOneStepOfDiffusionForSelectedSpace(true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
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
            GenerateOneStepOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam);
    }
    CATCH("generating n steps of simulation for whole cell space in one thread")
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        CellEngineConfigDataObject.MultiThreaded = true;

        LoggersManagerObject.Log(STREAM("MaxParticleIndex = " << MaxParticleIndex));

        SetZeroForAllParallelExecutionVariables();

        LoggersManagerObject.Log(STREAM("START MPI SIMULATION"));

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        if (CellEngineConfigDataObject.OpenGLGraphicsSwitchedOff == false && MPIProcessDataObject.CurrentMPIProcessIndex == 0)
        {
            int ValueToSend = 1;
            MPI_Bcast(&ValueToSend, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }

        GenerateNStepsOfSimulationForWholeCellSpaceInMPIProcess(NumberOfStepsOutside, NumberOfStepsInside, MPIProcessDataObject.CurrentMPIProcessIndex, 0, 0, 0);

        ExchangeParticlesBetweenMPIProcessesVer2();

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));

        LoggersManagerObject.Log(STREAM("END MPI SIMULATION"));

        GatherAllParallelExecutionVariables();
    }
    CATCH("generating n steps simulation for whole cell space in mpi process")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup1()
{
    try
    {
        MPI_Request MPIMessageRequest = MPI_REQUEST_NULL;
        for (UnsignedInt NeighbourProcessIndex = 0; NeighbourProcessIndex < NumberOfAllNeighbours; NeighbourProcessIndex++)
            if (NeighbourProcessesIndexes[NeighbourProcessIndex] != -1)
            {
                if (VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].empty() == true)
                {
                    VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].emplace_back(MPIParticleSenderStruct{ 0, 0, static_cast<int>(MPIProcessDataObject.CurrentMPIProcessIndex), static_cast<int>(NeighbourProcessesIndexes[NeighbourProcessIndex]), 0, 0, 0, 0, 0, 0 });
                    LoggersManagerObject.Log(STREAM("MPI Process Index to EMPTY send = " << NeighbourProcessesIndexes[NeighbourProcessIndex] << " " << VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex][0].ReceiverProcessIndex << " " << NeighbourProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                }
                else
                    LoggersManagerObject.Log(STREAM("MPI Process Index to send = " << NeighbourProcessesIndexes[NeighbourProcessIndex] << " " << VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex][0].ReceiverProcessIndex << " " << NeighbourProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

                LoggersManagerObject.LogUnconditional(STREAM("MPI Process Length Message SEND = " << VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].size() << " " << VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].size() * sizeof(MPIParticleSenderStruct) << " " << NeighbourProcessesIndexes[NeighbourProcessIndex] << " " << VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex][0].ReceiverProcessIndex << " " << NeighbourProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

                char BufferToSend[1024 * 1024];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].size();
                LoggersManagerObject.LogUnconditional(STREAM("MPI PACKED SIZE = " << NumberOfPackedStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                for (const auto& ParticleToSendElement : VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex])
                {
                    MPI_Pack(&ParticleToSendElement.ParticleIndex, 1, MPI_UNSIGNED, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ParticleKindId, 1, MPI_UNSIGNED, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SenderProcessIndex, 1, MPI_INT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.ReceiverProcessIndex, 1, MPI_INT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SectorPos.X, 1, MPI_UNSIGNED_SHORT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SectorPos.Y, 1, MPI_UNSIGNED_SHORT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.SectorPos.Z, 1, MPI_UNSIGNED_SHORT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewPosition.X, 1, MPI_FLOAT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewPosition.Y, 1, MPI_FLOAT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                    MPI_Pack(&ParticleToSendElement.NewPosition.Z, 1, MPI_FLOAT, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                }
                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, static_cast<int>(NeighbourProcessesIndexes[NeighbourProcessIndex]), 0, MPI_COMM_WORLD);

                VectorOfParticlesToSendToNeighbourProcesses[NeighbourProcessIndex].clear();
            }

        int NumberOfReceivedMessages = 0;
        while (NumberOfReceivedMessages < NumberOfActiveNeighbours)
        {
            MPI_Status MPIMessageStatus;

            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int NumberOfBytesReceived;
            MPI_Get_count(&MPIMessageStatus, MPI_UNSIGNED, &NumberOfBytesReceived);

            vector<UniqueIdInt> ReceivedConfirmationOfParticlesToRemove;

            char ReceivedParticlesToInsert1[1024 * 1024];
            MPI_Recv(&ReceivedParticlesToInsert1, 1024 * 1024, MPI_PACKED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int PositionInBuffer = 0;
            unsigned int NumberOfUnpackedParticleStructures;
            MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &NumberOfUnpackedParticleStructures, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
            LoggersManagerObject.LogUnconditional(STREAM("MPI UNPACKED SIZE TO REMOVE = " << NumberOfUnpackedParticleStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
            for (UnsignedInt UnpackedParticleStructureIndex = 0; UnpackedParticleStructureIndex < NumberOfUnpackedParticleStructures; UnpackedParticleStructureIndex++)
            {
                UniqueIdInt MPIParticleSenderStructElementLocalObject;
                MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                ReceivedConfirmationOfParticlesToRemove.emplace_back(MPIParticleSenderStructElementLocalObject);
            }

            if (ReceivedConfirmationOfParticlesToRemove[0] != 0)
                for (const auto& ParticleToRemoveConfirmedIndex : ReceivedConfirmationOfParticlesToRemove)
                    RemoveParticle(ParticleToRemoveConfirmedIndex, true);

            NumberOfReceivedMessages++;
        }
    }
    CATCH("exchange particles between mpi processes")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesGroup2Ver2()
{
    try
    {
        int NumberOfReceivedMessages = 0;

        vector<MPIParticleSenderStruct> ReceivedParticlesToInsertFromAllNeigbhours[NumberOfAllNeighbours];

        LoggersManagerObject.LogUnconditional(STREAM("NumberOfActiveNeighbours = " << NumberOfActiveNeighbours << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

        while (NumberOfReceivedMessages < NumberOfActiveNeighbours)
        {
            MPI_Status MPIMessageStatus;

            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

            int NumberOfBytesReceived;
            MPI_Get_count(&MPIMessageStatus, MPI_CHAR, &NumberOfBytesReceived);

            {
                LoggersManagerObject.LogUnconditional(STREAM("MPI Process Length Message RECEIVE = " << NumberOfBytesReceived << " " << MPIProcessDataObject.CurrentMPIProcessIndex << " " << NumberOfBytesReceived / sizeof(MPIParticleSenderStruct)));

                char ReceivedParticlesToInsert1[1024 * 1024];
                if (NumberOfBytesReceived > 1024 * 1024)
                    LoggersManagerObject.LogUnconditional(STREAM("ERROR NumberOfBytesReceived = " << NumberOfBytesReceived << " " << MPIProcessDataObject.CurrentMPIProcessIndex));

                MPI_Recv(&ReceivedParticlesToInsert1, 1024 * 1024, MPI_PACKED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &MPIMessageStatus);

                int PositionInBuffer = 0;
                unsigned int NumberOfUnpackedParticleStructures;
                MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &NumberOfUnpackedParticleStructures, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                LoggersManagerObject.LogUnconditional(STREAM("MPI UNPACKED SIZE = " << NumberOfUnpackedParticleStructures << " " << MPIProcessDataObject.CurrentMPIProcessIndex));
                for (UnsignedInt UnpackedParticleStructureIndex = 0; UnpackedParticleStructureIndex < NumberOfUnpackedParticleStructures; UnpackedParticleStructureIndex++)
                {
                    MPIParticleSenderStruct MPIParticleSenderStructElementLocalObject;
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleIndex, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ParticleKindId, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SenderProcessIndex, 1, MPI_INT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.ReceiverProcessIndex, 1, MPI_INT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SectorPos.X, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SectorPos.Y, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.SectorPos.Z, 1, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewPosition.X, 1, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewPosition.Y, 1, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(ReceivedParticlesToInsert1, 1024 * 1024, &PositionInBuffer, &MPIParticleSenderStructElementLocalObject.NewPosition.Z, 1, MPI_FLOAT, MPI_COMM_WORLD);

                    if (MPIParticleSenderStructElementLocalObject.ParticleIndex != 0)
                    {
                        bool Found = false;
                        UnsignedInt NeighbourProcessIndex;
                        for (NeighbourProcessIndex = 0; NeighbourProcessIndex < NumberOfAllNeighbours; NeighbourProcessIndex++)
                            if (NeighbourProcessesIndexes[NeighbourProcessIndex] == MPIParticleSenderStructElementLocalObject.SenderProcessIndex)
                            {
                                ReceivedParticlesToInsertFromAllNeigbhours[NeighbourProcessIndex].emplace_back(MPIParticleSenderStructElementLocalObject);
                                Found = true;
                                break;
                            }

                        if (Found == false)
                            cout << "GET PROCESS INDEX FROM = " << NeighbourProcessIndex << " " << NeighbourProcessesIndexes[NeighbourProcessIndex] << " " << MPIParticleSenderStructElementLocalObject.SenderProcessIndex << " TO = " << MPIParticleSenderStructElementLocalObject.ReceiverProcessIndex << " " << MPIProcessDataObject.CurrentMPIProcessIndex << endl;
                    }
                }

                NumberOfReceivedMessages++;
            }
        }

        for (const auto& ReceivedParticlesToInsert : ReceivedParticlesToInsertFromAllNeigbhours)
            if (ReceivedParticlesToInsert.empty() == false)
            {
                vector<UniqueIdInt> ConfirmationOfParticlesToRemoveToSent;

                for (const auto& ReceivedParticleIndexToInsert : ReceivedParticlesToInsert)
                    if (CheckInsertOfParticle(ReceivedParticleIndexToInsert) == true)
                        ConfirmationOfParticlesToRemoveToSent.emplace_back(ReceivedParticleIndexToInsert.ParticleIndex);

                if (ConfirmationOfParticlesToRemoveToSent.empty() == true)
                    ConfirmationOfParticlesToRemoveToSent.emplace_back(0);

                char BufferToSend[1024 * 1024];
                int PositionInBuffer = 0;
                unsigned int NumberOfPackedStructures = ConfirmationOfParticlesToRemoveToSent.size();

                MPI_Pack(&NumberOfPackedStructures, 1, MPI_UNSIGNED, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                for (const auto& ParticleToSendElement : ConfirmationOfParticlesToRemoveToSent)
                    MPI_Pack(&ParticleToSendElement, 1, MPI_UNSIGNED, BufferToSend, 1024 * 1024, &PositionInBuffer, MPI_COMM_WORLD);
                MPI_Send(BufferToSend, PositionInBuffer, MPI_PACKED, ReceivedParticlesToInsert[0].SenderProcessIndex, 0, MPI_COMM_WORLD);
            }
    }
    CATCH("exchange particles between mpi processes ver 2")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenMPIProcessesVer2()
{
    try
    {
        if (ProcessGroupNumber == 0)
        {
            ExchangeParticlesBetweenMPIProcessesGroup1();

            MPI_Barrier(MPI_COMM_WORLD);

            ExchangeParticlesBetweenMPIProcessesGroup2Ver2();

            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            ExchangeParticlesBetweenMPIProcessesGroup2Ver2();

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
        ErrorCounter = 0;
        NumberOfExecutedReactions = 0;
        NumberOfCancelledReactions = 0;
        NumberOfCancelledAReactions = 0;
        NumberOfCancelledBReactions = 0;
        AddedParticlesInReactions = 0;
        RemovedParticlesInReactions = 0;
        RestoredParticlesInCancelledReactions = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfExecutedReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledAReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledBReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RemovedParticlesInReactions = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RestoredParticlesInCancelledReactions = 0;
        }
    }
    CATCH("setting zero for all parallel execution variables")
}

void CellEngineSimulationParallelExecutionManager::GatherAllParallelExecutionVariables()
{
    try
    {
        FOR_EACH_THREAD_IN_XYZ
        {
            ErrorCounter += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter;
            NumberOfExecutedReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfExecutedReactions;
            NumberOfCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledReactions;
            NumberOfCancelledAReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledAReactions;
            NumberOfCancelledBReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->NumberOfCancelledBReactions;
            AddedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions;
            RemovedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RemovedParticlesInReactions;
            RestoredParticlesInCancelledReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->RestoredParticlesInCancelledReactions;
        }

        LoggersManagerObject.Log(STREAM("RestoredParticlesInCancelledReactions = " + to_string(RestoredParticlesInCancelledReactions)));
        LoggersManagerObject.Log(STREAM("RemovedParticlesInReactions = " + to_string(RemovedParticlesInReactions)));
        LoggersManagerObject.Log(STREAM("AddedParticlesInReactions = " + to_string(AddedParticlesInReactions)));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledAReactions = " + to_string(NumberOfCancelledAReactions)));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledBReactions = " + to_string(NumberOfCancelledBReactions)));
        LoggersManagerObject.Log(STREAM("NumberOfCancelledReactions = " + to_string(NumberOfCancelledReactions)));
        LoggersManagerObject.Log(STREAM("NumberOfExecutedReactions = " + to_string(NumberOfExecutedReactions)));
        LoggersManagerObject.Log(STREAM("ErrorCounter = " + to_string(ErrorCounter)));
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
                LoggersManagerObject.Log(STREAM("ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));

            if (ParticleObject.second.Center.X > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Y > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Z > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                LoggersManagerObject.Log(STREAM("Wrong Particle Center -> ParticleIndex = " << ParticleObject.first << " " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z  << " EntityId = " << ParticleObject.second.EntityId));
            if (ParticleObject.second.Center.X == 0 || ParticleObject.second.Center.Y == 0 || ParticleObject.second.Center.Z == 0)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ZERO -> ParticleIndex = " << ParticleObject.first << " " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z << " EntityId = " << ParticleObject.second.EntityId));
                NumberOfBadCenterParticles++;
            }
            if (ParticleObject.second.Center.X == 1 || ParticleObject.second.Center.Y == 1 || ParticleObject.second.Center.Z == 1)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ONE -> ParticleIndex = " << ParticleObject.first << " " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z << " EntityId = " << ParticleObject.second.EntityId));
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
