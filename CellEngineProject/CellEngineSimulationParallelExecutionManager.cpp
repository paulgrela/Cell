
#include <chrono>

#include "DateTimeUtils.h"

#include "CellEngineDataFile.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineSimulationParallelExecutionManager.h"

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

                    ThreadLocalParticlesInProximityZPos = std::make_shared<SimulationSpaceType>(Particles, false, ThreadIndexPos, CurrentThreadPosType{ ThreadXPos, ThreadYPos, ThreadZPos });
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

#define FOR_EACH_THREAD_IN_XYZ \
        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++) \
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++) \
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)

#ifdef SHORTER_CODE
void CellEngineSimulationParallelExecutionManager::JoinStatisticsFromThreads(vector<map<UnsignedInt, ReactionStatistics>>& SavedReactionsMap, const UnsignedInt SimulationStepNumber) const
{
    try
    {
        FOR_EACH_THREAD_IN_XYZ
            for (const auto& ReactionData : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap.back())
                SavedReactionsMap[SimulationStepNumber - 1][ReactionData.first].Counter += ReactionData.second.Counter;
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

        FOR_EACH_PARTICLE_IN_XYZ_CONST
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
            for (const auto& CancelledParticleIndex : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->CancelledParticlesIndexes)
                CancelledParticlesIndexes.insert(CancelledParticleIndex);
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

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads, barrier
    <>* SyncPoint)
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
                            // if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                            //     GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, false);
                            // if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                            //     GenerateOneRandomReactionForSelectedSpace(ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, true);
                        }

                SyncPoint->arrive_and_wait();

                for (UnsignedInt ParticleSectorXIndex = (ThreadXIndex - 1) * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex < ThreadXIndex * CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation; ParticleSectorXIndex++)
                    for (UnsignedInt ParticleSectorYIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex < ThreadYIndex * CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation; ParticleSectorYIndex++)
                        for (UnsignedInt ParticleSectorZIndex = (ThreadYIndex - 1) * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex < ThreadZIndex * CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation; ParticleSectorZIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("XStart = " << (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZStart = " << (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " XEnd = " << ThreadXIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " YEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace << " ZEnd = " << ThreadZIndex * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace << " PosX = " << ParticleSectorXIndex << " PosY = " << ParticleSectorYIndex << " PosZ = " << ParticleSectorZIndex));

                            // if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                            //     GenerateOneStepOfDiffusionForSelectedSpace(true, ParticleSectorXIndex, ParticleSectorYIndex, ParticleSectorZIndex, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
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

void CellEngineSimulationParallelExecutionManager::SetZeroForAllParallelExecutionVariables()
{
    try
    {
        ErrorCounter = 0;
        AddedParticlesInReactions = 0;

        FOR_EACH_THREAD_IN_XYZ
        {
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter = 0;
            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions = 0;
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
            AddedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions;
        }

        LoggersManagerObject.Log(STREAM("AddedParticlesInReactions = " + to_string(AddedParticlesInReactions)));
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

        FOR_EACH_PARTICLE_IN_XYZ_CONST
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
