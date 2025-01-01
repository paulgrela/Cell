
#include <chrono>

#include "DateTimeUtils.h"

#include "CellEngineDataFile.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineSimulationParallelExecutionManager.h"

using namespace std;

void LogCenterOfParticleWithThreadIndex(const Particle& ParticleObject, const ThreadIdType ThreadXIndex, const ThreadIdType ThreadYIndex, const ThreadIdType ThreadZIndex)
{
    LoggersManagerObject.Log(STREAM("Center: " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << endl));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex << endl));
    LoggersManagerObject.Log(STREAM(endl));
}

CellEngineSimulationParallelExecutionManager::CellEngineSimulationParallelExecutionManager() : SimulationSpaceDataForThreads(CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer)
{
}

void CellEngineSimulationParallelExecutionManager::FirstSendParticlesForThreads(const bool PrintCenterOfParticleWithThreadIndex, const bool PrintTime)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (auto& SimulationSpaceDataForThreadsXPos : SimulationSpaceDataForThreads)
            for (auto& SimulationSpaceDataForThreadsYPos : SimulationSpaceDataForThreadsXPos)
                for (auto& SimulationSpaceDataForThreadsZPos : SimulationSpaceDataForThreadsYPos)
                    SimulationSpaceDataForThreadsZPos->ParticlesForThreads.clear();

        for (const auto& ParticleObject : Particles)
        {
            UnsignedInt ThreadXIndex = floor(ParticleObject.second.Center.X / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
            UnsignedInt ThreadYIndex = floor(ParticleObject.second.Center.Y / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
            UnsignedInt ThreadZIndex = floor(ParticleObject.second.Center.Z / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

            if (PrintCenterOfParticleWithThreadIndex == true)
                LogCenterOfParticleWithThreadIndex(ParticleObject.second, ThreadXIndex, ThreadYIndex, ThreadZIndex);

            SimulationSpaceDataForThreads[ThreadXIndex][ThreadYIndex][ThreadZIndex]->ParticlesForThreads.insert(ParticleObject);
        }

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    InitiateFreeParticleIndexes(SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads, false);

                    if (PrintCenterOfParticleWithThreadIndex == true)
                        LoggersManagerObject.Log(STREAM("THREAD[" << ThreadXIndex << "," << ThreadYIndex << "," << ThreadZIndex << "] SIZE = " << SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.size()));
                }

        const auto stop_time = chrono::high_resolution_clock::now();

        if (PrintTime == true)
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "First sending particles to threads has taken time = ","Execution in threads")));
    }
    CATCH("first sending particles for threads")
}

void PrintThreadIndexes(const Particle& ParticleObject, const ThreadIdType CurrentThreadIndex, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, const UnsignedInt ThreadXIndexNew, const UnsignedInt ThreadYIndexNew, const UnsignedInt ThreadZIndexNew)
{
    LoggersManagerObject.Log(STREAM("Particle Center: " << ParticleObject.Index << " " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId).ParticleKindSpecialDataSector[0].ParticleType)));
    LoggersManagerObject.Log(STREAM("THREAD POS NEW = " << ThreadXIndexNew << ThreadYIndexNew << ThreadZIndexNew));
    LoggersManagerObject.Log(STREAM("THREAD POS = " << ThreadXIndex << ", " << ThreadYIndex << ", " << ThreadZIndex));
    LoggersManagerObject.Log(STREAM(endl));
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreads(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        UnsignedInt ExchangedParticleCount = 0;

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
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
                            SimulationSpaceDataForThreads[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew]->ParticlesForThreads.insert(SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.extract(ParticleIter++));
                            ExchangedParticleCount++;
                        }
                        else
                            ++ParticleIter;
                    }
                }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Exchanging particles between threads has taken time = ","Execution in threads")));

        if (PrintInfo == true)
        {
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(stop_time - start_time, "Exchanging particles between threads has taken time = ","Execution in threads")));
            LoggersManagerObject.Log(STREAM("Number of exchanged particles: " << ExchangedParticleCount));
            LoggersManagerObject.Log(STREAM("StateOfSimulationSpaceDivisionForThreads: " << StateOfSimulationSpaceDivisionForThreads));
        }
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsParallelInsert(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCount = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<unordered_map<UniqueIdInt, Particle>>>> ParticlesToExchange(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<unordered_map<UniqueIdInt, Particle>>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<unordered_map<UniqueIdInt, Particle>>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject };

            auto ParticleIter = SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.begin();

            while (ParticleIter != SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.end())
            {
                const UnsignedInt SizeOfXInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2));
                const UnsignedInt SizeOfYInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2));
                const UnsignedInt SizeOfZInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2));

                if (StateOfSimulationSpaceDivisionForThreads == true)
                    if (ParticleIter->second.Center.X < SizeOfXInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Y < SizeOfYInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Z < SizeOfZInOneThreadInSimulationSpaceDiv2)
                    {
                        ++ParticleIter;
                        continue;
                    }

                UnsignedInt ThreadXIndexNew = floor((ParticleIter->second.Center.X - SizeOfXInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
                UnsignedInt ThreadYIndexNew = floor((ParticleIter->second.Center.Y - SizeOfYInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
                UnsignedInt ThreadZIndexNew = floor((ParticleIter->second.Center.Z - SizeOfZInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if ((CurrentThreadPos.ThreadPosX - 1 != ThreadXIndexNew) || (CurrentThreadPos.ThreadPosY - 1 != ThreadYIndexNew) || (CurrentThreadPos.ThreadPosZ - 1 != ThreadZIndexNew))
                {
                    ParticlesToExchange[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew].insert(SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.extract(ParticleIter++));
                    ExchangedParticleCount++;
                }
                else
                    ++ParticleIter;
            }
        }

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                    if (ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].empty() == false)
                    {
                        lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->MainExchangeParticlesMutexObject };

                        for (const auto& ParticleToExchange : ParticlesToExchange[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1])
                            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.insert(ParticleToExchange);
                    }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Exchanging particles between threads has taken time = ","Execution in threads")));

        if (PrintInfo == true)
        {
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(stop_time - start_time, "Exchanging particles between threads has taken time = ","Execution in threads")));
            LoggersManagerObject.Log(STREAM("Number of exchanged particles: " << ExchangedParticleCount));
            LoggersManagerObject.Log(STREAM("StateOfSimulationSpaceDivisionForThreads: " << StateOfSimulationSpaceDivisionForThreads));
        }
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::ExchangeParticlesBetweenThreadsParallelExtract(const UnsignedInt StepOutside, const bool StateOfSimulationSpaceDivisionForThreads, const bool PrintInfo) const
{
    try
    {
        UnsignedInt ExchangedParticleCount = 0;

        const auto start_time = chrono::high_resolution_clock::now();

        vector<vector<vector<unordered_map<UniqueIdInt, Particle>>>> ParticlesToExchangeMap(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<unordered_map<UniqueIdInt, Particle>>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<unordered_map<UniqueIdInt, Particle>>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        {
            lock_guard LockGuardObject1{ SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject };

            auto ParticleIter = SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.begin();

            while (ParticleIter != SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.end())
            {
                const UnsignedInt SizeOfXInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2));
                const UnsignedInt SizeOfYInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2));
                const UnsignedInt SizeOfZInOneThreadInSimulationSpaceDiv2 = (StateOfSimulationSpaceDivisionForThreads == false ? 0 : (CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2));

                if (StateOfSimulationSpaceDivisionForThreads == true)
                    if (ParticleIter->second.Center.X < SizeOfXInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Y < SizeOfYInOneThreadInSimulationSpaceDiv2 || ParticleIter->second.Center.Z < SizeOfZInOneThreadInSimulationSpaceDiv2)
                    {
                        ++ParticleIter;
                        continue;
                    }

                UnsignedInt ThreadXIndexNew = floor((ParticleIter->second.Center.X - SizeOfXInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace);
                UnsignedInt ThreadYIndexNew = floor((ParticleIter->second.Center.Y - SizeOfYInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace);
                UnsignedInt ThreadZIndexNew = floor((ParticleIter->second.Center.Z - SizeOfZInOneThreadInSimulationSpaceDiv2) / CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace);

                if ((CurrentThreadPos.ThreadPosX - 1 != ThreadXIndexNew) || (CurrentThreadPos.ThreadPosY - 1 != ThreadYIndexNew) || (CurrentThreadPos.ThreadPosZ - 1 != ThreadZIndexNew))
                {
                    ParticlesToExchangeMap[ThreadXIndexNew][ThreadYIndexNew][ThreadZIndexNew].insert(SimulationSpaceDataForThreads[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->ParticlesForThreads.extract(ParticleIter++));
                    ExchangedParticleCount++;
                }
                else
                    ++ParticleIter;
            }
        }

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                    if (ParticlesToExchangeMap[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].empty() == false)
                    {
                        lock_guard LockGuardObject2{ SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->MainExchangeParticlesMutexObject };

                        auto ParticleIterLocalInside = ParticlesToExchangeMap[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].begin();
                        while (ParticleIterLocalInside != ParticlesToExchangeMap[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].end())
                            SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads.insert(ParticlesToExchangeMap[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1].extract(ParticleIterLocalInside++));
                    }

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Exchanging particles between threads has taken time = ","Execution in threads")));

        if (PrintInfo == true)
        {
            LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(stop_time - start_time, "Exchanging particles between threads has taken time = ","Execution in threads")));
            LoggersManagerObject.Log(STREAM("Number of exchanged particles: " << ExchangedParticleCount));
            LoggersManagerObject.Log(STREAM("StateOfSimulationSpaceDivisionForThreads: " << StateOfSimulationSpaceDivisionForThreads));
        }
    }
    CATCH("exchanging particles between threads")
}

void CellEngineSimulationParallelExecutionManager::GatherParticlesFromThreadsToParticlesInMainThread()
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        Particles.clear();

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                    for (const auto& ParticleObject : SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ParticlesForThreads)
                        Particles.insert(ParticleObject);

        InitiateFreeParticleIndexes(Particles, false);

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Gathering particles from threads to main thread has taken time = ","Execution in threads")));
    }
    CATCH("sending particles for threads")
}

void CellEngineSimulationParallelExecutionManager::GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(const UnsignedInt NumberOfStepsInside, const UnsignedInt StepOutside, const UnsignedInt ThreadXIndex, const UnsignedInt ThreadYIndex, const UnsignedInt ThreadZIndex, bool StateOfSimulationSpaceDivisionForThreads)
{
    try
    {
        for (UnsignedInt Step2 = 1; Step2 <= NumberOfStepsInside; Step2++)
        {
            LoggersManagerObject.Log(STREAM("STEP INSIDE = " << Step2 << " ThreadX = " << ThreadXIndex << " ThreadX = " << ThreadYIndex << " ThreadX = " << ThreadZIndex));

            UnsignedInt XStartParam = (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace;
            UnsignedInt XEndParam = (ThreadXIndex - 1) * CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace + CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace;
            UnsignedInt YStartParam = (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace;
            UnsignedInt YEndParam = (ThreadYIndex - 1) * CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace + CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace;
            UnsignedInt ZStartParam = (ThreadZIndex - 1) * CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace;
            UnsignedInt ZEndParam = (ThreadZIndex - 1) * CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace + CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace;

            if (StateOfSimulationSpaceDivisionForThreads == true)
            {
                XStartParam += CellEngineConfigDataObject.SizeOfXInOneThreadInSimulationSpace / 2;
                YStartParam += CellEngineConfigDataObject.SizeOfYInOneThreadInSimulationSpace / 2;
                ZStartParam += CellEngineConfigDataObject.SizeOfZInOneThreadInSimulationSpace / 2;
            }

            for (UnsignedInt PosX = XStartParam; PosX < XEndParam; PosX += CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace)
                for (UnsignedInt PosY = YStartParam; PosY < YEndParam; PosY += CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace)
                    for (UnsignedInt PosZ = ZStartParam; PosZ < ZEndParam; PosZ += CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace)
                    {
                        LoggersManagerObject.Log(STREAM("XStart = " << XStartParam << " YStart = " << YStartParam << " ZStart = " << ZStartParam << " XEnd = " << XEndParam << " YEnd = " << YEndParam << " ZEnd = " << ZEndParam << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ));

                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion, CellEngineConfigData::TypesOfSimulation::OnlyDiffusion }))
                            GenerateOneStepOfDiffusionForSelectedSpace(true, PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace);
                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::BothReactionsAndDiffusion }))
                            GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, false);
                        if (CellEngineUseful::IsIn(CellEngineConfigDataObject.TypeOfSimulation, { CellEngineConfigData::TypesOfSimulation::OnlyReactions }))
                            GenerateOneRandomReactionForSelectedSpace(PosX, PosY, PosZ, CellEngineConfigDataObject.SizeOfXInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfYInOneSectorInOneThreadInSimulationSpace, CellEngineConfigDataObject.SizeOfZInOneSectorInOneThreadInSimulationSpace, true);
                    }
        }
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}

inline UnsignedInt StepToChangeSimulationSpaceDivisionForThreads(const UnsignedInt StepOutside, bool StateOfSimulationSpaceDivisionForThreads)
{
    return ((StepOutside % CellEngineConfigDataObject.StepToChangeSpaceDivisionForThreads == 0) ? !StateOfSimulationSpaceDivisionForThreads : StateOfSimulationSpaceDivisionForThreads);
}

void CellEngineSimulationParallelExecutionManager::GenerateNStepsOfSimulationForWholeCellSpaceInThreads(const UnsignedInt NumberOfStepsOutside, const UnsignedInt NumberOfStepsInside)
{
    try
    {
        CheckParticlesCenters(false);

        LoggersManagerObject.Log(STREAM("MaxParticleIndex = " << MaxParticleIndex));

        ErrorCounter = 0;
        AddedParticlesInReactions = 0;

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter = 0;
                    SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions = 0;
                }

        bool StateOfSimulationSpaceDivisionForThreads = false;

        std::barrier SyncPoint(CellEngineConfigDataObject.NumberOfXThreadsInSimulation * CellEngineConfigDataObject.NumberOfYThreadsInSimulation * CellEngineConfigDataObject.NumberOfZThreadsInSimulation);

        vector<vector<vector<thread*>>> Threads(CellEngineConfigDataObject.NumberOfXThreadsInSimulation, vector<vector<thread*>>(CellEngineConfigDataObject.NumberOfYThreadsInSimulation, vector<thread*>(CellEngineConfigDataObject.NumberOfZThreadsInSimulation)));

        auto GenerateNStepsOfSimulationForWholeCellSpaceInThreadsLambda = [NumberOfStepsOutside, &StateOfSimulationSpaceDivisionForThreads, &SyncPoint, this](const UnsignedInt NumberOfStepsInside, const ThreadIdType CurrentThreadIndexParam, const UnsignedInt ThreadXIndexParam, const UnsignedInt ThreadYIndexParam, const UnsignedInt ThreadZIndexParam)
        {
            for (UnsignedInt StepOutside = 1; StepOutside <= NumberOfStepsOutside; StepOutside++)
            {
                SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->GenerateOneStepOfSimulationForWholeCellSpaceInOneThread(NumberOfStepsInside, StepOutside, ThreadXIndexParam, ThreadYIndexParam, ThreadZIndexParam, StateOfSimulationSpaceDivisionForThreads);

                SyncPoint.arrive_and_wait();

                if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert || CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                {
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert)
                        SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->ExchangeParticlesBetweenThreadsParallelInsert(StepOutside, StateOfSimulationSpaceDivisionForThreads, false);
                    if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::ParallelExtract)
                        SimulationSpaceDataForThreads[ThreadXIndexParam - 1][ThreadYIndexParam - 1][ThreadZIndexParam - 1]->ExchangeParticlesBetweenThreadsParallelExtract(StepOutside, StateOfSimulationSpaceDivisionForThreads, false);

                    SyncPoint.arrive_and_wait();

                    if (CurrentThreadIndexParam == 1)
                        StateOfSimulationSpaceDivisionForThreads = StepToChangeSimulationSpaceDivisionForThreads(StepOutside, StateOfSimulationSpaceDivisionForThreads);
                }
                else
                if (CellEngineConfigDataObject.TypeOfExchangeOfParticlesBetweenThreads == CellEngineConfigData::TypesOfExchangeOfParticlesBetweenThreads::InMainThread && CurrentThreadIndexParam == 1)
                {
                    ExchangeParticlesBetweenThreads(StepOutside, StateOfSimulationSpaceDivisionForThreads, false);

                                                                                                                        SyncPoint.arrive_and_wait();

                    StateOfSimulationSpaceDivisionForThreads = StepToChangeSimulationSpaceDivisionForThreads(StepOutside, StateOfSimulationSpaceDivisionForThreads);
                }

                                                                                                                        SyncPoint.arrive_and_wait();
            }
        };

        LoggersManagerObject.Log(STREAM("START THREADS"));

        const auto start_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOffLogs();

        ThreadIdType ThreadIndex = 1;
        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    Threads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1] = new thread(GenerateNStepsOfSimulationForWholeCellSpaceInThreadsLambda, NumberOfStepsInside, ThreadIndex, ThreadXIndex, ThreadYIndex, ThreadZIndex);
                    ThreadIndex++;
                }

        for (auto& ThreadX : Threads)
            for (auto& ThreadY : ThreadX)
                for (auto& ThreadZ : ThreadY)
                {
                    ThreadZ->join();
                    delete ThreadZ;
                }

        CellEngineUseful::SwitchOnLogs();

        const auto stop_time = chrono::high_resolution_clock::now();

        string ResultText = "Execution in threads for steps outside = " + to_string(NumberOfStepsOutside) + " and steps inside = " + to_string(NumberOfStepsInside) + " has taken time: ";
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, ResultText.c_str(),"Execution in threads")));

        LoggersManagerObject.Log(STREAM("END THREADS"));

        for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
            for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
                for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
                {
                    ErrorCounter += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->ErrorCounter;
                    AddedParticlesInReactions += SimulationSpaceDataForThreads[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->AddedParticlesInReactions;
                }

        LoggersManagerObject.Log(STREAM("AddedParticlesInReactions =  " + to_string(AddedParticlesInReactions)));
        LoggersManagerObject.Log(STREAM("ErrorCounter = " + to_string(ErrorCounter)));
    }
    CATCH("generating n steps simulation for whole cell space in threads")
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

void CellEngineSimulationParallelExecutionManager::CheckParticlesCenters(const bool PrintAllParticles) const
{
    try
    {
        UnsignedInt NumberOfZeroSizedParticles = 0;
        UnsignedInt NumberOfZeroCenterParticles = 0;

        LoggersManagerObject.Log(STREAM("Number of erased particles with bad center = " << erase_if(Particles, [](auto& P){ return P.second.Center.X == 0 || P.second.Center.Y == 0 || P.second.Center.Z == 0; })));

        for (const auto& ParticleObject : Particles)
        {
            if (PrintAllParticles ==  true)
                LoggersManagerObject.Log(STREAM("ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));

            if (ParticleObject.second.Center.X > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Y > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension || ParticleObject.second.Center.Z > CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                LoggersManagerObject.Log(STREAM("Wrong Particle Center -> ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));
            if (ParticleObject.second.Center.X == 0 || ParticleObject.second.Center.Y == 0 || ParticleObject.second.Center.Z == 0)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ZERO -> ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));
                NumberOfZeroCenterParticles++;
            }
            if (ParticleObject.second.Center.X == 1 || ParticleObject.second.Center.Y == 1 || ParticleObject.second.Center.Z == 1)
            {
                LoggersManagerObject.Log(STREAM("Wrong Particle Center ONE -> ParticleIndex = " << ParticleObject.second.Index << " " << ParticleObject.second.Center.X << " " << ParticleObject.second.Center.Y << " " << ParticleObject.second.Center.Z));
                NumberOfZeroCenterParticles++;
            }
            if (ParticleObject.second.ListOfVoxels.empty() == true)
                NumberOfZeroSizedParticles++;
        }
        LoggersManagerObject.Log(STREAM("All Particles Centers Checked. Number Of Zero Center Particles = " << NumberOfZeroCenterParticles));
        LoggersManagerObject.Log(STREAM("All Particles Size Checked. Number Of Zero Sized Particles = " << NumberOfZeroSizedParticles));
    }
    CATCH("generating n steps simulation for whole cell space in threads")
}
