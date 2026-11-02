
#include "CellEngineMacros.h"
#include "CellEngineDataFile.h"
#include "CellEngineMPIProcess.h"
#include "CellEngineParticlesFullAtomOperations.h"

constexpr bool PrintAdditionalInformation = false;
constexpr bool PrintAdditionalInformationToLogs = false;
constexpr bool PrintAdditionalInformationToFiles = true;
constexpr bool PrintAdditionalInformationToConsole = true;

static bool ExchangeParticleBetweenSectors(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt SectorPosX1, const SignedInt SectorPosY1, const SignedInt SectorPosZ1, const SignedInt SectorPosX2, const SignedInt SectorPosY2, const SignedInt SectorPosZ2)
{
    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.ParticleUniqueId) == false)
    {
        ListOfParticlesToChangeSectors.emplace_back(ParticleToBeMovedFromOneSectorToAnotherSector{ .ParticleUniqueId = ParticleObject.ParticleUniqueId, .SenderSectorPos = SectorPosType{ .PosX = SectorPosX1, .PosY = SectorPosY1, .PosZ = SectorPosZ1 }, .ReceiverSectorPos = SectorPosType{ .PosX = SectorPosX2, .PosY = SectorPosY2, .PosZ = SectorPosZ2 }});

        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        const auto ParticleFromSourceToMoveToTargetIterator = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.find(ListOfParticlesToChangeSectors.back().ParticleUniqueId);
        if (ParticleFromSourceToMoveToTargetIterator == ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.end())
            LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("ERROR2 = " << ListOfParticlesToChangeSectors.back().ParticleUniqueId << " " << ParticleObject.ParticleUniqueId));
        #endif
    }

    return false;
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ParallelTaskPosType& CurrentThreadPos)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
            return;

        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        const auto ParticleFromSourceToMoveToTargetIterator = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.find(ParticleObject.ParticleUniqueId);
        if (ParticleFromSourceToMoveToTargetIterator == ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.end())
            LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("ERROR1 = " << ParticleObject.ParticleUniqueId));
        #endif

        if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
        {
            if (CellEngineConfigDataObject.MultiThreaded == true)
            {
                if ((SectorPosX1 != SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 != SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 != SectorPosZ2))
                {
                    const auto Thread1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].ParallelTaskPos;
                    const auto Thread2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].ParallelTaskPos;

                    if (Thread1Pos != Thread2Pos)
                    {
                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                            if (CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX][Thread2Pos.PosY][Thread2Pos.PosZ]->NeighborParallelTasksIndexes[NeighborProcessIndex] != -1)
                            {
                                DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("NeighborProcessIndex = " << NeighborProcessIndex << " ThreadPos = " << Thread2Pos.PosX << "'" << Thread2Pos.PosY << "'" << Thread2Pos.PosZ << " ProcessIndex = " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX - 1][Thread2Pos.PosY - 1][Thread2Pos.PosZ - 1]->MPIProcessIndex - 1));)

                                const ThreadIdType NeighborThreadIndexFoundInSourceThreadNeighborThreadIndexes = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[CurrentThreadPos.PosX][CurrentThreadPos.PosY][CurrentThreadPos.PosZ]->NeighborParallelTasksIndexes[NeighborProcessIndex];
                                const ThreadIdType TargetThreadIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX][Thread2Pos.PosY][Thread2Pos.PosZ]->CurrentThreadIndex;

                                if (NeighborThreadIndexFoundInSourceThreadNeighborThreadIndexes == TargetThreadIndex)
                                {
                                    #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                                    if (Thread1Pos != CurrentThreadPos)
                                    {
                                        if constexpr(PrintAdditionalInformationToLogs == true)
                                            LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << "CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
                                        if constexpr(PrintAdditionalInformationToFiles == true)
                                            LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
                                        if constexpr(PrintAdditionalInformationToConsole == true)
                                            LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = " << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
                                    }
                                    else
                                    {
                                        if constexpr(PrintAdditionalInformationToLogs == true)
                                            LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = "));
                                        if constexpr(PrintAdditionalInformationToFiles == true)
                                            LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", "  << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
                                        if constexpr(PrintAdditionalInformationToConsole == true)
                                            LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.PosX << "," << Thread2Pos.PosY << "," << Thread2Pos.PosZ << ") FROM THREAD = (" << Thread1Pos.PosX << "," << Thread1Pos.PosY << "," << Thread1Pos.PosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", "  << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.ParticleUniqueId));
                                    }
                                    #endif

                                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("VectorOfParticlesToSendToNeighborParallelTasks = " << NeighborProcessIndex << " " << ParticleObject.ParticleUniqueId << " " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX - 1][Thread2Pos.PosY - 1][Thread2Pos.PosZ - 1]->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].size()));)

                                    CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX][Thread2Pos.PosY][Thread2Pos.PosZ]->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].emplace_back(ParticleSenderStructForParallelMultiTaskingExecution{ .ParticleUniqueId = ParticleObject.ParticleUniqueId, .ParticleKindId = ParticleObject.ParticleKindId, .SenderParallelTaskIndex = 0, .ReceiverParallelTaskIndex = 0, .SenderParallelTaskPos = { .PosX = Thread1Pos.PosX, .PosY = Thread1Pos.PosY, .PosZ = Thread1Pos.PosZ }, .ReceiverParallelTaskPos = { .PosX = Thread2Pos.PosX, .PosY = Thread2Pos.PosY, .PosZ = Thread2Pos.PosZ }, .SenderSectorPos{ .X = static_cast<uint16_t>(SectorPosX1), .Y = static_cast<uint16_t>(SectorPosY1), .Z = static_cast<uint16_t>(SectorPosZ1) }, .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewParticleCenter = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }, .Vector = { .X = VectorX, .Y = VectorY, .Z = VectorZ }, .ParticleObject = ParticleObject });

                                    CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread1Pos.PosX][Thread1Pos.PosY][Thread1Pos.PosZ]->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask++;

                                    DEBUG_LOG_DIFFUSION(LoggersManagerObject.LogUnconditional(STREAM("VectorOfParticlesToSendToNeighborParallelTasks = " << NeighborProcessIndex << " " << ParticleObject.ParticleUniqueId << " " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX - 1][Thread2Pos.PosY - 1][Thread2Pos.PosZ - 1]->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].back().ParticleUniqueId << " " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.PosX - 1][Thread2Pos.PosY - 1][Thread2Pos.PosZ - 1]->VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].size()));)

                                    break;
                                }
                            }
                    }
                    else
                    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.ParticleUniqueId) == false)
                    {
                        CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
                        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

                        ListOfParticlesToChangeSectors.emplace_back(ParticleToBeMovedFromOneSectorToAnotherSector{ .ParticleUniqueId = ParticleObject.ParticleUniqueId, .SenderSectorPos = SectorPosType{ .PosX = SectorPosX1, .PosY = SectorPosY1, .PosZ = SectorPosZ1 }, .ReceiverSectorPos = SectorPosType{ .PosX = SectorPosX2, .PosY = SectorPosY2, .PosZ = SectorPosZ2 }});
                    }
                }
            }
            else
            {
                CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

                ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
            }
        }
    }
    CATCH_AND_THROW("moving particle by vector for threads")
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForMPIProcesses(Particle& ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStructForParallelMultiTaskingExecution>* VectorOfParticlesToSendToNeighborParallelTasks, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
            return;

        #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
        const auto ParticleFromSourceToMoveToTargetIterator = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.find(ParticleObject.ParticleUniqueId);
        if (ParticleFromSourceToMoveToTargetIterator == ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.end())
            LoggersManagerObject.Log(STREAM("ERROR1 = " << ParticleObject.ParticleUniqueId));
        #endif

        if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
        {
            if ((SectorPosX1 != SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 != SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 != SectorPosZ2))
            {
                const auto Process1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].MPIProcessIndex;
                const auto Process2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].MPIProcessIndex;

                if (Process1Pos != Process2Pos)
                {
                    for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                        if (NeighborProcessesIndexes[NeighborProcessIndex] == Process2Pos)
                        {
                            #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
                            if (Process1Pos != MPIProcessDataObject.CurrentMPIProcessIndex)
                            {
                                if constexpr(PrintAdditionalInformationToLogs == true)
                                    LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE BAD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                                if constexpr(PrintAdditionalInformationToFiles == true)
                                    LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                                if constexpr(PrintAdditionalInformationToConsole == true)
                                    LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                            }
                            else
                            {
                                if constexpr(PrintAdditionalInformationToLogs == true)
                                    LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                                if constexpr(PrintAdditionalInformationToFiles == true)
                                    LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                                if constexpr(PrintAdditionalInformationToConsole == true)
                                    LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << " S2 = " << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << " P = " << ParticleObject.Center.X << "," << ParticleObject.Center.Y << "," << ParticleObject.Center.Z << " V = " << VectorX << "," << VectorY << "," << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << "," << ParticleObject.Center.Y + VectorY << "," << ParticleObject.Center.Z + VectorZ));
                            }
                            #endif

                            VectorOfParticlesToSendToNeighborParallelTasks[NeighborProcessIndex].emplace_back(ParticleSenderStructForParallelMultiTaskingExecution{ .ParticleUniqueId = ParticleObject.ParticleUniqueId, .ParticleKindId = ParticleObject.ParticleKindId, .SenderParallelTaskIndex = static_cast<int>(Process1Pos), .ReceiverParallelTaskIndex = static_cast<int>(Process2Pos), .SenderSectorPos{ .X = static_cast<uint16_t>(SectorPosX1), .Y = static_cast<uint16_t>(SectorPosY1), .Z = static_cast<uint16_t>(SectorPosZ1) }, .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewParticleCenter = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }, .Vector = { .X = VectorX, .Y = VectorY, .Z = VectorZ }, .ParticleObject = ParticleObject });

                            CellEngineDataFileObjectPointer->CellEngineFullAtomSimulationSpaceObjectPointer->NumberOfParticlesToBeMovedInDiffusionSentToCheckInAnotherParallelTask++;

                            break;
                        }
                }
                else
                {
                    CellEngineUseful::MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
                    ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

                    ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
                }
            }
        }
    }
    CATCH_AND_THROW("moving particle by vector for mpi processes")
}
