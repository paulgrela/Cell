
#include "CellEngineMacros.h"

#include "CellEngineDataFile.h"
#include "CellEngineParticlesFullAtomOperations.h"

#include "CellEngineImGuiMenu.h"

constexpr bool PrintAdditionalInformation = false;
constexpr bool PrintAdditionalInformationToLogs = false;
constexpr bool PrintAdditionalInformationToFiles = true;
constexpr bool PrintAdditionalInformationToConsole = true;

void CellEngineParticlesFullAtomOperations::SetProperThreadIndexForEveryParticlesSector(ParticlesContainer<Particle>& ParticlesSectors)
{
    try
    {
        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            const ThreadPosType ThreadPos = { .ThreadPosX = static_cast<SignedInt>(ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation + 1), .ThreadPosY = static_cast<SignedInt>(ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation + 1), .ThreadPosZ = static_cast<SignedInt>(ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation + 1) };
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].ThreadPos = ThreadPos;
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ThreadPos.ThreadPosX - 1][ThreadPos.ThreadPosY - 1][ThreadPos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;

            #ifdef SIMULATION_DETAILED_DEBUG_LOG
            if constexpr (PrintAdditionalInformation == true)
                LoggersManagerObject.Log(STREAM("ThreadPos = " << ThreadPos.ThreadPosX << "'" << ThreadPos.ThreadPosY << "'" << ThreadPos.ThreadPosZ << " ProcessIndex = " << ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex));
            #endif
        }
    }
    CATCH("setting proper thread index for every particles sector")
}

static bool ExchangeParticleBetweenSectors(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt SectorPosX1, const SignedInt SectorPosY1, const SignedInt SectorPosZ1, const SignedInt SectorPosX2, const SignedInt SectorPosY2, const SignedInt SectorPosZ2)
{
    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
    {
        ListOfParticlesToChangeSectors.emplace_back(ParticleToBeMovedFromOneSectorToAnotherSector{ .ParticleIndex = ParticleObject.Index, .SenderSectorPos = SectorPosType{ .SectorPosX = SectorPosX1, .SectorPosY = SectorPosY1, .SectorPosZ = SectorPosZ1 }, .ReceiverSectorPos = SectorPosType{ .SectorPosX = SectorPosX2, .SectorPosY = SectorPosY2, .SectorPosZ = SectorPosZ2 }});

        #ifdef SIMULATION_DETAILED_DEBUG_LOG
        const auto ParticleFromSourceToMoveToTargetIterator = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.find(ListOfParticlesToChangeSectors.back().ParticleIndex);
        if (ParticleFromSourceToMoveToTargetIterator == ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.end())
            LoggersManagerObject.Log(STREAM("ERROR2 = " << ListOfParticlesToChangeSectors.back().ParticleIndex << " " << ParticleObject.Index));
        #endif
    }

    return false;
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighborProcessesIndexes, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ThreadPosType& CurrentThreadPos)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
            return;

        #ifdef SIMULATION_DETAILED_DEBUG_LOG
        const auto ParticleFromSourceToMoveToTargetIterator = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.find(ParticleObject.Index);
        if (ParticleFromSourceToMoveToTargetIterator == ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.end())
            LoggersManagerObject.Log(STREAM("ERROR1 = " << ParticleObject.Index));
        #endif

        MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
        {
            if (CellEngineConfigDataObject.MultiThreaded == true)
            {
                bool NewSectorNeighborThreadFound = false;

                if ((SectorPosX1 != SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 != SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 != SectorPosZ2))
                {
                    const auto Thread1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].ThreadPos;
                    const auto Thread2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].ThreadPos;

                    if (Thread1Pos != Thread2Pos)
                    {
                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                        {
                            DEBUGLOG(LoggersManagerObject.Log(STREAM("NeighborProcessIndex = " << NeighborProcessesIndexes[NeighborProcessIndex] << "ThreadPos = " << Thread2Pos.ThreadPosX << "'" << Thread2Pos.ThreadPosY << "'" << Thread2Pos.ThreadPosZ << " ProcessIndex = " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1));)

                            const ThreadIdType NeighborThreadIndexFoundInSourceThreadNeighborThreadIndexes = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->NeighborProcessesIndexes[NeighborProcessIndex];
                            const ThreadIdType TargetThreadIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->CurrentThreadIndex - 1;
                            if (NeighborThreadIndexFoundInSourceThreadNeighborThreadIndexes == TargetThreadIndex)
                            {
                                #ifdef SIMULATION_DETAILED_DEBUG_LOG
                                if (Thread1Pos != CurrentThreadPos)
                                {
                                    if constexpr(PrintAdditionalInformationToLogs == true)
                                        LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << "C urrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.Index));
                                    if constexpr(PrintAdditionalInformationToFiles == true)
                                        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.Index));
                                    if constexpr(PrintAdditionalInformationToConsole == true)
                                        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = " << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.Index));
                                }
                                else
                                {
                                    if constexpr(PrintAdditionalInformationToLogs == true)
                                        LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = "));
                                    if constexpr(PrintAdditionalInformationToFiles == true)
                                        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", " << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ","  << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.Index));
                                    if constexpr(PrintAdditionalInformationToConsole == true)
                                        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD - TO THREAD = (" << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << ") FROM THREAD = (" << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << ") CurrentThreadIndex = " << TargetThreadIndex << " CurrentProcessIndex = " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = (" << SectorPosX1 << "," << SectorPosY1 << "," << SectorPosZ1 << ") S2 = (" << SectorPosX2 << "," << SectorPosY2 << "," << SectorPosZ2 << ") P = (" << ParticleObject.Center.X << ", " << ParticleObject.Center.Y << ", " << ParticleObject.Center.Z << ") V = (" << VectorX << ", " << VectorY << ", "  << VectorZ << ") PSHIFT = (" << ParticleObject.Center.X + VectorX << ", " << ParticleObject.Center.Y + VectorY << ", " << ParticleObject.Center.Z + VectorZ << ") PARTCLE_INDEX = " << ParticleObject.Index));
                                }
                                #endif

                                CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].emplace_back(ParticleSenderStructMultiThreaded{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = 0, .ReceiverProcessIndex = 0, .SenderThreadPos = { .ThreadPosX = Thread1Pos.ThreadPosX, .ThreadPosY = Thread1Pos.ThreadPosY, .ThreadPosZ = Thread1Pos.ThreadPosZ }, .ReceiverThreadPos = { .ThreadPosX = Thread2Pos.ThreadPosX, .ThreadPosY = Thread2Pos.ThreadPosY, .ThreadPosZ = Thread2Pos.ThreadPosZ }, .SenderSectorPos{ .X = static_cast<uint16_t>(SectorPosX1), .Y = static_cast<uint16_t>(SectorPosY1), .Z = static_cast<uint16_t>(SectorPosZ1) }, .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }, .ParticleObject = ParticleObject });
                                //CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->VectorOfParticlesToSendToNeighborProcessesOrThreads[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = 0, .ReceiverProcessIndex = 0, .SenderThreadPos = { .ThreadPosX = Thread1Pos.ThreadPosX, .ThreadPosY = Thread1Pos.ThreadPosY, .ThreadPosZ = Thread1Pos.ThreadPosZ }, .ReceiverThreadPos = { .ThreadPosX = Thread2Pos.ThreadPosX, .ThreadPosY = Thread2Pos.ThreadPosY, .ThreadPosZ = Thread2Pos.ThreadPosZ }, .SenderSectorPos{ .X = static_cast<uint16_t>(SectorPosX1), .Y = static_cast<uint16_t>(SectorPosY1), .Z = static_cast<uint16_t>(SectorPosZ1) }, .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }, ParticleObject });
                                //CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->VectorOfParticlesToSendToNeighborProcessesOrThreads[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = 0, .ReceiverProcessIndex = 0, .SenderThreadPos = { .ThreadPosX = Thread1Pos.ThreadPosX, .ThreadPosY = Thread1Pos.ThreadPosY, .ThreadPosZ = Thread1Pos.ThreadPosZ }, .ReceiverThreadPos = { .ThreadPosX = Thread2Pos.ThreadPosX, .ThreadPosY = Thread2Pos.ThreadPosY, .ThreadPosZ = Thread2Pos.ThreadPosZ }, .SenderSectorPos{ .X = static_cast<uint16_t>(SectorPosX1), .Y = static_cast<uint16_t>(SectorPosY1), .Z = static_cast<uint16_t>(SectorPosZ1) }, .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }});
                                // CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->VectorOfWholeParticlesToSendToNeighborProcessesOrThreads[NeighborProcessIndex].emplace_back(ParticleObject);

                                NewSectorNeighborThreadFound = true;
                                break;
                            }
                        }
                    }
                    else
                    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
                    {
                        ListOfParticlesToChangeSectors.emplace_back(ParticleToBeMovedFromOneSectorToAnotherSector{ .ParticleIndex = ParticleObject.Index, .SenderSectorPos = SectorPosType{ .SectorPosX = SectorPosX1, .SectorPosY = SectorPosY1, .SectorPosZ = SectorPosZ1 }, .ReceiverSectorPos = SectorPosType{ .SectorPosX = SectorPosX2, .SectorPosY = SectorPosY2, .SectorPosZ = SectorPosZ2 }});
                        NewSectorNeighborThreadFound = true;
                    }
                }

                if (NewSectorNeighborThreadFound == false)
                {
                    MoveAllAtomsInParticleAtomsListByVector(ParticleObject, -VectorX, -VectorY, -VectorZ);
                    ParticleObject.SetCenterCoordinates(ParticleObject.Center.X - VectorX, ParticleObject.Center.Y - VectorY, ParticleObject.Center.Z - VectorZ);
                }
            }
            else
                ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ParticleObjectIter, ListOfParticlesToChangeSectors, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
        }
    }
    CATCH_AND_THROW("moving particle by vector for threads")
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForMPIProcesses(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStructMPIMultiProcess>* VectorOfParticlesToSendToNeighborProcesses, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ThreadPosType CurrentThreadPos)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
            return;

        MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
        {
            bool NewSectorNeighborProcessFound = false;

            if ((SectorPosX1 != SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 != SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 != SectorPosZ2))
            {
                const auto Process1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].MPIProcessIndex;
                const auto Process2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].MPIProcessIndex;

                if (Process1Pos != Process2Pos)
                {
                    for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                        if (NeighborProcessesIndexes[NeighborProcessIndex] == Process2Pos)
                        {
                            #ifdef SIMULATION_DETAILED_DEBUG_LOG
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

                            VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].emplace_back(ParticleSenderStructMPIMultiProcess{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = static_cast<int>(Process1Pos), .ReceiverProcessIndex = static_cast<int>(Process2Pos), .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }});
                            //VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = static_cast<int>(Process1Pos), .ReceiverProcessIndex = static_cast<int>(Process2Pos), .ReceiverSectorPos = { .X = static_cast<uint16_t>(SectorPosX2), .Y = static_cast<uint16_t>(SectorPosY2), .Z = static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { .X = ParticleObject.Center.X, .Y = ParticleObject.Center.Y, .Z = ParticleObject.Center.Z }});

                            NewSectorNeighborProcessFound = true;
                            break;
                        }
                }
                else
                {
                    ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ParticleObjectIter, ListOfParticlesToChangeSectors, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
                    NewSectorNeighborProcessFound = true;
                }
            }

            if (NewSectorNeighborProcessFound == false)
            {
                MoveAllAtomsInParticleAtomsListByVector(ParticleObject, -VectorX, -VectorY, -VectorZ);
                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X - VectorX, ParticleObject.Center.Y - VectorY, ParticleObject.Center.Z - VectorZ);
            }
        }
    }
    CATCH_AND_THROW("moving particle by vector for mpi processes")
}
