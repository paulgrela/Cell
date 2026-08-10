
#include "CellEngineMacros.h"

#include "CellEngineDataFile.h"
#include "CellEngineParticlesFullAtomOperations.h"

#include "CellEngineImGuiMenu.h"

constexpr bool PrintAdditionalInformation = false;
constexpr bool PrintAdditionalInformationToLogs = false;
constexpr bool PrintAdditionalInformationToFiles = true;
constexpr bool PrintAdditionalInformationToConsole = false;

void CellEngineParticlesFullAtomOperations::SetProperThreadIndexForEveryParticlesSector(ParticlesContainer<Particle>& ParticlesSectors)
{
    try
    {
        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            const ThreadPosType ThreadPos = { static_cast<SignedInt>(ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation + 1), static_cast<SignedInt>(ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation + 1), static_cast<SignedInt>(ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation + 1) };
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].ThreadPos = ThreadPos;
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ThreadPos.ThreadPosX - 1][ThreadPos.ThreadPosY - 1][ThreadPos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;

            if constexpr (PrintAdditionalInformation == true)
                LoggersManagerObject.Log(STREAM("ThreadPos = " << ThreadPos.ThreadPosX << "'" << ThreadPos.ThreadPosY << "'" << ThreadPos.ThreadPosZ << " ProcessIndex = " << ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex));
        }
    }
    CATCH("setting proper thread index for every particles sector")
}

static bool ExchangeParticleBetweenSectors(const Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, UnsignedInt const SectorPosX1, UnsignedInt const SectorPosY1, const UnsignedInt SectorPosZ1, const UnsignedInt SectorPosX2, const UnsignedInt SectorPosY2, const UnsignedInt SectorPosZ2)
{
    #ifdef CONTAINERS_FOR_SPEED
    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
    #else
    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObjectIter->second.Index) == false)
    #endif
    {
        #ifdef CONTAINERS_FOR_SPEED
        ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
        #else
        if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObjectIter->first) == true)
        {
            //ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObjectIter++));
            const auto ParticleIndexCopiedObject = ParticleObjectIter->first;
            const auto ParticleCopiedObject = ParticleObjectIter->second;
            ParticleObjectIter = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.erase(ParticleObjectIter);
            ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert({ ParticleIndexCopiedObject, ParticleCopiedObject });
        }
        else
            ++ParticleObjectIter;
        #endif
    }
    else
    {
        #ifdef CONTAINERS_FOR_SPEED
        #else
        if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObjectIter->second.Index) == false)
            if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObjectIter->second.Index)->second == &ParticleObject)
            {
                ++ParticleObjectIter;

                if constexpr(PrintAdditionalInformation == true)
                    std::cout << "Earlier correct sector for particle index: " << ParticleObject.EntityId << " " << ParticleObject.Index << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                return true;
            }
        ++ParticleObjectIter;
        #endif
    }

    return false;
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighborThreads, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ThreadPosType& CurrentThreadPos)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
        {
            #ifndef CONTAINERS_FOR_SPEED
            ++ParticleObjectIter;
            #endif
            return;
        }

        MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);
        ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
        {
            if (CellEngineConfigDataObject.MultiThreaded == true)
            {
                bool NewSectorNeighborThreadFound = false;

                if ((SectorPosX1 != SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 != SectorPosY2 && SectorPosZ1 == SectorPosZ2) || (SectorPosX1 == SectorPosX2 && SectorPosY1 == SectorPosY2 && SectorPosZ1 != SectorPosZ2))
                {
                    #ifdef CONTAINERS_FOR_SPEED
                    const auto Thread1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].ThreadPos;
                    const auto Thread2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].ThreadPos;

                    if (Thread1Pos != Thread2Pos)
                    {
                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                        {
                            LoggersManagerObject.Log(STREAM("NeighborProcessIndex = " << NeighborProcessesIndexes[NeighborProcessIndex] << "ThreadPos = " << Thread2Pos.ThreadPosX << "'" << Thread2Pos.ThreadPosY << "'" << Thread2Pos.ThreadPosZ << " ProcessIndex = " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1));

                            if (CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[CurrentThreadPos.ThreadPosX - 1][CurrentThreadPos.ThreadPosY - 1][CurrentThreadPos.ThreadPosZ - 1]->NeighborProcessesIndexes[NeighborProcessIndex] == CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1)
                            {
                                if (Thread1Pos != CurrentThreadPos)
                                {
                                    if constexpr(PrintAdditionalInformationToLogs == true)
                                        LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE BAD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                    if constexpr(PrintAdditionalInformationToFiles == true)
                                        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                    if constexpr(PrintAdditionalInformationToConsole == true)
                                        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                }
                                else
                                {
                                    if constexpr(PrintAdditionalInformationToLogs == true)
                                        LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                    if constexpr(PrintAdditionalInformationToFiles == true)
                                        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                    if constexpr(PrintAdditionalInformationToConsole == true)
                                        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                }

                                VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ .ParticleIndex = ParticleObject.Index, .ParticleKindId = ParticleObject.EntityId, .SenderProcessIndex = 0, .ReceiverProcessIndex = 0, .SectorPos = { static_cast<uint16_t>(SectorPosX2), static_cast<uint16_t>(SectorPosY2), static_cast<uint16_t>(SectorPosZ2) }, .NewPosition = { ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z }});

                                NewSectorNeighborThreadFound = true;
                                break;
                            }
                        }
                    }
                    else
                    if (Thread1Pos == CurrentThreadPos)
                        if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
                            ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
                    #endif
                }

                if (NewSectorNeighborThreadFound == false)
                {
                    MoveAllAtomsInParticleAtomsListByVector(ParticleObject, -VectorX, -VectorY, -VectorZ);
                    ParticleObject.SetCenterCoordinates(ParticleObject.Center.X - VectorX, ParticleObject.Center.Y - VectorY, ParticleObject.Center.Z - VectorZ);
                }
            }
            else
                ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ParticleObjectIter, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
        }
        #ifndef CONTAINERS_FOR_SPEED
        else
            ++ParticleObjectIter;
        #endif
    }
    CATCH_AND_THROW("moving particle by vector for threads")
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForMPIProcesses(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighborProcesses, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ThreadPosType CurrentThreadPos)
{
    try
    {
        auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
        auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

        if (SectorPosX2 == -1 || SectorPosY2 == -1 || SectorPosZ2 == -1)
        {
            #ifndef CONTAINERS_FOR_SPEED
            ++ParticleObjectIter;
            #endif
            return;
        }

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
                            if (Process1Pos != MPIProcessDataObject.CurrentMPIProcessIndex)
                                LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE BAD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                            else
                                LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Process2Pos << " FROM " << Process1Pos << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));

                            VectorOfParticlesToSendToNeighborProcesses[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ ParticleObject.Index, ParticleObject.EntityId, static_cast<int>(Process1Pos), static_cast<int>(Process2Pos), { static_cast<uint16_t>(SectorPosX2), static_cast<uint16_t>(SectorPosY2), static_cast<uint16_t>(SectorPosZ2) }, { ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z } });

                            NewSectorNeighborProcessFound = true;
                            break;
                        }
                }
                else
                {
                    ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, ParticleObjectIter, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
                    NewSectorNeighborProcessFound = true;
                }
            }

            if (NewSectorNeighborProcessFound == false)
            {
                MoveAllAtomsInParticleAtomsListByVector(ParticleObject, -VectorX, -VectorY, -VectorZ);
                ParticleObject.SetCenterCoordinates(ParticleObject.Center.X - VectorX, ParticleObject.Center.Y - VectorY, ParticleObject.Center.Z - VectorZ);
            }
        }
        #ifndef CONTAINERS_FOR_SPEED
        else
            ++ParticleObjectIter;
        #endif
    }
    CATCH_AND_THROW("moving particle by vector for mpi processes")
}
