
#include "CellEngineMacros.h"

#include "CellEngineDataFile.h"
#include "CellEngineParticlesFullAtomOperations.h"

#include "CellEngineImGuiMenu.h"

constexpr bool PrintInfoWarning = false;

void CellEngineParticlesFullAtomOperations::SetProperThreadIndexForEveryParticlesSector(ParticlesContainer<Particle>& ParticlesSectors)
{
    try
    {
        FOR_EACH_SECTOR_IN_XYZ_ONLY
        {
            const ThreadPosType ThreadPos = { static_cast<SignedInt>(ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation + 1), static_cast<SignedInt>(ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation + 1), static_cast<SignedInt>(ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation + 1) };
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].ThreadPos = ThreadPos;
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ThreadPos.ThreadPosX - 1][ThreadPos.ThreadPosY - 1][ThreadPos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;

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
        if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObject.Index) == false)
            if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObject.Index)->second == &ParticleObject)
            {
                if constexpr(PrintInfoWarning == true)
                    std::cout << "Earlier correct sector for particle index: " << ParticleObject.EntityId << " " << ParticleObject.Index << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                return true;
            }
        #else
        if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObjectIter->second.Index) == false)
            if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObjectIter->second.Index)->second == &ParticleObject)
            {
                ++ParticleObjectIter;

                if constexpr(PrintInfoWarning == true)
                    std::cout << "Earlier correct sector for particle index: " << ParticleObject.EntityId << " " << ParticleObject.Index << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                return true;
            }
        ++ParticleObjectIter;
        #endif
    }

    return false;
}

//void CellEngineParticlesFullAtomOperations::MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, const ThreadPosType* NeighborThreadsIndexes, const unique_ptr<ThreadsIndexesType>& ThreadsIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighborThreads, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const ThreadPosType& CurrentThreadPos)
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
                        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("C5"));

                        for (UnsignedInt NeighborProcessIndex = 0; NeighborProcessIndex < NumberOfAllNeighbors; NeighborProcessIndex++)
                        {
                            //LoggersManagerObject.Log(STREAM("ThreadPos = " << Thread2Pos.ThreadPosX << "'" << Thread2Pos.ThreadPosY << "'" << Thread2Pos.ThreadPosZ << " ProcessIndex = " << ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex));
                            LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("NeighborProcessIndex = " << NeighborProcessesIndexes[NeighborProcessIndex] << "ThreadPos = " << Thread2Pos.ThreadPosX << "'" << Thread2Pos.ThreadPosY << "'" << Thread2Pos.ThreadPosZ << " ProcessIndex = " << CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1));
                            //if (NeighborProcessesIndexes[NeighborProcessIndex] == CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1)
                            if (NeighborProcessesIndexes[NeighborProcessIndex] == CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1)
                                //if (NeighborThreadsIndexes[NeighborProcessIndex] == Thread2Pos)
                                    //if (NeighborThreadsIndexes[NeighborProcessIndex] == Thread2Pos)
                                        //if (NeighborThreadsIndexes[NeighborProcessIndex].ThreadPosX == Thread2Pos.ThreadPosX - 1 && NeighborThreadsIndexes[NeighborProcessIndex].ThreadPosY == Thread2Pos.ThreadPosY - 1 && NeighborThreadsIndexes[NeighborProcessIndex].ThreadPosZ == Thread2Pos.ThreadPosZ - 1)
                            {
                                LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("C6"));

                                if (Thread1Pos != CurrentThreadPos)
                                    LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE BAD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex << " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));
                                else
                                    LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("PROCESS TO SEND PARTICLE GOOD = " << Thread2Pos.ThreadPosX << "," << Thread2Pos.ThreadPosY << "," << Thread2Pos.ThreadPosZ << " FROM " << Thread1Pos.ThreadPosX << "," << Thread1Pos.ThreadPosY << "," << Thread1Pos.ThreadPosZ << " Current Process " << MPIProcessDataObject.CurrentMPIProcessIndex <<  " S1 = " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " S2 = " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << " P = " << ParticleObject.Center.X << " " << ParticleObject.Center.Y << " " << ParticleObject.Center.Z << " V = " << VectorX << " " << VectorY << " " << VectorZ << " PSHIFT = " << ParticleObject.Center.X + VectorX << " " << ParticleObject.Center.Y + VectorY << " " << ParticleObject.Center.Z + VectorZ));

                                // const auto Thread1Index = (*ThreadsIndexes)[Thread1Pos.ThreadPosX][Thread1Pos.ThreadPosY][Thread1Pos.ThreadPosZ];
                                // const auto Thread2Index = (*ThreadsIndexes)[Thread2Pos.ThreadPosX][Thread2Pos.ThreadPosY][Thread2Pos.ThreadPosZ];
                                const auto Thread1Index = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread1Pos.ThreadPosX - 1][Thread1Pos.ThreadPosY - 1][Thread1Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;
                                const auto Thread2Index = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;;
                                VectorOfParticlesToSendToNeighborThreads[NeighborProcessIndex].emplace_back(ParticleSenderStruct{ ParticleObject.Index, ParticleObject.EntityId, static_cast<int>(Thread1Index), static_cast<int>(Thread2Index), { static_cast<uint16_t>(SectorPosX2), static_cast<uint16_t>(SectorPosY2), static_cast<uint16_t>(SectorPosZ2) }, { ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z } });

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
