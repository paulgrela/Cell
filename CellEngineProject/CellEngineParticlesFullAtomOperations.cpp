
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
            const ThreadPosType ThreadPos = { ParticleSectorXIndex / CellEngineConfigDataObject.NumberOfXSectorsInOneThreadInSimulation + 1, ParticleSectorYIndex / CellEngineConfigDataObject.NumberOfYSectorsInOneThreadInSimulation + 1, ParticleSectorZIndex / CellEngineConfigDataObject.NumberOfZSectorsInOneThreadInSimulation + 1 };
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].ThreadPos = ThreadPos;
            ParticlesSectors[ParticleSectorXIndex][ParticleSectorYIndex][ParticleSectorZIndex].MPIProcessIndex = CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ThreadPos.ThreadPosX - 1][ThreadPos.ThreadPosY - 1][ThreadPos.ThreadPosZ - 1]->GetMPIProcessIndex() - 1;
        }
    }
    CATCH("setting proper thread index for every particles sector")
}

bool ExchangeParticleBetweenSectors(const Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, UnsignedInt const SectorPosX1, UnsignedInt const SectorPosY1, const UnsignedInt SectorPosZ1, const UnsignedInt SectorPosX2, const UnsignedInt SectorPosY2, const UnsignedInt SectorPosZ2)
{
    if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
        ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
    else
    {
        if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObject.Index) == false)
            if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObject.Index)->second == &ParticleObject)
            {
                if constexpr(PrintInfoWarning == true)
                    std::cout << "Earlier correct sector for particle index: " << ParticleObject.EntityId << " " << ParticleObject.Index << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                return true;
            }
    }

    return false;
}

void CellEngineParticlesFullAtomOperations::MoveParticleByVector(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SignedInt* NeighbourProcessesIndexes, std::vector<MPIParticleSenderStruct>* VectorOfParticlesToSendToNeighbourProcesses, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
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
            if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false)
            {
                const auto Thread1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].ThreadPos;
                const auto Thread2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].ThreadPos;

                if (CellEngineConfigDataObject.MultiThreaded == true && Thread1Pos != Thread2Pos)
                {
                    std::scoped_lock LockGuardScopedLock{ CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread1Pos.ThreadPosX - 1][Thread1Pos.ThreadPosY - 1][Thread1Pos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject, CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[Thread2Pos.ThreadPosX - 1][Thread2Pos.ThreadPosY - 1][Thread2Pos.ThreadPosZ - 1]->MainExchangeParticlesMutexObject };

                    if (ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2) == true)
                        return;
                }
            }
            else
            if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true && NeighbourProcessesIndexes != nullptr && VectorOfParticlesToSendToNeighbourProcesses != nullptr)
            {
                bool NewSectorNeighbourProcessFound = false;
                const auto Process1Pos = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].MPIProcessIndex;
                const auto Process2Pos = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].MPIProcessIndex;
                if (Process1Pos != Process2Pos)
                {
                    for (UnsignedInt NeigbourhProcessIndex = 0; NeigbourhProcessIndex < NumberOfAllNeighbours; NeigbourhProcessIndex++)
                        if (NeighbourProcessesIndexes[NeigbourhProcessIndex] == Process2Pos)
                        {
                            LoggersManagerObject.Log(STREAM("PROCESS TO SEND PARTICLE = " << Process2Pos));

                            //VectorOfParticlesToSendToNeighbourProcesses[NeigbourhProcessIndex].emplace_back(MPIParticleSenderStruct{ ParticleObject.Index, ParticleObject.EntityId, static_cast<int>(Process1Pos), static_cast<int>(Process2Pos), { static_cast<uint16_t>(SectorPosX2), static_cast<uint16_t>(SectorPosY2), static_cast<uint16_t>(SectorPosZ2) }, { ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z } });
                            VectorOfParticlesToSendToNeighbourProcesses[NeigbourhProcessIndex].emplace_back(MPIParticleSenderStruct{ ParticleObject.Index, ParticleObject.EntityId, static_cast<int>(MPIProcessDataObject.CurrentMPIProcessIndex), static_cast<int>(Process2Pos), { static_cast<uint16_t>(SectorPosX2), static_cast<uint16_t>(SectorPosY2), static_cast<uint16_t>(SectorPosZ2) }, { ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z } });
                                                                                                                        //ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].FreeIndexesOfParticles.push(ParticleObject.Index);
                                                                                                                        // const UniqueIdInt ParticleIndex = ParticleObject.Index;
                                                                                                                        // ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.erase(ParticleIndex);

                            NewSectorNeighbourProcessFound = true;
                            break;
                        }
                }
                else
                if (Process1Pos == Process2Pos)
                {
                    ExchangeParticleBetweenSectors(ParticleObject, ParticlesInSector, SectorPosX1, SectorPosY1, SectorPosZ1, SectorPosX2, SectorPosY2, SectorPosZ2);
                    NewSectorNeighbourProcessFound = true;
                }

                if (NewSectorNeighbourProcessFound == false)
                {
                    MoveAllAtomsInParticleAtomsListByVector(ParticleObject, -VectorX, -VectorY, -VectorZ);
                    ParticleObject.SetCenterCoordinates(ParticleObject.Center.X - VectorX, ParticleObject.Center.Y - VectorY, ParticleObject.Center.Z - VectorZ);
                }
            }
        }
    }
    CATCH_AND_THROW("moving particle by vector")
}
