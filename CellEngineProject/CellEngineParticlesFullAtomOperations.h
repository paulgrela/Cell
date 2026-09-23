
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

#include <cmath>
#include <memory>
#include <ranges>

#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"
#include "CellEngineBasicParallelExecutionData.h"

class CellEngineParticlesFullAtomOperations
{
public:
    static void SetProperThreadIndexForEveryParticlesSector(ParticlesContainer<Particle>& ParticlesSectors);
protected:
    static void MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighborProcesses, RealType VectorX, RealType VectorY, RealType VectorZ, const ThreadPosType& CurrentThreadPos);
    static void MoveParticleByVectorForMPIProcesses(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighbourProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighbourProcesses, RealType VectorX, RealType VectorY, RealType VectorZ, ThreadPosType CurrentThreadPos);
protected:
    static inline void MoveAllAtomsInParticleAtomsListByVector(Particle &ParticleObject, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        try
        {
            for (auto& AtomObject : ParticleObject.ListOfAtoms)
            {
                AtomObject.X += VectorX;
                AtomObject.Y += VectorY;
                AtomObject.Z += VectorZ;
            }
        }
        CATCH_AND_THROW("moving all atoms from atoms list by vector")
    }
protected:
    static bool CheckSectorPos(const UnsignedInt SectorPosX, UnsignedInt SectorPosY, UnsignedInt SectorPosZ)
    {
        return (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
    }
protected:
    static bool CheckDistanceOfParticlesInSector(const RealType Radius, const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles)
            if (ParticleInSectorObject.second.Index != Index)
                if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { .X = TestedPosX, .Y = TestedPosY, .Z = TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                    return false;

        return true;
    }
protected:
    static bool CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles | std::ranges::views::values)
            if (ParticleInSectorObject.Index != Index)
                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.ListOfAtoms)
                    if (DistanceOfPoints({ .X = AtomParticleInSectorObject.X, .Y = AtomParticleInSectorObject.Y, .Z = AtomParticleInSectorObject.Z }, { .X = TestedPosX, .Y = TestedPosY, .Z = TestedPosZ }) < 2 * AtomRadius)
                        return false;

        return true;
    }
protected:
    static PosType GetNewPosMovedByVector(const RealType X, const RealType Y, const RealType Z, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        return { .PosX = X + VectorX, .PosY = Y + VectorY, .PosZ = Z + VectorZ };
    }

protected:
    static inline bool CheckBoundsForSpace(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (!(TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX && TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY && TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ))
            return false;

        return true;
    }



protected:
    static inline bool CheckBoundsForSectorReactions(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        // std::cout << "R2 = (" << TestedPosX << " " << TestedPosY << " " << TestedPosZ << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")" << std::endl;

        if (!((TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX) && (TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY) && (TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ)))
            return false;

        return true;
    }
protected:
    static inline bool CheckBoundsBySectorAndBySpaceReactions(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckBounds == true)
        {
            if (CompareBoundsBySectorsBounds == true)
            {
                if (CheckBoundsForSectorReactions(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
            if (CompareBoundsBySpaceBounds == true)
            {
                if (CheckBoundsForSpace(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
        }

        return true;
    }
protected:
    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorReactions(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CheckBoundsBySectorAndBySpaceReactions(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                return false;

            if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                return false;

            if (CheckDistanceOfParticlesInSector(Radius, Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                return false;
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }



protected:
    static inline bool CheckBoundsForSectorDiffusion(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        //std::cout << "R2 = (" << SectorPosX << " " << SectorPosY << " " << SectorPosZ << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")" << std::endl;

        if (!(SectorPosX >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos) && SectorPosX < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX) && SectorPosY >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos) && SectorPosY < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY) && SectorPosZ >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos) && SectorPosZ < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ)))
            return false;

        return true;
    }
protected:
    static inline bool CheckBoundsBySectorAndBySpaceDiffusion(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckBounds == true)
        {
            if (CompareBoundsBySectorsBounds == true)
            {
                if (CheckBoundsForSectorDiffusion(SectorPosX, SectorPosY, SectorPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
            if (CompareBoundsBySpaceBounds == true)
            {
                if (CheckBoundsForSpace(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;
            }
        }

        return true;
    }
protected:
    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorDiffusion(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CheckBoundsBySectorAndBySpaceDiffusion(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                return false;

            if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                return false;

            if (CheckDistanceOfParticlesInSector(Radius, Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                return false;
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }



protected:
    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorReactions(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            // auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            // auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (CheckBoundsBySectorAndBySpaceReactions(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }
protected:
    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorDiffusion(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            // auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            // auto [SectorPosXC, SectorPosYC, SectorPosCZ] = CellEngineUseful::GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (CheckBoundsBySectorAndBySpaceDiffusion(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
                    return false;

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }



protected:
    //CCC - w sprawdzeniu czy nowa czastka w reakcji
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorReactions(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorReactions(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
    }
protected:
    static inline bool CheckFreeSpaceAndBoundsForListOfAtomsReactions(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
    {
         return CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(ListOfAtoms, Radius, 0, { .X = 0, .Y = 0, .Z = 0 }, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, true);
    }





protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorDiffusion(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorDiffusion(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
    }
protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForMPIProcesses(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighbourProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighbourProcesses, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const ThreadPosType CurrentThreadPos)
    {
        try
        {
            //#ifdef CONTAINERS_FOR_SPEED
            if (CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
            //#else
            //if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObjectIter->second.ListOfAtoms, ParticleObjectIter->second.Radius, ParticleObjectIter->second.Index, ParticleObjectIter->second.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
            //#endif
                MoveParticleByVectorForMPIProcesses(ParticleObject, ParticlesInSector, ParticleObjectIter, ListOfParticlesToChangeSectors, NeighbourProcessesIndexes, VectorOfParticlesToSendToNeighbourProcesses, VectorX, VectorY, VectorZ, CurrentThreadPos);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }

    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForThreads(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, ParticlesDetailedContainer<Particle>::iterator& ParticleObjectIter, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighbourProcessesIndexes, std::vector<ParticleSenderStruct>* VectorOfParticlesToSendToNeighbourProcesses, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const ThreadPosType& CurrentThreadPos)
    {
        try
        {
            //#ifdef CONTAINERS_FOR_SPEED
            if (CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
            //#else
            //if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObjectIter->second.ListOfAtoms, ParticleObjectIter->second.Radius, ParticleObjectIter->second.Index, ParticleObjectIter->second.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
            //#endif
                MoveParticleByVectorForThreads(ParticleObject, ParticlesInSector, ParticleObjectIter, ListOfParticlesToChangeSectors, NeighbourProcessesIndexes, VectorOfParticlesToSendToNeighbourProcesses, VectorX, VectorY, VectorZ, CurrentThreadPos);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    static inline bool CheckFreeSpaceAndBoundsForListOfAtomsDiffusion(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ListOfAtoms, Radius, 0, { .X = 0, .Y = 0, .Z = 0 }, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, false, true);
    }




protected:
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const bool CheckOnlyParticlesCenters)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ 0, 0, 0, 0, 0, 0, 0, 0, 0}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
protected:
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SectorPosType& CurrentSectorPos, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ, const ThreadPosType CurrentThreadPos)
    {
        try
        {
            bool FoundFreeSpace = false;

            const RealType VecX = NewPositionParticleObject.ListOfAtoms[0].X - ParticleObject.ListOfAtoms[0].X;
            const RealType VecY = NewPositionParticleObject.ListOfAtoms[0].Y - ParticleObject.ListOfAtoms[0].Y;
            const RealType VecZ = NewPositionParticleObject.ListOfAtoms[0].Z - ParticleObject.ListOfAtoms[0].Z;

            for (RealType PosX = VecX - AddX; PosX < VecX + AddX; PosX += 1.0)
                for (RealType PosY = VecY - AddY; PosY < VecY + AddY; PosY += 1.0)
                    for (RealType PosZ = VecZ - AddZ; PosZ < VecZ + AddZ; PosZ += 1.0)
                        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, CellEngineConfigDataObject.CheckOnlyParticlesCenters) == true)
                        {
                            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "FREE SPACE FOUND " << VecX << " " << VecY << " " << VecZ << " " << PosX << " " << PosY << " " << PosZ << terminal_colors_utils::white));
                            ParticlesDetailedContainer<Particle>::iterator ParticleIter;
                            MoveParticleByVectorForMPIProcesses(ParticleObject, ParticlesInSector, ParticleIter, ListOfParticlesToChangeSectors, nullptr, nullptr, PosX, PosY, PosZ, CurrentThreadPos);
                            FoundFreeSpace = true;
                            goto Outside;
                        }
            Outside:

            if (FoundFreeSpace == false)
                LoggersManagerObject.Log(STREAM(terminal_colors_utils::red << "FREE SPACE NOT FOUND " << VecX << " " << VecY << " " << VecZ << " " << terminal_colors_utils::white));
        }
        CATCH("moving particle near other particles")
    }
};

#endif
