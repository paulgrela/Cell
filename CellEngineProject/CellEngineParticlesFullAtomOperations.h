
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

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
    static void MoveParticleByVectorForThreads(Particle& ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, RealType VectorX, RealType VectorY, RealType VectorZ, const ThreadPosType& CurrentThreadPos);
    static void MoveParticleByVectorForMPIProcesses(Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighborProcessesIndexes, std::vector<ParticleSenderStructMPIMultiProcess>* VectorOfParticlesToSendToNeighborProcesses, RealType VectorX, RealType VectorY, RealType VectorZ);
protected:
    static bool CheckSectorPos(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ)
    {
        return (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
    }

    static bool CheckDistanceOfParticlesInSector(const RealType Radius, const UniqueIdUInt64Type Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles | std::views::values)
            if (ParticleInSectorObject.ParticleUniqueId != Index)
                if (DistanceOfParticleFromPoint(ParticleInSectorObject, { .X = TestedPosX, .Y = TestedPosY, .Z = TestedPosZ }) < ParticleInSectorObject.Radius + Radius)
                    return false;

        return true;
    }

    static bool CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(const UniqueIdUInt64Type Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles | std::views::values)
            if (ParticleInSectorObject.ParticleUniqueId != Index)
                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.ListOfAtoms)
                    if (DistanceOfPoints({ .X = AtomParticleInSectorObject.X, .Y = AtomParticleInSectorObject.Y, .Z = AtomParticleInSectorObject.Z }, { .X = TestedPosX, .Y = TestedPosY, .Z = TestedPosZ }) < 2 * AtomRadius)
                        return false;

        return true;
    }

    static PosType GetNewPosMovedByVector(const RealType X, const RealType Y, const RealType Z, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        return { .PosX = X + VectorX, .PosY = Y + VectorY, .PosZ = Z + VectorZ };
    }

    static inline bool CheckBoundsForSpaceV2(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
    {
        if (!(TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX && TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY && TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ))
            return false;

        return true;
    }

protected:
    static inline bool CheckBoundsForSpaceV1Reactions(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
    {
        DEBUG_LOG_CHEMICAL_REACTION(LoggersManagerObject.Log(STREAM("R2 = (" << TestedPosX << " " << TestedPosY << " " << TestedPosZ << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")"));)

        if (!((TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX) && (TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY) && (TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ)))
            return false;

        return true;
    }

    static inline bool CheckBoundsBySpaceV1V2Reactions(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySpaceBounds1, const bool CompareBoundsBySpaceBounds2)
    {
        if (CheckBounds == true)
        {
            if (CompareBoundsBySpaceBounds1 == true)
            {
                if (CheckBoundsForSpaceV1Reactions(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam) == false)
                    return false;
            }
            if (CompareBoundsBySpaceBounds2 == true)
            {
                if (CheckBoundsForSpaceV2(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam) == false)
                    return false;
            }
        }

        return true;
    }

    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorReactions(const RealType Radius, const UniqueIdUInt64Type Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySpaceBounds1, const bool CompareBoundsBySpaceBounds2)
    {
        try
        {
            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CheckBoundsBySpaceV1V2Reactions(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySpaceBounds1, CompareBoundsBySpaceBounds2) == false)
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
    static inline bool CheckBoundsForSpaceV1InDiffusion(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
    {
        DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("R2 = (" << TestedPosX << " " << TestedPosY << " " << TestedPosZ << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")"));)

        if (!(
            (TestedPosX >= -CellEngineConfigDataObject.ShiftCenterX && TestedPosX < CellEngineConfigDataObject.ShiftCenterX && TestedPosY >= -CellEngineConfigDataObject.ShiftCenterY && TestedPosY < CellEngineConfigDataObject.ShiftCenterY && TestedPosZ >= -CellEngineConfigDataObject.ShiftCenterZ && TestedPosZ < CellEngineConfigDataObject.ShiftCenterZ)
            &&
            (TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos - ParticleMoveSizeInDiffusionRealType && TestedPosX < SimulationSpaceSectorBoundsObjectParam.EndXPos + ParticleMoveSizeInDiffusionRealType && TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos - ParticleMoveSizeInDiffusionRealType && TestedPosY < SimulationSpaceSectorBoundsObjectParam.EndYPos + ParticleMoveSizeInDiffusionRealType && TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos - ParticleMoveSizeInDiffusionRealType && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.EndZPos + ParticleMoveSizeInDiffusionRealType)
            ))
        {
            DEBUG_LOG_DIFFUSION(LoggersManagerObject.Log(STREAM("R3 = (" << TestedPosX << " " << TestedPosY << " " << TestedPosZ << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")"));)

            return false;
        }

        return true;
    }

protected:
    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorDiffusion(const RealType Radius, const UniqueIdUInt64Type Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CompareBoundsBySpaceBounds1)
    {
        try
        {
            #ifdef SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION
            auto [TestedPosX1, TestedPosY1, TestedPosZ1] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, 0, 0, 0);
            LoggersManagerObject.Log(STREAM("R2S = (" << TestedPosX1 << " " << TestedPosY1 << " " << TestedPosZ1 << ") (" << SimulationSpaceSectorBoundsObjectParam.StartXPos << "," << SimulationSpaceSectorBoundsObjectParam.EndXPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartYPos << "," << SimulationSpaceSectorBoundsObjectParam.EndYPos << ") (" << SimulationSpaceSectorBoundsObjectParam.StartZPos << "," << SimulationSpaceSectorBoundsObjectParam.EndZPos << ")"));
            #endif

            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CompareBoundsBySpaceBounds1 == true)
            {
                if (CheckBoundsForSpaceV1InDiffusion(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam) == false)
                    return false;
            }

            if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                return false;

            // if (CheckDistanceOfParticlesInSector(Radius, Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
            //     return false;
        }
        CATCH("checking free space only by center and bounds for particle moved by vector")

        return true;
    }

protected:
    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorReactions(const ListOfAtomsType& ListOfAtoms, const UniqueIdUInt64Type Index, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySpaceBounds1, const bool CompareBoundsBySpaceBounds2)
    {
        try
        {
            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (CheckBoundsBySpaceV1V2Reactions(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySpaceBounds1, CompareBoundsBySpaceBounds2) == false)
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

    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorDiffusion(const ListOfAtomsType& ListOfAtoms, const UniqueIdUInt64Type Index, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CompareBoundsBySpaceBounds1)
    {
        try
        {
            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (CompareBoundsBySpaceBounds1 == true)
                {
                    if (CheckBoundsForSpaceV1InDiffusion(TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam) == false)
                        return false;
                }

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
    static inline bool CheckFreeSpaceAndBoundsForListOfAtomsReactions(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(ListOfAtoms, Radius, 0, { .X = 0, .Y = 0, .Z = 0 }, ParticlesInSector, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, true);
    }

    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdUInt64Type Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CheckBounds, const bool CompareBoundsBySpaceBounds1, const bool CompareBoundsBySpaceBounds2)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorReactions(Radius, Index, Center, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySpaceBounds1, CompareBoundsBySpaceBounds2);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorReactions(ListOfAtoms, Index, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySpaceBounds1, CompareBoundsBySpaceBounds2);
    }

    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdUInt64Type Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CompareBoundsBySpaceBounds1)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVectorDiffusion(Radius, Index, Center, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CompareBoundsBySpaceBounds1);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVectorDiffusion(ListOfAtoms, Index, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CompareBoundsBySpaceBounds1);
    }

protected:
    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForMPIProcessesDiffusion(Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const SignedInt* NeighbourProcessesIndexes, std::vector<ParticleSenderStructMPIMultiProcess>* VectorOfParticlesToSendToNeighbourProcesses, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.ParticleUniqueId, ParticleObject.Center, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ .StartXPos = StartXPosParam, .StartYPos = StartYPosParam, .StartZPos = StartZPosParam, .SizeX = SizeXParam, .SizeY = SizeYParam, .SizeZ = SizeZParam, .EndXPos = SizeXParam, .EndYPos = SizeYParam, .EndZPos = SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true) == true)
                MoveParticleByVectorForMPIProcesses(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, NeighbourProcessesIndexes, VectorOfParticlesToSendToNeighbourProcesses, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }

    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBoundsForThreadsDiffusion(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const ThreadPosType& CurrentThreadPos)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.ParticleUniqueId, ParticleObject.Center, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ .StartXPos = StartXPosParam, .StartYPos = StartYPosParam, .StartZPos = StartZPosParam, .SizeX = SizeXParam, .SizeY = SizeYParam, .SizeZ = SizeZParam, .EndXPos = SizeXParam, .EndYPos = SizeYParam, .EndZPos = SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true) == true)
                MoveParticleByVectorForThreads(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, VectorX, VectorY, VectorZ, CurrentThreadPos);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }

protected:
    static inline bool CheckFreeSpaceAndBoundsForListOfAtomsDiffusion(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVectorDiffusion(ListOfAtoms, Radius, 0, { .X = 0, .Y = 0, .Z = 0 }, ParticlesInSector, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true);
    }

protected:
    static inline bool CheckFreeSpaceForParticleMovedByVectorReactionsOfTranscriptionTranslation(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVectorReactions(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.ParticleUniqueId, ParticleObject.Center, ParticlesInSector, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ .StartXPos = 0, .StartYPos = 0, .StartZPos = 0, .SizeX = 0, .SizeY = 0, .SizeZ = 0, .EndXPos = 0, .EndYPos = 0, .EndZPos = 0}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }

protected:
    static inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpaceReactionsOfTranscriptionTranslation(Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, std::vector<ParticleToBeMovedFromOneSectorToAnotherSector>& ListOfParticlesToChangeSectors, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ)
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
                        if (CheckFreeSpaceForParticleMovedByVectorReactionsOfTranscriptionTranslation(ParticleObject, ParticlesInSector, PosX, PosY, PosZ) == true)
                        {
                            DEBUG_LOG(LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "FREE SPACE FOUND " << VecX << " " << VecY << " " << VecZ << " " << PosX << " " << PosY << " " << PosZ << terminal_colors_utils::white));)

                            MoveParticleByVectorForMPIProcesses(ParticleObject, ParticlesInSector, ListOfParticlesToChangeSectors, nullptr, nullptr, PosX, PosY, PosZ);
                            FoundFreeSpace = true;
                            goto Outside;
                        }
            Outside:

            #ifdef SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION
            if (FoundFreeSpace == false)
                LoggersManagerObject.Log(STREAM(terminal_colors_utils::red << "FREE SPACE NOT FOUND " << VecX << " " << VecY << " " << VecZ << " " << terminal_colors_utils::white));
            #endif
        }
        CATCH("moving particle near other particles")
    }
};

#endif
