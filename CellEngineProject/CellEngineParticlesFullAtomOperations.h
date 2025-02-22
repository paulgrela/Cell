
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

#include <cmath>

#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"

class CellEngineParticlesFullAtomOperations
{
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
        CATCH_AND_THROW("moving all voxels in particle voxel list by vector")
    }
protected:
    static inline void MoveParticleByVector(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        try
        {
            auto [SectorPosX1, SectorPosY1, SectorPosZ1] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
            auto [SectorPosX2, SectorPosY2, SectorPosZ2] = CellEngineUseful::GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
            {
                if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
                    ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
                else
                {
                    if (ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.contains(ParticleObject.Index) == false)
                        if (&ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.find(ParticleObject.Index)->second == &ParticleObject)
                        {
                            std::cout << "Bad sector for particle index: " << ParticleObject.Index << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                            return;
                        }

                    const UniqueIdInt FormerParticleIndex = ParticleObject.Index;

                    ParticleObject.Index = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].FreeIndexesOfParticles.top();
                    ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].FreeIndexesOfParticles.pop();

                    ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles[ParticleObject.Index] = ParticleObject;

                    ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.erase(ParticleObject.Index);

                    std::cout << "New particle index: " << ParticleObject.Index << " Former particle index: " << FormerParticleIndex << " " << SectorPosX1 << " " << SectorPosY1 << " " << SectorPosZ1 << " " << SectorPosX2 << " " << SectorPosY2 << " " << SectorPosZ2 << std::endl;
                }

                // if (ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.contains(ParticleObject.Index) == false)
                //     ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(ParticleObject.Index));
                // else
                // {
                //     const UniqueIdInt FormerParticleIndex = ParticleObject.Index;
                //
                //     //UniqueIdInt NewParticleIndex = ParticleObject.Index = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].FreeIndexesOfParticles.top();
                //     const UniqueIdInt NewParticleIndex = ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].FreeIndexesOfParticles.top();
                //     ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].FreeIndexesOfParticles.pop();
                //
                //     //wyglada na to ze wczesniej czastka przypisana do zlego sektora ale czemu i niema jej w sektorze1
                //     auto ParticleExtract = ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].Particles.extract(FormerParticleIndex);
                //     if (ParticleExtract.empty() == false)
                //     {
                //         ParticleExtract.mapped().Index = NewParticleIndex;
                //         ParticleExtract.key() = NewParticleIndex;
                //         std::cout << "Correct New particle index: " << NewParticleIndex << " Former particle index: " << FormerParticleIndex << std::endl;
                //     }
                //     else
                //         std::cout << "New particle index: " << NewParticleIndex << " Former particle index: " << FormerParticleIndex << std::endl;
                //
                //     ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].Particles.insert(std::move(ParticleExtract));
                // }
            }
        }
        CATCH_AND_THROW("moving particle by vector")
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
                if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                    return false;

        return true;
    }
protected:
    static bool CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ].Particles)
            if (ParticleInSectorObject.second.Index != Index)
                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                    if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                        return false;

        return true;
    }
protected:
    static PosType GetNewPosMovedByVector(const RealType X, const RealType Y, const RealType Z, const RealType VectorX, const RealType VectorY, const RealType VectorZ)
    {
        return { X + VectorX, Y + VectorY, Z + VectorZ };
    }
protected:
    static inline bool CheckBoundsForSector(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (!(SectorPosX >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos) && SectorPosX < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX) && SectorPosY >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos) && SectorPosY < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY) && SectorPosZ >= static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos) && SectorPosZ < static_cast<UnsignedInt>(SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ)))
            return false;

        return true;
    }
    static inline bool CheckBoundsForSpace(const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (!(TestedPosX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && TestedPosX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX && TestedPosY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && TestedPosY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY && TestedPosZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && TestedPosZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ))
            return false;

        return true;
    }
protected:
    static inline bool CheckBoundsBySectorAndBySpace(const UnsignedInt SectorPosX, const UnsignedInt SectorPosY, const UnsignedInt SectorPosZ, const RealType TestedPosX, const RealType TestedPosY, const RealType TestedPosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckBounds == true)
        {
            if (CompareBoundsBySectorsBounds == true)
            {
                if (CheckBoundsForSector(SectorPosX, SectorPosY, SectorPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
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
    static inline bool CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

            if (CheckBoundsBySectorAndBySpace(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
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
    static inline bool CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        try
        {
            auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(Center.X, Center.Y, Center.Z, VectorX, VectorY, VectorZ);
            auto [SectorPosX, SectorPosY, SectorPosZ] = CellEngineUseful::GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

            for (auto &AtomParticleObject : ListOfAtoms)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);

                if (CheckBoundsBySectorAndBySpace(SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds) == false)
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
    //BBB - uzyte w normalnej dyfuzji
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(const ListOfAtomsType& ListOfAtoms, const RealType Radius, const UniqueIdInt Index, const vector3_Real32 Center, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters, const bool CheckBounds, const bool CompareBoundsBySectorsBounds, const bool CompareBoundsBySpaceBounds)
    {
        if (CheckOnlyParticlesCenters == true)
            return CheckFreeSpaceOnlyByTestingCenterAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
        else
            return CheckFreeSpaceByTestingAllAtomsAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, Index, Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBoundsObjectParam, CheckBounds, CompareBoundsBySectorsBounds, CompareBoundsBySpaceBounds);
    }
protected:
    //CCC - w sprawdzeniu czy nowa czastka w reakcji
    static inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType Radius, const RealType PosX, const RealType PosY, const RealType PosZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
    {
         return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, { 0, 0, 0 }, ParticlesInSector, CurrentSectorPos, PosX, PosY, PosZ, SimulationSpaceSectorBoundsObjectParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, false, true);
    }
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const bool CheckOnlyParticlesCenters)
    {
        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ 0, 0, 0, 0, 0, 0, 0, 0, 0}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const RealType VectorX, const RealType VectorY, const RealType VectorZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject.ListOfAtoms, ParticleObject.Radius, ParticleObject.Index, ParticleObject.Center, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, SimulationSpaceSectorBounds{ StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, StartXPosParam + SizeXParam, StartYPosParam + SizeYParam, StartZPosParam + SizeZParam }, CellEngineConfigDataObject.CheckOnlyParticlesCenters, true, true, false) == true)
                MoveParticleByVector(ParticleObject, ParticlesInSector, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const SectorPosType& CurrentSectorPos, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ)
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
                            MoveParticleByVector(ParticleObject, ParticlesInSector, PosX, PosY, PosZ);
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
