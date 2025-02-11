
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

#include <cmath>

#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"

class CellEngineParticlesFullAtomOperations
{
protected:
    static inline void MoveAllAtomsInParticleAtomsListByVector(Particle &ParticleObject, const float VectorX, const float VectorY, const float VectorZ)
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
    static CurrentSectorPosType GetSectorPos(const float X, const float Y, const float Z)
    {
        const UnsignedInt SectorPosX = std::floor((X + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
        const UnsignedInt SectorPosY = std::floor((Y + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
        const UnsignedInt SectorPosZ = std::floor((Z + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

        return { SectorPosX, SectorPosY, SectorPosZ };
    }

    static inline void MoveParticleByVector(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ)
    {
        try
        {
            MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            auto [SectorPosX1, SectorPosY1, SectorPosZ1] = GetSectorPos(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z);
            auto [SectorPosX2, SectorPosY2, SectorPosZ2] = GetSectorPos(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
                ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].extract(ParticleObject.Index));
        }
        CATCH_AND_THROW("moving particle by vector")
    }
protected:
    static bool CheckSectorPos(const UnsignedInt SectorPosX, UnsignedInt SectorPosY, UnsignedInt SectorPosZ)
    {
        return (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
    }
protected:
    static bool CheckDistanceOfParticlesInSector(const float Radius, const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, UnsignedInt SectorPosY, UnsignedInt SectorPosZ, const float TestedPosX, const float TestedPosY, const float TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
            if (ParticleInSectorObject.second.Index != Index)
                if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                    return false;

        return true;
    }
protected:
    static bool CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(const UniqueIdInt Index, const ParticlesContainer<Particle>& ParticlesInSector, const UnsignedInt SectorPosX, UnsignedInt SectorPosY, UnsignedInt SectorPosZ, const float TestedPosX, const float TestedPosY, const float TestedPosZ)
    {
        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
            if (ParticleInSectorObject.second.Index != Index)
                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                    if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                        return false;

        return true;
    }
protected:
    static PosType GetNewPosMovedByVector(const float X, const float Y, const float Z, const float VectorX, const float VectorY, const float VectorZ)
    {
        const float TestedPosX = X + VectorX;
        const float TestedPosY = Y + VectorY;
        const float TestedPosZ = Z + VectorZ;

        return { TestedPosX, TestedPosY, TestedPosZ };
    }

    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const bool CheckOnlyParticlesCenters)
    {
        try
        {
            if (CheckOnlyParticlesCenters == true)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z, VectorX, VectorY, VectorZ);

                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSector(ParticleObject.Radius, ParticleObject.Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
            else
            {
                auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

                for (auto &AtomParticleObject : ParticleObject.ListOfAtoms)
                {
                    auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);

                    if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                        return false;

                    if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(ParticleObject.Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                        return false;
                }
            }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }
protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam, const bool CheckOnlyParticlesCenters)
    {
        try
        {
            if (CheckOnlyParticlesCenters == true)
            {
                auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                //if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                if (!(SectorPosX >= StartXPosParam && SectorPosX < StartXPosParam + SizeXParam && SectorPosY >= StartYPosParam && SectorPosY < StartYPosParam + SizeYParam && SectorPosZ >= StartZPosParam && SectorPosZ < StartZPosParam + SizeZParam))
                    return false;

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSector(ParticleObject.Radius, ParticleObject.Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
            else
            {
                auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(ParticleObject.Center.X, ParticleObject.Center.Y, ParticleObject.Center.Z, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

                for (auto &AtomParticleObject : ParticleObject.ListOfAtoms)
                {
                    auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);

                    //if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    if (!(SectorPosX >= StartXPosParam && SectorPosX < StartXPosParam + SizeXParam && SectorPosY >= StartYPosParam && SectorPosY < StartYPosParam + SizeYParam && SectorPosZ >= StartZPosParam && SectorPosZ < StartZPosParam + SizeZParam))
                        return false;

                    if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                        return false;

                    if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(ParticleObject.Index, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                        return false;
                }
            }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }
protected:
    //CCC - w sprawdzeniu czy nowa czastka w reakcji
    static inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, const ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float Radius, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
    {
        try
        {
            const float StartXPosParam = SimulationSpaceSectorBoundsObjectParam.StartXPos;
            const float StartYPosParam = SimulationSpaceSectorBoundsObjectParam.StartYPos;
            const float StartZPosParam = SimulationSpaceSectorBoundsObjectParam.StartZPos;
            const float SizeXParam = SimulationSpaceSectorBoundsObjectParam.SizeX;
            const float SizeYParam = SimulationSpaceSectorBoundsObjectParam.SizeY;
            const float SizeZParam = SimulationSpaceSectorBoundsObjectParam.SizeZ;

            if (CheckOnlyParticlesCenters == true)
            {
                const float TestedPosX = VectorX;
                const float TestedPosY = VectorY;
                const float TestedPosZ = VectorZ;

                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosX, TestedPosY, TestedPosZ);

                if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;

                if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                    return false;

                if (CheckDistanceOfParticlesInSector(Radius, 0, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                    return false;
            }
            else
            {
                auto [TestedPosXC, TestedPosYC, TestedPosZC] = GetNewPosMovedByVector(VectorX, VectorY, VectorZ, VectorX, VectorY, VectorZ);
                auto [SectorPosX, SectorPosY, SectorPosZ] = GetSectorPos(TestedPosXC, TestedPosYC, TestedPosZC);

                for (auto &AtomParticleObject : ListOfAtoms)
                {
                    auto [TestedPosX, TestedPosY, TestedPosZ] = GetNewPosMovedByVector(AtomParticleObject.X, AtomParticleObject.Y, AtomParticleObject.Z, VectorX, VectorY, VectorZ);

                    if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                        return false;

                    if (CheckSectorPos(SectorPosX, SectorPosY, SectorPosZ) == false)
                        return false;

                    if (CheckDistanceOfParticlesInSectorByComparingAllAtomsDistances(0, ParticlesInSector, SectorPosX, SectorPosY, SectorPosZ, TestedPosX, TestedPosY, TestedPosZ) == false)
                        return false;
                }
            }
        }
        CATCH("checking free space for list of voxels")

        return true;
    }

protected:
    //BBB - uzyte w normalnej dyfuzji
    static inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject, ParticlesInSector, CurrentSectorPos, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, CellEngineConfigDataObject.CheckOnlyParticlesCenters) == true)
                MoveParticleByVector(ParticleObject, ParticlesInSector, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    //AAA - uzyte w reakcjach transkrypcji i translacji
    static inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const Particle &NewPositionParticleObject, const float AddX, const float AddY, const float AddZ)
    {
        try
        {
            bool FoundFreeSpace = false;

            const float VecX = NewPositionParticleObject.ListOfAtoms[0].X - ParticleObject.ListOfAtoms[0].X;
            const float VecY = NewPositionParticleObject.ListOfAtoms[0].Y - ParticleObject.ListOfAtoms[0].Y;
            const float VecZ = NewPositionParticleObject.ListOfAtoms[0].Z - ParticleObject.ListOfAtoms[0].Z;

            for (float PosX = VecX - AddX; PosX < VecX + AddX; PosX += 1.0)
                for (float PosY = VecY - AddY; PosY < VecY + AddY; PosY += 1.0)
                    for (float PosZ = VecZ - AddZ; PosZ < VecZ + AddZ; PosZ += 1.0)
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
