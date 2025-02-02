
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
            for (auto& AtomObject: ParticleObject.ListOfAtoms)
            {
                AtomObject.X += VectorX;
                AtomObject.Y += VectorY;
                AtomObject.Z += VectorZ;
            }
        }
        CATCH_AND_THROW("moving all voxels in particle voxel list by vector")
    }
protected:
    static inline void MoveParticleByVector(Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ)
    {
        try
        {
            //SetAllAtomsInListOfAtomsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());

            MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);

            //SetAllAtomsInListOfAtomsToValue(ParticleObject.ListOfVoxels, ParticleObject.Index);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            const UnsignedInt SectorPosX1 = std::floor((ParticleObject.Center.X + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
            const UnsignedInt SectorPosY1 = std::floor((ParticleObject.Center.X + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
            const UnsignedInt SectorPosZ1 = std::floor((ParticleObject.Center.X + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

            const UnsignedInt SectorPosX2 = std::floor((ParticleObject.Center.X + VectorX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
            const UnsignedInt SectorPosY2 = std::floor((ParticleObject.Center.X + VectorX + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
            const UnsignedInt SectorPosZ2 = std::floor((ParticleObject.Center.X + VectorX + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

            if (SectorPosX1 != SectorPosX2 || SectorPosY1 != SectorPosY2 || SectorPosZ1 != SectorPosZ2)
                ParticlesInSector[SectorPosX2][SectorPosY2][SectorPosZ2].insert(ParticlesInSector[SectorPosX1][SectorPosY1][SectorPosZ1].extract(ParticleObject.Index));
        }
        CATCH_AND_THROW("moving particle by vector")
    }
protected:

    static inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const bool CheckOnlyParticlesCenters)
    {
        try
        {
            if (CheckOnlyParticlesCenters == true)
            {
                const float TestedPosX = ParticleObject.Center.X + VectorX;
                const float TestedPosY = ParticleObject.Center.Y + VectorY;
                const float TestedPosZ = ParticleObject.Center.Z + VectorZ;

                const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                    for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
                        if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                            if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + ParticleObject.Radius)
                                return false;
            }
            else
                for (auto &AtomParticleObject : ParticleObject.ListOfAtoms)
                {
                    const float TestedPosX = AtomParticleObject.X + VectorX;
                    const float TestedPosY = AtomParticleObject.Y + VectorY;
                    const float TestedPosZ = AtomParticleObject.Z + VectorZ;

                    const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                    const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                    const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                    if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
                            if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                                    if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                        return false;
                }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }
protected:
    static inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(const Particle &ParticleObject, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam, const bool CheckOnlyParticlesCenters)
    {
        try
        {
            if (CheckOnlyParticlesCenters == true)
            {
                const float TestedPosX = ParticleObject.Center.X + VectorX;
                const float TestedPosY = ParticleObject.Center.Y + VectorY;
                const float TestedPosZ = ParticleObject.Center.Z + VectorZ;

                const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;

                if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                    for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
                        if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                            if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + ParticleObject.Radius)
                                return false;
            }
            else
                for (auto &AtomParticleObject : ParticleObject.ListOfAtoms)
                {
                    const float TestedPosX = AtomParticleObject.X + VectorX;
                    const float TestedPosY = AtomParticleObject.Y + VectorY;
                    const float TestedPosZ = AtomParticleObject.Z + VectorZ;

                    const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                    const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                    const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                    if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                        return false;

                    if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                        for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
                            if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                                for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                                    if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                        return false;
                }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }
protected:
    static inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, ParticlesContainer<Particle>& ParticlesInSector, const CurrentSectorPosType& CurrentSectorPos, const float Radius, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyParticlesCenters)
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

                const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;

                if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                    for (const auto& ParticleInSectorObject : ParticlesInSector[SectorPosX][SectorPosY][SectorPosZ])
                        if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                            return false;
            }
            else
                for (auto &AtomParticleObject : ListOfAtoms)
                {
                    const float TestedPosX = AtomParticleObject.X + VectorX;
                    const float TestedPosY = AtomParticleObject.Y + VectorY;
                    const float TestedPosZ = AtomParticleObject.Z + VectorZ;

                    const UnsignedInt SectorPosX = std::floor((TestedPosX + CellEngineConfigDataObject.ShiftCenterX) / CellEngineConfigDataObject.SizeOfParticlesSectorX);
                    const UnsignedInt SectorPosY = std::floor((TestedPosY + CellEngineConfigDataObject.ShiftCenterY) / CellEngineConfigDataObject.SizeOfParticlesSectorY);
                    const UnsignedInt SectorPosZ = std::floor((TestedPosZ + CellEngineConfigDataObject.ShiftCenterZ) / CellEngineConfigDataObject.SizeOfParticlesSectorZ);

                    if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                        return false;

                    for (const auto& ParticleInSectorObject : ParticlesInSector[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ])
                        if (SectorPosX < CellEngineConfigDataObject.NumberOfParticlesSectorsInX && SectorPosY < CellEngineConfigDataObject.NumberOfParticlesSectorsInY && SectorPosZ < CellEngineConfigDataObject.NumberOfParticlesSectorsInZ)
                            for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                                if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                    return false;
                }
        }
        CATCH("checking free space for list of voxels")

        return true;
    }

protected:
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
