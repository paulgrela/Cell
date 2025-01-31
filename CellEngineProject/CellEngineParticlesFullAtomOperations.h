
#ifndef CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_FULL_ATOM_OPERATIONS_H

#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"

class CellEngineParticlesFullAtomOperations
{
public:

    using SimulationSpaceVoxel = UniqueIdInt;

    // void SetAllVoxelsInListOfVoxelsToValueForOuterClass(ListOfVoxelsType& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
    // {
    //     SetAllAtomsInListOfAtomsToValue(ListOfVoxels, SimulationSpaceVoxelValue);
    // }
public:
    static inline SimulationSpaceVoxel GetZeroSimulationSpaceVoxel()
    {
        return 0;
    }
protected:

    using Space_2048_2048_2048 = SimulationSpaceVoxel[2048][2048][2048];
    void* SpacePointer = nullptr;
    [[nodiscard]] inline SimulationSpaceVoxel& GetSpaceVoxel(const float x, const float y, const float z) const
    {
        return (*static_cast<Space_2048_2048_2048*>(SpacePointer))[1][1][1];
    }



    // inline void SetAllAtomsInListOfAtomsToValue(ListOfVoxelsType& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
    // {
    //     try
    //     {
    //         for (auto &Voxel: ListOfVoxels)
    //             GetSpaceVoxel(Voxel.X, Voxel.Y, Voxel.Z) = SimulationSpaceVoxelValue;
    //     }
    //     CATCH("making all zero voxels in list of voxels")
    // }
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
    static inline void MoveParticleByVector(Particle &ParticleObject, const float VectorX, const float VectorY, const float VectorZ)
    {
        try
        {
            //SetAllAtomsInListOfAtomsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());

            MoveAllAtomsInParticleAtomsListByVector(ParticleObject, VectorX, VectorY, VectorZ);

            //SetAllAtomsInListOfAtomsToValue(ParticleObject.ListOfVoxels, ParticleObject.Index);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);

            //PRZENIES DO INNEGO SEKTORA - I TU KLOPS BO NIE MAM DOSTEPU DO CALEJ
        }
        CATCH_AND_THROW("moving particle by vector")
    }
protected:

    inline bool CheckFreeSpaceForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ, const bool CheckOnlyCenter)
    {
        try
        {
            if (CheckOnlyCenter == true)
            {
                const float TestedPosX = ParticleObject.Center.X + VectorX;
                const float TestedPosY = ParticleObject.Center.Y + VectorY;
                const float TestedPosZ = ParticleObject.Center.Z + VectorZ;

                for (const auto& ParticleInSectorObject : ParticlesInSector)
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

                    for (const auto& ParticleInSectorObject : ParticlesInSector)
                        if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                            for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                                if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                    return false;
                }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }











// protected:
//     inline bool CheckFreeSpaceForListOfVoxels(const ListOfAtomsType& ListOfAtoms, const float VectorX, const float VectorY, const float VectorZ)
//     {
//         try
//         {
//             for (auto &VoxelOfParticle: ListOfAtoms)
//                 if (GetSpaceVoxel(VoxelOfParticle.X + VectorX, VoxelOfParticle.Y + VectorY, VoxelOfParticle.Z + VectorZ) != GetZeroSimulationSpaceVoxel())
//                     return false;
//         }
//         CATCH("checking free space for list of voxels")
//
//         return true;
//     }

// protected:
//     static inline bool CheckBoundsForParticleMovedByVector(Particle &ParticleObject, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
//     {
//         try
//         {
//             for (auto &VoxelOfParticle: ParticleObject.ListOfAtoms)
//                 if (!(VoxelOfParticle.X + VectorX >= StartXPosParam && VoxelOfParticle.X + VectorX < StartXPosParam + SizeXParam && VoxelOfParticle.Y + VectorY >= StartYPosParam && VoxelOfParticle.Y + VectorY < StartYPosParam + SizeYParam && VoxelOfParticle.Z + VectorZ >= StartZPosParam && VoxelOfParticle.Z + VectorZ < StartZPosParam + SizeZParam))
//                     return false;
//         }
//         CATCH("checking bounds for particle moved by vector")
//
//         return true;
//     }












    inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(const Particle &ParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam, const bool CheckOnlyCenter)
    {
        try
        {
            if (CheckOnlyCenter == true)
            {
                const float TestedPosX = ParticleObject.Center.X + VectorX;
                const float TestedPosY = ParticleObject.Center.Y + VectorY;
                const float TestedPosZ = ParticleObject.Center.Z + VectorZ;

                if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;

                for (const auto& ParticleInSectorObject : ParticlesInSector)
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

                    if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                        return false;

                    for (const auto& ParticleInSectorObject : ParticlesInSector)
                        if (ParticleInSectorObject.second.Index != ParticleObject.Index)
                            for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                            {
                                if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                    return false;
                            }
                }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }


// protected:
//     inline bool CheckFreeSpaceAndBoundsForParticleMovedByVector(Particle &ParticleObject, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
//     {
//         try
//         {
//             for (auto &VoxelOfParticle: ParticleObject.ListOfAtoms)
//             {
//                 float TestedPosX = VoxelOfParticle.X + VectorX;
//                 float TestedPosY = VoxelOfParticle.Y + VectorY;
//                 float TestedPosZ = VoxelOfParticle.Z + VectorZ;
//                 if ((GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != GetZeroSimulationSpaceVoxel() && GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != ParticleObject.Index) || !(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
//                     return false;
//             }
//         }
//         CATCH("checking bounds and free space for particle moved by vector")
//
//         return true;
//     }




protected:
    //inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyCenter)
    inline bool CheckFreeSpaceAndBoundsForListOfAtoms(const ListOfAtomsType& ListOfAtoms, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float Radius, const float VectorX, const float VectorY, const float VectorZ, const SimulationSpaceSectorBounds& SimulationSpaceSectorBoundsObjectParam, const bool CheckOnlyCenter)
    {
        try
        {
            // for (auto &VoxelOfParticle: ListOfAtoms)
            //     if (GetSpaceVoxel(VoxelOfParticle.X + VectorX, VoxelOfParticle.Y + VectorY, VoxelOfParticle.Z + VectorZ) != GetZeroSimulationSpaceVoxel() || !(VoxelOfParticle.X + VectorX >= SimulationSpaceSectorBoundsObjectParam.StartXPos && VoxelOfParticle.X + VectorX < SimulationSpaceSectorBoundsObjectParam.StartXPos + SimulationSpaceSectorBoundsObjectParam.SizeX && VoxelOfParticle.Y + VectorY >= SimulationSpaceSectorBoundsObjectParam.StartYPos && VoxelOfParticle.Y + VectorY < SimulationSpaceSectorBoundsObjectParam.StartYPos + SimulationSpaceSectorBoundsObjectParam.SizeY && VoxelOfParticle.Z + VectorZ >= SimulationSpaceSectorBoundsObjectParam.StartZPos && VoxelOfParticle.Z + VectorZ < SimulationSpaceSectorBoundsObjectParam.StartZPos + SimulationSpaceSectorBoundsObjectParam.SizeZ))
            //         return false;

            const float StartXPosParam = SimulationSpaceSectorBoundsObjectParam.StartXPos;
            const float StartYPosParam = SimulationSpaceSectorBoundsObjectParam.StartYPos;
            const float StartZPosParam = SimulationSpaceSectorBoundsObjectParam.StartZPos;
            const float SizeXParam = SimulationSpaceSectorBoundsObjectParam.SizeX;
            const float SizeYParam = SimulationSpaceSectorBoundsObjectParam.SizeY;
            const float SizeZParam = SimulationSpaceSectorBoundsObjectParam.SizeZ;

            if (CheckOnlyCenter == true)
            {
                const float TestedPosX = VectorX;
                const float TestedPosY = VectorY;
                const float TestedPosZ = VectorZ;

                if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;

                for (const auto& ParticleInSectorObject : ParticlesInSector)
                    if (DistanceOfParticleFromPoint(ParticleInSectorObject.second, { TestedPosX, TestedPosY, TestedPosZ }) < ParticleInSectorObject.second.Radius + Radius)
                        return false;
            }
            else
                for (auto &AtomParticleObject : ListOfAtoms)
                {
                    const float TestedPosX = AtomParticleObject.X + VectorX;
                    const float TestedPosY = AtomParticleObject.Y + VectorY;
                    const float TestedPosZ = AtomParticleObject.Z + VectorZ;

                    if (!(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                        return false;

                    for (const auto& ParticleInSectorObject : ParticlesInSector)
                        for (const auto &AtomParticleInSectorObject : ParticleInSectorObject.second.ListOfAtoms)
                        {
                            if (DistanceOfPoints({ AtomParticleInSectorObject.X, AtomParticleInSectorObject.Y, AtomParticleInSectorObject.Z }, { TestedPosX, TestedPosY, TestedPosZ }) < 2 * AtomRadius)
                                return false;
                        }
                }

        }
        CATCH("checking free space for list of voxels")

        return true;
    }






























protected:
    inline bool MoveParticleNearOtherParticleIfFullAtomSpaceIsEmpty(Particle &ParticleObject, const Particle &NewPositionParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float AddX, const float AddY, const float AddZ)
    {
        return MoveParticleByVectorIfFullAtomSpaceIsEmpty(ParticleObject, ParticlesInSector, NewPositionParticleObject.ListOfAtoms[0].X - ParticleObject.ListOfAtoms[0].X + AddX, NewPositionParticleObject.ListOfAtoms[0].Y - ParticleObject.ListOfAtoms[0].Y + AddY, NewPositionParticleObject.ListOfAtoms[0].Z - ParticleObject.ListOfAtoms[0].Z + AddZ);
    }
protected:
    inline bool MoveParticleByVectorIfFullAtomSpaceIsEmpty(Particle &ParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ)
    {
        try
        {
            if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, ParticlesInSector, VectorX, VectorY, VectorZ, true) == true)
                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH("moving particle by vector if voxel space is empty")

        return true;
    }
protected:
    inline bool MoveParticleByVectorIfFullAtomSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float VectorX, const float VectorY, const float VectorZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
    {
        try
        {
            if (CheckFreeSpaceAndBoundsForParticleMovedByVector(ParticleObject, ParticlesInSector, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, true) == true)
                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH_AND_THROW("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    inline void MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const ParticlesContainerInternal<Particle>& ParticlesInSector, const float AddX, const float AddY, const float AddZ)
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
                        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, ParticlesInSector, PosX, PosY, PosZ, true) == true)
                        {
                            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "FREE SPACE FOUND " << VecX << " " << VecY << " " << VecZ << " " << PosX << " " << PosY << " " << PosZ << terminal_colors_utils::white));
                            MoveParticleByVector(ParticleObject, PosX, PosY, PosZ);
                            FoundFreeSpace = true;
                            goto Outside;
                        }
            Outside:

            if (FoundFreeSpace == false)
                LoggersManagerObject.Log(STREAM(terminal_colors_utils::red << "FREE SPACE NOT FOUND " << VecX << " " << VecY << " " << VecZ << " " << terminal_colors_utils::white));
        }
        CATCH("moving particle near other particles")
    }
protected:
    inline void MoveParticleNearOtherParticle(Particle &ParticleObject, const Particle &NewPositionParticleObject, const float AddX, const float AddY, const float AddZ)
    {
        try
        {
            MoveParticleByVector(ParticleObject, NewPositionParticleObject.ListOfAtoms[0].X - ParticleObject.ListOfAtoms[0].X + AddX, NewPositionParticleObject.ListOfAtoms[0].Y - ParticleObject.ListOfAtoms[0].Y + AddY, NewPositionParticleObject.ListOfAtoms[0].Z - ParticleObject.ListOfAtoms[0].Z + AddZ);
        }
        CATCH("moving particle near other particles")
    }
};

#endif
