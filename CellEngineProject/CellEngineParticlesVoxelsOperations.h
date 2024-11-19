
#ifndef CELL_ENGINE_PARTICLES_VOXELS_OPERATIONS_H
#define CELL_ENGINE_PARTICLES_VOXELS_OPERATIONS_H

#include <shared_mutex>

#include "DestinationPlatform.h"
#include "TerminalColorsUtils.h"

#include "CellEngineParticle.h"
#include "CellEngineBasicVoxelsOperations.h"

class CellEngineParticlesVoxelsOperations : public CellEngineBasicVoxelsOperations
{
public:
    void SetAllVoxelsInListOfVoxelsToValueForOuterClass(std::vector<vector3_16> &ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
    {
        SetAllVoxelsInListOfVoxelsToValue(ListOfVoxels, SimulationSpaceVoxelValue);
    }
protected:
    static inline SimulationSpaceVoxel GetZeroSimulationSpaceVoxel()
    {
        return 0;
    }
protected:
    inline void SetAllVoxelsInListOfVoxelsToValue(std::vector<vector3_16> &ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
    {
        try
        {
            for (auto &Voxel: ListOfVoxels)
                GetSpaceVoxel(Voxel.X, Voxel.Y, Voxel.Z) = SimulationSpaceVoxelValue;
        }
        CATCH("making all zero voxels in list of voxels")
    }
protected:
    inline void MoveParticleByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
    {
        try
        {
            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
            {
                VoxelOfParticle.X += VectorX;
                VoxelOfParticle.Y += VectorY;
                VoxelOfParticle.Z += VectorZ;
            }
            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, ParticleObject.Index);

            ParticleObject.SetCenterCoordinates(ParticleObject.Center.X + VectorX, ParticleObject.Center.Y + VectorY, ParticleObject.Center.Z + VectorZ);
        }
        CATCH("moving particle by vector")
    }
protected:
    inline void MoveParticleNearOtherParticle(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
    {
        try
        {
            MoveParticleByVector(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
        }
        CATCH("moving particle near other particles")
    }
protected:
    inline bool CheckFreeSpaceForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
    {
        try
        {
            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
            {
                UnsignedInt TestedPosX = VoxelOfParticle.X + VectorX;
                UnsignedInt TestedPosY = VoxelOfParticle.Y + VectorY;
                UnsignedInt TestedPosZ = VoxelOfParticle.Z + VectorZ;
                if (GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != GetZeroSimulationSpaceVoxel() && GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != ParticleObject.Index)
                    return false;
            }
        }
        CATCH("checking free space for particle moved by vector")

        return true;
    }
protected:
    static inline bool CheckBoundsForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
    {
        try
        {
            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
                if (!(VoxelOfParticle.X + VectorX >= StartXPosParam && VoxelOfParticle.X + VectorX < StartXPosParam + SizeXParam && VoxelOfParticle.Y + VectorY >= StartYPosParam && VoxelOfParticle.Y + VectorY < StartYPosParam + SizeYParam && VoxelOfParticle.Z + VectorZ >= StartZPosParam && VoxelOfParticle.Z + VectorZ < StartZPosParam + SizeZParam))
                    return false;
        }
        CATCH("checking bounds for particle moved by vector")

        return true;
    }
protected:
    inline bool CheckBoundsAndFreeSpaceForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
    {
        try
        {
            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
            {
                UnsignedInt TestedPosX = VoxelOfParticle.X + VectorX;
                UnsignedInt TestedPosY = VoxelOfParticle.Y + VectorY;
                UnsignedInt TestedPosZ = VoxelOfParticle.Z + VectorZ;
                if ((GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != GetZeroSimulationSpaceVoxel() && GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != ParticleObject.Index) || !(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
                    return false;
            }
        }
        CATCH("checking bounds and free space for particle moved by vector")

        return true;
    }
protected:
    inline bool MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
    {
        return MoveParticleByVectorIfVoxelSpaceIsEmpty(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
    }
protected:
    inline bool MoveParticleByVectorIfVoxelSpaceIsEmpty(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
    {
        try
        {
            if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, VectorX, VectorY, VectorZ) == true)
                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH("moving particle by vector if voxel space is empty")

        return true;
    }
protected:
    inline bool MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
    {
        try
        {
            if (CheckBoundsAndFreeSpaceForParticleMovedByVector(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
            else
                return false;
        }
        CATCH("moving particle by vector if voxel space is empty and is in bounds")

        return true;
    }
protected:
    inline void MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
    {
        try
        {
            bool FoundFreeSpace = false;

            SignedInt VecX = NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X;
            SignedInt VecY = NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y;
            SignedInt VecZ = NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z;

            for (SignedInt PosX = VecX - AddX; PosX < VecX + AddX; PosX++)
                for (SignedInt PosY = VecY - AddY; PosY < VecY + AddY; PosY++)
                    for (SignedInt PosZ = VecZ - AddZ; PosZ < VecZ + AddZ; PosZ++)
                        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, PosX, PosY, PosZ) == true)
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
};

#endif
