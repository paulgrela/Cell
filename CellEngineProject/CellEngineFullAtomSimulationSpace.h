
#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_H

#include "CellEngineAtom.h"
#include "CellEngineTypes.h"
#include "CellEngineParticle.h"

#include "CellEngineSimulationSpace.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

class CellEngineParticlesFullAtomOperations //: public CellEngineBasicVoxelsOperations
{
//public:
//    void SetAllVoxelsInListOfVoxelsToValueForOuterClass(std::vector<vector3_16> &ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
//    {
//        SetAllVoxelsInListOfVoxelsToValue(ListOfVoxels, SimulationSpaceVoxelValue);
//    }
//protected:
//    static inline SimulationSpaceVoxel GetZeroSimulationSpaceVoxel()
//    {
//        return 0;
//    }
//protected:
//    inline void SetAllVoxelsInListOfVoxelsToValue(std::vector<vector3_16> &ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
//    {
//        try
//        {
//            for (auto &Voxel: ListOfVoxels)
//                GetSpaceVoxel(Voxel.X, Voxel.Y, Voxel.Z) = SimulationSpaceVoxelValue;
//        }
//        CATCH("making all zero voxels in list of voxels")
//    }
//protected:
//    inline void MoveParticleByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
//    {
//        try
//        {
//            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
//            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
//            {
//                VoxelOfParticle.X += VectorX;
//                VoxelOfParticle.Y += VectorY;
//                VoxelOfParticle.Z += VectorZ;
//            }
//            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, ParticleObject.Index);
//        }
//        CATCH("moving particle by vector")
//    }
//protected:
//    inline void MoveParticleNearOtherParticle(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
//    {
//        try
//        {
//            MoveParticleByVector(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
//        }
//        CATCH("moving particle near other particles")
//    }
//protected:
//    inline bool CheckFreeSpaceForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
//    {
//        try
//        {
//            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
//                if (GetSpaceVoxel(VoxelOfParticle.X + VectorX, VoxelOfParticle.Y + VectorY, VoxelOfParticle.Z + VectorZ) != GetZeroSimulationSpaceVoxel())
//                    return false;
//        }
//        CATCH("checking free space for particle moved by vector")
//
//        return true;
//    }
//protected:
//    static inline bool CheckBoundsForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
//    {
//        try
//        {
//            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
//                if (!(VoxelOfParticle.X + VectorX >= StartXPosParam && VoxelOfParticle.X + VectorX < StartXPosParam + SizeXParam && VoxelOfParticle.Y + VectorY >= StartYPosParam && VoxelOfParticle.Y + VectorY < StartYPosParam + SizeYParam && VoxelOfParticle.Z + VectorZ >= StartZPosParam && VoxelOfParticle.Z + VectorZ < StartZPosParam + SizeZParam))
//                    return false;
//        }
//        CATCH("checking bounds for particle moved by vector")
//
//        return true;
//    }
//protected:
//    inline bool CheckBoundsAndFreeSpaceForParticleMovedByVector(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
//    {
//        try
//        {
//            for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
//            {
//                UnsignedInt TestedPosX = VoxelOfParticle.X + VectorX;
//                UnsignedInt TestedPosY = VoxelOfParticle.Y + VectorY;
//                UnsignedInt TestedPosZ = VoxelOfParticle.Z + VectorZ;
//                if ((GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != GetZeroSimulationSpaceVoxel() && GetSpaceVoxel(TestedPosX, TestedPosY, TestedPosZ) != ParticleObject.Index) || !(TestedPosX >= StartXPosParam && TestedPosX < StartXPosParam + SizeXParam && TestedPosY >= StartYPosParam && TestedPosY < StartYPosParam + SizeYParam && TestedPosZ >= StartZPosParam && TestedPosZ < StartZPosParam + SizeZParam))
//                    return false;
//            }
//        }
//        CATCH("checking bounds and free space for particle moved by vector")
//
//        return true;
//    }
//protected:
//    inline bool MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
//    {
//        return MoveParticleByVectorIfVoxelSpaceIsEmpty(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
//    }
//protected:
//    inline bool MoveParticleByVectorIfVoxelSpaceIsEmpty(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ)
//    {
//        try
//        {
//            if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, VectorX, VectorY, VectorZ) == true)
//                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
//            else
//                return false;
//        }
//        CATCH("moving particle by vector if voxel space is empty")
//
//        return true;
//    }
//protected:
//    inline bool MoveParticleByVectorIfVoxelSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, const SignedInt VectorX, const SignedInt VectorY, const SignedInt VectorZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
//    {
//        try
//        {
//            if (CheckBoundsAndFreeSpaceForParticleMovedByVector(ParticleObject, VectorX, VectorY, VectorZ, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
//                MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
//            else
//                return false;
//        }
//        CATCH("moving particle by vector if voxel space is empty and is in bounds")
//
//        return true;
//    }
//protected:
//    inline void MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
//    {
//        try
//        {
//            bool FoundFreeSpace = false;
//
//            SignedInt VecX = NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X;
//            SignedInt VecY = NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y;
//            SignedInt VecZ = NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z;
//
//            for (SignedInt PosX = VecX - AddX; PosX < VecX + AddX; PosX++)
//                for (SignedInt PosY = VecY - AddY; PosY < VecY + AddY; PosY++)
//                    for (SignedInt PosZ = VecZ - AddZ; PosZ < VecZ + AddZ; PosZ++)
//                        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, PosX, PosY, PosZ) == true)
//                        {
//                            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "FREE SPACE FOUND " << VecX << " " << VecY << " " << VecZ << " " << PosX << " " << PosY << " " << PosZ << terminal_colors_utils::white));
//                            MoveParticleByVector(ParticleObject, PosX, PosY, PosZ);
//                            FoundFreeSpace = true;
//                            goto Outside;
//                        }
//            Outside:
//
//            if (FoundFreeSpace == false)
//                LoggersManagerObject.Log(STREAM(terminal_colors_utils::red << "FREE SPACE NOT FOUND " << VecX << " " << VecY << " " << VecZ << " " << terminal_colors_utils::white));
//        }
//        CATCH("moving particle near other particles")
//    }
};

class CellEngineParticlesFullAtomShapesGenerator : virtual public CellEngineParticlesFullAtomOperations
{
protected:
    virtual Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) = 0;
    virtual void ClearVoxelSpaceAndParticles() = 0;
protected:
    typedef bool (CellEngineParticlesFullAtomShapesGenerator::*CheckFreeSpaceForSelectedSpaceType)(UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UniqueIdInt );
    typedef void (CellEngineParticlesFullAtomShapesGenerator::*SetValueToVoxelsForSelectedSpaceType)(std::vector<vector3_16>*, UniqueIdInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt);
protected:
    bool GenerateParticleVoxelsWhenSelectedSpaceIsFree(UnsignedInt LocalNewParticleIndex, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToVoxelsForSelectedSpaceType SetValueToVoxelsForSelectedSpace);
protected:
    void SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(std::vector <vector3_16> *FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ);
public:
    bool CheckFreeSpaceInCuboidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt SizeOfParticleX, UnsignedInt SizeOfParticleY, UnsignedInt SizeOfParticleZ, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForSphereSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
    bool CheckFreeSpaceForEllipsoidSelectedSpace(UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam, UniqueIdInt ValueToCheck);
public:
    void SetValueToVoxelsForCuboidSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
    void SetValueToVoxelsForSphereSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
    void SetValueToVoxelsForEllipsoidSelectedSpace(std::vector<vector3_16>* FilledSpaceVoxels, UniqueIdInt VoxelValue, UnsignedInt PosXStart, UnsignedInt PosYStart, UnsignedInt PosZStart, UnsignedInt StepX, UnsignedInt StepY, UnsignedInt StepZ, UnsignedInt RadiusXParam, UnsignedInt RadiusYParam, UnsignedInt RadiusZParam);
protected:
    explicit CellEngineParticlesFullAtomShapesGenerator(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam)
    {
    }
};

class CellEngineChemicalReactionsInFullAtomSimulationSpace : virtual public CellEngineChemicalReactionsInSimulationSpace, virtual public CellEngineParticlesFullAtomOperations
{
protected:
    void ClearSpaceForParticle(Particle& ParticleObject, bool ClearVoxels) override;
    void MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, SignedInt AddX, SignedInt AddY, SignedInt AddZ) override;
    void FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::set<UnsignedInt> &FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
protected:
    explicit CellEngineChemicalReactionsInFullAtomSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : CellEngineChemicalReactionsInSimulationSpace(ParticlesParam), CellEngineBasicParticlesOperations(ParticlesParam)
    {
    }
};

class CellEngineFullAtomSimulationSpace : public CellEngineSimulationSpace, public CellEngineChemicalReactionsInFullAtomSimulationSpace, public CellEngineVoxelSimulationSpaceStatistics
{
public:
    SimulationSpaceVoxel GetSpaceVoxelForOuterClass(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
    Particle& GetParticleFromIndexForOuterClass(UniqueIdInt ParticleIndex);
private:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    void SetAtomInFullAtomSimulationSpace(UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom);
    //Particle& GetParticleFromIndexForGenerator(UniqueIdInt ParticleIndex) override;
public:
    //void ClearVoxelSpaceAndParticles() override;
public:
//    void ClearWholeVoxelSpace();
//    void ClearSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam);
protected:
    void FillParticleElementInSpace(UniqueIdInt ParticleIndex, vector3_64 NewVoxel) override;
    bool MoveParticleByVectorIfSpaceIsEmptyAndIsInBounds(Particle &ParticleObject, SignedInt VectorX, SignedInt VectorY, SignedInt VectorZ, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) override;
public:
    explicit CellEngineFullAtomSimulationSpace(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam);
    ~CellEngineFullAtomSimulationSpace();
};


#endif
