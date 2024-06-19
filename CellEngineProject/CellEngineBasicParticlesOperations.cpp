
#include <vector>

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineBasicParticlesOperations.h"

using namespace std;

void  CellEngineBasicParticlesOperations::InitiateFreeParticleIndexes()
{
    try
    {
        FreeIndexesOfParticles = {};

        for (UnsignedInt FreeIndex = MaxParticleIndex + 100'000'000; FreeIndex >= MaxParticleIndex + 1; FreeIndex--)
            if (Particles.find(FreeIndex) == Particles.end())
                FreeIndexesOfParticles.push(FreeIndex);
    }
    CATCH("initiating free particle indexes")
}

void CellEngineBasicParticlesOperations::PreprocessData(const bool UpdateParticleKindListOfVoxelsBool)
{
    try
    {
        InitiateFreeParticleIndexes();
        GetMinMaxCoordinatesForAllParticles(UpdateParticleKindListOfVoxelsBool);
    }
    CATCH("preprocessing data for voxel simulation space")
}

void CellEngineBasicParticlesOperations::SetStartValuesForSpaceMinMax()
{
    XMin = YMin = ZMin = 10000;
    XMax = YMax = ZMax = 0;
}

void CellEngineBasicParticlesOperations::GetMinMaxOfCoordinates(const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam)
{
    try
    {
        XMinParam = min(PosX, XMinParam);
        XMaxParam = max(PosX, XMaxParam);
        YMinParam = min(PosY, YMinParam);
        YMaxParam = max(PosY, YMaxParam);
        ZMinParam = min(PosZ, ZMinParam);
        ZMaxParam = max(PosZ, ZMaxParam);
    }
    CATCH("getting min max of coordinates")
}

void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForAllParticles(const bool UpdateParticleKindListOfVoxelsBool)
{
    try
    {
        for (auto& ParticleObject : Particles)
                                                        if (ParticleObject.second.EntityId != 0)
            GetMinMaxCoordinatesForParticle(ParticleObject.second, UpdateParticleKindListOfVoxelsBool);
    }
    CATCH("getting min max coordinates for all particles")
}

void CellEngineBasicParticlesOperations::UpdateParticleKindListOfVoxels(Particle& ParticleObject, const UnsignedInt ParticleXMin, const UnsignedInt ParticleXMax, const UnsignedInt ParticleYMin, const UnsignedInt ParticleYMax, const UnsignedInt ParticleZMin, const UnsignedInt ParticleZMax)
{
    try
    {
        auto& ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        if (ParticleKindObject.ListOfVoxels.empty() == true)
        {
            for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
                ParticleKindObject.ListOfVoxels.emplace_back(VoxelCoordinates.X - ParticleXMin, VoxelCoordinates.Y - ParticleYMin, VoxelCoordinates.Z - ParticleZMin);

            ParticleKindObject.XSizeDiv2 = (ParticleXMax - ParticleXMin) / 2;
            ParticleKindObject.YSizeDiv2 = (ParticleYMax - ParticleYMin) / 2;
            ParticleKindObject.ZSizeDiv2 = (ParticleZMax - ParticleZMin) / 2;
        }
    }
    CATCH("updating particle kind list of voxels")
}

void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle(Particle& ParticleObject, const bool UpdateParticleKindListOfVoxelsBool)
{
    try
    {
        UnsignedInt ParticleXMin{}, ParticleXMax{}, ParticleYMin{}, ParticleYMax{}, ParticleZMin{}, ParticleZMax{};
        ParticleXMin = ParticleYMin = ParticleZMin = 10000;
        ParticleXMax = ParticleYMax = ParticleZMax = 0;

        for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
            GetMinMaxOfCoordinates(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);

        ParticleObject.SetCenterCoordinates(ParticleXMin + (ParticleXMax - ParticleXMin) / 2, ParticleYMin + (ParticleYMax - ParticleYMin) / 2, ParticleZMin + (ParticleZMax - ParticleZMin) / 2);

        if (UpdateParticleKindListOfVoxelsBool == true)
            UpdateParticleKindListOfVoxels(ParticleObject, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);
    }
    CATCH("getting min max coordinates for one particle")
}

vector<UniqueIdInt> CellEngineBasicParticlesOperations::GetAllParticlesWithChosenEntityId(const UniqueIdInt EntityId)
{
    vector<UniqueIdInt> ListOfParticlesIndexes;

    try
    {
        ListOfParticlesIndexes.reserve(1024 * 1024);
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == EntityId)
                ListOfParticlesIndexes.emplace_back(ParticleObject.first);
    }
    CATCH("getting all particles with chosen entity id")

    return ListOfParticlesIndexes;
}

UnsignedInt CellEngineBasicParticlesOperations::GetNumberOfParticlesWithChosenEntityId(const UniqueIdInt EntityId)
{
    UnsignedInt ParticleCounter = 0;

    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == EntityId)
                ParticleCounter++;
    }
    CATCH("getting number of all particles with chosen entity id")

    return ParticleCounter;
}
