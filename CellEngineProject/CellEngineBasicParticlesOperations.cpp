
#include <vector>

#include "CellEngineMacros.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineBasicParticlesOperations.h"

using namespace std;

void  CellEngineBasicParticlesOperations::InitiateFreeParticleIndexes(const ParticlesContainerInternal<Particle>& LocalParticles, const bool PrintInfo)
{
    try
    {
        if (PrintInfo == true)
            LoggersManagerObject.Log(STREAM("Scope of particle indexes for current thread = (" << to_string((CurrentThreadIndex + 1) * ParticleIndexesCreatorFactor - 1) << " , " << to_string(CurrentThreadIndex * ParticleIndexesCreatorFactor) << ") LocalParticles.size() = " << LocalParticles.size()));

        FreeIndexesOfParticles = {};

        for (UniqueIdInt FreeIndex = (CurrentThreadIndex + 1) * ParticleIndexesCreatorFactor - 1; FreeIndex > CurrentThreadIndex * ParticleIndexesCreatorFactor; FreeIndex--)
            if (!LocalParticles.contains(FreeIndex))
                FreeIndexesOfParticles.push(FreeIndex);

        if (PrintInfo == true)
            LoggersManagerObject.Log(STREAM("FreeIndexesOfParticles.size() = " << FreeIndexesOfParticles.size()));
    }
    CATCH("initiating free particle indexes")
}

template <class T, class A>
void CellEngineBasicParticlesOperations::PreprocessData(vector<A> Particle::*ListOfElements, std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool)
{
    try
    {
        LoggersManagerObject.Log(STREAM("Preprocess data"));

        InitiateFreeParticleIndexes(Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ], true);
        GetMinMaxCoordinatesForAllParticles<T, A>(ListOfElements, ListOfElementsOfParticleKind, UpdateParticleKindListOfElementsBool);
    }
    CATCH("preprocessing data for voxel simulation space")
}

template <class T>
void CellEngineBasicParticlesOperations::GetMinMaxOfCoordinates(const T PosX, const T PosY, const T PosZ, T& XMinParam, T& XMaxParam, T& YMinParam, T& YMaxParam, T& ZMinParam, T& ZMaxParam)
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

template <class T, class A>
void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForAllParticles(vector<A> Particle::*ListOfElements, std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool) const
{
    try
    {
        LoggersManagerObject.Log(STREAM("GET MIN MAX FOR ALL PARTICLES"));

        FOR_EACH_PARTICLE_IN_XYZ
            if (ParticleObject.second.EntityId != 0)
                GetMinMaxCoordinatesForParticle<T, A>(ParticleObject.second, ListOfElements, ListOfElementsOfParticleKind, UpdateParticleKindListOfElementsBool);
    }
    CATCH("getting min max coordinates for all particles")
}

template <class T, class A>
void CellEngineBasicParticlesOperations::UpdateParticleKindListOfElements(const Particle& ParticleObject, vector<A> Particle::*ListOfElements, vector<A> ParticleKind::*ListOfElementsOfParticleKind, const T ParticleXMin, const T ParticleXMax, const T ParticleYMin, const T ParticleYMax, const T ParticleZMin, const T ParticleZMax)
{
    try
    {
        auto& ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        if (static_cast<vector<A>>(ParticleKindObject.*ListOfElementsOfParticleKind).empty() == true)
        {
            for (const auto& VoxelCoordinates : ParticleObject.*ListOfElements)
                static_cast<vector<A>>(ParticleKindObject.*ListOfElementsOfParticleKind).emplace_back(VoxelCoordinates.X - ParticleXMin, VoxelCoordinates.Y - ParticleYMin, VoxelCoordinates.Z - ParticleZMin);

            ParticleKindObject.XSizeDiv2 = static_cast<float>(ParticleXMax - ParticleXMin) / 2;
            ParticleKindObject.YSizeDiv2 = static_cast<float>(ParticleYMax - ParticleYMin) / 2;
            ParticleKindObject.ZSizeDiv2 = static_cast<float>(ParticleZMax - ParticleZMin) / 2;
        }
    }
    CATCH("updating particle kind list of Elements")
}

template <class T, class A>
void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle(Particle& ParticleObject, vector<A> Particle::*ListOfElements, vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool)
{
    try
    {
        T ParticleXMin{}, ParticleXMax{}, ParticleYMin{}, ParticleYMax{}, ParticleZMin{}, ParticleZMax{};
        ParticleXMin = ParticleYMin = ParticleZMin = 10000;
        if constexpr(std::is_same_v<T, UnsignedInt>)
            ParticleXMax = ParticleYMax = ParticleZMax = 0;
        else
            ParticleXMax = ParticleYMax = ParticleZMax = -10000;

        for (const auto& VoxelCoordinates : ParticleObject.*ListOfElements)
            GetMinMaxOfCoordinates<T>(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);

        ParticleObject.SetCenterCoordinates(ParticleXMin + static_cast<float>(ParticleXMax - ParticleXMin) / 2, ParticleYMin + static_cast<float>(ParticleYMax - ParticleYMin) / 2, ParticleZMin + static_cast<float>(ParticleZMax - ParticleZMin) / 2);

        if (UpdateParticleKindListOfElementsBool == true)
            UpdateParticleKindListOfElements<T, A>(ParticleObject, ListOfElements, ListOfElementsOfParticleKind, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);
    }
    CATCH("getting min max coordinates for one particle")
}


vector<UniqueIdInt> CellEngineBasicParticlesOperations::GetAllParticlesWithChosenParticleType(const ParticlesTypes ParticleTypeParam) const
{
    vector<UniqueIdInt> ListOfParticlesIndexes;

    try
    {
        ListOfParticlesIndexes.reserve(1024 * 1024);
        FOR_EACH_PARTICLE_IN_XYZ_CONST
            if (ParticleObject.second.EntityId != 0)
            {
                auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.second.EntityId);
                if (ParticleKindObject.ParticleKindSpecialDataSector.empty() == false)
                    for (const auto& ParticleKindSpecialDataObject : ParticleKindObject.ParticleKindSpecialDataSector)
                        if (ParticleKindSpecialDataObject.ParticleType == ParticleTypeParam)
                            ListOfParticlesIndexes.emplace_back(ParticleObject.first);
            }
    }
    CATCH("getting all particles with chosen entity id")

    return ListOfParticlesIndexes;
}

vector<UniqueIdInt> CellEngineBasicParticlesOperations::GetAllParticlesWithChosenEntityId(const UniqueIdInt EntityId) const
{
    vector<UniqueIdInt> ListOfParticlesIndexes;

    try
    {
        ListOfParticlesIndexes.reserve(1024 * 1024);
        FOR_EACH_PARTICLE_IN_XYZ_CONST
            if (ParticleObject.second.EntityId == EntityId)
                ListOfParticlesIndexes.emplace_back(ParticleObject.first);
    }
    CATCH("getting all particles with chosen entity id")

    return ListOfParticlesIndexes;
}

UnsignedInt CellEngineBasicParticlesOperations::GetNumberOfParticlesWithChosenEntityId(const UniqueIdInt EntityId) const
{
    UnsignedInt ParticleCounter = 0;

    try
    {
        FOR_EACH_PARTICLE_IN_XYZ_CONST
            if (ParticleObject.second.EntityId == EntityId)
                ParticleCounter++;
    }
    CATCH("getting number of all particles with chosen entity id")

    return ParticleCounter;
}

template void CellEngineBasicParticlesOperations::PreprocessData<float, CellEngineAtom>(std::vector<CellEngineAtom> Particle::*ListOfElements, std::vector<CellEngineAtom> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::PreprocessData<UnsignedInt, vector3_16>(std::vector<vector3_16> Particle::*ListOfElements, std::vector<vector3_16> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::GetMinMaxOfCoordinates<UnsignedInt>(UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam);
template void CellEngineBasicParticlesOperations::UpdateParticleKindListOfElements<UnsignedInt, vector3_16>(const Particle& ParticleObject, vector<vector3_16> Particle::*ListOfElements, vector<vector3_16> ParticleKind::*ListOfElementsOfParticleKind, const UnsignedInt ParticleXMin, const UnsignedInt ParticleXMax, const UnsignedInt ParticleYMin, const UnsignedInt ParticleYMax, const UnsignedInt ParticleZMin, const UnsignedInt ParticleZMax);
template void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(Particle& ParticleObject, vector<vector3_16> Particle::*ListOfElements, vector<vector3_16> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<float, CellEngineAtom>(Particle& ParticleObject, vector<CellEngineAtom> Particle::*ListOfElements, vector<CellEngineAtom> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool);
