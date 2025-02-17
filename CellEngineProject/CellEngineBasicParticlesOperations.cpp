
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
void CellEngineBasicParticlesOperations::PreprocessData(const vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool)
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
void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForAllParticles(const vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool) const
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
void CellEngineBasicParticlesOperations::UpdateParticleKindListOfElements(const Particle& ParticleObject, const vector<A> Particle::*ListOfElements, const vector<A> ParticleKind::*ListOfElementsOfParticleKind, const T ParticleXMin, const T ParticleXMax, const T ParticleYMin, const T ParticleYMax, const T ParticleZMin, const T ParticleZMax, const T XSizeDiv2, const T YSizeDiv2, const T ZSizeDiv2)
{
    try
    {
        auto& ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        if (static_cast<vector<A>>(ParticleKindObject.*ListOfElementsOfParticleKind).empty() == true)
        {
            for (const auto& ElementCoordinates : ParticleObject.*ListOfElements)
                static_cast<vector<A>>(ParticleKindObject.*ListOfElementsOfParticleKind).emplace_back(ElementCoordinates.X - ParticleXMin, ElementCoordinates.Y - ParticleYMin, ElementCoordinates.Z - ParticleZMin);

            ParticleKindObject.XSizeDiv2 = XSizeDiv2;
            ParticleKindObject.YSizeDiv2 = YSizeDiv2;
            ParticleKindObject.ZSizeDiv2 = ZSizeDiv2;

            ParticleKindObject.Radius = ParticleObject.Radius;
        }
    }
    CATCH("updating particle kind list of Elements")
}

template <class T, class A>
void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle(Particle& ParticleObject, const vector<A> Particle::*ListOfElements, const vector<A> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool)
{
    try
    {
        T ParticleXMin{}, ParticleXMax{}, ParticleYMin{}, ParticleYMax{}, ParticleZMin{}, ParticleZMax{};
        ParticleXMin = ParticleYMin = ParticleZMin = 10000;
        if constexpr(std::is_same_v<T, UnsignedInt>)
            ParticleXMax = ParticleYMax = ParticleZMax = 0;
        else
            ParticleXMax = ParticleYMax = ParticleZMax = -10000;

        UnsignedInt ElementsCounter = 0;
        for (const auto& ElementCoordinates : ParticleObject.*ListOfElements)
        {
            GetMinMaxOfCoordinates<T>(ElementCoordinates.X, ElementCoordinates.Y, ElementCoordinates.Z, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);
            ElementsCounter++;
        }

        T XSizeDiv2 = (ParticleXMax - ParticleXMin) / 2;
        T YSizeDiv2 = (ParticleYMax - ParticleYMin) / 2;
        T ZSizeDiv2 = (ParticleZMax - ParticleZMin) / 2;

        ParticleObject.Radius = max({ XSizeDiv2, YSizeDiv2, ZSizeDiv2 });

        ParticleObject.SetCenterCoordinates(ParticleXMin + XSizeDiv2, ParticleYMin + YSizeDiv2, ParticleZMin + ZSizeDiv2);
        
        if (UpdateParticleKindListOfElementsBool == true)
            UpdateParticleKindListOfElements<T, A>(ParticleObject, ListOfElements, ListOfElementsOfParticleKind, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax, XSizeDiv2, YSizeDiv2, ZSizeDiv2);
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

template void CellEngineBasicParticlesOperations::PreprocessData<UnsignedInt, vector3_16>(const std::vector<vector3_16> Particle::*ListOfElements, const std::vector<vector3_16> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::PreprocessData<RealType, CellEngineAtom>(const std::vector<CellEngineAtom> Particle::*ListOfElements, const std::vector<CellEngineAtom> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::GetMinMaxOfCoordinates<UnsignedInt>(UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam);
template void CellEngineBasicParticlesOperations::GetMinMaxOfCoordinates<RealType>(RealType PosX, RealType PosY, RealType PosZ, RealType& XMinParam, RealType& XMaxParam, RealType& YMinParam, RealType& YMaxParam, RealType& ZMinParam, RealType& ZMaxParam);
template void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(Particle& ParticleObject, const vector<vector3_16> Particle::*ListOfElements, const vector<vector3_16> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool);
template void CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<RealType, CellEngineAtom>(Particle& ParticleObject, const vector<CellEngineAtom> Particle::*ListOfElements, const vector<CellEngineAtom> ParticleKind::*ListOfElementsOfParticleKind, const bool UpdateParticleKindListOfElementsBool);
