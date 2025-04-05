
#ifndef CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H
#define CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H

#include <stack>
#include <shared_mutex>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineParticleKind.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParallelExecutionData.h"

class CellEngineBasicParticlesOperations : public CellEngineBasicParallelExecutionData
{
protected:
    UnsignedInt MaxParticleIndex{};
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
protected:
    ParticlesContainer<Particle>& Particles;
protected:
    inline ParticlesDetailedContainer<Particle>& GetParticles()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].Particles;
        else
        {
            if (CurrentThreadIndex == 0)
                return Particles[0][0][0].Particles;
            else
                return ParticlesForThreads;
        }
    }
protected:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return GetParticles()[ParticleIndex];
    }
protected:
    inline std::stack<UniqueIdInt>& GetFreeIndexes()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].FreeIndexesOfParticles;
        else
            return FreeIndexesOfParticles;
    }
public:
    [[nodiscard]] UniqueIdInt GetFreeIndexesOfParticleSize() const
    {
        return FreeIndexesOfParticles.size();
    }
protected:
    void InitiateFreeParticleIndexes(const ParticlesDetailedContainer<Particle>& LocalParticles, bool PrintInfo);
protected:
    inline UniqueIdInt GetNewFreeIndexOfParticle()
    {
        if (GetFreeIndexes().empty() == false)
        {
            const UniqueIdInt FreeIndexOfParticle = GetFreeIndexes().top();
            GetFreeIndexes().pop();
            return FreeIndexOfParticle;
        }
        else
        {
            LoggersManagerObject.Log(STREAM("Lack of new free indexes of particles"));
            return MaxParticleIndex + 1;
        }
    }
public:
    void SetCurrentSectorPos(const SectorPosType& CurrentSectorPosParam)
    {
        CurrentSectorPos = CurrentSectorPosParam;
    }
public:
    UniqueIdInt AddNewParticle(const Particle& ParticleParam)
    {
        GetParticles()[ParticleParam.Index] = ParticleParam;
        return MaxParticleIndex = ParticleParam.Index;
    }
protected:
    virtual void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearElements) = 0;
public:
    template <class T, class A>
    void PreprocessData(const std::vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool);
protected:
    template <class T, class A>
    void GetMinMaxCoordinatesForAllParticles(const std::vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElementsBool) const;
    template <class T>
    static void GetMinMaxOfCoordinates(T PosX, T PosY, T PosZ, T& XMinParam, T& XMaxParam, T& YMinParam, T& YMaxParam, T& ZMinParam, T& ZMaxParam);
    template <class T, class A>
    static void UpdateParticleKindListOfElements(const Particle& ParticleObject, const std::vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, T ParticleXMin, T ParticleXMax, T ParticleYMin, T ParticleYMax, T ParticleZMin, T ParticleZMax, T XSizeDiv2, T YSizeDiv2, T ZSizeDiv2);
public:
    template <class T, class A>
    static void GetMinMaxCoordinatesForParticle(Particle& ParticleObject, const std::vector<A> Particle::*ListOfElements, const std::vector<A> ParticleKind::*ListOfElementsOfParticleKind, bool UpdateParticleKindListOfElements);
protected:
    std::vector<UniqueIdInt> GetAllParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam) const;
    std::vector<UniqueIdInt> GetAllParticlesWithChosenEntityId(UniqueIdInt EntityId) const;
    UnsignedInt GetNumberOfParticlesWithChosenEntityId(UniqueIdInt EntityId) const;
protected:
    explicit CellEngineBasicParticlesOperations(ParticlesContainer<Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
    }
public:
    virtual ~CellEngineBasicParticlesOperations() = default;
};

#endif
