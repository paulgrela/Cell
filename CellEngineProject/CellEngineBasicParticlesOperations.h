
#ifndef CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H
#define CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H

#include <stack>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineParticleKind.h"
#include "CellEngineParticleUniqueIdGenerator.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParallelExecutionData.h"

class CellEngineBasicParticlesOperations : public CellEngineBasicParallelExecutionData
{
protected:
    CellEngineParticleUniqueIdGenerator ParticleUniqueIdGenerator;
protected:
    UnsignedInt MaxParticleIndex{};//USUNAC
    //std::stack<UniqueIdUnsignedInt> FreeIndexesOfParticlesGlobal;
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
    inline Particle& GetParticleFromIndex(const UniqueIdUnsignedInt ParticleIndex)
    {
        return GetParticles()[ParticleIndex];
    }
protected:

    // inline std::stack<UniqueIdUnsignedInt>& GetFreeIndexes()
    // {
    //     if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
    //         return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].FreeIndexesOfParticles;
    //     else
    //         return FreeIndexesOfParticlesGlobal;
    // }
    inline CellEngineParticleUniqueIdGenerator& GetFreeParticleIndexes()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator;
        else
            return ParticleUniqueIdGenerator;
    }
public:
    // [[nodiscard]] UniqueIdUnsignedInt GetFreeIndexesOfParticleSize() const
    // {
    //     return FreeIndexesOfParticlesGlobal.size();
    // }
protected:
    void InitiateFreeParticleIndexesForAllSectors();
    // void InitiateFreeParticleIndexes(const ParticlesDetailedContainer<Particle>& LocalParticles, bool PrintInfo);
protected:


    inline UniqueIdUnsignedInt GetNewFreeIndexOfParticleFinal()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();
        else
        {
            if (CurrentThreadIndex == 0)
                return Particles[0][0][0].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();
            else
                return ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();
        }
    }
    inline UniqueIdUnsignedInt GetNewFreeIndexOfParticle()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
        {
            UniqueIdUnsignedInt NewUniqueParticleIndex = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

            while (GetParticles().contains(NewUniqueParticleIndex) == true)
                NewUniqueParticleIndex = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

            return NewUniqueParticleIndex;
        }
        else
        {
            if (CurrentThreadIndex == 0)
            {
                UniqueIdUnsignedInt NewUniqueParticleIndex = Particles[0][0][0].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

                while (GetParticles().contains(NewUniqueParticleIndex) == true)
                    NewUniqueParticleIndex = Particles[0][0][0].ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

                return NewUniqueParticleIndex;
            }
            else
            {
                UniqueIdUnsignedInt NewUniqueParticleIndex = ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

                while (GetParticles().contains(NewUniqueParticleIndex) == true)
                    NewUniqueParticleIndex = ParticleUniqueIdGenerator.GetNewUniqueParticleIndex();

                return NewUniqueParticleIndex;
            }
        }
    }
    // inline UniqueIdUnsignedInt GetNewFreeIndexOfParticleReuse()
    // {
    //     if (GetFreeIndexes().empty() == false)
    //     {
    //         const UniqueIdUnsignedInt FreeIndexOfParticle = GetFreeIndexes().top();
    //         GetFreeIndexes().pop();
    //         return FreeIndexOfParticle;
    //     }
    //     else
    //     {
    //         LoggersManagerObject.Log(STREAM("Lack of new free indexes of particles"));
    //         return MaxParticleIndex + 1;
    //     }
    // }


public:
    void SetCurrentSectorPos(const SectorPosType& CurrentSectorPosParam)
    {
        CurrentSectorPos = CurrentSectorPosParam;
    }
public:
    UniqueIdUnsignedInt AddNewParticle(const Particle& ParticleParam)
    {
        GetParticles()[ParticleParam.Index] = ParticleParam;
        return MaxParticleIndex = ParticleParam.Index;
    }
protected:
    virtual void RemoveParticle(UniqueIdUnsignedInt ParticleIndex, bool ClearElements) = 0;
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
    std::vector<UniqueIdUnsignedInt> GetAllParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam) const;
    std::vector<UniqueIdUnsignedInt> GetAllParticlesWithChosenEntityId(UniqueIdUnsignedInt EntityId) const;
    UnsignedInt GetNumberOfParticlesWithChosenEntityId(UniqueIdUnsignedInt EntityId) const;
protected:
    explicit CellEngineBasicParticlesOperations(ParticlesContainer<Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
    }
public:
    virtual ~CellEngineBasicParticlesOperations() = default;
};

#endif
