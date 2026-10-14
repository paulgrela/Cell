
#ifndef CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H
#define CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineParticleKind.h"
#include "CellEngineParticleUniqueIdGenerator.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineBasicParallelExecutionData.h"

class CellEngineBasicParticlesOperations : public CellEngineBasicParallelExecutionData
{
protected:
    CellEngineParticleUniqueIdGenerator ParticleUniqueIdGenerator{};
protected:
    UnsignedInt MaxParticleUniqueId{};
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
public:
    inline Particle& GetParticleFromUniqueId(const UniqueIdUInt64Type ParticleUniqueId)
    {
        return GetParticles()[ParticleUniqueId];
    }
// protected:
//     inline CellEngineParticleUniqueIdGenerator& GetFreeParticleUniqueIdes()
//     {
//         if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
//             return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator;
//         else
//             return ParticleUniqueIdGenerator;
//     }
protected:
    void InitiateFreeParticlesUniqueIdsForAllSectors();
protected:
    inline UniqueIdUInt64Type GetNewParticleUniqueIdFinal()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewParticleUniqueId();
        else
        {
            if (CurrentThreadIndex == 0)
                return Particles[0][0][0].ParticleUniqueIdGenerator.GetNewParticleUniqueId();
            else
                return ParticleUniqueIdGenerator.GetNewParticleUniqueId();
        }
    }
    inline UniqueIdUInt64Type GetNewFreeParticleUniqueId()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
        {
            UniqueIdUInt64Type NewParticleUniqueId = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewParticleUniqueId();

            while (GetParticles().contains(NewParticleUniqueId) == true)
                NewParticleUniqueId = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].ParticleUniqueIdGenerator.GetNewParticleUniqueId();

            return NewParticleUniqueId;
        }
        else
        {
            if (CurrentThreadIndex == 0)
            {
                UniqueIdUInt64Type NewParticleUniqueId = Particles[0][0][0].ParticleUniqueIdGenerator.GetNewParticleUniqueId();

                while (GetParticles().contains(NewParticleUniqueId) == true)
                    NewParticleUniqueId = Particles[0][0][0].ParticleUniqueIdGenerator.GetNewParticleUniqueId();

                return NewParticleUniqueId;
            }
            else
            {
                UniqueIdUInt64Type NewParticleUniqueId = ParticleUniqueIdGenerator.GetNewParticleUniqueId();

                while (GetParticles().contains(NewParticleUniqueId) == true)
                    NewParticleUniqueId = ParticleUniqueIdGenerator.GetNewParticleUniqueId();

                return NewParticleUniqueId;
            }
        }
    }
public:
    void SetCurrentSectorPos(const SectorPosType& CurrentSectorPosParam)
    {
        CurrentSectorPos = CurrentSectorPosParam;
    }
public:
    UniqueIdUInt64Type AddNewParticle(const Particle& ParticleParam)
    {
        GetParticles()[ParticleParam.ParticleUniqueId] = ParticleParam;
        return MaxParticleUniqueId = ParticleParam.ParticleUniqueId;
    }
protected:
    virtual void RemoveParticle(UniqueIdUInt64Type ParticleUniqueId, bool ClearElements, bool CountAllParticlesBool) = 0;
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
    std::vector<UniqueIdUInt64Type> GetAllParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam) const;
    std::vector<UniqueIdUInt64Type> GetAllParticlesWithChosenParticleKindId(UniqueIdUInt64Type ParticleKindId) const;
    UnsignedInt GetNumberOfParticlesWithChosenParticleKindId(UniqueIdUInt64Type ParticleKindId) const;
protected:
    explicit CellEngineBasicParticlesOperations(ParticlesContainer<Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
    }
public:
    virtual ~CellEngineBasicParticlesOperations() = default;
};

#endif
