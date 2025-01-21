
#ifndef CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H
#define CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H

#include <imgui.h>
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
    inline ParticlesContainerInternal<Particle>& GetParticles()
    {
        if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSpace)
            return Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ];
        else
        {
            if (CurrentThreadIndex == 0)
                return Particles[0][0][0];
            else
                return ParticlesForThreads;
        }
    }
protected:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return GetParticles()[ParticleIndex];
    }
public:
    [[nodiscard]] UniqueIdInt GetFreeIndexesOfParticleSize() const
    {
        return FreeIndexesOfParticles.size();
    }
protected:
    void InitiateFreeParticleIndexes(const ParticlesContainerInternal<Particle>& LocalParticles, bool PrintInfo);
protected:
    inline UniqueIdInt GetNewFreeIndexOfParticle()
    {
        if (FreeIndexesOfParticles.empty() == false)
        {
            UniqueIdInt FreeIndexOfParticle = FreeIndexesOfParticles.top();
            FreeIndexesOfParticles.pop();
            return FreeIndexOfParticle;
        }
        else
        {
            LoggersManagerObject.Log(STREAM("Lack of new free indexes of particles"));
            return MaxParticleIndex + 1;
        }
    }
public:
    void SetCurrentSectorPos(const CurrentSectorPosType& CurrentSectorPosParam)
    {
        CurrentSectorPos = CurrentSectorPosParam;
    }
public:
    UniqueIdInt AddNewParticle(const Particle& ParticleParam)
    {
                                                                                                                        //if (ParticleParam.Center.X == 0 || ParticleParam.Center.Y == 0 || ParticleParam.Center.Z == 0)
                                                                                                                        if (ParticleParam.Index == 0)
                                                                                                                        {
                                                                                                                            std::cout << "AAAYYY" << std::endl;
                                                                                                                            //getchar();
                                                                                                                            //throw std::runtime_error("TTT");
                                                                                                                        }
        GetParticles()[ParticleParam.Index] = ParticleParam;
        return MaxParticleIndex = ParticleParam.Index;
    }
protected:
    virtual void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) = 0;
public:
    template <class T>
    void PreprocessData(bool UpdateParticleKindListOfVoxelsBool);
protected:
    template <class T>
    void GetMinMaxCoordinatesForAllParticles(bool UpdateParticleKindListOfVoxelsBool) const;
    template <class T>
    static void GetMinMaxOfCoordinates(T PosX, T PosY, T PosZ, T& XMinParam, T& XMaxParam, T& YMinParam, T& YMaxParam, T& ZMinParam, T& ZMaxParam);
    template <class T>
    static void UpdateParticleKindListOfVoxels(const Particle& ParticleObject, T ParticleXMin, T ParticleXMax, T ParticleYMin, T ParticleYMax, T ParticleZMin, T ParticleZMax);
public:
    template <class T>
    static void GetMinMaxCoordinatesForParticle(Particle& ParticleObject, bool UpdateParticleKindListOfVoxels);
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
