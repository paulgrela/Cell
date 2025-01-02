
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
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
protected:
    UnsignedInt MaxParticleIndex{};
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
protected:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
protected:
    inline std::unordered_map<UniqueIdInt, Particle>& GetParticles()
    {
        if (CurrentThreadIndex == 0)
            return Particles;
        else
            return ParticlesForThreads;
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
    void InitiateFreeParticleIndexes(const std::unordered_map<UniqueIdInt, Particle>& LocalParticles, bool PrintInfo);
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
    UniqueIdInt AddNewParticle(const Particle& ParticleParam)
    {
        GetParticles()[ParticleParam.Index] = ParticleParam;
        return MaxParticleIndex = ParticleParam.Index;
    }
protected:
    virtual void RemoveParticle(UniqueIdInt ParticleIndex, bool ClearVoxels) = 0;
public:
    void PreprocessData(bool UpdateParticleKindListOfVoxelsBool);
protected:
    void SetStartValuesForSpaceMinMax();
    void GetMinMaxCoordinatesForAllParticles(bool UpdateParticleKindListOfVoxelsBool);
    static void GetMinMaxOfCoordinates(UnsignedInt PosX, UnsignedInt PosY, UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam);
    static void UpdateParticleKindListOfVoxels(const Particle& ParticleObject, UnsignedInt ParticleXMin, UnsignedInt ParticleXMax, UnsignedInt ParticleYMin, UnsignedInt ParticleYMax, UnsignedInt ParticleZMin, UnsignedInt ParticleZMax);
public:
    static void GetMinMaxCoordinatesForParticle(Particle& ParticleObject, bool UpdateParticleKindListOfVoxels);
protected:
    std::vector<UniqueIdInt> GetAllParticlesWithChosenParticleType(ParticlesTypes ParticleTypeParam);
    std::vector<UniqueIdInt> GetAllParticlesWithChosenEntityId(UniqueIdInt EntityId);
    UnsignedInt GetNumberOfParticlesWithChosenEntityId(UniqueIdInt EntityId);
protected:
    explicit CellEngineBasicParticlesOperations(std::unordered_map<UniqueIdInt, Particle>& ParticlesParam) : Particles(ParticlesParam)
    {
    }
public:
    virtual ~CellEngineBasicParticlesOperations() = default;
};

#endif
