
#ifndef CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H
#define CELL_ENGINE_BASIC_PARTICLES_OPERATIONS_H

#include <stack>
#include <shared_mutex>

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineParticleKind.h"
#include "CellEngineParticlesVoxelsOperations.h"

class CellEngineBasicParticlesOperations
{
protected:
    static inline std::shared_mutex MainParticlesSharedMutexObject;
    static inline std::mutex MainParticlesMutexObject;
    static inline std::mutex MainParticlesIndexesMutexObject;
protected:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
protected:
    UnsignedInt MaxParticleIndex{};
    std::stack<UniqueIdInt> FreeIndexesOfParticles;
protected:
    std::unordered_map<UniqueIdInt, Particle>& Particles;
protected:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        std::shared_lock<std::shared_mutex> LockGuardObject{ MainParticlesSharedMutexObject };

        return Particles[ParticleIndex];
    }
public:
    [[nodiscard]] UniqueIdInt GetFreeIndexesOfParticleSize() const
    {
        return FreeIndexesOfParticles.size();
    }
protected:
    void InitiateFreeParticleIndexes();
protected:
    inline UniqueIdInt GetNewFreeIndexOfParticle()
    {
        std::lock_guard<std::mutex> LockGuardObject{ MainParticlesIndexesMutexObject };

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
        std::lock_guard<std::mutex> LockGuardObject{ MainParticlesMutexObject };

        Particles[ParticleParam.Index] = ParticleParam;
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
    static void UpdateParticleKindListOfVoxels(Particle& ParticleObject, UnsignedInt ParticleXMin, UnsignedInt ParticleXMax, UnsignedInt ParticleYMin, UnsignedInt ParticleYMax, UnsignedInt ParticleZMin, UnsignedInt ParticleZMax);
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
