
#ifndef CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H
#define CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H

#include "ExceptionsMacro.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineBasicParticlesOperations.h"

class ReactionStatistics
{
public:
    UnsignedInt ReactionId{};
public:
    UnsignedInt Counter = 0;
};

class ParticleKindStatistics
{
public:
    EntityIdInt EntityId{};
public:
    UnsignedInt Counter = 0;
};

class CellEngineSimulationSpaceStatistics : virtual public CellEngineBasicParticlesOperations
{
protected:
    UnsignedInt SimulationStepNumber = 0;
    UnsignedInt ModuloStepToSaveStatistics = 10;
protected:
    static inline std::mutex MainStatisticsMutexObject;
    static inline std::shared_mutex MainStatisticsSharedMutexObject;
protected:
    bool SaveReactionsStatisticsBool = true;
    bool SaveParticlesAsCopiedMapBool = true;
    bool SaveParticlesAsVectorElementsBool = true;
    bool SortParticlesAsSortedVectorElementsBool = true;
protected:
    std::vector<std::vector<Particle>> ParticlesSnapshots;
    std::vector<std::unordered_map<UniqueIdInt, Particle>> ParticlesSnapshotsCopiedUnorderedMap;
    std::vector<std::vector<ParticleKindStatistics>> ParticlesKindsSnapshotsVectorSortedByCounter;
protected:
    std::vector<std::map<EntityIdInt, ParticleKindStatistics>> ParticlesKindsSnapshotsCopiedMap;
protected:
    std::vector<std::map<UnsignedInt, ReactionStatistics>> SavedReactionsMap;
protected:
    void MakeSimulationStepNumberZeroForStatistics();
    void IncSimulationStepNumberForStatistics();
    void GenerateNewEmptyElementsForContainersForStatistics();
    void SaveParticlesStatistics();
    void CheckConditionsToIncSimulationStepNumberForStatistics();
protected:
    void SaveParticlesAsCopiedMap();
    void SaveParticlesAsVectorElements();
    void SaveParticlesAsSortedVectorElements();
    void SaveReactionForStatistics(const ChemicalReaction& ReactionParam);
protected:
    void GetNumberOfParticlesFromParticleKind(EntityIdInt ParticleKindId);
public:
    CellEngineSimulationSpaceStatistics() noexcept
    {
    }
    virtual ~CellEngineSimulationSpaceStatistics() = default;
};

#endif
