
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

class ParticleKindHistogramComparison
{
public:
    EntityIdInt EntityId{};
public:
    UnsignedInt Counter1 = 0;
    UnsignedInt Counter2 = 0;
    SignedInt Difference = 0;
};

class CellEngineSimulationSpaceStatistics : virtual public CellEngineBasicParticlesOperations
{
    friend class CellEngineImGuiMenu;
    friend class CellEngineSimulationParallelExecutionManager;
protected:
    UnsignedInt SimulationStepNumber = 0;
    UnsignedInt ModuloStepToSaveStatistics = 10;
protected:
    bool SaveReactionsStatisticsBool = true;
    bool SaveParticlesAsCopiedMapBool = true;
    bool SaveParticlesAsVectorElementsBool = true;
    bool SortParticlesAsSortedVectorElementsBool = true;
    bool SortHistogramOfParticlesAsSortedVectorElementsBool = true;
protected:
    std::vector<std::vector<Particle>> ParticlesSnapshots;
        std::vector<ParticlesDetailedContainer<Particle>> ParticlesSnapshotsCopiedUnorderedMap;
        std::vector<std::vector<EntityIdInt>> ParticlesSnapshotsCopiedVectorForMPI;
    std::vector<std::vector<ParticleKindStatistics>> ParticlesKindsSnapshotsVectorSortedByCounter;
    std::vector<std::vector<ParticleKindHistogramComparison>> ParticlesKindsHistogramComparisons;
protected:
    std::vector<std::map<EntityIdInt, ParticleKindStatistics>> ParticlesKindsSnapshotsCopiedMap;
protected:
        std::vector<std::map<UnsignedInt, ReactionStatistics>> SavedReactionsMap;
        std::vector<std::vector<ReactionStatistics>> SavedReactionsMapForMPI;
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
protected:
    void CompareHistogramsOfParticles(UnsignedInt SimulationStepNumber1, UnsignedInt SimulationStepNumber2);
protected:
    void SaveReactionForStatistics(const ChemicalReaction& ReactionParam);
protected:
    void GetNumberOfParticlesFromParticleKind(EntityIdInt ParticleKindId);
public:
    CellEngineSimulationSpaceStatistics();
    ~CellEngineSimulationSpaceStatistics() override = default;
};

#endif
