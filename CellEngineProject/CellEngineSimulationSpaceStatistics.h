
#ifndef CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H
#define CELL_ENGINE_SIMULATION_SPACE_STATISTICS_H

#include "ExceptionsMacro.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
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
protected:
    std::vector<std::vector<Particle>> ParticlesSnapshots;
    std::vector<std::unordered_map<UniqueIdInt, Particle>> ParticlesSnapshotsCopiedUnorderedMap;
    std::vector<std::vector<ParticleKindStatistics>> ParticlesKindsSnapshotsVectorSortedByCounter;
protected:
    std::vector<std::map<EntityIdInt, ParticleKindStatistics>> ParticlesSnapshotsCopiedMap;
protected:
    std::vector<std::map<UnsignedInt, ReactionStatistics>> SavedReactionsMap;
protected:
    void MakeSimulationStepNumberZeroForStatistics();
    void IncSimulationStepNumberForStatistics();
protected:
    void SaveParticlesAsCopiedMad();
    void SaveParticlesAsVectorElements();
    void SaveParticlesAsSortedVectorElements();
    void SaveReactionForStatistics(const Reaction& ReactionParam);
};

#endif
