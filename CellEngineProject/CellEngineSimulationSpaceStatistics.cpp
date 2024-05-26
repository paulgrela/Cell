
#include "CellEngineSimulationSpaceStatistics.h"

void CellEngineSimulationSpaceStatistics::SaveParticlesAsCopiedMad()
{
    ParticlesSnapshotsCopiedUnorderedMap.emplace_back(Particles);
}

void CellEngineSimulationSpaceStatistics::SaveParticlesAsVectorElements()
{
    ParticlesSnapshots[SimulationStepNumber].reserve(Particles.size());
    std::transform(Particles.begin(), Particles.end(), std::back_inserter(ParticlesSnapshots[SimulationStepNumber]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );
}

void CellEngineSimulationSpaceStatistics::SaveParticlesAsSortedVectorElements()
{
    ParticlesSnapshotsCopiedMap.clear();

    for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber])
        ParticlesSnapshotsCopiedMap[SimulationStepNumber][ParticlesSnapshotsCopiedUnorderedMapElement.first].Counter++;

    std::transform(ParticlesSnapshotsCopiedMap[SimulationStepNumber].begin(), ParticlesSnapshotsCopiedMap[SimulationStepNumber].end(), std::back_inserter(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );

    std::sort(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber].end(), [](const auto& P1, const auto& P2){ return P1.Counter > P2.Counter; } );
}

void CellEngineSimulationSpaceStatistics::SaveReactionForStatistics(const Reaction& ReactionParam)
{
    SavedReactionsMap[SimulationStepNumber][ReactionParam.Id].Counter++;
}