
#include "ExceptionsMacro.h"

#include "CellEngineSimulationSpaceStatistics.h"

void CellEngineSimulationSpaceStatistics::SaveParticlesAsCopiedMad()
{
    try
    {
        ParticlesSnapshotsCopiedUnorderedMap.emplace_back(Particles);
    }
    CATCH("saving particles as copied map")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesAsVectorElements()
{
    try
    {
        ParticlesSnapshots[SimulationStepNumber].reserve(Particles.size());

        transform(Particles.begin(), Particles.end(), std::back_inserter(ParticlesSnapshots[SimulationStepNumber]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );
    }
    CATCH("saving particles as vector elements")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesAsSortedVectorElements()
{
    try
    {
        ParticlesSnapshotsCopiedMap.clear();

        for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber])
            ParticlesSnapshotsCopiedMap[SimulationStepNumber][ParticlesSnapshotsCopiedUnorderedMapElement.first].Counter++;

        transform(ParticlesSnapshotsCopiedMap[SimulationStepNumber].begin(), ParticlesSnapshotsCopiedMap[SimulationStepNumber].end(), back_inserter(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );

        sort(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber].end(), [](const auto& P1, const auto& P2){ return P1.Counter > P2.Counter; } );
    }
    CATCH("saving particles as sorted vector elements")
}

void CellEngineSimulationSpaceStatistics::SaveReactionForStatistics(const Reaction& ReactionParam)
{
    try
    {
        SavedReactionsMap[SimulationStepNumber][ReactionParam.ReactionIdNum].Counter++;
    }
    CATCH("saving reaction for statistics")
}