
#include "ExceptionsMacro.h"

#include "CellEngineSimulationSpaceStatistics.h"

void CellEngineSimulationSpaceStatistics::MakeSimulationStepNumberZeroForStatistics()
{
    try
    {
        SimulationStepNumber = 0;

        ParticlesSnapshots.clear();
        ParticlesSnapshotsCopiedUnorderedMap.clear();
        ParticlesKindsSnapshotsVectorSortedByCounter.clear();
        ParticlesSnapshotsCopiedMap.clear();
        SavedReactionsMap.clear();
    }
    CATCH("making simulation step number zero for statistics")
}

void CellEngineSimulationSpaceStatistics::IncSimulationStepNumberForStatistics()
{
    try
    {
        SimulationStepNumber++;
    }
    CATCH("incrementing simulation step number for statistics")
}

void CellEngineSimulationSpaceStatistics::GenerateNewEmptyElementsForContainersForStatistics()
{
    try
    {
        ParticlesSnapshots.emplace_back();
        ParticlesSnapshotsCopiedUnorderedMap.emplace_back();
        ParticlesKindsSnapshotsVectorSortedByCounter.emplace_back();
        ParticlesSnapshotsCopiedMap.emplace_back();

        SavedReactionsMap.emplace_back();
    }
    CATCH("incrementing simulation step number for statistics")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesStatistics()
{
    try
    {
        if (SaveParticlesAsCopiedMapBool == true)
            SaveParticlesAsCopiedMad();
        if (SaveParticlesAsVectorElementsBool == true)
            SaveParticlesAsVectorElements();
        if (SortParticlesAsSortedVectorElementsBool == true)
            SaveParticlesAsSortedVectorElements();
    }
    CATCH("saving particles statistics")
}

void CellEngineSimulationSpaceStatistics::CheckConditionsToIncSimulationStepNumberForStatistics()
{
    try
    {
        if (SimulationStepNumber % ModuloStepToSaveStatistics == 0)
        {
            GenerateNewEmptyElementsForContainersForStatistics();
            SaveParticlesStatistics();
        }
    }
    CATCH("checking conditions to increment simulation step number for statistics")
}

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
        ParticlesSnapshots[SimulationStepNumber - 1].reserve(Particles.size());

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
            ParticlesSnapshotsCopiedMap[SimulationStepNumber - 1][ParticlesSnapshotsCopiedUnorderedMapElement.first].Counter++;

        transform(ParticlesSnapshotsCopiedMap[SimulationStepNumber - 1].begin(), ParticlesSnapshotsCopiedMap[SimulationStepNumber - 1].end(), back_inserter(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );

        sort(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].end(), [](const auto& P1, const auto& P2){ return P1.Counter > P2.Counter; } );
    }
    CATCH("saving particles as sorted vector elements")
}

void CellEngineSimulationSpaceStatistics::SaveReactionForStatistics(const Reaction& ReactionParam)
{
    try
    {
        auto ReactionIter = SavedReactionsMap[SimulationStepNumber - 1].find(ReactionParam.ReactionIdNum);

        if (ReactionIter != SavedReactionsMap[SimulationStepNumber - 1].end())
            ReactionIter->second.Counter++;
        else
            SavedReactionsMap[SimulationStepNumber - 1][ReactionParam.ReactionIdNum] = ReactionStatistics{ ReactionParam.ReactionIdNum, 0 };
    }
    CATCH("saving reaction for statistics")
}