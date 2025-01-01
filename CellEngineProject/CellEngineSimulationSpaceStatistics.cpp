
#include <shared_mutex>

#include "ExceptionsMacro.h"
#include "CellEngineDataFile.h"
#include "CellEngineSimulationSpaceStatistics.h"
#include "CellEngineParticlesKindsManager.h"

#define LONGER_CODE_

using namespace std;

void CellEngineSimulationSpaceStatistics::MakeSimulationStepNumberZeroForStatistics()
{
    try
    {
        SimulationStepNumber = 0;

        ParticlesSnapshots.clear();
        ParticlesSnapshotsCopiedUnorderedMap.clear();
        ParticlesKindsSnapshotsVectorSortedByCounter.clear();
        ParticlesKindsSnapshotsCopiedMap.clear();
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

        ParticlesKindsSnapshotsVectorSortedByCounter.emplace_back();
        ParticlesKindsSnapshotsCopiedMap.emplace_back();

        SavedReactionsMap.emplace_back();
    }
    CATCH("incrementing simulation step number for statistics")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesStatistics()
{
    try
    {
        if (SaveParticlesAsCopiedMapBool == true)
            SaveParticlesAsCopiedMap();
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

void CellEngineSimulationSpaceStatistics::SaveParticlesAsCopiedMap()
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

        transform(Particles.begin(), Particles.end(), std::back_inserter(ParticlesSnapshots[SimulationStepNumber - 1]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );
    }
    CATCH("saving particles as vector elements")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesAsSortedVectorElements()
{
    try
    {
        ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].clear();

        LoggersManagerObject.LogStatistics(STREAM("Size of all particles copied for statistics = " << ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber - 1].size()));

        #ifdef LONGER_CODE
        for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber - 1])
        {
            UniqueIdInt ParticleKind = GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId;
            ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][ParticleKind] = { ParticleKind, ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][ParticleKind].Counter++ };
        }
        #else
        for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber - 1])
            if (auto FoundResult = ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].find(GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId); FoundResult != ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].end())
                FoundResult->second.Counter++;
            else
                ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId] = { GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId, 0 };
        #endif

        LoggersManagerObject.LogStatistics(STREAM("Size of all particle kinds with number of particles for particle kind = " << ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].size()));

        transform(ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].end(), back_inserter(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );

        sort(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].end(), [](const auto& P1, const auto& P2){ return P1.Counter > P2.Counter; } );
    }
    CATCH("saving particles as sorted vector elements")
}

void CellEngineSimulationSpaceStatistics::SaveReactionForStatistics(const ChemicalReaction& ReactionParam)
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

void CellEngineSimulationSpaceStatistics::GetNumberOfParticlesFromParticleKind(const EntityIdInt ParticleKindId)
{
    try
    {
        auto FoundResult = find_if(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].end(), [ParticleKindId](const auto& P){ return P.EntityId == ParticleKindId; });
        if (FoundResult != ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].end())
            LoggersManagerObject.LogStatistics(STREAM("Particle Kind Name = " << ParticlesKindsManagerObject.GetParticleKind(FoundResult->EntityId).IdStr << " Particle Kind Id = " << FoundResult->EntityId << " Number of Particles = " << FoundResult->Counter));
    }
    CATCH("getting number of particles from particle kind")
}

// void CellEngineSimulationSpaceStatistics::JoinStatisticsFromThreads(std::map<UnsignedInt, ReactionStatistics>& SavedReactionsMap) const
// {
//     try
//     {
//         for (UnsignedInt ThreadXIndex = 1; ThreadXIndex <= CellEngineConfigDataObject.NumberOfXThreadsInSimulation; ThreadXIndex++)
//             for (UnsignedInt ThreadYIndex = 1; ThreadYIndex <= CellEngineConfigDataObject.NumberOfYThreadsInSimulation; ThreadYIndex++)
//                 for (UnsignedInt ThreadZIndex = 1; ThreadZIndex <= CellEngineConfigDataObject.NumberOfZThreadsInSimulation; ThreadZIndex++)
//                     for (const auto& ReactionData : CellEngineDataFileObjectPointer->CellEngineSimulationSpaceForThreadsObjectsPointer[ThreadXIndex - 1][ThreadYIndex - 1][ThreadZIndex - 1]->SavedReactionsMap[SimulationStepNumber - 1])
//                         if (SavedReactionsMap.contains(ReactionData.second.ReactionId))
//                             SavedReactionsMap.find(ReactionData.second.ReactionId)->second.Counter += ReactionData.second.Counter;
//                         else
//                             SavedReactionsMap.insert(ReactionData);
//     }
//     CATCH("joining statistics from threads")
// }
