
#include <shared_mutex>

#include "ExceptionsMacro.h"
#include "CellEngineDataFile.h"
#include "CellEngineSimulationSpaceStatistics.h"
#include "CellEngineParticlesKindsManager.h"

using namespace std;

CellEngineSimulationSpaceStatistics::CellEngineSimulationSpaceStatistics()
{
}

void CellEngineSimulationSpaceStatistics::MakeSimulationStepNumberZeroForStatistics()
{
    try
    {
        SimulationStepNumber = 0;

        ParticlesSnapshots.clear();
        ParticlesSnapshotsCopiedUnorderedMap.clear();
        ParticlesKindsSnapshotsVectorSortedByCounter.clear();
        ParticlesKindsSnapshotsCopiedMap.clear();
        ParticlesKindsHistogramComparisons.clear();

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
        ParticlesKindsHistogramComparisons.emplace_back();

        SavedReactionsMap.emplace_back();

        LoggersManagerObject.LogStatistics(STREAM("Size of arrays to statistics = " << ParticlesSnapshots.size() << " " << ParticlesKindsSnapshotsVectorSortedByCounter.size() << " " << ParticlesKindsSnapshotsCopiedMap.size() << " " << ParticlesKindsHistogramComparisons.size() << " " << SavedReactionsMap.size()));
    }
    CATCH("incrementing simulation step number for statistics")
}

void CellEngineSimulationSpaceStatistics::SaveParticlesStatistics()
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("SAVE PARTICLES ARRAYS"));

        if (SaveParticlesAsCopiedMapBool == true)
            SaveParticlesAsCopiedMap();
        if (SaveParticlesAsVectorElementsBool == true)
            SaveParticlesAsVectorElements();
        if (SortParticlesAsSortedVectorElementsBool == true)
            SaveParticlesAsSortedVectorElements();
        if (SortHistogramOfParticlesAsSortedVectorElementsBool == true)
            CompareHistogramsOfParticles(SimulationStepNumber - 1, SimulationStepNumber);
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

        #ifdef SHORTER_CODE
        for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber - 1])
        {
            UniqueIdInt ParticleKindId = GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId;
            ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][ParticleKindId] = { ParticleKindId, ++ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][ParticleKindId].Counter };
        }
        #else
        for (const auto& ParticlesSnapshotsCopiedUnorderedMapElement : ParticlesSnapshotsCopiedUnorderedMap[SimulationStepNumber - 1])
        {
            UniqueIdInt ParticleKindId = GetParticleFromIndex(ParticlesSnapshotsCopiedUnorderedMapElement.first).EntityId;
            if (auto FoundResult = ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].find(ParticleKindId); FoundResult != ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].end())
                FoundResult->second.Counter++;
            else
                ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1][ParticleKindId] = { ParticleKindId, 1 };
        }
        #endif

        LoggersManagerObject.LogStatistics(STREAM("Size of all particle kinds with number of particles for particle kind = " << ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].size()));

        transform(ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsCopiedMap[SimulationStepNumber - 1].end(), back_inserter(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1]), [](const auto& ParticlesMapElement){ return ParticlesMapElement.second; } );

        sort(ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].begin(), ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber - 1].end(), [](const auto& P1, const auto& P2){ return P1.Counter > P2.Counter; } );
    }
    CATCH("saving particles as sorted vector elements")
}

void CellEngineSimulationSpaceStatistics::CompareHistogramsOfParticles(const UnsignedInt SimulationStepNumber1, const UnsignedInt SimulationStepNumber2)
{
    try
    {
        LoggersManagerObject.LogStatistics(STREAM("Size of particles vector sorted by number = " << ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber1 - 1].size() << " " << ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber2 - 1].size() << " " << SimulationStepNumber1 << " " << SimulationStepNumber2));

        for (const auto& ParticleKindStatisticsObject1 : ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber1 - 1])
        {
            bool Found = false;
            for (const auto& ParticleKindStatisticsObject2 : ParticlesKindsSnapshotsVectorSortedByCounter[SimulationStepNumber2 - 1])
                if (ParticleKindStatisticsObject2.EntityId == ParticleKindStatisticsObject1.EntityId)
                {
                    ParticlesKindsHistogramComparisons[SimulationStepNumber - 1].emplace_back(ParticleKindHistogramComparison{ ParticleKindStatisticsObject1.EntityId, ParticleKindStatisticsObject1.Counter, ParticleKindStatisticsObject2.Counter, static_cast<SignedInt>(ParticleKindStatisticsObject2.Counter - ParticleKindStatisticsObject1.Counter) });
                    Found = true;
                    break;
                }
            if (Found == false)
                ParticlesKindsHistogramComparisons[SimulationStepNumber - 1].emplace_back(ParticleKindHistogramComparison{ ParticleKindStatisticsObject1.EntityId, ParticleKindStatisticsObject1.Counter, 0, static_cast<SignedInt>(0 - ParticleKindStatisticsObject1.Counter) });
        }
    }
    CATCH("comparing histograms of particles")
}

#ifdef SHORTER_CODE
void CellEngineSimulationSpaceStatistics::SaveReactionForStatistics(const ChemicalReaction& ReactionParam)
{
    try
    {
        SavedReactionsMap[SimulationStepNumber - 1][ReactionParam.ReactionIdNum].Counter++;
    }
    CATCH("saving reaction for statistics")
}
#else
void CellEngineSimulationSpaceStatistics::SaveReactionForStatisticsExtended(const ChemicalReaction& ReactionParam)
{
    try
    {
        auto ReactionIter = SavedReactionsMap[SimulationStepNumber - 1].find(ReactionParam.ReactionIdNum);
        if (ReactionIter != SavedReactionsMap[SimulationStepNumber - 1].end())
            ReactionIter->second.Counter++;
        else
            SavedReactionsMap[SimulationStepNumber - 1][ReactionParam.ReactionIdNum] = ReactionStatistics{ ReactionParam.ReactionIdNum, 1 };
    }
    CATCH("saving reaction for statistics extended")
}
#endif

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