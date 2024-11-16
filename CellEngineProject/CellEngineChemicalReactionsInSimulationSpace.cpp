
#include <unordered_set>

#include "CellEngineConstants.h"
#include "CellEngineParticle.h"
#include "CellEngineAminoAcids.h"

#include "CellEngineExecutionTimeStatistics.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

using namespace std;

void CellEngineChemicalReactionsInSimulationSpace::RemoveParticle(const UniqueIdInt ParticleIndex, const bool ClearVoxels)
{
    try
    {
        Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

        CutDNAPrev(&ParticleObject);
        CutDNANext(&ParticleObject);
        SeparateTwoPairedDNANucleotides(&ParticleObject);
        DeleteLinkedParticlesPointersList(ParticleObject);
        ClearSpaceForParticle(ParticleObject, ClearVoxels);

        lock_guard LockGuardObject{ MainParticlesIndexesMutexObject };
        FreeIndexesOfParticles.push(ParticleIndex);

        if (CurrentThreadIndex == 0)
            Particles.erase(ParticleIndex);
        else
            ParticlesForThreads.erase(ParticleIndex);
    }
    CATCH("removing particle")
}

void CellEngineChemicalReactionsInSimulationSpace::MakingZeroSizeForContainersForFoundParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesKindsFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesSortedByCapacityFoundInProximity.clear();

        GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesWithFreeNextEndingsFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesWithFreePrevEndingsFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesWithFreePrevEndingsFoundInProximity.clear();

        GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesFreeFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFreeFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFoundInProximity.clear();

        GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesFullFreeFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFullFreeFoundInProximity.clear();

        GetThreadsLocalParticlesInProximity(CurrentThreadPos).tRNAUnchargedFoundInProximity.clear();
        GetThreadsLocalParticlesInProximity(CurrentThreadPos).tRNAChargedFoundInProximity.clear();
    }
    CATCH("making zero size for containers for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::UpdateFoundNucleotidesForFoundParticlesInProximity(const UnsignedInt ParticleIndex, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineConfigDataObject.RNAInOneParticle == false)
        {
            if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
            {
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            }
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }
        else
        {
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId) && ParticleObject.SequenceStr.substr(0, RNAStartSequence.length()) == RNAStartSequence)
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }

        if (CellEngineUseful::IsFreeDNANucleotide(ParticleObject.EntityId) == true)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).DNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsFreeRNANucleotide(ParticleObject.EntityId) == true)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).RNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineAminoAcidsManagerObject.IstRNAUncharged(ParticleObject.EntityId) == true)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).tRNAUnchargedFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineAminoAcidsManagerObject.IstRNACharged(ParticleObject.EntityId) == true)
            GetThreadsLocalParticlesInProximity(CurrentThreadPos).tRNAChargedFoundInProximity.emplace_back(ParticleIndex);

        GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesKindsFoundInProximity[ParticleObject.EntityId]++;
    }
    CATCH("updating found nucleotides for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::SaveParticleFoundInProximity(const UniqueIdInt ParticleIndex, unordered_set<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        if (FoundParticleIndexes.contains(ParticleIndex) == false)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleIndex);

            Particle& ParticleFromIndex = GetParticleFromIndex(ParticleIndex);

            if (UpdateNucleotides == true)
                UpdateFoundNucleotidesForFoundParticlesInProximity(ParticleIndex, CurrentThreadPos);

            GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesKindsFoundInProximity[ParticleFromIndex.EntityId]++;

            FoundParticleIndexes.insert(ParticleIndex);

            const auto stop_time = chrono::high_resolution_clock::now();

            CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForSavingFoundParticles += chrono::duration(stop_time - start_time);
        }
    }
    CATCH("saving particle found in proximity")
}

bool CellEngineChemicalReactionsInSimulationSpace::FindParticlesInProximityOfSimulationSpaceForSelectedSpace(const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        unordered_set<UnsignedInt> FoundParticleIndexes;

        MakingZeroSizeForContainersForFoundParticlesInProximity(CurrentThreadPos);

        FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(FoundParticleIndexes, UpdateNucleotides, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, CurrentThreadPos);

        LoggersManagerObject.Log(STREAM(endl << "Number of Particles Kinds Found In Proximity = " << GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesKindsFoundInProximity.size()));

        if (GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesSortedByCapacityFoundInProximity.empty() == false)
        {
            //sort(ParticlesSortedByCapacityFoundInProximity.begin(), ParticlesSortedByCapacityFoundInProximity.end(), [this](const UnsignedInt PK1, const UnsignedInt PK2) { return GetParticleFromIndex(PK1).ListOfVoxels.size() > GetParticleFromIndex(PK2).ListOfVoxels.size(); });
            PrintInformationAboutFoundParticlesInProximity(CurrentThreadPos);
        }
        else
        {
            LoggersManagerObject.Log(STREAM(endl << "No particle found in proximity"));
            return false;
        }

        LoggersManagerObject.Log(STREAM("Looking for particles in proximity done"));
    }
    CATCH("finding particles in proximity of simulation space for selected space")

    return true;
}

void CellEngineChemicalReactionsInSimulationSpace::PrintInformationAboutFoundParticlesInProximity(const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
        for (const auto& LocalParticleIndexObjectToWrite : GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesSortedByCapacityFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).ChainId) : '0') << " GENOME INDEX = " << GetParticleFromIndex(LocalParticleIndexObjectToWrite).GenomeIndex));
        LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInProximity List"));
        for (const auto& LocalParticleKindObjectToWrite : GetThreadsLocalParticlesInProximity(CurrentThreadPos).ParticlesKindsFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleKind EntityId = " << to_string(LocalParticleKindObjectToWrite.first) << " in quantity = " << to_string(LocalParticleKindObjectToWrite.second)));
    }
    CATCH("printing information found particles in proximity")
}

bool CellEngineChemicalReactionsInSimulationSpace::FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, const UnsignedInt AdditionalBoundFactor, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId)));

        auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalBoundFactor, CurrentThreadPos);
    }
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

