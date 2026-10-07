
#include <unordered_set>

#include "CellEngineDataFile.h"
#include "CellEngineMacros.h"
#include "CellEngineConstants.h"
#include "CellEngineParticle.h"
#include "CellEngineAminoAcids.h"

#include "CellEngineExecutionTimeStatistics.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

using namespace std;

void CellEngineChemicalReactionsInSimulationSpace::RemoveParticle(const UniqueIdUnsignedInt ParticleIndex, const bool ClearVoxels, const bool CountAllParticlesBool)
{
    try
    {
        #ifdef SIMULATION_DETAILED_DEBUG_LOG
        const bool NotFoundParticle2 = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].Particles.contains(ParticleIndex);
        #endif

        Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

        #ifdef SIMULATION_DETAILED_DEBUG_LOG
        const bool NotFoundParticle1 = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].Particles.contains(ParticleIndex);
        UnsignedInt NumberOfAllParticlesBeforeRemoval;
        if (CountAllParticlesBool == true)
            NumberOfAllParticlesBeforeRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        else
            NumberOfAllParticlesBeforeRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(CurrentSectorPos);
        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("NUMBER OF ALL PARTICLES BEFORE REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " Particle to be removed has index = " << ParticleObject.Index << " " << ParticleIndex << " NOT FOUND PARTICLE = " << NotFoundParticle1 << " " << NotFoundParticle2));
        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("NUMBER OF ALL PARTICLES BEFORE REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " Particle to be removed has index = " << ParticleObject.Index << " " << ParticleIndex << " NOT FOUND PARTICLE = " << NotFoundParticle1 << " " << NotFoundParticle2));
        #endif

        CutDNAPrev(&ParticleObject);
        CutDNANext(&ParticleObject);
        SeparateTwoPairedDNANucleotides(&ParticleObject);
        DeleteLinkedParticlesPointersList(ParticleObject);
        ClearSpaceForParticle(ParticleObject, ClearVoxels);

        const UnsignedInt NumberOfRemovedParticelAfterErase = GetParticles().erase(ParticleIndex);

        #ifdef SIMULATION_DETAILED_DEBUG_LOG
        UnsignedInt NumberOfAllParticlesAfterRemoval;
        if (CountAllParticlesBool == true)
            NumberOfAllParticlesAfterRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        else
            NumberOfAllParticlesAfterRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(CurrentSectorPos);
        const UnsignedInt DifferenceOfAllParticlesBeforeAndAfterRemoval = NumberOfAllParticlesBeforeRemoval - NumberOfAllParticlesAfterRemoval;
        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("NUMBER OF ALL PARTICLES AFTER REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " " << NumberOfAllParticlesAfterRemoval << " DIFFERENCE = " << DifferenceOfAllParticlesBeforeAndAfterRemoval << " " << NumberOfRemovedParticelAfterErase << " Particle to be removed has index = " << ParticleObject.Index << " " << ParticleIndex));
        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("NUMBER OF ALL PARTICLES AFTER REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " " << NumberOfAllParticlesAfterRemoval << " DIFFERENCE = " << DifferenceOfAllParticlesBeforeAndAfterRemoval << " " << NumberOfRemovedParticelAfterErase <<  " Particle to be removed has index = " << ParticleObject.Index << " " << ParticleIndex));
        #endif
    }
    CATCH("removing particle")
}

void CellEngineChemicalReactionsInSimulationSpace::MakingZeroSizeForContainersForFoundParticlesInProximity(const ThreadPosType& CurrentThreadPos)
{
    try
    {
        LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.clear();

        LocalThreadParticlesInProximityObject.NucleotidesWithFreeNextEndingsFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.NucleotidesWithFreePrevEndingsFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.DNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.DNANucleotidesWithFreePrevEndingsFoundInProximity.clear();

        LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.clear();

        LocalThreadParticlesInProximityObject.DNANucleotidesFullFreeFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.RNANucleotidesFullFreeFoundInProximity.clear();

        LocalThreadParticlesInProximityObject.tRNAUnchargedFoundInProximity.clear();
        LocalThreadParticlesInProximityObject.tRNAChargedFoundInProximity.clear();
    }
    CATCH("making zero size for containers for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::UpdateFoundNucleotidesForFoundParticlesInProximity(const UnsignedInt ParticleIndex)
{
    try
    {
        const Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            LocalThreadParticlesInProximityObject.NucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            LocalThreadParticlesInProximityObject.NucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            LocalThreadParticlesInProximityObject.DNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            LocalThreadParticlesInProximityObject.DNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineConfigDataObject.RNAInOneParticle == false)
        {
            if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
            {
                LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
                LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            }
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }
        else
        {
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId) && ParticleObject.SequenceStr.substr(0, RNAStartSequence.length()) == RNAStartSequence)
                LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
                LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }

        if (CellEngineUseful::IsFreeDNANucleotide(ParticleObject.EntityId) == true)
            LocalThreadParticlesInProximityObject.DNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsFreeRNANucleotide(ParticleObject.EntityId) == true)
            LocalThreadParticlesInProximityObject.RNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineAminoAcidsManagerObject.IstRNAUncharged(ParticleObject.EntityId) == true)
            LocalThreadParticlesInProximityObject.tRNAUnchargedFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineAminoAcidsManagerObject.IstRNACharged(ParticleObject.EntityId) == true)
            LocalThreadParticlesInProximityObject.tRNAChargedFoundInProximity.emplace_back(ParticleIndex);

        LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ParticleObject.EntityId]++;
    }
    CATCH("updating found nucleotides for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::SaveParticleFoundInProximity(const UniqueIdUnsignedInt ParticleIndex, MainSetType<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides)
{
    try
    {
        if (FoundParticleIndexes.contains(ParticleIndex) == false)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleIndex);

            const Particle& ParticleFromIndex = GetParticleFromIndex(ParticleIndex);

            if (UpdateNucleotides == true)
                UpdateFoundNucleotidesForFoundParticlesInProximity(ParticleIndex);

            LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ParticleFromIndex.EntityId]++;

            FoundParticleIndexes.insert(ParticleIndex);

            const auto stop_time = chrono::high_resolution_clock::now();

            CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForSavingFoundParticles += chrono::duration(stop_time - start_time);
        }
    }
    CATCH("saving particle found in proximity")
}

bool CellEngineChemicalReactionsInSimulationSpace::FindParticlesInProximityOfSimulationSpaceForSelectedSpace(const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        MainSetType<UnsignedInt> FoundParticleIndexes;

        MakingZeroSizeForContainersForFoundParticlesInProximity(CurrentThreadPos);

        FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(FoundParticleIndexes, UpdateNucleotides, StartXPosParam, StartYPosParam, StartZPosParam, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType);

        DEBUGLOG(LoggersManagerObject.Log(STREAM(endl << "Number of Particles Kinds Found In Proximity = " << LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size()));)

        if (LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.empty() == false)
        {
            if constexpr (AdditionalSortParticlesInProximityByCapacity == true)
                ranges::sort(LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity, [this](const UnsignedInt PK1, const UnsignedInt PK2) { return GetParticleFromIndex(PK1).ListOfVoxels.size() > GetParticleFromIndex(PK2).ListOfVoxels.size(); });

            DEBUGLOG(PrintInformationAboutFoundParticlesInProximity(CurrentThreadPos);)
        }
        else
        {
            DEBUGLOG(LoggersManagerObject.Log(STREAM(endl << "No particle found in proximity"));)
            return false;
        }

        DEBUGLOG(LoggersManagerObject.Log(STREAM("Looking for particles in proximity done"));)
    }
    CATCH("finding particles in proximity of simulation space for selected space")

    return true;
}

void CellEngineChemicalReactionsInSimulationSpace::PrintInformationAboutFoundParticlesInProximity(const ThreadPosType& CurrentThreadPos)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
        for (const auto& LocalParticleIndexObjectToWrite : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).ChainId) : '0') << " GENOME INDEX = " << GetParticleFromIndex(LocalParticleIndexObjectToWrite).GenomeIndex));
        LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInProximity List"));
        for (const auto& [LocalParticleKindObjectToWriteFirst, LocalParticleKindObjectToWriteSecond] : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleKind EntityId = " << to_string(LocalParticleKindObjectToWriteFirst) << " in quantity = " << to_string(LocalParticleKindObjectToWriteSecond)));
    }
    CATCH("printing information found particles in proximity")
}

