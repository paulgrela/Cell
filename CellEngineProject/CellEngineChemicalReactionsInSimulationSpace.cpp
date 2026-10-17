
#include <unordered_set>

#include "CellEngineMacros.h"
#include "CellEngineDataFile.h"
#include "CellEngineConstants.h"
#include "CellEngineParticle.h"
#include "CellEngineAminoAcids.h"

#include "CellEngineSimulationExecutionMeasureTimeStatisticsData.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsInSimulationSpace.h"

using namespace std;

void CellEngineChemicalReactionsInSimulationSpace::RemoveParticle(const UniqueIdUInt64Type ParticleUniqueId, const bool ClearVoxels, const bool CountAllParticlesBool)
{
    try
    {
        #if defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION) || defined(SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION)
        const bool NotFoundParticle2 = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].Particles.contains(ParticleUniqueId);
        #endif

        Particle& ParticleObject = GetParticleFromUniqueId(ParticleUniqueId);

        #if defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION) || defined(SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION)
        const bool NotFoundParticle1 = Particles[CurrentSectorPos.SectorPosX][CurrentSectorPos.SectorPosY][CurrentSectorPos.SectorPosZ].Particles.contains(ParticleUniqueId);
        UnsignedInt NumberOfAllParticlesBeforeRemoval;
        if (CountAllParticlesBool == true)
            NumberOfAllParticlesBeforeRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        else
            NumberOfAllParticlesBeforeRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(CurrentSectorPos);
        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("NUMBER OF ALL PARTICLES BEFORE REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " Particle to be removed has index = " << ParticleObject.ParticleUniqueId << " " << ParticleUniqueId << " NOT FOUND PARTICLE = " << NotFoundParticle1 << " " << NotFoundParticle2));
        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("NUMBER OF ALL PARTICLES BEFORE REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " Particle to be removed has index = " << ParticleObject.ParticleUniqueId << " " << ParticleUniqueId << " NOT FOUND PARTICLE = " << NotFoundParticle1 << " " << NotFoundParticle2));
        #endif

        CutDNAPrev(&ParticleObject);
        CutDNANext(&ParticleObject);
        SeparateTwoPairedDNANucleotides(&ParticleObject);
        DeleteLinkedParticlesPointersList(ParticleObject);
        ClearSpaceForParticle(ParticleObject, ClearVoxels);

        const UnsignedInt NumberOfRemovedParticelAfterErase = GetParticles().erase(ParticleUniqueId);

        #if defined(SIMULATION_DETAILED_DEBUG_LOG_DIFFUSION) || defined(SIMULATION_DETAILED_DEBUG_LOG_CHEMICAL_REACTION)
        UnsignedInt NumberOfAllParticlesAfterRemoval;
        if (CountAllParticlesBool == true)
            NumberOfAllParticlesAfterRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticles();
        else
            NumberOfAllParticlesAfterRemoval = CellEngineDataFileObjectPointer->GetNumberOfAllParticlesInSector(CurrentSectorPos);
        const UnsignedInt DifferenceOfAllParticlesBeforeAndAfterRemoval = NumberOfAllParticlesBeforeRemoval - NumberOfAllParticlesAfterRemoval;
        LoggersManagerObject.LogOnlyToFilesUnconditional(STREAM("NUMBER OF ALL PARTICLES AFTER REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " " << NumberOfAllParticlesAfterRemoval << " DIFFERENCE = " << DifferenceOfAllParticlesBeforeAndAfterRemoval << " " << NumberOfRemovedParticelAfterErase << " Particle to be removed has index = " << ParticleObject.ParticleUniqueId << " " << ParticleUniqueId));
        LoggersManagerObject.LogOnlyToConsoleUnconditional(STREAM("NUMBER OF ALL PARTICLES AFTER REMOVAL = " << NumberOfAllParticlesBeforeRemoval << " " << NumberOfAllParticlesAfterRemoval << " DIFFERENCE = " << DifferenceOfAllParticlesBeforeAndAfterRemoval << " " << NumberOfRemovedParticelAfterErase <<  " Particle to be removed has index = " << ParticleObject.ParticleUniqueId << " " << ParticleUniqueId));
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

void CellEngineChemicalReactionsInSimulationSpace::UpdateFoundNucleotidesForFoundParticlesInProximity(const UnsignedInt ParticleUniqueId)
{
    try
    {
        const Particle& ParticleObject = GetParticleFromUniqueId(ParticleUniqueId);

        if (CellEngineUseful::IsDNAorRNA(ParticleObject.ParticleKindId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            LocalThreadParticlesInProximityObject.NucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleUniqueId);
        if (CellEngineUseful::IsDNAorRNA(ParticleObject.ParticleKindId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            LocalThreadParticlesInProximityObject.NucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleUniqueId);

        if (CellEngineUseful::IsDNA(ParticleObject.ParticleKindId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            LocalThreadParticlesInProximityObject.DNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleUniqueId);
        if (CellEngineUseful::IsDNA(ParticleObject.ParticleKindId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            LocalThreadParticlesInProximityObject.DNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleUniqueId);

        if (CellEngineConfigDataObject.RNAInOneParticle == false)
        {
            if (CellEngineUseful::IsDNAorRNA(ParticleObject.ParticleKindId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
            {
                LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.emplace_back(ParticleUniqueId);
                LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.emplace_back(ParticleUniqueId);
            }
            if (CellEngineUseful::IsRNA(ParticleObject.ParticleKindId))
                LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.emplace_back(ParticleUniqueId);
        }
        else
        {
            if (CellEngineUseful::IsRNA(ParticleObject.ParticleKindId) && ParticleObject.SequenceStr.substr(0, RNAStartSequence.length()) == RNAStartSequence)
                LocalThreadParticlesInProximityObject.RNANucleotidesFreeFoundInProximity.emplace_back(ParticleUniqueId);
            if (CellEngineUseful::IsDNA(ParticleObject.ParticleKindId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
                LocalThreadParticlesInProximityObject.NucleotidesFreeFoundInProximity.emplace_back(ParticleUniqueId);
            if (CellEngineUseful::IsRNA(ParticleObject.ParticleKindId))
                LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.emplace_back(ParticleUniqueId);
        }

        if (CellEngineUseful::IsFreeDNANucleotide(ParticleObject.ParticleKindId) == true)
            LocalThreadParticlesInProximityObject.DNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleUniqueId);
        if (CellEngineUseful::IsFreeRNANucleotide(ParticleObject.ParticleKindId) == true)
            LocalThreadParticlesInProximityObject.RNANucleotidesFullFreeFoundInProximity.emplace_back(ParticleUniqueId);

        if (CellEngineAminoAcidsManagerObject.IstRNAUncharged(ParticleObject.ParticleKindId) == true)
            LocalThreadParticlesInProximityObject.tRNAUnchargedFoundInProximity.emplace_back(ParticleUniqueId);
        if (CellEngineAminoAcidsManagerObject.IstRNACharged(ParticleObject.ParticleKindId) == true)
            LocalThreadParticlesInProximityObject.tRNAChargedFoundInProximity.emplace_back(ParticleUniqueId);

        LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ParticleObject.ParticleKindId]++;
    }
    CATCH("updating found nucleotides for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::SaveParticleFoundInProximity(const UniqueIdUInt64Type ParticleUniqueId, MainSetType<UnsignedInt>& FoundParticleUniqueIdes, const bool UpdateNucleotides)
{
    try
    {
        if (FoundParticleUniqueIdes.contains(ParticleUniqueId) == false)
        {
            const auto start_time = chrono::high_resolution_clock::now();

            LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleUniqueId);

            const Particle& ParticleFromIndex = GetParticleFromUniqueId(ParticleUniqueId);

            if (UpdateNucleotides == true)
                UpdateFoundNucleotidesForFoundParticlesInProximity(ParticleUniqueId);

            LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity[ParticleFromIndex.ParticleKindId]++;

            FoundParticleUniqueIdes.insert(ParticleUniqueId);

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
        MainSetType<UnsignedInt> FoundParticleUniqueIdes;

        MakingZeroSizeForContainersForFoundParticlesInProximity(CurrentThreadPos);

        FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(FoundParticleUniqueIdes, UpdateNucleotides, StartXPosParam, StartYPosParam, StartZPosParam, UnusedValueRealType, UnusedValueRealType, UnusedValueRealType);

        DEBUG_LOG(LoggersManagerObject.Log(STREAM(endl << "Number of Particles Kinds Found In Proximity = " << LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity.size()));)

        if (LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity.empty() == false)
        {
            if constexpr (AdditionalSortParticlesInProximityByCapacity == true)
                ranges::sort(LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity, [this](const UnsignedInt PK1, const UnsignedInt PK2) { return GetParticleFromUniqueId(PK1).ListOfVoxels.size() > GetParticleFromUniqueId(PK2).ListOfVoxels.size(); });

            DEBUG_LOG(PrintInformationAboutFoundParticlesInProximity(CurrentThreadPos);)
        }
        else
        {
            DEBUG_LOG(LoggersManagerObject.Log(STREAM(endl << "No particle found in proximity"));)
            return false;
        }

        DEBUG_LOG(LoggersManagerObject.Log(STREAM("Looking for particles in proximity done"));)
    }
    CATCH("finding particles in proximity of simulation space for selected space")

    return true;
}

void CellEngineChemicalReactionsInSimulationSpace::PrintInformationAboutFoundParticlesInProximity(const ThreadPosType& CurrentThreadPos)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
        for (const auto& LocalParticleUniqueIdObjectToWrite : LocalThreadParticlesInProximityObject.ParticlesSortedByCapacityFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleUniqueId = " << to_string(LocalParticleUniqueIdObjectToWrite) << " ParticleKindId = " << to_string(GetParticleFromUniqueId(LocalParticleUniqueIdObjectToWrite).ParticleKindId) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromUniqueId(LocalParticleUniqueIdObjectToWrite).ParticleKindId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromUniqueId(LocalParticleUniqueIdObjectToWrite).ChainId) : '0') << " GENOME INDEX = " << GetParticleFromUniqueId(LocalParticleUniqueIdObjectToWrite).GenomeIndex));
        LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInProximity List"));
        for (const auto& [LocalParticleKindObjectToWriteFirst, LocalParticleKindObjectToWriteSecond] : LocalThreadParticlesInProximityObject.ParticlesKindsFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleKind ParticleKindId = " << to_string(LocalParticleKindObjectToWriteFirst) << " in quantity = " << to_string(LocalParticleKindObjectToWriteSecond)));
    }
    CATCH("printing information found particles in proximity")
}

