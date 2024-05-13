
#include "CellEngineConstants.h"
#include "CellEngineParticle.h"

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
        FreeIndexesOfParticles.push(ParticleIndex);
        Particles.erase(ParticleIndex);
    }
    CATCH("removing particle")
}

void CellEngineChemicalReactionsInSimulationSpace::MakingZeroSizeForContainersForFoundParticlesInProximity()
{
    try
    {
        ParticlesKindsFoundInProximity.clear();
        ParticlesSortedByCapacityFoundInProximity.clear();
        NucleotidesWithFreeNextEndingsFoundInProximity.clear();
        NucleotidesWithFreePrevEndingsFoundInProximity.clear();
        DNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        DNANucleotidesWithFreePrevEndingsFoundInProximity.clear();
        NucleotidesFreeFoundInProximity.clear();
        RNANucleotidesFreeFoundInProximity.clear();
        RNANucleotidesFoundInProximity.clear();
    }
    CATCH("making zero size for containers for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::UpdateFoundNucleotidesForFoundParticlesInProximity(const UnsignedInt ParticleIndex)
{
    try
    {
        Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);

        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            NucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            NucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            DNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            DNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineConfigDataObject.RNAInOneParticle == false)
        {
            if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
            {
                NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
                RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            }
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }
        else
        {
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId) && ParticleObject.SequenceStr.length() == 1)
                RNANucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsDNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
                NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
            if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
                RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);
        }

        ParticlesKindsFoundInProximity[ParticleObject.EntityId]++;
    }
    CATCH("updating found nucleotides for found particles in proximity")
}

void CellEngineChemicalReactionsInSimulationSpace::SaveParticleFoundInProximity(const UniqueIdInt ParticleIndex, set<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides)
{
    try
    {
        if (FoundParticleIndexes.find(ParticleIndex) == FoundParticleIndexes.end())
        {
            ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleIndex);

            Particle& ParticleFromIndex = GetParticleFromIndex(ParticleIndex);

            if (UpdateNucleotides == true)
                UpdateFoundNucleotidesForFoundParticlesInProximity(ParticleIndex);

            ParticlesKindsFoundInProximity[ParticleFromIndex.EntityId]++;

            FoundParticleIndexes.insert(ParticleIndex);
        }
    }
    CATCH("saving particle found in proximity")
}

bool CellEngineChemicalReactionsInSimulationSpace::FindParticlesInProximityOfSimulationSpaceForSelectedSpace(const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        set<UnsignedInt> FoundParticleIndexes;

        MakingZeroSizeForContainersForFoundParticlesInProximity();

        FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(FoundParticleIndexes, UpdateNucleotides, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam);

        if (ParticlesSortedByCapacityFoundInProximity.empty() == false)
        {
            sort(ParticlesSortedByCapacityFoundInProximity.begin(), ParticlesSortedByCapacityFoundInProximity.end(), [this](UnsignedInt PK1, UnsignedInt PK2) { return GetParticleFromIndex(PK1).ListOfVoxels.size() > GetParticleFromIndex(PK2).ListOfVoxels.size(); });
            PrintInformationAboutFoundParticlesInProximity();
        }
        else
        {
            LoggersManagerObject.Log(STREAM(endl << "No particle found in proximity "));
            return false;
        }
    }
    CATCH("finding particles in proximity of simulation space for selected space")

    return true;
}

void CellEngineChemicalReactionsInSimulationSpace::PrintInformationAboutFoundParticlesInProximity()
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
        for (const auto& LocalParticleIndexObjectToWrite : ParticlesSortedByCapacityFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).ChainId) : '0') << " GENOME INDEX = " << GetParticleFromIndex(LocalParticleIndexObjectToWrite).GenomeIndex));
        LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInProximity List"));
        for (const auto& LocalParticleKindObjectToWrite : ParticlesKindsFoundInProximity)
            LoggersManagerObject.Log(STREAM("ParticleKind EntityId = " << to_string(LocalParticleKindObjectToWrite.first) << " in quantity = " << to_string(LocalParticleKindObjectToWrite.second)));
    }
    CATCH("printing information found particles in proximity")
}

bool CellEngineChemicalReactionsInSimulationSpace::FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, const UnsignedInt AdditionalBoundFactor)
{
    try
    {
        LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId)));

        auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

        FindParticlesInProximityOfSimulationSpaceForSelectedSpace(true, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalBoundFactor);
    }
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

