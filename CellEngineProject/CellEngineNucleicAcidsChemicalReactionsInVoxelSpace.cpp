
#include <set>

#include "CellEngineParticle.h"

#include "CellEngineNucleicAcidsChemicalReactionsInVoxelSpace.h"

using namespace std;

void CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::RemoveParticle(const UniqueIdInt ParticleIndex, const bool ClearVoxels)
{
    try
    {
        Particle& ParticleObject = GetParticleFromIndex(ParticleIndex);
        CutDNAPrev(&ParticleObject);
        CutDNANext(&ParticleObject);
        SeparateTwoPairedDNANucleotides(&ParticleObject);
        DeleteLinkedParticlesPointersList(ParticleObject);
        if (ClearVoxels == true)
            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
        FreeIndexesOfParticles.push(ParticleIndex);
        Particles.erase(ParticleIndex);
    }
    CATCH("removing particle")
}

void CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::MakingZeroSizeForContainersForFoundParticlesInProximity()
{
    try
    {
        ParticlesKindsFoundInProximity.clear();
        ParticlesSortedByCapacityFoundInProximity.clear();
        NucleotidesWithFreeNextEndingsFoundInProximity.clear();
        NucleotidesWithFreePrevEndingsFoundInProximity.clear();
        DNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        DNANucleotidesWithFreePrevEndingsFoundInProximity.clear();
        RNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        RNANucleotidesWithFreePrevEndingsFoundInProximity.clear();
        NucleotidesFreeFoundInProximity.clear();
        RNANucleotidesFoundInProximity.clear();
    }
    CATCH("making zero size for containers for found particles in proximity")
}

void CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::UpdateFoundNucleotidesForFoundParticlesInProximity(const UnsignedInt ParticleIndex)
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

        if (CellEngineUseful::IsRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev != nullptr)
            RNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsRNA(ParticleObject.EntityId) && ParticleObject.Prev == nullptr && ParticleObject.Next != nullptr)
            RNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

        if (CellEngineUseful::IsDNAorRNA(ParticleObject.EntityId) && ParticleObject.Next == nullptr && ParticleObject.Prev == nullptr)
            NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
        if (CellEngineUseful::IsRNA(ParticleObject.EntityId))
            RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);

        ParticlesKindsFoundInProximity[ParticleObject.EntityId]++;
    }
    CATCH("updating found nucleotides for found particles in proximity")
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        set<UnsignedInt> FoundParticleIndexes;

        MakingZeroSizeForContainersForFoundParticlesInProximity();

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX >= 0 && PosY >= 0 && PosZ >= 0 && PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                        {
                            UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosY, PosZ);

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
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

void CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::PrintInformationAboutFoundParticlesInProximity()
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

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, const UnsignedInt AdditionalBoundFactor)
{
    try
    {
        LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId)));

        auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

        FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(true, ParticleObject.Center.X - ParticleKindObject.XSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Y - ParticleKindObject.YSizeDiv2 - AdditionalBoundFactor, ParticleObject.Center.Z - ParticleKindObject.ZSizeDiv2 - AdditionalBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalBoundFactor);
    }
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT 10 or 40 inside 1"));

            auto NucleotidePtr1 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }));

            if (NucleotidePtr1 != nullptr && NucleotidePtr1->Prev != nullptr)
            {
                CutDNANext(NucleotidePtr1->Prev);

                if (ReactionObject.Id == 40 || ReactionObject.Id == 41 || ReactionObject.Id == 42)
                {
                    auto NucleotidePtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter2, *GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).PairedNucleotide, false, false, [](const Particle*){ return true; }));
                    if (NucleotidePtr2 != nullptr && NucleotidePtr2->Prev != nullptr)
                    {
                        CutDNANext(NucleotidePtr2->Prev);

                        SeparateDNAStrands(&Particle::Next, ReactionObject.AdditionalParameter1 < ReactionObject.AdditionalParameter2 ? NucleotidePtr1 : NucleotidePtr2, abs(static_cast<SignedInt>(ReactionObject.AdditionalParameter2) - static_cast<SignedInt>(ReactionObject.AdditionalParameter1)));
                    }
                }
            }

            return true;
        }
    }
    CATCH("cutting dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 2)
        {
            LoggersManagerObject.Log(STREAM("LINK 1 inside 1"));

            Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first);
            Particle* NucleotideObjectForReactionPtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[1].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first), false, false, [](const Particle*){ return true; }));

            Particle* NucleotideObjectForReactionPtr1Inv = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first);
            Particle* NucleotideObjectForReactionPtr2Inv = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }));

            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2INV GENOME INDEX = " << NucleotideObjectForReactionPtr2Inv->GenomeIndex));

            if (NucleotideObjectForReactionPtr1->Prev == nullptr && NucleotideObjectForReactionPtr2->Next == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
            else
            if (NucleotideObjectForReactionPtr2Inv->Prev == nullptr && NucleotideObjectForReactionPtr1Inv->Next == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr2, NucleotideObjectForReactionPtr1);
            else
            if (NucleotideObjectForReactionPtr2Inv->Next == nullptr && NucleotideObjectForReactionPtr1Inv->Prev == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr1Inv, NucleotideObjectForReactionPtr2Inv);

            return true;
        }
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesWithFreePrevEndingsFoundInProximity.empty() == false && NucleotidesWithFreeNextEndingsFoundInProximity.empty() == false)
        {
            LoggersManagerObject.Log(STREAM("LINK 1 or 2 ANY inside 1"));

            if (DistanceOfParticles(GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]), GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0])) <= 2.0)
            {
                LoggersManagerObject.Log(STREAM("LINK 1 or 2 ANY CLOSE ENOUGH"));

                LinkDNA(&GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]), &GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0]));

                if (ReactionObject.Id == 80)
                    LinkDNA(GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]).PairedNucleotide, GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0]).PairedNucleotide);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT CRISPER inside 1"));

            Particle* FirstStrand = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.AdditionalParameter1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }))->Prev;
            if (FirstStrand != nullptr && FirstStrand->Prev != nullptr)
            {
                FirstStrand = FirstStrand->Prev;
                CutDNANext(FirstStrand);
                if (ReactionObject.Id == 110)
                    CutDNANext(FirstStrand->PairedNucleotide);
            }
        }
        else
            return true;
    }
    CATCH("cutting dna crisper in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesWithFreePrevEndingsFoundInProximity.empty() == false && NucleotidesWithFreeNextEndingsFoundInProximity.empty() == false)
        {
            LoggersManagerObject.Log(STREAM("LINK 2 ANY COMPATIBLE inside 1"));

            if (DistanceOfParticles(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]), GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0])) <= 2.0 && (DistanceOfParticles(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]), GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1])) <= 2.0))
            {
                LoggersManagerObject.Log(STREAM("LINK 2 ANY COMPATIBLE CLOSE ENOUGH - " << to_string(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]).GenomeIndex) << " "));

                LinkDNA(&GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]), &GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0]));
                LinkDNA(&GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]), &GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]));

                JoinDNAStrands(&Particle::Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).PairedNucleotide->Prev);
            }
        }
        else
            return true;
    }
    CATCH("linking dna ligase in any place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("LINK 2 REACTION 60"));

        auto NucleotidesIndexesChosenForReactionCopy(NucleotidesIndexesChosenForReaction);

        Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first);
        Particle* NucleotideObjectForReactionPtr1Paired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first).PairedNucleotide;
        Particle* NucleotideObjectForReactionPtr1Inv = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first);
        Particle* NucleotideObjectForReactionPtr1InvPaired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotide;

        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));

        if ((NucleotideObjectForReactionPtr1Inv->Next == nullptr || NucleotideObjectForReactionPtr1Inv->Prev == nullptr) && (NucleotideObjectForReactionPtr1->Next != nullptr && NucleotideObjectForReactionPtr1->Prev != nullptr))
        {
            swap(NucleotidesIndexesChosenForReactionCopy[0], NucleotidesIndexesChosenForReactionCopy[1]);

            swap(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr1Inv);
            swap(NucleotideObjectForReactionPtr1Paired, NucleotideObjectForReactionPtr1InvPaired);

            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 AFTER SWAP GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV AFTER SWAP GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));
        }

        Particle* NucleotideObjectForReactionPtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReactionCopy[1].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first), false, false, [](const Particle*){ return true; }));
        Particle* NucleotideObjectForReactionPtr2Paired = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReactionCopy[1].second].SequenceStr.length() - 1, *GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotide, false, false, [](const Particle*){ return true; }));

        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));

        if (NucleotideObjectForReactionPtr1Paired == nullptr && NucleotideObjectForReactionPtr2Paired != nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X1"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Next, 32, *NucleotideObjectForReactionPtr1, true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Next, 32, *(NucleotideObjectForReactionPtr2Paired->Next), true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtr1->GenomeIndex << " " << ParticlePtr1->PairedNucleotide->GenomeIndex << " " << ParticlePtrPrev2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtr1->PairedNucleotide, ParticlePtrPrev2);
                JoinDNAStrands(&Particle::Next, NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2Paired->Next);
            }
        }
        else
        if (NucleotideObjectForReactionPtr1Paired != nullptr && NucleotideObjectForReactionPtr2Paired == nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X2"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Prev, 32, *NucleotideObjectForReactionPtr1Paired->Prev, true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Prev, 32, *(NucleotideObjectForReactionPtr2), true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtrPrev1->GenomeIndex << " " << ParticlePtr2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtrPrev1, ParticlePtr2->PairedNucleotide);
                JoinDNAStrands(&Particle::Prev, NucleotideObjectForReactionPtr1Paired->Prev, NucleotideObjectForReactionPtr2);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("POLYMERASE DNA START REACTION"));

        auto& ParticleObject = GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first);

        if (NucleotidesFreeFoundInProximity.empty() == false && ParticleObject.LinkedParticlesPointersList.empty() == true)
        {
            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(ParticlesIndexesChosenForReaction[0].first) << " Nucleotide = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.ChainId) << " Nucleotide Index = " << ParticleObject.GenomeIndex));

            ParticleObject.AddNewLinkToParticle(&GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]));
            ParticleObject.AddNewLinkToParticle(GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).Next);

            MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
            MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]), ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing polymerase start dna special reaction function")

    return false;
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("POLYMERASE DNA CONTINUE REACTION"));

        auto& ParticleObject = GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first);

        for (auto& NucleotideFree : NucleotidesFreeFoundInProximity)
            LoggersManagerObject.Log(STREAM("FOUND FREE NUCLEOTIDES = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(NucleotideFree).ChainId)));

        LoggersManagerObject.Log(STREAM("Letter to find = " <<  CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->ChainId) << " NEXT " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->Next->ChainId) << " Particle Object Index = " << ParticleObject.Index << " " << ParticlesIndexesChosenForReaction[0].first ));

        auto ChosenNucleotideIterator = find_if(NucleotidesFreeFoundInProximity.cbegin(), NucleotidesFreeFoundInProximity.cend(), [this, ParticleObject](const UniqueIdInt& NucleotideParticleIndex){ return &GetParticleFromIndex(NucleotideParticleIndex) != ParticleObject.LinkedParticlesPointersList[0] && GetParticleFromIndex(NucleotideParticleIndex).ChainId == ParticleObject.LinkedParticlesPointersList[1]->ChainId; });
        if (ChosenNucleotideIterator != NucleotidesFreeFoundInProximity.end())
        {
            Particle* ChosenNucleotide = &GetParticleFromIndex(*ChosenNucleotideIterator);

            LoggersManagerObject.Log(STREAM("Letter to compare = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ChosenNucleotide->ChainId) << " " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->ChainId) << " NEXT " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->Next->ChainId)));

            LinkDNA(ChosenNucleotide, ParticleObject.LinkedParticlesPointersList[0]);
            ParticleObject.LinkedParticlesPointersList[0] = ChosenNucleotide;
            ParticleObject.LinkedParticlesPointersList[1] = ParticleObject.LinkedParticlesPointersList[1]->Next;

            MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
            MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(*ChosenNucleotide, ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing polymerase continue dna special reaction function")

    return false;
}



tuple<vector<ChainIdInt>, string> CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, const UnsignedInt SizeOfLoop)
{
    string TemplateSequenceStr;

    vector<ChainIdInt> TemplateSequence;

    try
    {
        UnsignedInt LengthOfTemplateForRNA =  32;

        UnsignedInt NucleotidesFoundInProximityCounter = 0;

        while (TemplateSequenceStr.empty() == true && NucleotidesFoundInProximityCounter < SizeOfLoop)
        {
            auto [ParticlePtrNext, ParticlePtrPrevNext, NucleotidesCounterNext, NucleotidesSequenceToCompareStringNext, NucleotidesSequenceToCompareVectorNext] = GetNucleotidesSequence(&Particle::Next, LengthOfTemplateForRNA, GetParticleFromIndex(NucleotidesFoundInProximity[NucleotidesFoundInProximityCounter]), true, true, [](const Particle *){ return true; });
            auto [ParticlePtrBack, ParticlePtrBackNext, NucleotidesCounterBack, NucleotidesSequenceToCompareStringBack, NucleotidesSequenceToCompareVectorBack] = GetNucleotidesSequence(&Particle::Prev, LengthOfTemplateForRNA, GetParticleFromIndex(NucleotidesFoundInProximity[NucleotidesFoundInProximityCounter]), true, true, [](const Particle *){ return true; });

            reverse(NucleotidesSequenceToCompareVectorBack.begin(), NucleotidesSequenceToCompareVectorBack.end());
            TemplateSequence = NucleotidesSequenceToCompareVectorBack;
            TemplateSequence.insert(end(TemplateSequence), begin(NucleotidesSequenceToCompareVectorNext), end(NucleotidesSequenceToCompareVectorNext));
            for (auto &Nucleotide: TemplateSequence)
                Nucleotide = CellEngineUseful::GetPairedChainIdForDNAorRNA(Nucleotide);

            reverse(NucleotidesSequenceToCompareStringBack.begin(), NucleotidesSequenceToCompareStringBack.end());
            TemplateSequenceStr = NucleotidesSequenceToCompareStringBack + NucleotidesSequenceToCompareStringNext;
            string OriginalTemplateRNASequenceStr = TemplateSequenceStr;
            for (auto &TemplateSequenceNucleotideChar: TemplateSequenceStr)
                TemplateSequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(TemplateSequenceNucleotideChar)));

            NucleotidesFoundInProximityCounter++;
        }
    }
    CATCH("getting nucleotides sequence in both directions")

    return { TemplateSequence, TemplateSequenceStr };
}

bool CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    bool FoundSequenceNotFit = false;

    try
    {
        string TemplateSequenceStr = ParticleKindForReactionObject.SequenceStr;

        if (TemplateSequenceStr == "ANY")
        {
            LoggersManagerObject.Log(STREAM("DNA SEQUENCE FOUND = ANY => ParticleKindForReactionObject.EntityId = " << to_string(ParticleKindForReactionObject.EntityId)));
            return true;
        }

        vector<ChainIdInt> TemplateSequence = ParticleKindForReactionObject.Sequence;

        string OriginalTemplateRNASequenceStr;

        if (TemplateSequenceStr == "RNA")
            if (RNANucleotidesFoundInProximity.empty() == false)
            {
                tie(TemplateSequence, TemplateSequenceStr) = GetNucleotidesSequenceInBothDirections(RNANucleotidesFoundInProximity, RNANucleotidesFoundInProximity.size());

                if (TemplateSequenceStr.empty() == true)
                    return false;
            }

        auto [ParticlePtr, ParticlePtrPrev, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector] = GetNucleotidesSequence(&Particle::Next, TemplateSequenceStr.length(), ParticleObjectForReaction, true, true, [](const Particle*){ return true; });

        LoggersManagerObject.Log(STREAM("DNA SEQUENCE COMPARE = #" << NucleotidesSequenceToCompareString << "#" << TemplateSequenceStr << "#" << OriginalTemplateRNASequenceStr << "#" << to_string(ParticleKindForReactionObject.EntityId)));

        if (TypeOfComparison == ComparisonType::ByVectorLoop)
        {
            if (NucleotidesSequenceToCompareVector.size() >= TemplateSequence.size())
            {
                LoggersManagerObject.Log(STREAM("LOOP COMPARISON SIZE = " << to_string(NucleotidesSequenceToCompareVector.size()) << " " << to_string(TemplateSequence.size())));

                for (UnsignedInt NucleotideNum = 0; NucleotideNum < TemplateSequence.size(); NucleotideNum++)
                    if (TemplateSequence[NucleotideNum] != NucleotidesSequenceToCompareVector[NucleotideNum])
                    {
                        LoggersManagerObject.Log(STREAM("LOOP COMPARISON BREAK = " << to_string(NucleotideNum) << "#"));

                        FoundSequenceNotFit = true;
                        break;
                    }
            }
            else
                FoundSequenceNotFit = true;
        }
        else
        if (TypeOfComparison == ComparisonType::ByString)
            FoundSequenceNotFit = !(NucleotidesSequenceToCompareString == TemplateSequenceStr);

        if (FoundSequenceNotFit == false)
            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "DNA SEQUENCE FOUND FIT" << terminal_colors_utils::white));
    }
    CATCH("comparing fitness of dna sequence by nucleotides loop")

    return !FoundSequenceNotFit;
}

tuple<vector<pair<UniqueIdInt, UnsignedInt>>, bool> CellEngineNucleicAcidsChemicalReactionsInVoxelSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject)
{
    bool AllAreZero = false;

    vector<pair<UniqueIdInt, UnsignedInt>> NucleotidesIndexesChosenForReaction, ParticlesIndexesChosenForReaction, AllParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : ParticlesSortedByCapacityFoundInProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z)));

            vector<ParticleKindForReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsSpecialDNA(ParticleKindForReactionObjectParam.EntityId) && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByVectorLoop, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

            auto PositionInReactants = ReactantIterator - ReactionObject.Reactants.begin();

            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == true)
                if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0)
            {
                AllParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);
                LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.Center.X) << " Y = " << to_string(ParticleObjectTestedForReaction.Center.Y) << " Z = " << to_string(ParticleObjectTestedForReaction.Center.Z) << endl));
                ReactantsCounters[PositionInReactants]--;
            }

            AllAreZero = all_of(ReactantsCounters.begin(), ReactantsCounters.end(), [this](const UnsignedInt& Counter){ return Counter == 0; });
            if (AllAreZero == true)
            {
                LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));
                break;
            }
            LoggersManagerObject.Log(STREAM(""));
        }

        if (AllAreZero == true || (AllAreZero == false && (ReactionObject.Id == 30 || ReactionObject.Id == 80 || ReactionObject.Id == 70)))
            ReactionObject.SpecialReactionFunction(this, AllParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
    }
    CATCH("choosing particles for reaction from all particles in proximity")

    if (AllAreZero == true)
    {
        LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size())));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { vector<pair<UniqueIdInt, UnsignedInt>>(), false };
}
