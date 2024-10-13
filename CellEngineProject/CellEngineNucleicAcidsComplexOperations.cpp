
#include "CellEngineConstants.h"
#include "CellEngineParticle.h"
#include "CellEngineAminoAcids.h"
#include "CellEngineNucleicAcidsComplexOperations.h"

bool CellEngineNucleicAcidsComplexOperations::CutDNAInChosenPlaceSpecialReactionFunction(const bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT inside 1"));

            auto NucleotidePtr1 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }));

            if (NucleotidePtr1 != nullptr && NucleotidePtr1->Prev != nullptr)
            {
                CutDNANext(NucleotidePtr1->Prev);

                if (BothStrandsBool == true)
                {
                    auto NucleotidePtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter2, *GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).PairedNucleotidePtr, false, false, [](const Particle*){ return true; }));
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

bool CellEngineNucleicAcidsComplexOperations::CutDNASingleStrandInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return CutDNAInChosenPlaceSpecialReactionFunction(false, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::CutDNABothStrandsInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return CutDNAInChosenPlaceSpecialReactionFunction(true, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
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


bool CellEngineNucleicAcidsComplexOperations::LinkDNAInAnyPlaceSpecialReactionFunction(const bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
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

                if (BothStrandsBool == true)
                    LinkDNA(GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]).PairedNucleotidePtr, GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0]).PairedNucleotidePtr);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::LinkDNASingleStrandInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return LinkDNAInAnyPlaceSpecialReactionFunction(false, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::LinkDNABothStrandsInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return LinkDNAInAnyPlaceSpecialReactionFunction(true, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::CutDNACrisperInChosenPlaceSpecialReactionFunction(const bool BothStrandsBool, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
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
                if (BothStrandsBool == true)
                    CutDNANext(FirstStrand->PairedNucleotidePtr);
            }
        }
        else
            return true;
    }
    CATCH("cutting dna crisper in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::CutDNASingleStrandCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return CutDNACrisperInChosenPlaceSpecialReactionFunction(false, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::CutDNABothStrandsCrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    return CutDNACrisperInChosenPlaceSpecialReactionFunction(true, ParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
}

bool CellEngineNucleicAcidsComplexOperations::LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
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

                if (GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]).GenomeIndex < GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).GenomeIndex)
                    JoinDNAStrands(&Particle::Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).PairedNucleotidePtr->Prev);
                else
                    JoinDNAStrands(&Particle::Next, GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]).Next, GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]).PairedNucleotidePtr->Next);
            }
        }
        else
            return true;
    }
    CATCH("linking dna ligase in any place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("LINK 2 REACTION 60"));

        auto NucleotidesIndexesChosenForReactionCopy(NucleotidesIndexesChosenForReaction);

        Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first);
        Particle* NucleotideObjectForReactionPtr1Paired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first).PairedNucleotidePtr;
        Particle* NucleotideObjectForReactionPtr1Inv = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first);
        Particle* NucleotideObjectForReactionPtr1InvPaired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotidePtr;

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
        Particle* NucleotideObjectForReactionPtr2Paired = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReactionCopy[1].second].SequenceStr.length() - 1, *GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotidePtr, false, false, [](const Particle*){ return true; }));

        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));

        if (NucleotideObjectForReactionPtr1Paired == nullptr && NucleotideObjectForReactionPtr2Paired != nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X1"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Next, 32, *NucleotideObjectForReactionPtr1, true, false, [](const Particle* P){ return P->PairedNucleotidePtr == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Next, 32, *(NucleotideObjectForReactionPtr2Paired->Next), true, false, [](const Particle* P){ return P->PairedNucleotidePtr == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtr1->GenomeIndex << " " << ParticlePtr1->PairedNucleotidePtr->GenomeIndex << " " << ParticlePtrPrev2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtr1->PairedNucleotidePtr, ParticlePtrPrev2);
                JoinDNAStrands(&Particle::Next, NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2Paired->Next);
            }
        }
        else
        if (NucleotideObjectForReactionPtr1Paired != nullptr && NucleotideObjectForReactionPtr2Paired == nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X2"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Prev, 32, *NucleotideObjectForReactionPtr1Paired->Prev, true, false, [](const Particle* P){ return P->PairedNucleotidePtr == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Prev, 32, *(NucleotideObjectForReactionPtr2), true, false, [](const Particle* P){ return P->PairedNucleotidePtr == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtrPrev1->GenomeIndex << " " << ParticlePtr2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtrPrev1, ParticlePtr2->PairedNucleotidePtr);
                JoinDNAStrands(&Particle::Prev, NucleotideObjectForReactionPtr1Paired->Prev, NucleotideObjectForReactionPtr2);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::PolymeraseRNATranscriptionStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("POLYMERASE RNA TRANSCRIPTION START REACTION"));

        auto LocalParticlesIndexesChosenForReaction(ParticlesIndexesChosenForReaction);
        for (const auto& ParticleChosenForReactionObject : LocalParticlesIndexesChosenForReaction)
            LoggersManagerObject.Log(STREAM("ParticleIndexA = " << to_string(ParticleChosenForReactionObject.first) << " " << to_string(ParticleChosenForReactionObject.second) << " Particle Type = " << ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(ParticleChosenForReactionObject.first).EntityId).IdStr));

        sort(LocalParticlesIndexesChosenForReaction.begin(), LocalParticlesIndexesChosenForReaction.end(), [](const auto& PK1, const auto& PK2) { return PK1.second < PK2.second; });

        if (auto& ParticleObject = GetParticleFromIndex(LocalParticlesIndexesChosenForReaction[0].first); RNANucleotidesFreeFoundInProximity.empty() == false && ParticleObject.LinkedParticlesPointersList.empty() == true)
        {
            for (const auto& ParticleChosenForReactionObject : LocalParticlesIndexesChosenForReaction)
                LoggersManagerObject.Log(STREAM("ParticleIndexB = " << to_string(ParticleChosenForReactionObject.first) << " " << to_string(ParticleChosenForReactionObject.second) << " Particle Type = " << ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(ParticleChosenForReactionObject.first).EntityId).IdStr));

            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticlesIndexesChosenForReaction[0].first) << " Nucleotide = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.ChainId) << " Nucleotide Index = " << ParticleObject.GenomeIndex));

            ParticleObject.AddNewLinkToParticle(&GetParticleFromIndex(RNANucleotidesFreeFoundInProximity[0]));
            ParticleObject.AddNewLinkToParticle(GoSomeNucleotides(&Particle::Next, ReactionObject.Reactants[1].SequenceStr.length(), *GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).Next));

            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ParticleObject.LinkedParticlesPointersList[0], ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing polymerase start dna transcription special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::PolymeraseRNATranscriptionContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        string SequenceOfLettersToCheckFinishSequence;

        LoggersManagerObject.Log(STREAM("POLYMERASE RNA TRANSCRIPTION CONTINUE REACTION"));

        auto& ParticleObject = GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first);

        for (const auto &NucleotideFree: RNANucleotidesFullFreeFoundInProximity)
            LoggersManagerObject.Log(STREAM("FOUND FREE NUCLEOTIDES = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(NucleotideFree).ChainId)));

        LoggersManagerObject.Log(STREAM("Letter to find = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->ChainId) << " NEXT " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->Next->ChainId) << " Particle Object Index = " << ParticleObject.Index << " " << ParticlesIndexesChosenForReaction[0].first));

        auto ChosenNucleotideIterator = find_if(RNANucleotidesFullFreeFoundInProximity.cbegin(), RNANucleotidesFullFreeFoundInProximity.cend(), [this, ParticleObject](const UniqueIdInt &NucleotideParticleIndex){ return &GetParticleFromIndex(NucleotideParticleIndex) != ParticleObject.LinkedParticlesPointersList[0] && CellEngineUseful::IsRNANucleotidePairedForRNAEqual(GetParticleFromIndex(NucleotideParticleIndex).EntityId, ParticleObject.LinkedParticlesPointersList[1]->ChainId) == true; });
        if (ChosenNucleotideIterator != RNANucleotidesFullFreeFoundInProximity.end())
        {
            Particle *ChosenNucleotide = &GetParticleFromIndex(*ChosenNucleotideIterator);

            LoggersManagerObject.Log(STREAM("Letter to compare = " << CellEngineUseful::GetLetterFromNucleotideForDNAorRNA(ChosenNucleotide->EntityId) << " " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->ChainId) << " NEXT " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticleObject.LinkedParticlesPointersList[1]->Next->ChainId)));

            if (CellEngineConfigDataObject.RNAInOneParticle == false)
            {
                LinkDNA(ChosenNucleotide, ParticleObject.LinkedParticlesPointersList[0]);
                ParticleObject.LinkedParticlesPointersList[0] = ChosenNucleotide;
                ParticleObject.LinkedParticlesPointersList[1] = ParticleObject.LinkedParticlesPointersList[1]->Next;

                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ChosenNucleotide, ParticleObject, 2, 2, 2);

                SequenceOfLettersToCheckFinishSequence = get<3>(GetNucleotidesSequence(&Particle::Prev, MaxLengthOfGene, *(ParticleObject.LinkedParticlesPointersList[1]), true, false, [](const Particle* P){ return true; }));
            }
            else
            {
                ParticleObject.LinkedParticlesPointersList[0]->SequenceStr += CellEngineUseful::GetLetterFromNucleotideForDNAorRNA(ChosenNucleotide->EntityId);
                LoggersManagerObject.Log(STREAM("SequenceStr = " << ParticleObject.LinkedParticlesPointersList[0]->SequenceStr));

                ParticleObject.LinkedParticlesPointersList[1] = ParticleObject.LinkedParticlesPointersList[1]->Next;

                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ParticleObject.LinkedParticlesPointersList[0], *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
                vector<vector3_16> Centers;
                EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(*ChosenNucleotideIterator, Centers);

                SequenceOfLettersToCheckFinishSequence = ParticleObject.LinkedParticlesPointersList[0]->SequenceStr;
            }
        }
        else
            LoggersManagerObject.Log(STREAM("Nucleotide not found - " << NucleotidesFreeFoundInProximity.size()));

        if (SequenceOfLettersToCheckFinishSequence.length() > 0 && SequenceOfLettersToCheckFinishSequence.length() % 3 == 0)
            if (CellEngineUseful::IsIn(SequenceOfLettersToCheckFinishSequence.substr(SequenceOfLettersToCheckFinishSequence.length() - 3, 3), { "UAG", "UAA", "UGA" }))
            {
                ParticleObject.LinkedParticlesPointersList.clear();

                LoggersManagerObject.Log(STREAM("DETACHED FROM DNA AND RNA AFTER PROCESS OF TRANSCRIPTION"));
            }
    }
    CATCH("executing polymerase continue dna transcription special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::RibosomeTranslationStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("RIBOSOME TRANSLATION START REACTION"));

        auto LocalParticlesIndexesChosenForReaction(ParticlesIndexesChosenForReaction);

        for (const auto& ParticleChosenForReactionObject : LocalParticlesIndexesChosenForReaction)
            LoggersManagerObject.Log(STREAM("ParticleIndex before sort = " << to_string(ParticleChosenForReactionObject.first) << " " << to_string(ParticleChosenForReactionObject.second) << " Particle Type = " << ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(ParticleChosenForReactionObject.first).EntityId).IdStr));

        sort(LocalParticlesIndexesChosenForReaction.begin(), LocalParticlesIndexesChosenForReaction.end(), [](const auto& PK1, const auto& PK2) { return PK1.second < PK2.second; });

        if (auto& ParticleObject = GetParticleFromIndex(LocalParticlesIndexesChosenForReaction[0].first); ParticleObject.LinkedParticlesPointersList.empty() == true)
        {
            for (const auto& ParticleChosenForReactionObject : LocalParticlesIndexesChosenForReaction)
                LoggersManagerObject.Log(STREAM("ParticleIndex after sort = " << to_string(ParticleChosenForReactionObject.first) << " " << to_string(ParticleChosenForReactionObject.second) << " Particle Type = " << ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(ParticleChosenForReactionObject.first).EntityId).IdStr));

            ParticleObject.AddNewLinkToParticle(&GetParticleFromIndex(LocalParticlesIndexesChosenForReaction[1].first));
            ParticleObject.AddNewLinkToParticle(&GetParticleFromIndex(RNANucleotidesFreeFoundInProximity[0]));

            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ParticleObject.LinkedParticlesPointersList[0], ParticleObject, 2, 2, 2);
            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ParticleObject.LinkedParticlesPointersList[1], ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing ribosome start dna translation special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::RibosomeTranslationContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>> &NucleotidesIndexesChosenForReaction, const ChemicalReaction& ReactionObject)
{
    try
    {
        string SequenceOfLettersToCheckFinishSequence;

        LoggersManagerObject.Log(STREAM("RIBOSOME TRANSLATION CONTINUE REACTION"));

        LoggersManagerObject.Log(STREAM("SIZE OF tRNA charged found = " << tRNAChargedFoundInProximity.size()));

        for (const auto &tRNAChargedFree: tRNAChargedFoundInProximity)
            LoggersManagerObject.Log(STREAM("FOUND FREE CHARGED TRNA = " << ParticlesKindsManagerObject.GetParticleKind(GetParticleFromIndex(tRNAChargedFree).EntityId).IdStr));

        auto& ParticleObject = GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first);

        if (CellEngineConfigDataObject.RNAInOneParticle == false)
        {
        }
        else
        {
            auto Codon = ParticleObject.LinkedParticlesPointersList[1]->SequenceStr.substr(RNAStartSequence.length() + ParticleObject.LinkedParticlesPointersList[1]->PositionInSequence, 3);

            LoggersManagerObject.Log(STREAM("Codon = " << Codon << " " << ParticleObject.LinkedParticlesPointersList[1]->PositionInSequence));

            auto ChosentRNAChargedIterator = find_if(tRNAChargedFoundInProximity.begin(), tRNAChargedFoundInProximity.end(), [this, Codon](const UniqueIdInt &tRNAParticleIndex){ return CellEngineAminoAcidsManagerObject.IstRNAChargedWithAminoAcidForCodon(GetParticleFromIndex(tRNAParticleIndex).EntityId, Codon); });

            if (ChosentRNAChargedIterator != tRNAChargedFoundInProximity.end())
            {
                ParticleObject.LinkedParticlesPointersList[0]->SequenceStr += Codon;
                ParticleObject.LinkedParticlesPointersList[1]->PositionInSequence += 3;
                LoggersManagerObject.Log(STREAM("SequenceStr = " << ParticleObject.LinkedParticlesPointersList[0]->SequenceStr));

                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
                MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ParticleObject.LinkedParticlesPointersList[0], *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
                vector<vector3_16> Centers;
                EraseParticleChosenForReactionAndGetCentersForNewProductsOfReaction(*ChosentRNAChargedIterator, Centers);

                SequenceOfLettersToCheckFinishSequence = ParticleObject.LinkedParticlesPointersList[0]->SequenceStr;
            }
            else
                LoggersManagerObject.Log(STREAM("trna not found - " << NucleotidesFreeFoundInProximity.size()));
        }

        if (SequenceOfLettersToCheckFinishSequence.length() > 0 && SequenceOfLettersToCheckFinishSequence.size() % 3 == 0)
            if (CellEngineUseful::IsIn(SequenceOfLettersToCheckFinishSequence.substr(SequenceOfLettersToCheckFinishSequence.length() - 3, 3), { "UAG", "UAA", "UGA" }))
            {
                ParticleObject.LinkedParticlesPointersList.clear();

                LoggersManagerObject.Log(STREAM("DETACHED FROM RNA AND PROTEIN AFTER PROCESS OF TRANSLATION"));
            }
    }
    CATCH("executing ribosome continue dna translation special reaction function")

    return false;
}

