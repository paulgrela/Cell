
#include "CellEngineParticle.h"
#include "CellEngineNucleicAcidsComplexOperations.h"

bool CellEngineNucleicAcidsComplexOperations::CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

bool CellEngineNucleicAcidsComplexOperations::PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]), ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing polymerase start dna special reaction function")

    return false;
}

bool CellEngineNucleicAcidsComplexOperations::PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
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

            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(ParticleObject, *ParticleObject.LinkedParticlesPointersList[1], 2, 2, 2);
            MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(*ChosenNucleotide, ParticleObject, 2, 2, 2);
        }
    }
    CATCH("executing polymerase continue dna special reaction function")

    return false;
}