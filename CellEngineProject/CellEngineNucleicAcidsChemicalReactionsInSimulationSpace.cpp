
#include "CellEngineNucleicAcidsChemicalReactionsInSimulationSpace.h"

tuple<vector<ChainIdInt>, string> CellEngineNucleicAcidsChemicalReactionsInSimulationSpace::GetNucleotidesSequenceInBothDirections(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, const UnsignedInt SizeOfLoop)
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
            for (auto &TemplateSequenceNucleotideChar: TemplateSequenceStr)
                TemplateSequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(TemplateSequenceNucleotideChar)));

            NucleotidesFoundInProximityCounter++;
        }
    }
    CATCH("getting nucleotides sequence in both directions")

    return { TemplateSequence, TemplateSequenceStr };
}

tuple<vector<ChainIdInt>, string> CellEngineNucleicAcidsChemicalReactionsInSimulationSpace::GetNucleotidesSequenceFromRNAInOneParticle(const std::vector<UniqueIdInt>& NucleotidesFoundInProximity, const UnsignedInt SizeOfLoop)
{
    string TemplateSequenceStr;

    vector<ChainIdInt> TemplateSequence;

    try
    {
        TemplateSequenceStr = GetParticleFromIndex(*LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.begin()).SequenceStr;

        for (auto &TemplateSequenceNucleotideChar: TemplateSequenceStr)
            TemplateSequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(TemplateSequenceNucleotideChar)));

        for (auto &TemplateSequenceNucleotideChar: TemplateSequenceStr)
            TemplateSequence.emplace_back(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(TemplateSequenceNucleotideChar));
    }
    CATCH("getting nucleotides sequence in both directions")

    return { TemplateSequence, TemplateSequenceStr };
}

bool CellEngineNucleicAcidsChemicalReactionsInSimulationSpace::CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForChemicalReaction& ParticleKindForReactionObject, Particle& ParticleObjectTestedForReaction)
{
    bool FoundSequenceNotFit = false;
    bool SpecialCompareFunctionResult = true;

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
            if (LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.empty() == false)
            {
                if (CellEngineConfigDataObject.RNAInOneParticle == false)
                    tie(TemplateSequence, TemplateSequenceStr) = GetNucleotidesSequenceInBothDirections(LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity, LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.size());
                else
                    tie(TemplateSequence, TemplateSequenceStr) = GetNucleotidesSequenceFromRNAInOneParticle(LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity, LocalThreadParticlesInProximityObject.RNANucleotidesFoundInProximity.size());

                if (TemplateSequenceStr.empty() == true)
                    return false;
            }

        auto [ParticlePtr, ParticlePtrPrev, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector] = GetNucleotidesSequence(&Particle::Next, TemplateSequenceStr.length(), ParticleObjectTestedForReaction, true, true, [](const Particle*){ return true; });

        LoggersManagerObject.Log(STREAM("DNA SEQUENCE COMPARE = #" << NucleotidesSequenceToCompareString << "#" << TemplateSequenceStr << "#" << OriginalTemplateRNASequenceStr << "#" << to_string(ParticleKindForReactionObject.EntityId)));

        if (TypeOfComparison == ComparisonType::ByVectorLoop)
            CellEngineUseful::CompareSequences(TemplateSequence, NucleotidesSequenceToCompareVector, FoundSequenceNotFit, false);
        else
        if (TypeOfComparison == ComparisonType::ByString)
            FoundSequenceNotFit = !(NucleotidesSequenceToCompareString == TemplateSequenceStr);

        if (FoundSequenceNotFit == false)
            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "DNA SEQUENCE FOUND FIT" << terminal_colors_utils::white));

        if (FoundSequenceNotFit == false)
            if (ParticleKindForReactionObject.SpecialCompareFunction != nullptr)
                SpecialCompareFunctionResult = ParticleKindForReactionObject.SpecialCompareFunction(ParticleObjectTestedForReaction.GenomeIndex);
    }
    CATCH("comparing fitness of dna sequence by nucleotides loop")

    return !FoundSequenceNotFit && SpecialCompareFunctionResult;
}