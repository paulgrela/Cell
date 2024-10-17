
#include <string>
#include <sstream>

#include "CellEngineTypes.h"
#include "CellEngineUseful.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineGenesPromotersAndGenesStartCodonsFinder.h"

using namespace std;

void FindInterGenesSequencesFromGenesData()
{
    try
    {
        vector<Gene> LocalGenesVector;
        transform(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), back_inserter(LocalGenesVector), [](const auto& Gene){ return Gene.second; });
        sort(LocalGenesVector.begin(), LocalGenesVector.end(), [](const Gene& G1, const Gene& G2){ return G1.StartPosInGenome < G2.StartPosInGenome; });

        for (const auto& Gene : LocalGenesVector)
        {
            UnsignedInt StartPos = Gene.StartPosInGenome;
            UnsignedInt EndPos = Gene.EndPosInGenome;
            LoggersManagerObject.Log(STREAM("Gene = " << Gene.NumId << " " << StartPos << " " << EndPos));
        }
    }
    CATCH("finding inter genes sequences from genes")
}

inline void FindCodons(const string& GenomeStr, const vector<string>& StartCodonsList, map<GeneIdInt, UnsignedInt>& FoundGenesCounter, UnsignedInt& NumberOfFoundPromoters, UnsignedInt& NumberOfFoundPromotersConfirmed, const UnsignedInt ScopeForStartGeneIndex, const UnsignedInt PosInGenomeToFindCodonStart, const UnsignedInt PosInGenomeToFind35Box, const UnsignedInt PosInGenomeToFind10Box, const bool SwitchLogsBool)
{
    try
    {
        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOffLogs();

        for (const string& StartCodon : StartCodonsList)
            if (GenomeStr.substr(PosInGenomeToFindCodonStart, StartCodon.size()) == StartCodon)
            {
                LoggersManagerObject.Log(STREAM("Promoter found at -35: " << PosInGenomeToFind35Box << " and -10: " << PosInGenomeToFind10Box));
                LoggersManagerObject.Log(STREAM("Start codon (" << StartCodon << ") found at: " << PosInGenomeToFindCodonStart << std::endl));

                auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&PosInGenomeToFindCodonStart, &ScopeForStartGeneIndex](const auto& G1) { return PosInGenomeToFindCodonStart - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= PosInGenomeToFindCodonStart + ScopeForStartGeneIndex; });
                if (FindResult != ParticlesKindsManagerObject.Genes.end())
                {
                    LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));

                    FoundGenesCounter[FindResult->second.NumId]++;
                    NumberOfFoundPromotersConfirmed++;
                }
                NumberOfFoundPromoters++;

                break;
            }

        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOnLogs();
    }
    CATCH("finding codons")
}

void FindPromotersAndStartCodons1(const string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        const string Promoter35BoxesList = "TTGACA";
        const string Promoter10BoxesList = "TATAA";
        const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

        const SignedInt MaxDistanceBetween35and10Boxes = 50;
        const SignedInt MaxDistanceFrom10BoxToStartCodon = 50;

        UnsignedInt NumberOfFoundPromoters = 0;
        UnsignedInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 5;

        for (UnsignedInt PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - Promoter35BoxesList.size(); ++PosInGenomeToFind35Box)
            if (GenomeStr.substr(PosInGenomeToFind35Box, Promoter35BoxesList.size()) == Promoter35BoxesList)
                for (UnsignedInt PosInGenomeToFind10Box = PosInGenomeToFind35Box + Promoter35BoxesList.size(); PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size(); ++PosInGenomeToFind10Box)
                    if (GenomeStr.substr(PosInGenomeToFind10Box, Promoter10BoxesList.size()) == Promoter10BoxesList)
                    {
                        for (UnsignedInt PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + Promoter10BoxesList.size(); PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + Promoter10BoxesList.size() + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size(); ++PosInGenomeToFindCodonStart)
                            FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                        break;
                    }

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Found ConfirmedPromoters = " << NumberOfFoundPromotersConfirmed << endl));
    }
    CATCH("finding promoters and start codons 1")
}

inline SignedInt CalculateDifference(const std::string& Sequence1, const std::string& Sequence2)
{
    SignedInt Mismatches = 0;

    for (UnsignedInt Pos = 0; Pos < Sequence1.size(); ++Pos)
        if (Sequence1[Pos] != Sequence2[Pos])
            Mismatches++;

    return Mismatches;
}

void FindPromotersAndStartCodons2(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
        const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
        const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

        const SignedInt MaxDistanceBetween35and10Boxes = 50;
        const SignedInt MaxDistanceFrom10BoxToStartCodon = 100;
        const SignedInt MaxNumberOfMismatchesForBox35 = 1;
        const SignedInt MaxNumberOfMismatchesForBox10 = 1;

        UnsignedInt NumberOfFoundPromoters = 0;
        UnsignedInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 5;

        for (UnsignedInt PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - 6; ++PosInGenomeToFind35Box)
            for (const string& Promoter35Box : Promoter35BoxesList)
                if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) <= MaxNumberOfMismatchesForBox35)
                    for (UnsignedInt PosInGenomeToFind10Box = PosInGenomeToFind35Box + 6; PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                        for (const string& Promoter10Box : Promoter10BoxesList)
                            if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) <= MaxNumberOfMismatchesForBox10)
                            {
                                for (UnsignedInt PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + 6; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + 6 + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                                    FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                                break;
                            }

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
        LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
    }
    CATCH("finding promoters and start codons 2")
}

inline SignedInt CalculateSimilarity(const string& Sequence1, const string& Sequence2)
{
    SignedInt Matches = 0;

    for (UnsignedInt Pos = 0; Pos < Sequence1.size(); ++Pos)
        if (Sequence1[Pos] == Sequence2[Pos])
            Matches++;

    return Matches;
}

void FindPromotersAndStartCodons3(const string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
        const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
        const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

        const SignedInt MaxDistanceBetween35and10Boxes = 19;
        const SignedInt MaxDistanceFrom10BoxToStartCodon = 50;
        const SignedInt MinDistanceFrom10BoxToStartCodon = 10;
        const SignedInt MinimumSimilarityScoreFor35Box = 5;
        const SignedInt MinimumSimilarityScoreFor10Box = 5;

        UnsignedInt NumberOfFoundPromoters = 0;
        UnsignedInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 3;

        for (UnsignedInt PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - 6; ++PosInGenomeToFind35Box)
            for (const string& Promoter35Box : Promoter35BoxesList)
                if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
                    for (UnsignedInt PosInGenomeToFind10Box = PosInGenomeToFind35Box + 16; PosInGenomeToFind10Box <= PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                        for (const string& Promoter10Box : Promoter10BoxesList)
                            if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                            {
                                for (UnsignedInt PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + MinDistanceFrom10BoxToStartCodon; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                                    FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                                break;
                            }

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
        LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
    }
    CATCH("finding promoters and start codons 3")
}

void FindPromotersForGenesFromGeneStartPos(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOffLogs();

        std::vector<pair<GeneIdInt, UnsignedInt>> GenesStarts;
        for (const auto& Gene : ParticlesKindsManagerObject.Genes)
            GenesStarts.emplace_back(Gene.second.NumId, Gene.second.StartPosInGenome);

        const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
        const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };

        const UnsignedInt MinDistanceFrom10BoxToStartCodon = 10;
        const UnsignedInt MaxDistanceFrom10BoxToStartCodon = 50;
        const UnsignedInt MaxDistanceBetween35and10Boxes = 19;
        const UnsignedInt MinimumSimilarityScoreFor35Box = 5;
        const UnsignedInt MinimumSimilarityScoreFor10Box = 5;

        UnsignedInt NumberOfFoundPromoters = 0;
        map<GeneIdInt, UnsignedInt> FoundPromotersCounter;

        for (const auto& GeneStartPos : GenesStarts)
        {
            bool PromoterFound = false;

            for (SignedInt PosInGenomeToFind10Box = static_cast<SignedInt>(GeneStartPos.second) - static_cast<SignedInt>(MinDistanceFrom10BoxToStartCodon); PosInGenomeToFind10Box >= std::max(0, static_cast<int>(GeneStartPos.second - MaxDistanceFrom10BoxToStartCodon)); --PosInGenomeToFind10Box)
                for (const string& Promoter10Box : Promoter10BoxesList)
                    if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                        for (SignedInt PosInGenomeToFind35Box = PosInGenomeToFind10Box - static_cast<SignedInt>(MaxDistanceBetween35and10Boxes); PosInGenomeToFind35Box >= std::max(0, static_cast<int>(PosInGenomeToFind10Box) - static_cast<int>(MaxDistanceBetween35and10Boxes) - 5); --PosInGenomeToFind35Box)
                            for (const string& Promoter35Box : Promoter35BoxesList)
                                if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
                                {
                                    LoggersManagerObject.Log(STREAM("Promoter found for gene starting at " << GeneStartPos.second));
                                    LoggersManagerObject.Log(STREAM("-35 box at position: " << PosInGenomeToFind35Box));
                                    LoggersManagerObject.Log(STREAM("-10 box at position: " << PosInGenomeToFind10Box << endl));

                                    ParticlesKindsManagerObject.Promoters.emplace(PosInGenomeToFind10Box, Promoter{ GeneStartPos.first, static_cast<UnsignedInt>(PosInGenomeToFind10Box), static_cast<UnsignedInt>(PosInGenomeToFind35Box), GeneStartPos.second });
                                    FoundPromotersCounter[PosInGenomeToFind10Box]++;
                                    NumberOfFoundPromoters++;
                                    PromoterFound = true;

                                    goto Outside;
                                }
            Outside:

            if (PromoterFound == false)
                LoggersManagerObject.Log(STREAM("No promoter found for gene starting at " << GeneStartPos.second << std::endl));
        }

        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundPromotersCounter.size() << endl));
    }
    CATCH("finding promoters for genes from gene start position")
}

void FindPromoters(const std::vector<std::string>& GenomesLines, const std::vector<std::vector<UniqueIdInt>>& Genomes, const bool SwitchLogsBool)
{
    try
    {
        vector<Gene> LocalGenesVector;
        transform(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), back_inserter(LocalGenesVector), [](const auto& Gene){ return Gene.second; });
        sort(LocalGenesVector.begin(), LocalGenesVector.end(), [](const Gene& G1, const Gene& G2){ return G1.StartPosInGenome < G2.StartPosInGenome; });

        for (const auto& Gene : LocalGenesVector)
        {
            UnsignedInt StartPos = Gene.StartPosInGenome;
            UnsignedInt EndPos = Gene.EndPosInGenome;
            LoggersManagerObject.Log(STREAM("Gene = " << Gene.NumId << " " << StartPos << " " << EndPos));
        }

        UnsignedInt NumberOfFoundPromoterSequences = 0;
        string AttachPolymeraseToDNAStartSequenceStr = "TATAAT";

        UnsignedInt NumberOfFoundPromoters = 0;
        UnsignedInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 5;
        const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };
        const SignedInt MaxDistanceToStart = 100;

        auto AttachPolymeraseToDNAStartSequence = CellEngineUseful::ConvertStringSequenceToChainIdSequence(AttachPolymeraseToDNAStartSequenceStr);

        for (UnsignedInt StartGenomeIndex = 0; StartGenomeIndex < Genomes[0].size() - AttachPolymeraseToDNAStartSequenceStr.size(); StartGenomeIndex++)
        {
            auto NucleotidesSequenceToCompareVector = CellEngineUseful::ConvertStringSequenceToChainIdSequence(GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()));

            bool FoundSequenceNotFit = false;
            CellEngineUseful::CompareSequences(AttachPolymeraseToDNAStartSequence, NucleotidesSequenceToCompareVector, FoundSequenceNotFit, true);
            if (FoundSequenceNotFit == false)
            {
                if (SwitchLogsBool == true)
                    CellEngineUseful::SwitchOffLogs();

                LoggersManagerObject.Log(STREAM("Promoter sequence found = " << GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()) << " StartGenomeIndex = " << StartGenomeIndex));

                 for (UnsignedInt GeneIndex = 1; GeneIndex < LocalGenesVector.size(); GeneIndex++)
                     if (StartGenomeIndex < LocalGenesVector[GeneIndex].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length() && StartGenomeIndex > LocalGenesVector[GeneIndex - 1].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length())
                     {
                         LoggersManagerObject.Log(STREAM("Promoter sequence for gene = " << GeneIndex << " StartGenomeIndex for sequence " << AttachPolymeraseToDNAStartSequenceStr << " = " << StartGenomeIndex << " Gene.StartPosInGenome " << LocalGenesVector[GeneIndex].StartPosInGenome << " Diff = " << LocalGenesVector[GeneIndex].StartPosInGenome - StartGenomeIndex));
                         break;
                     }

                if (SwitchLogsBool == true)
                    CellEngineUseful::SwitchOnLogs();

                for (UnsignedInt PosInGenomeToFindCodonStart = StartGenomeIndex + 6; PosInGenomeToFindCodonStart < StartGenomeIndex + 6 + MaxDistanceToStart && PosInGenomeToFindCodonStart < Genomes[0].size() - 3; ++PosInGenomeToFindCodonStart)
                    FindCodons(GenomesLines[0], StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, 1, 1, SwitchLogsBool);

                NumberOfFoundPromoterSequences++;
            }
        }

        LoggersManagerObject.Log(STREAM("Number Of Found Promoter Sequences = " << NumberOfFoundPromoterSequences << endl));

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
        LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
    }
    CATCH("checking genome promoters")
}

bool IsPotentialHairpin(const std::string& SequenceStr)
{
    const UnsignedInt Length = SequenceStr.size();
    SignedInt Matches = 0;

    for (SignedInt Pos = 0; Pos < Length / 2; ++Pos)
    {
        const char Left = SequenceStr[Pos];
        const char Right = SequenceStr[Length - 1 - Pos];

        if ((Left == 'G' && Right == 'C') || (Left == 'C' && Right == 'G') || (Left == 'A' && Right == 'T') || (Left == 'T' && Right == 'A'))
            Matches++;
    }

    return Matches >= Length / 2;
}

void FindTranscriptionTerminators1(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOffLogs();

        const SignedInt HairpinLength = 8;
        const SignedInt MinimumLengthOfADownstreamOfHairpin = 4;
        const SignedInt MaximumLengthOfADownstreamOfHairpin = 8;

        UnsignedInt NumberOfFoundTerminators = 0;
        map<GeneIdInt, UnsignedInt> FoundTerminatorsCounter;

        for (UnsignedInt PosInGenome = 0; PosInGenome < GenomeStr.size() - HairpinLength - MaximumLengthOfADownstreamOfHairpin; ++PosInGenome)
        {
            std::string HairpinCandidate = GenomeStr.substr(PosInGenome, HairpinLength);

            if (IsPotentialHairpin(HairpinCandidate))
            {
                bool FoundAs = false;
                for (UnsignedInt PosInGenomeHairpin = PosInGenome + HairpinLength; PosInGenomeHairpin < PosInGenome + HairpinLength + MaximumLengthOfADownstreamOfHairpin; ++PosInGenomeHairpin)
                {
                    SignedInt LengthOfRunning = 0;
                    while (GenomeStr[PosInGenomeHairpin + LengthOfRunning] == 'A' && LengthOfRunning < MaximumLengthOfADownstreamOfHairpin)
                        LengthOfRunning++;

                    if (LengthOfRunning >= MinimumLengthOfADownstreamOfHairpin)
                    {
                        LoggersManagerObject.Log(STREAM("Potential terminator found at position: " << PosInGenome << endl));
                        LoggersManagerObject.Log(STREAM("Hairpin sequence: " << HairpinCandidate << endl));
                        LoggersManagerObject.Log(STREAM("Run of A's at position: " << PosInGenomeHairpin << " (length: " << LengthOfRunning << ")" << endl));

                        FoundTerminatorsCounter[PosInGenome]++;
                        NumberOfFoundTerminators++;
                        FoundAs = true;

                        break;
                    }
                }

                if (FoundAs)
                    PosInGenome += HairpinLength + MaximumLengthOfADownstreamOfHairpin;
            }
        }

        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("Number Of Found Terminators = " << NumberOfFoundTerminators));
        LoggersManagerObject.Log(STREAM("Number Of Found Different Terminators = " << FoundTerminatorsCounter.size()));
    }
    CATCH("finding transcription terminators 1")
}

SignedInt CalculateHairpinScore1(const std::string& LeftStem, const std::string& RightStem)
{
    SignedInt Score = 0;
    const UnsignedInt Length = LeftStem.size();

    for (SignedInt Pos = 0; Pos < Length; ++Pos)
    {
        const char Left = LeftStem[Pos];
        const char Right = RightStem[Length - 1 - Pos];

        if ((Left == 'G' && Right == 'C') || (Left == 'C' && Right == 'G'))
            Score += 3;
        else
        if ((Left == 'A' && Right == 'T') || (Left == 'T' && Right == 'A'))
            Score += 1;
        else
            Score -= 2;
    }

    return Score;
}

void FindTerminatorsForGenes1(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOffLogs();

        std::vector<pair<GeneIdInt, pair<UnsignedInt, UnsignedInt>>> GenesStartsAndEnds;
        for (const auto& Gene : ParticlesKindsManagerObject.Genes)
            GenesStartsAndEnds.emplace_back(Gene.second.NumId, make_pair(Gene.second.StartPosInGenome, Gene.second.EndPosInGenome));

        const SignedInt MinHairpinLength = 8;
        const SignedInt MaxHairpinLength = 15;
        const SignedInt MinimumLengthOfADownstreamOfHairpin = 4;
        const SignedInt MaximumLengthOfADownstreamOfHairpin = 10;
        const SignedInt MaxDistanceFromGeneEnd = 150;
        const SignedInt MaxDistanceBetweenHairpinAndAs = 40;
        const SignedInt MinHairpinScore = 5;

        SignedInt NumberOfFoundTerminators = 0;
        map<GeneIdInt, UnsignedInt> FoundTerminatorsCounter;

        for (UnsignedInt GeneIndex = 0; GeneIndex < GenesStartsAndEnds.size(); ++GeneIndex)
        {
            UnsignedInt GeneEnd = GenesStartsAndEnds[GeneIndex].second.second;
            UnsignedInt NextGeneStart = (GeneIndex < GenesStartsAndEnds.size() - 1) ? GenesStartsAndEnds[GeneIndex + 1].second.first : GenomeStr.size();

            bool TerminatorFound = false;

            UnsignedInt SearchRegionStart = GeneEnd + 1;
            UnsignedInt SearchRegionEnd = std::min(GeneEnd + MaxDistanceFromGeneEnd, NextGeneStart - 1);

            for (UnsignedInt PosInGenome = SearchRegionStart; PosInGenome < SearchRegionEnd; ++PosInGenome)
            {
                for (SignedInt HairpinLength = MinHairpinLength; HairpinLength <= MaxHairpinLength; ++HairpinLength)
                {
                    if (PosInGenome + HairpinLength * 2 >= GenomeStr.size())
                        continue;

                    std::string LeftStem = GenomeStr.substr(PosInGenome, HairpinLength);
                    std::string RightStem = GenomeStr.substr(PosInGenome + HairpinLength, HairpinLength);

                    SignedInt HairpinScore = CalculateHairpinScore1(LeftStem, RightStem);

                    if (HairpinScore >= MinHairpinScore)
                    {
                        for (UnsignedInt PosInGenomeHairpin = PosInGenome + HairpinLength * 2; PosInGenomeHairpin < PosInGenome + HairpinLength * 2 + MaxDistanceBetweenHairpinAndAs; ++PosInGenomeHairpin)
                        {
                            SignedInt LengthOfRunning = 0;
                            while (GenomeStr[PosInGenomeHairpin + LengthOfRunning] == 'A' && LengthOfRunning < MaximumLengthOfADownstreamOfHairpin)
                                LengthOfRunning++;

                            if (LengthOfRunning >= MinimumLengthOfADownstreamOfHairpin)
                            {
                                if (PosInGenomeHairpin + LengthOfRunning < NextGeneStart)
                                {
                                    LoggersManagerObject.Log(STREAM("Terminator found downstream of gene ending at position: " << GeneEnd << endl));
                                    LoggersManagerObject.Log(STREAM("Hairpin left stem: " << LeftStem << " right stem: " << RightStem << " (Score: " << HairpinScore << ")" << endl));
                                    LoggersManagerObject.Log(STREAM("Run of A's at position: " << PosInGenomeHairpin << " (length: " << LengthOfRunning << ")" << endl));

                                    FoundTerminatorsCounter[PosInGenome]++;
                                    NumberOfFoundTerminators++;
                                    TerminatorFound = true;
                                }
                                break;
                            }
                        }
                    }
                    if (TerminatorFound)
                        break;
                }
                if (TerminatorFound)
                    break;
            }
        }

        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("Number Of Found Terminators = " << NumberOfFoundTerminators));
        LoggersManagerObject.Log(STREAM("Number Of Found Different Terminators = " << FoundTerminatorsCounter.size()));
    }
    CATCH("finding termantors for genes 1")
}

SignedInt CalculateHairpinScore2(const std::string& LeftStem, const std::string& RightStem)
{
    SignedInt Score = 0;
    const SignedInt Length = LeftStem.size();

    for (SignedInt Index = 0; Index < Length; ++Index)
    {
        char Left = LeftStem[Index];
        char Right = RightStem[Length - 1 - Index];

        if ((Left == 'G' && Right == 'C') || (Left == 'C' && Right == 'G'))
            Score += 3;
        else
        if ((Left == 'A' && Right == 'T') || (Left == 'T' && Right == 'A'))
            Score += 1;
        else
            Score -= 1;
    }

    return Score;
}

void FindTerminatorsForGenes2(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    try
    {
        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOffLogs();

        std::vector<pair<GeneIdInt, UnsignedInt>> GenesEnds;
        for (const auto& Gene : ParticlesKindsManagerObject.Genes)
            GenesEnds.emplace_back(Gene.second.NumId, Gene.second.EndPosInGenome);

        const SignedInt MinHairpinLength = 8;
        const SignedInt MaxHairpinLength = 15;
        const SignedInt MinimumLengthOfADownstreamOfHairpin = 4;
        const SignedInt MaximumLengthOfADownstreamOfHairpin = 10;
        const SignedInt MaxDistanceFromGeneEnd = 150;
        const SignedInt MaxDistanceBetweenHairpinAndAs = 40;
        const SignedInt MinHairpinScore = 5;

        SignedInt NumberOfFoundTerminators = 0;
        map<GeneIdInt, UnsignedInt> FoundTerminatorsCounter;

        for (auto GeneEnd : GenesEnds)
        {
            bool TerminatorFound = false;

            UnsignedInt SearchRegionStart = GeneEnd.second + 1;
            UnsignedInt SearchRegionEnd = std::min(GeneEnd.second + MaxDistanceFromGeneEnd, GenomeStr.size() - 1);

            for (UnsignedInt PosInGenome = SearchRegionStart; PosInGenome < SearchRegionEnd; ++PosInGenome)
            {
                for (SignedInt HairpinLength = MinHairpinLength; HairpinLength <= MaxHairpinLength; ++HairpinLength)
                {
                    if (PosInGenome + HairpinLength * 2 >= GenomeStr.size())
                        continue;

                    std::string LeftStem = GenomeStr.substr(PosInGenome, HairpinLength);
                    std::string RightStem = GenomeStr.substr(PosInGenome + HairpinLength, HairpinLength);

                    SignedInt HairpinScore = CalculateHairpinScore2(LeftStem, RightStem);

                    if (HairpinScore >= MinHairpinScore)
                        for (UnsignedInt PosInGenomeHairpin = PosInGenome + HairpinLength * 2; PosInGenomeHairpin < PosInGenome + HairpinLength * 2 + MaxDistanceBetweenHairpinAndAs; ++PosInGenomeHairpin)
                        {
                            SignedInt LengthOfRunning = 0;
                            while (GenomeStr[PosInGenomeHairpin + LengthOfRunning] == 'A' && LengthOfRunning < MaximumLengthOfADownstreamOfHairpin)
                                LengthOfRunning++;

                            if (LengthOfRunning >= MinimumLengthOfADownstreamOfHairpin)
                            {
                                LoggersManagerObject.Log(STREAM("Terminator found downstream of gene ending at position: " << GeneEnd.second << endl));
                                LoggersManagerObject.Log(STREAM("Hairpin left stem: " << LeftStem << " right stem: " << RightStem << " (Score: " << HairpinScore << ")" << endl));
                                LoggersManagerObject.Log(STREAM("Run of A's at position: " << PosInGenomeHairpin << " (length: " << LengthOfRunning << ")" << endl));

                                FoundTerminatorsCounter[PosInGenome]++;
                                NumberOfFoundTerminators++;
                                TerminatorFound = true;
                                break;
                            }
                        }
                    if (TerminatorFound)
                        break;
                }
                if (TerminatorFound)
                    break;
            }
        }

        if (SwitchLogsBool == true)
            CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("Number Of Found Terminators = " << NumberOfFoundTerminators));
        LoggersManagerObject.Log(STREAM("Number Of Found Different Terminators = " << FoundTerminatorsCounter.size()));
    }
    CATCH("finding terminators for genes 2")
}

void TestSeveralDifferentKindsOfPromotersFindingsAndTerminatorFindingsAlgorithms(const std::vector<std::string>& GenomesLines, const std::vector<std::vector<UniqueIdInt>>& Genomes)
{
    try
    {
        FindPromoters(GenomesLines, Genomes, true);

        FindPromotersAndStartCodons1(GenomesLines[0], true);
        FindPromotersAndStartCodons2(GenomesLines[0], true);
        FindPromotersAndStartCodons3(GenomesLines[0], true);

        FindPromotersForGenesFromGeneStartPos(GenomesLines[0], true);

        FindTranscriptionTerminators1(GenomesLines[0], true);
        FindTerminatorsForGenes1(GenomesLines[0], true);
        FindTerminatorsForGenes2(GenomesLines[0], true);
    }
    CATCH("testing several different kinds of promoters finding algorithms")
}