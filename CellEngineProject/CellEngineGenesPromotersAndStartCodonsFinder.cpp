
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

void FindCodons(const string& GenomeStr, const vector<string>& StartCodonsList, map<GeneIdInt, UnsignedInt>& FoundGenesCounter, UniqueIdInt& NumberOfFoundPromoters, UniqueIdInt& NumberOfFoundPromotersConfirmed, const UnsignedInt ScopeForStartGeneIndex, const size_t PosInGenomeToFindCodonStart, const size_t PosInGenomeToFind35Box, const size_t PosInGenomeToFind10Box, const bool SwitchLogsBool)
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
    const string Promoter35BoxesList = "TTGACA";
    const string Promoter10BoxesList = "TATAA";
    const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

    const int MaxDistanceBetween35and10Boxes = 50;
    const int MaxDistanceFrom10BoxToStartCodon = 50;

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - Promoter35BoxesList.size(); ++PosInGenomeToFind35Box)
        if (GenomeStr.substr(PosInGenomeToFind35Box, Promoter35BoxesList.size()) == Promoter35BoxesList)
            for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + Promoter35BoxesList.size(); PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size(); ++PosInGenomeToFind10Box)
                if (GenomeStr.substr(PosInGenomeToFind10Box, Promoter10BoxesList.size()) == Promoter10BoxesList)
                {
                    for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + Promoter10BoxesList.size(); PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + Promoter10BoxesList.size() + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size(); ++PosInGenomeToFindCodonStart)
                        FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                    break;
                }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found ConfirmedPromoters = " << NumberOfFoundPromotersConfirmed << endl));
}

int CalculateDifference(const std::string& Sequence1, const std::string& Sequence2)
{
    int Mismatches = 0;
    for (size_t i = 0; i < Sequence1.size(); ++i)
    {
        if (Sequence1[i] != Sequence2[i])
        {
            Mismatches++;
        }
    }
    return Mismatches;
}

void FindPromotersAndStartCodons2(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
    const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

    const int MaxDistanceBetween35and10Boxes = 50;
    const int MaxDistanceFrom10BoxToStartCodon = 100;
    const int MaxNumberOfMismatchesForBox35 = 1;
    const int MaxNumberOfMismatchesForBox10 = 1;

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - 6; ++PosInGenomeToFind35Box)
        for (const string& Promoter35Box : Promoter35BoxesList)
            if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) <= MaxNumberOfMismatchesForBox35)
                for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + 6; PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                    for (const string& Promoter10Box : Promoter10BoxesList)
                        if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) <= MaxNumberOfMismatchesForBox10)
                        {
                            for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + 6; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + 6 + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                                FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                            break;
                        }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

int CalculateSimilarity(const string& Sequence1, const string& Sequence2)
{
    int Matches = 0;
    for (size_t i = 0; i < Sequence1.size(); ++i)
    {
        if (Sequence1[i] == Sequence2[i])
        {
            Matches++;
        }
    }
    return Matches;
}

void FindPromotersAndStartCodons3(const string& GenomeStr, const bool SwitchLogsBool)
{
    const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
    const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };

    const int MaxDistanceBetween35and10Boxes = 19;
    const int MaxDistanceFrom10BoxToStartCodon = 50;
    const int MinDistanceFrom10BoxToStartCodon = 10;
    const int MinimumSimilarityScoreFor35Box = 5;
    const int MinimumSimilarityScoreFor10Box = 5;

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 3;

    for (size_t PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - 6; ++PosInGenomeToFind35Box)
        for (const string& Promoter35Box : Promoter35BoxesList)
            if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
                for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + 16; PosInGenomeToFind10Box <= PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                    for (const string& Promoter10Box : Promoter10BoxesList)
                        if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                        {
                            for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + MinDistanceFrom10BoxToStartCodon; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                                FindCodons(GenomeStr, StartCodonsList, FoundGenesCounter, NumberOfFoundPromoters, NumberOfFoundPromotersConfirmed, ScopeForStartGeneIndex, PosInGenomeToFindCodonStart, PosInGenomeToFind35Box, PosInGenomeToFind10Box, SwitchLogsBool);
                            break;
                        }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

void FindPromotersForGenesFromGeneStartPos(const std::string& GenomeStr, const bool SwitchLogsBool)
{
    if (SwitchLogsBool == true)
        CellEngineUseful::SwitchOffLogs();

    std::vector<size_t> GenesStarts;
    for (const auto& Gene : ParticlesKindsManagerObject.Genes)
        GenesStarts.emplace_back(Gene.second.StartPosInGenome);

    const vector<string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const vector<string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };

    const int MinDistanceFrom10BoxToStartCodon = 10;
    const int MaxDistanceFrom10BoxToStartCodon = 50;
    const int MaxDistanceBetween35and10Boxes = 19;
    const int MinimumSimilarityScoreFor35Box = 5;
    const int MinimumSimilarityScoreFor10Box = 5;

    UniqueIdInt NumberOfFoundPromoters = 0;
    map<GeneIdInt, UnsignedInt> FoundPromotersCounter;

    for (size_t GeneStartPos : GenesStarts)
    {
        bool PromoterFound = false;

        for (int PosInGenomeToFind10Box = GeneStartPos - MinDistanceFrom10BoxToStartCodon; PosInGenomeToFind10Box >= std::max(0, static_cast<int>(GeneStartPos - MaxDistanceFrom10BoxToStartCodon)); --PosInGenomeToFind10Box)
            for (const string& Promoter10Box : Promoter10BoxesList)
                if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                    for (int PosInGenomeToFind35Box = PosInGenomeToFind10Box - MaxDistanceBetween35and10Boxes; PosInGenomeToFind35Box >= std::max(0, PosInGenomeToFind10Box - MaxDistanceBetween35and10Boxes - 5); --PosInGenomeToFind35Box)
                        for (const string& Promoter35Box : Promoter35BoxesList)
                            if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
                            {
                                LoggersManagerObject.Log(STREAM("Promoter found for gene starting at " << GeneStartPos));
                                LoggersManagerObject.Log(STREAM("-35 box at position: " << PosInGenomeToFind35Box));
                                LoggersManagerObject.Log(STREAM("-10 box at position: " << PosInGenomeToFind10Box << endl));

                                FoundPromotersCounter[PosInGenomeToFind10Box]++;
                                NumberOfFoundPromoters++;
                                PromoterFound = true;

                                goto Outside;
                            }
        Outside:

        if (PromoterFound == false)
            LoggersManagerObject.Log(STREAM("No promoter found for gene starting at " << GeneStartPos << std::endl));
    }

    if (SwitchLogsBool == true)
        CellEngineUseful::SwitchOnLogs();

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundPromotersCounter.size() << endl));
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

        UniqueIdInt NumberOfFoundPromoters = 0;
        UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 5;
        const vector<string> StartCodonsList = { "ATG", "GTG", "TTG" };
        const int MaxDistanceToStart = 100;

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

                for (size_t PosInGenomeToFindCodonStart = StartGenomeIndex + 6; PosInGenomeToFindCodonStart < StartGenomeIndex + 6 + MaxDistanceToStart && PosInGenomeToFindCodonStart < Genomes[0].size() - 3; ++PosInGenomeToFindCodonStart)
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

void TestSeveralDifferentKindsOfPromotersFindingAlgorithms(const std::vector<std::string>& GenomesLines, const std::vector<std::vector<UniqueIdInt>>& Genomes)
{
    FindPromoters(GenomesLines, Genomes, true);

    FindPromotersAndStartCodons1(GenomesLines[0], true);
    FindPromotersAndStartCodons2(GenomesLines[0], true);
    FindPromotersAndStartCodons3(GenomesLines[0], true);

    FindPromotersForGenesFromGeneStartPos(GenomesLines[0], true);
}