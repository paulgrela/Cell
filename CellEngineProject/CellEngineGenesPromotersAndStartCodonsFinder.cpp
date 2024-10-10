
#include <string>
#include <sstream>

#include "CellEngineTypes.h"

#include "CellEngineParticlesKindsManager.h"
#include "CellEngineGenesPromotersAndGenesStartCodonsFinder.h"

using namespace std;

void FindPromotersAndStartCodons1(const std::string& GenomeStr)
{
    const std::string Promoter35BoxesList = "TTGACA";
    const std::string Promoter10BoxesList = "TATAA";
    const std::vector<std::string> StartCodonsList = { "ATG", "GTG", "TTG" };

    const int MaxDistanceBetween35and10Boxes = 50;
    const int MaxDistanceFrom10BoxToStartCodon = 50;

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - Promoter35BoxesList.size(); ++PosInGenomeToFind35Box)
    {
        if (GenomeStr.substr(PosInGenomeToFind35Box, Promoter35BoxesList.size()) == Promoter35BoxesList)
        {
            for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + Promoter35BoxesList.size(); PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size(); ++PosInGenomeToFind10Box)
            {
                if (GenomeStr.substr(PosInGenomeToFind10Box, Promoter10BoxesList.size()) == Promoter10BoxesList)
                {
                    for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + Promoter10BoxesList.size(); PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + Promoter10BoxesList.size() + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size(); ++PosInGenomeToFindCodonStart)
                    {
                        for (const std::string& StartCodon : StartCodonsList)
                        {
                            if (GenomeStr.substr(PosInGenomeToFindCodonStart, StartCodon.size()) == StartCodon)
                            {
                                // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                                // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                                auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&PosInGenomeToFindCodonStart, &ScopeForStartGeneIndex](const auto& G1) { return PosInGenomeToFindCodonStart - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= PosInGenomeToFindCodonStart + ScopeForStartGeneIndex; });
                                if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                {
                                    //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                    NumberOfFoundPromotersConfirmed++;
                                }
                                NumberOfFoundPromoters++;

                                break;
                            }
                        }
                    }
                    break;
                }
            }
        }
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

void FindPromotersAndStartCodons2(const std::string& GenomeStr)
{
    const std::vector<std::string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const std::vector<std::string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
    const std::vector<std::string> StartCodonsList = { "ATG", "GTG", "TTG" };

    const int MaxDistanceBetween35and10Boxes = 50;
    const int MaxDistanceFrom10BoxToStartCodon = 100;
    const int MaxNumberOfMismatchesForBox35 = 1;
    const int MaxNumberOfMismatchesForBox10 = 1;

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t PosInGenomeToFind35Box = 0; PosInGenomeToFind35Box < GenomeStr.size() - 6; ++PosInGenomeToFind35Box)
    {
        for (const std::string& Promoter35Box : Promoter35BoxesList)
        {
            if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) <= MaxNumberOfMismatchesForBox35)
            {
                for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + 6; PosInGenomeToFind10Box < PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                {
                    for (const std::string& Promoter10Box : Promoter10BoxesList)
                    {
                        if (CalculateDifference(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) <= MaxNumberOfMismatchesForBox10)
                        {
                            for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + 6; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + 6 + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                            {
                                for (const std::string& start_codon : StartCodonsList)
                                {
                                    if (GenomeStr.substr(PosInGenomeToFindCodonStart, start_codon.size()) == start_codon)
                                    {
                                        // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                                        // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                                        auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&PosInGenomeToFindCodonStart, &ScopeForStartGeneIndex](const auto& G1) { return PosInGenomeToFindCodonStart - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= PosInGenomeToFindCodonStart + ScopeForStartGeneIndex; });
                                        if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                        {
                                            //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                            FoundGenesCounter[FindResult->second.NumId]++;
                                            NumberOfFoundPromotersConfirmed++;
                                        }
                                        NumberOfFoundPromoters++;

                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

int CalculateSimilarity(const std::string& Sequence1, const std::string& Sequence2)
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

void FindPromotersAndStartCodons3(const std::string& GenomeStr)
{
    const std::vector<std::string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const std::vector<std::string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };
    const std::vector<std::string> StartCodonsList = { "ATG", "GTG", "TTG" };

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
    {
        for (const std::string& Promoter35Box : Promoter35BoxesList)
        {
            if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
            {
                for (size_t PosInGenomeToFind10Box = PosInGenomeToFind35Box + 16; PosInGenomeToFind10Box <= PosInGenomeToFind35Box + MaxDistanceBetween35and10Boxes && PosInGenomeToFind10Box < GenomeStr.size() - 6; ++PosInGenomeToFind10Box)
                {
                    for (const std::string& Promoter10Box : Promoter10BoxesList)
                    {
                        if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                        {
                            for (size_t PosInGenomeToFindCodonStart = PosInGenomeToFind10Box + MinDistanceFrom10BoxToStartCodon; PosInGenomeToFindCodonStart < PosInGenomeToFind10Box + MaxDistanceFrom10BoxToStartCodon && PosInGenomeToFindCodonStart < GenomeStr.size() - 3; ++PosInGenomeToFindCodonStart)
                            {
                                for (const std::string& StartCodon : StartCodonsList)
                                {
                                    if (GenomeStr.substr(PosInGenomeToFindCodonStart, 3) == StartCodon)
                                    {
                                        // std::cout << "Promoter found: -35 at " << i << " and -10 at " << j << std::endl;
                                        // std::cout << "Start codon (" << start_codon << ") found at: " << k << std::endl;
                                        // std::cout << "----\n";

                                        auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&PosInGenomeToFindCodonStart, &ScopeForStartGeneIndex](const auto& G1) { return PosInGenomeToFindCodonStart - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= PosInGenomeToFindCodonStart + ScopeForStartGeneIndex; });
                                        if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                        {
                                            //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                            FoundGenesCounter[FindResult->second.NumId]++;
                                            NumberOfFoundPromotersConfirmed++;
                                        }
                                        NumberOfFoundPromoters++;

                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

void FindPromotersForGenesFromGeneStartPos(const std::string& GenomeStr, const std::vector<size_t>& GenesStarts)
{
    const std::vector<std::string> Promoter35BoxesList = { "TTGACA", "TTTACA", "CTGACA", "TTGATA" };
    const std::vector<std::string> Promoter10BoxesList = { "TATAAT", "TATAAA", "TATGAT", "TATATT" };

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
        {
            for (const std::string& Promoter10Box : Promoter10BoxesList)
            {
                if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind10Box, 6), Promoter10Box) >= MinimumSimilarityScoreFor10Box)
                {
                    for (int PosInGenomeToFind35Box = PosInGenomeToFind10Box - MaxDistanceBetween35and10Boxes; PosInGenomeToFind35Box >= std::max(0, PosInGenomeToFind10Box - MaxDistanceBetween35and10Boxes - 5); --PosInGenomeToFind35Box)
                    {
                        for (const std::string& Promoter35Box : Promoter35BoxesList)
                        {
                            if (CalculateSimilarity(GenomeStr.substr(PosInGenomeToFind35Box, 6), Promoter35Box) >= MinimumSimilarityScoreFor35Box)
                            {
                                LoggersManagerObject.Log(STREAM("Promoter found for gene starting at " << GeneStartPos << std::endl));
                                LoggersManagerObject.Log(STREAM("-35 box at position: " << PosInGenomeToFind35Box << std::endl));
                                LoggersManagerObject.Log(STREAM("-10 box at position: " << PosInGenomeToFind10Box << std::endl));

                                FoundPromotersCounter[PosInGenomeToFind10Box]++;
                                NumberOfFoundPromoters++;
                                PromoterFound = true;
                                break;
                            }
                        }
                        if (PromoterFound)
                            break;
                    }
                }
                if (PromoterFound)
                    break;
            }
            if (PromoterFound)
                break;
        }

        if (!PromoterFound)
        {
            LoggersManagerObject.Log(STREAM("No promoter found for gene starting at " << GeneStartPos << std::endl));
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundPromotersCounter.size() << endl));
}


inline void CompareSequences1(const std::vector<ChainIdInt>& TemplateSequence, const std::vector<ChainIdInt>& NucleotidesSequenceToCompareVector, bool& FoundSequenceNotFit)
{
    if (NucleotidesSequenceToCompareVector.size() >= TemplateSequence.size())
    {
        //LoggersManagerObject.Log(STREAM("LOOP COMPARISON SIZE = " << std::to_string(NucleotidesSequenceToCompareVector.size()) << " " << std::to_string(TemplateSequence.size())));

        for (UnsignedInt NucleotideNum = 0; NucleotideNum < TemplateSequence.size(); NucleotideNum++)
            if (CellEngineUseful::CompareIUPACNucleotideCode(TemplateSequence[NucleotideNum], NucleotidesSequenceToCompareVector[NucleotideNum]) == false)
            {
                //LoggersManagerObject.Log(STREAM("LOOP COMPARISON BREAK = " << std::to_string(NucleotideNum) << "#"));

                FoundSequenceNotFit = true;
                break;
            }
    }
    else
        FoundSequenceNotFit = true;
}

void FindInterGenesSequencesFromGenes()
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

void FindPromoters(const std::vector<std::string>& GenomesLines, const std::vector<std::vector<UniqueIdInt>>& Genomes)
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
        const std::vector<std::string> start_codons = {"ATG", "GTG", "TTG"};
        const int max_distance_to_start = 100;

        auto AttachPolymeraseToDNAStartSequence = CellEngineUseful::ConvertStringSequenceToChainIdSequence(AttachPolymeraseToDNAStartSequenceStr);

        for (UnsignedInt StartGenomeIndex = 0; StartGenomeIndex < Genomes[0].size() - AttachPolymeraseToDNAStartSequenceStr.size(); StartGenomeIndex++)
        {
            auto NucleotidesSequenceToCompareVector = CellEngineUseful::ConvertStringSequenceToChainIdSequence(GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()));

            bool FoundSequenceNotFit = false;
            CompareSequences1(AttachPolymeraseToDNAStartSequence, NucleotidesSequenceToCompareVector, FoundSequenceNotFit);
            if (FoundSequenceNotFit == false)
            {
                //LoggersManagerObject.Log(STREAM("Promoter sequence found = " << GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()) << " StartGenomeIndex = " << StartGenomeIndex));

                // for (UnsignedInt GeneIndex = 1; GeneIndex < LocalGenesVector.size(); GeneIndex++)
                //     if (StartGenomeIndex < LocalGenesVector[GeneIndex].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length() && StartGenomeIndex > LocalGenesVector[GeneIndex - 1].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length())
                //     {
                //         LoggersManagerObject.Log(STREAM("Promoter sequence for gene = " << GeneIndex << " StartGenomeIndex for sequence " << AttachPolymeraseToDNAStartSequenceStr << " = " << StartGenomeIndex << " Gene.StartPosInGenome " << LocalGenesVector[GeneIndex].StartPosInGenome << " Diff = " << LocalGenesVector[GeneIndex].StartPosInGenome - StartGenomeIndex));
                //         break;
                //     }

                for (size_t k = StartGenomeIndex + 6; k < StartGenomeIndex + 6 + max_distance_to_start && k < Genomes[0].size() - 3; ++k)
                {
                    for (const std::string& start_codon : start_codons)
                    {
                        if (GenomesLines[0].substr(k, 3) == start_codon)
                        //if (genome.substr(k, start_codon.size()) == start_codon)
                        {
                            // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                            // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                            //auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k](const auto& G1) { return G1.second.StartPosInGenome == k - 1 || G1.second.StartPosInGenome == k || G1.second.StartPosInGenome == k + 1; });
                            auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k, &ScopeForStartGeneIndex](const auto& G1) { return k - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= k + ScopeForStartGeneIndex; });
                            if (FindResult != ParticlesKindsManagerObject.Genes.end())
                            {
                                //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                FoundGenesCounter[FindResult->second.NumId]++;
                                NumberOfFoundPromotersConfirmed++;
                            }
                            NumberOfFoundPromoters++;

                            break;
                        }
                    }
                }

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

void TestSeveralDifferentKindsOfPromotersFindingAlgorithms(const std::vector<std::string>& GenomesLines)
{
    FindPromotersAndStartCodons1(GenomesLines[0]);
    FindPromotersAndStartCodons2(GenomesLines[0]);
    FindPromotersAndStartCodons3(GenomesLines[0]);

    std::vector<size_t> GenesStarts;
    for (const auto& Gene : ParticlesKindsManagerObject.Genes)
        GenesStarts.emplace_back(Gene.second.StartPosInGenome);
    FindPromotersForGenesFromGeneStartPos(GenomesLines[0], GenesStarts);
}


