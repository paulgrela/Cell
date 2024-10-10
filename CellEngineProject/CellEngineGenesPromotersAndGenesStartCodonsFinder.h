
#ifndef CELL_ENGINE_GENES_PROMOTERS_AND_GENES_START_CODONS_FINDER_H
#define CELL_ENGINE_GENES_PROMOTERS_AND_GENES_START_CODONS_FINDER_H

void FindPromotersAndStartCodons1(const std::string& GenomeStr);
void FindPromotersAndStartCodons2(const std::string& GenomeStr);
void FindPromotersAndStartCodons3(const std::string& GenomeStr);
void FindPromotersForGenesFromGeneStartPos(const std::string& GenomeStr, const std::vector<size_t>& GenesStarts);

void FindInterGenesSequencesFromGenes();
void FindPromoters(const std::vector<std::string>& GenomesLines, std::vector<std::vector<UniqueIdInt>>& Genomes);

void TestSeveralDifferentKindsOfPromotersFindingAlgorithms(const std::vector<std::string>& GenomesLines);

#endif
