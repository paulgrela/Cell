
#ifndef CELL_ENGINE_GENES_PROMOTERS_AND_GENES_START_CODONS_FINDER_H
#define CELL_ENGINE_GENES_PROMOTERS_AND_GENES_START_CODONS_FINDER_H

void FindPromoters(const std::vector<std::string>& GenomesLines, std::vector<std::vector<UniqueIdInt>>& Genomes, bool SwitchOffLogsBool);

void FindPromotersAndStartCodons1(const std::string& GenomeStr, bool SwitchLogsBool);
void FindPromotersAndStartCodons2(const std::string& GenomeStr, bool SwitchLogsBool);
void FindPromotersAndStartCodons3(const std::string& GenomeStr, bool SwitchLogsBool);
void FindPromotersForGenesFromGeneStartPos(const std::string& GenomeStr, bool SwitchLogsBoolBool);

void FindInterGenesSequencesFromGenesData();

void TestSeveralDifferentKindsOfPromotersFindingAlgorithms(const std::vector<std::string>& GenomesLines, const std::vector<std::vector<UniqueIdInt>>& Genomes);

#endif
