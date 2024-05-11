
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_CREATOR_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_CREATOR_H

#include "CellEngineChemicalReactions.h"

class CellEngineChemicalReactionsCreator : virtual public CellEngineChemicalReactions
{
private:
//    EntityIdInt ParticleKindId = 100000;
//    EntityIdInt ReactionId = 10000;
private:
    //std::map<std::string, std::string> MappedNamesOfProteins;
    //std::multimap<std::string, ParticleKindSpecialData> ParticlesDataForGenerator;
public:
    static void AddParticlesKinds();
    void ReadChemicalReactionsFromFile();
    void AddChemicalReactions();
public:
    static void ReadReactionsFromJSONFile(const std::string& FileName);
    void ReadReactionsFromXMLFile(const std::string& FileName);
public:
    static void PrintGenesFile();
    void PrintAllParticlesData();
    static void PrintAllParticleKinds();
public:
    void CheckHowManyParticleDataForGeneratorIsNotInParticleKindsAndAddThem(bool UpdateParticleKinds);
    static void ReadAndParseGenesFile(const std::string& FileName);
    void RemapProteinsNames(const std::string& ParticlesDirectory);
    void GetRemappingNamesForProteins(const std::string& ParticlesDirectory);
    void ParticlesDataFromParsedCSVStructure(const std::vector<std::vector<std::string>>& ParsedCSVFileStructure, UnsignedInt StartRow, UnsignedInt EndRow, UnsignedInt NameCol, SignedInt GeneCol, SignedInt AddedParticleCol, SignedInt CleanTranscriptionProductCol, SignedInt CounterCol, bool FromConcentration, bool IsProtein, UnsignedInt CounterParam, const std::string& NamePrefix, const std::string& Description, UnsignedInt s1, UnsignedInt s2, ParticlesTypes ParticleType);
public:
    void ReadCSVFiles(bool Read, const std::string& ParticlesDirectory);
    void ReadTSVFiles(bool Read, const std::string& ParticlesDirectory);
};

#endif
