
#ifndef CELL_ENGINE_ILLINOIS_DATA_CREATOR_H
#define CELL_ENGINE_ILLINOIS_DATA_CREATOR_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "CellEngineConstants.h"
#include "CellEngineParticleKind.h"
#include "CellEngineChemicalReactionsEngine.h"

class CellEngineIllinoisDataCreator
{
private:
    EntityIdInt ParticleKindId = StartParticleKindId;
    EntityIdInt ReactionId = StartReactionId;
private:
    std::map<std::string, std::string> MappedNamesOfProteins;
    std::multimap<std::string, ParticleKindSpecialData> ParticlesDataForGenerator;
public:
    void ReadAllIllinoisDataFromFiles();
public:
    static void ReadReactionsFromJSONFile(const std::string& FileName, bool ReadFromFile);
    void ReadReactionsFromXMLFile(const std::string& FileName);
    void GetParticlesFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement);
    void GetProteinsFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement);
    void GetProperReactionsListFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement);
    void AddXMLChemicalReaction(ChemicalReaction& ReactionObject);
public:
    static void AddSingleParticleKind(ParticlesTypes ParticlesTypesObject, EntityIdInt& ParticleKindIdParam, const std::string& IdStrParam, const std::string& NameParam, const std::string& FormulaParam, SignedInt GeneIdParam, ElectricChargeType ElectricChargeParam, const std::string& CompartmentParam, UnsignedInt CounterParam);
public:
    static void PrintGenesFile();
    void PrintAllParticlesData();
    static void PrintAllParticleKinds();
public:
    void CheckHowManyParticleDataForGeneratorIsNotInParticleKindsAndAddThem(bool UpdateParticleKinds);
    static void CheckHowManyParticlesKindsHasCounterAtStartOfSimulationEquZeroAndAddThem(bool UpdateParticleKinds);
    static void ReadAndParseGenesFile(const std::string& FileName);
    void RemapProteinsNames(const std::string& ParticlesDirectory);
    void GetRemappingNamesForProteins(const std::string& ParticlesDirectory);
    void ParticlesDataFromParsedCSVStructure(const std::vector<std::vector<std::string>>& ParsedCSVFileStructure, UnsignedInt StartRow, UnsignedInt EndRow, UnsignedInt NameCol, SignedInt GeneCol, SignedInt AddedParticleCol, SignedInt CleanTranscriptionProductCol, SignedInt CounterCol, bool FromConcentration, bool IsProtein, UnsignedInt CounterParam, const std::string& NamePrefix, const std::string& Description, UnsignedInt s1, UnsignedInt s2, ParticlesTypes ParticleType);
public:
    void ReadCSVFiles(bool Read, const std::string& ParticlesDirectory);
    void ReadTSVFiles(bool Read, const std::string& ParticlesDirectory);
};

#endif
