
#include <cmath>
#include <regex>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Logger.h"
#include "StringUtils.h"
#include "ExceptionsMacro.h"
#include "DestinationPlatform.h"

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineChemicalReactionsCreator.h"

using namespace std;

std::map<std::string, std::string> MappedNamesOfProteins;
std::multimap<std::string, UnsignedInt> ParticleCountersContainer;

EntityIdInt ParticleKindId = 100000;

//wczytuje geny

//potem wczytuje czastki wszystkie

//potem trzeba je szukac po ID tekstowym - liniowo by dodac:
// ilosc
// identyfikator genu
// oznaczenie JAKIEGO JEST TYPU

// ksztalt - list of voxels



void ReadReactionsFromJSONFile(const string& FileName)
{
    try
    {
        boost::property_tree::ptree ReactionsPropertyTreeJSON;

        read_json(FileName, ReactionsPropertyTreeJSON);

        LoggersManagerObject.Log(STREAM("JSON NUMBER OF TYPES OF PARTICLES = " << ReactionsPropertyTreeJSON.get_child("metabolites").size()));
        UnsignedInt ParticleNumber = 1;
        for (const auto& ReactionsPropertyTreeJSONTreeElementParticle : ReactionsPropertyTreeJSON.get_child("metabolites"))
        {
            LoggersManagerObject.Log(STREAM("PARTICLE ID = " << ParticleNumber << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get<string>("id")));
            LoggersManagerObject.Log(STREAM("PARTICLE NAME = " << ParticleNumber << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get<string>("name")));
            LoggersManagerObject.Log(STREAM("PARTICLE FORMULA = " << ParticleNumber << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get<string>("formula")));
            LoggersManagerObject.Log(STREAM("PARTICLE CHARGE = " << ParticleNumber << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get<int>("charge")));
            LoggersManagerObject.Log(STREAM("PARTICLE COMPARTMENT = " << ParticleNumber << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get<string>("compartment")));
            LoggersManagerObject.Log(STREAM(""));
            ParticleNumber++;
        }

        LoggersManagerObject.Log(STREAM("JSON NUMBER OF TYPES OF REACTIONS = " << ReactionsPropertyTreeJSON.get_child("reactions").size()));
        UnsignedInt ReactionNumber = 1;
        for (const auto& ReactionsPropertyTreeJSONTreeElementReaction : ReactionsPropertyTreeJSON.get_child("reactions"))
        {
            LoggersManagerObject.Log(STREAM("REACTION ID = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("id")));
            LoggersManagerObject.Log(STREAM("REACTION NAME = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("name")));
            for (const auto& ReactionsPropertyTreeJSONTreeElementParticle : ReactionsPropertyTreeJSONTreeElementReaction.second.get_child("metabolites"))
            {
                LoggersManagerObject.Log(STREAM("PARTICLE = " << ReactionsPropertyTreeJSONTreeElementParticle.first << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get_value<string>()));
                LoggersManagerObject.Log(STREAM("PARTICLE = " << ReactionsPropertyTreeJSONTreeElementParticle.first << " " << ReactionsPropertyTreeJSONTreeElementParticle.second.get_value<double>()));
            }
            LoggersManagerObject.Log(STREAM("LOWER BOUND = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("lower_bound")));
            LoggersManagerObject.Log(STREAM("UPPER BOUND = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("upper_bound")));
            LoggersManagerObject.Log(STREAM("GENE REACTION RULE = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("gene_reaction_rule")));

            LoggersManagerObject.Log(STREAM(""));
            ReactionNumber++;
        }
    }
    CATCH("reading reactions from json file")
}

void ReadReactionsFromXMLFile(const string& FileName)
{
    try
    {
        boost::property_tree::ptree ReactionsPropertyTreeXML;

        read_xml(FileName, ReactionsPropertyTreeXML, boost::property_tree::xml_parser::trim_whitespace);

        for (const auto& ReactionsPropertyTreeXMLTreeElement : ReactionsPropertyTreeXML.get_child("sbml").get_child("model"))
        {
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfSpecies")
            {
                LoggersManagerObject.Log(STREAM("XML NUMBER OF TYPES OF PARTICLES = " << ReactionsPropertyTreeXMLTreeElement.second.size()));
                for (const auto& ReactionsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement.second)
                {
                    ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id"), ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name"), ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.formula"), 0, 0, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment"), 0 });
                    ParticleKindId++;
                    ParticleCountersContainer.insert(make_pair(ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id"), 10));
                    LoggersManagerObject.Log(STREAM("PARTICLE = " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:charge") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula")));
                }
            }
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "fbc:listOfGeneProducts")
            {
                LoggersManagerObject.Log(STREAM("XML NUMBER OF GENE PRODUCTS PROTEINS = " << ReactionsPropertyTreeXMLTreeElement.second.size()));
                for (const auto& ProteinsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement.second)
                {
                    string ProteinName = "JCVISYN3A_" + ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id").substr(9, 4);
                    ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ProteinName, ProteinName, ProteinName, 0, 0, "c", 0 });
                    ParticleKindId++;
                    ParticleCountersContainer.insert(make_pair(ProteinName, 10));
                    LoggersManagerObject.Log(STREAM("PARTICLE = " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:name") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:label")));
                }
            }
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfReactions")
            {
                LoggersManagerObject.Log(STREAM("XML NUMBER OF TYPES OF REACTIONS = " << ReactionsPropertyTreeXMLTreeElement.second.size()));
                for (const auto& ReactionsPropertyTreeXMLTreeElementReaction : ReactionsPropertyTreeXMLTreeElement.second)
                {
                    LoggersManagerObject.Log(STREAM("REACTION = " << ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.id") << " " << ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.reversible") << " " << ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.name") << " " << ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.fbc:upperFluxBound") << " " << ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.fbc:lowerFluxBound")));

                    if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("listOfReactants"))
                        for (const auto& ReactionsPropertyTreeXMLTreeElementReactant : ReactionsPropertyTreeXMLTreeElementReaction.second.get_child("listOfReactants"))
                            LoggersManagerObject.Log(STREAM("REACTANT = " << ReactionsPropertyTreeXMLTreeElementReactant.second.get<string>("<xmlattr>.constant") << " " << ReactionsPropertyTreeXMLTreeElementReactant.second.get<string>("<xmlattr>.species") << " " << ReactionsPropertyTreeXMLTreeElementReactant.second.get<string>("<xmlattr>.stoichiometry")));
                    if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("listOfProducts"))
                        for (const auto& ReactionsPropertyTreeXMLTreeElementProduct : ReactionsPropertyTreeXMLTreeElementReaction.second.get_child("listOfProducts"))
                            LoggersManagerObject.Log(STREAM("REACTANT = " << ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.constant") << " " << ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.species") << " " << ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.stoichiometry")));

                    if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("fbc:geneProductAssociation"))
                    {
                        auto& ReactionsPropertyTreeXMLTreeElementReactionGeneProduct = ReactionsPropertyTreeXMLTreeElementReaction.second.get_child("fbc:geneProductAssociation");
                        if (ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child_optional("fbc:or"))
                        {
                            for (const auto& ReactionsPropertyTreeXMLTreeElementGeneProduct : ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child("fbc:or"))
                                LoggersManagerObject.Log(STREAM("GENE PRODUCT OR = " << ReactionsPropertyTreeXMLTreeElementGeneProduct.second.get<string>("<xmlattr>.fbc:geneProduct")));
                        }
                        else
                        if (ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child_optional("fbc:and"))
                        {
                            for (const auto& ReactionsPropertyTreeXMLTreeElementGeneProduct : ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child("fbc:and"))
                                LoggersManagerObject.Log(STREAM("GENE PRODUCT AND = " << ReactionsPropertyTreeXMLTreeElementGeneProduct.second.get<string>("<xmlattr>.fbc:geneProduct")));
                        }
                        else
                            for (const auto& ReactionsPropertyTreeXMLTreeElementGeneProduct : ReactionsPropertyTreeXMLTreeElementReaction.second.get_child("fbc:geneProductAssociation"))
                                LoggersManagerObject.Log(STREAM("GENE PRODUCT ONE = " << ReactionsPropertyTreeXMLTreeElementGeneProduct.second.get<string>("<xmlattr>.fbc:geneProduct")));

                    }
                }
                LoggersManagerObject.Log(STREAM(""));
            }
        }
    }
    CATCH("reading reactions from xml file")
}

void ReadAndParseGenesFile(const string& FileName)
{
    smatch SMatchObject;
    regex GeneRegexObject(R"(=[\w. ()]+)");

    try
    {
        ifstream Data(FileName);
        string Line;
        Gene GeneObject;
        string SequenceStr;
        UnsignedInt MainLineNumber = 0;

        while (std::getline(Data, Line))
        {
            if (Line.substr(0, 4) == ">lcl")
            {
                MainLineNumber++;

                if (SequenceStr.empty() == false)
                {
                    GeneObject.Sequence = SequenceStr;
                    ParticlesKindsManagerObject.Genes[GeneObject.NumId] = GeneObject;
                    SequenceStr = "";
                }

                auto pos = Line.cbegin();
                auto end = Line.cend();
                vector<string> MatchSave1, MatchSave2;
                for ( ; regex_search(pos, end, SMatchObject, GeneRegexObject); pos = SMatchObject.suffix().first)
                    MatchSave1.emplace_back(SMatchObject[0].str().substr(1, SMatchObject[0].length() - 1));

                if (Line.find("gene=") == string::npos)
                    MatchSave1.insert(MatchSave1.begin() + 0, "");
                if (Line.find("pseudo=") != string::npos)
                    MatchSave1.insert(MatchSave1.begin() + 3, "");

                auto pos1 = MatchSave1[4].cbegin();
                auto end1 = MatchSave1[4].cend();
                regex GeneScopeRegexObject(R"((\d)+)");
                for ( ; regex_search(pos1, end1, SMatchObject, GeneScopeRegexObject); pos1 = SMatchObject.suffix().first)
                    MatchSave2.emplace_back(SMatchObject[0]);

                GeneIdInt GeneId = stoi(MatchSave1[1].substr(10, 4));
                UnsignedInt StartPos = stoi(MatchSave2[0]);
                UnsignedInt EndPos = stoi(MatchSave2[1]);
                GeneObject = { GeneId, MatchSave1[0], MatchSave1[1], MatchSave1[2], StartPos, EndPos, "" };
            }
            else
                SequenceStr += Line.substr(0, Line.length() - 1);
        }

        LoggersManagerObject.Log(STREAM("Num of genes = " << MainLineNumber));
    }
    CATCH("reading and parsing csv file")
};

void PrintGenesFile()
{
    try
    {
        for (const auto& GeneObject : ParticlesKindsManagerObject.Genes)
        {
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.NumId));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.StrId));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.Description));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.ProteinId));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.StartPosInGenome));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.EndPosInGenome));
            LoggersManagerObject.Log(STREAM("Gene = " << GeneObject.second.Sequence));
        }
        LoggersManagerObject.Log(STREAM("Gene Size = " << ParticlesKindsManagerObject.Genes.size()));
    }
    CATCH("printing gene file")
}




constexpr long double constexpr_power(long double num, unsigned int pow)
{
    return (pow == 0 ? 1 : num * constexpr_power(num, pow - 1));
}

constexpr long double CapacityOfCell = (4.0 / 3.0) * M_PI * constexpr_power((200 * 1E-09), 3);

constexpr long double AvogardoConstant = 6.022 * 1E+23;

vector<vector<string>> ReadAndParseCSVFile(const string& FileName, const char Separator)
{
    vector<vector<string>> ParsedCSVFileStructure;

    try
    {
        ifstream Data(FileName);
        string Line;
        while (std::getline(Data, Line))
        {
            stringstream LineStream(Line);
            string Cell;
            vector<std::string> ParsedRow;
            while (std::getline(LineStream, Cell, Separator))
                ParsedRow.push_back(Cell);
            ParsedCSVFileStructure.push_back(ParsedRow);
        }
    }
    CATCH("reading and parsing csv file")

    return ParsedCSVFileStructure;
};

void GetNumberOfParticlesFromParsedCSVStructure(const vector<vector<string>>& ParsedCSVFileStructure, const UnsignedInt StartRow, const UnsignedInt EndRow, const UnsignedInt NameCol, const UnsignedInt NumberCol, const bool FromConcentration, const string& NamePrefix, const UnsignedInt s1, const UnsignedInt s2)
{
    try
    {
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= StartRow && Row <= EndRow)
                ParticleCountersContainer.insert(make_pair(NamePrefix + (s1 != 0 ? ParsedCSVFileStructure[Row][NameCol].substr(s1, s2) : ParsedCSVFileStructure[Row][NameCol]), (FromConcentration == false ? stoi(ParsedCSVFileStructure[Row][NumberCol]) : UnsignedInt(stold(ParsedCSVFileStructure[Row][NumberCol]) * AvogardoConstant * CapacityOfCell))));
    }
    CATCH("getting number of particles from parsed csv structure")
}

void SetNumberOfParticlesForParticlesTypesFromParsedStructure(const vector<vector<string>>& ParsedCSVFileStructure, const UnsignedInt StartRow, const UnsignedInt EndRow, const UnsignedInt NameCol, UnsignedInt NumberOfParticles, const string& NamePrefix, const UnsignedInt s1, const UnsignedInt s2)
{
    try
    {
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= StartRow && Row <= EndRow)
                ParticleCountersContainer.insert(make_pair(NamePrefix + (s1 != 0 ? ParsedCSVFileStructure[Row][NameCol].substr(s1, s2) : ParsedCSVFileStructure[Row][NameCol]), NumberOfParticles));
    }
    CATCH("setting number of particles for particles types from parsed structure")
}

void PrintCountersForAllParticlesTypes()
{
    try
    {
        for (auto ParticleIterator = ParticleCountersContainer.begin(); ParticleIterator != ParticleCountersContainer.end(); ParticleIterator = ParticleCountersContainer.upper_bound(ParticleIterator->first))
        {
            string CounterStr;

            auto Range = ParticleCountersContainer.equal_range(ParticleIterator->first);
            for (auto& CounterIterator = Range.first; CounterIterator != Range.second; CounterIterator++)
                CounterStr += (to_string(CounterIterator->second) + " ");

            LoggersManagerObject.Log(STREAM(ParticleIterator->first << " " << CounterStr));
        }
    }
    CATCH("printing counters for all particles")
}

void RemapProteinsNames(const string& ParticlesDirectory)
{
    try
    {
        auto ParsedCSVFileStructure = ReadAndParseCSVFile(ParticlesDirectory + string("JCVI-syn3A quantitative proteomics.csv"), ',');
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= 2 && Row <= 430)
                MappedNamesOfProteins.insert(make_pair(ParsedCSVFileStructure[Row][0], ParsedCSVFileStructure[Row][1]));
        LoggersManagerObject.Log(STREAM("MP SIZE = " << MappedNamesOfProteins.size()));

        vector<pair<string, UnsignedInt>> TempVectorForParticleContainerKeys;
        copy(ParticleCountersContainer.begin(), ParticleCountersContainer.end(), back_inserter(TempVectorForParticleContainerKeys));

        for (auto& ParticleObject : TempVectorForParticleContainerKeys)
        {
            auto ParticleObjectMapIterator = MappedNamesOfProteins.find(ParticleObject.first);
            if (ParticleObjectMapIterator != MappedNamesOfProteins.end())
            {
                auto ParticleContainerNode = ParticleCountersContainer.extract(ParticleObject.first);
                ParticleContainerNode.key() = ParticleObjectMapIterator->second;
                ParticleCountersContainer.insert(std::move(ParticleContainerNode));
            }
        }
    }
    CATCH("remapping protein names")
}

void ReadCSVFiles(bool Read, const string& ParticlesDirectory)
{
    try
    {
        if (Read == true)
        {
            GetNumberOfParticlesFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("proteomics.csv"), ','), 2, 429, 0, 21, false, "", 0, 0);

            RemapProteinsNames(ParticlesDirectory);

            GetNumberOfParticlesFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("mRNA_Counts.csv"), ','), 1, 455, 1, 2, true, "mrna_", 0, 0); //mRNA generowac JCVISYN3A_0932 zamienic na MMSYN1_0932
            GetNumberOfParticlesFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("Escher_metData.csv"), ','), 1, 240, 0, 1, true, "M_", 0, 0); //male czastki - ile

            ReadAndParseCSVFile(ParticlesDirectory + string("ribo_protein_metabolites.csv"), ','); //pierwsza kolumna nazwa bialka nalezacego do rybosomow = TYLKO DO LISTY GENERATORA RYOSOMOW
            ReadAndParseCSVFile(ParticlesDirectory + string("membrane_protein_metabolites.csv"), ','); //lista bialek blonowych - TYLKO DO GENERATORA BIALEK RYBOSOMOW
            ReadAndParseCSVFile(ParticlesDirectory + string("trna_metabolites_synthase.csv"), ','); //lista trna z bialkami do ktorych sie wiaza - Z NIEGO POBRAC GEN CO GO KODUJE DO GENERATORA

            GetNumberOfParticlesFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("protein_metabolites_frac.csv"), ','), 1, 15, 1, 3, true, "JCVISYN3A_", 7, 4); //lista bialek w formie aktywnej i nieaktywnej - ale u mnie to niewazne bo kazde jako osobna czastka
            //wazna tylko informacja z jakiego genu pochodzi dana czastka
        }
    }
    CATCH("reading tsv files")
}

void ReadTSVFiles(bool Read, const string& ParticlesDirectory)
{
    try
    {
        if (Read == true)
        {
            SetNumberOfParticlesForParticlesTypesFromParsedStructure(ReadAndParseCSVFile(ParticlesDirectory + string("Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"), ' '), 8, 310, 0, 10, "", 0, 0);
            SetNumberOfParticlesForParticlesTypesFromParsedStructure(ReadAndParseCSVFile(ParticlesDirectory + string("lipid_NoH2O_balanced_model.tsv"), ' '), 8, 42, 0, 10, "", 0, 0);
            SetNumberOfParticlesForParticlesTypesFromParsedStructure(ReadAndParseCSVFile(ParticlesDirectory + string("transport_NoH2O_Zane-TB-DB.tsv"), ' '), 8, 120, 0, 10, "", 0, 0);
        }
    }
    CATCH("reading tsv files")
}

void CellEngineChemicalReactionsCreator::ReadChemicalReactionsFromFile()
{
    try
    {
        ReadAndParseGenesFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENES.txt"));
        PrintGenesFile();

        string ReactionsDirectory = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP;

        ReadReactionsFromXMLFile(ReactionsDirectory + string("iMB155.xml"));
        ReadReactionsFromJSONFile(ReactionsDirectory + string("iMB155.json"));

        string ParticlesDirectory = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("particles") + OS_DIR_SEP;

        ReadCSVFiles(true, ParticlesDirectory);
        ReadTSVFiles(false, ParticlesDirectory + string("tsv") + OS_DIR_SEP);

        PrintCountersForAllParticlesTypes();

        LoggersManagerObject.Log(STREAM("REACTIONS READ FROM FILE"));
        getchar();
    }
    CATCH("reading chemical reactions from file")
}

void CellEngineChemicalReactionsCreator::AddParticlesKinds()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", "H2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", "C6H12O6", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", "02", 0, 0, "c", 0  });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", "CO2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", "CH2CH2", 0, 0, "c", 0   });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", "CH3CH2(OH)", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", "CH3CHCH2", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", "HX", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", "CH3CHXCH3", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", "TEST", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", "CH2CH2O", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", "0", 0, 0, "c", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 12, "Polymerase", "POL", "POL", 0, 0, "c", 0 });

        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", "DNA", 0, 0, "c", 0 });

        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", "?", "?", 0, 0, "c", 0 });

        LoggersManagerObject.Log(STREAM("ADDED PARTICLES KINDS"));
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineChemicalReactionsCreator::AddChemicalReactions()
{
    try
    {
        const string DNASequenceForTestFindingDNA = "TACAAAAAAAGAGGTGTTAGC";

        const string DNASequence1ForTestCutLink1 = "TACAAAAAAAGAGGTGTT";
        const string DNASequence2ForTestCutLink1 = "AGCTCTTATTA";

        const string DNASequence1ForTestCutLink1Any = "ANY";

        const string DNASequence1ForTestCutLink2 = "TACAAAAAAAGAGGTGTT";

        const string DNASequenceForTestCutCrisper = "RNA";
        const string RNASequenceForTestCutCrisper = "ANY";

        AddChemicalReaction(Reaction(1101, "STD ONLY WITH SEQ", "CH3CH2(OH) + DNA + ", { { 5, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 10, 1, "", true } }, nullptr, "standard normal reaction example"));

        AddChemicalReaction(Reaction(10, "CUT 1 SEQ", "CH3CH2(OH) + DNA + ", { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "only in presence of chosen sequence of DNA"));

        AddChemicalReaction(Reaction(20, "LINK 1 SEQ", "CH3CH2(OH) + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink1, false }, { 10001, 1, DNASequence2ForTestCutLink1, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInChosenPlaceSpecialReactionFunction, "links one strand of DNA with 2 endings with chosen sequences"));
        AddChemicalReaction(Reaction(30, "LINK 1 ANY", "CH3CH2(OH) + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 2, DNASequence1ForTestCutLink1Any, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "link one strand of DNA with 2 endings with any sequence"));

        AddChemicalReaction(Reaction(40, "CUT 2 SEQ SHIFT 3 10", "CH3CHCH2 + DNA + ", 3, 10, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=3 shift=10"));
        AddChemicalReaction(Reaction(41, "CUT 2 SEQ SHIFT 7 3", "CH3CHCH2 + DNA + ", 7, 3, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=7 shift2=3"));
        AddChemicalReaction(Reaction(42, "CUT 2 SEQ SHIFT 10 3", "CH3CHCH2 + DNA + ", 10, 3, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineChemicalReactionsInSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=10 shift2=3"));

        AddChemicalReaction(Reaction(60, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGC", false }, { 10001, 1, "TCTTATT", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));
        AddChemicalReaction(Reaction(61, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGCTCTT", false }, { 10001, 1, "ATTATGA", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));

        AddChemicalReaction(Reaction(70, "LINK 2 ANY COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "ANY", false }, { 10001, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction, "Links both DNA strands with any sequence if opposite joining strands when complementary"));

        AddChemicalReaction(Reaction(80, "LINK 2 ANY EQU SAME", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "ANY", false }, { 10001, 1, "ANY", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "links both strands of DNA if they are cut equally in the same place"));

        AddChemicalReaction(Reaction(100, "CUT CRISPER 1", "CH3CHCH2 + RNA + DNA + ", 3, 7, { { 5, 1, "", false }, { 10001, 1, DNASequenceForTestCutCrisper, false }, { 10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of one strand of DNA with RNA as template"));
        AddChemicalReaction(Reaction(110, "CUT CRISPER 2", "CH3CHCH2 + RNA + DNA + ", 3, 7, { { 5, 1, "", false }, { 10001, 1, DNASequenceForTestCutCrisper, false }, { 10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of two strands of DNA with RNA as template"));

        AddChemicalReaction(Reaction(150, "POLYMERASE DNA START SEQ", "CH3CHCH2 + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGC", false } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseDNAStartSpecialReactionFunction, "links particle with DNA when found sequence and joins first nucleotide "));
        AddChemicalReaction(Reaction(160, "POLYMERASE DNA CONTINUE", "CH3CHCH2 + DNA + ", { { 5, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, {}, &CellEngineChemicalReactionsInSimulationSpace::PolymeraseDNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        AddChemicalReaction(Reaction(1001, "STD", "C6H12O6 + O2 + ", { { 1, 1, "", true }, { 2, 6, "", true } }, { { 3, 6, "", true }, { 0, 6, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1002, "STD", "CH2CH2 + H2O + ", { { 4, 1, "", true }, { 0, 1, "", true } }, { { 5, 1, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1003, "STD", "CH3CHCH2 + HX + ", { { 6, 1, "", true }, { 7, 1, "", true } }, { { 8, 1, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1004, "STD", "CH2CH2 + O + ", { { 4,  1, "", true }, { 11, 1, "", true } }, { { 10, 1, "", true } }, nullptr));

        PreprocessChemicalReactions();

        LoggersManagerObject.Log(STREAM("ADDED CHEMICAL REACTIONS"));
    }
    CATCH("adding particles kinds and reactions")
};
