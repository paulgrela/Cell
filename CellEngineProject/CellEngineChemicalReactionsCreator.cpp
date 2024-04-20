
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
                    LoggersManagerObject.Log(STREAM("PARTICLE = " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:charge") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula")));
                }
            }
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "fbc:listOfGeneProducts")
                LoggersManagerObject.Log(STREAM("XML NUMBER OF GENE PRODUCTS PROTEINS = " << ReactionsPropertyTreeXMLTreeElement.second.size()));
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

struct LocalParticle
{
    std::string Name;
    UnsignedInt Counter;
};

std::multimap<std::string, UnsignedInt> Container;

vector<vector<string>> ReadAndParseCSVFile(const string& FileName)
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
            while (std::getline(LineStream, Cell, ','))
                ParsedRow.push_back(Cell);
            ParsedCSVFileStructure.push_back(ParsedRow);
        }
    }
    CATCH("reading and parsing csv file")

    return ParsedCSVFileStructure;
};

void GetNumberOfParticlesFromParsedCSVStructure(const vector<vector<string>>& ParsedCSVFileStructure, const UnsignedInt StartRow, const UnsignedInt EndRow, const UnsignedInt NameCol, const UnsignedInt NumberCol)
{
    try
    {
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= StartRow && Row <= EndRow)
                Container.insert(make_pair(ParsedCSVFileStructure[Row][NameCol], stoi(ParsedCSVFileStructure[Row][NumberCol])));
    }
    CATCH("")
}

void PrintCountersForAllParticlesTypes()
{
    try
    {
        for (auto ParticleIterator = Container.begin(); ParticleIterator != Container.end(); ParticleIterator = Container.upper_bound(ParticleIterator->first))
        {
            string CounterStr;
            auto Range = Container.equal_range(ParticleIterator->first);
            for (auto& CounterIterator = Range.first; CounterIterator != Range.second; CounterIterator++)
                CounterStr += (to_string(CounterIterator->second) + " ");
            LoggersManagerObject.Log(STREAM(ParticleIterator->first << " " << CounterStr));
        }
    }
    CATCH("printing counters for all particles")
}

void CellEngineChemicalReactionsCreator::ReadChemicalReactionsFromFile()
{
    try
    {
        ReadReactionsFromXMLFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP + string("iMB155.xml"));
        ReadReactionsFromJSONFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP + string("iMB155.json"));

        GetNumberOfParticlesFromParsedCSVStructure(ReadAndParseCSVFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("particles") + OS_DIR_SEP + string("proteomics.csv")), 2, 430, 0, 21);

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
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", 0, 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 12, "Polymerase", "POL", 0, 0 });

        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", 0, 0 });

        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", "?", 0, 0 });

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
