
#include <cmath>
#include <regex>

#include "Logger.h"
#include "StringUtils.h"
#include "ExceptionsMacro.h"
#include "DestinationPlatform.h"

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineSimulationSpace.h"
#include "CellEngineIllinoisDataCreator.h"
#include "CellEngineChemicalReactionsManager.h"

using namespace std;

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

void CellEngineIllinoisDataCreator::ReadReactionsFromJSONFile(const string& FileName, const bool ReadFromFile)
{
    try
    {
        if (ReadFromFile == true)
        {
            boost::property_tree::ptree ReactionsPropertyTreeJSON;

            read_json(FileName, ReactionsPropertyTreeJSON);

            LoggersManagerObject.Log(STREAM("JSON NUMBER OF TYPES OF PARTICLES = " << ReactionsPropertyTreeJSON.get_child("metabolites").size()));
            UnsignedInt ParticleNumber = 1;
            for (const auto &ReactionsPropertyTreeJSONTreeElementParticle: ReactionsPropertyTreeJSON.get_child("metabolites"))
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
            for (const auto &ReactionsPropertyTreeJSONTreeElementReaction: ReactionsPropertyTreeJSON.get_child("reactions"))
            {
                LoggersManagerObject.Log(STREAM("REACTION ID = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("id")));
                LoggersManagerObject.Log(STREAM("REACTION NAME = " << ReactionNumber << " " << ReactionsPropertyTreeJSONTreeElementReaction.second.get<string>("name")));
                for (const auto &ReactionsPropertyTreeJSONTreeElementParticle: ReactionsPropertyTreeJSONTreeElementReaction.second.get_child("metabolites"))
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
    }
    CATCH("reading reactions from json file")
}

void CellEngineIllinoisDataCreator::GetParticlesFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement)
{
    try
    {
        LoggersManagerObject.Log(STREAM("XML NUMBER OF TYPES OF PARTICLES = " << ReactionsPropertyTreeXMLTreeElement.size()));
        for (const auto& ReactionsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement)
        {
            auto IdStr = ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id");
            ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, IdStr, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name"), ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula"), 0, 0, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment"), 0 });

            auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(IdStr);
            for (auto& ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

            LoggersManagerObject.Log(STREAM("PARTICLE = " << IdStr << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:charge") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula")));
            LoggersManagerObject.Log(STREAM(""));

            ParticleKindId++;
        }
    }
    CATCH("getting particles from reactions xml file")
}

void CellEngineIllinoisDataCreator::GetNewParticlesFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement)
{
    try
    {
        UnsignedInt NumberOfNewParticles = 0;

        LoggersManagerObject.Log(STREAM("XML NUMBER OF TYPES OF PARTICLES = " << ReactionsPropertyTreeXMLTreeElement.size()));
        for (const auto& ReactionsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement)
            if (auto IdStr = ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.id"); ParticlesKindsManagerObject.GetParticleKindFromStrId(IdStr).has_value() == false)
            {
                ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, IdStr, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name"), ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula"), 0, 0, ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment"), 0 });

                auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(IdStr);
                for (auto& ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                    ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

                LoggersManagerObject.Log(STREAM("NEW PARTICLE = " << IdStr << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.name") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.compartment") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:charge") << " " << ReactionsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:chemicalFormula")));
                LoggersManagerObject.Log(STREAM(""));

                NumberOfNewParticles++;
                ParticleKindId++;
            }
            else
                LoggersManagerObject.Log(STREAM("ALREADY EXISTS = " << IdStr));

        LoggersManagerObject.Log(STREAM("Number of new particles = " << NumberOfNewParticles));
    }
    CATCH("getting particles from reactions xml file")
}

void CellEngineIllinoisDataCreator::GetProteinsFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement)
{
    try
    {
        LoggersManagerObject.Log(STREAM("XML NUMBER OF GENE PRODUCTS PROTEINS = " << ReactionsPropertyTreeXMLTreeElement.size()));
        for (const auto& ProteinsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement)
        {
            GeneIdInt GeneId = stoi(ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id").substr(9, 4));
            string ProteinName = JCVISYN3APredStr + ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id").substr(9, 4);
            ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ProteinName, "", ProteinName, GeneId, 0, "c", 10 });

            auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(ProteinName);
            for (auto& ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

            LoggersManagerObject.Log(STREAM("PARTICLE PROTEIN = " << GeneId << " " << ProteinName << " " << ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).GeneId << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:name") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:label")));
            LoggersManagerObject.Log(STREAM(""));

            ParticleKindId++;
        }
    }
    CATCH("getting proteins from reactions xml")
}

void CellEngineIllinoisDataCreator::GetNewProteinsFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement)
{
    try
    {
        UnsignedInt NumberOfNewProteins = 0;

        LoggersManagerObject.Log(STREAM("XML NUMBER OF GENE PRODUCTS PROTEINS = " << ReactionsPropertyTreeXMLTreeElement.size()));
        for (const auto& ProteinsPropertyTreeXMLTreeElementParticle : ReactionsPropertyTreeXMLTreeElement)
        {
            GeneIdInt GeneId = stoi(ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id").substr(9, 4));

            if (string ProteinName = JCVISYN3APredStr + ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id").substr(9, 4); ParticlesKindsManagerObject.GetParticleKindFromStrId(ProteinName).has_value() == false)
            {
                ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ProteinName, "", ProteinName, GeneId, 0, "c", 10 });

                auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(ProteinName);
                for (auto& ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                    ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

                LoggersManagerObject.Log(STREAM("NEW PARTICLE PROTEIN = " << GeneId << " " << ProteinName << " " << ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).GeneId << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:id") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:name") << " " << ProteinsPropertyTreeXMLTreeElementParticle.second.get<string>("<xmlattr>.fbc:label")));
                LoggersManagerObject.Log(STREAM(""));

                NumberOfNewProteins++;
                ParticleKindId++;
            }
        }
        LoggersManagerObject.Log(STREAM("Number of new proteins = " << NumberOfNewProteins));
    }
    CATCH("getting proteins from reactions xml")
}

void GetDataForGeneProductsForReactionFromXMLFile(const boost::property_tree::ptree& PropertyTreeXMLTreeElementReactionGeneProduct, vector<ParticleKindForChemicalReaction>& LocalVectorOfElements, const string& Str0, const string& Str1, const string& Str2, const string& Str3)
{
    try
    {
        for (const auto &ReactionsPropertyTreeXMLTreeElementGeneProduct : PropertyTreeXMLTreeElementReactionGeneProduct.get_child(Str0))
        {
            SignedInt GeneId = stoi(ReactionsPropertyTreeXMLTreeElementGeneProduct.second.get<string>("<xmlattr>.fbc:geneProduct").substr(9, 4));
            if (auto ParticleKindResult = ParticlesKindsManagerObject.GetParticleKindFromGeneId(GeneId); ParticleKindResult.has_value() == true)
            {
                LocalVectorOfElements.emplace_back(ParticleKindResult->EntityId, 1, "", true);
                LoggersManagerObject.Log(STREAM(Str1 << JCVISYN3APredStr + to_string(GeneId) << " " << 1));
            }
            else
                LoggersManagerObject.Log(STREAM(Str2));

            LoggersManagerObject.Log(STREAM(Str3 << ReactionsPropertyTreeXMLTreeElementGeneProduct.second.get<string>("<xmlattr>.fbc:geneProduct")));
        }
    }
    CATCH("getting data for gene products from xml file")
}

void GetDataForReactionFromXMLFile(const boost::property_tree::ptree& PropertyTreeXMLTreeElementReaction, vector<ParticleKindForChemicalReaction>& LocalVectorOfElements, const string& Str0, const string& Str1, const string& Str2)
{
    try
    {
        for (const auto& ReactionsPropertyTreeXMLTreeElementProduct : PropertyTreeXMLTreeElementReaction.get_child(Str0))
        {
            auto ParticleKindIdStr = ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.species");
            UnsignedInt Counter = stoi(ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.stoichiometry"));

            if (auto LocalParticleKindObject = ParticlesKindsManagerObject.GetParticleKindFromStrId(ParticleKindIdStr); LocalParticleKindObject.has_value() == true)
            {
                 LocalVectorOfElements.emplace_back(LocalParticleKindObject->EntityId, Counter, "", true);
                 LoggersManagerObject.Log(STREAM(Str1 << ReactionsPropertyTreeXMLTreeElementProduct.second.get<string>("<xmlattr>.constant") << " " << ParticleKindIdStr << " " << Counter));
            }
            else
                 LoggersManagerObject.Log(STREAM(Str2 << " : " << ParticleKindIdStr));
        }
    }
    CATCH("getting data for reaction from xml file")
}

void CellEngineIllinoisDataCreator::AddXMLChemicalReaction(ChemicalReaction& ReactionObject)
{
    try
    {
        if (ChemicalReactionsManagerObject.ChemicalReactionsPosFromIdStr.contains(ReactionObject.ReactionIdStr) == false)
        {
            LoggersManagerObject.Log(STREAM("ADDED NEW TYPE OF REACTION ID_STR = " << ReactionObject.ReactionIdStr));

            ReactionObject.ReactionIdNum = ReactionId;
            ReactionObject.ReactantsStr = ChemicalReactionsManager::GetStringOfSortedParticlesDataNames(ReactionObject.Reactants);
            ReactionObject.SpecialReactionFunction = nullptr;
            ChemicalReactionsManagerObject.AddChemicalReaction(ReactionObject);

            if (ReactionObject.Reversible == true)
            {
                ReactionObject.ReactantsStr = ChemicalReactionsManager::GetStringOfSortedParticlesDataNames(ReactionObject.Products);
                swap(ReactionObject.Products, ReactionObject.Reactants);
                ChemicalReactionsManagerObject.AddChemicalReaction(ReactionObject);
            }

            ReactionId++;
        }
    }
    CATCH("adding xml chemical reaction")
}

void CellEngineIllinoisDataCreator::GetProperReactionsListFromXMLFile(const boost::property_tree::ptree& ReactionsPropertyTreeXMLTreeElement)
{
    try
    {
        LoggersManagerObject.Log(STREAM("XML NUMBER OF TYPES OF REACTIONS = " << ReactionsPropertyTreeXMLTreeElement.size()));
        for (const auto& ReactionsPropertyTreeXMLTreeElementReaction : ReactionsPropertyTreeXMLTreeElement)
        {
            ChemicalReaction ReactionObject;
            ReactionObject.ReactionIdStr = ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.id");
            ReactionObject.ReactionName = ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.name");
            ReactionObject.Reversible = ReactionsPropertyTreeXMLTreeElementReaction.second.get<bool>("<xmlattr>.reversible");
            ReactionObject.UpperFluxBound = ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.fbc:upperFluxBound");
            ReactionObject.LowerFluxBound = ReactionsPropertyTreeXMLTreeElementReaction.second.get<string>("<xmlattr>.fbc:lowerFluxBound");

            LoggersManagerObject.Log(STREAM("REACTION = " << ReactionObject.ReactionIdStr << " " << ReactionObject.Reversible << " " << ReactionObject.ReactionName << " " << ReactionObject.UpperFluxBound << " " << ReactionObject.LowerFluxBound));

            std::vector<ParticleKindForChemicalReaction> LocalReactantsOr;

            if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("listOfReactants"))
                GetDataForReactionFromXMLFile(ReactionsPropertyTreeXMLTreeElementReaction.second, ReactionObject.Reactants, "listOfReactants", "PRODUCT = ", "ERROR: PRODUCT in REACTION NOT FOUND");

            if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("listOfProducts"))
                GetDataForReactionFromXMLFile(ReactionsPropertyTreeXMLTreeElementReaction.second, ReactionObject.Products, "listOfProducts", "REACTANT = ", "ERROR: REACTANT in REACTION NOT FOUND");

            if (ReactionsPropertyTreeXMLTreeElementReaction.second.get_child_optional("fbc:geneProductAssociation"))
            {
                auto& ReactionsPropertyTreeXMLTreeElementReactionGeneProduct = ReactionsPropertyTreeXMLTreeElementReaction.second.get_child("fbc:geneProductAssociation");
                if (ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child_optional("fbc:or"))
                    GetDataForGeneProductsForReactionFromXMLFile(ReactionsPropertyTreeXMLTreeElementReactionGeneProduct, LocalReactantsOr, "fbc:or", "GENE PRODUCT OR = ", "ERROR: GENE PRODUCT OR in REACTION NOT FOUND", "GENE PRODUCT OR = ");
                else
                if (ReactionsPropertyTreeXMLTreeElementReactionGeneProduct.get_child_optional("fbc:and"))
                    GetDataForGeneProductsForReactionFromXMLFile(ReactionsPropertyTreeXMLTreeElementReactionGeneProduct, ReactionObject.Reactants, "fbc:and", "GENE PRODUCT AND = ", "ERROR: GENE PRODUCT AND in REACTION NOT FOUND", "GENE PRODUCT AND = ");
                else
                    GetDataForGeneProductsForReactionFromXMLFile(ReactionsPropertyTreeXMLTreeElementReaction.second, ReactionObject.Reactants, "fbc:geneProductAssociation", "GENE PRODUCT = ", "ERROR: GENE PRODUCT in REACTION NOT FOUND", "GENE PRODUCT ONE = ");
            }

            if (LocalReactantsOr.empty() == false)
                for (const auto& LocalReactantsOrObject : LocalReactantsOr)
                {
                    ReactionObject.Reactants.push_back(LocalReactantsOrObject);

                    AddXMLChemicalReaction(ReactionObject);

                    ReactionObject.Reactants.pop_back();
                }
            else
                AddXMLChemicalReaction(ReactionObject);

            LoggersManagerObject.Log(STREAM(""));
        }
        LoggersManagerObject.Log(STREAM(""));
    }
    CATCH("get proper reactions list from xml file")
}

void CellEngineIllinoisDataCreator::ReadReactionsFromXMLFile(const string& FileName)
{
    try
    {
        boost::property_tree::ptree ReactionsPropertyTreeXML;

        read_xml(FileName, ReactionsPropertyTreeXML, boost::property_tree::xml_parser::trim_whitespace);

        for (const auto& ReactionsPropertyTreeXMLTreeElement : ReactionsPropertyTreeXML.get_child("sbml").get_child("model"))
        {
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfSpecies")
                GetParticlesFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "fbc:listOfGeneProducts")
                GetProteinsFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfReactions")
                GetProperReactionsListFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
        }
    }
    CATCH("reading reactions from xml file")
}

void CellEngineIllinoisDataCreator::ReadNewReactionsFromXMLFile(const string& FileName)
{
    try
    {
        ParticleKindId = max_element(ParticlesKindsManagerObject.ParticlesKinds.begin(), ParticlesKindsManagerObject.ParticlesKinds.end(), [](const pair<EntityIdInt, ParticleKind>& lhs, const pair<EntityIdInt, ParticleKind>& rhs){ return lhs.first < rhs.first; })->first + 1;
        ReactionId = max_element(ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.begin(), ChemicalReactionsManagerObject.ChemicalReactionsPosFromId.end(), [](const pair<UnsignedInt, UnsignedInt>& lhs, const pair<UnsignedInt, UnsignedInt>& rhs){ return lhs.first < rhs.first; })->first + 1;

        boost::property_tree::ptree ReactionsPropertyTreeXML;

        read_xml(FileName, ReactionsPropertyTreeXML, boost::property_tree::xml_parser::trim_whitespace);

        for (const auto& ReactionsPropertyTreeXMLTreeElement : ReactionsPropertyTreeXML.get_child("sbml").get_child("model"))
        {
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfSpecies")
                GetNewParticlesFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "fbc:listOfGeneProducts")
                GetNewProteinsFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
            else
            if (ReactionsPropertyTreeXMLTreeElement.first == "listOfReactions")
                GetProperReactionsListFromXMLFile(ReactionsPropertyTreeXMLTreeElement.second);
        }
    }
    CATCH("reading new reactions from xml file")
}

void CellEngineIllinoisDataCreator::CheckHowManyParticleDataForGeneratorIsNotInParticleKindsAndAddThem(const bool UpdateParticleKinds)
{
    try
    {
        for (auto ParticleIterator = ParticlesDataForGenerator.begin(); ParticleIterator != ParticlesDataForGenerator.end(); ParticleIterator = ParticlesDataForGenerator.upper_bound(ParticleIterator->first))
            if (ParticlesKindsManagerObject.GetParticleKindFromStrId(ParticleIterator->first).has_value() == false)
                if (ParticlesKindsManagerObject.GetParticleKindFromGeneId(ParticleIterator->second.GeneId).has_value() == false)
                {
                    LoggersManagerObject.Log(STREAM("LACKING " << ParticleIterator->first << " EVEN BY GENE " << ParticleIterator->second.GeneId << " in ParticlesKinds"));
                    if (UpdateParticleKinds == true)
                    {
                        ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ParticleIterator->first, "ADDED IN SECOND ROUND", ParticleIterator->first, static_cast<UnsignedInt>(ParticleIterator->second.GeneId), 0, "c", 10 });
                        auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(ParticleIterator->first);
                        for (auto &ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                            ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

                        ParticleKindId++;
                    }
                }

        for (auto ParticleIterator = ParticlesDataForGenerator.begin(); ParticleIterator != ParticlesDataForGenerator.end(); ParticleIterator = ParticlesDataForGenerator.upper_bound(ParticleIterator->first))
            if (ParticleIterator->second.ParticleType == ParticlesTypes::mRNA)
                if (ParticlesKindsManagerObject.GetParticleKindFromStrId(ParticleIterator->first).has_value() == false)
                {
                    LoggersManagerObject.Log(STREAM("LACKING mRNA " << ParticleIterator->first << " EVEN BY GENE " << ParticleIterator->second.GeneId << " in ParticlesKinds"));
                    if (UpdateParticleKinds == true)
                    {
                        ParticlesKindsManagerObject.AddParticleKind({ ParticleKindId, ParticleIterator->first, "ADDED IN SECOND ROUND", ParticleIterator->first, static_cast<UnsignedInt>(ParticleIterator->second.GeneId), 0, "c", 10 });
                        auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(ParticleIterator->first);
                        for (auto &ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                            ParticlesKindsManagerObject.GetParticleKind(ParticleKindId).ParticleKindSpecialDataSector.emplace_back(ParticleDataForGeneratorIterator->second);

                        ParticleKindId++;
                    }
                }

        if (UpdateParticleKinds == true)
            for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
                if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                    if (ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType == ParticlesTypes::OtherProtein)
                        for (auto ParticleIterator = ParticlesDataForGenerator.begin(); ParticleIterator != ParticlesDataForGenerator.end(); ParticleIterator = ParticlesDataForGenerator.upper_bound(ParticleIterator->first))
                        {
                            auto ParticlesDataForGeneratorRange = ParticlesDataForGenerator.equal_range(ParticleIterator->first);
                            for (auto &ParticleDataForGeneratorIterator = ParticlesDataForGeneratorRange.first; ParticleDataForGeneratorIterator != ParticlesDataForGeneratorRange.second; ParticleDataForGeneratorIterator++)
                                if (ParticleDataForGeneratorIterator->second.GeneId == ParticleKindObject.second.GeneId && ParticleDataForGeneratorIterator->second.ParticleType != ParticlesTypes::OtherProtein && ParticleDataForGeneratorIterator->second.IsProtein == true)
                                {
                                    ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType = ParticleDataForGeneratorIterator->second.ParticleType;
                                    ParticleKindObject.second.Name = ParticleDataForGeneratorIterator->first;
                                    LoggersManagerObject.Log(STREAM("UPDATED TYPE = " << ParticleDataForGeneratorIterator->second.GeneId << " " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleDataForGeneratorIterator->second.ParticleType)));
                                }
                        }
    }
    CATCH("checking how many particle data for generator in not in particle kinds")
};

void CellEngineIllinoisDataCreator::CheckHowManyParticlesKindsHasCounterAtStartOfSimulationEquZeroAndAddThem(const bool UpdateParticleKinds)
{
    try
    {
        if (UpdateParticleKinds == true)
            for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
                if (ParticleKindObject.second.EntityId >= StartParticleKindId)
                    if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == true)
                    {
                        ParticleKindObject.second.ParticleKindSpecialDataSector.emplace_back(ParticleKindSpecialData{ -1, "", "", false, ParticlesTypes::Basic, false, 100 });
                        LoggersManagerObject.Log(STREAM("UPDATED COUNTER FOR PARTICLE = " << ParticleKindObject.second.EntityId << " " << ParticleKindObject.second.IdStr << " " << ParticleKindObject.second.Name << " " << ParticleKindObject.second.Formula << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType)));
                    }
    }
    CATCH("checking how many particles kinds has counter at start of simulation equ zero")
};

void CellEngineIllinoisDataCreator::ReadAndParseGenesFile(const string& FileName)
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
                GeneObject = { GeneId, MatchSave1[0], MatchSave1[1], MatchSave1[2], "", StartPos, EndPos };

                if (MatchSave1[2].starts_with("30S"))
                    ParticlesKindsManagerObject.Ribosomes30SProteinsList.emplace_back(GeneId);
                else
                if (MatchSave1[2].starts_with("50S"))
                    ParticlesKindsManagerObject.Ribosomes50SProteinsList.emplace_back(GeneId);
            }
            else
                SequenceStr += Line.substr(0, Line.length() - 1);
        }

        LoggersManagerObject.Log(STREAM("Num of genes = " << MainLineNumber));
    }
    CATCH("reading and parsing csv file")
};

void CellEngineIllinoisDataCreator::PrintGenesFile()
{
    try
    {
        for (const auto& GeneObject : ParticlesKindsManagerObject.Genes)
        {
            LoggersManagerObject.Log(STREAM("Gene NumId = " << GeneObject.second.NumId));
            LoggersManagerObject.Log(STREAM("Gene StrId = " << GeneObject.second.StrId));
            LoggersManagerObject.Log(STREAM("Gene Description = " << GeneObject.second.Description));
            LoggersManagerObject.Log(STREAM("Gene ProteinId = " << GeneObject.second.ProteinId));
            LoggersManagerObject.Log(STREAM("Gene Start Position in Genome= " << GeneObject.second.StartPosInGenome));
            LoggersManagerObject.Log(STREAM("Gene End Position in Genome = " << GeneObject.second.EndPosInGenome));
            LoggersManagerObject.Log(STREAM("Gene Sequence = " << GeneObject.second.Sequence));
            LoggersManagerObject.Log(STREAM("Gene END" << endl));
        }
        LoggersManagerObject.Log(STREAM("Gene Size = " << ParticlesKindsManagerObject.Genes.size()));
    }
    CATCH("printing gene file")
}

void CellEngineIllinoisDataCreator::RemapProteinsNames(const string& ParticlesDirectory)
{
    try
    {
        auto ParsedCSVFileStructure = ReadAndParseCSVFile(ParticlesDirectory + string("JCVI-syn3A quantitative proteomics_corected.csv"), ',');
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= 2 && Row <= 460)
                MappedNamesOfProteins.insert(make_pair(ParsedCSVFileStructure[Row][0], ParsedCSVFileStructure[Row][1]));
        LoggersManagerObject.Log(STREAM("MP SIZE = " << MappedNamesOfProteins.size()));

        vector<pair<string, ParticleKindSpecialData>> TempVectorForParticleContainerKeys;
        copy(ParticlesDataForGenerator.begin(), ParticlesDataForGenerator.end(), back_inserter(TempVectorForParticleContainerKeys));

        for (auto& ParticleObject : TempVectorForParticleContainerKeys)
        {
            auto ParticleObjectMapIterator = MappedNamesOfProteins.find(ParticleObject.first);
            if (ParticleObjectMapIterator != MappedNamesOfProteins.end())
            {
                auto ParticleContainerNode = ParticlesDataForGenerator.extract(ParticleObject.first);
                ParticleContainerNode.key() = ParticleObjectMapIterator->second;
                ParticlesDataForGenerator.insert(std::move(ParticleContainerNode));
            }
        }
    }
    CATCH("remapping protein names")
}

void CellEngineIllinoisDataCreator::GetRemappingNamesForProteins(const string& ParticlesDirectory)
{
    try
    {
        auto ParsedCSVFileStructure = ReadAndParseCSVFile(ParticlesDirectory + string("JCVI-syn3A quantitative proteomics_corrected.csv"), ',');
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= 2 && Row <= 460)
                MappedNamesOfProteins.insert(make_pair(ParsedCSVFileStructure[Row][0], ParsedCSVFileStructure[Row][1]));

        LoggersManagerObject.Log(STREAM("MP SIZE = " << MappedNamesOfProteins.size()));
    }
    CATCH("get remapping names for proteins")
}

void CellEngineIllinoisDataCreator::ParticlesDataFromParsedCSVStructure(const vector<vector<string>>& ParsedCSVFileStructure, const UnsignedInt StartRow, const UnsignedInt EndRow, const UnsignedInt NameCol, const SignedInt GeneCol, const SignedInt AddedParticleCol, const SignedInt CleanTranscriptionProductCol, const SignedInt CounterCol, const bool FromConcentration, bool IsProtein, const UnsignedInt CounterParam, const string& NamePrefix, const string& Description, const UnsignedInt s1, const UnsignedInt s2, const ParticlesTypes ParticleType)
{
    try
    {
        for (UnsignedInt Row = 0; Row <= ParsedCSVFileStructure.size(); Row++)
            if (Row >= StartRow && Row <= EndRow)
            {
                string Name;
                UnsignedInt CounterAtStartOfSimulation;
                SignedInt GeneId = -1;
                string AddedParticleStr;
                SignedInt CleanTranscriptionProduct = -1;
                if (GeneCol != -1)
                    GeneId = stoi(ParsedCSVFileStructure[Row][GeneCol].substr(s1, s2));
                if (AddedParticleCol != -1)
                    AddedParticleStr = ParsedCSVFileStructure[Row][AddedParticleCol];
                if (CleanTranscriptionProductCol != -1)
                    CleanTranscriptionProduct = stoi(ParsedCSVFileStructure[Row][CleanTranscriptionProductCol]);
                if(CounterCol != -1)
                    CounterAtStartOfSimulation = (FromConcentration == false ? stoi(ParsedCSVFileStructure[Row][CounterCol]) : UnsignedInt(stold(ParsedCSVFileStructure[Row][CounterCol]) * AvogardoConstant * CapacityOfCell));
                else
                    CounterAtStartOfSimulation = CounterParam;
                Name = ParsedCSVFileStructure[Row][NameCol];
                auto NameIter = MappedNamesOfProteins.find(Name);
                if (NameIter != MappedNamesOfProteins.end())
                {
                    Name = NameIter->second;
                    GeneId = stoi(Name.substr(10, 4));
                }
                ParticlesDataForGenerator.insert(make_pair(NamePrefix + Name, ParticleKindSpecialData{ GeneId, Description, AddedParticleStr, CleanTranscriptionProduct, ParticleType, IsProtein, CounterAtStartOfSimulation }));
            }
    }
    CATCH("getting particles data from parsed csv structure")
}

void CellEngineIllinoisDataCreator::PrintAllParticlesData()
{
    try
    {
        for (auto ParticleIterator = ParticlesDataForGenerator.begin(); ParticleIterator != ParticlesDataForGenerator.end(); ParticleIterator = ParticlesDataForGenerator.upper_bound(ParticleIterator->first))
        {
            LoggersManagerObject.Log(STREAM("P NAME = " << ParticleIterator->first));
            auto Range = ParticlesDataForGenerator.equal_range(ParticleIterator->first);
            for (auto& ParticleDataIterator = Range.first; ParticleDataIterator != Range.second; ParticleDataIterator++)
                LoggersManagerObject.Log(STREAM(string("P GENE = " + string(ParticleDataIterator->second.GeneId != -1 ? JCVISYN3APredStr + to_string(ParticleDataIterator->second.GeneId) : "NoGene") + " TYPE = " + ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleDataIterator->second.ParticleType) + " D = #" + ParticleDataIterator->second.Description + "# Added = #" + ParticleDataIterator->second.AddedParticle + "# CLEAN PRODUCT = #" + to_string(ParticleDataIterator->second.CleanProductOfTranscription) + "# COUNTER = " + to_string(ParticleDataIterator->second.CounterAtStartOfSimulation))));

            LoggersManagerObject.Log(STREAM(""));
        }
    }
    CATCH("printing all particles data")
}

void CellEngineIllinoisDataCreator::ReadCSVFiles(bool Read, const string& ParticlesDirectory)
{
    try
    {
        if (Read == true)
        {
            GetRemappingNamesForProteins(ParticlesDirectory);

            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("proteomics.csv"), ','), 2, 429, 0, -1, -1, -1, 21, false, true, 0,"", "protein_from_gene", 0, 0, ParticlesTypes::OtherProtein);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("proteomics_additional.csv"), ','), 0, 29, 0, -1, -1, -1, 2, false, true, 0, "", "protein_from_gene", 0, 0, ParticlesTypes::OtherProtein);

            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("Escher_metData.csv"), ','), 1, 240, 0, -1, -1, -1, 1, true, false, 0,"M_", "basic", 0, 0, ParticlesTypes::Basic);

            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("trna_metabolites_synthase.csv"), ','), 1, 29, 0, 2, 3, -1, -1, false, false, 100, "", "_uncharged_trna_", 7, 4, ParticlesTypes::tRNA_uncharged);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("trna_metabolites_synthase.csv"), ','), 1, 29, 1, 2, 3, -1, -1, false, false, 100, "", "_charged_trna_", 7, 4, ParticlesTypes::tRNA_charged);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("mrna_counts.csv"), ','), 1, 454, 1, 1, -1, -1, -1, true, false, 1, "mrna_", "", 10, 4, ParticlesTypes::mRNA);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("rrna_metabolites.csv"), ','), 1, 6, 0, 1, -1, -1, -1, false, false, 10, "rrna_", "", 10, 4, ParticlesTypes::rRNA);

            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("polymerase_rna_proteins.csv"), ','), 1, 4, 0, 1, -1, -1, -1, false, true, 10, "", "", 7, 4, ParticlesTypes::RNAPolymeraseProtein);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("ribosomes_proteins_metabolites.csv"), ','), 1, 48, 0, 1, -1, -1, -1, false, true, 10, "", "", 7, 4, ParticlesTypes::RibosomeProtein);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("membrane_proteins_metabolites.csv"), ','), 1, 94, 0, 1, -1, -1, -1, false, true, 10, "", "", 7, 4, ParticlesTypes::MembraneProtein);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("proteins_metabolites_frac.csv"), ','), 1, 15, 0, 1, -1, 2, -1, false, true, 10, "", "", 7, 4, ParticlesTypes::ProteinFrac);
        }
    }
    CATCH("reading tsv files")
}

void CellEngineIllinoisDataCreator::ReadTSVFiles(bool Read, const string& ParticlesDirectory)
{
    try
    {
        if (Read == true)
        {
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"), ' '), 8, 310, 0, -1, -1, -1, 5, true, false, 0,"M_", "basic", 0, 0, ParticlesTypes::Basic);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("lipid_NoH2O_balanced_model.tsv"), ' '), 8, 42, 0, -1, -1, -1, 5, true, false, 0,"M_", "basic", 0, 0, ParticlesTypes::Basic);
            ParticlesDataFromParsedCSVStructure(ReadAndParseCSVFile(ParticlesDirectory + string("transport_NoH2O_Zane-TB-DB.tsv"), ' '), 8, 120, 0, -1, -1, -1, 5, true, false, 0,"M_", "basic", 0, 0, ParticlesTypes::Basic);
        }
    }
    CATCH("reading tsv files")
}

void CellEngineIllinoisDataCreator::ReadChemicalReactionsFromFiles()
{
    try
    {
        string ReactionsDirectory = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP;

        ReadReactionsFromXMLFile(ReactionsDirectory + string("iMB155.xml"));
        ReadReactionsFromJSONFile(ReactionsDirectory + string("iMB155.json"), false);
    }
    CATCH("reading chemical reactions from file")
}

void CellEngineIllinoisDataCreator::ReadNewChemicalReactionsFromFiles()
{
    try
    {
        ReadNewReactionsFromXMLFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP + string("Syn3A_updated.xml"));
    }
    CATCH("reading new chemical reactions from file")
}

void CellEngineIllinoisDataCreator::ReadAllIllinoisDataFromFiles()
{
    try
    {
        ParticlesKindsManagerObject.Genes.clear();
        ParticlesKindsManagerObject.ParticlesKinds.clear();
        ChemicalReactionsManagerObject.ChemicalReactions.clear();

        EntityIdInt LocalDNAIdentifier = CellEngineConfigDataObject.DNAIdentifier;
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::DNANucleotide, LocalDNAIdentifier, "DNANucleotide", "DNANucleotide", "DNA", -1, 0, "c", 1);
        EntityIdInt LocalRNAIdentifier = CellEngineConfigDataObject.RNAIdentifier;
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::DNANucleotide, LocalRNAIdentifier, "RNANucleotide", "RNANucleotide", "RNA", -1, 0, "c", 1);

        ReadAndParseGenesFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENES.txt"));
        PrintGenesFile();

        string ParticlesDirectory = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("particles") + OS_DIR_SEP;

        ReadCSVFiles(true, ParticlesDirectory);
        PrintAllParticlesData();

        ReadTSVFiles(false, ParticlesDirectory + string("tsv") + OS_DIR_SEP);

        ReadChemicalReactionsFromFiles();

        CheckHowManyParticleDataForGeneratorIsNotInParticleKindsAndAddThem(true);

        CheckHowManyParticlesKindsHasCounterAtStartOfSimulationEquZeroAndAddThem(true);

        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::Ribosome, ParticleKindId, "particle_RIBOSOME", "Ribosome70S", "R70S", -1, 0, "c", 100);
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::DNAPolymerase, ParticleKindId, "particle_DNAPolymerase", "DNAPolymerase", "dnapol", -1, 0, "c", 10);
        ParticlesKindsManagerObject.AddSingleParticleKind(ParticlesTypes::RNAPolymerase, ParticleKindId, "particle_RNAPolymerase", "RNAPolymerase", "rnapol", -1, 0, "c", 400);

        ParticlesKindsManagerObject.PrintAllParticleKinds();

        CheckHowManyParticleDataForGeneratorIsNotInParticleKindsAndAddThem(false);

        LoggersManagerObject.Log(STREAM("ALL DATA READ FROM FILE"));
    }
    CATCH("reading all illinois data from file")
}

void CellEngineIllinoisDataCreator::ReadGenes()
{
    try
    {
        ParticlesKindsManagerObject.Genes.clear();
        ReadAndParseGenesFile(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENES.txt"));
        PrintGenesFile();
    }
    CATCH("reading all illinois data from file")
}