
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Logger.h"
#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"
#include "DestinationPlatform.h"
#include "CellEngineParticle.h"
#include "CellEngineChemicalReactions.h"

using namespace std;

void CellEngineChemicalReactions::PreprocessChemicalReactions()
{
    try
    {
        for (auto& ReactionObject : Reactions)
            sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForReaction& PK1, ParticleKindForReaction& PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); } );
    }
    CATCH("executing constructor CellEngineChemicalReactions")
}

bool CellEngineChemicalReactions::TryToMakeRandomChemicalReaction(UnsignedInt NumberOfReactants)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

        vector<UnsignedInt> RandomParticlesTypes = GetRandomParticles(NumberOfReactants);

        sort(begin(RandomParticlesTypes), end(RandomParticlesTypes));
        auto IteratorUnique = unique(begin(RandomParticlesTypes), end(RandomParticlesTypes));
        if (IteratorUnique == RandomParticlesTypes.end())
        {
            vector<string> ParticlesSymbolsForReactionToSort;
            ParticlesSymbolsForReactionToSort.reserve(10);
            for (auto& RandomParticleType : RandomParticlesTypes)
                ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).Symbol);
            std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());

            string ReactionSymbolsStr;
            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                ReactionSymbolsStr += (ParticleSymbolForReaction + " + ");
            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            auto ReactionIter = ReactionsIdByString.find(ReactionSymbolsStr);
            if (ReactionIter != ReactionsIdByString.end())
            {
                auto& ReactionObject = Reactions[ReactionIter->second];

                bool IsPossible = IsChemicalReactionPossible(ReactionObject);
                if (IsPossible == true)
                {
                    LoggersManagerObject.Log(STREAM("REACTION POSSIBLE" << endl));

                    MakeChemicalReaction(ReactionObject);
                    return true;
                }
                else
                    LoggersManagerObject.Log(STREAM("Reaction impossible!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Reaction for particles does not exist!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Particles types for reaction are not unique!"));

        LoggersManagerObject.Log(STREAM("END OF REACTION" << endl));
    }
    CATCH("making random reaction")

    return false;
}

void ReadChemicalReactionsFromFile()
{
    try
    {
        namespace pt = boost::property_tree;

        boost::property_tree::ptree root;

        boost::property_tree::read_json(string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("reactions") + OS_DIR_SEP + string("iMB155.json"), root);

        LoggersManagerObject.Log(STREAM("REACTIONS READ FROM FILE"));
    }
    CATCH("reading chemical reactions from file")
}