
#include "Logger.h"
#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

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
                ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).Formula);
            std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());

            string ReactionSymbolsStr;
            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                ReactionSymbolsStr += (ParticleSymbolForReaction + " + ");

            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            auto NumberOfElementsForKey = ReactionsPosFromString.count(ReactionSymbolsStr);
            if (NumberOfElementsForKey > 0)
            {
                auto ReactionIter = ReactionsPosFromString.equal_range(ReactionSymbolsStr).first;
                if (NumberOfElementsForKey > 1)
                {
                    uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(0, NumberOfElementsForKey - 1);
                    UnsignedInt RandomShift = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

                    LoggersManagerObject.Log(STREAM("Random Shift = [" << to_string(RandomShift) << "]" << endl));

                    for (UnsignedInt Step = 1; Step <= RandomShift; Step++)
                        ++ReactionIter;
                }
                auto& ReactionObject = Reactions[ReactionIter->second];

                LoggersManagerObject.Log(STREAM("Reaction Position In Array = [" << to_string(ReactionIter->second) << "]" << endl));
                LoggersManagerObject.Log(STREAM("Reaction Id = [" << to_string(ReactionObject.Id) << "]" << endl));

                bool IsPossible = IsChemicalReactionPossible(ReactionObject);
                if (IsPossible == true)
                {
                    LoggersManagerObject.Log(STREAM("REACTION POSSIBLE" << endl));

                    return MakeChemicalReaction(ReactionObject);
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