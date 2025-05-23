
#include "Logger.h"
#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

#include "CellEngineParticle.h"
#include "CellEngineChemicalReactionsEngine.h"
#include "CellEngineChemicalReactionsManager.h"

using namespace std;

bool CellEngineChemicalReactionsEngine::TryToMakeRandomChemicalReaction(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

        vector<UnsignedInt> RandomParticlesTypes = GetRandomParticles(NumberOfReactants, MaxNumberOfReactants);

        sort(begin(RandomParticlesTypes), end(RandomParticlesTypes));

        if (auto IteratorUnique = unique(begin(RandomParticlesTypes), end(RandomParticlesTypes)); IteratorUnique == RandomParticlesTypes.end())
        {
            vector<string> ParticlesSymbolsForReactionToSort;
            ParticlesSymbolsForReactionToSort.reserve(10);
            for (auto& RandomParticleType : RandomParticlesTypes)
                ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).IdStr);
            sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());

            string ReactionSymbolsStr;
            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                ReactionSymbolsStr += (ParticleSymbolForReaction + "+");

            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            auto NumberOfElementsForKey = ChemicalReactionsManagerObject.ChemicalReactionsPosFromString.count(ReactionSymbolsStr);
            if (NumberOfElementsForKey > 0)
            {
                auto ReactionIter = ChemicalReactionsManagerObject.ChemicalReactionsPosFromString.equal_range(ReactionSymbolsStr).first;
                if (NumberOfElementsForKey > 1)
                {
                    uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(0, NumberOfElementsForKey - 1);

                    const UnsignedInt RandomShift = GetRandomValue<uniform_int_distribution, UnsignedInt>(UniformDistributionObjectSizeOfParticle_Uint64t);

                    LoggersManagerObject.Log(STREAM("Random Shift = [" << to_string(RandomShift) << "]" << endl));

                    for (UnsignedInt Step = 1; Step <= RandomShift; Step++)
                        ++ReactionIter;
                }
                auto& ReactionObject = ChemicalReactionsManagerObject.ChemicalReactions[ReactionIter->second];

                LoggersManagerObject.Log(STREAM("Reaction Position In Array = [" << to_string(ReactionIter->second) << "]" << endl));
                LoggersManagerObject.Log(STREAM("Reaction Id Num = [" << to_string(ReactionObject.ReactionIdNum) << "]" << endl));

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