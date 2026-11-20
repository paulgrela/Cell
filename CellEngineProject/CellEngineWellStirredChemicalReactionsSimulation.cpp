
#include "ExceptionsMacro.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineChemicalReactionsManager.h"
#include "CellEngineWellStirredChemicalReactionsSimulation.h"

std::vector<UnsignedInt> CellEngineWellStirredChemicalReactionsSimulation::GetRandomParticles(const UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(1, ParticlesKindsManagerObject.ParticlesKinds.size());

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(UniformDistributionObjectUint64t(mt64X));
            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
            RandomParticlesTypes.back()--;
            if (ParticlesKindsManagerObject.ParticlesKinds[RandomParticlesTypes.back()].Counter == 0)
            {
                LoggersManagerObject.Log(STREAM("Number of random reactant is zero!"));
                return {};
            }
        }
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

bool CellEngineWellStirredChemicalReactionsSimulation::MakeChemicalReaction(ChemicalReaction& ReactionObject)
{
    try
    {
        for (const auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.ParticleKindId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Counter));
        for (const auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.ParticleKindId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.ParticleKindId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.ParticleKindId].Counter));
        LoggersManagerObject.Log(STREAM(""));

        for (const auto& ReactionReactant : ReactionObject.Reactants)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Counter -=  ReactionReactant.Counter;

        for (const auto& ReactionProduct : ReactionObject.Products)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.ParticleKindId].Counter +=  ReactionProduct.Counter;

        for (const auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.ParticleKindId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Counter));
        for (const auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.ParticleKindId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.ParticleKindId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.ParticleKindId].Counter));
        LoggersManagerObject.Log(STREAM(""));
    }
    CATCH("making reaction")

    return true;
}

bool CellEngineWellStirredChemicalReactionsSimulation::IsChemicalReactionPossible(const ChemicalReaction& ReactionObject)
{
    return ranges::all_of(ReactionObject.Reactants, [](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.ParticleKindId].Counter; });
};

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToMakeRandomChemicalReaction(2, 2);
}

bool CellEngineWellStirredChemicalReactionsSimulation::TryToMakeRandomChemicalReaction(const UnsignedInt NumberOfReactants, const UnsignedInt MaxNumberOfReactants)
{
    try
    {
        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

        vector<UnsignedInt> RandomParticlesTypes = GetRandomParticles(NumberOfReactants, MaxNumberOfReactants);

        ranges::sort(RandomParticlesTypes);

        if (const auto IteratorUnique = ranges::unique(RandomParticlesTypes).begin(); IteratorUnique == RandomParticlesTypes.end())
        {
            vector<string> ParticlesSymbolsForReactionToSort;
            ParticlesSymbolsForReactionToSort.reserve(10);
            for (const auto& RandomParticleType : RandomParticlesTypes)
                ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).IdStr);
            ranges::sort(ParticlesSymbolsForReactionToSort);

            string ReactionSymbolsStr;
            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                ReactionSymbolsStr += (ParticleSymbolForReaction + "+");

            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            const auto NumberOfElementsForKey = ChemicalReactionsManagerObject.ChemicalReactionsPosFromString.count(ReactionSymbolsStr);
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
                LoggersManagerObject.Log(STREAM("Reaction Id Num = [" << to_string(ReactionObject.ChemicalReactionIdNum) << "]" << endl));

                const bool IsPossible = IsChemicalReactionPossible(ReactionObject);
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