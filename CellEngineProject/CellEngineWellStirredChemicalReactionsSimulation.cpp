
#include "ExceptionsMacro.h"
#include "CellEngineWellStirredChemicalReactionsSimulation.h"

std::vector<UnsignedInt> CellEngineWellStirredChemicalReactionsSimulation::GetRandomParticles(const UnsignedInt NumberOfReactants)
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

void CellEngineWellStirredChemicalReactionsSimulation::MakeReaction(Reaction& ReactionObject)
{
    try
    {
        for (auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Symbol << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter));
        for (auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Symbol << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter));
        LoggersManagerObject.Log(STREAM(""));

        for (auto& ReactionReactant : ReactionObject.Reactants)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter -=  ReactionReactant.Counter;

        for (auto& ReactionProduct : ReactionObject.Products)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter +=  ReactionProduct.Counter;

        for (auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Symbol << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter));
        for (auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Symbol << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter));
        LoggersManagerObject.Log(STREAM(""));
    }
    CATCH("making reaction")
}

bool CellEngineWellStirredChemicalReactionsSimulation::IsReactionPossible(const Reaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter; });
};

//bool CellEngineWellStirredChemicalReactionsSimulation::TryToMakeRandomReaction(const UnsignedInt NumberOfReactants)
//{
//    try
//    {
//        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));
//
//        vector<UnsignedInt> RandomParticlesTypes = GetRandomParticles(NumberOfReactants);
//
//        sort(begin(RandomParticlesTypes), end(RandomParticlesTypes));
//        auto IteratorUnique = unique(begin(RandomParticlesTypes), end(RandomParticlesTypes));
//        if (IteratorUnique == RandomParticlesTypes.end())
//        {
//            vector<string> ParticlesSymbolsForReactionToSort;
//            ParticlesSymbolsForReactionToSort.reserve(10);
//            for (auto& RandomParticleType : RandomParticlesTypes)
//                ParticlesSymbolsForReactionToSort.emplace_back(ParticlesKindsManagerObject.GetParticleKind(RandomParticleType).Symbol);
//            std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());
//
//            string ReactionSymbolsStr;
//            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
//                ReactionSymbolsStr += (ParticleSymbolForReaction + " + ");
//            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));
//
//            auto ReactionIter = ReactionsIdByString.find(ReactionSymbolsStr);
//            if (ReactionIter != ReactionsIdByString.end())
//            {
//                auto& ReactionObject = Reactions[ReactionIter->second];
//
//                bool IsPossible = IsReactionPossible(ReactionObject);
//                if (IsPossible == true)
//                {
//                    MakeReaction(ReactionObject);
//                    return true;
//                }
//                else
//                    LoggersManagerObject.Log(STREAM("Particles types are the same!"));
//            }
//            else
//                LoggersManagerObject.Log(STREAM("Reaction for particles does not exist!"));
//        }
//        else
//            LoggersManagerObject.Log(STREAM("Particles types for reaction are not unique!"));
//
//        LoggersManagerObject.Log(STREAM("END OF REACTION" << endl));
//    }
//    CATCH("making random reaction")
//
//    return false;
//}

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToMakeRandomReaction(2);
}