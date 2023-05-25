
#include "ExceptionsMacro.h"
#include "CellEngineWellStirredChemicalReactionsSimulation.h"

void CellEngineWellStirredChemicalReactionsSimulation::AddParticleKind(const ParticleKind& ParticleParam)
{
    Particles.emplace_back(ParticleParam);
}

void CellEngineWellStirredChemicalReactionsSimulation::AddReaction(const Reaction& ReactionParam)
{
    Reactions.emplace_back(ReactionParam);
    ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
}

std::vector<UnsignedInt> CellEngineWellStirredChemicalReactionsSimulation::GetRandomParticles(const UnsignedInt NumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(1, Particles.size());

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(UniformDistributionObjectUint64t(mt64X));
            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
            RandomParticlesTypes.back()--;
            if (Particles[RandomParticlesTypes.back()].Counter == 0)
            {
                LoggersManagerObject.Log(STREAM("Number of random reactant is zero!"));
                return {};
            }
        }
    }
    CATCH("getting random particles")

    return RandomParticlesTypes;
}

void CellEngineWellStirredChemicalReactionsSimulation::TryToDoRandomReaction(const UnsignedInt NumberOfReactants)
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
                ParticlesSymbolsForReactionToSort.emplace_back(Particles[RandomParticleType].Symbol);
            std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());

            string ReactionSymbolsStr;
            for (auto& ParticleSymbolForReaction : ParticlesSymbolsForReactionToSort)
                ReactionSymbolsStr += (ParticleSymbolForReaction + " + ");
            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            auto ReactionIter = ReactionsIdByString.find(ReactionSymbolsStr);
            if (ReactionIter != ReactionsIdByString.end())
            {
                Reaction* ReactionObject = &Reactions[ReactionIter->second];

                bool IsPossible = all_of(ReactionObject->Reactants.begin(), ReactionObject->Reactants.end(), [this](const ParticleKind& ReactionReactant){ return ReactionReactant.Counter <= Particles[ReactionReactant.EntityId].Counter; });

                if (IsPossible == true)
                {
                    for (auto& ReactionReactant : ReactionObject->Reactants)
                        LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.EntityId << " " << Particles[ReactionReactant.EntityId].Symbol << " " << Particles[ReactionReactant.EntityId].Counter));
                    for (auto& ReactionProduct : ReactionObject->Products)
                        LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.EntityId << " " << Particles[ReactionProduct.EntityId].Symbol << " " << Particles[ReactionProduct.EntityId].Counter));
                    LoggersManagerObject.Log(STREAM(""));

                    for (auto& ReactionReactant : ReactionObject->Reactants)
                        Particles[ReactionReactant.EntityId].Counter -=  ReactionReactant.Counter;

                    for (auto& ReactionProduct : ReactionObject->Products)
                        Particles[ReactionProduct.EntityId].Counter +=  ReactionProduct.Counter;

                    for (auto& ReactionReactant : ReactionObject->Reactants)
                        LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.EntityId << " " << Particles[ReactionReactant.EntityId].Symbol << " " << Particles[ReactionReactant.EntityId].Counter));
                    for (auto& ReactionProduct : ReactionObject->Products)
                        LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.EntityId << " " << Particles[ReactionProduct.EntityId].Symbol << " " << Particles[ReactionProduct.EntityId].Counter));
                    LoggersManagerObject.Log(STREAM(""));
                }
                else
                    LoggersManagerObject.Log(STREAM("Particles types are the same!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Reaction for particles does not exist!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Particles types for reaction are not unique!"));

        LoggersManagerObject.Log(STREAM("END OF REACTION" << endl));
    }
    CATCH("making random reaction")
}

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToDoRandomReaction(2);
}