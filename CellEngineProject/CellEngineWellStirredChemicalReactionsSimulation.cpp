
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

bool CellEngineWellStirredChemicalReactionsSimulation::MakeChemicalReaction(Reaction& ReactionObject)
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

    return true;
}

bool CellEngineWellStirredChemicalReactionsSimulation::IsChemicalReactionPossible(const Reaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter; });
};

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToMakeRandomChemicalReaction(2);
}