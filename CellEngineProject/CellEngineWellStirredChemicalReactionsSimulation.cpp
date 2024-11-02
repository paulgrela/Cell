
#include "ExceptionsMacro.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineWellStirredChemicalReactionsSimulation.h"

std::vector<UnsignedInt> CellEngineWellStirredChemicalReactionsSimulation::GetRandomParticles(const UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos)
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

bool CellEngineWellStirredChemicalReactionsSimulation::MakeChemicalReaction(ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        for (const auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter));
        for (const auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter));
        LoggersManagerObject.Log(STREAM(""));

        for (const auto& ReactionReactant : ReactionObject.Reactants)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter -=  ReactionReactant.Counter;

        for (const auto& ReactionProduct : ReactionObject.Products)
            ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter +=  ReactionProduct.Counter;

        for (const auto& ReactionReactant : ReactionObject.Reactants)
            LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter));
        for (const auto& ReactionProduct : ReactionObject.Products)
            LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.EntityId << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Formula << " " << ParticlesKindsManagerObject.ParticlesKinds[ReactionProduct.EntityId].Counter));
        LoggersManagerObject.Log(STREAM(""));
    }
    CATCH("making reaction")

    return true;
}

bool CellEngineWellStirredChemicalReactionsSimulation::IsChemicalReactionPossible(const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos)
{
    return all_of(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [this](const ParticleKindForChemicalReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsManagerObject.ParticlesKinds[ReactionReactant.EntityId].Counter; });
};

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToMakeRandomChemicalReaction(2, 2, { 0, 0, 0 });
}