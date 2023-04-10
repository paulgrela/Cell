
#include "CellEngineWellStirredChemicalReactionsSimulation.h"

void CellEngineWellStirredChemicalReactionsSimulation::TryToDoRandomReaction(const UnsignedIntType NumberOfReactants)
{
    LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

    std::uniform_int_distribution<uint64_t> UniformDistributionObject1X_Uint64t(1, Particles.size());

    vector<UnsignedIntType> RandomParticlesTypes;
    for (UnsignedIntType ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
    {
        RandomParticlesTypes.emplace_back(UniformDistributionObject1X_Uint64t(mt64X));
        LoggersManagerObject.Log(STREAM("Particle Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        RandomParticlesTypes.back()--;
        if (Particles[RandomParticlesTypes.back()].Counter == 0)
        {
            LoggersManagerObject.Log(STREAM("Number of random reactant is zero!"));
            return;
        }
    }

    //RandomParticlesTypes.clear(); RandomParticlesTypes = {9, 0};

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

        auto ReactionIter = Reactions.find(ReactionSymbolsStr);
        if (ReactionIter != Reactions.end())
        {
            bool IsPossible = true;

            for (auto& ReactionReactant : ReactionIter->second.Reactants)
                if (ReactionReactant.Counter > Particles[ReactionReactant.Identifier].Counter)
                {
                    IsPossible = false;
                    break;
                }

            if (IsPossible == true)
            {
                for (auto& ReactionReactant : ReactionIter->second.Reactants)
                    LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.Identifier << " "  << Particles[ReactionReactant.Identifier].Symbol << " " << Particles[ReactionReactant.Identifier].Counter));
                for (auto& ReactionProduct : ReactionIter->second.Products)
                    LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.Identifier << " "  << Particles[ReactionProduct.Identifier].Symbol << " " << Particles[ReactionProduct.Identifier].Counter));
                LoggersManagerObject.Log(STREAM(""));

                for (auto& ReactionReactant : ReactionIter->second.Reactants)
                    Particles[ReactionReactant.Identifier].Counter -=  ReactionReactant.Counter;

                for (auto& ReactionProduct : ReactionIter->second.Products)
                    Particles[ReactionProduct.Identifier].Counter +=  ReactionProduct.Counter;

                for (auto& ReactionReactant : ReactionIter->second.Reactants)
                    LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.Identifier << " "  << Particles[ReactionReactant.Identifier].Symbol << " " << Particles[ReactionReactant.Identifier].Counter));
                for (auto& ReactionProduct : ReactionIter->second.Products)
                    LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.Identifier << " " << Particles[ReactionProduct.Identifier].Symbol << " " << Particles[ReactionProduct.Identifier].Counter));
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

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedIntType NumberOfSteps)
{
    for (UnsignedIntType Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToDoRandomReaction(2);
}

void TestWellStirredChemicalReactionsSimulation()
{
    CellEngineWellStirredChemicalReactionsSimulation sim;

    sim.AddParticles({ 0, "Water", "H2O", false, 100 });
    sim.AddParticles({ 1, "Glucose", "C6H12O6", false, 50 });
    sim.AddParticles({ 2, "Oxygen", "0", false, 10 });
    sim.AddParticles({ 3, "Carbon dioxide", "CO2" , false, 5 });
    sim.AddParticles({ 4, "Eten", "CH2CH2", false, 15 });
    sim.AddParticles({ 5, "Ethanol", "CH3CH2(OH)", false, 25 });
    sim.AddParticles({ 6, "Propen", "CH3CHCH2", false, 5 });
    sim.AddParticles({ 7, "HX", "HX", false, 10 });
    sim.AddParticles({ 8, "2Halogenopropan", "CH3CHXCH3",false, 10 });
    sim.AddParticles({ 9, "Eten", "CH2CH2", false, 10 });
    sim.AddParticles({ 10, "Ethylene", "CH2CH2O", false, 10 });

    sim.AddReaction(ReactionType("C6H12O6 + O6 + ", { {1, 1}, {2, 6} }, { {3, 6}, {2, 6} }));
    sim.AddReaction(ReactionType("CH2CH2 + H2O + ", { {4, 1}, {0, 1} }, { {5, 1} }));
    sim.AddReaction(ReactionType("CH3CHCH2 + HX + ", { {6, 1}, {7, 1} }, { {8, 1} }));
    sim.AddReaction(ReactionType("CH2CH2 + O + ", { {9, 1}, {2, 1} }, { {10, 1} }));

    sim.Run(2);
}
