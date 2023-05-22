
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

void CellEngineWellStirredChemicalReactionsSimulation::TryToDoRandomReaction(const UnsignedInt NumberOfReactants)
{
    LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

    std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(1, Particles.size());

    vector<UnsignedInt> RandomParticlesTypes;
    for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
    {
        RandomParticlesTypes.emplace_back(UniformDistributionObjectUint64t(mt64X));
        LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        RandomParticlesTypes.back()--;
        if (Particles[RandomParticlesTypes.back()].Counter == 0)
        {
            LoggersManagerObject.Log(STREAM("Number of random reactant is zero!"));
            return;
        }
    }

                                                                                                                        RandomParticlesTypes.clear(); RandomParticlesTypes = {9, 0};// zapisac jako testy jednostkowe gdzie rand jest zapisane

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

        //SPRAWDZIC CZY DNA SIE POKRYWA KOD DNA Z ZADANYM W REAKCJI - W STRINGU REAKCJI MUSI BYC WIEC TO ZADZIALA

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

void CellEngineWellStirredChemicalReactionsSimulation::Run(UnsignedInt NumberOfSteps)
{
    for (UnsignedInt Steps = 1; Steps <= NumberOfSteps; Steps++)
        TryToDoRandomReaction(2);
}

void TestWellStirredChemicalReactionsSimulation()
{
    CellEngineWellStirredChemicalReactionsSimulation sim;

    sim.AddParticleKind({ 0, "Water", "H2O", 100 });
    sim.AddParticleKind({ 1, "Glucose", "C6H12O6", 50 });
    sim.AddParticleKind({ 2, "Oxygen", "0", 10 });
    sim.AddParticleKind({ 3, "Carbon dioxide", "CO2", 5 });
    sim.AddParticleKind({ 4, "Eten", "CH2CH2", 15 });
    sim.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 25 });
    sim.AddParticleKind({ 6, "Propen", "CH3CHCH2", 5 });
    sim.AddParticleKind({ 7, "HX", "HX", 10 });
    sim.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 10 });
    sim.AddParticleKind({ 9, "Eten", "CH2CH2", 10 });
    sim.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 10 });
    sim.AddParticleKind({ 11, "DNA", "CGATATTAAATAGGGCCT", 10 });

    sim.AddReaction(Reaction("C6H12O6 + O6 + ", { { 1, 1 }, { 2, 6 } }, { { 3, 6 }, { 2, 6 } }));
    sim.AddReaction(Reaction("CH2CH2 + H2O + ", { { 4, 1 }, { 0, 1 } }, { { 5, 1 } }));
    sim.AddReaction(Reaction("CH3CHCH2 + HX + ", { { 6, 1 }, { 7, 1 } }, { { 8, 1 } }));
    sim.AddReaction(Reaction("CH2CH2 + O + ", { { 9, 1 }, { 2, 1 } }, { { 10, 1 } }));

    sim.Run(2);
}
