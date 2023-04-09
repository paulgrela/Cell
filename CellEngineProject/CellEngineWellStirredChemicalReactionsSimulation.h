
#ifndef CELL_ENGINE_WELL_STIRRED_CHEMICAL_REACTIONS_SIMULATION_H
#define CELL_ENGINE_WELL_STIRRED_CHEMICAL_REACTIONS_SIMULATION_H

#include <list>
#include <string>
#include <vector>
#include <random>
#include <iostream>
#include <unordered_map>

#include "CellEngineTypes.h"

using namespace std;

class Particle
{
public:
    UnsignedIntType Identifier{};
    std::string Name;
    std::string Symbol;
    bool SelectedForReaction;
    UnsignedIntType Counter;
public:
    Particle(UnsignedIntType IdentifierParam, std::string NameParam, std::string SymbolParam, bool SelectedForReactionParam, UnsignedIntType CounterParam) : Identifier(IdentifierParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), SelectedForReaction(SelectedForReactionParam), Counter(CounterParam)
    {}
    Particle(UnsignedIntType IdentifierParam, UnsignedIntType CounterParam) : Identifier(IdentifierParam), SelectedForReaction(false), Counter(CounterParam)
    {}
};

class ReactionType
{
private:
    UnsignedIntType Identifier{};
public:
    std::string ReactantsStr;
public:
    std::vector<Particle> Reactants;
    std::vector<Particle> Products;
public:
    ReactionType() = delete;
    ReactionType(std::string ReactantsStrParam, std::vector<Particle> ReactantsParam, std::vector<Particle> ProductsParam) : ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam))
    {
    }
};

class CellEngineWellStirredChemicalReactionsSimulation
{
public:
    std::vector<Particle> Particles;
    std::unordered_map<std::string, ReactionType> Reactions;
public:
    void AddParticles(const Particle& ParticleParam)
    {
        Particles.emplace_back(ParticleParam);
    }
public:
    void AddReaction(const ReactionType& ReactionParam)
    {
        Reactions.insert(std::make_pair(ReactionParam.ReactantsStr, ReactionParam));
    }
public:
    void TryDoReaction(const UnsignedIntType NumberOfReactants)
    {
        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

        std::mt19937_64 mt64X{ std::random_device{}() };
        std::uniform_int_distribution<uint64_t> UniformDistributionObject1X_Uint64t(1, Particles.size());

        vector<UnsignedIntType> RandomParticlesTypes;
        for (UnsignedIntType ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(UniformDistributionObject1X_Uint64t(mt64X));
            LoggersManagerObject.Log(STREAM("Particle Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
            RandomParticlesTypes.back()--;
        }
        
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
                        LoggersManagerObject.Log(STREAM("BEFORE REACTANT = " << ReactionReactant.Identifier << " " << Particles[ReactionReactant.Identifier].Counter));
                    for (auto& ReactionProduct : ReactionIter->second.Products)
                        LoggersManagerObject.Log(STREAM("BEFORE PRODUCT = " << ReactionProduct.Identifier << " " << Particles[ReactionProduct.Identifier].Counter));
                    LoggersManagerObject.Log(STREAM(""));

                    for (auto& ReactionReactant : ReactionIter->second.Reactants)
                        Particles[ReactionReactant.Identifier].Counter -=  ReactionReactant.Counter;

                    for (auto& ReactionProduct : ReactionIter->second.Products)
                        Particles[ReactionProduct.Identifier].Counter +=  ReactionProduct.Counter;

                    for (auto& ReactionReactant : ReactionIter->second.Reactants)
                        LoggersManagerObject.Log(STREAM("AFTER REACTANT = " << ReactionReactant.Identifier << " " << Particles[ReactionReactant.Identifier].Counter));
                    for (auto& ReactionProduct : ReactionIter->second.Products)
                        LoggersManagerObject.Log(STREAM("AFTER PRODUCT = " << ReactionProduct.Identifier << " " << Particles[ReactionProduct.Identifier].Counter));
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
public:
    void Run(UnsignedIntType NumberOfSteps)
    {
        for (UnsignedIntType Steps = 1; Steps <= NumberOfSteps; Steps++)
            TryDoReaction(2);
    }
};

void TestWellStirredChemicalReactionsSimulation()
{
    CellEngineWellStirredChemicalReactionsSimulation sim;

    sim.AddParticles({ 1, "Water", "H2O", false, 100 });
    sim.AddParticles({ 2, "Glucose", "C6H12O6", false, 50 });
    sim.AddParticles({ 3, "Oxygen", "0", false, 10 });
    sim.AddParticles({ 4, "Carbon dioxide", "CO2" , false, 5 });
    sim.AddParticles({ 5, "Eten", "CH2CH2", false, 15 });
    sim.AddParticles({ 6, "Ethanol", "CH3CH2(OH)", false, 25 });
    sim.AddParticles({ 7, "Propen", "CH3CHCH2", false, 5 });
    sim.AddParticles({ 8, "HX", "HX", false, 10 });
    sim.AddParticles({ 9, "2Halogenopropan", "CH3CHXCH3",false, 10 });
    sim.AddParticles({ 10, "Eten", "CH2CH2", false, 10 });
    sim.AddParticles({ 11, "Ethylene", "CH2CH2O", false, 10 });

    sim.AddReaction(ReactionType("C6H12O6 + O6 + ", { {2, 1}, {3, 6} }, { {4, 6}, {1, 6} }));
    sim.AddReaction(ReactionType("CH2CH2 + H2O + ", { {5 - 1, 1}, {1 - 1, 1} }, { {6 - 1, 1} }));
    sim.AddReaction(ReactionType("CH3CHCH2 + HX + ", { {7, 1}, {8, 1} }, { {9, 1} }));
    sim.AddReaction(ReactionType("CH2CH2 + O + ", { {10, 1}, {3, 1} }, { {11, 1} }));

    sim.Run(1);
}

#endif
