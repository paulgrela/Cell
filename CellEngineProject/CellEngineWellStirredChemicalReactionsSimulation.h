
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
    void TryDoReaction()
    {
        LoggersManagerObject.Log(STREAM(endl << "REACTION" << endl));

        std::mt19937_64 mt64X{ std::random_device{}() };
        std::uniform_int_distribution<uint64_t> UniformDistributionObject1X_Uint64t(1, Particles.size());
        UnsignedIntType RandomParticleType1 = UniformDistributionObject1X_Uint64t(mt64X);
        UnsignedIntType RandomParticleType2 = UniformDistributionObject1X_Uint64t(mt64X);
        LoggersManagerObject.Log(STREAM("Particle 1 = (" << to_string(RandomParticleType1) << " | " << Particles[RandomParticleType1].Symbol << ") Particle 2 = (" << to_string(RandomParticleType2) << " | " << Particles[RandomParticleType2].Symbol << ")"));
        RandomParticleType1--;
        RandomParticleType2--;

        if (RandomParticleType1 != RandomParticleType2)
        {
            vector<string> ParticlesSymbolsForReactionToSort{ Particles[RandomParticleType1].Symbol, Particles[RandomParticleType2].Symbol };
            std::sort(ParticlesSymbolsForReactionToSort.begin(), ParticlesSymbolsForReactionToSort.end());
            string ReactionSymbolsStr = ParticlesSymbolsForReactionToSort[0] + " + " + ParticlesSymbolsForReactionToSort[1];

            LoggersManagerObject.Log(STREAM("Reaction Symbols = [" << ReactionSymbolsStr << "]" << endl));

            auto ReactionIter = Reactions.find(ReactionSymbolsStr);
            if (ReactionIter != Reactions.end())
            {
                LoggersManagerObject.Log(STREAM("Reaction parameters [ Reactant 0 Counter = " << ReactionIter->second.Reactants[0].Counter << " Particle Counter = " << Particles[RandomParticleType1].Counter << "]"));
                LoggersManagerObject.Log(STREAM("Reaction parameters [ Reactant 1 Counter = " << ReactionIter->second.Reactants[1].Counter << " Particle Counter = " << Particles[RandomParticleType2].Counter << "]"));

                bool IsPossible = true;
                if (ReactionIter->second.Reactants[0].Counter > Particles[RandomParticleType1].Counter)
                    IsPossible = false;
                if (ReactionIter->second.Reactants[1].Counter > Particles[RandomParticleType2].Counter)
                    IsPossible = false;
                if (IsPossible == true)
                {
                    LoggersManagerObject.Log(STREAM("BEFORE = (" << Particles[RandomParticleType1].Counter << " " << Particles[RandomParticleType2].Counter << ") "));

                    Particles[RandomParticleType1].Counter -=  ReactionIter->second.Reactants[0].Counter;
                    Particles[RandomParticleType2].Counter -=  ReactionIter->second.Reactants[1].Counter;

                    for (auto& ReactionProduct : ReactionIter->second.Products)
                        Particles[ReactionProduct.Identifier].Counter +=  ReactionProduct.Counter;

                    LoggersManagerObject.Log(STREAM("AFTER R1 = " << ReactionIter->second.Products[0].Identifier << " " << Particles[ReactionIter->second.Products[0].Identifier].Counter));
                    LoggersManagerObject.Log(STREAM("AFTER R2 = " << ReactionIter->second.Products[1].Identifier << " " << Particles[ReactionIter->second.Products[1].Identifier].Counter));

                    getchar();
                }
                else
                    LoggersManagerObject.Log(STREAM("Particles types are the same!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Reaction for particles does not exist!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Particles types for reaction are the same!"));

        LoggersManagerObject.Log(STREAM("END OF REACTION" << endl));
    }
public:
    void Run(UnsignedIntType NumberOfSteps)
    {
        for (UnsignedIntType Steps = 1; Steps <= NumberOfSteps; Steps++)
            TryDoReaction();
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

    sim.AddReaction(ReactionType("C6H12O6 + O6", { {2, 1}, {3, 6} }, { {4, 6}, {1, 6} }));
    sim.AddReaction(ReactionType("CH2CH2 + H2O", { {5, 1}, {1, 1} }, { {6, 1} }));
    sim.AddReaction(ReactionType("CH3CHCH2 + HX", { {7, 1}, {8, 1} }, { {9, 1} }));
    sim.AddReaction(ReactionType("CH2CH2 + O", { {10, 1}, {3, 1} }, { {11, 1} }));

    sim.Run(15);
}

#endif
