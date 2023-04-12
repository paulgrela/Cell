
#ifndef CELL_ENGINE_WELL_STIRRED_CHEMICAL_REACTIONS_SIMULATION_H
#define CELL_ENGINE_WELL_STIRRED_CHEMICAL_REACTIONS_SIMULATION_H

#include <list>
#include <string>
#include <vector>
#include <random>
#include <iostream>
#include <unordered_map>

#include "Logger.h"
#include "CellEngineTypes.h"

using namespace std;

class Particle
{
public:
    UnsignedInt Identifier{};
    std::string Name;
    std::string Symbol;
    bool SelectedForReaction;
    UnsignedInt Counter;
public:
    Particle(UnsignedInt IdentifierParam, std::string NameParam, std::string SymbolParam, bool SelectedForReactionParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), SelectedForReaction(SelectedForReactionParam), Counter(CounterParam)
    {}
    Particle(UnsignedInt IdentifierParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), SelectedForReaction(false), Counter(CounterParam)
    {}
};

class ReactionType
{
private:
    UnsignedInt Identifier{};
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
    std::mt19937_64 mt64X{ std::random_device{}() };
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
    void TryToDoRandomReaction(UnsignedInt NumberOfReactants);
public:
    void Run(UnsignedInt NumberOfSteps);
};

void TestWellStirredChemicalReactionsSimulation();

#endif
