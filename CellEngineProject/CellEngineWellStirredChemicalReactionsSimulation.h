
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
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"

using namespace std;

class CellEngineWellStirredChemicalReactionsSimulation
{
public:
    std::mt19937_64 mt64X{ std::random_device{}() };
public:
    std::vector<ParticleKind> Particles;
    std::unordered_map<std::string, Reaction> Reactions;
public:
    void AddParticleKind(const ParticleKind& ParticleParam)
    {
        Particles.emplace_back(ParticleParam);
    }
public:
    void AddReaction(const Reaction& ReactionParam)
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
