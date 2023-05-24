
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
    std::vector<Reaction> Reactions;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    void AddParticleKind(const ParticleKind& ParticleParam);
    void AddReaction(const Reaction& ReactionParam);
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants);
    void TryToDoRandomReaction(UnsignedInt NumberOfReactants);
public:
    void Run(UnsignedInt NumberOfSteps);
};

#endif
