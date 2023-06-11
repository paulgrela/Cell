
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
#include "CellEngineChemicalReactions.h"

using namespace std;

class CellEngineWellStirredChemicalReactionsSimulation : public CellEngineChemicalReactions
{
public:
    std::mt19937_64 mt64X{ std::random_device{}() };
public:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) override;
    bool IsChemicalReactionPossible(const Reaction& ReactionObject)override;
    void MakeChemicalReaction(Reaction& ReactionObject) override;
public:
    void Run(UnsignedInt NumberOfSteps);
};

#endif
