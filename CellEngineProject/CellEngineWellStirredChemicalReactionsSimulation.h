
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
#include "CellEngineChemicalReaction.h"
#include "CellEngineChemicalReactionsEngine.h"

using namespace std;

class CellEngineWellStirredChemicalReactionsSimulation : public CellEngineChemicalReactionsEngine
{
public:
    std::mt19937_64 mt64X{ std::random_device{}() };
public:
    std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants, const CurrentThreadPosType& CurrentThreadPos) override;
    bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos)override;
    bool MakeChemicalReaction(ChemicalReaction& ReactionObject, const CurrentThreadPosType& CurrentThreadPos) override;
public:
    void Run(UnsignedInt NumberOfSteps);
};

#endif
