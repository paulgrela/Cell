
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_H

#include <string>

#include "ExceptionsMacro.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineRandomDeviceEngine.h"

class CellEngineChemicalReactionsEngine : virtual public CellEngineRandomDeviceEngine
{
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants) = 0;
    virtual bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) = 0;
    virtual bool MakeChemicalReaction(ChemicalReaction& ReactionObject) = 0;
public:
    bool TryToMakeRandomChemicalReaction(UnsignedInt NumberOfReactants, UnsignedInt MaxNumberOfReactants);
};

#endif
