
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_H

#include "CellEngineChemicalReaction.h"
#include "CellEngineRandomDeviceEngine.h"

class CellEngineChemicalReactionsEngine : virtual public CellEngineRandomDeviceEngine
{
public:
    virtual bool IsChemicalReactionPossible(const ChemicalReaction& ReactionObject) = 0;
    virtual bool MakeChemicalReaction(ChemicalReaction& ReactionObject) = 0;
};

#endif
