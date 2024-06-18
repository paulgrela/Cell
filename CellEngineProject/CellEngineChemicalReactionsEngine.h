
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_H

#include <string>

#include "ExceptionsMacro.h"
#include "CellEngineReaction.h"
#include "CellEngineRandomDeviceEngine.h"

class CellEngineChemicalReactionsEngine : virtual public CellEngineRandomDeviceEngine
{
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) = 0;
    virtual bool IsChemicalReactionPossible(const Reaction& ReactionObject) = 0;
    virtual bool MakeChemicalReaction(Reaction& ReactionObject) = 0;
public:
    bool TryToMakeRandomChemicalReaction(UnsignedInt NumberOfReactants);
};

#endif
