
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_CREATOR_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_CREATOR_H

#include "CellEngineChemicalReactions.h"

class CellEngineChemicalReactionsCreator : virtual public CellEngineChemicalReactions
{
public:
    static void AddParticlesKinds();
    static void ReadChemicalReactionsFromFile();
    void AddChemicalReactions();
};

#endif
