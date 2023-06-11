
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_H

#include <string>

#include "ExceptionsMacro.h"
#include "CellEngineReaction.h"

class CellEngineChemicalReactions
{
protected:
    std::vector<Reaction> Reactions;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    CellEngineChemicalReactions();
public:
    void AddReaction(const Reaction& ReactionParam)
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) = 0;
    virtual bool IsReactionPossible(const Reaction& ReactionObject) = 0;
    virtual void MakeReaction(Reaction& ReactionObject) = 0;
public:
    bool TryToMakeRandomReaction(UnsignedInt NumberOfReactants);
};

void ReadChemicalReactionsFromFile();

#endif
