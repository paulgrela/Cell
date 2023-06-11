
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
    void PreprocessReactions();
public:
    void AddChemicalReaction(const Reaction& ReactionParam)
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) = 0;
    virtual bool IsChemicalReactionPossible(const Reaction& ReactionObject) = 0;
    virtual void MakeChemicalReaction(Reaction& ReactionObject) = 0;
public:
    bool TryToMakeRandomChemicalReaction(UnsignedInt NumberOfReactants);
};

void ReadChemicalReactionsFromFile();

#endif
