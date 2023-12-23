
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_H

#include <string>

#include "ExceptionsMacro.h"
#include "CellEngineReaction.h"

class CellEngineChemicalReactions
{
private:
    std::mt19937_64 mt64R{ std::random_device{}() };
protected:
    std::vector<Reaction> Reactions;
    std::unordered_map<UnsignedInt, UnsignedInt> ReactionsPosFromId;
    //std::unordered_multimap<std::string, UnsignedInt> ReactionsPosFromString;
    std::multimap<std::string, UnsignedInt> ReactionsPosFromString;
public:
    void PreprocessChemicalReactions();
public:
    void AddChemicalReaction(const Reaction& ReactionParam)
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsPosFromId.insert(std::make_pair(ReactionParam.Id, Reactions.size() - 1));
        ReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
public:
    virtual std::vector<UnsignedInt> GetRandomParticles(UnsignedInt NumberOfReactants) = 0;
    virtual bool IsChemicalReactionPossible(const Reaction& ReactionObject) = 0;
    virtual bool MakeChemicalReaction(Reaction& ReactionObject) = 0;
public:
    bool TryToMakeRandomChemicalReaction(UnsignedInt NumberOfReactants);
};

void ReadChemicalReactionsFromFile();

#endif
