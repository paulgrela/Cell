
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H

#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesKindsManager.h"

class CellEngineChemicalReactionsManager : virtual public CellEngineRandomDeviceEngine
{
public:
    std::vector<Reaction> Reactions;
    std::unordered_map<UnsignedInt, UnsignedInt> ReactionsPosFromId;
    std::unordered_multimap<std::string, UnsignedInt> ReactionsPosFromString;
public:
    Reaction& GetReactionFromNumId(UnsignedInt ReactionId)
    {
        return Reactions[ReactionsPosFromId.find(ReactionId)->second];
    }
public:
    void PreprocessChemicalReactions()
    {
        for (auto& ReactionObject : Reactions)
            sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForReaction& PK1, ParticleKindForReaction& PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); } );
    }
public:
    void AddChemicalReaction(const Reaction& ReactionParam)
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsPosFromId.insert(std::make_pair(ReactionParam.ReactionIdNum, Reactions.size() - 1));
        ReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
};

inline CellEngineChemicalReactionsManager CellEngineChemicalReactionsManagerObject;

#endif
