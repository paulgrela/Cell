
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H

#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesKindsManager.h"

class CellEngineChemicalReactionsManager : virtual public CellEngineRandomDeviceEngine
{
public:
    std::vector<ChemicalReaction> Reactions;
    std::unordered_map<UnsignedInt, UnsignedInt> ReactionsPosFromId;
    std::unordered_multimap<std::string, UnsignedInt> ReactionsPosFromString;
public:
    ChemicalReaction& GetReactionFromNumId(UnsignedInt ReactionId)
    {
        return Reactions[ReactionsPosFromId.find(ReactionId)->second];
    }
public:
    void PreprocessChemicalReactions()
    {
        for (auto& ReactionObject : Reactions)
            sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForChemicalReaction& PK1, ParticleKindForChemicalReaction& PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); } );
    }
public:
    void AddChemicalReaction(const ChemicalReaction& ReactionParam)
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsPosFromId.insert(std::make_pair(ReactionParam.ReactionIdNum, Reactions.size() - 1));
        ReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
};

inline CellEngineChemicalReactionsManager CellEngineChemicalReactionsManagerObject;

#endif
