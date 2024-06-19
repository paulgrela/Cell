
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H

#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesKindsManager.h"

class ChemicalReactionsManager : virtual public CellEngineRandomDeviceEngine
{
public:
    std::vector<ChemicalReaction> ChemicalReactions;
    std::unordered_map<UnsignedInt, UnsignedInt> ChemicalReactionsPosFromId;
    std::unordered_multimap<std::string, UnsignedInt> ChemicalReactionsPosFromString;
public:
    ChemicalReaction& GetReactionFromNumId(UnsignedInt ReactionId)
    {
        return ChemicalReactions[ChemicalReactionsPosFromId.find(ReactionId)->second];
    }
public:
    void PreprocessChemicalReactions()
    {
        for (auto& ReactionObject : ChemicalReactions)
            sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForChemicalReaction& PK1, ParticleKindForChemicalReaction& PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); } );
    }
public:
    void AddChemicalReaction(const ChemicalReaction& ReactionParam)
    {
        ChemicalReactions.emplace_back(ReactionParam);
        ChemicalReactionsPosFromId.insert(std::make_pair(ReactionParam.ReactionIdNum, ChemicalReactions.size() - 1));
        ChemicalReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, ChemicalReactions.size() - 1));
    }
};

inline ChemicalReactionsManager ChemicalReactionsManagerObject;

#endif
