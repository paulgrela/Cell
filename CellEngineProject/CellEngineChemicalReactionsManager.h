
#ifndef CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H
#define CELL_ENGINE_CHEMICAL_REACTIONS_MANAGER_H

#include "CellEngineParticle.h"
#include "CellEngineChemicalReaction.h"
#include "CellEngineRandomDeviceEngine.h"
#include "CellEngineParticlesKindsManager.h"

class ChemicalReactionsManager : virtual public CellEngineRandomDeviceEngine
{
public:
    UnsignedInt MaxNumberOfReactants{};
public:
    std::vector<ChemicalReaction> ChemicalReactions;
    std::unordered_map<UnsignedInt, UnsignedInt> ChemicalReactionsPosFromId;
    std::unordered_multimap<std::string, UnsignedInt> ChemicalReactionsPosFromString;
public:
    void AddChemicalReaction(const ChemicalReaction& ReactionParam)
    {
        ChemicalReactions.emplace_back(ReactionParam);
        ChemicalReactionsPosFromId.insert(std::make_pair(ReactionParam.ReactionIdNum, ChemicalReactions.size() - 1));
        ChemicalReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, ChemicalReactions.size() - 1));
    }
public:
    ChemicalReaction& GetReactionFromNumId(UnsignedInt ReactionId)
    {
        return ChemicalReactions[ChemicalReactionsPosFromId.find(ReactionId)->second];
    }
public:
    void PreprocessChemicalReactions()
    {
        for (auto& ReactionObject : ChemicalReactions)
        {
            sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForChemicalReaction &PK1, ParticleKindForChemicalReaction &PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); });
            ReactionObject.ReactantsStr = GetStringOfSortedParticlesDataNames(ReactionObject.Reactants);
        }
        MaxNumberOfReactants = std::max_element(ChemicalReactions.begin(), ChemicalReactions.end(), [](const ChemicalReaction& R1, const ChemicalReaction& R2){ return R1.Reactants.size() < R2.Reactants.size(); })->Reactants.size();
    }
public:
    static std::string GetStringOfSortedParticlesDataNames(const std::vector<ParticleKindForChemicalReaction>& ReactionArguments)
    // static std::string GetStringOfSortedParticlesDataNames(std::vector<ParticleKindForChemicalReaction>& LocalData)
    {
        std::string KeyStringOfReaction;

        try
        {
            // sort(LocalData.begin(), LocalData.end(), [](const ParticleKindForChemicalReaction& P1, const ParticleKindForChemicalReaction& P2){ return ParticlesKindsManagerObject.GetParticleKind(P1.EntityId).Formula < ParticlesKindsManagerObject.GetParticleKind(P2.EntityId).Formula; } );
            //  for (const auto& LocalReactant : LocalData)
            //      KeyStringOfReaction += ParticlesKindsManagerObject.GetParticleKind(LocalReactant.EntityId).Formula + "+";

            auto LocalReactionArguments(ReactionArguments);
            sort(LocalReactionArguments.begin(), LocalReactionArguments.end(), [](const ParticleKindForChemicalReaction& P1, const ParticleKindForChemicalReaction& P2){ return ParticlesKindsManagerObject.GetParticleKind(P1.EntityId).IdStr < ParticlesKindsManagerObject.GetParticleKind(P2.EntityId).IdStr; } );
            for (const auto& LocalReactant : LocalReactionArguments)
                KeyStringOfReaction += ParticlesKindsManagerObject.GetParticleKind(LocalReactant.EntityId).IdStr + "+";
        }
        CATCH("getting string of sorted particles data names")

        return KeyStringOfReaction;
    }
public:
    void PrintChemicalReactions()
    {
        LoggersManagerObject.Log(STREAM("Number of chemical reactions = " << ChemicalReactions.size()));

        LoggersManagerObject.Log(STREAM("Maximal number of reactants = " << MaxNumberOfReactants));

        for (const auto& ChemicalReactionObject : ChemicalReactions)
            LoggersManagerObject.Log(STREAM("ReactionId = " << ChemicalReactionObject.ReactionIdNum << " ReactantsStr = |" << ChemicalReactionObject.ReactantsStr << "| ReactionIdStr = " << ChemicalReactionObject.ReactionIdStr << " Reactants Size " << ChemicalReactionObject.Reactants.size() << " Products Size " << ChemicalReactionObject.Products.size()));
    }
};

inline ChemicalReactionsManager ChemicalReactionsManagerObject;

#endif
