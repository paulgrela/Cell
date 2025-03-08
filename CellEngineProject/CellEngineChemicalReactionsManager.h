
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
    std::unordered_map<std::string, UnsignedInt> ChemicalReactionsPosFromIdStr;
    std::unordered_multimap<std::string, UnsignedInt> ChemicalReactionsPosFromString;
public:
    void ClearAllChemicalReactions()
    {
        ChemicalReactions.clear();
        ChemicalReactionsPosFromId.clear();
        ChemicalReactionsPosFromIdStr.clear();
        ChemicalReactionsPosFromString.clear();
    }
public:
    void AddChemicalReaction(const ChemicalReaction& ReactionParam)
    {
        ChemicalReactions.emplace_back(ReactionParam);
        ChemicalReactionsPosFromId.insert(std::make_pair(ReactionParam.ReactionIdNum, ChemicalReactions.size() - 1));
        ChemicalReactionsPosFromIdStr.insert(std::make_pair(ReactionParam.ReactionIdStr, ChemicalReactions.size() - 1));
        ChemicalReactionsPosFromString.insert(std::make_pair(ReactionParam.ReactantsStr, ChemicalReactions.size() - 1));
    }
public:
    ChemicalReaction& GetReactionFromNumId(const UnsignedInt ReactionId)
    {
        return ChemicalReactions[ChemicalReactionsPosFromId.find(ReactionId)->second];
    }
private:
    void ReverseReactantsAndProductsBecauseOfFormerError()
    {
        for (auto& ChemicalReactionObject : ChemicalReactions)
        {
            for (auto& ReactantObject : ChemicalReactionObject.Reactants)
                if (ParticlesKindsManagerObject.GetParticleKind(ReactantObject.EntityId).IdStr.substr(0, 10) == JCVISYN3APredStr)
                {
                    ReactantObject.ToRemoveInReaction = false;
                    ChemicalReactionObject.Products.emplace_back(ReactantObject);
                }
            erase_if(ChemicalReactionObject.Reactants, [](auto& R){ return ParticlesKindsManagerObject.GetParticleKind(R.EntityId).IdStr.substr(0, 10) == JCVISYN3APredStr; });

            swap(ChemicalReactionObject.Products, ChemicalReactionObject.Reactants);
        }
    }
    void SortReactantsAndProductsForEveryChemicalReaction()
    {
        for (auto& ReactionObject : ChemicalReactions)
        {
            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::VoxelSimulationSpace)
                sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForChemicalReaction &PK1, ParticleKindForChemicalReaction &PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfVoxels.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfVoxels.size(); });
            else
            if (CellEngineConfigDataObject.TypeOfSpace == CellEngineConfigData::TypesOfSpace::FullAtomSimulationSpace)
                sort(ReactionObject.Products.begin(), ReactionObject.Products.end(), [](ParticleKindForChemicalReaction &PK1, ParticleKindForChemicalReaction &PK2){ return ParticlesKindsManagerObject.GetParticleKind(PK1.EntityId).ListOfAtoms.size() > ParticlesKindsManagerObject.GetParticleKind(PK2.EntityId).ListOfAtoms.size(); });

            ReactionObject.ReactantsStr = GetStringOfSortedParticlesDataNames(ReactionObject.Reactants);
        }
    }
public:
    void PreprocessAllChemicalReactions()
    {
        if (CellEngineConfigDataObject.ReverseReactantsAndProductsBecauseOfFormerErrorBool == true)
            ReverseReactantsAndProductsBecauseOfFormerError();

        SortReactantsAndProductsForEveryChemicalReaction();

        MaxNumberOfReactants = std::max_element(ChemicalReactions.begin(), ChemicalReactions.end(), [](const ChemicalReaction& R1, const ChemicalReaction& R2){ return R1.Reactants.size() < R2.Reactants.size(); })->Reactants.size();
    }
public:
    static std::string GetStringOfSortedParticlesDataNames(const std::vector<ParticleKindForChemicalReaction>& ReactionArguments)
    {
        std::string KeyStringOfReaction;

        try
        {
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
