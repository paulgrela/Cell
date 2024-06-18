
#ifndef CELL_ENGINE_REACTION_H
#define CELL_ENGINE_REACTION_H

#include <string>
#include <vector>

#include <functional>

#include "CellEngineTypes.h"
#include "CellEngineUseful.h"

class ParticleKind;

class CellEngineChemicalReactionsInSimulationSpace;

class ParticleKindForChemicalReaction
{
public:
    EntityIdInt EntityId{};
    UnsignedInt Counter{};
    bool ToRemoveInReaction;
    std::vector<UniqueIdInt> LinkedParticleTypes;
public:
    std::string SequenceStr;
    std::vector<ChainIdInt> Sequence;
public:
    ParticleKindForChemicalReaction(EntityIdInt EntityIdParam, UnsignedInt CounterParam, std::string SequenceStrParam, bool ToRemoveInReactionParam) : EntityId(EntityIdParam), Counter(CounterParam), SequenceStr(std::move(SequenceStrParam)), ToRemoveInReaction(ToRemoveInReactionParam)
    {
        for (auto& NucleotideLetter : SequenceStr)
            Sequence.emplace_back(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(NucleotideLetter));
    }
    ParticleKindForChemicalReaction(EntityIdInt EntityIdParam, UnsignedInt CounterParam, std::string SequenceStrParam, bool ToRemoveInReactionParam, std::vector<UniqueIdInt> LinkedParticlesTypesParam) : ParticleKindForChemicalReaction(EntityIdParam, CounterParam, std::move(SequenceStrParam), ToRemoveInReactionParam)
    {
        LinkedParticleTypes = std::move(LinkedParticlesTypesParam);
    }
};

class ChemicalReaction
{
    using SpecialReactionFunctionType = std::function<bool (CellEngineChemicalReactionsInSimulationSpace*, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>&, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>&, const ChemicalReaction&)>;
public:
    UnsignedInt ReactionIdNum{};
    std::string ReactionIdStr;
    std::string ReactionName;
public:
    std::string ReactantsStr;
public:
    std::uint32_t Duration{};
public:
    bool Reversible{};
    std::string UpperFluxBound;
    std::string LowerFluxBound;
public:
    UnsignedInt AdditionalParameter1 = 0;
    UnsignedInt AdditionalParameter2 = 0;
public:
    SpecialReactionFunctionType SpecialReactionFunction;
public:
    std::vector<ParticleKindForChemicalReaction> Reactants;
    std::vector<ParticleKindForChemicalReaction> Products;
public:
    std::string Comment;
public:
    ChemicalReaction() = default;
public:
    ChemicalReaction(UnsignedInt IdNumParam, std::string ReactionNameParam, std::string ReactantsStrParam, std::vector<ParticleKindForChemicalReaction> ReactantsParam, std::vector<ParticleKindForChemicalReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam, std::string CommentParam = "") : ReactionIdNum(IdNumParam), ReactionName(std::move(ReactionNameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), SpecialReactionFunction(std::move(SpecialReactionFunctionParam)), Comment(std::move(CommentParam))
    {
    }
    ChemicalReaction(UnsignedInt IdNumParam, std::string ReactionNameParam, std::string ReactantsStrParam, UnsignedInt AdditionalParameterParam1, UnsignedInt AdditionalParameterParam2, std::vector<ParticleKindForChemicalReaction> ReactantsParam, std::vector<ParticleKindForChemicalReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam, std::string CommentParam = "") : ReactionIdNum(IdNumParam), ReactionName(std::move(ReactionNameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), AdditionalParameter1(AdditionalParameterParam1), AdditionalParameter2(AdditionalParameterParam2), SpecialReactionFunction(std::move(SpecialReactionFunctionParam)), Comment(std::move(CommentParam))
    {
    }
};

#endif
