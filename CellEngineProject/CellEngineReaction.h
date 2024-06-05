
#ifndef CELL_ENGINE_REACTION_H
#define CELL_ENGINE_REACTION_H

#include <string>
#include <vector>

#include <functional>

#include "CellEngineTypes.h"

class ParticleKind;

class CellEngineChemicalReactionsInSimulationSpace;

class Reaction
{
    using SpecialReactionFunctionType = std::function<bool (CellEngineChemicalReactionsInSimulationSpace*, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>&, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>&, const Reaction&)>;
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
    std::vector<ParticleKindForReaction> Reactants;
    std::vector<ParticleKindForReaction> Products;
public:
    std::string Comment;
public:
    Reaction() = default;
public:
    Reaction(UnsignedInt IdNumParam, std::string ReactionNameParam, std::string ReactantsStrParam, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam, std::string CommentParam = "") : ReactionIdNum(IdNumParam), ReactionName(std::move(ReactionNameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), SpecialReactionFunction(std::move(SpecialReactionFunctionParam)), Comment(std::move(CommentParam))
    {
    }
    Reaction(UnsignedInt IdNumParam, std::string ReactionNameParam, std::string ReactantsStrParam, UnsignedInt AdditionalParameterParam1, UnsignedInt AdditionalParameterParam2, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam, std::string CommentParam = "") : ReactionIdNum(IdNumParam), ReactionName(std::move(ReactionNameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), AdditionalParameter1(AdditionalParameterParam1), AdditionalParameter2(AdditionalParameterParam2), SpecialReactionFunction(std::move(SpecialReactionFunctionParam)), Comment(std::move(CommentParam))
    {
    }
};

#endif
