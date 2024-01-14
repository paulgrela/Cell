
#ifndef CELL_ENGINE_REACTION_H
#define CELL_ENGINE_REACTION_H

#include <string>
#include <vector>

#include <functional>

#include "CellEngineTypes.h"

class ParticleKind;

class CellEngineVoxelSimulationSpace;

class Reaction
{
    using SpecialReactionFunctionType = std::function<bool (CellEngineVoxelSimulationSpace*, const std::vector<std::pair<UniqueIdInt, UnsignedInt>>&, const Reaction&)>;
public:
    UnsignedInt Id{};
    std::string Name;
public:
    std::uint32_t Duration{};
public:
    UnsignedInt AdditionalParameter1 = 0;
    UnsignedInt AdditionalParameter2 = 0;
public:
    std::string ReactantsStr;
public:
    SpecialReactionFunctionType SpecialReactionFunction;
public:
    std::vector<ParticleKindForReaction> Reactants;
    std::vector<ParticleKindForReaction> Products;
public:
    Reaction() = delete;
    Reaction(UnsignedInt IdParam, std::string NameParam, std::string ReactantsStrParam, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam) : Id(IdParam), Name(std::move(NameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), SpecialReactionFunction(std::move(SpecialReactionFunctionParam))
    {
    }
    Reaction(UnsignedInt IdParam, std::string NameParam, std::string ReactantsStrParam, UnsignedInt AdditionalParameterParam1, UnsignedInt AdditionalParameterParam2, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam, SpecialReactionFunctionType SpecialReactionFunctionParam) : Id(IdParam), Name(std::move(NameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), AdditionalParameter1(AdditionalParameterParam1), AdditionalParameter2(AdditionalParameterParam2), SpecialReactionFunction(std::move(SpecialReactionFunctionParam))
    {
    }
};

#endif
