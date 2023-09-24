
#ifndef CELL_ENGINE_REACTION_H
#define CELL_ENGINE_REACTION_H

#include <string>
#include <vector>

#include "CellEngineTypes.h"

class ParticleKind;

class Reaction
{
public:
    UnsignedInt Id{};
    std::string Name;
public:
    std::uint32_t Duration{};
public:
    UnsignedInt AdditionalParameter = 0;
public:
    std::string ReactantsStr;
public:
    std::vector<ParticleKindForReaction> Reactants;
    std::vector<ParticleKindForReaction> Products;
public:
    Reaction() = delete;
    Reaction(UnsignedInt IdParam, std::string NameParam, std::string ReactantsStrParam, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam) : Id(IdParam), Name(std::move(NameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam))
    {
    }
    Reaction(UnsignedInt IdParam, std::string NameParam, std::string ReactantsStrParam, UnsignedInt AdditionalParameterParam, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam) : Id(IdParam), Name(std::move(NameParam)), ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam)), AdditionalParameter(AdditionalParameterParam)
    {
    }
};

#endif
