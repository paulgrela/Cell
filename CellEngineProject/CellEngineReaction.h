
#ifndef CELL_ENGINE_REACTION_H
#define CELL_ENGINE_REACTION_H

#include <string>
#include <vector>

#include "CellEngineTypes.h"

class ParticleKind;

class Reaction
{
private:
    UnsignedInt Identifier{};
public:
    std::uint32_t Duration{};
public:
    std::string ReactantsStr;
public:
    std::vector<ParticleKindForReaction> Reactants;
    std::vector<ParticleKindForReaction> Products;
public:
    Reaction() = delete;
    Reaction(std::string ReactantsStrParam, std::vector<ParticleKindForReaction> ReactantsParam, std::vector<ParticleKindForReaction> ProductsParam) : ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam))
    {
    }
};

#endif
