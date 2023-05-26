
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
    std::string ReactantsStr;
public:
    std::vector<ParticleKind> Reactants;
    std::vector<ParticleKind> Products;
public:
    Reaction() = delete;
    Reaction(std::string ReactantsStrParam, std::vector<ParticleKind> ReactantsParam, std::vector<ParticleKind> ProductsParam) : ReactantsStr(std::move(ReactantsStrParam)), Reactants(std::move(ReactantsParam)), Products(std::move(ProductsParam))
    {
    }
};

#endif
