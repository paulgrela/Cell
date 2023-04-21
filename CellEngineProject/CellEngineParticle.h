
#ifndef CELL_ENGINE_PARTICLE_H
#define CELL_ENGINE_PARTICLE_H

#include <list>
#include <vector>
#include <string>

#include "CellEngineTypes.h"

class Particle
{
public:
    bool SelectedForReaction;
public:
    std::vector<vector3_64> ListOfVoxels;
public:
    explicit Particle(std::vector<vector3_64> ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {}
};

class ParticleKind
{
public:
    UnsignedInt Identifier{};
    std::string Name;
    std::string Symbol;
    UnsignedInt Counter;
    std::list<Particle> ParticlesObjects;
public:
    ParticleKind(UnsignedInt IdentifierParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam)
    {}
    ParticleKind(UnsignedInt IdentifierParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), Counter(CounterParam)
    {}
};

#endif
