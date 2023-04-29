
#ifndef CELL_ENGINE_PARTICLE_H
#define CELL_ENGINE_PARTICLE_H

#include <list>
#include <vector>
#include <string>
#include <unordered_map>

#include "CellEngineTypes.h"

class Particle
{
public:
    bool SelectedForReaction{};
public:
    UnsignedInt UniqueIdentifier{};
public:
    std::vector<vector3_64> ListOfVoxels;
public:
    explicit Particle(std::vector<vector3_64>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {}
public:
    Particle() = default;
};

class ParticleKind
{
public:
    UnsignedInt Identifier{};
    ChainIdInt ChainId{};
    std::string Name;
    std::string Symbol;
    UnsignedInt Counter;
public:
    std::list<Particle> ParticlesObjects;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    ParticleKind(UnsignedInt IdentifierParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam, std::unordered_map<std::string, UnsignedInt> ReactionsIdByStringParam) : Identifier(IdentifierParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam), ReactionsIdByString(std::move(ReactionsIdByStringParam))
    {}
    ParticleKind(UnsignedInt IdentifierParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam)
    {}
    ParticleKind(UnsignedInt IdentifierParam, UnsignedInt CounterParam) : Identifier(IdentifierParam), Counter(CounterParam)
    {}
};

#endif
