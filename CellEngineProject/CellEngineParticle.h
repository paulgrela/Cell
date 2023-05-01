
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
    EntityIdInt EntityId{};
    ChainIdInt ChainId{};
    vector3_16 UniqueColor{};
public:
    std::vector<vector3_64> ListOfVoxels;
public:
    explicit Particle(std::vector<vector3_64>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {
    }
    explicit Particle(EntityIdInt EntityIdParam, ChainIdInt ChainIdParam, vector3_16 UniqueColorParam) : EntityId(EntityIdParam), ChainId(ChainIdParam), UniqueColor(UniqueColorParam)
    {
    }
public:
    Particle() = default;
};

class ParticleKind
{
public:
    UnsignedInt EntityId{};
    std::string Name;
    std::string Symbol;
    UnsignedInt Counter;
public:
    std::list<Particle> ParticlesObjects;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    ParticleKind(UnsignedInt EntityIdParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam, std::unordered_map<std::string, UnsignedInt> ReactionsIdByStringParam) : EntityId(EntityIdParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam), ReactionsIdByString(std::move(ReactionsIdByStringParam))
    {
    }
    ParticleKind(UnsignedInt EntityIdParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam)
    {
    }
    ParticleKind(UnsignedInt EntityIdParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Counter(CounterParam)
    {
    }
};

#endif
