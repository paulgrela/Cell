
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
    UniqueIdInt Index{};
    UniqueIdInt GenomeIndex{};
    vector3_16 UniqueColor{};
public:
    std::vector<vector3_16> ListOfVoxels;
public:
    explicit Particle(std::vector<vector3_16>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {
    }
    explicit Particle(UniqueIdInt IndexParam, EntityIdInt EntityIdParam, ChainIdInt ChainIdParam, UniqueIdInt GenomeIndexParam, vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeIndex(GenomeIndexParam), UniqueColor(UniqueColorParam)
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
    std::string SequenceStr;
    std::vector<ChainIdInt> Sequence;
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
