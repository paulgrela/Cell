
#ifndef CELL_ENGINE_PARTICLE_H
#define CELL_ENGINE_PARTICLE_H

#include <list>
#include <vector>
#include <string>
#include <unordered_map>

#include "CellEngineTypes.h"
#include "CellEngineReaction.h"

#define PARTICLES_IN_VECTOR_

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
    std::vector<vector3_16> ListOfVoxels;
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

struct GraphicParticleKind
{
    UnsignedInt Identifier;
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3_16 AtomColor;
    vector3_16 ParticleColor;
    vector3_16 RandomParticleColor;
    std::string NameFromXML;
    std::string NameFromDataFile;
};

struct GraphicAtomKind
{
    std::string Name;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3_16 Color;
    vmath::vec3 ColorVmathVec3;
};

inline bool operator==(const GraphicAtomKind& AtomKindParameter, const std::string& NameStr)
{
    return AtomKindParameter.Name == NameStr;
}

class CellEngineSimulationManager
{
public:
    std::string GenomeLine;
public:
    std::vector<UniqueIdInt> Genome1;
    std::vector<UniqueIdInt> Genome2;
public:
    std::vector<GraphicAtomKind> AtomsKinds;
    std::vector<GraphicParticleKind> ParticlesKinds;
    std::unordered_map<UnsignedInt, GraphicParticleKind> ParticlesKindsXML;
    std::unordered_map<UnsignedInt, UnsignedInt> ParticlesKindsPos;
public:
    std::vector<ParticleKind> ParticlesKinds1;
    UnsignedInt MaxParticleIndex{};
public:
#ifdef PARTICLES_IN_VECTOR
    std::vector<Particle> Particles;
#else
    std::unordered_map<UniqueIdInt, Particle> Particles;
#endif
private:
    std::vector<Reaction> Reactions;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    void AddParticleKind(const ParticleKind& ParticleParam)
    {
        ParticlesKinds1.emplace_back(ParticleParam);
    }

    UniqueIdInt AddNewParticle(UniqueIdInt ParticleIndex, const Particle& ParticleParam)
    {
        #ifdef PARTICLES_IN_VECTOR
        Particles.emplace_back(ParticleParam);
        return MaxParticleIndex = Particles.size() - 1;
        #else
        Particles[ParticleIndex] = ParticleParam;
        return MaxParticleIndex = ParticleIndex;
        #endif
    }

    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return Particles[ParticleIndex];
    }

    void AddReaction(const Reaction& ReactionParam)
    {
        try
        {
            Reactions.emplace_back(ReactionParam);
            ReactionsIdByString.insert(make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
        }
        CATCH("adding reaction")
    }

    std::vector<GraphicAtomKind>::iterator GetAtomKindDataForAtom(char Name)
    {
        std::vector<GraphicAtomKind>::iterator AtomKindObjectIterator;

        try
        {
            AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, Name));
            if (AtomKindObjectIterator == AtomsKinds.end())
                AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, 'E'));
        }
        CATCH("getting atom kind data for atom")

        return AtomKindObjectIterator;
    }
};

inline CellEngineSimulationManager CellEngineSimulationManagerObject;

#endif
