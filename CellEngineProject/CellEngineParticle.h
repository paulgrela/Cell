
#ifndef CELL_ENGINE_PARTICLE_H
#define CELL_ENGINE_PARTICLE_H

#include <optional>
#include <algorithm>

#include <list>
#include <utility>
#include <vector>
#include <string>
#include <unordered_map>

#include "DoublyLinkedList.h"

#include "CellEngineAtom.h"
#include "CellEngineTypes.h"
#include "CellEngineUseful.h"
#include "CellEngineConstants.h"

template <class T>
class PairedNucleotide
{
public:
    T* PairedNucleotidePtr = nullptr;
public:
    UniqueIdInt PairedNucleotideTemporary = 0;
public:
    static void LinkPairedNucleotides(T* PairedNucleotide1, T* PairedNucleotide2)
    {
        PairedNucleotide1->PairedNucleotidePtr = PairedNucleotide2;
        PairedNucleotide2->PairedNucleotidePtr = PairedNucleotide1;
    }
};

template <class T>
class LinkedParticles
{
public:
    std::vector<T*> LinkedParticlesPointersList;
public:
    std::vector<UniqueIdInt> LinkedParticlesPointersListTemporary;
public:
    void AddNewLinkToParticle(T* NewLinkToParticle)
    {
        LinkedParticlesPointersList.emplace_back(NewLinkToParticle);
    }
    void RemoveLastLinkToParticle()
    {
        LinkedParticlesPointersList.pop_back();
    }
    void SetLinkToParticle(const UnsignedInt Position, const T* NewLinkToParticle)
    {
        LinkedParticlesPointersList[Position] = NewLinkToParticle;
    }
    void SetNullToLinkToParticle(const UnsignedInt Position)
    {
        LinkedParticlesPointersList[Position] = nullptr;
    }
    void SetNullToAllLinksToParticles()
    {
        fill(begin(LinkedParticlesPointersList), end(LinkedParticlesPointersList), nullptr);
    }
};

using ListOfVoxelsType = std::vector<vector3_16>;
using ListOfAtomsType = std::vector<CellEngineAtom>;
using ListOfCentersType = std::vector<vector3_Real32>;

class Particle : public DoublyLinkedListNode<Particle, UniqueIdInt>, public PairedNucleotide<Particle>, public LinkedParticles<Particle>
{
public:
    bool SelectedForReaction{};
public:
    EntityIdInt EntityId{};
    ChainIdInt ChainId{};
    UniqueIdInt Index{};
    UniqueIdInt GenomeThread{};
    UniqueIdInt GenomeIndex{};
    UniqueIdInt GenomeIndexPrev{};
    UniqueIdInt GenomeIndexNext{};
    ElectricChargeType ElectricCharge{};
    SectorPosType ParticleSectorPos{};
public:
    vector3_Real32 Center{};
    RealType Radius{};
public:
    ListOfVoxelsType ListOfVoxels;
    ListOfAtomsType ListOfAtoms;
    std::vector<std::pair<UnsignedInt, UnsignedInt>> BondsBetweenAtomsToDraw;
public:
    std::string SequenceStr;
    UnsignedInt PositionInSequence{};
public:
    char Nucleotide = '?';
public:
    vector3_16 ParticleColor{};
    vector3_16 UniqueParticleColor{};
    vector3_16 RandomParticleKindColor{};
public:
    UnsignedInt OrderInParticleIndex{};
public:
    void SetCenterCoordinates(const RealType XCenterParam, const RealType YCenterParam, const RealType ZCenterParam)
    {
        Center.X = XCenterParam;
        Center.Y = YCenterParam;
        Center.Z = ZCenterParam;
    }
public:
    explicit Particle(const UniqueIdInt IndexParam, const EntityIdInt EntityIdParam, const ChainIdInt ChainIdParam, const UniqueIdInt GenomeThreadParam, const UniqueIdInt GenomeIndexParam, const ElectricChargeType ElectricChargeParam, const vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), ElectricCharge(ElectricChargeParam), UniqueParticleColor(UniqueColorParam)
    {
    }
    explicit Particle(const UniqueIdInt IndexParam, const EntityIdInt EntityIdParam, const ChainIdInt ChainIdParam, const UniqueIdInt GenomeThreadParam, const UniqueIdInt GenomeIndexParam, const ElectricChargeType ElectricChargeParam, std::string SequenceStrParam, const vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), ElectricCharge(ElectricChargeParam), SequenceStr(std::move(SequenceStrParam)), UniqueParticleColor(UniqueColorParam)
    {
    }
public:
    Particle() = default;
};

inline double DistanceOfParticles(const Particle& Particle1, const Particle& Particle2)
{
    return sqrt(pow(static_cast<double>(Particle1.Center.X) - static_cast<double>(Particle2.Center.X), 2.0) + pow(static_cast<double>(Particle1.Center.Y) - static_cast<double>(Particle2.Center.Y), 2.0) + pow(static_cast<double>(Particle1.Center.Z) - static_cast<double>(Particle2.Center.Z), 2.0));
}

inline double DistanceOfParticleFromPoint(const Particle& Particle1, const vector3_Real32& Point)
{
    return sqrt(pow(static_cast<double>(Particle1.Center.X) - static_cast<double>(Point.X), 2.0) + pow(static_cast<double>(Particle1.Center.Y) - static_cast<double>(Point.Y), 2.0) + pow(static_cast<double>(Particle1.Center.Z) - static_cast<double>(Point.Z), 2.0));
}

inline double DistanceOfPoints(const vector3_Real32& Point1, const vector3_Real32& Point2)
{
    return sqrt(pow(static_cast<double>(Point1.X) - static_cast<double>(Point2.X), 2.0) + pow(static_cast<double>(Point1.Y) - static_cast<double>(Point2.Y), 2.0) + pow(static_cast<double>(Point1.Z) - static_cast<double>(Point2.Z), 2.0));
}

#endif
