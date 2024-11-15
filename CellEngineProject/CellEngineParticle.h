
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
    vector3_16 UniqueColor{};
    ElectricChargeType ElectricCharge{};
public:
    std::vector<vector3_16> ListOfVoxels;
    vector3_16 Center{};
public:
    std::string SequenceStr;
    UnsignedInt PositionInSequence{};
public:
    void SetCenterCoordinates(const PositionInt XCenterParam, const PositionInt YCenterParam, const PositionInt ZCenterParam)
    {
        if (XCenterParam < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
            Center.X = XCenterParam;
        if (YCenterParam < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
            Center.Y = YCenterParam;
        if (ZCenterParam < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
            Center.Z = ZCenterParam;
    }
public:
    explicit Particle(std::vector<vector3_16>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {
    }
    explicit Particle(const UniqueIdInt IndexParam, const EntityIdInt EntityIdParam, const ChainIdInt ChainIdParam, const UniqueIdInt GenomeThreadParam, const UniqueIdInt GenomeIndexParam, const ElectricChargeType ElectricChargeParam, const vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), ElectricCharge(ElectricChargeParam), UniqueColor(UniqueColorParam)
    {
    }
    explicit Particle(const UniqueIdInt IndexParam, const EntityIdInt EntityIdParam, const ChainIdInt ChainIdParam, const UniqueIdInt GenomeThreadParam, const UniqueIdInt GenomeIndexParam, const ElectricChargeType ElectricChargeParam, std::string SequenceStrParam, const vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), ElectricCharge(ElectricChargeParam), SequenceStr(std::move(SequenceStrParam)), UniqueColor(UniqueColorParam)
    {
    }
public:
    Particle() = default;
};

inline double DistanceOfParticles(const Particle& Particle1, const Particle& Particle2)
{
    return sqrt(pow(static_cast<double>(Particle1.Center.X) - static_cast<double>(Particle2.Center.X), 2.0) + pow(static_cast<double>(Particle1.Center.Y) - static_cast<double>(Particle2.Center.Y), 2.0) + pow(static_cast<double>(Particle1.Center.Z) - static_cast<double>(Particle2.Center.Z), 2.0));
}

#endif
