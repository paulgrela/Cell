
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
#include "CellEngineColors.h"

template <class T>
class PairedNucleotide
{
public:
    T* PairedNucleotide = nullptr;
public:
    UniqueIdInt PairedNucleotideTemporary = 0;
public:
    static void LinkPairedNucleotides(T* PairedNucleotide1, T* PairedNucleotide2)
    {
        PairedNucleotide1->PairedNucleotide = PairedNucleotide2;
        PairedNucleotide2->PairedNucleotide = PairedNucleotide1;
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
public:
    void SetCenterCoordinates(PositionInt XCenterParam, PositionInt YCenterParam, PositionInt ZCenterParam)
    {
        Center.X = XCenterParam;
        Center.Y = YCenterParam;
        Center.Z = ZCenterParam;
    }
public:
    explicit Particle(std::vector<vector3_16>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {
    }
    explicit Particle(UniqueIdInt IndexParam, EntityIdInt EntityIdParam, ChainIdInt ChainIdParam, UniqueIdInt GenomeThreadParam, UniqueIdInt GenomeIndexParam, ElectricChargeType ElectricChargeParam, vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), ElectricCharge(ElectricChargeParam), UniqueColor(UniqueColorParam)
    {
    }
public:
    Particle() = default;
};

inline double DistanceOfParticles(const Particle& Particle1, const Particle& Particle2)
{
    return sqrt(pow(static_cast<double>(Particle1.Center.X) - static_cast<double>(Particle2.Center.X), 2.0) + pow(static_cast<double>(Particle1.Center.Y) - static_cast<double>(Particle2.Center.Y), 2.0) + pow(static_cast<double>(Particle1.Center.Z) - static_cast<double>(Particle2.Center.Z), 2.0));
}

struct ParticleKindGraphicData
{
    EntityIdInt EntityId{};
    bool Visible{};
    float SizeX{};
    float SizeY{};
    float SizeZ{};
    vector3_16 AtomColor{};
    vector3_16 ParticleColor{};
    vector3_16 RandomParticleColor{};
    std::string NameFromXML;
    std::string NameFromDataFile;
};

enum class ParticleDestinationTypes : UnsignedInt
{
    Empty = 0,
    Basic = 1,
    Lipid = 2,
    tRNA = 3,
    mRNA = 4,
    rRNA = 5,
    PolymeraseProtein = 6,
    RibosomesProtein = 7,
    MembraneProtein = 8,
    Other = 9
};

class Gene
{
public:
    GeneIdInt NumId;
    std::string StrId;
    std::string Description;
    std::string ProteinId;
    UnsignedInt StartPosInGenome;
    UnsignedInt EndPosInGenome;
    std::string Sequence;
};

class ParticleKind
{
public:
    EntityIdInt EntityId{};
    std::string Id;
    std::string Name;
    std::string Formula;
    GeneIdInt GeneId{};
    UnsignedInt Counter{};
    ElectricChargeType ElectricCharge{};
    std::string Compartment{};
public:
    ParticleDestinationTypes ProteinDestinationType = ParticleDestinationTypes::Empty;
public:
    std::vector<vector3_16> ListOfVoxels;
    UnsignedInt XSizeDiv2{}, YSizeDiv2{}, ZSizeDiv2{};
public:
    ParticleKindGraphicData GraphicData;
public:
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    void SetReactionsIdByString(std::unordered_map<std::string, UnsignedInt>& ReactionsIdByStringParam)
    {
        ReactionsIdByString = ReactionsIdByStringParam;
    }
public:
    ParticleKind() = default;
public:
    ParticleKind(UnsignedInt EntityIdParam, std::string IdParam, std::string NameParam, std::string FormulaParam, GeneIdInt GeneIdParam, ElectricChargeType ElectricChargeParam, std::string CompartmentParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Id(std::move(IdParam)), Name(std::move(NameParam)), Formula(std::move(FormulaParam)), GeneId(GeneIdParam), ElectricCharge(ElectricChargeParam), Compartment(std::move(CompartmentParam)), Counter(CounterParam)
    {
    }
    ParticleKind(UnsignedInt EntityIdParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Counter(CounterParam)
    {
    }
    explicit ParticleKind(const ParticleKindGraphicData& ParticleKindGraphicDataObjectParam) : EntityId(ParticleKindGraphicDataObjectParam.EntityId), GraphicData(ParticleKindGraphicDataObjectParam)
    {
    }
};

class ParticleKindForReaction
{
public:
    EntityIdInt EntityId{};
    UnsignedInt Counter{};
    bool ToRemoveInReaction;
    std::vector<UniqueIdInt> LinkedParticleTypes;
public:
    std::string SequenceStr;
    std::vector<ChainIdInt> Sequence;
public:
    ParticleKindForReaction(EntityIdInt EntityIdParam, UnsignedInt CounterParam, std::string SequenceStrParam, bool ToRemoveInReactionParam) : EntityId(EntityIdParam), Counter(CounterParam), SequenceStr(std::move(SequenceStrParam)), ToRemoveInReaction(ToRemoveInReactionParam)
    {
        for (auto& NucleotideLetter : SequenceStr)
            Sequence.emplace_back(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(NucleotideLetter));
    }
    ParticleKindForReaction(EntityIdInt EntityIdParam, UnsignedInt CounterParam, std::string SequenceStrParam, bool ToRemoveInReactionParam, std::vector<UniqueIdInt> LinkedParticlesTypesParam) : ParticleKindForReaction(EntityIdParam, CounterParam, std::move(SequenceStrParam), ToRemoveInReactionParam)
    {
        LinkedParticleTypes = std::move(LinkedParticlesTypesParam);
    }
};

struct AtomKindGraphicData
{
    std::string Name;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3_16 Color;
    vmath::vec3 ColorVmathVec3;
};

inline bool operator==(const AtomKindGraphicData& AtomKindParameter, const std::string& NameStr)
{
    return AtomKindParameter.Name == NameStr;
}

class ParticlesKindsManager
{
public:
    std::unordered_map<GeneIdInt, Gene> Genes;
    std::unordered_map<EntityIdInt, ParticleKind> ParticlesKinds;
    std::unordered_map<EntityIdInt, ParticleKindGraphicData> GraphicParticlesKindsFromConfigXML;
public:
    std::vector<AtomKindGraphicData> AtomsKindsGraphicData;
public:
    ParticleKind& GetParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds.find(EntityId)->second;
    }
    std::optional<ParticleKind> GetParticleKindFromStrId(const std::string& StrId)
    {
        for (auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.Id == StrId)
                return ParticleKindObject.second;

        return {};
    }
    ParticleKindGraphicData& GetGraphicParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds.find(EntityId)->second.GraphicData;
    }
    void AddParticleKind(const ParticleKind& ParticleKindObjectParam)
    {
        ParticlesKinds[ParticleKindObjectParam.EntityId] = ParticleKindObjectParam;
        ParticlesKinds[ParticleKindObjectParam.EntityId].GraphicData = ParticleKindGraphicData{ParticleKindObjectParam.EntityId, true, 1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectParam.Name, ParticleKindObjectParam.Formula };
    }
public:
    std::vector<AtomKindGraphicData>::iterator GetGraphicAtomKindDataFromAtomName(char Name)
    {
        std::vector<AtomKindGraphicData>::iterator AtomKindObjectIterator;

        try
        {
            AtomKindObjectIterator = std::find(AtomsKindsGraphicData.begin(), AtomsKindsGraphicData.end(), std::string(1, Name));
            if (AtomKindObjectIterator == AtomsKindsGraphicData.end())
                AtomKindObjectIterator = std::find(AtomsKindsGraphicData.begin(), AtomsKindsGraphicData.end(), std::string(1, 'E'));
        }
        CATCH("getting atom kind data for atom")

        return AtomKindObjectIterator;
    }
};

inline ParticlesKindsManager ParticlesKindsManagerObject;

#endif
