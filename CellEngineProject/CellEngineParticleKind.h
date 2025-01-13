
#ifndef CELL_ENGINE_PARTICLE_KIND_H
#define CELL_ENGINE_PARTICLE_KIND_H

#include <set>
#include <string>
#include <vector>
#include <unordered_map>

#include "CellEngineParticle.h"
#include "CellEngineTypes.h"

struct ParticleKindGraphicData
{
    EntityIdInt EntityId{};
    bool Visible{};
    bool Selected{};
    float SizeX{};
    float SizeY{};
    float SizeZ{};
    vector3_16 AtomColor{};
    vector3_16 ParticleColor{};
    vector3_16 RandomParticleColor{};
    std::string NameFromXML;
    std::string NameFromDataFile;
};

enum class ParticlesTypes : UnsignedInt
{
    DNANucleotide = 0,
    Basic = 1,
    Lipid = 2,
    mRNA = 3,
    rRNA = 4,
    tRNA_uncharged = 5,
    tRNA_charged = 6,
    RNAPolymeraseProtein = 7,
    PolymeraseProtein = 8,
    RibosomeProtein = 9,
    MembraneProtein = 10,
    OtherProtein = 11,
    ProteinFrac = 12,
    RNAPolymerase = 13,
    DNAPolymerase = 14,
    Ribosome = 15,
    Other = 16,
    Empty = 17,
    ProteinInBuildingProcess = 18
};

class InterGeneSequence
{
public:
    UnsignedInt StartPosInGenome;
    UnsignedInt EndPosInGenome;
    std::string Sequence;
public:
    InterGeneSequence(const UnsignedInt StartPosInGenomeParam, const UnsignedInt EndPosInGenomeParam, std::string SequenceParam) : StartPosInGenome(StartPosInGenomeParam), EndPosInGenome(EndPosInGenomeParam), Sequence(std::move(SequenceParam))
    {
    }
};

class Promoter
{
public:
    GeneIdInt GeneId;
    UnsignedInt Box10Position;
    UnsignedInt Box35Position;
    UnsignedInt StartCodonPosition;
};

class Terminator
{
public:
    std::string LeftStem;
    std::string RightStem;
    UnsignedInt GeneEnd;
    UnsignedInt PosInGenomeHairpin;
    SignedInt LengthOfRunning;
};

class Gene
{
public:
    GeneIdInt NumId;
    std::string StrId;
    std::string Description;
    std::string ProteinId;
    std::string Sequence;
    UnsignedInt StartPosInGenome;
    UnsignedInt EndPosInGenome;
};

class ParticleKindSpecialData
{
public:
    SignedInt GeneId{ -1 };
    std::string Description;
    std::string AddedParticle;
    SignedInt CleanProductOfTranscription{ false };
    ParticlesTypes ParticleType = ParticlesTypes::Empty;
    bool IsProtein{ false };
    UnsignedInt CounterAtStartOfSimulation{ 0 };
};

class ParticleKind
{
public:
    EntityIdInt EntityId{};
    std::string IdStr;
    std::string Name;
    std::string Formula;
    std::string Compartment{};
    ElectricChargeType ElectricCharge{};
public:
    std::set<UnsignedInt> AssociatedChemicalReactions;
    std::vector<ParticleKindSpecialData> ParticleKindSpecialDataSector;
    GeneIdInt GeneId{};
public:
    UnsignedInt Counter{};
public:
    //std::vector<vector3_16> ListOfVoxels;
    ListOfElements ListOfVoxels;
    UnsignedInt XSizeDiv2{}, YSizeDiv2{}, ZSizeDiv2{};
public:
    ParticleKindGraphicData GraphicData;
public:
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
public:
    void SetReactionsIdByString(const std::unordered_map<std::string, UnsignedInt>& ReactionsIdByStringParam)
    {
        ReactionsIdByString = ReactionsIdByStringParam;
    }
public:
    ParticleKind() = default;
public:
    ParticleKind(const UnsignedInt EntityIdParam, std::string IdStrParam, std::string NameParam, std::string FormulaParam, const GeneIdInt GeneIdParam, ElectricChargeType const ElectricChargeParam, std::string CompartmentParam, const UnsignedInt CounterParam) : EntityId(EntityIdParam), IdStr(std::move(IdStrParam)), Name(std::move(NameParam)), Formula(std::move(FormulaParam)), GeneId(GeneIdParam), ElectricCharge(ElectricChargeParam), Compartment(std::move(CompartmentParam)), Counter(CounterParam)
    {
    }
    ParticleKind(const UnsignedInt EntityIdParam, const UnsignedInt CounterParam) : EntityId(EntityIdParam), Counter(CounterParam)
    {
    }
    explicit ParticleKind(const ParticleKindGraphicData& ParticleKindGraphicDataObjectParam) : EntityId(ParticleKindGraphicDataObjectParam.EntityId), GraphicData(ParticleKindGraphicDataObjectParam)
    {
    }
};

#endif
