
#ifndef CELL_ENGINE_PARTICLE_KIND_H
#define CELL_ENGINE_PARTICLE_KIND_H

#include <string>
#include <vector>
#include <unordered_map>

#include "CellEngineTypes.h"

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
    Empty = 17
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
    std::vector<ParticleKindSpecialData> ParticleKindSpecialDataSector;
    GeneIdInt GeneId{};
public:
    UnsignedInt Counter{};
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
    ParticleKind(UnsignedInt EntityIdParam, std::string IdStrParam, std::string NameParam, std::string FormulaParam, GeneIdInt GeneIdParam, ElectricChargeType ElectricChargeParam, std::string CompartmentParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), IdStr(std::move(IdStrParam)), Name(std::move(NameParam)), Formula(std::move(FormulaParam)), GeneId(GeneIdParam), ElectricCharge(ElectricChargeParam), Compartment(std::move(CompartmentParam)), Counter(CounterParam)
    {
    }
    ParticleKind(UnsignedInt EntityIdParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Counter(CounterParam)
    {
    }
    explicit ParticleKind(const ParticleKindGraphicData& ParticleKindGraphicDataObjectParam) : EntityId(ParticleKindGraphicDataObjectParam.EntityId), GraphicData(ParticleKindGraphicDataObjectParam)
    {
    }
};

#endif
