
#ifndef CELL_ENGINE_PARTICLE_H
#define CELL_ENGINE_PARTICLE_H

#include <list>
#include <vector>
#include <string>
#include <unordered_map>

#include "CellEngineTypes.h"
#include "CellEngineUseful.h"
#include "CellEngineColors.h"

class Particle
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
    ElectricChargeIdInt ElectricCharge{};
public:
    std::vector<vector3_16> ListOfVoxels;
    UnsignedInt XCenter{}, YCenter{}, ZCenter{};
public:
    void SetCenterCoordinates(UnsignedInt XCenterParam, UnsignedInt YCenterParam, UnsignedInt ZCenterParam)
    {
        XCenter = XCenterParam;
        YCenter = YCenterParam;
        ZCenter = ZCenterParam;
    }
public:
    explicit Particle(std::vector<vector3_16>& ListOfVoxelsParam) : ListOfVoxels(std::move(ListOfVoxelsParam)), SelectedForReaction(false)
    {
    }
    explicit Particle(UniqueIdInt IndexParam, EntityIdInt EntityIdParam, ChainIdInt ChainIdParam, UniqueIdInt GenomeThreadParam, UniqueIdInt GenomeIndexParam, vector3_16 UniqueColorParam) : Index(IndexParam), EntityId(EntityIdParam), ChainId(ChainIdParam), GenomeThread(GenomeThreadParam), GenomeIndex(GenomeIndexParam), UniqueColor(UniqueColorParam)
    {
    }
public:
    Particle() = default;
};

struct ParticleKindGraphicData
{
    UnsignedInt EntityId{};
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

class ParticleKind
{
public:
    UnsignedInt EntityId{};
    std::string Name;
    std::string Symbol;
    UnsignedInt Counter{};
    ElectricChargeIdInt ElectricCharge;
public:
    std::string SelfSequenceStr;
    std::vector<ChainIdInt> SelfSequence;
    std::string SequenceStr;
    std::vector<ChainIdInt> Sequence;
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
    ParticleKind(UnsignedInt EntityIdParam, std::string NameParam, std::string SymbolParam, UnsignedInt CounterParam) : EntityId(EntityIdParam), Name(std::move(NameParam)), Symbol(std::move(SymbolParam)), Counter(CounterParam)
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
    UniqueIdInt EntityId{};
    UnsignedInt Counter{};
    bool ToRemoveInReaction;
public:
    ParticleKindForReaction(UnsignedInt EntityIdParam, UnsignedInt CounterParam, bool ToRemoveInReactionParam) : EntityId(EntityIdParam), Counter(CounterParam), ToRemoveInReaction(ToRemoveInReactionParam)
    {
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
    std::vector<ParticleKind> ParticlesKinds;
    std::unordered_map<UnsignedInt, UnsignedInt> ParticlesKindsPos;
    std::unordered_map<UnsignedInt, ParticleKindGraphicData> GraphicParticlesKindsFromConfigXML;
public:
    std::vector<AtomKindGraphicData> AtomsKindsGraphicData;
public:
    ParticleKind& GetParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds[ParticlesKindsPos.find(EntityId)->second];
    }
    ParticleKindGraphicData& GetGraphicParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds[ParticlesKindsPos.find(EntityId)->second].GraphicData;
    }
    void AddParticleKind(const ParticleKind& ParticleParam)
    {
        ParticlesKinds.emplace_back(ParticleParam);
        ParticlesKinds.back().GraphicData = ParticleKindGraphicData{ ParticleParam.EntityId, true, 1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleParam.Name, ParticleParam.Symbol };
        ParticlesKindsPos[ParticleParam.EntityId] = ParticlesKinds.size() - 1;
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
