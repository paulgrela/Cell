
#ifndef CELL_ENGINE_PARTICLES_KINDS_MANAGER_H
#define CELL_ENGINE_PARTICLES_KINDS_MANAGER_H

#include <optional>

#include "ExceptionsMacro.h"
#include "CellEngineUseful.h"
#include "CellEngineColors.h"
#include "CellEngineParticleKind.h"

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
    std::vector<GeneIdInt> Ribosomes30SProteinsList;
    std::vector<GeneIdInt> Ribosomes50SProteinsList;
public:
    std::vector<AtomKindGraphicData> AtomsKindsGraphicData;
public:
    ParticleKind& GetParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds.find(EntityId)->second;
    }
public:
    std::optional<ParticleKind> GetParticleKindFromStrId(const std::string& StrId)
    {
        for (auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.IdStr == StrId)
                return ParticleKindObject.second;

        return {};
    }
public:
    std::optional<ParticleKind> GetParticleKindFromGeneId(const SignedInt GeneId)
    {
        for (auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.GeneId == GeneId)
                return ParticleKindObject.second;

        for (auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                if (ParticleKindObject.second.ParticleKindSpecialDataSector[0].GeneId == GeneId)
                    return ParticleKindObject.second;

        return {};
    }
public:
    ParticleKindGraphicData& GetGraphicParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds.find(EntityId)->second.GraphicData;
    }
public:
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
public:
    static std::string ConvertParticleTypeToString(ParticlesTypes ParticleType)
    {
        switch (ParticleType)
        {
            case ParticlesTypes::Empty : return "Empty";
            case ParticlesTypes::Basic : return "Basic";
            case ParticlesTypes::Lipid : return "Lipid";
            case ParticlesTypes::mRNA : return "mRNA";
            case ParticlesTypes::rRNA : return "rRNA";
            case ParticlesTypes::tRNA_charged : return "tRNA_charged";
            case ParticlesTypes::tRNA_uncharged : return "tRNA_uncharged";
            case ParticlesTypes::RNAPolymeraseProtein : return "RNAPolymeraseProtein";
            case ParticlesTypes::PolymeraseProtein : return "PolymeraseProtein";
            case ParticlesTypes::RibosomeProtein : return "RibosomeProtein";
            case ParticlesTypes::MembraneProtein : return "MembraneProtein";
            case ParticlesTypes::OtherProtein : return "OtherProtein";
            case ParticlesTypes::ProteinFrac : return "ProteinFrac";
            case ParticlesTypes::Other : return "Other";
            case ParticlesTypes::RNAPolymerase : return "RNAPolymerase";
            case ParticlesTypes::DNAPolymerase : return "DNAPolymerase";
            case ParticlesTypes::Ribosome : return "Ribosomes";
            default : return "NONE";
        }
    };
};

inline ParticlesKindsManager ParticlesKindsManagerObject;

#endif
