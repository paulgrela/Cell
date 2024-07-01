
#ifndef CELL_ENGINE_PARTICLES_KINDS_MANAGER_H
#define CELL_ENGINE_PARTICLES_KINDS_MANAGER_H

#include <string>
#include <optional>

#include "ExceptionsMacro.h"
#include "CellEngineUseful.h"
#include "CellEngineColors.h"
#include "CellEngineConstants.h"
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
        ParticlesKinds[ParticleKindObjectParam.EntityId].GraphicData = ParticleKindGraphicData{ParticleKindObjectParam.EntityId, true, false, 1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectParam.Name, ParticleKindObjectParam.Formula };
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
    void AddSingleParticleKind(const ParticlesTypes ParticlesTypesObject, EntityIdInt& ParticleKindIdParam, const std::string& IdStrParam, const std::string& NameParam, const std::string& FormulaParam, const SignedInt GeneIdParam, const ElectricChargeType ElectricChargeParam, const std::string& CompartmentParam, const UnsignedInt CounterParam)
    {
        try
        {
            AddParticleKind({ ParticleKindIdParam, IdStrParam, NameParam, FormulaParam, static_cast<UnsignedInt>(GeneIdParam == -1 ? 0 : GeneIdParam), ElectricChargeParam, CompartmentParam, CounterParam });
            GetParticleKind(ParticleKindIdParam).ParticleKindSpecialDataSector.emplace_back(ParticleKindSpecialData{ GeneIdParam, "", "", false, ParticlesTypesObject, false, CounterParam });

            ParticleKindIdParam++;

            LoggersManagerObject.Log(STREAM("PARTICLE ADDED"));
        }
        CATCH("adding single particle kind")
    }
public:
    void PrintAllParticleKinds()
    {
        try
        {
            for (const auto& ParticleKindObject : ParticlesKinds)
            {
                LoggersManagerObject.Log(STREAM("KIND PARTICLE = " << ParticleKindObject.second.EntityId << " " << ParticleKindObject.second.IdStr << " " << ParticleKindObject.second.Name << " " << ParticleKindObject.second.Formula << " " << ParticleKindObject.second.ElectricCharge << " GeneId = " << ParticleKindObject.second.GeneId));
                for (const auto& ParticleKindParticleKindSpecialDataSectorObject : ParticleKindObject.second.ParticleKindSpecialDataSector)
                    LoggersManagerObject.Log(STREAM(std::string("KIND GENE = " + std::string(ParticleKindParticleKindSpecialDataSectorObject.GeneId != -1 ? "JCVISYN3A_" + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.GeneId) : "NoGene") + " TYPE = " + ParticlesKindsManager::ConvertParticleTypeToString(ParticleKindParticleKindSpecialDataSectorObject.ParticleType) + " D = #" + ParticleKindParticleKindSpecialDataSectorObject.Description + "# Added = #" + ParticleKindParticleKindSpecialDataSectorObject.AddedParticle + "# CLEAN PRODUCT = #" + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.CleanProductOfTranscription) + "# COUNTER = " + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation))));                ;

                LoggersManagerObject.Log(STREAM(""));
            }
        }
        CATCH("printing all particle kinds")
    };
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
