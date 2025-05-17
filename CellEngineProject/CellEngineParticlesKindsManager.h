
#ifndef CELL_ENGINE_PARTICLES_KINDS_MANAGER_H
#define CELL_ENGINE_PARTICLES_KINDS_MANAGER_H

#include <string>
#include <optional>

#include "ExceptionsMacro.h"
#include "CellEngineUseful.h"
#include "CellEngineParticleKind.h"

#ifndef USING_MODULES
#include "CellEngineColors.h"
#endif

struct AtomKindGraphicData
{
    std::string Name;
    RealType SizeX;
    RealType SizeY;
    RealType SizeZ;
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
    MainMapType<GeneIdInt, Gene> Genes;
    MainMapType<UnsignedInt, Promoter> Promoters;
    MainMapType<std::string, Terminator> Terminators;
    MainMapType<EntityIdInt, ParticleKind> ParticlesKinds;
    MainMapType<EntityIdInt, ParticleKindGraphicData> GraphicParticlesKindsFromConfigXML;
public:
    std::vector<InterGeneSequence> InterGenesSequences;
public:
    std::vector<GeneIdInt> Ribosomes30SProteinsList;
    std::vector<GeneIdInt> Ribosomes50SProteinsList;
public:
    std::vector<AtomKindGraphicData> AtomsKindsGraphicData;
public:
    Promoter& GetPromoterFromBox10Position(const UnsignedInt Box10Position)
    {
        return Promoters.find(Box10Position)->second;
    }
public:
    ParticleKind& GetParticleKind(const EntityIdInt EntityId)
    {
        return ParticlesKinds.find(EntityId)->second;
    }
public:
    std::optional<ParticleKind> GetParticleKindFromStrId(const std::string& StrId) const
    {
        for (const auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.IdStr == StrId)
                return ParticleKindObject.second;

        return {};
    }
public:
    std::optional<ParticleKind> GetParticleKindFromGeneId(const SignedInt GeneId) const
    {
        for (const auto& ParticleKindObject : ParticlesKinds)
            if (ParticleKindObject.second.GeneId == GeneId)
                return ParticleKindObject.second;

        for (const auto& ParticleKindObject : ParticlesKinds)
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
    vector3_16 GetParticleKindGraphicDataFromConfigXMLData(const EntityIdInt EntityId)
    {
        if (auto ParticleKindObjectIterator = GraphicParticlesKindsFromConfigXML.find(EntityId); ParticleKindObjectIterator == GraphicParticlesKindsFromConfigXML.end())
            return GraphicParticlesKindsFromConfigXML.find(10000)->second.ParticleColor;
        else
            return ParticleKindObjectIterator->second.ParticleColor;
    }
public:
    void AddParticleKind(const ParticleKind& ParticleKindObjectParam)
    {
        ParticlesKinds[ParticleKindObjectParam.EntityId] = ParticleKindObjectParam;
        #ifndef USING_MODULES
        ParticlesKinds[ParticleKindObjectParam.EntityId].GraphicData = ParticleKindGraphicData{ParticleKindObjectParam.EntityId, true, false, 1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), ParticleKindObjectParam.Name, ParticleKindObjectParam.Formula };
        #endif
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
    void PrintAllParticleKinds() const
    {
        try
        {
            for (const auto& ParticleKindObject : ParticlesKinds)
            {
                LoggersManagerObject.Log(STREAM("KIND PARTICLE = " << ParticleKindObject.second.EntityId << " " << ParticleKindObject.second.IdStr << " " << ParticleKindObject.second.Name << " " << ParticleKindObject.second.Formula << " " << ParticleKindObject.second.ElectricCharge << " GeneId = " << ParticleKindObject.second.GeneId));
                for (const auto& ParticleKindParticleKindSpecialDataSectorObject : ParticleKindObject.second.ParticleKindSpecialDataSector)
                    LoggersManagerObject.Log(STREAM(std::string("KIND GENE = " + std::string(ParticleKindParticleKindSpecialDataSectorObject.GeneId != -1 ? JCVISYN3APredStr + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.GeneId) : "NoGene") + " TYPE = " + ParticlesKindsManager::ConvertParticleTypeToString(ParticleKindParticleKindSpecialDataSectorObject.ParticleType) + " D = #" + ParticleKindParticleKindSpecialDataSectorObject.Description + "# Added = #" + ParticleKindParticleKindSpecialDataSectorObject.AddedParticle + "# CLEAN PRODUCT = #" + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.CleanProductOfTranscription) + "# COUNTER = " + std::to_string(ParticleKindParticleKindSpecialDataSectorObject.CounterAtStartOfSimulation))));                ;

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
            case ParticlesTypes::Ribosome : return "Ribosome";
            default : return "NONE";
        }
    };
};

inline ParticlesKindsManager ParticlesKindsManagerObject;

#endif
