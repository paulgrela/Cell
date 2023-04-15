
#ifndef CELL_ENGINE_CONFIG_DATA_H
#define CELL_ENGINE_CONFIG_DATA_H

#include "vmath.h"
#include <random>

#include "sb7color.h"

#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

struct ParticleKind
{
    UnsignedInt Identifier;
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3 AtomColor;
    vector3 ParticleColor;
    vector3 RandomParticleColor;
    std::string NameFromXML;
    std::string NameFromDataFile;
};

struct AtomKind
{
    std::string Name;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3 Color;
    vmath::vec3 ColorVmathVec3;
};

inline bool operator==(const AtomKind& AtomKindParameter, const std::string& NameStr)
{
    return AtomKindParameter.Name == NameStr;
}

class CellEngineConfigData
{
public:
    int XTopMainWindow{};
    int YTopMainWindow{};
    int WidthMainWindow{};
    int HeightMainWindow{};
    int XTopMenuWindow{};
    int YTopMenuWindow{};
    int WidthMenuWindow{};
    int HeightMenuWindow{};
    int XTopSecondWindow{};
    int YTopSecondWindow{};
    int WidthSecondWindow{};
    int HeightSecondWindow{};
public:
    bool PrintLogToConsole = true;
    bool PrintLogToFiles = true;

    bool PrintLogLineNumberToConsole = false;
    bool PrintLogDateTimeToConsole = false;
    bool PrintLogProcessIdToConsole = false;
    bool PrintLogProcessPriorityLevelToConsole = false;
    bool PrintLogThreadIdToConsole = false;

    bool PrintLogLineNumberToFile = true;
    bool PrintLogDateTimeToFile = true;
    bool PrintLogProcessIdToFile = false;
    bool PrintLogProcessPriorityLevelToFile = false;
    bool PrintLogThreadIdToFile = false;

    UnsignedInt MaximalNumberOfLinesInOneFile = 100000;
public:
    bool MultiThreaded{};
    bool SetProcessPriorityHighest{};
public:
    enum class TypesOfSpace : UnsignedInt
    {
        FullAtomSpace = 1,
        VoxelSimulationSpace = 2,
    };
    TypesOfSpace TypeOfSpace = TypesOfSpace::FullAtomSpace;
public:
    UnsignedInt NumberOfVoxelSimulationSpaceInEachDimension;
public:
    enum class RandomColorEngineTypes : UnsignedInt
    {
        Rand = 1,
        mt19937 = 2,
        DefaultRandomEngine = 3
    };
    RandomColorEngineTypes RandomColorEngineObject = RandomColorEngineTypes::Rand;
public:
    std::string CellStateFileName;
public:
    float SpecularPower = 8.0f;
    float SpecularAlbedo = 0.2222f;
public:
    enum class MakeColorsType : UnsignedInt
    {
        DrawColorForEveryAtom = 1,
        DrawColorForEveryParticle = 2,
        DrawRandomColorForEveryParticle = 3
    };
    MakeColorsType MakeColorsTypeObject = MakeColorsType::DrawColorForEveryParticle;
public:
    enum class StencilForDrawingObjectsTypes : UnsignedInt
    {
        StencilForDrawingOnlyParticlesCenters = 1,
        StencilForDrawingOnlyInAtomScale = 2
    };
    StencilForDrawingObjectsTypes StencilForDrawingObjectsTypesObject = StencilForDrawingObjectsTypes::StencilForDrawingOnlyInAtomScale;
    UnsignedInt NumberOfStencilBufferLoops{};
    bool UseStencilBuffer = false;
public:
    bool DrawBondsBetweenParticlesCenters{};
    bool DrawBondsBetweenAtoms{};
public:
    bool ShowDetailsInAtomScale = false;
    bool CheckAtomVisibility{};
    float CutZ{};
    float Distance{};
    UnsignedInt LoadOfAtomsStep{};
    bool ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside = false;
public:
    float XLowToDrawInAtomScale{};
    float XHighToDrawInAtomScale{};
    float YLowToDrawInAtomScale{};
    float YHighToDrawInAtomScale{};
    float ZLowToDrawInAtomScale{};
    float ZHighToDrawInAtomScale{};
public:
    enum class SizeOfAtomsDrawingTypes
    {
        AtomSize = 1,
        ParticleSize = 2,
        AutomaticChangeSize = 3
    };
    SizeOfAtomsDrawingTypes SizeOfAtomsDrawingTypesObject = SizeOfAtomsDrawingTypes::AutomaticChangeSize;
public:
    float SizeOfAtomChangeStep{};
    float SizeOfAtomX{};
    float SizeOfAtomY{};
    float SizeOfAtomZ{};
    float RotationAngle1{};
    float RotationAngle2{};
    float RotationAngle3{};
    float CameraXPosition{};
    float CameraYPosition{};
    float CameraZPosition{};
    float CameraXMoveShortStep{};
    float CameraYMoveShortStep{};
    float CameraZMoveShortStep{};
    float CameraXMoveLongStep{};
    float CameraYMoveLongStep{};
    float CameraZMoveLongStep{};
    float ViewPositionX{};
    float ViewPositionY{};
    float ViewPositionZ{};
    float ViewXMoveShortStep{};
    float ViewYMoveShortStep{};
    float ViewZMoveShortStep{};
    float ViewXMoveLongStep{};
    float ViewYMoveLongStep{};
    float ViewZMoveLongStep{};
public:
    bool ViewChangeUsingLongStep{};
    bool AutomaticChangeOfSizeOfAtom{};
    bool AutomaticChangeOfLoadAtomsStep{};
public:
    UnsignedInt ChosenStructureIndex = 0;
public:
    vmath::vec3 BackgroundColors[4] = { sb7::FromVec4ToVec3(sb7::color::White), sb7::FromVec4ToVec3(sb7::color::White), sb7::FromVec4ToVec3(sb7::color::White), sb7::FromVec4ToVec3(sb7::color::White) };
    UnsignedInt ChosenBackgroundColor{};
public:
    bool PrintAtomDescriptionOnScreen{};
    bool LogParametersOfRenderingToFile{};
public:
    std::string AtomDescriptionStr1;
    std::string AtomDescriptionStr2;
    std::string AtomDescriptionStr3;
    std::string AtomDescriptionStr4;
    std::string TimeParametersOfRenderingStr;
    std::string NumberOfRenderedAtomsParametersOfRenderingStr;
public:
    std::vector<AtomKind> AtomsKinds;
    std::vector<ParticleKind> ParticlesKinds;
    std::unordered_map<UnsignedInt, ParticleKind> ParticlesKindsXML;
    std::unordered_map<UnsignedInt, UnsignedInt> ParticlesKindsPos;
public:
    bool ImGuiDemoWindowMenu = false;
    bool ImGuiLightVersion = false;
public:
    int ChosenShapeOfAtoms = 1;
    UnsignedInt DNAIdentifier{};
    UnsignedInt RNAIdentifier{};
public:
    std::vector<AtomKind>::iterator GetAtomKindDataForAtom(char Name)
    {
        std::vector<AtomKind>::iterator AtomKindObjectIterator;

        try
        {
            AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, Name));
            if (AtomKindObjectIterator == AtomsKinds.end())
                AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, 'E'));
        }
        CATCH("getting atom kind data for atom")

        return AtomKindObjectIterator;
    }
public:
    std::mt19937_64 mt64{ std::random_device{}() };
    std::default_random_engine DefaultRandomEngineObject{ static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()) };
public:
    std::uniform_real_distribution<float> UniformDistributionObject;
public:
    void SelectRandomEngineForColors()
    {
        try
        {
            switch (RandomColorEngineObject)
            {
                case RandomColorEngineTypes::Rand : srand((unsigned int)time(nullptr)); break;
                case RandomColorEngineTypes::mt19937 : mt64.seed(std::random_device{}()); break;
                case RandomColorEngineTypes::DefaultRandomEngine: DefaultRandomEngineObject.seed(static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count())); break;
                default:  break;
            }
        }
        CATCH("selecting random engines for colors");
    }
public:
    vmath::vec3 GetRandomColor()
    {
        vmath::vec3 ReturnRandomColor;

        try
        {
            switch (RandomColorEngineObject)
            {
                case CellEngineConfigData::RandomColorEngineTypes::Rand : ReturnRandomColor = vmath::vec3((float)rand() / static_cast<float>(RAND_MAX), (float)rand() / static_cast<float>(RAND_MAX), (float)rand() / static_cast<float>(RAND_MAX)); break;
                case CellEngineConfigData::RandomColorEngineTypes::mt19937 : ReturnRandomColor = vmath::vec3(UniformDistributionObject(mt64), UniformDistributionObject(mt64), UniformDistributionObject(mt64)); break;
                case CellEngineConfigData::RandomColorEngineTypes::DefaultRandomEngine : ReturnRandomColor = vmath::vec3(UniformDistributionObject(DefaultRandomEngineObject), UniformDistributionObject(DefaultRandomEngineObject), UniformDistributionObject(DefaultRandomEngineObject));break;
                default:  break;
            }
        }
        CATCH("getting random color")

        return ReturnRandomColor;
    }
public:
    inline bool IsDNAorRNA(const EntityIdInt EntityId) const
    {
        return EntityId == DNAIdentifier || EntityId == RNAIdentifier;
    }
public:
    vector3 DNANR1Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 DNANR2Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 DNANR3Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 DNANR4Color = GetVector3FormVMathVec3(GetRandomColor());

    vector3 RNANR1Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 RNANR2Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 RNANR3Color = GetVector3FormVMathVec3(GetRandomColor());
    vector3 RNANR4Color = GetVector3FormVMathVec3(GetRandomColor());

    vector3 GetDNAorRNAColor(EntityIdInt EntityId, ChainIdInt ChainId) const
    {
        if (EntityId == DNAIdentifier)
        {
            switch (ChainId)
            {
                case 01:
                case 11: return DNANR1Color; break;
                case 02:
                case 12: return DNANR2Color; break;
                case 03:
                case 13: return DNANR3Color; break;
                case 04:
                case 14: return DNANR4Color; break;
                default : break;
            }
        }
        else
        if (EntityId == RNAIdentifier)
        {
            switch (ChainId)
            {
                case 01: return RNANR1Color; break;
                case 02: return RNANR2Color; break;
                case 03: return RNANR3Color; break;
                case 04: return RNANR4Color; break;
                default : break;
            }
        }

        return { 0, 0, 0 };
    }

};

inline CellEngineConfigData CellEngineConfigDataObject;

#endif
