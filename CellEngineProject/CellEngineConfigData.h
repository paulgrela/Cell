
#ifndef CELL_ENGINE_CONFIG_DATA_H
#define CELL_ENGINE_CONFIG_DATA_H

#include "vmath.h"
#include <random>

#include "sb7color.h"

#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

struct GraphicParticleKind
{
    UnsignedInt Identifier;
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3_16 AtomColor;
    vector3_16 ParticleColor;
    vector3_16 RandomParticleColor;
    std::string NameFromXML;
    std::string NameFromDataFile;
};

struct GraphicAtomKind
{
    std::string Name;
    float SizeX;
    float SizeY;
    float SizeZ;
    vector3_16 Color;
    vmath::vec3 ColorVmathVec3;
};

inline bool operator==(const GraphicAtomKind& AtomKindParameter, const std::string& NameStr)
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
    UnsignedInt NumberOfVoxelsInVoxelSimulationSpaceInEachDimension{};
    float DivisionFactorForVoxelSimulationSpace{};
public:
    UnsignedInt VoxelSimulationSpaceSelectionStartXPos{}, VoxelSimulationSpaceSelectionStartYPos{}, VoxelSimulationSpaceSelectionStartZPos{};
    UnsignedInt VoxelSimulationSpaceSelectionStepX{}, VoxelSimulationSpaceSelectionStepY{}, VoxelSimulationSpaceSelectionStepZ{};
    UnsignedInt VoxelSimulationSpaceSelectionSizeX{}, VoxelSimulationSpaceSelectionSizeY{}, VoxelSimulationSpaceSelectionSizeZ{};
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
        DrawRandomColorForEveryParticleKind = 3,
        DrawRandomColorForEveryUniqueParticle = 4
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
    vmath::vec3 BackgroundColors[4] = { vmath::FromVec4ToVec3(sb7::color::White), vmath::FromVec4ToVec3(sb7::color::White), vmath::FromVec4ToVec3(sb7::color::White), vmath::FromVec4ToVec3(sb7::color::White) };
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
    std::vector<GraphicAtomKind> AtomsKinds;
    std::vector<GraphicParticleKind> ParticlesKinds;
    std::unordered_map<UnsignedInt, GraphicParticleKind> ParticlesKindsXML;
    std::unordered_map<UnsignedInt, UnsignedInt> ParticlesKindsPos;
public:
    bool ImGuiDemoWindowMenu = false;
    bool ImGuiLightVersion = false;
public:
    int ChosenShapeOfAtoms = 1;
    UnsignedInt DNAIdentifier{};
    UnsignedInt RNAIdentifier{};
public:
    std::vector<GraphicAtomKind>::iterator GetAtomKindDataForAtom(char Name)
    {
        std::vector<GraphicAtomKind>::iterator AtomKindObjectIterator;

        try
        {
            AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, Name));
            if (AtomKindObjectIterator == AtomsKinds.end())
                AtomKindObjectIterator = std::find(AtomsKinds.begin(), AtomsKinds.end(), std::string(1, 'E'));
        }
        CATCH("getting atom kind data for atom")

        return AtomKindObjectIterator;
    }
    inline bool IsDNAorRNA(const EntityIdInt EntityId) const
    {
        return EntityId == DNAIdentifier || EntityId == RNAIdentifier;
    }
    static bool IsNucleotide(const std::string_view ChainName)
    {
        return (ChainName == "NU01" || ChainName == "NU11" || ChainName == "NU02" || ChainName == "NU12" || ChainName == "NU03" || ChainName == "NU13" || ChainName == "NU04" || ChainName == "NU14");
    }
};

inline CellEngineConfigData CellEngineConfigDataObject;

#endif
