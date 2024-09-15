
#ifndef CELL_ENGINE_CONFIG_DATA_H
#define CELL_ENGINE_CONFIG_DATA_H

#include "vmath.h"
#include <map>
#include <random>

#include "sb7color.h"

#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

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

    bool PrintLogToCommonFileWhenPrintLogToSpecialFile = true;
public:
    bool PrintAtomDescriptionOnScreen{};
    bool LogParametersOfRenderingToFile{};
public:
    std::string TimeParametersOfRenderingStr;
    std::string NumberOfRenderedAtomsParametersOfRenderingStr;
public:
    bool GenomeReadFromFile = false;
public:
    bool MultiThreaded{};
    bool SetProcessPriorityHighest{};
public:
    enum class TypesOfFileToRead : UnsignedInt
    {
        BinaryFile = 0,
        CIFFile = 1,
        PDBFile = 2
    };
    TypesOfFileToRead TypeOfFileToRead = TypesOfFileToRead::BinaryFile;
public:
    enum class TypesOfSpace : UnsignedInt
    {
        FullAtomSpace = 1,
        VoxelSimulationSpace = 2
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
    double RadiusOfCellForDNA{};
public:
    enum class SelectedSpaceStartParametersDrawTypes : UnsignedInt
    {
        DrawFromCenter = 1,
        DrawFromCorner = 2
    };
    SelectedSpaceStartParametersDrawTypes SelectedSpaceStartParametersDrawTypesObject =  SelectedSpaceStartParametersDrawTypes::DrawFromCenter;
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
    std::string CellStateFileNameBackup;
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
    bool ImGuiDemoWindowMenu = false;
    bool ImGuiLightVersion = false;
public:
    int ChosenShapeOfAtoms = 1;
    EntityIdInt DNAIdentifier{};
    EntityIdInt RNAIdentifier{};
    EntityIdInt ProteinInBuildingProcessIdentifier{};
public:
    bool RNAInOneParticle = true;
};

inline CellEngineConfigData CellEngineConfigDataObject;

#endif
