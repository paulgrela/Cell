
#ifndef CELL_ENGINE_CONFIG_DATA_H
#define CELL_ENGINE_CONFIG_DATA_H

#include "vmath.h"
#include <map>
#include <random>
#include <boost/property_tree/ptree.hpp>

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
    bool MixedFullAtomWithVoxelSpace = false;
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
        FullAtomSimulationSpace = 1,
        VoxelSimulationSpace = 2,
    };
    TypesOfSpace TypeOfSpace = TypesOfSpace::FullAtomSimulationSpace;
public:
    bool MultiThreaded{};
    bool SetProcessPriorityHighest{};
public:
    bool FullAtomMPIParallelProcessesExecution = false;
public:
    bool OpenGLGraphicsSwitchedOff = false;
public:
    UnsignedInt NumberOfParticlesSectorsInX = 1;
    UnsignedInt NumberOfParticlesSectorsInY = 1;
    UnsignedInt NumberOfParticlesSectorsInZ = 1;
public:
    RealType ShiftCenterX = 2560.0f;
    RealType ShiftCenterY = 2560.0f;
    RealType ShiftCenterZ = 2560.0f;
public:
    RealType SizeOfParticlesSectorX = 128.0f;
    RealType SizeOfParticlesSectorY = 128.0f;
    RealType SizeOfParticlesSectorZ = 128.0f;
public:
    UnsignedInt SizeOfSimulationSpaceInEachDimension{};
public:
    float DivisionFactorForSimulationSpace{};
    float DivisionFactorForReadingPositionsOfParticles{ 1.0 };
    float DivisionFactorForGeneratingPositionsOfParticles{ 1.0 };
public:
    UnsignedInt Radius1ForGenerationOfParticles{ 400 };
    UnsignedInt Radius1SizeForGenerationOfParticles{ 400 };
    UnsignedInt Radius2ForGenerationOfParticles{ 420 };
    UnsignedInt Radius2SizeForGenerationOfParticles{ 45 };
public:
    bool CompareBoundsBySectorsBounds = true;
    bool CompareBoundsBySpaceBounds = true;
public:
    bool DNAPaired = true;
    UnsignedInt SizeOfAdvanceDuringGenerationOfRandomDNA = 50;
public:
    UnsignedInt SimulationSpaceSelectionStartXPos{}, SimulationSpaceSelectionStartYPos{}, SimulationSpaceSelectionStartZPos{};
    UnsignedInt SimulationSpaceSelectionStepX{}, SimulationSpaceSelectionStepY{}, SimulationSpaceSelectionStepZ{};
    UnsignedInt SimulationSpaceSelectionSizeX{}, SimulationSpaceSelectionSizeY{}, SimulationSpaceSelectionSizeZ{};
public:
    UnsignedInt SizeOfBigPartOfTheCellMultiplyFactor = 1;
public:
    int NumberOfStepsInSimulationOutside = 1;
    int NumberOfStepsInSimulationInside = 1;
    enum class TypesOfSimulation : UnsignedInt
    {
        BothReactionsAndDiffusion = 1,
        OnlyReactions = 2,
        OnlyDiffusion = 3
    };
    TypesOfSimulation TypeOfSimulation = TypesOfSimulation::BothReactionsAndDiffusion;
public:
    UnsignedInt NumberOfXThreadsInSimulation = 1;
    UnsignedInt NumberOfYThreadsInSimulation = 1;
    UnsignedInt NumberOfZThreadsInSimulation = 1;
    UnsignedInt NumberOfXSectorsInOneThreadInSimulation = 20;
    UnsignedInt NumberOfYSectorsInOneThreadInSimulation = 20;
    UnsignedInt NumberOfZSectorsInOneThreadInSimulation = 20;
    UnsignedInt SizeOfXInOneThreadInSimulationSpace = 1024;
    UnsignedInt SizeOfYInOneThreadInSimulationSpace = 1024;
    UnsignedInt SizeOfZInOneThreadInSimulationSpace = 1024;
    UnsignedInt SizeOfXInOneSectorInOneThreadInSimulationSpace = 32;
    UnsignedInt SizeOfYInOneSectorInOneThreadInSimulationSpace = 32;
    UnsignedInt SizeOfZInOneSectorInOneThreadInSimulationSpace = 32;
public:
    UnsignedInt StepToChangeSpaceDivisionForThreads = 2;
public:
    enum class TypesOfExchangeOfParticlesBetweenThreads : UnsignedInt
    {
        InMainThread = 0,
        ParallelInsert = 1,
        ParallelExtract = 2,
        NoExchangePhase = 3
    };
    TypesOfExchangeOfParticlesBetweenThreads TypeOfExchangeOfParticlesBetweenThreads = TypesOfExchangeOfParticlesBetweenThreads::ParallelInsert;
public:
    bool CheckOnlyParticlesCenters = true;
public:
    bool UseMutexBetweenMainScreenThreadAndMenuThreads = true;
public:
    double RadiusOfCellForDNA{};
public:
    enum class SelectedSpaceStartParametersDrawTypes : UnsignedInt
    {
        DrawFromCenter = 1,
        DrawFromCorner = 2
    };
    SelectedSpaceStartParametersDrawTypes SelectedSpaceStartParametersDrawTypesObject = SelectedSpaceStartParametersDrawTypes::DrawFromCenter;
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
    std::string CellGenomePositionsFileName;
    std::string CellGenomeSequenceFileName;
    std::string CellGenomeSequenceFastaFileName;
    std::string CellGenomeSequenceFastaOneFileName;
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
    UnsignedInt NumberOfStencilBufferLoops{};
    bool UseStencilBuffer = false;
public:
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
    EntityIdInt ATP_ID{};
    EntityIdInt CTP_ID{};
    EntityIdInt GTP_ID{};
    EntityIdInt TTP_ID{};
    EntityIdInt UTP_ID{};
public:
    bool RNAInOneParticle = true;
public:
    bool ReverseReactantsAndProductsBecauseOfFormerErrorBool = true;
};

inline CellEngineConfigData CellEngineConfigDataObject;

#endif
