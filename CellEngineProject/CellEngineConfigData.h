
#ifndef CELL_ENGINE_CONFIG_DATA_H
#define CELL_ENGINE_CONFIG_DATA_H

#include "vmath.h"

#include "ExceptionsMacro.h"
#include "CellEngineTypes.h"

struct ParticleKind
{
    int Identifier;
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vmath::vec3 Color;
    std::string NameFromXML;
    std::string NameFromDataFile;
};

struct AtomKind
{
    std::string Name;
    float SizeX;
    float SizeY;
    float SizeZ;
    vmath::vec3 Color;
};

inline bool operator==(const AtomKind& AtomKindParameter, const std::string& NameStr)
{
    return AtomKindParameter.Name == NameStr;
}

class CellEngineConfigData
{
public:
    int XTopMainWindow;
    int YTopMainWindow;
    int WidthMainWindow;
    int HeightMainWindow;
    int XTopMenuWindow;
    int YTopMenuWindow;
    int WidthMenuWindow;
    int HeightMenuWindow;
    int XTopSecondWindow;
    int YTopSecondWindow;
    int WidthSecondWindow;
    int HeightSecondWindow;
public:
    bool MultiThreaded;
    bool SetProcessPriorityHighest;
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

    uint64_t MaximalNumberOfLinesInOneFile = 100000;
public:
    std::string CellStateFileName;
public:
    float SpecularPower = 8.0f;
    float SpecularAlbedo = 0.2222f;
public:
    enum class MakeColorsType
    {
        DrawColorForEveryAtom = 1,
        DrawColorForEveryParticle = 2,
        DrawRandomColorForEveryParticle = 3
    };
    MakeColorsType MakeColorsTypeObject;
public:
    enum class StencilForDrawingObjectsTypes
    {
        StencilForDrawingOnlyParticlesCenters = 1,
        StencilForDrawingOnlyInAtomScale = 2
    };
    StencilForDrawingObjectsTypes StencilForDrawingObjectsTypesObject;
    std::uint64_t NumberOfStencilBufferLoops;
    bool UseStencilBuffer = false;
public:
    bool DrawBondsBetweenParticlesCenters;
    bool DrawBondsBetweenAtoms;
public:
    bool ShowDetailsInAtomScale = false;
    bool CheckAtomVisibility;
    float CutZ;
    float Distance;
    std::uint64_t LoadOfAtomsStep;
    bool ShowAtomsInEachPartOfTheCellWhenObserverIsFromOutside = false;
public:
    float XLowToDrawInAtomScale;
    float XHighToDrawInAtomScale;
    float YLowToDrawInAtomScale;
    float YHighToDrawInAtomScale;
    float ZLowToDrawInAtomScale;
    float ZHighToDrawInAtomScale;
public:
    enum class SizeOfAtomsDrawingTypes
    {
        AtomSize = 1,
        ParticleSize = 2,
        AutomaticChangeSize = 3
    };
    SizeOfAtomsDrawingTypes SizeOfAtomsDrawingTypesObject;
public:
    float SizeOfAtomChangeStep;
    float SizeOfAtomX;
    float SizeOfAtomY;
    float SizeOfAtomZ;
    float RotationAngle1;
    float RotationAngle2;
    float RotationAngle3;
    float CameraXPosition;
    float CameraYPosition;
    float CameraZPosition;
    float CameraXMoveShortStep;
    float CameraYMoveShortStep;
    float CameraZMoveShortStep;
    float CameraXMoveLongStep;
    float CameraYMoveLongStep;
    float CameraZMoveLongStep;
    float ViewPositionX;
    float ViewPositionY;
    float ViewPositionZ;
    float ViewXMoveShortStep;
    float ViewYMoveShortStep;
    float ViewZMoveShortStep;
    float ViewXMoveLongStep;
    float ViewYMoveLongStep;
    float ViewZMoveLongStep;
public:
    bool ViewChangeUsingLongStep;
    bool AutomaticChangeOfSizeOfAtom;
    bool AutomaticChangeOfLoadAtomsStep;
public:
    UnsignedIntType ChosenStructureIndex = 0;
public:
    vmath::vec3 BackgroundColors[4];
    std::uint64_t ChosenBackgroundColor;
public:
    bool PrintAtomDescriptionOnScreen;
    bool LogParametersOfRenderingToFile;
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
    std::unordered_map<UnsignedIntType, ParticleKind> ParticlesKindsXML;
    std::unordered_map<UnsignedIntType, UnsignedIntType> ParticlesKindsPos;
public:
    bool ImGuiDemoWindowMenu = false;
    bool ImGuiLightVersion = false;
public:
    int ChosenShapeOfAtoms = 1;
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
};

inline CellEngineConfigData CellEngineConfigDataObject;

#endif
