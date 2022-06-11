
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

CellEnginePDBDataFile::CellEnginePDBDataFile(const string_view FileName)
{
    ChosenStructureIndex = 0;

    NumberOfStencilBufferLoops = 3;
    StencilForDrawingObjectsTypesObject = StencilForDrawingObjectsTypes::StencilForDrawingOnlyParticlesCenters;

    MakeColorsTypeObject = MakeColorsType::DrawColorForEveryAtom;

    DrawBondsBetweenParticlesCenters = false;
    DrawBondsBetweenAtoms = false;

    ShowDetailsInAtomScale = false;
    CheckAtomVisibility = false;
    CutZ = 1000;
    Distance = 1000;
    LoadOfAtomsStep = 0;

    XLowToDrawInAtomScale = 0;
    XHighToDrawInAtomScale = 0;
    YLowToDrawInAtomScale = 0;
    YHighToDrawInAtomScale = 0;
    ZLowToDrawInAtomScale = 0;
    ZHighToDrawInAtomScale = 0;

    SizeOfAtomsDrawingTypesObject = CellEngineDataFile::SizeOfAtomsDrawingTypes::AutomaticChangeSize;

    SizeStep = 0.01;

    SizeOfAtomX = 1;
    SizeOfAtomY = 1;
    SizeOfAtomZ = 1;

    CameraXPosition = 0.0;
    CameraYPosition = 0.0;
    CameraZPosition = 0.0;

    CameraXMoveShortStep = 1;
    CameraYMoveShortStep = 1;
    CameraZMoveShortStep = 1;

    CameraXMoveLongStep = 1;
    CameraYMoveLongStep = 1;
    CameraZMoveLongStep = 1;

    ViewXMoveShortStep = 1;
    ViewYMoveShortStep = 1;
    ViewZMoveShortStep = 1;

    ViewXMoveLongStep = 1;
    ViewYMoveLongStep = 1;
    ViewZMoveLongStep = 1;

    ViewChangeUsingLongStep = false;
    AutomaticChangeOfSizeOfAtom = false;
    AutomaticChangeOfLoadAtomsStep = false;

    BackgroundColors[1] = FromVec4ToVec3(sb7::color::Cyan);
    BackgroundColors[2] = FromVec4ToVec3(sb7::color::White);
    BackgroundColors[3] = FromVec4ToVec3(sb7::color::Black);

    ChosenBackgroundColor = 3;
    //BackgroundColorTypeObject = BackgroundColorType::Background1Color;

	ReadDataFromFile(FileName);
}

CellEngineAtom CellEnginePDBDataFile::ParseRecord(const char* LocalPDBRecord)
{
    CellEngineAtom CellEngineAtomObject;

    try
    {
        string RecordStr = LocalPDBRecord;

        CellEngineAtomObject.EntityId = 0;
        CellEngineAtomObject.Serial = stoi(RecordStr.substr(6, 5));
        string NameStr = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        string ResNameStr = trim_whitespace_surrounding(RecordStr.substr(17, 3));
        strncpy(CellEngineAtomObject.Name, NameStr.c_str(), NameStr.length() + 1);
        strncpy(CellEngineAtomObject.ResName, ResNameStr.c_str(), ResNameStr.length() + 1);

        CellEngineAtomObject.X = stof(RecordStr.substr(30, 8));
        CellEngineAtomObject.Y = stof(RecordStr.substr(38, 8));
        CellEngineAtomObject.Z = stof(RecordStr.substr(46, 8));

        CellEngineAtomObject.AtomColor = CellEngineAtomObject.ParticleColor = CellEngineAtomObject.RandomParticleColor = ChooseColorForAtom(CellEngineAtomObject);
        CellEngineAtomObject.SizeXAtom = CellEngineAtomObject.SizeYAtom = CellEngineAtomObject.SizeZAtom =  CellEngineAtomObject.SizeXParticle = CellEngineAtomObject.SizeYParticle = CellEngineAtomObject.SizeZParticle = 1;
    }
    CATCH("parsing atom record")

    return CellEngineAtomObject;
}

void CellEnginePDBDataFile::ReadDataFromFile(const string_view FileName)
{
    try
    {
        string Line;
        std::vector<CellEngineAtom> StructureObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(FileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING OF PDB FILE"));

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == "ATOM" || Line.substr(0, 6) == "HETATM")
                StructureObject.emplace_back(ParseRecord(Line.c_str()));
            else
            if (Line.substr(0, 3) == "END" )
            {
                StructureObject.emplace_back(CellEngineAtom());
                ParticlesCenters.push_back(StructureObject);
                StructureObject.clear();
            }
        }

        LoggersManagerObject.Log(STREAM("FINISHED READING OF PDB FILE"));

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Reading data from pdb file has taken time: ","executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from PDB file")
}