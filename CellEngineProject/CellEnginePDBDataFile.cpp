
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

CellEnginePDBDataFile::CellEnginePDBDataFile(const string_view FileName)
{
    NumberOfStencilBufferLoop = 3;
    StencilForParticlesCenters = true;

    MakeColorsTypeObject = MakeColorsType::DrawColorForEveryAtom;

    DrawBondsBetweenParticlesCenters = false;
    DrawBondsBetweenAtoms = false;

    CheckAtomVisibility = false;
    CutZ = 1000;
    Distance = 1000;

    SizeX = 1;
    SizeY = 1;
    SizeZ = 1;

    SizeStep = 0.01;

    ViewStep = 3;

    CameraXMoveStep = 1;
    CameraYMoveStep = 1;
    CameraZMoveStep = 1;

    ViewXMoveStep = 1;
    ViewYMoveStep = 1;
    ViewZMoveStep = 1;

    ShowDetailsInAtomScale = false;

    ChosenStructureIndex = 0;

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