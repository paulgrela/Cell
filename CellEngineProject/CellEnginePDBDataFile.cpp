
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

CellEngineAtom CellEnginePDBDataFile::ParseRecord(const char* LocalPDBRecord)
{
    CellEngineAtom CellEngineAtomObject;

    try
    {
        string RecordStr = LocalPDBRecord;

        CellEngineAtomObject.Serial = stoi(RecordStr.substr(6, 5));
        string NameStr = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        string ResNameStr = trim_whitespace_surrounding(RecordStr.substr(17, 3));
        strncpy(CellEngineAtomObject.Name, NameStr.c_str(), NameStr.length() + 1);
        strncpy(CellEngineAtomObject.ResName, ResNameStr.c_str(), ResNameStr.length() + 1);

        CellEngineAtomObject.X = stof(RecordStr.substr(30, 8));
        CellEngineAtomObject.Y = stof(RecordStr.substr(38, 8));
        CellEngineAtomObject.Z = stof(RecordStr.substr(46, 8));
    }
    CATCH("parsing element record")

    return CellEngineAtomObject;
}

CellEnginePDBDataFile::CellEnginePDBDataFile(const string_view FileName)
{
    ChosenStructureIndex = 0;

	ReadDataFromFile(FileName);
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
                Atoms.push_back(StructureObject);
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

IntType CellEnginePDBDataFile::GetNumberOfElements() const
{
	return Atoms[ChosenStructureIndex].size();
}

const CellEngineAtom& CellEnginePDBDataFile::GetElement(IntType DataRawIndex) const
{
	return Atoms[ChosenStructureIndex][DataRawIndex];
}

FloatVectorType CellEnginePDBDataFile::MassCenter()
{
	FloatVectorType MassCenter(0.0, 0.0, 0.0);

	try
    {
        for (const CellEngineAtom& AtomObject : Atoms[ChosenStructureIndex])
            MassCenter += AtomObject.Position();

        MassCenter /= Atoms[ChosenStructureIndex].size();
	}
    CATCH_AND_THROW("counting mass center for PDB file")

	return MassCenter;
}
