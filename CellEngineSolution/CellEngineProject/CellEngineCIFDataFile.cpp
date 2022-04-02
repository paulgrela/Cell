
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEngineCIFDataFile.h"

using namespace std;
using namespace string_utils;

Atom::Atom(const char* CIFRecord, UnsignedIntType AtomIndex)
{
    try
    {
        ParseRecord(CIFRecord);
        this->AtomIndex = AtomIndex;
    }
    CATCH("initation of Atom")
}

void Atom::ParseRecord(const char* LocalCIFRecord)
{
    try
    {
        string RecordStr = LocalCIFRecord;

        Serial = stoi(RecordStr.substr(6, 5));
        Name = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        ResName = trim_whitespace_surrounding(RecordStr.substr(17, 3));

        X = stod(RecordStr.substr(30, 8));
        Y = stod(RecordStr.substr(38, 8));
        Z = stod(RecordStr.substr(46, 8));
    }
    CATCH("parsing Atom record")
}

DoubleVectorType Atom::Position() const
{
    return DoubleVectorType(X, Y, Z);
}

CIFDataFile::CIFDataFile(const string_view FileName)
{
    ChosenStructureIndex = 0;

    ReadDataFromFile(FileName);
}

void CIFDataFile::ReadDataFromFile(const std::string_view FileName)
{
    try
    {
        string Line;
        std::vector<Atom> StructureObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(FileName).c_str(), std::ios_base::in);

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == "ATOM" || Line.substr(0, 6) == "HETATM")
                StructureObject.emplace_back(Atom(Line.c_str(), StructureObject.size()));
            else
            if (Line.substr(0, 3) == "END" )
            {
                Atoms.push_back(StructureObject);
                StructureObject.clear();
            }
        }

        LoggersManagerObject.Log(STREAM("FINISHED READING OF CIF FILE"));

        const auto stop_time = chrono::high_resolution_clock::now();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of reading data from CIF file has taken time: ","executing printing duration_time")));

        File.close();
    }
    CATCH("reading data from CIF file")
}

IntType CIFDataFile::GetNumberOfAtoms() const
{
    return Atoms[ChosenStructureIndex].size();
}

const Atom& CIFDataFile::GetAtom(IntType DataRawIndex) const
{
    return Atoms[ChosenStructureIndex][DataRawIndex];
}

DoubleVectorType CIFDataFile::MassCenter() const
{
    DoubleVectorType MassCenter(0.0, 0.0, 0.0);

    try
    {
        for (const Atom& AtomObject : Atoms[ChosenStructureIndex])
            MassCenter += AtomObject.Position();

        MassCenter /= Atoms[ChosenStructureIndex].size();
    }
    CATCH_AND_THROW("counting mass center")

    return MassCenter;
}
