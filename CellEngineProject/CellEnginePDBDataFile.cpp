
#include <fstream>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

Element::Element(const char* PDBRecord, UnsignedIntType ElementIndex)
{
    try
    {
        ParseRecord(PDBRecord);
        this->ElementIndex = ElementIndex;
    }
    CATCH("initation of element")
}

void Element::ParseRecord(const char* LocalPDBRecord)
{
    try
    {
        string RecordStr = LocalPDBRecord;

        Serial = stoi(RecordStr.substr(6, 5));
        Name = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        ResName = trim_whitespace_surrounding(RecordStr.substr(17, 3));

        X = stof(RecordStr.substr(30, 8));
        Y = stof(RecordStr.substr(38, 8));
        Z = stof(RecordStr.substr(46, 8));
    }
    CATCH("parsing element record")
}

FloatVectorType Element::Position() const
{
	return FloatVectorType(X, Y, Z);
}

PDBDataFile::PDBDataFile(const string_view FileName)
{
    ChosenStructureIndex = 0;

	ReadDataFromFile(FileName);
}

void PDBDataFile::ReadDataFromFile(const string_view FileName)
{
    try
    {
        string Line;
        std::vector<Element> StructureObject;

        const auto start_time = chrono::high_resolution_clock::now();

        std::ifstream File(string(FileName).c_str(), std::ios_base::in);

        LoggersManagerObject.Log(STREAM("STARTED READING OF PDB FILE"));

        while (getline(File, Line, '\n'))
        {
            if (Line.substr(0, 4) == "ATOM" || Line.substr(0, 6) == "HETATM")
                StructureObject.emplace_back(Element(Line.c_str(), StructureObject.size()));
            else
            if (Line.substr(0, 3) == "END" )
            {
                StructureObject.emplace_back(Element());
                Elements.push_back(StructureObject);
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

IntType PDBDataFile::GetNumberOfElements() const
{
	return Elements[ChosenStructureIndex].size();
}

const Element& PDBDataFile::GetElement(IntType DataRawIndex) const
{
	return Elements[ChosenStructureIndex][DataRawIndex];
}

FloatVectorType PDBDataFile::MassCenter() const
{
	FloatVectorType MassCenter(0.0, 0.0, 0.0);

	try
    {
        for (const Element& ElementObject : Elements[ChosenStructureIndex])
            MassCenter += ElementObject.Position();

        MassCenter /= Elements[ChosenStructureIndex].size();
	}
    CATCH_AND_THROW("counting mass center for PDB file")

	return MassCenter;
}
