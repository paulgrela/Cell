
#include <cmath>
#include <fstream>
#include <cstring>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

Atom::Atom(const char* PDBRecord, UnsignedIntType AtomIndex)
{
    try
    {
        ParseRecord(PDBRecord);
        this->AtomIndex = AtomIndex;
    }
    CATCH("initation of atom")
}

void Atom::ParseRecord(const char* LocalPDBRecord)
{
    try
    {
        string RecordStr = LocalPDBRecord;

        Serial = stoi(RecordStr.substr(6, 5));
        Name = trim_whitespace_surrounding(RecordStr.substr(12, 4));
        ResName = trim_whitespace_surrounding(RecordStr.substr(17, 3));

        X = stod(RecordStr.substr(30, 8));
        Y = stod(RecordStr.substr(38, 8));
        Z = stod(RecordStr.substr(46, 8));
    }
    CATCH("parsing arom record")
}

DoubleVectorType Atom::Position() const
{
	return DoubleVectorType(X, Y, Z);
}

PDBDataFile::PDBDataFile(const char* FileName) : FileName(FileName)
{
    ChosenStructureIndex = 0;

	ReadDataFromFile();
}

void PDBDataFile::ReadDataFromFile()
{
    try
    {
        std::vector<Atom> StructureObject;

        const IntType MaxLineSize = 256;
        char Line[MaxLineSize] = "\0";

        std::ifstream File(FileName, std::ios_base::in);

        while (File.getline(Line, MaxLineSize))
        {
            if (strstr(Line, "ATOM") == Line)
                StructureObject.emplace_back(Atom(Line, StructureObject.size()));
            else
            if (strstr(Line, "END") == Line)
            {
                Atoms.push_back(StructureObject);
                StructureObject.clear();
            }
        }
        File.close();
    }
    CATCH("reading data from PDB file")
}

IntType PDBDataFile::GetNumberOfAtoms() const
{
	return Atoms[ChosenStructureIndex].size();
}

const Atom& PDBDataFile::GetAtom(IntType DataRawIndex) const
{
	return Atoms[ChosenStructureIndex][DataRawIndex];
}

DoubleVectorType PDBDataFile::MassCenter() const
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
