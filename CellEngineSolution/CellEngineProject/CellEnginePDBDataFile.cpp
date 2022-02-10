
#include <cmath>
#include <fstream>
#include <cstring>

#include "ExceptionsMacro.h"

#include "StringUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

char* Atom::GetAtomSymbol(const char* Name, char* AtomSymbol)
{
    try
    {
        strcpy(AtomSymbol, "\0\0\0");

        IntType CharIndex = 0;
        do
        {
            AtomSymbol[0] = Name[CharIndex];
            CharIndex++;
        }
        while (AtomSymbol[0] < 'A' || AtomSymbol[0] > 'Z');

        if (Name[CharIndex] >= 'a' && Name[CharIndex] <= 'z')
            AtomSymbol[1] = Name[CharIndex];
        else
            AtomSymbol[1] = '\0';

        AtomSymbol[2] = '\0';
    }
    CATCH("getting atom symbol")

	return AtomSymbol;
}

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
        char SSerial[5];
        char SX[8], SY[8], SZ[8];
        char* EndPtr;

        strncpy(RecordName, LocalPDBRecord, 6);
        RecordName[6] = '\0';

        strncpy(SSerial, LocalPDBRecord + 6, 5);
        Serial = strtol(SSerial, &EndPtr, 10);

        strncpy(Name, LocalPDBRecord + 12, 4);
        Name[4] = '\0';

        strncpy(ResName, LocalPDBRecord + 17, 3);
        ResName[3] = '\0';

        strncpy(SX, LocalPDBRecord + 30, 8);
        X = strtod(SX, &EndPtr);
        strncpy(SY, LocalPDBRecord + 38, 8);
        Y = strtod(SY, &EndPtr);
        strncpy(SZ, LocalPDBRecord + 46, 8);
        Z = strtod(SZ, &EndPtr);

        strncpy(Element, LocalPDBRecord + 76, 2);
        Element[2] = '\0';

        string RecordNameStr = trim_whitespace_surrounding(RecordName);
        strncpy(RecordName, RecordNameStr.c_str(), RecordNameStr.length());

        string NameStr = trim_whitespace_surrounding(Name);
        strncpy(Name, NameStr.c_str(), NameStr.length());

        string ResNameStr = trim_whitespace_surrounding(ResName);
        strncpy(ResName, ResNameStr.c_str(), ResNameStr.length());

        string ElementStr = trim_whitespace_surrounding(Element);
        strncpy(Element, ElementStr.c_str(), ElementStr.length());
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
