
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

        //const IntType MaxLineSize = 256;
        //char Line[MaxLineSize] = "\0";
        string Line;

        //std::ifstream File(FileName, std::ios_base::in);
        std::ifstream File(FileName, std::ios_base::in);
        //int a = 0;
        //while (File.getline(Line, MaxLineSize))
//        bool EndOfModelFound = false;

        while (getline(File, Line, '\n'))
        //while (getline(File, Line, '\n') && a < 100000)
        {
            //a++;
            //if (a % 10000 == 0 )
            //    cout << a << endl;
            //if (strstr(Line, "ATOM") == Line || strstr(Line, "HETATM") == Line)
            if (Line.substr(0, 4) == "ATOM" || Line.substr(0, 6) == "HETATM")
            {
                StructureObject.emplace_back(Atom(Line.c_str(), StructureObject.size()));
//                EndOfModelFound = false;
            }
            else
            if (Line.substr(0, 3) == "END" )//lub koniec pliku
//            //if (strstr(Line, "END") == Line)
            {
                Atoms.push_back(StructureObject);
                StructureObject.clear();
//                EndOfModelFound = true;
            }
        }

//        if (EndOfModelFound == false)
//        {
//            Atoms.push_back(StructureObject);
//            StructureObject.clear();
//        }
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
