
#include <cmath>
#include <fstream>
#include <cstring>

#include "StringUtils.h"
#include "CellEnginePDBDataFile.h"

using namespace std;
using namespace string_utils;

bool Atom::IsDataRecord(const char* PDBRecord)
{
	return strstr(PDBRecord, "ATOM") == PDBRecord;
}

char* Atom::GetAtomSymbol(const char* Name, char* AtomSymbol)
{
	strcpy(AtomSymbol, "\0\0\0");

	IntType Index = 0;
	do
	{
		AtomSymbol[0] = Name[Index];
		Index++;
	} 
	while (AtomSymbol[0] < 'A' || AtomSymbol[0] > 'Z');

	if (Name[Index] >= 'a' && Name[Index] <= 'z')
	{
		AtomSymbol[1] = Name[Index];
	}
	else 
		AtomSymbol[1] = '\0';

	AtomSymbol[2] = '\0';

	return AtomSymbol;
}

Atom::Atom(const char* PDBRecord, UnsignedIntType AtomIndex)
{
	this->PDBRecord = new char[strlen(PDBRecord) + 1];
	strcpy(this->PDBRecord, PDBRecord);
	ParseRecord(this->PDBRecord);
	this->AtomIndex = AtomIndex;
}

Atom::~Atom()
{
	delete[] PDBRecord;
}

void Atom::ParseRecord(const char* LocalPDBRecord)
{
	char SSerial[5], SResSeq[4];
	char SX[8], SY[8], SZ[8];
	char SOccupancy[6], STempFactor[6], SCharge[3];
	char* EndPtr;

	strncpy(RecordName, LocalPDBRecord, 6);
	RecordName[6] = '\0';
	
	strncpy(SSerial, LocalPDBRecord + 6, 5);
	Serial = strtol(SSerial, &EndPtr, 10);
	
	strncpy(Name, LocalPDBRecord + 12, 4);
	Name[4] = '\0';

	AltLoc = LocalPDBRecord[16];
	
	strncpy(ResName, LocalPDBRecord + 17, 3);
	ResName[3] = '\0';
	
	ChainID = LocalPDBRecord[21];
	strncpy(SResSeq, LocalPDBRecord + 22, 4);
	ResSeq = strtol(SResSeq, &EndPtr, 10);
	ICode = LocalPDBRecord[26];
	strncpy(SX, LocalPDBRecord + 30, 8);
	X = strtod(SX, &EndPtr);
	strncpy(SY, LocalPDBRecord + 38, 8);
	Y = strtod(SY, &EndPtr);
	strncpy(SZ, LocalPDBRecord + 46, 8);
	Z = strtod(SZ, &EndPtr);
	strncpy(SOccupancy, LocalPDBRecord + 54, 6);
	Occupancy = strtod(SOccupancy, &EndPtr);
	strncpy(STempFactor, LocalPDBRecord + 60, 6);
	TempFactor = strtod(STempFactor, &EndPtr);
	
	strncpy(Element, LocalPDBRecord + 76, 2);
	Element[2] = '\0';

	strncpy(SCharge, LocalPDBRecord + 78, 2);
	Charge = strtod(SCharge, &EndPtr);

	string RecordNameStr = trim_whitespace_surrounding(RecordName);
	strncpy(RecordName, RecordNameStr.c_str(), RecordNameStr.length());
	
	string NameStr = trim_whitespace_surrounding(Name);
	strncpy(Name, NameStr.c_str(), NameStr.length());
	
	string ResNameStr = trim_whitespace_surrounding(ResName);
	strncpy(ResName, ResNameStr.c_str(), ResNameStr.length());

	string ElementStr = trim_whitespace_surrounding(Element);
	strncpy(Element, ElementStr.c_str(), ElementStr.length());
}

DoubleVectorType Atom::Position() const
{
	return DoubleVectorType(X, Y, Z);
}

PDBDataFile::PDBDataFile(const char* FileName) : FileName(FileName), NumberOfDataRaws(0)
{
	AnalyzeFile();
	ReadDataFromFile();
}

bool PDBDataFile::IsMarkerOfLastFrame(const char* PDBRecord)
{
	return strstr(PDBRecord, "END") == PDBRecord;
}

void PDBDataFile::AnalyzeFile()
{
	const IntType MaxLineSize = 256;
	char Line[MaxLineSize] = "\0";

	NumberOfDataRaws = 0;

	std::ifstream File(FileName, std::ios_base::in);

	while (!IsMarkerOfLastFrame(Line))
	{
		File.getline(Line, MaxLineSize);
		if (Atom::IsDataRecord(Line))
			NumberOfDataRaws++;
	}

	File.close();
}

void PDBDataFile::ReadDataFromFile()
{
	const IntType MaxLineSize = 256;
	char Line[MaxLineSize] = "\0";

	std::ifstream File(FileName, std::ios_base::in);

	for (IntType DataRawIndex = 0; DataRawIndex < NumberOfDataRaws; DataRawIndex++)
	{
		File.getline(Line, MaxLineSize);

		if (!Atom::IsDataRecord(Line))
			DataRawIndex--;
		else
		{
			Atoms.push_back(make_unique<Atom>(Line, DataRawIndex));

			if (IsMarkerOfLastFrame(Line))
			{
				throw std::logic_error("Too early end of frame");
				File.close();
			}
		}
	}

	File.close();
}

IntType PDBDataFile::GetNumberOfAtoms() const
{
	return NumberOfDataRaws;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-local-addr"
const unique_ptr<Atom>& PDBDataFile::GetAtom(IntType DataRawIndex) const
{
	if (DataRawIndex < 0 || DataRawIndex >= NumberOfDataRaws)
		return std::unique_ptr<Atom>(nullptr);
	return Atoms[DataRawIndex];
}
#pragma GCC diagnostic pop

DoubleVectorType PDBDataFile::MassCenter() const
{
	DoubleVectorType MassCenter(0.0, 0.0, 0.0);
	for (IntType DataRawIndex = 0; DataRawIndex < NumberOfDataRaws; DataRawIndex++)
	{
		MassCenter += Atoms[DataRawIndex]->Position();
	}
	MassCenter /= (double)NumberOfDataRaws;
	return MassCenter;
}
