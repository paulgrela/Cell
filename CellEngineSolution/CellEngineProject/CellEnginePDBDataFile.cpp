
#include <math.h>
#include <fstream>
#include <string.h>

#include "CellEnginePDBDataFile.h"

using namespace std;

// ATOM

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
		//atomSymbol[0]=toupper(Name[DataRawIndex]);
		AtomSymbol[0] = Name[Index];
		Index++;
	} 
	while (AtomSymbol[0] < 'A' || AtomSymbol[0] > 'Z');

	if (Name[Index] >= 'a' && Name[Index] <= 'z')
	{
		AtomSymbol[1] = Name[Index];
	}
	else AtomSymbol[1] = '\0';
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

									string trim_whitespace_surrounding(const string& s)
									{
										const char whitespace[]{ " \t\n" };
										const size_t first(s.find_first_not_of(whitespace));
										if (string::npos == first)
											return {};
										const size_t last(s.find_last_not_of(whitespace));
										return s.substr(first, (last - first + 1));
									}

void Atom::ParseRecord(const char* PDBRecord)
{
	char SSerial[5], SResSeq[4];
	char SX[8], SY[8], SZ[8];
	char SOccupancy[6], STempFactor[6], SCharge[3];
	char* EndPtr;

	strncpy(RecordName, PDBRecord, 6); //0-5
	RecordName[6] = '\0';
	
	strncpy(SSerial, PDBRecord + 6, 5); //6-10, 11=spacja
	Serial = strtol(SSerial, &EndPtr, 10);
	
	strncpy(Name, PDBRecord + 12, 4); //12-15	
	Name[4] = '\0';

	AltLoc = PDBRecord[16]; //16
	
	strncpy(ResName, PDBRecord + 17, 3); //17-19, 20=spacja
	ResName[3] = '\0';
	
	ChainID = PDBRecord[21]; //21
	strncpy(SResSeq, PDBRecord + 22, 4); //22-25
	ResSeq = strtol(SResSeq, &EndPtr, 10);
	ICode = PDBRecord[26]; //26
	strncpy(SX, PDBRecord + 30, 8); //30-37
	X = strtod(SX, &EndPtr);
	strncpy(SY, PDBRecord + 38, 8); //38-45
	Y = strtod(SY, &EndPtr);
	strncpy(SZ, PDBRecord + 46, 8); //46-53
	Z = strtod(SZ, &EndPtr);
	strncpy(SOccupancy, PDBRecord + 54, 6); //54-59
	Occupancy = strtod(SOccupancy, &EndPtr);
	strncpy(STempFactor, PDBRecord + 60, 6); //60-65, 66-75=spacje
	TempFactor = strtod(STempFactor, &EndPtr);
	
	strncpy(Element, PDBRecord + 76, 2); //76-77
	Element[2] = '\0';

	strncpy(SCharge, PDBRecord + 78, 2); //78-79
	Charge = strtod(SCharge, &EndPtr);

	//usuwanie spacji
	/*
	//XXX
	TrimStr(RecordName);
	TrimStr(Name);
	TrimStr(ResName);
	TrimStr(Element);
	*/
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

//PDB

PDBDataFile::PDBDataFile(const char* FileName) : FileName(FileName), NumberOfDataRaws(0)
{
	AnalyzeFile(); //moze zglosic wyjatek
	
	//Atoms.resize(NumberOfDataRaws);	
	//for (unique_ptr<Atom>& AtomObjectPtr : Atoms)
	//	AtomObjectPtr = nullptr;



	//Atoms = new Atom * [NumberOfDataRaws];
	//for (IntType DataRawIndex = 0; DataRawIndex < NumberOfDataRaws; DataRawIndex++)
	//	Atoms[DataRawIndex] = NULL; //potrzebne do usuwania
	ReadDataFromFile(); //moze zglosic wyjatek
}

PDBDataFile::~PDBDataFile()
{
	//for (IntType DataRawIndex = 0; DataRawIndex < NumberOfDataRaws; DataRawIndex++)
	//	delete Atoms[DataRawIndex];
	//delete[] Atoms;
}

bool PDBDataFile::IsMarkerOfLastFrame(const char* PDBRecord)
{
	return strstr(PDBRecord, "END") == PDBRecord;
}

void PDBDataFile::AnalyzeFile() //z analizy wyciagam informacje o wielkosci tablicy atomow, ktora nalezy zaalokowac
{
	//XXX
	//if (!TestFile(FileName))
	//	throw std::logic_error("PDB File does not exist");

	const IntType MaxLineSize = 256; //zgodnie z dokumentacja PDB wystarczyloby 80
	char Line[MaxLineSize] = "\0";

	NumberOfDataRaws = 0;

	std::ifstream File(FileName, std::ios_base::in);

	//dane
	while (!IsMarkerOfLastFrame(Line)) //do konca pierwszej ramki, plik PDB moze miec ich wiele
	{
		File.getline(Line, MaxLineSize);
		if (Atom::IsDataRecord(Line))
			NumberOfDataRaws++;
	}

	File.close();
}


void PDBDataFile::ReadDataFromFile() //jedna klatka
{
	const IntType MaxLineSize = 256; //zgodnie z dokumentacja PDB wystarczyloby 80
	char Line[MaxLineSize] = "\0";

	std::ifstream File(FileName, std::ios_base::in);

	for (IntType DataRawIndex = 0; DataRawIndex < NumberOfDataRaws; DataRawIndex++)
	{
		File.getline(Line, MaxLineSize);

		if (!Atom::IsDataRecord(Line))
			DataRawIndex--;
		else
		{
			//if (Atoms[DataRawIndex] != NULL)
			//	delete Atoms[DataRawIndex]; //kasowanie poprzedniego obiektu

			//Atoms[DataRawIndex] = new Atom(Line);
			//Atoms.push_back(new Atom(Line));
			Atoms.push_back(make_unique<Atom>(Line, DataRawIndex));

			//if (Atoms[DataRawIndex] != nullptr)
			//	Atoms[DataRawIndex].reset();
			//Atoms[DataRawIndex] = make_unique<Atom>(Line);

			/*
			if(Atoms[DataRawIndex]->Serial-1!=DataRawIndex)
			{
				throw std::logic_error("Numer atomu w pliku PDB nie zgadza sie z indeksem");
				File.close();
			}
			*/
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

const unique_ptr<Atom>& PDBDataFile::GetAtom(IntType DataRawIndex) const
{
	if (DataRawIndex < 0 || DataRawIndex >= NumberOfDataRaws)
		return nullptr;
	return Atoms[DataRawIndex];
}

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
