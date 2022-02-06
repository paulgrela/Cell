#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>

#include <memory>

#include "VectorType.h"

#include "CellEngineTypes.h"

struct Atom
{
public:
	static bool IsDataRecord(const char* PDBRecord);
	static char* GetAtomSymbol(const char* Name, char* elementSymbol);

	char* PDBRecord;	
	char RecordName[7];
	UnsignedIntType AtomIndex;	
	IntType Serial;	
	char Name[5];	
	char AltLoc;
	char ResName[4];	
	char ChainID;
	IntType ResSeq;
	char ICode;

	double X;
	double Y;
	double Z;
	
	double Occupancy;
	double TempFactor;
	
	char Element[3];
	double Charge;

	Atom(const char* PDBRecord, UnsignedIntType AtomIndex);
	~Atom();
	[[nodiscard]] DoubleVectorType Position() const;

private:
	void ParseRecord(const char* LocalPDBRecord);
};

class PDBDataFile
{
private:
	static bool IsMarkerOfLastFrame(const char* PDBRecord);

	const char* FileName;
	IntType NumberOfDataRaws;

	std::vector<std::unique_ptr<Atom>> Atoms;

	void ReadDataFromFile();
public:
	explicit PDBDataFile(const char* FileName);
	~PDBDataFile() = default;
public:
	[[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const std::vector<std::unique_ptr<Atom>>& GetAtoms() const
	{
		return Atoms;
	}
    [[nodiscard]] const std::unique_ptr<Atom>& GetAtom(IntType Index) const;
    [[nodiscard]] DoubleVectorType MassCenter() const;
};

#endif
