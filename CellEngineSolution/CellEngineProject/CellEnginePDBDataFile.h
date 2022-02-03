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

//odpowiada linii z pliku PDB
struct Atom
{
public:
	static bool IsDataRecord(const char* PDBRecord);
	static char* GetAtomSymbol(const char* Name, char* elementSymbol);

	//pola, nazwy pol zgodne z nazwami z dokumentacji pliku PDB
	char* PDBRecord;
	
	//std::string RecordName; //dodatkowy znak na \0
	char RecordName[7]; //dodatkowy znak na \0

	UnsignedIntType AtomIndex;
	
	IntType Serial;
	
	//std::string Name;
	char Name[5];
	
	char AltLoc;
	
	//std::string ResName;
	char ResName[4];
	
	char ChainID;
	IntType ResSeq;
	char ICode;
	double X;
	double Y;
	double Z;
	double Occupancy;
	double TempFactor;
	
	//std::string Element;
	char Element[3];

	double Charge;

	Atom(const char* PDBRecord, UnsignedIntType AtomIndex);
	~Atom();
	DoubleVectorType Position() const;

private:
	void ParseRecord(const char* PDBRecord);

};

class PDBDataFile
{
private:
	static bool IsMarkerOfLastFrame(const char* PDBRecord);

	const char* FileName;
	IntType NumberOfDataRaws; //ilosc atomow

	//Atom** Atoms; //Atoms z jednej klatki
	//std::vector<Atom> Atoms;

	std::vector<std::unique_ptr<Atom>> Atoms;

	void AnalyzeFile();
	void ReadDataFromFile();

public:
	PDBDataFile(const char* FileName);
	~PDBDataFile();
public:
	IntType GetNumberOfAtoms() const;
	const std::vector<std::unique_ptr<Atom>>& GetAtoms() const
	{
		return Atoms;
	}
	const std::unique_ptr<Atom>& GetAtom(IntType Index) const;
	DoubleVectorType MassCenter() const;
};

#endif
