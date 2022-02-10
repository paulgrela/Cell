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
	char RecordName[7];
	UnsignedIntType AtomIndex;	
	IntType Serial;	
	char Name[5];	
	char ResName[4];
	double X;
	double Y;
	double Z;
	char Element[3];
public:
	Atom(const char* PDBRecord, UnsignedIntType AtomIndex);
	~Atom() = default;
public:
	[[nodiscard]] DoubleVectorType Position() const;
public:
    static char* GetAtomSymbol(const char* Name, char* elementSymbol);
private:
	void ParseRecord(const char* LocalPDBRecord);
};

class PDBDataFile
{
private:
	const char* FileName;
    std::vector<std::vector<Atom>> Atoms;
private:
	void ReadDataFromFile();
public:
    UnsignedIntType ChosenStructureIndex;
public:
	explicit PDBDataFile(const char* FileName);
	~PDBDataFile() = default;
public:
    [[nodiscard]] const std::vector<Atom>& GetAtoms() const
    {
        return Atoms[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures() const
    {
        return Atoms.size();
    }
    [[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const Atom& GetAtom(IntType Index) const;
    [[nodiscard]] DoubleVectorType MassCenter() const;
};

#endif
