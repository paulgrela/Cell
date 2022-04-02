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

struct Element
{
public:
	UnsignedIntType ElementIndex;
	IntType Serial;
	std::string Name;
	std::string ResName;
	double X;
	double Y;
	double Z;
public:
	Element(const char* PDBRecord, UnsignedIntType ElementIndex);
	~Element() = default;
public:
	[[nodiscard]] DoubleVectorType Position() const;
private:
	void ParseRecord(const char* LocalPDBRecord);
};

class PDBDataFile
{
private:
    std::vector<std::vector<Element>> Elements;
private:
	void ReadDataFromFile(const std::string_view FileName);
public:
    UnsignedIntType ChosenStructureIndex;
public:
	explicit PDBDataFile(const std::string_view FileName);
	~PDBDataFile() = default;
public:
    [[nodiscard]] const std::vector<Element>& GetElements() const
    {
        return Elements[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures() const
    {
        return Elements.size();
    }
    [[nodiscard]] IntType GetNumberOfElements() const;
    [[nodiscard]] const Element& GetElement(IntType Index) const;
    [[nodiscard]] DoubleVectorType MassCenter() const;
};

#endif
