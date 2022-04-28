#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include <vector>
#include <string.h>

#include <fstream>
#include <stdexcept>

#include <memory>

#include "VectorType.h"

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"

/*
struct Element
{
public:
	UnsignedIntType ElementIndex;
	IntType Serial;
	std::string Name;
	std::string ResName;
	float X;
	float Y;
	float Z;
public:
	Element(const char* PDBRecord, UnsignedIntType ElementIndex);
    Element() = default;
	~Element() = default;
public:
	[[nodiscard]] FloatVectorType Position() const;
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
    [[nodiscard]] FloatVectorType MassCenter() const;
};
*/
#include <cstdio>
#include <cstring>
/*
struct Atom
{
public:
    UnsignedIntType AtomIndex;
    IntType Serial;
    char Name[2];
    char ResName[4];
    float X;
    float Y;
    float Z;
    char Chain[6];
    UnsignedIntType EntityId;
public:
    Atom() = default;
    Atom(const char* CIFRecord, UnsignedIntType AtomIndex);
    Atom(float XParam, float YParam, float ZParam,  UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6]) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    ~Atom() = default;
public:
    [[nodiscard]] FloatVectorType Position() const;
private:
    void ParseRecord(const char* LocalPDBRecord);
};
*/
/*
class AtomPDB: public AtomBase
{
public:
    AtomPDB() = default;
    AtomPDB(const char* PDBRecord, UnsignedIntType AtomIndex);
public:
    [[nodiscard]] FloatVectorType Position() const;
    void ParseRecord(const char* LocalPDBRecord);
};
*/

/*
class AtomPDB: public AtomBase
{
public:
    AtomPDB() = default;
    AtomPDB(const char* PDBRecord, UnsignedIntType AtomIndex);
public:
    [[nodiscard]] FloatVectorType Position() const;
    void ParseRecord(const char* LocalPDBRecord);
};
*/

class PDBDataFile
{
private:
    std::vector<std::vector<AtomBase>> Atoms;
private:
    void ReadDataFromFile(const std::string_view FileName);
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit PDBDataFile(const std::string_view FileName);
    ~PDBDataFile() = default;
public:
    static AtomBase ParseRecord(const char* LocalPDBRecord);
    [[nodiscard]] std::vector<AtomBase>& GetAtoms()
    {
        return Atoms[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures() const
    {
        return Atoms.size();
    }
    [[nodiscard]] IntType GetNumberOfElements() const;
    [[nodiscard]] const AtomBase& GetElement(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() const;
};

#endif
