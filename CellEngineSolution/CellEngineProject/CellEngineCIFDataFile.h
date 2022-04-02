#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>

#include <memory>
#include <unordered_map>

#include "VectorType.h"

#include "CellEngineTypes.h"

struct Atom
{
public:
    UnsignedIntType AtomIndex;
    IntType Serial;
    std::string Name;
    std::string ResName;
    double X;
    double Y;
    double Z;
    std::string Chain;
    UnsignedIntType EntityId;
public:
    Atom(const char* CIFRecord, UnsignedIntType AtomIndex);
    ~Atom() = default;
public:
    [[nodiscard]] DoubleVectorType Position() const;
private:
    void ParseRecord(const char* LocalPDBRecord);
};

struct Entity
{
public:
    UnsignedIntType  EntityId;
    std::vector<Atom> Atoms;
};

struct Matrix
{
public:
    double Matrix[3][3];
};

class CIFDataFile
{
private:
    std::vector<Matrix> Matrixes;
    std::unordered_map<std::string, Entity> Entities; // Created during reading atoms for temporary purpose
private:
    std::vector<Atom> AtomsFinal;
    std::vector<std::vector<Atom>> Atoms;
private:
    void ReadDataFromFile(const std::string_view);
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit CIFDataFile(const std::string_view FileName);
    ~CIFDataFile() = default;
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