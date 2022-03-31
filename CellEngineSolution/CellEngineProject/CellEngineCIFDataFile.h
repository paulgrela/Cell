#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

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
    UnsignedIntType AtomIndex;
    IntType Serial;
    std::string Name;
    std::string ResName;
    double X;
    double Y;
    double Z;
public:
    Atom(const char* CIFRecord, UnsignedIntType AtomIndex);
    ~Atom() = default;
public:
    [[nodiscard]] DoubleVectorType Position() const;
private:
    void ParseRecord(const char* LocalPDBRecord);
};

class CIFDataFile
{
private:
    const char* FileName;
    std::vector<std::vector<Atom>> Atoms;
private:
    void ReadDataFromFile();
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit CIFDataFile(const char* FileName);
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