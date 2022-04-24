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
    char Name[2];
    char ResName[4];
    double X;
    double Y;
    double Z;
    char Chain[6];
    UnsignedIntType EntityId;
public:
    Atom(const char* CIFRecord, UnsignedIntType AtomIndex);
    Atom(double XParam, double YParam, double ZParam,  UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6]) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    ~Atom() = default;
public:
    [[nodiscard]] DoubleVectorType Position() const;
private:
    void ParseRecord(const char* LocalPDBRecord);
};

struct Matrix
{
public:
    UnsignedIntType MatrixId;
    double Matrix[3][4];
};

class CIFDataFile
{
private:
    std::vector<Matrix> Matrixes;
    std::unordered_map<std::string, std::vector<Atom>> ChainsNames;
private:
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