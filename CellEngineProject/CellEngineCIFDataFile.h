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
#include "CellEngineAtom.h"

struct Matrix
{
public:
    UnsignedIntType MatrixId;
    float Matrix[3][4];
};

class CIFDataFile
{
private:
    std::vector<Matrix> Matrixes;
    std::unordered_map<std::string, std::vector<AtomBase>> ChainsNames;
private:
    std::vector<std::vector<AtomBase>> Atoms;
private:
    void ReadDataFromFile(const std::string_view LocalCIFRecord);
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit CIFDataFile(const std::string_view FileName);
    ~CIFDataFile() = default;
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
    [[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const AtomBase& GetAtom(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() const;
};

#endif