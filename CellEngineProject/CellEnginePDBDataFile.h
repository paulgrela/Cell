#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include <vector>
#include <cstring>

#include <fstream>
#include <stdexcept>

#include <memory>

#include "VectorType.h"

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"

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
