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

#include "CellEngineDataFile.h"

class CellEnginePDBDataFile: public CellEngineDataFile
{
private:
    std::vector<std::vector<CellEngineAtom>> Atoms;
private:
    void ReadDataFromFile(const std::string_view FileName);
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit CellEnginePDBDataFile(const std::string_view FileName);
    ~CellEnginePDBDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
    [[nodiscard]] std::vector<CellEngineAtom>& GetAtoms() override
    {
        return Atoms[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures() const
    {
        return Atoms.size();
    }
    [[nodiscard]] IntType GetNumberOfElements() const;
    [[nodiscard]] const CellEngineAtom& GetElement(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() override;
};

#endif
