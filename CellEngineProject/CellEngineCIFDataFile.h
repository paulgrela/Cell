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

#include "CellEngineDataFile.h"

struct AtomsPositionMatrix3x4
{
public:
    float Matrix[3][4];
};

class CellEngineCIFDataFile: public CellEngineDataFile
{
private:
    std::vector<AtomsPositionMatrix3x4> AtomsPositionsMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
private:
    void ReadDataFromFile(std::string_view LocalCIFRecord);
public:
    explicit CellEngineCIFDataFile(std::string_view FileName);
    ~CellEngineCIFDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
    [[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const CellEngineAtom& GetAtom(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() override;
};

#endif