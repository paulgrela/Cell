#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>

#include <memory>
#include <unordered_map>

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"

#include "CellEngineDataFile.h"

struct TransformationMatrix3x4
{
    float Matrix[3][4];
};

class CellEngineCIFDataFile: public CellEngineDataFile
{
private:
    std::unordered_map<std::uint64_t, TransformationMatrix3x4> TransformationsMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
private:
    void ReadDataFromFile(std::string_view LocalCIFRecord) override;
public:
    explicit CellEngineCIFDataFile(std::string_view FileName);
    ~CellEngineCIFDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif