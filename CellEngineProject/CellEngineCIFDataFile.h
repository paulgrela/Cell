#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

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
    void ReadDataFromFile() override;
public:
    explicit CellEngineCIFDataFile() = default;
    ~CellEngineCIFDataFile() = default;
public:
    CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif