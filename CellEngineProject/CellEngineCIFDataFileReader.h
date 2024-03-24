#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include "CellEngineDataFile.h"

struct TransformationMatrix3x4
{
    float Matrix[3][4];
};

class CellEngineCIFDataFileReader : virtual public CellEngineDataFile
{
private:
    std::unordered_map<UnsignedInt, TransformationMatrix3x4> TransformationsMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
protected:
    void ReadDataFromFile() override;
public:
    explicit CellEngineCIFDataFileReader() = default;
    ~CellEngineCIFDataFileReader() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

class CellEngineParticlesDataFileReader : virtual public CellEngineDataFile
{
protected:
    void ReadDataFromFile() override;
public:
    explicit CellEngineParticlesDataFileReader() = default;
    ~CellEngineParticlesDataFileReader() = default;
public:
    void ReadParticlesFromFile();
    void ReadParticlesFromFileAndPrepareData();
};

#endif