#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include "CellEngineDataFile.h"
#include "CellEngineBuildParticlesDataOperations.h"
#include "CellEngineParticlesBinaryDataFileReaderWriter.h"

struct TransformationMatrix3x4
{
    float Matrix[3][4];
};

class CellEngineCIFDataFileReader : virtual public CellEngineDataFile, virtual public CellEngineBuildParticlesDataOperations
{
private:
    std::unordered_map<UnsignedInt, TransformationMatrix3x4> TransformationsMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
public:
    explicit CellEngineCIFDataFileReader() = default;
    ~CellEngineCIFDataFileReader() = default;
public:
    void ReadDataFromCIFFile();
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif