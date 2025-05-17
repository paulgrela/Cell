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
    MainMapType<UnsignedInt, TransformationMatrix3x4> TransformationsMatrixes;
    MainMapType<std::string, std::vector<CellEngineAtom>> ChainsNames;
public:
    explicit CellEngineCIFDataFileReader() = default;
    ~CellEngineCIFDataFileReader() override = default;
public:
    void ReadDataFromCIFFile(bool SetStartValuesBool);
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif