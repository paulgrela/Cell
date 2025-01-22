#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include "CellEngineDataFile.h"

class CellEnginePDBDataFileReader: public CellEngineDataFile
{
private:
    void ReadDataFromFile(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type) override;
public:
    explicit CellEnginePDBDataFileReader() = default;
    ~CellEnginePDBDataFileReader() override = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif
