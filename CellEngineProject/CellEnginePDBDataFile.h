#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include "CellEngineDataFile.h"

class CellEnginePDBDataFile: public CellEngineDataFile
{
private:
    void ReadDataFromFile(bool StartValuesBool, CellEngineConfigData::TypesOfFileToRead Type) override;
public:
    explicit CellEnginePDBDataFile() = default;
    ~CellEnginePDBDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif
