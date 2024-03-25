#pragma once

#ifndef CELL_ENGINE_PDB_DATA_FILE_H
#define CELL_ENGINE_PDB_DATA_FILE_H

#include "CellEngineDataFile.h"

class CellEnginePDBDataFile: public CellEngineDataFile
{
private:
    void ReadDataFromFile(bool StartValuesBool) override;
private:
    void SetStartValues() override
    {
    }
    UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) override
    {
        return 0;
    }
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, const UniqueIdInt ParticleIndex) override
    {
    }
    void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) override
    {
    }
    void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) override
    {
    }
    void PreprocessData() override
    {
    }
    void PrintStatistics() override
    {
    }
public:
    explicit CellEnginePDBDataFile() = default;
    ~CellEnginePDBDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

#endif
