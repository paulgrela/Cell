#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include "CellEngineDataFile.h"

struct TransformationMatrix3x4
{
    float Matrix[3][4];
};

class CellEngineCIFDataFile : public CellEngineDataFile
{
private:
    std::unordered_map<UnsignedInt, TransformationMatrix3x4> TransformationsMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
private:
    void ReadDataFromFile() override;
public:
    explicit CellEngineCIFDataFile() = default;
    ~CellEngineCIFDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
protected:
    virtual void SetStartValues() = 0;
    virtual void SetParticleKindData(const EntityIdInt EntityId, const ChainIdInt ChainId) = 0;
    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom) = 0;
    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
    virtual void PrintStatistics() = 0;
};

#endif