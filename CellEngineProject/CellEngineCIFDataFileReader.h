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
    void ReadDataFromFile(bool StartValuesBool) override;
public:
    explicit CellEngineCIFDataFileReader() = default;
    ~CellEngineCIFDataFileReader() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
};

class CellEngineParticlesDataFileReader : virtual public CellEngineDataFile
{
protected:
    void ReadDataFromFile(bool StartValuesBool) override;
public:
    explicit CellEngineParticlesDataFileReader() = default;
    ~CellEngineParticlesDataFileReader() = default;
private:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return (*Particles)[ParticleIndex];
    }
public:
    void ReadParticlesFromFile();
    void PrepareParticlesAfterReadingFromFile();
    void ReadParticlesFromFileAndPrepareData(bool StartValuesBool);
};

#endif