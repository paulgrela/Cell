#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include "CellEngineDataFile.h"

class CellEngineBuildParticlesDataOperations
{
protected:
    virtual void SetStartValues() = 0;
    virtual UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) = 0;
    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, UniqueIdInt ParticleIndex) = 0;
    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
    virtual void PreprocessData() = 0;
    virtual void PrintStatistics() = 0;
};

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

class CellEngineParticlesBinaryDataFileReaderWriter : virtual public CellEngineDataFile, virtual public CellEngineBuildParticlesDataOperations
{
public:
    explicit CellEngineParticlesBinaryDataFileReaderWriter() = default;
    ~CellEngineParticlesBinaryDataFileReaderWriter() = default;
private:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return (*Particles)[ParticleIndex];
    }
public:
    void SaveDataToFile() override;
public:
    void SaveParticlesToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveParticlesKindsToBinaryFile(std::ofstream& ParticlesDataFile);
    void SaveParticlesKindsAndParticlesToBinaryFile();
public:
    void ReadParticlesFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadParticlesKindsFromBinaryFile(std::ifstream& ParticlesDataFile);
    void ReadParticlesKindsAndParticlesFromBinaryFile();
    void PrepareParticlesAfterReadingFromBinaryFile();
    void ReadParticlesFromBinaryFileAndPrepareData(bool StartValuesBool);
};

class CellEngineParticlesDataFile : public CellEngineCIFDataFileReader, public CellEngineParticlesBinaryDataFileReaderWriter
{
public:
    void ReadDataFromFile(bool StartValuesBool, CellEngineConfigData::TypesOfFileToRead Type) override;
public:
    explicit CellEngineParticlesDataFile() = default;
    ~CellEngineParticlesDataFile() = default;
};

#endif