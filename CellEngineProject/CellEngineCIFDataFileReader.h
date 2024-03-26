#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include "CellEngineDataFile.h"

struct TransformationMatrix3x4
{
    float Matrix[3][4];
};

class AAA
{
protected:
    virtual UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) = 0;
    virtual void SetStartValues() = 0;
    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, UniqueIdInt ParticleIndex) = 0;
    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
    virtual void PreprocessData() = 0;
    virtual void PrintStatistics() = 0;
};

class CellEngineCIFDataFileReader : virtual public CellEngineDataFile, virtual public AAA
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


class CellEngineParticlesBinaryDataFileReader : virtual public CellEngineDataFile, virtual public AAA
{
public:
    explicit CellEngineParticlesBinaryDataFileReader() = default;
    ~CellEngineParticlesBinaryDataFileReader() = default;
private:
    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
    {
        return (*Particles)[ParticleIndex];
    }
public:
    void ReadParticlesFromFile();
    void PrepareParticlesAfterReadingFromFile();
    void ReadParticlesFromBinaryFileAndPrepareData(bool StartValuesBool);
};

class CellEngineParticlesDataFileReader : public CellEngineCIFDataFileReader, public CellEngineParticlesBinaryDataFileReader
{
public:
    void ReadDataFromFile(bool StartValuesBool, CellEngineConfigData::TypesOfFileToRead Type) override;
public:
    explicit CellEngineParticlesDataFileReader() = default;
    ~CellEngineParticlesDataFileReader() = default;
};


//class CellEngineParticlesDataFileReader : virtual public CellEngineDataFile
//class CellEngineParticlesDataFileReader : public CellEngineDataFile
//{
//public:
//    void ReadDataFromFile(bool StartValuesBool, CellEngineConfigData::TypesOfFileToRead Type) override;
//public:
//    explicit CellEngineParticlesDataFileReader() = default;
//    ~CellEngineParticlesDataFileReader() = default;
//private:
//    inline Particle& GetParticleFromIndex(const UniqueIdInt ParticleIndex)
//    {
//        return (*Particles)[ParticleIndex];
//    }
//public:
//    void ReadDataFromCIFFile();
//    void ReadParticlesFromFile();
//    void PrepareParticlesAfterReadingFromFile();
//    void ReadParticlesFromBinaryFileAndPrepareData(bool StartValuesBool);
//
//private:
//    std::unordered_map<UnsignedInt, TransformationMatrix3x4> TransformationsMatrixes;
//    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
//public:
//    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
//
//protected:
//    virtual UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) = 0;
//    virtual void SetStartValues() = 0;
//    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, UniqueIdInt ParticleIndex) = 0;
//    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
//    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
//    virtual void PreprocessData() = 0;
//    virtual void PrintStatistics() = 0;
//};

#endif