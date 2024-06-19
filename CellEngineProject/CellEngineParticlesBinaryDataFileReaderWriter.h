
#ifndef CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H
#define CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H

#include "CellEngineDataFile.h"
#include "CellEngineBuildParticlesDataOperations.h"

class CellEngineParticlesBinaryDataFileReaderWriter : virtual public CellEngineDataFile, virtual public CellEngineBuildParticlesDataOperations//, virtual public CellEngineChemicalReactions
{
public:
    explicit CellEngineParticlesBinaryDataFileReaderWriter() = default;
    ~CellEngineParticlesBinaryDataFileReaderWriter() = default;
public:
    void SaveDataToFile() override;
public:
    void SaveParticlesToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveParticlesKindsToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveGenesToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveChemicalReactionsToBinaryFile(std::ofstream& ParticlesDataFile);
    void SaveParticlesKindsAndParticlesAndChemicalReactionsAndGenesToBinaryFile();
public:
    void ReadParticlesFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadParticlesKindsFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadGenesFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadChemicalReactionsFromBinaryFile(std::ifstream& ParticlesDataFile);
    void ReadParticlesKindsAndParticlesAndChemicalReactionsAndGenesFromBinaryFile(CellEngineConfigData::TypesOfFileToRead Type);
    void PrepareParticlesAfterReadingFromBinaryFile();
    void ReadAllDataFromBinaryFileAndPrepareData(const bool StartValuesBool, const bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type);
};


#endif
