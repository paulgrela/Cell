
#ifndef CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H
#define CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H

#include "CellEngineDataFile.h"
#include "CellEngineBuildParticlesDataOperations.h"

class CellEngineParticlesBinaryDataFileReaderWriter : virtual public CellEngineDataFile, virtual public CellEngineBuildParticlesDataOperations
{
public:
    explicit CellEngineParticlesBinaryDataFileReaderWriter() = default;
    ~CellEngineParticlesBinaryDataFileReaderWriter() override = default;
public:
    void SaveDataToFile() override;
public:
    void SaveParticlesToBinaryFile(std::ofstream& ParticlesDataFile) const;
    static void SaveParticlesKindsToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveGenesToBinaryFile(std::ofstream& ParticlesDataFile);
    static void SaveChemicalReactionsToBinaryFile(std::ofstream& ParticlesDataFile);
    void SaveParticlesKindsAndParticlesAndChemicalReactionsAndGenesToBinaryFile() const;
public:
    void ReadParticlesFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadParticlesKindsFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadGenesFromBinaryFile(std::ifstream& ParticlesDataFile);
    static void ReadChemicalReactionsFromBinaryFile(std::ifstream& ParticlesDataFile);
    void PrepareParticlesAfterReadingFromBinaryFile();
    static void PreprocessLinkAndAssociateEveryParticleKindWithProperChemicalReaction();
    void ReadParticlesKindsAndParticlesAndChemicalReactionsAndGenesFromBinaryFile(CellEngineConfigData::TypesOfFileToRead Type);
    void ReadAllDataFromBinaryFileAndPrepareData(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type);
public:
    static void FindNucleotidesIdentifiers();
    void FindGenomeParameters() const;
};


#endif
