
#ifndef CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H
#define CELL_ENGINE_PARTICLES_BINARY_DATA_FILE_READER_WRITER_H

#include "CellEngineDataFile.h"
#include "CellEngineBuildParticlesDataOperations.h"

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
    void ReadParticlesKindsAndParticlesFromBinaryFile(CellEngineConfigData::TypesOfFileToRead Type);
    void PrepareParticlesAfterReadingFromBinaryFile();
    void ReadParticlesFromBinaryFileAndPrepareData(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type);
};


#endif
