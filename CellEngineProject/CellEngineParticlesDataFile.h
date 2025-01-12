
#ifndef CELL_ENGINE_PARTICLES_DATA_FILE_H
#define CELL_ENGINE_PARTICLES_DATA_FILE_H

#include "CellEngineCIFDataFileReader.h"
#include "CellEngineParticlesBinaryDataFileReaderWriter.h"

class CellEngineParticlesDataFile : public CellEngineCIFDataFileReader, public CellEngineParticlesBinaryDataFileReaderWriter
{
public:
    void ReadDataFromFile(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type) override;
    void ReadCIFFileWithBinaryData(bool StartValuesBool, bool UpdateParticleKindListOfVoxelsBool, CellEngineConfigData::TypesOfFileToRead Type);
public:
    explicit CellEngineParticlesDataFile() = default;
    ~CellEngineParticlesDataFile() override = default;
};

#endif
