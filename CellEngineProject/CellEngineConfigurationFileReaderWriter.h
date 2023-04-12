#pragma once

#ifndef CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H
#define CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H

class CellEngineConfigurationFileReaderWriter
{
public:
    std::string ConfigFileName;
public:
    CellEngineConfigurationFileReaderWriter() = default;
public:
    void ReadCellConfigurationFile(const char* ConfigFileNameParameter, UnsignedInt ExecuteCellStateId);
    void SaveTestStatisticsToFile(UnsignedInt ExecuteCellStateId) const;
};

inline CellEngineConfigurationFileReaderWriter CellEngineConfigurationFileReaderWriterObject;

#endif