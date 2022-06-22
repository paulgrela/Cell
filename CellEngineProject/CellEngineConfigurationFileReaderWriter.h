#pragma once

#ifndef _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_
#define _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_

class CellEngineConfigurationFileReaderWriter
{
public:
    std::string ConfigFileName;
public:
    CellEngineConfigurationFileReaderWriter() = default;
public:
    void ReadChessConfigurationFile(const char* ConfigFileNameParameter, const std::uint64_t ExecuteCellStateId);
    void SaveTestStatisticsToFile(const uint64_t ExecuteCellStateId) const;
};

inline CellEngineConfigurationFileReaderWriter CellEngineConfigurationFileReaderWriterObject;

#endif