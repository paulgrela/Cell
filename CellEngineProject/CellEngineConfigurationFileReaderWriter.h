#pragma once

#ifndef _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_
#define _CELL_ENGINE_CONFIGURATION_FILE_READER_WRITER_H_

#include <string>
#include <vector>
#include <atomic>

#include "CellEngineTypes.h"

class CellEngineConfigurationFileReaderWriter
{
public:
    bool MultiThreaded;
    bool SetProcessPriorityHighest;
public:
    bool PrintLogToConsole = true;
    bool PrintLogToFiles = true;

    bool PrintLogLineNumberToConsole = false;
    bool PrintLogDateTimeToConsole = false;
    bool PrintLogProcessIdToConsole = false;
    bool PrintLogProcessPriorityLevelToConsole = false;
    bool PrintLogThreadIdToConsole = false;

    bool PrintLogLineNumberToFile = true;
    bool PrintLogDateTimeToFile = true;
    bool PrintLogProcessIdToFile = false;
    bool PrintLogProcessPriorityLevelToFile = false;
    bool PrintLogThreadIdToFile = false;

    uint64_t MaximalNumberOfLinesInOneFile = 100000;
public:
    std::string ConfigFileName;
public:
    CellEngineConfigurationFileReaderWriter() = default;
public:
    void ReadChessConfigurationFile(const char* ConfigFileNameParameter, std::unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer, const std::uint64_t ExecuteCellStateId);
    void SaveTestStatisticsToFile(std::unique_ptr<CellEngineDataFile>& CellEngineDataFileObjectPointer, const uint64_t ExecuteCellStateId) const;
};

inline CellEngineConfigurationFileReaderWriter CellEngineConfigurationFileReaderWriterObject;

#endif