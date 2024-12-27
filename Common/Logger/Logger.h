#pragma once

#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <mutex>
#include <thread>

#include <sstream>
#include <fstream>

#include <vector>
#include <unordered_map>

#include <functional>
#include <initializer_list>

using ThreadIdType = std::uint64_t;

class Logger
{
public:
	ThreadIdType ThisThreadId;
private:
	std::uint64_t FileNumberInLog = 1;
	std::uint64_t LineNumberInLog = 1;
private:
    std::string LogDirectory;
	std::string MainDirectoryName;
	std::string LoggerName;
	std::string TaskName;
private:
	std::vector<std::ofstream> UserLogFiles;
	std::vector<std::ofstream> SpecialLogFiles;
private:
	static inline std::mutex LogMessageCoutMutexObject;
public:
	Logger(const char* LogDirectoryParameter, const char* MainDirectoryNameParameter, const char* LoggerNameParameter, const char* TaskNameParameter, std::uint64_t ThisThreadIdParameter);
	~Logger();
private:
	void AllocResourcesForFiles();
	void CreateDirectories();
	void OpenLogFiles();
	void CloseLogFiles();
private:
	void CloseOldLogFilesAndOpenNewLogFiles();
	void LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(ThreadIdType CurrentThreadId);
private:
	[[nodiscard]] std::string CreateLogString(const std::string& MessageStr, bool LogLineInfo, ThreadIdType CurrentThreadId, std::uint64_t LineNumberInCommonLog, bool PrintLogLineNumber, bool PrintLogDateTime, bool PrintLogProcessId, bool PrintLogProcessPriorityLevel, bool PrintLogThreadId) const;
private:
	static void WriteToCommonLogFromThread(bool Condition, const std::string& MessageStr, std::ostream& StreamObject, ThreadIdType CurrentThreadId, std::uint64_t FileNumber);
	void WriteToLogsFromThread(const std::string& MessageStrToFile, ThreadIdType CurrentThreadId, std::int64_t SpecialLogFileIndex);
public:
	void LogMessageBool(const std::string& MessageStr, bool LogLineInfo, ThreadIdType CurrentThreadId, bool PrintLogToConsoleUnconditional, bool PrintLogToFilesUnconditional, bool PrintLogToConsole, bool PrintLogToFiles, std::int64_t SpecialLogFileIndex);
};

class LoggersManager
{
	friend class Logger;
private:
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

	std::uint64_t MaximalNumberOfLinesInOneFile = 100000;

    bool PrintLogToCommonFileWhenPrintLogToSpecialFile = false;
private:
    std::vector<bool> UseSpecialLogFiles = { false, false, false, false, false, false, false, false, false };

    const std::int64_t LogWarningsFileIndex = 0;
    const std::int64_t LogErrorsFileIndex = 1;
    const std::int64_t LogExceptionsFileIndex = 2;
    const std::int64_t LogErrorsAndExceptionsFileIndex = 3;
    const std::int64_t LogCriticalFileIndex = 4;
    const std::int64_t LogInformationFileIndex = 5;
    const std::int64_t LogImportantFileIndex = 6;
    const std::int64_t LogStatisticsFileIndex = 7;
    const std::int64_t LogDebugFileIndex = 8;
private:
	std::string TaskName;
    std::string LogDirectory;
	std::string ActualDateTimeStr;
    bool LogThreadsToSeparatedFiles{};
	std::uint64_t FileNumberToIncreaseLineNumber{};
private:
	std::function<void(const ThreadIdType CurrentThreadId, const std::uint64_t FileNumber, const std::string& MessageStr)> DrawMessageFunctionObject;
private:
	std::unique_ptr<Logger> LoggerMainObjectPointer;
	std::unordered_map<ThreadIdType, std::unique_ptr<Logger>> LoggersThreadsObjectsPointersMap;
private:
	static inline std::mutex CreateNewLoggerForThreadMutexObject;
private:
	void LogMessageBool(const std::string& MessageStr, bool LogLineInfo, bool PrintLogToConsoleUnconditional, bool PrintLogToFilesUnconditional, bool PrintToConsole, bool PrintLogToFiles, std::int64_t SpecialLogFileIndex);
public:
	void Log(const std::stringstream& Message);

	void LogOnlyToFiles(const std::stringstream& Message);
	void LogOnlyToConsole(const std::stringstream& Message);
	void LogWithoutLineInfo(const std::stringstream& Message);
	void LogWithoutLineInfoOnlyToFiles(const std::stringstream& Message);

	void LogUnconditional(const std::stringstream& Message);
	void LogOnlyToConsoleUnconditional(const std::stringstream& Message);
	void LogOnlyToFilesUnconditional(const std::stringstream& Message);

    void LogInColorTerminal(std::ostream& color(std::ostream& s), const std::stringstream& Message);
public:
    [[maybe_unused]] void LogWarning(const std::stringstream& Message);
    [[maybe_unused]] void LogError(const std::stringstream& Message);
    [[maybe_unused]] void LogException(const std::stringstream& Message);
    [[maybe_unused]] void LogErrorAndException(const std::stringstream& Message);
    [[maybe_unused]] void LogCritical(const std::stringstream& Message);
    [[maybe_unused]] void LogInformation(const std::stringstream& Message);
    [[maybe_unused]] void LogImportant(const std::stringstream& Message);
    [[maybe_unused]] void LogStatistics(const std::stringstream& Message);
    [[maybe_unused]] void LogDebug(const std::stringstream& Message);
public:
    [[maybe_unused]] void LogWarn(const std::stringstream& Message);
    [[maybe_unused]] void LogErr(const std::stringstream& Message);
    [[maybe_unused]] void LogExc(const std::stringstream& Message);
    [[maybe_unused]] void LogErrAndExc(const std::stringstream& Message);
    [[maybe_unused]] void LogCrit(const std::stringstream& Message);
    [[maybe_unused]] void LogInfo(const std::stringstream& Message);
    [[maybe_unused]] void LogImp(const std::stringstream& Message);
    [[maybe_unused]] void LogStat(const std::stringstream& Message);
    [[maybe_unused]] void LogDeb(const std::stringstream& Message);
private:
    std::vector<std::string> SpecialLogFilesNames = { "Warnings", "Errors", "Exceptions", "ErrorsAndExceptions", "Critical", "Information", "Important", "Statistics", "Debug" };
	std::vector<std::string> UserLogFilesNames;
	std::vector<std::string> SelectiveWordsForUserLogFiles;
	std::vector<std::function<const bool(const std::string&)>> SelectiveWordsFunctions;
public:	
	void InitializeFilesNames(std::initializer_list<const std::string> InitialFilesNames);
	void InitializeSelectiveWordsFunctions(std::initializer_list<std::function<bool(const std::string&)>> InitialSelectiveWordsFunctions);
	void InitializePrintingParameters(bool PrintLogToConsoleParam, bool PrintLogToFilesParam, bool PrintLogLineNumberToConsoleParam, bool PrintLogDateTimeToConsoleParam, bool PrintLogProcessIdToConsoleParam, bool PrintLogProcessPriorityLevelToConsoleParam, bool PrintLogThreadIdToConsoleParam, bool PrintLogLineNumberToFileParam, bool PrintLogDateTimeToFileParam, bool PrintLogProcessIdToFileParam, bool PrintLogProcessPriorityLevelToFileParam, bool PrintLogThreadIdToFileParam, uint64_t MaximalNumberOfLinesInOneFileParam, bool PrintLogToCommonFileWhenPrintLogToSpecialFile);
	void InitializeSpecialLogFiles(bool CreateLogWarningsFileParam, bool CreateLogErrorsFileParam, bool CreateLogExceptionsFileParam, bool CreateLogErrorsAndExceptionsFileParam, bool CreateLogCriticalFileParam, bool CreateLogInformationFileParam, bool CreateLogImportantFileParam, bool CreateLogStatisticsFileParam, bool CreateLogDebugFileParam);
	void InitializeLoggerManagerDataForTask(const std::string& TaskNameParameter, const std::string& LogDirectoryParameter, const std::string& ActualDateTimeStrParameter, bool LogThreadsToSeparatedFilesParameter, std::uint64_t FileNumberToIncreaseLineNumberParameter, std::function<void(const ThreadIdType CurrentThreadId, const std::uint64_t FileNumber, const std::string& MessageStr)> DrawMessageFunctionObjectParameter);
};

inline LoggersManager LoggersManagerObject;

#define STREAM(text) (std::stringstream&)(std::stringstream() << text)

#endif
