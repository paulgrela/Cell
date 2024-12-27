
#include "DestinationPlatform.h"

#ifdef WINDOWS_PLATFORM
#include <windows.h>
#endif

#include <limits>
#include <string>
#include <iostream>

#ifdef WINDOWS_PLATFORM
#include <dir.h>
#endif

#ifdef UNIX_PLATFORM
#include <unistd.h>
#include <sys/resource.h>
#include <filesystem>
#endif

#include "ExceptionsMacro.h"

#include "Logger.h"
#include "StringUtils.h"
#include "DateTimeUtils.h"
#include "TerminalColorsUtils.h"

using namespace std;
using namespace terminal_colors_utils;

constexpr char EndLineChar = '\n';

Logger::Logger(const char* LogDirectoryParameter, const char* MainDirectoryNameParameter, const char* LoggerNameParameter, const char* TaskNameParameter, const uint64_t ThisThreadIdParameter)
{
	try
	{
		ThisThreadId = ThisThreadIdParameter;
        LogDirectory = LogDirectoryParameter;
		MainDirectoryName = MainDirectoryNameParameter;
		LoggerName = LoggerNameParameter;
		TaskName = TaskNameParameter;

		CreateDirectories();
		AllocResourcesForFiles();
		OpenLogFiles();
	}
	CATCH_COUT("logger constructor")
}

Logger::~Logger()
{
	CloseLogFiles();
}

void Logger::CreateDirectories()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
		{
		    #ifdef WINDOWS_PLATFORM
			mkdir(string(string(LogDirectory) + string("logs")).c_str());
			mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName).c_str());
			mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName).c_str());
			mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName).c_str());

			for (const string& FileName : LoggersManagerObject.UserLogFilesNames)
				mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + FileName).c_str());

			for (const string& FileName : LoggersManagerObject.SpecialLogFilesNames)
                if (LoggersManagerObject.UseSpecialLogFiles[FileNumber] == true)
				    mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + FileName).c_str());
            #endif
            #ifdef UNIX_PLATFORM
            filesystem::create_directory(string(LogDirectory) + string("logs"));
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName);
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName);
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName);

            for (const string& FileName : LoggersManagerObject.UserLogFilesNames)
                filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + FileName);

            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.SpecialLogFilesNames.size(); FileNumber++)
                if (LoggersManagerObject.UseSpecialLogFiles[FileNumber] == true)
                    filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + LoggersManagerObject.SpecialLogFilesNames[FileNumber]);
            #endif
        }
	}
	CATCH_AND_THROW_COUT("creating directories in logger")
}

void Logger::AllocResourcesForFiles()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
		{
			UserLogFiles.clear();
			UserLogFiles.resize(LoggersManagerObject.UserLogFilesNames.size());

            SpecialLogFiles.clear();
            SpecialLogFiles.resize(LoggersManagerObject.UseSpecialLogFiles.size());
		}
	}
	CATCH_AND_THROW_COUT("allocation of  resources for file in logger")
}

void Logger::OpenLogFiles()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
        {
            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.UserLogFilesNames.size(); FileNumber++)
                UserLogFiles[FileNumber].open(string(LogDirectory) + OS_DIR_SEP + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + LoggersManagerObject.UserLogFilesNames[FileNumber] + OS_DIR_SEP + LoggersManagerObject.UserLogFilesNames[FileNumber] + string_utils::align_str(to_string(FileNumberInLog), '0', 5) + ".log.txt");

            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.SpecialLogFilesNames.size(); FileNumber++)
                if (LoggersManagerObject.UseSpecialLogFiles[FileNumber] == true)
                    SpecialLogFiles[FileNumber].open(string(LogDirectory) + OS_DIR_SEP + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + LoggersManagerObject.SpecialLogFilesNames[FileNumber] + OS_DIR_SEP + LoggersManagerObject.SpecialLogFilesNames[FileNumber] + string_utils::align_str(to_string(FileNumberInLog), '0', 5) + ".log.txt");
        }
	}
	CATCH_AND_THROW_COUT("opening log files in logger")
}

void Logger::CloseLogFiles()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
        {
			for (auto& FileNumber : UserLogFiles)
				FileNumber.close();

            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.SpecialLogFilesNames.size(); FileNumber++)
                if (LoggersManagerObject.UseSpecialLogFiles[FileNumber] == true)
                    SpecialLogFiles[FileNumber].close();
        }
	}
	CATCH_AND_THROW_COUT("closing log files in logger")
}

void Logger::CloseOldLogFilesAndOpenNewLogFiles()
{
	try
	{
		CloseLogFiles();

		FileNumberInLog++;

		OpenLogFiles();
	}
	CATCH_AND_THROW_COUT("closing old log files and opening new log files after maximal limit of lines in file is exceeded in logger")
}

void Logger::LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(const ThreadIdType CurrentThreadId)
{
	try
	{
		if (LoggersManagerObject.MaximalNumberOfLinesInOneFile != 0 && LineNumberInLog != 0 && LineNumberInLog % LoggersManagerObject.MaximalNumberOfLinesInOneFile == 0)
			CloseOldLogFilesAndOpenNewLogFiles();
	}
	CATCH_COUT("limiting log size by closing old log files and opening new log files after maximal limit of lines in old file is exceeded in logger")
}

string Logger::CreateLogString(const string& MessageStr, const bool LogLineInfo, const ThreadIdType CurrentThreadId, uint64_t LineNumberInCommonLog, const bool PrintLogLineNumber, const bool PrintLogDateTime, const bool PrintLogProcessId, const bool PrintLogProcessPriorityLevel, const bool PrintLogThreadId) const
{
	using namespace string_utils;

	string LocalMessageStr;

	try
	{
		if (LogLineInfo == true)
		{
			if (PrintLogLineNumber == true)
			{
				if (LoggersManagerObject.LoggerMainObjectPointer.get() == this)
					LocalMessageStr = LocalMessageStr + "[" + align_str(to_string(LineNumberInCommonLog), '0', 5) + "] ";
				else
					LocalMessageStr = LocalMessageStr + "[" + align_str(to_string(LineNumberInLog), '0', 5) + "] ";

				if (LoggersManagerObject.LoggerMainObjectPointer)
					if (LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId != CurrentThreadId && LoggersManagerObject.LoggerMainObjectPointer.get() != this)
						LocalMessageStr = LocalMessageStr + "[" + align_str(to_string(LineNumberInCommonLog), '0', 5) + "] ";
			}

            #ifdef WINDOWS_PLATFORM
            if (PrintLogDateTime == true)
				LocalMessageStr = LocalMessageStr + "[" + GetActualDateTimeWindows("-", "-", " ", ":", ":", ":") + "] ";

			if (PrintLogProcessId == true)
				LocalMessageStr = LocalMessageStr + "[ProcessId = " + to_string(int64_t(GetCurrentProcessId())) + "] ";

			if (PrintLogProcessPriorityLevel == true)
				LocalMessageStr = LocalMessageStr + "[ProcessPriorityLevel = " + to_string(int64_t(GetPriorityClass(GetCurrentProcess()))) + "] ";
            #endif

            #ifdef UNIX_PLATFORM
            if (PrintLogDateTime == true)
                LocalMessageStr = LocalMessageStr + "[" + GetActualDateTimeStandardCPP("-", "-", " ", ":", ":") + "] ";

            if (PrintLogProcessId == true)
                LocalMessageStr = LocalMessageStr + "[ProcessId = " + to_string(int64_t(getpid())) + "] ";

            if (PrintLogProcessPriorityLevel == true)
                LocalMessageStr = LocalMessageStr + "[ProcessPriorityLevel = " + to_string(int64_t(getpriority(PRIO_PROCESS, getpid()))) + "] ";
            #endif

			if (PrintLogThreadId == true)
				LocalMessageStr = LocalMessageStr + "[ThreadId = " + to_string(CurrentThreadId) + "] ";
		}

		LocalMessageStr = LocalMessageStr + MessageStr + EndLineChar;
	}
	CATCH_AND_THROW_COUT("creating log string in logger");

	return LocalMessageStr;
}

void Logger::WriteToCommonLogFromThread(const bool Condition, const string& MessageStr, ostream& StreamObject, const ThreadIdType CurrentThreadId, const uint64_t FileNumber)
{
	try
	{
		if (Condition == true)
		{
			StreamObject << MessageStr << flush;

			if (FileNumber != numeric_limits<uint64_t>::max())
				if (LoggersManagerObject.DrawMessageFunctionObject)
					LoggersManagerObject.DrawMessageFunctionObject(CurrentThreadId, FileNumber, MessageStr);
		}
	}
	CATCH_AND_THROW_COUT("writing to common log from thread in logger")
}

void Logger::WriteToLogsFromThread(const string& MessageStrToFile, const ThreadIdType CurrentThreadId, const std::int64_t SpecialLogFileIndex)
{
	try
	{
        if (SpecialLogFileIndex == -1)
            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.UserLogFilesNames.size(); FileNumber++)
                WriteToCommonLogFromThread(UserLogFiles[FileNumber].is_open() == true && LoggersManagerObject.SelectiveWordsFunctions[FileNumber](MessageStrToFile), MessageStrToFile, UserLogFiles[FileNumber], CurrentThreadId, FileNumber);
        else
        if (LoggersManagerObject.UseSpecialLogFiles[SpecialLogFileIndex] == true)
        {
            if (LoggersManagerObject.PrintLogToCommonFileWhenPrintLogToSpecialFile == true)
                WriteToCommonLogFromThread(UserLogFiles[0].is_open() == true && LoggersManagerObject.SelectiveWordsFunctions[0](MessageStrToFile), MessageStrToFile, UserLogFiles[0], CurrentThreadId, 0);
            WriteToCommonLogFromThread(SpecialLogFiles[SpecialLogFileIndex].is_open() == true, MessageStrToFile, SpecialLogFiles[SpecialLogFileIndex], CurrentThreadId, SpecialLogFileIndex);
        }
	}
	CATCH_AND_THROW_COUT("writing to logs from thread in logger")
}

void Logger::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const ThreadIdType CurrentThreadId, const bool PrintLogToConsoleUnconditional, const bool PrintLogToFilesUnconditional, const bool PrintLogToConsole, const bool PrintLogToFiles, const std::int64_t SpecialLogFileIndex)
{
	try
	{
		string LocalMessageStr;

		lock_guard LockGuardObject{ LogMessageCoutMutexObject };

		uint64_t LineNumberInCommonLog = LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog;
		LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog++;

		if ((LoggersManagerObject.PrintLogToConsole == true && PrintLogToConsole == true) || PrintLogToConsoleUnconditional == true)
		{
			LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToConsole, LoggersManagerObject.PrintLogDateTimeToConsole, LoggersManagerObject.PrintLogProcessIdToConsole, LoggersManagerObject.PrintLogProcessPriorityLevelToConsole, LoggersManagerObject.PrintLogThreadIdToConsole);
			WriteToCommonLogFromThread(true, LocalMessageStr, cout, CurrentThreadId, numeric_limits<uint64_t>::max());
		}

		if ((LoggersManagerObject.PrintLogToFiles == true && PrintLogToFiles == true) || PrintLogToFilesUnconditional == true)
			if (LoggersManagerObject.LoggerMainObjectPointer)
			{
				if (LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId != CurrentThreadId)
				{
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					LoggersManagerObject.LoggerMainObjectPointer->WriteToLogsFromThread(LocalMessageStr, LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId, SpecialLogFileIndex);

					LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(CurrentThreadId);

					LineNumberInLog++;
					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId, SpecialLogFileIndex);
				}
				else
				{
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId, SpecialLogFileIndex);
				}
			}
	}
	CATCH_COUT("logging message in logger")
}

void  LoggersManager::InitializeFilesNames(const initializer_list<const string> InitialFilesNames)
{
	try
	{
		for (const string& InitialFileName : InitialFilesNames)
			UserLogFilesNames.push_back(InitialFileName);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager files names")
}

void  LoggersManager::InitializeSelectiveWordsFunctions(const initializer_list<function<bool(const string&)>> InitialSelectiveWordsFunctions)
{
	try
	{
		for (const auto& InitialSelectiveWordsFunction : InitialSelectiveWordsFunctions)
			SelectiveWordsFunctions.emplace_back(InitialSelectiveWordsFunction);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager selective words functions")
}

void LoggersManager::InitializePrintingParameters(bool PrintLogToConsoleParam, bool PrintLogToFilesParam, bool PrintLogLineNumberToConsoleParam, bool PrintLogDateTimeToConsoleParam, bool PrintLogProcessIdToConsoleParam, bool PrintLogProcessPriorityLevelToConsoleParam, bool PrintLogThreadIdToConsoleParam, bool PrintLogLineNumberToFileParam, bool PrintLogDateTimeToFileParam, bool PrintLogProcessIdToFileParam, bool PrintLogProcessPriorityLevelToFileParam, bool PrintLogThreadIdToFileParam, uint64_t MaximalNumberOfLinesInOneFileParam, bool PrintLogToCommonFileWhenPrintLogToSpecialFile)
{
	try
	{
		this->PrintLogToConsole = PrintLogToConsoleParam;
		this->PrintLogToFiles = PrintLogToFilesParam;
	
		this->PrintLogLineNumberToConsole = PrintLogLineNumberToConsoleParam;
		this->PrintLogDateTimeToConsole = PrintLogDateTimeToConsoleParam;
		this->PrintLogProcessIdToConsole = PrintLogProcessIdToConsoleParam;
		this->PrintLogProcessPriorityLevelToConsole = PrintLogProcessPriorityLevelToConsoleParam;
		this->PrintLogThreadIdToConsole = PrintLogThreadIdToConsoleParam;

		this->PrintLogLineNumberToFile = PrintLogLineNumberToFileParam;
		this->PrintLogDateTimeToFile = PrintLogDateTimeToFileParam;
		this->PrintLogProcessIdToFile = PrintLogProcessIdToFileParam;
		this->PrintLogProcessPriorityLevelToFile = PrintLogProcessPriorityLevelToFileParam;
		this->PrintLogThreadIdToFile = PrintLogThreadIdToFileParam;

		this->MaximalNumberOfLinesInOneFile = MaximalNumberOfLinesInOneFileParam;

        this->PrintLogToCommonFileWhenPrintLogToSpecialFile = PrintLogToCommonFileWhenPrintLogToSpecialFile;
	}
	CATCH_AND_THROW_COUT("initializing loggers manager printing conditions")
};

void LoggersManager::InitializeSpecialLogFiles(bool CreateLogWarningsFileParam, bool CreateLogErrorsFileParam, bool CreateLogExceptionsFileParam, bool CreateLogErrorsAndExceptionsFileParam, bool CreateLogCriticalFileParam, bool CreateLogInformationFileParam, bool CreateLogImportantFileParam, bool CreateLogStatisticsFileParam, bool CreateLogDebugFileParam)
{
    try
    {
        UseSpecialLogFiles[LogWarningsFileIndex] = CreateLogWarningsFileParam;
        UseSpecialLogFiles[LogErrorsFileIndex] = CreateLogErrorsFileParam;
        UseSpecialLogFiles[LogExceptionsFileIndex] = CreateLogExceptionsFileParam;
        UseSpecialLogFiles[LogErrorsAndExceptionsFileIndex] = CreateLogErrorsAndExceptionsFileParam;
        UseSpecialLogFiles[LogCriticalFileIndex] = CreateLogCriticalFileParam;
        UseSpecialLogFiles[LogInformationFileIndex] = CreateLogInformationFileParam;
        UseSpecialLogFiles[LogImportantFileIndex] = CreateLogImportantFileParam;
        UseSpecialLogFiles[LogStatisticsFileIndex] = CreateLogStatisticsFileParam;
        UseSpecialLogFiles[LogDebugFileIndex] = CreateLogDebugFileParam;
    }
    CATCH("initializing special log files")
}

void LoggersManager::InitializeLoggerManagerDataForTask(const string& TaskNameParameter, const std::string& LogDirectoryParameter, const string& ActualDateTimeStrParameter, const bool LogThreadsToSeparatedFilesParameter, const uint64_t FileNumberToIncreaseLineNumberParameter, function<void(const ThreadIdType CurrentThreadId, const uint64_t FileNumber, const string& MessageStr)> DrawMessageFunctionObjectParameter)
{
	try
	{
		LoggerMainObjectPointer.reset();
		LoggersThreadsObjectsPointersMap.clear();

		TaskName = TaskNameParameter;
        LogDirectory = LogDirectoryParameter;
		ActualDateTimeStr = ActualDateTimeStrParameter;
        LogThreadsToSeparatedFiles = LogThreadsToSeparatedFilesParameter;
		FileNumberToIncreaseLineNumber = FileNumberToIncreaseLineNumberParameter;

		if (LoggersManagerObject.FileNumberToIncreaseLineNumber >= LoggersManagerObject.UserLogFilesNames.size())
			throw runtime_error("FileNumberToIncreaseLineNumber is out of range of user defined files.");

		DrawMessageFunctionObject = std::move(DrawMessageFunctionObjectParameter);

        ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str());
		LoggerMainObjectPointer = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("LOGGER_COMMON").c_str(), TaskName.c_str(), CurrentThreadId);
		LoggerMainObjectPointer->LogMessageBool("START MAIN LOGGER_COMMON\n", true, CurrentThreadId, false, false, true, true, -1);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager data for task")
}

void LoggersManager::Log(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, false, false, true, true, -1);
};

void LoggersManager::LogWithoutLineInfo(const stringstream& Message)
{
	LogMessageBool(Message.str(), false, false, false, true, true, -1);
}

void LoggersManager::LogOnlyToFiles(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, false, false, false, true, -1);
};

void LoggersManager::LogOnlyToConsole(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, false, false, true, false, -1);
};

void LoggersManager::LogWithoutLineInfoOnlyToFiles(const stringstream& Message)
{
	LogMessageBool(Message.str(), false, false, false, false, true, -1);
}

void LoggersManager::LogUnconditional(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, true, true, false, false, -1);
};

void LoggersManager::LogOnlyToConsoleUnconditional(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, true, false, false, false, -1);
};

void LoggersManager::LogOnlyToFilesUnconditional(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, false, true, false, false, -1);
};

void LoggersManager::LogInColorTerminal(ostream& color(ostream& s), const stringstream& Message)
{
    this->LogOnlyToFiles(Message);
    cout << color << Message.str() << white << endl;
}

void LoggersManager::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const bool PrintLogToConsoleUnconditional, const bool PrintLogToFilesUnconditional, const bool PrintToConsole, const bool PrintLogToFiles, const std::int64_t SpecialLogFileIndex)
{
	try
	{
        ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str());

        if (LoggerMainObjectPointer)
        {
            if (LogThreadsToSeparatedFiles == true)
            {
                if (CurrentThreadId == LoggerMainObjectPointer->ThisThreadId)
                    LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsole, PrintLogToFiles, SpecialLogFileIndex);
                else
                {
                    {
                        lock_guard CreateNewLoggerForThreadLockGuardMutexObject{ CreateNewLoggerForThreadMutexObject };

                        auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId);

                        if (FoundLoggerIterator == LoggersThreadsObjectsPointersMap.end())
                            LoggersThreadsObjectsPointersMap[CurrentThreadId] = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("THREAD_" + to_string(LoggersThreadsObjectsPointersMap.size() + 1) + "_" + (stringstream() << CurrentThreadId).str()).c_str(), TaskName.c_str(), CurrentThreadId);
                    }

                    LoggersThreadsObjectsPointersMap[CurrentThreadId]->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsole, PrintLogToFiles, SpecialLogFileIndex);
                }
            }
            else
                LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsole, PrintLogToFiles, SpecialLogFileIndex);
        }
	}
	CATCH_COUT("logging message in loggers manager")
}

[[maybe_unused]] void LoggersManager::LogWarning(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogWarningsFileIndex);
}

[[maybe_unused]] void LoggersManager::LogError(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogErrorsFileIndex);
}

[[maybe_unused]] void LoggersManager::LogException(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogExceptionsFileIndex);
}

[[maybe_unused]] void LoggersManager::LogErrorAndException(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogErrorsAndExceptionsFileIndex);
}

[[maybe_unused]] void LoggersManager::LogCritical(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogCriticalFileIndex);
}

[[maybe_unused]] void LoggersManager::LogInformation(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogInformationFileIndex);
}

[[maybe_unused]] void LoggersManager::LogImportant(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogImportantFileIndex);
}

[[maybe_unused]] void LoggersManager::LogStatistics(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogStatisticsFileIndex);
}

[[maybe_unused]] void LoggersManager::LogDebug(const stringstream& Message)
{
    LogMessageBool(Message.str(), true, false, false, true, true, LogDebugFileIndex);
}

[[maybe_unused]] void LoggersManager::LogWarn(const stringstream& Message)
{
    LogWarning(Message);
}

[[maybe_unused]] void LoggersManager::LogErr(const stringstream& Message)
{
    LogError(Message);
}

[[maybe_unused]] void LoggersManager::LogExc(const stringstream& Message)
{
    LogException(Message);
}

[[maybe_unused]] void LoggersManager::LogErrAndExc(const stringstream& Message)
{
    LogErrorAndException(Message);
}

[[maybe_unused]] void LoggersManager::LogCrit(const stringstream& Message)
{
    LogCritical(Message);
}

[[maybe_unused]] void LoggersManager::LogInfo(const stringstream& Message)
{
    LogInformation(Message);
}

[[maybe_unused]] void LoggersManager::LogImp(const stringstream& Message)
{
    LogImportant(Message);
}

[[maybe_unused]] void LoggersManager::LogStat(const stringstream& Message)
{
    LogStatistics(Message);
}

[[maybe_unused]] void LoggersManager::LogDeb(const stringstream& Message)
{
    LogDebug(Message);
}