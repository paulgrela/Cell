
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

			for (const string& FileName : LoggersManagerObject.FilesNames)
				mkdir(string(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + FileName).c_str());
            #endif
            #ifdef UNIX_PLATFORM
            filesystem::create_directory(string(LogDirectory) + string("logs"));
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName);
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName);
            filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName);

            for (const string& FileName : LoggersManagerObject.FilesNames)
                filesystem::create_directory(string(LogDirectory) + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + FileName);
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
			Files.clear();
			Files.resize(LoggersManagerObject.FilesNames.size());
		}
	}
	CATCH_AND_THROW_COUT("aloccing resources for file in logger")
}

void Logger::OpenLogFiles()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
            for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.FilesNames.size(); FileNumber++)
                Files[FileNumber].open(string(LogDirectory) + OS_DIR_SEP + string("logs") + OS_DIR_SEP + MainDirectoryName + OS_DIR_SEP + TaskName + OS_DIR_SEP + LoggerName + OS_DIR_SEP + LoggersManagerObject.FilesNames[FileNumber] + OS_DIR_SEP + LoggersManagerObject.FilesNames[FileNumber] + string_utils::align_str(to_string(FileNumberInLog), '0', 5) + ".log.txt");
	}
	CATCH_AND_THROW_COUT("opening log files in logger")
}

void Logger::CloseLogFiles()
{
	try
	{
		if (LoggersManagerObject.PrintLogToFiles == true)
			for (auto& FileNumber : Files)
				FileNumber.close();
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

void Logger::WriteToLogsFromThread(const string& MessageStrToFile, const ThreadIdType CurrentThreadId)
{
	try
	{
		for (uint64_t FileNumber = 0; FileNumber < LoggersManagerObject.FilesNames.size(); FileNumber++)
			WriteToCommonLogFromThread(Files[FileNumber].is_open() == true && LoggersManagerObject.SelectiveWordsFunctions[FileNumber](MessageStrToFile), MessageStrToFile, Files[FileNumber], CurrentThreadId, FileNumber);
	}
	CATCH_AND_THROW_COUT("writing to logs from thread in logger")
}

void Logger::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const ThreadIdType CurrentThreadId, const bool PrintToConsole)
{
	try
	{
		string LocalMessageStr;
		uint64_t LineNumberInCommonLog;

		lock_guard<mutex> LockGuardObject{ LogMessageCoutMutexObject };

		LineNumberInCommonLog = LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog;
		LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog++;

		if (LoggersManagerObject.PrintLogToConsole == true && PrintToConsole == true)
		{
			LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToConsole, LoggersManagerObject.PrintLogDateTimeToConsole, LoggersManagerObject.PrintLogProcessIdToConsole, LoggersManagerObject.PrintLogProcessPriorityLevelToConsole, LoggersManagerObject.PrintLogThreadIdToConsole);
			WriteToCommonLogFromThread(true, LocalMessageStr, cout, CurrentThreadId, numeric_limits<uint64_t>::max());
		}

		if (LoggersManagerObject.PrintLogToFiles == true)
			if (LoggersManagerObject.LoggerMainObjectPointer)
			{
				if (LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId != CurrentThreadId)
				{
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					LoggersManagerObject.LoggerMainObjectPointer->WriteToLogsFromThread(LocalMessageStr, LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(CurrentThreadId);

					LineNumberInLog++;
					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId);
				}
				else
				{
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId);
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
			FilesNames.push_back(InitialFileName);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager files names")
}

void  LoggersManager::InitializeSelectiveWordsFunctions(const initializer_list<function<bool(const string&)>> InitialSelectiveWordsFunctions)
{
	try
	{
		for (const auto& InitialSelectiveWordsFunction : InitialSelectiveWordsFunctions)
			SelectiveWordsFunctions.push_back(InitialSelectiveWordsFunction);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager selective words functions")
}

void LoggersManager::InitializePrintingParameters(bool PrintLogToConsole, bool PrintLogToFiles, bool PrintLogLineNumberToConsole, bool PrintLogDateTimeToConsole, bool PrintLogProcessIdToConsole, bool PrintLogProcessPriorityLevelToConsole, bool PrintLogThreadIdToConsole, bool PrintLogLineNumberToFile, bool PrintLogDateTimeToFile, bool PrintLogProcessIdToFile, bool PrintLogProcessPriorityLevelToFile, bool PrintLogThreadIdToFile, uint64_t MaximalNumberOfLinesInOneFile)
{
	try
	{
		this->PrintLogToConsole = PrintLogToConsole;
		this->PrintLogToFiles = PrintLogToFiles;
	
		this->PrintLogLineNumberToConsole = PrintLogLineNumberToConsole;
		this->PrintLogDateTimeToConsole = PrintLogDateTimeToConsole;
		this->PrintLogProcessIdToConsole = PrintLogProcessIdToConsole;
		this->PrintLogProcessPriorityLevelToConsole = PrintLogProcessPriorityLevelToConsole;
		this->PrintLogThreadIdToConsole = PrintLogThreadIdToConsole;

		this->PrintLogLineNumberToFile = PrintLogLineNumberToFile;
		this->PrintLogDateTimeToFile = PrintLogDateTimeToFile;
		this->PrintLogProcessIdToFile = PrintLogProcessIdToFile;
		this->PrintLogProcessPriorityLevelToFile = PrintLogProcessPriorityLevelToFile;
		this->PrintLogThreadIdToFile = PrintLogThreadIdToFile;

		this->MaximalNumberOfLinesInOneFile = MaximalNumberOfLinesInOneFile;
	}
	CATCH_AND_THROW_COUT("initializing loggers manager printing conditions")
};

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

		if (LoggersManagerObject.FileNumberToIncreaseLineNumber >= LoggersManagerObject.FilesNames.size())
			throw runtime_error("FileNumberToIncreaseLineNumber is out of range of user defined files.");

		DrawMessageFunctionObject = DrawMessageFunctionObjectParameter;

		ThreadIdType CurrentThreadId = stoll(static_cast<stringstream&>(stringstream() << this_thread::get_id()).str());
		LoggerMainObjectPointer = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("LOGGER_COMMON").c_str(), TaskName.c_str(), CurrentThreadId);
		LoggerMainObjectPointer->LogMessageBool("START MAIN LOGGER_COMMON\n", true, CurrentThreadId, true);
	}
	CATCH_AND_THROW_COUT("initializing loggers manager data for task")
}

void LoggersManager::Log(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, true);
};

void LoggersManager::LogWithoutLineInfo(const stringstream& Message)
{
	LogMessageBool(Message.str(), false, true);
}

void LoggersManager::LogOnlyToFiles(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, false);
};

void LoggersManager::LogWithoutLineInfoOnlyToFiles(const stringstream& Message)
{
	LogMessageBool(Message.str(), false, false);
}

void LoggersManager::LogInColorTerminal(ostream& color(ostream& s), const stringstream& Message)
{
    this->LogOnlyToFiles(Message);
    cout << color << Message.str() << white << endl;
}

void LoggersManager::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const bool PrintToConsole)
{
	try
	{
		ThreadIdType CurrentThreadId = stoll(static_cast<stringstream&>(stringstream() << this_thread::get_id()).str());

        if (LoggerMainObjectPointer)
        {
            if (LogThreadsToSeparatedFiles == true)
            {
                if (CurrentThreadId == LoggerMainObjectPointer->ThisThreadId)
                    LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintToConsole);
                else
                {
                    {
                        lock_guard<mutex> CreateNewLoggerForThreadLockGuardMutexObject{ CreateNewLoggerForThreadMutexObject };

                        auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId);

                        if (FoundLoggerIterator == LoggersThreadsObjectsPointersMap.end())
                            LoggersThreadsObjectsPointersMap[CurrentThreadId] = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("THREAD_" + to_string(LoggersThreadsObjectsPointersMap.size() + 1) + "_" + static_cast<stringstream&>(stringstream() << CurrentThreadId).str()).c_str(), TaskName.c_str(), CurrentThreadId);
                    }

                    LoggersThreadsObjectsPointersMap[CurrentThreadId]->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintToConsole);
                }
            }
            else
                LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintToConsole);
        }
	}
	CATCH_COUT("logging message in loggers manager")
}