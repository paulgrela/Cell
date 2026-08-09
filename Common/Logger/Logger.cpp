
#include <syncstream>

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

Logger::Logger(const char* LogDirectoryParameter, const char* MainDirectoryNameParameter, const char* LoggerNameParameter, const char* TaskNameParameter, const uint64_t ThisThreadIdParameter, const bool PrintLogToFilesUnconditionalParam)
{
	try
	{
		PrintLogToFilesUnconditional = PrintLogToFilesUnconditionalParam;
		//
		// if (PrintLogToFilesUnconditional == true)
		// {
		// 	cout << "XXX555 = " << PrintLogToFilesUnconditional << endl;
		// 	getchar();
		// }

		ThisThreadId = ThisThreadIdParameter;
        LogDirectory = LogDirectoryParameter;
		MainDirectoryName = MainDirectoryNameParameter;
		LoggerName = LoggerNameParameter;
		TaskName = TaskNameParameter;

		CreateDirectories();
		AllocResourcesForFiles();
		OpenLogFiles();

		PrintLogToFilesUnconditional = false;
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
		//if (LoggersManagerObject.PrintLogToFiles == true)
		if (LoggersManagerObject.PrintLogToFiles == true || PrintLogToFilesUnconditional == true)
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
		//if (LoggersManagerObject.PrintLogToFiles == true)
		if (LoggersManagerObject.PrintLogToFiles == true || PrintLogToFilesUnconditional == true)
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
		// if (PrintLogToFilesUnconditional == true)
		// cout << "XXX444 = " << PrintLogToFilesUnconditional << endl; //getchar();
		//if (LoggersManagerObject.PrintLogToFiles == true)
		// if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != LoggerMainObjectPointer->ThisThreadId)
		// 	if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator == LoggersThreadsObjectsPointersMap.end())
		//	if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != ThisThreadId)
		//		cout << "XXX6 = " << PrintLogToFilesUnconditional << " " << CurrentThreadId << " " << ThisThreadId << endl;
				//if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator == LoggersThreadsObjectsPointersMap.end())

		if (LoggersManagerObject.PrintLogToFiles == true || PrintLogToFilesUnconditional == true)
		{
			//cout << "XXX5 = " << PrintLogToFilesUnconditional << endl;

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
		if (LoggersManagerObject.PrintLogToFiles == true || PrintLogToFilesUnconditional == true)
		//if (LoggersManagerObject.PrintLogToFiles == true || PrintLogToFilesUnconditional == true)
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

string Logger::CreateLogString(const string& MessageStr, const bool LogLineInfo, const ThreadIdType CurrentThreadId, const uint64_t LineNumberInCommonLog, const bool PrintLogLineNumber, const bool PrintLogDateTime, const bool PrintLogProcessId, const bool PrintLogProcessPriorityLevel, const bool PrintLogThreadId) const
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
			//std::osyncstream(StreamObject) << MessageStr << flush;

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

void Logger::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const ThreadIdType CurrentThreadId, const bool PrintLogToConsoleUnconditional, const bool PrintLogToFilesUnconditionalNOT, const bool PrintLogToConsole, const bool PrintLogToFiles, const std::int64_t SpecialLogFileIndex, const bool PrintLogToFilesUnconditionalParam)
{
	try
	{
		string LocalMessageStr;

		lock_guard LockGuardObject{ LogMessageCoutMutexObject };

		const uint64_t LineNumberInCommonLog = LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog;
		LoggersManagerObject.LoggerMainObjectPointer->LineNumberInLog++;

		if ((LoggersManagerObject.PrintLogToConsole == true && PrintLogToConsole == true) || PrintLogToConsoleUnconditional == true)
		{
			//cout << "XXX3" << endl;
			LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToConsole, LoggersManagerObject.PrintLogDateTimeToConsole, LoggersManagerObject.PrintLogProcessIdToConsole, LoggersManagerObject.PrintLogProcessPriorityLevelToConsole, LoggersManagerObject.PrintLogThreadIdToConsole);
			WriteToCommonLogFromThread(true, LocalMessageStr, cout, CurrentThreadId, numeric_limits<uint64_t>::max());
		}

		//cout << "XXX4" << endl;
		//if ((LoggersManagerObject.PrintLogToFiles == true && PrintLogToFiles == true) || PrintLogToFilesUnconditional == true)
		if ((LoggersManagerObject.PrintLogToFiles == true && PrintLogToFiles == true) || PrintLogToFilesUnconditionalParam == true)
			if (LoggersManagerObject.LoggerMainObjectPointer)
			{
				//cout << "XXX5" << endl;
				if (LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId != CurrentThreadId)
				{
					//cout << "XXX6" << endl;
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = LoggersManagerObject.LoggerMainObjectPointer->CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					LoggersManagerObject.LoggerMainObjectPointer->WriteToLogsFromThread(LocalMessageStr, LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId, SpecialLogFileIndex);

					LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(CurrentThreadId);

					//cout << "XXX60" << endl;
					LineNumberInLog++;
					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId, SpecialLogFileIndex);
					//cout << "XXX61" << endl;
				}
				else
				{
					//cout << "XXX7" << endl;
					LoggersManagerObject.LoggerMainObjectPointer->LimitLogSizeByClosingOldLogFilesAndOpeningNewLogFilesAfterMaximalLimitOfLinesInOldFileIsExceeded(LoggersManagerObject.LoggerMainObjectPointer->ThisThreadId);

					LocalMessageStr = CreateLogString(MessageStr, LogLineInfo, CurrentThreadId, LineNumberInCommonLog, LoggersManagerObject.PrintLogLineNumberToFile, LoggersManagerObject.PrintLogDateTimeToFile, LoggersManagerObject.PrintLogProcessIdToFile, LoggersManagerObject.PrintLogProcessPriorityLevelToFile, LoggersManagerObject.PrintLogThreadIdToFile);
					WriteToLogsFromThread(LocalMessageStr, CurrentThreadId, SpecialLogFileIndex);
					//cout << "XXX71" << endl;
				}
			}
	}
	CATCH_COUT("logging message in logger")
}

void LoggersManager::InitializeFilesNames(const initializer_list<const string> InitialFilesNames)
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

void LoggersManager::InitializePrintingParameters(const bool PrintLogToConsoleParam, const bool PrintLogToFilesParam, const bool PrintLogLineNumberToConsoleParam, const bool PrintLogDateTimeToConsoleParam, const bool PrintLogProcessIdToConsoleParam, const bool PrintLogProcessPriorityLevelToConsoleParam, const bool PrintLogThreadIdToConsoleParam, const bool PrintLogLineNumberToFileParam, const bool PrintLogDateTimeToFileParam, const bool PrintLogProcessIdToFileParam, const bool PrintLogProcessPriorityLevelToFileParam, const bool PrintLogThreadIdToFileParam, const uint64_t MaximalNumberOfLinesInOneFileParam, const bool PrintLogToCommonFileWhenPrintLogToSpecialFileParam)
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

        this->PrintLogToCommonFileWhenPrintLogToSpecialFile = PrintLogToCommonFileWhenPrintLogToSpecialFileParam;
	}
	CATCH_AND_THROW_COUT("initializing loggers manager printing conditions")
};

void LoggersManager::InitializeSpecialLogFiles(const bool CreateLogWarningsFileParam, const bool CreateLogErrorsFileParam, const bool CreateLogExceptionsFileParam, const bool CreateLogErrorsAndExceptionsFileParam, const bool CreateLogCriticalFileParam, const bool CreateLogInformationFileParam, const bool CreateLogImportantFileParam, const bool CreateLogStatisticsFileParam, const bool CreateLogDebugFileParam)
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
	// if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != LoggerMainObjectPointer->ThisThreadId)
	// 	if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator != LoggersThreadsObjectsPointersMap.end())
	// 	{
	// 		LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = true;
	// 		LogMessageBool(Message.str(), true, false, true, false, false, -1);
	// 		LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = false;
	// 	}

	LogMessageBool(Message.str(), true, true, true, false, false, -1);
};

void LoggersManager::LogOnlyToConsoleUnconditional(const stringstream& Message)
{
	LogMessageBool(Message.str(), true, true, false, false, false, -1);
};

void LoggersManager::LogOnlyToFilesUnconditional(const stringstream& Message)
{
	//cout << "XXX0" << endl; //getchar();
	//if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != LoggerMainObjectPointer->ThisThreadId)
		//if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator != LoggersThreadsObjectsPointersMap.end())
		{
			//cout << "XXX1" << endl;
			//LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = true; //NIE TEN NUMBER WATKU BO TO WATEK PRZED POWOLANIEM OBIEKTU WIEC MUSI BYC CARSH
			LogMessageBool(Message.str(), true, false, true, false, false, -1, true);

			//if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != LoggerMainObjectPointer->ThisThreadId)
			//	if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator != LoggersThreadsObjectsPointersMap.end())
			//LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = false;

			//cout << "XXX2" << endl;
			//return;
		}

	//LoggersManagerObject.PrintLogToFiles = true;
	//
	//LogMessageBool(Message.str(), true, false, true, false, false, -1);

	//cout << "XXX0" << endl;

	// const auto RememberPrintLogToFiles = LoggersManagerObject.PrintLogToFiles;
	// LoggersManagerObject.PrintLogToFiles = true;
	//
	//LogMessageBool(Message.str(), true, false, true, false, false, -1);
	//
	// LoggersManagerObject.PrintLogToFiles = RememberPrintLogToFiles;
};

void LoggersManager::LogInColorTerminal(ostream& color(ostream& s), const stringstream& Message)
{
    this->LogOnlyToFiles(Message);
    cout << color << Message.str() << white << endl;
}

void LoggersManager::LogMessageBool(const string& MessageStr, const bool LogLineInfo, const bool PrintLogToConsoleUnconditional, const bool PrintLogToFilesUnconditional, const bool PrintToConsoleParam, const bool PrintLogToFilesParam, const std::int64_t SpecialLogFileIndex, const bool PrintLogToFilesUnconditionalParam)
{
	try
	{
        ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str());

        if (LoggerMainObjectPointer)
        {
        								// if (PrintLogToFilesUnconditionalParam == true)
        								// cout << "XXX3SSS-" << LoggerMainObjectPointer << " " << PrintLogToFilesUnconditionalParam << endl;

            if (LogThreadsToSeparatedFiles == true)
            {
            							// if (PrintLogToFilesUnconditionalParam == true)
            							// cout << "XXX4KKK-" << LogThreadsToSeparatedFiles << " " << PrintLogToFilesUnconditionalParam << endl;

                if (CurrentThreadId == LoggerMainObjectPointer->ThisThreadId)
                    LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsoleParam, PrintLogToFilesParam, SpecialLogFileIndex);
                else
                {
                	// if (PrintLogToFilesUnconditionalParam == true)
                	// cout << "XXX5RRR-" << LogThreadsToSeparatedFiles << " " << PrintLogToFilesUnconditionalParam << endl;

	                {
		                lock_guard CreateNewLoggerForThreadLockGuardMutexObject{ CreateNewLoggerForThreadMutexObject };

	                	// cout << "XXX2" << endl;
                		// if (PrintLogToFilesUnconditionalParam == true)
                		// cout << "XXX6RRR-" << LogThreadsToSeparatedFiles << " " << PrintLogToFilesUnconditionalParam << endl;

	                	if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator == LoggersThreadsObjectsPointersMap.end())
	                	{
	                		//cout << "XXX3" << endl; getchar();
	                		//if (const ThreadIdType CurrentThreadId = stoll((stringstream() << this_thread::get_id()).str()); CurrentThreadId != LoggerMainObjectPointer->ThisThreadId)
	                		//if (const auto FoundLoggerIterator = LoggersThreadsObjectsPointersMap.find(CurrentThreadId); FoundLoggerIterator != LoggersThreadsObjectsPointersMap.end())
	                			//LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = true; //NIE TEN NUMBER WATKU BO TO WATEK PRZED POWOLANIEM OBIEKTU WIEC MUSI BYC CRASH

	                																									//if (PrintLogToFilesUnconditionalParam == true)
	                																									//cout << "XXX3DDD = " << PrintLogToFilesUnconditionalParam << endl; //getchar();

	                		//else faktycznie nie otwiera plików - problem logiki - co z otwieraniem plików - kiedy to zrobić - czy w momencie pierwszego wejscia i tworzenia obiektu czy
	                		//gdy mam pisac ale pliki nie otwarte - ale wtedy już jest zrobione
	                		//CZY ZAWSZE OTWIERAC PLIKI DLA WATKOW - chyba tak - ale wtedy za wczesniej pisze dane do plikow

	                		//LoggersThreadsObjectsPointersMap[CurrentThreadId] = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("THREAD_" + to_string(LoggersThreadsObjectsPointersMap.size() + 1) + "_" + (stringstream() << CurrentThreadId).str()).c_str(), TaskName.c_str(), CurrentThreadId, PrintLogToFilesUnconditionalParam);
	                		LoggersThreadsObjectsPointersMap[CurrentThreadId] = make_unique<Logger>(LogDirectory.c_str(), ActualDateTimeStr.c_str(), string("THREAD_" + to_string(LoggersThreadsObjectsPointersMap.size() + 1) + "_" + (stringstream() << CurrentThreadId).str()).c_str(), TaskName.c_str(), CurrentThreadId, true);
	                	}
                    }

                    //LoggersThreadsObjectsPointersMap[CurrentThreadId]->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsoleParam, PrintLogToFilesParam, SpecialLogFileIndex);
                	LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = PrintLogToFilesUnconditionalParam;
                    LoggersThreadsObjectsPointersMap[CurrentThreadId]->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsoleParam, PrintLogToFilesParam, SpecialLogFileIndex, PrintLogToFilesUnconditionalParam);
                	LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = false;
					//LoggersThreadsObjectsPointersMap[CurrentThreadId]->PrintLogToFilesUnconditional = false;
                }
            }
            else
                LoggerMainObjectPointer->LogMessageBool(MessageStr, LogLineInfo, CurrentThreadId, PrintLogToConsoleUnconditional, PrintLogToFilesUnconditional, PrintToConsoleParam, PrintLogToFilesParam, SpecialLogFileIndex);
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