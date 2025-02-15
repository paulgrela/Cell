
#ifndef CELL_ENGINE_EXECUTION_TIME_STATISTICS_H
#define CELL_ENGINE_EXECUTION_TIME_STATISTICS_H

#include <chrono>

#include "Logger.h"
#include "DateTimeUtils.h"

class CellEngineExecutionTimeStatistics
{
public:
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForFindingParticles { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForSavingFoundParticles { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForMakingChemicalReactions { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForMakingCancelledChemicalReactions { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions { 0 };
public:
    void ZeroMeasureTime()
    {
        ExecutionDurationTimeForFindingParticles = std::chrono::seconds::zero();
        ExecutionDurationTimeForSavingFoundParticles = std::chrono::seconds::zero();
        ExecutionDurationTimeForMakingChemicalReactions = std::chrono::seconds::zero();
        ExecutionDurationTimeForMakingCancelledChemicalReactions = std::chrono::seconds::zero();
        ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions = std::chrono::seconds::zero();
        ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions = std::chrono::seconds::zero();
    }
public:
    void PrintMeasureTime() const
    {
        LoggersManagerObject.Log(STREAM(""));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForFindingParticles, "Execution of finding particles has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForSavingFoundParticles, "Execution of saving found particles has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForMakingChemicalReactions, "Execution of making chemical reactions has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForMakingCancelledChemicalReactions, "Execution of making cancelled chemical reactions has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions, "Execution of making chemical reactions special functions has taken time: ","Execution in threads")));
        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLine(ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions, "Execution of saving choosing particles for making chemical reactions has taken time: ","Execution in threads")));
    }
};

inline CellEngineExecutionTimeStatistics CellEngineExecutionTimeStatisticsObject;

#endif
