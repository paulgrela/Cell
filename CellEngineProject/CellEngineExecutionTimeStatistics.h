
#ifndef CELL_ENGINE_EXECUTION_TIME_STATISTICS_H
#define CELL_ENGINE_EXECUTION_TIME_STATISTICS_H

#include <chrono>

class CellEngineExecutionTimeStatistics
{
public:
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForFindingParticles { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForSavingFoundParticles { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForMakingChemicalReactions { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForMakingChemicalReactionsSpecialFunctions { 0 };
    std::common_type<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>>::type ExecutionDurationTimeForChoosingParticlesForMakingChemicalReactions { 0 };
};

inline CellEngineExecutionTimeStatistics CellEngineExecutionTimeStatisticsObject;

#endif
