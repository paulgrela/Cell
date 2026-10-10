
#ifndef CELL_ENGINE_MEASURE_TIME_DATA_H
#define CELL_ENGINE_MEASURE_TIME_DATA_H

#include <chrono>

namespace CellEngineMeasureTimeData
{
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForPreparingParticlesForDrawingByCPUComputationsPhase1{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForPreparingParticlesForDrawingByCPUComputationsPhase2{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForComputingViewAndModelMatrixesByCPUForDrawingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForComputingViewAndModelMatrixesAndVisibilityByCPUForDrawingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesAndAtomsToGPUMemoryForComputations{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForComputationsOfParticlesInGPUInComputeShader{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForSwapingGraphicBuffersForDrawingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForDrawingParticlesWithSwapingGraphicBuffers{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForDrawingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForDrawingImGuiMenu{ 0 };
};

#endif
