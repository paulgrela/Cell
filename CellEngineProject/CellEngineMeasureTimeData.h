
#ifndef CELL_ENGINE_MEASURE_TIME_DATA_H
#define CELL_ENGINE_MEASURE_TIME_DATA_H

#include <chrono>

namespace CellEngineMeasureTimeData
{
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForDrawingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForPreparingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForTotalPreparingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCheckingPreparingParticles{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCheckingPreparingParticles2{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCheckingPreparingParticles3{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory0{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory1{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory2{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory21{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory22{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory23{ 0 };
    inline std::common_type_t<std::chrono::duration<long, std::ratio<1, 1000000000>>, std::chrono::duration<long, std::ratio<1, 1000000000>>> ExecutionDurationTimeForCopyingParticlesToGraphicMemory3{ 0 };
};

#endif
