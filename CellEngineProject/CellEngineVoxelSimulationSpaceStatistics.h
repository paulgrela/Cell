
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_STATISTICS_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_STATISTICS_H

#include "CellEngineBasicParticlesOperations.h"

class CellEngineVoxelSimulationSpaceStatistics : virtual public CellEngineBasicParticlesOperations, virtual public CellEngineParticlesVoxelsOperations
{
protected:
    UnsignedInt SumOfNotEmptyVoxels{};
public:
    void CountStatisticsOfVoxelSimulationSpace();
    [[nodiscard]] UnsignedInt GetSumOfNotEmptyVoxels() const
    {
        return SumOfNotEmptyVoxels;
    }
};

#endif
