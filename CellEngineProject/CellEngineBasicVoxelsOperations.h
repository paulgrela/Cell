
#ifndef CELL_ENGINE_BASIC_VOXELS_OPERATIONS_H
#define CELL_ENGINE_BASIC_VOXELS_OPERATIONS_H

#include "CellEngineTypes.h"

using SimulationSpaceVoxel = UniqueIdInt;

constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048 = 2048;

using Space_2048_2048_2048 = SimulationSpaceVoxel[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst2048];

class CellEngineBasicVoxelsOperations
{
protected:
    void* SpacePointer = nullptr;
protected:
    inline SimulationSpaceVoxel& GetSpaceVoxel(UnsignedInt x, UnsignedInt y, UnsignedInt z)
    {
        return (*static_cast<Space_2048_2048_2048*>(SpacePointer))[x][y][z];
    }
};

#endif
