
#include "CellEngineVoxelSimulationSpaceStatistics.h"

void CellEngineVoxelSimulationSpaceStatistics::CountStatisticsOfVoxelSimulationSpace()
{
    try
    {
        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension; PosZ++)
                {
                    if (GetParticleFromUniqueId(GetSpaceVoxel(PosX, PosY, PosZ)).ParticleKindId != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
    CATCH("counting statistics of voxel simulation space")
}