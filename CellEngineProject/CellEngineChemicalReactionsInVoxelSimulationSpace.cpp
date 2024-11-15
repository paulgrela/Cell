
#include <unordered_set>

#include "CellEngineExecutionTimeStatistics.h"
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"

void CellEngineChemicalReactionsInVoxelSimulationSpace::FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam, const CurrentThreadPosType& CurrentThreadPos)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX >= 0 && PosY >= 0 && PosZ >= 0 && PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        if (GetSpaceVoxel(PosX, PosY, PosZ) != GetZeroSimulationSpaceVoxel())
                        {
                            if (CurrentThreadIndex == 0)
                            {
                                SaveParticleFoundInProximity(GetSpaceVoxel(PosX, PosY, PosZ), FoundParticleIndexes, UpdateNucleotides, CurrentThreadPos);
                            }
                            else
                            {
                                if (ParticlesForThreads.contains(GetSpaceVoxel(PosX, PosY, PosZ)))
                                    SaveParticleFoundInProximity(GetSpaceVoxel(PosX, PosY, PosZ), FoundParticleIndexes, UpdateNucleotides, CurrentThreadPos);
                            }
                        }

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles += chrono::duration(stop_time - start_time);
    }
    CATCH("finding particles in proximity of simulation space for selected local space")
};

void CellEngineChemicalReactionsInVoxelSimulationSpace::MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const SignedInt AddX, const SignedInt AddY, const SignedInt AddZ)
{
    try
    {
        MoveParticleNearOtherParticleIfVoxelSpaceIsEmptyOrNearSpace(ParticleObject, NewPositionParticleObject, AddX, AddY, AddZ);
    }
    CATCH("moving particle near other particle if space is empty or to near space")
}

void CellEngineChemicalReactionsInVoxelSimulationSpace::ClearSpaceForParticle(Particle& ParticleObject, const bool ClearVoxels)
{
    try
    {
        if (ClearVoxels == true)
            SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
    }
    CATCH("clearing space for voxel")
}
