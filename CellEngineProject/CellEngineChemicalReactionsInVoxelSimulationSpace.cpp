
#include <unordered_set>

#include "CellEngineExecutionTimeStatistics.h"
#include "CellEngineParticlesVoxelsOperations.h"
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"

void CellEngineChemicalReactionsInVoxelSimulationSpace::FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, const bool UpdateNucleotides, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
                        if (const UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosY, PosZ); ParticleIndex != GetZeroSimulationSpaceVoxel())
                            if (!FoundParticleIndexes.contains(ParticleIndex))
                                if (GetParticles().contains(ParticleIndex))
                                    SaveParticleFoundInProximity(ParticleIndex, FoundParticleIndexes, UpdateNucleotides);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles += chrono::duration(stop_time - start_time);
    }
    CATCH("finding particles in proximity of simulation space for selected local space")
};

void CellEngineChemicalReactionsInVoxelSimulationSpace::MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const RealType AddX, const RealType AddY, const RealType AddZ)
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
