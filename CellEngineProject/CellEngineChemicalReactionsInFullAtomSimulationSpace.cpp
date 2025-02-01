
#include "CellEngineExecutionTimeStatistics.h"
#include "CellEngineChemicalReactionsInFullAtomSimulationSpace.h"

void CellEngineChemicalReactionsInFullAtomSimulationSpace::FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::unordered_set<UnsignedInt>& FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        const auto start_time = chrono::high_resolution_clock::now();

        // for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
        //     for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
        //         for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
        //             if (PosX < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension)
        //             {
        //                 UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosY, PosZ);
        //                 if (ParticleIndex != GetZeroSimulationSpaceVoxel())
        //                     if (!FoundParticleIndexes.contains(ParticleIndex))
        //                         if (GetParticles().contains(ParticleIndex))
        //                             SaveParticleFoundInProximity(ParticleIndex, FoundParticleIndexes, UpdateNucleotides);
        //             }


        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineExecutionTimeStatisticsObject.ExecutionDurationTimeForFindingParticles += chrono::duration(stop_time - start_time);
    }
    CATCH("finding particles in proximity of simulation space for selected local space")
};

void CellEngineChemicalReactionsInFullAtomSimulationSpace::MoveParticleNearOtherParticleIfSpaceIsEmptyOrNearSpace(Particle &ParticleObject, const Particle &NewPositionParticleObject, const float AddX, const float AddY, const float AddZ)
{
    try
    {
        MoveParticleNearOtherParticleIfFullAtomSpaceIsEmptyOrNearSpace(ParticleObject, Particles, CurrentSectorPos, NewPositionParticleObject, AddX, AddY, AddZ);
    }
    CATCH("moving particle near other particle if space is empty or to near space")
}

void CellEngineChemicalReactionsInFullAtomSimulationSpace::ClearSpaceForParticle(Particle& ParticleObject, const bool ClearFullAtoms)
{
    // try
    // {
    //     if (ClearFullAtoms == true)
    //         SetAllAtomsInListOfAtomsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
    // }
    // CATCH("clearing space for voxel")
}