
#include "CellEngineChemicalReactionsInVoxelSimulationSpace.h"

void CellEngineChemicalReactionsInVoxelSimulationSpace::FindParticlesInProximityInSimulationSpaceForSelectedLocalSpace(std::set<UnsignedInt> &FoundParticleIndexes, bool UpdateNucleotides, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX >= 0 && PosY >= 0 && PosZ >= 0 && PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        if (GetSpaceVoxel(PosX, PosY, PosZ) != GetZeroSimulationSpaceVoxel())
                            SaveParticleFoundInProximity(GetSpaceVoxel(PosX, PosY, PosZ), FoundParticleIndexes, UpdateNucleotides);
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
