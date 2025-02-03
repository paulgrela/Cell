
#include "CellEngineParticlesFullAtomShapesGenerator.h"

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace(const float PosXStart, const float PosYStart, const float PosZStart, const float StepX, const float StepY, const float StepZ, const float SizeOfParticleX, const float SizeOfParticleY, const float SizeOfParticleZ, const UniqueIdInt ValueToCheck)
{
    try
    {
        for (float PosX = PosXStart; PosX < PosXStart + SizeOfParticleX; PosX += StepX)
            for (float PosY = PosYStart; PosY < PosYStart + SizeOfParticleY; PosY += StepY)
                for (float PosZ = PosZStart; PosZ < PosZStart + SizeOfParticleZ; PosZ += StepZ)
                    //if (GetSpaceFullAtom(PosX, PosY, PosZ) != ValueToCheck)
                        ;return false;
    }
    CATCH("checking free space in cuboid selected space")

    return true;
}

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, const UniqueIdInt FullAtomValue, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float StepXParam, const float StepYParam, const float StepZParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
{
    try
    {
        for (float PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (float PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (float PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    ;//SetValueToSpaceFullAtomWithFillingListOfFullAtomsOfParticle(FilledSpaceFullAtoms, FullAtomValue, PosX, PosY, PosZ);
    }
    CATCH("setting value to full atoms for cuboid selected space")
}

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace(const float PosXStart, const float PosYStart, const float PosZStart, const float StepX, const float StepY, const float StepZ, const float RadiusXParam, const float RadiusYParam, const float RadiusZParam, const UniqueIdInt ValueToCheck)
{
    try
    {
        for (SignedInt x = 0; x < RadiusXParam * 2; x += static_cast<SignedInt>(StepX))
            for (SignedInt y = 0; y < RadiusYParam * 2; y += static_cast<SignedInt>(StepY))
                for (SignedInt z = 0; z < RadiusZParam * 2; z+= static_cast<SignedInt>(StepZ))
                {
                    SignedInt dx = static_cast<SignedInt>(RadiusXParam) - x;
                    SignedInt dy = static_cast<SignedInt>(RadiusYParam) - y;
                    SignedInt dz = static_cast<SignedInt>(RadiusZParam) - z;
                    if ((dx * dx * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam + dy * dy * RadiusXParam * RadiusXParam * RadiusZParam * RadiusZParam + dz * dz * RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam) <= (RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam))
                        ;//if (GetSpaceVoxel(PosXStart + dx, PosYStart + dy, PosZStart + dz) != ValueToCheck)
                        //    return false;
                }
    }
    CATCH("checking free space in ellipsoid selected space")

    return true;
};

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, UniqueIdInt VoxelValue, float PosXStart, float PosYStart, float PosZStart, const float StepX, const float StepY, const float StepZ, float RadiusXParam, float RadiusYParam, float RadiusZParam)
{
    try
    {
        for (SignedInt x = 0; x < RadiusXParam * 2; x += static_cast<SignedInt>(StepX))
            for (SignedInt y = 0; y < RadiusYParam * 2; y += static_cast<SignedInt>(StepY))
                for (SignedInt z = 0; z < RadiusZParam * 2; z += static_cast<SignedInt>(StepZ))
                {
                    SignedInt dx = static_cast<SignedInt>(RadiusXParam) - x;
                    SignedInt dy = static_cast<SignedInt>(RadiusYParam) - y;
                    SignedInt dz = static_cast<SignedInt>(RadiusZParam) - z;
                    if ((dx * dx * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam + dy * dy * RadiusXParam * RadiusXParam * RadiusZParam * RadiusZParam + dz * dz * RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam) <= (RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam))
                        ;//SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(FilledSpaceVoxels, VoxelValue, PosXStart + dx, PosYStart + dy, PosZStart + dz);
                }
    }
    CATCH("setting value to voxels for ellipsoid selected space")
};

bool CellEngineParticlesFullAtomShapesGenerator::GenerateParticleAtomsWhenSelectedSpaceIsFree(float LocalNewParticleIndex, float PosXStart, float PosYStart, float PosZStart, float SizeOfParticleX, float SizeOfParticleY, float SizeOfParticleZ, float StartXPosParam, float StartYPosParam, float StartZPosParam, float SizeXParam, float SizeYParam, float SizeZParam, CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, SetValueToAtomsForSelectedSpaceType SetValueToAtomsForSelectedSpace)
{
    try
    {
        ListOfAtomsType FilledAtomsForRandomParticle;

        if ((this->*CheckFreeSpaceForSelectedSpace)(PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ, 0) == true)
        {
            if (PosXStart + SizeOfParticleX < StartXPosParam + SizeXParam && PosYStart + SizeOfParticleY < StartYPosParam + SizeYParam && PosZStart + SizeOfParticleZ < StartZPosParam + SizeZParam)
                (this->*SetValueToAtomsForSelectedSpace)(&FilledAtomsForRandomParticle, LocalNewParticleIndex, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

            if (FilledAtomsForRandomParticle.empty() == false)
                GetParticleFromIndexForGenerator(LocalNewParticleIndex).ListOfAtoms = FilledAtomsForRandomParticle;

            CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<float, CellEngineAtom>(GetParticleFromIndexForGenerator(LocalNewParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, true);

            return true;
        }
    }
    CATCH_AND_THROW("generate particle atoms when selected space is free")

    return false;
}
