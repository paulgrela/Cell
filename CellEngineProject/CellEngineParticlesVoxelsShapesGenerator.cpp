
#include "CellEngineParticlesVoxelsShapesGenerator.h"

using namespace std;

void CellEngineParticlesVoxelsShapesGenerator::SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(ListOfVoxelsType* FilledSpaceVoxels, const UniqueIdInt VoxelValue, const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ) const
{
    try
    {
        if (FilledSpaceVoxels != nullptr)
            FilledSpaceVoxels->emplace_back(PosX, PosY, PosZ);
        GetSpaceVoxel(PosX, PosY, PosZ) = VoxelValue;
    }
    CATCH("setting value to voxel")
}

bool CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace(const UnsignedInt PosXStart, const UnsignedInt PosYStart, const UnsignedInt PosZStart, const UnsignedInt StepX, const UnsignedInt StepY, const UnsignedInt StepZ, const UnsignedInt SizeOfParticleX, const UnsignedInt SizeOfParticleY, const UnsignedInt SizeOfParticleZ, const UniqueIdInt ValueToCheck)
{
    try
    {
        for (UnsignedInt PosX = PosXStart; PosX < PosXStart + SizeOfParticleX; PosX += StepX)
            for (UnsignedInt PosY = PosYStart; PosY < PosYStart + SizeOfParticleY; PosY += StepY)
                for (UnsignedInt PosZ = PosZStart; PosZ < PosZStart + SizeOfParticleZ; PosZ += StepZ)
                    if (GetSpaceVoxel(PosX, PosY, PosZ) != ValueToCheck)
                        return false;
    }
    CATCH("checking free space in cuboid selected space")

    return true;
}

void CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace(ListOfVoxelsType* FilledSpaceVoxels, const UniqueIdInt VoxelValue, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(FilledSpaceVoxels, VoxelValue, PosX, PosY, PosZ);
    }
    CATCH("setting value to voxels for cuboid selected space")
}

bool CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace(const UnsignedInt PosXStart, const UnsignedInt PosYStart, const UnsignedInt PosZStart, const UnsignedInt StepX, const UnsignedInt StepY, const UnsignedInt StepZ, const UnsignedInt RadiusXParam, const UnsignedInt RadiusYParam, const UnsignedInt RadiusZParam, const UniqueIdInt ValueToCheck)
{
    try
    {
        UnsignedInt RadiusParam = RadiusXParam;

        for (SignedInt x = 0; x < RadiusParam * 2; x += static_cast<SignedInt>(StepX))
            for (SignedInt y = 0; y < RadiusParam * 2; y += static_cast<SignedInt>(StepY))
                for (SignedInt z = 0; z < RadiusParam * 2; z += static_cast<SignedInt>(StepZ))
                {
                    SignedInt dx = static_cast<SignedInt>(RadiusParam) - x;
                    SignedInt dy = static_cast<SignedInt>(RadiusParam) - y;
                    SignedInt dz = static_cast<SignedInt>(RadiusParam) - z;
                    if ((dx * dx + dy * dy + dz * dz) <= (RadiusParam * RadiusParam))
                        if (GetSpaceVoxel(PosXStart + dx, PosYStart + dy, PosZStart + dz) != ValueToCheck)
                            return false;
                }
    }
    CATCH("checking free space in sphere selected space")

    return true;
};

void CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace(ListOfVoxelsType* FilledSpaceVoxels, UniqueIdInt VoxelValue, const UnsignedInt PosXStart, const UnsignedInt PosYStart, const  UnsignedInt PosZStart, const UnsignedInt StepX, const UnsignedInt StepY, const UnsignedInt StepZ, const UnsignedInt RadiusXParam, const UnsignedInt RadiusYParam, const UnsignedInt RadiusZParam)
{
    try
    {
        UnsignedInt RadiusParam = RadiusXParam;

        for (SignedInt x = 0; x < RadiusParam * 2; x += static_cast<SignedInt>(StepX))
            for (SignedInt y = 0; y < RadiusParam * 2; y += static_cast<SignedInt>(StepY))
                for (SignedInt z = 0; z < RadiusParam * 2; z += static_cast<SignedInt>(StepZ))
                {
                    SignedInt dx = static_cast<SignedInt>(RadiusParam) - x;
                    SignedInt dy = static_cast<SignedInt>(RadiusParam) - y;
                    SignedInt dz = static_cast<SignedInt>(RadiusParam) - z;
                    if ((dx * dx + dy * dy + dz * dz) <= (RadiusParam * RadiusParam))
                        SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(FilledSpaceVoxels, VoxelValue, PosXStart + dx, PosYStart + dy, PosZStart + dz);
                }
    }
    CATCH("setting value to voxels for sphere selected space")
};

bool CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace(const UnsignedInt PosXStart, const UnsignedInt PosYStart, const UnsignedInt PosZStart, const UnsignedInt StepX, const UnsignedInt StepY, const UnsignedInt StepZ, const UnsignedInt RadiusXParam, const UnsignedInt RadiusYParam, const UnsignedInt RadiusZParam, const UniqueIdInt ValueToCheck)
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
                        if (GetSpaceVoxel(PosXStart + dx, PosYStart + dy, PosZStart + dz) != ValueToCheck)
                            return false;
                }
    }
    CATCH("checking free space in ellipsoid selected space")

    return true;
};

void CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForEllipsoidSelectedSpace(ListOfVoxelsType* FilledSpaceVoxels, const UniqueIdInt VoxelValue, const UnsignedInt PosXStart, const UnsignedInt PosYStart, const UnsignedInt PosZStart, const UnsignedInt StepX, const UnsignedInt StepY, const UnsignedInt StepZ, const UnsignedInt RadiusXParam, const UnsignedInt RadiusYParam, const UnsignedInt RadiusZParam)
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
                        SetValueToSpaceVoxelWithFillingListOfVoxelsOfParticle(FilledSpaceVoxels, VoxelValue, PosXStart + dx, PosYStart + dy, PosZStart + dz);
                }
    }
    CATCH("setting value to voxels for ellipsoid selected space")
};

bool CellEngineParticlesVoxelsShapesGenerator::GenerateParticleVoxelsWhenSelectedSpaceIsFree(const UnsignedInt LocalNewParticleIndex, const UnsignedInt PosXStart, const UnsignedInt PosYStart, const UnsignedInt PosZStart, const UnsignedInt SizeOfParticleX, const UnsignedInt SizeOfParticleY, const UnsignedInt SizeOfParticleZ, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam, const CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, const SetValueToVoxelsForSelectedSpaceType SetValueToVoxelsForSelectedSpace)
{
    try
    {
        ListOfVoxelsType FilledVoxelsForRandomParticle;

        if ((this->*CheckFreeSpaceForSelectedSpace)(PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ, GetZeroSimulationSpaceVoxel()) == true)
        {
            if (PosXStart + SizeOfParticleX < StartXPosParam + SizeXParam && PosYStart + SizeOfParticleY < StartYPosParam + SizeYParam && PosZStart + SizeOfParticleZ < StartZPosParam + SizeZParam)
                (this->*SetValueToVoxelsForSelectedSpace)(&FilledVoxelsForRandomParticle, LocalNewParticleIndex, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

            if (FilledVoxelsForRandomParticle.empty() == false)
                GetParticleFromIndexForGenerator(LocalNewParticleIndex).ListOfVoxels = FilledVoxelsForRandomParticle;

            CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(GetParticleFromIndexForGenerator(LocalNewParticleIndex), &Particle::ListOfVoxels, &ParticleKind::ListOfVoxels, true);

            return true;
        }
    }
    CATCH_AND_THROW("generate particle in selected space")

    return false;
}
