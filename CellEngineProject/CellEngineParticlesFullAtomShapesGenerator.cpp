
#include "CellEngineParticlesFullAtomShapesGenerator.h"

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, const float PosXStart, const float PosYStart, const float PosZStart, const float StepX, const float StepY, const float StepZ, const float SizeOfParticleX, const float SizeOfParticleY, const float SizeOfParticleZ, const UniqueIdInt ValueToCheck)
{
    try
    {
        vector3_float32 Center{};
        const float Radius = (SizeOfParticleX + SizeOfParticleY + SizeOfParticleZ) / 3.0f;
        ListOfAtomsType ListOfAtoms;

        if (CellEngineConfigDataObject.CheckOnlyParticlesCenters == true)
            Center = { PosXStart + SizeOfParticleX / 2, PosYStart + SizeOfParticleY / 2, PosZStart  + SizeOfParticleZ / 2 };
        else
            SetValueToAtomsForCuboidSelectedSpace(&ListOfAtoms, PosXStart, PosYStart, PosZStart, StepX, StepY, StepZ, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, Center, ParticlesParam, GetSectorPos(Center.X, Center.Y, Center.Z), 0, 0, 0, SimulationSpaceSectorBounds{}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
    CATCH("checking free space in cuboid selected space")

    return false;
}

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float StepXParam, const float StepYParam, const float StepZParam, const float SizeXParam, const float SizeYParam, const float SizeZParam)
{
    try
    {
        for (float PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (float PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (float PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    FilledSpaceAtoms->emplace_back(PosX, PosY, PosZ);
    }
    CATCH("setting value to full atoms for cuboid selected space")
}

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, const float PosXStart, const float PosYStart, const float PosZStart, const float StepX, const float StepY, const float StepZ, const float RadiusXParam, const float RadiusYParam, const float RadiusZParam, const UniqueIdInt ValueToCheck)
{
    try
    {
        vector3_float32 Center{};
        const float Radius = (RadiusXParam + RadiusYParam + RadiusZParam) / 3.0f;
        ListOfAtomsType ListOfAtoms;

        if (CellEngineConfigDataObject.CheckOnlyParticlesCenters == true)
            Center = { PosXStart, PosYStart, PosZStart };
        else
            SetValueToAtomsForEllipsoidSelectedSpace(&ListOfAtoms, PosXStart, PosYStart, PosZStart, StepX, StepY, StepZ, RadiusXParam, RadiusYParam, RadiusZParam);

        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, Center, ParticlesParam, GetSectorPos(Center.X, Center.Y, Center.Z), 0, 0, 0, SimulationSpaceSectorBounds{}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
    CATCH("checking free space in ellipsoid selected space")

    return false;
};

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType* FilledSpaceAtoms, const float PosXStart, const float PosYStart, const float PosZStart, const float StepX, const float StepY, const float StepZ, const float RadiusXParam, const float RadiusYParam, const float RadiusZParam)
{
    try
    {
        for (float x = 0; x < RadiusXParam * 2; x += StepX)
            for (float y = 0; y < RadiusYParam * 2; y += StepY)
                for (float z = 0; z < RadiusZParam * 2; z += StepZ)
                {
                    const float dx = RadiusXParam - x;
                    const float dy = RadiusYParam - y;
                    const float dz = RadiusZParam - z;
                    if ((dx * dx * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam + dy * dy * RadiusXParam * RadiusXParam * RadiusZParam * RadiusZParam + dz * dz * RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam) <= (RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam))
                        FilledSpaceAtoms->emplace_back(PosXStart + dx, PosYStart + dy, PosZStart + dz);
                }
    }
    CATCH("setting value to voxels for ellipsoid selected space")
};

bool CellEngineParticlesFullAtomShapesGenerator::GenerateParticleAtomsWhenSelectedSpaceIsFree(const ParticlesContainer<Particle>& ParticlesParam, const UnsignedInt LocalNewParticleIndex, const float PosXStart, const float PosYStart, const float PosZStart, const float SizeOfParticleX, const float SizeOfParticleY, const float SizeOfParticleZ, const float StartXPosParam, const float StartYPosParam, const float StartZPosParam, const float SizeXParam, const float SizeYParam, const float SizeZParam, const CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, const SetValueToAtomsForSelectedSpaceType SetValueToAtomsForSelectedSpace)
{
    try
    {
        ListOfAtomsType FilledAtomsForRandomParticle;

        if ((this->*CheckFreeSpaceForSelectedSpace)(ParticlesParam, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ, 0) == true)
        {
            if (PosXStart + SizeOfParticleX < StartXPosParam + SizeXParam && PosYStart + SizeOfParticleY < StartYPosParam + SizeYParam && PosZStart + SizeOfParticleZ < StartZPosParam + SizeZParam)
                (this->*SetValueToAtomsForSelectedSpace)(&FilledAtomsForRandomParticle, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

            if (FilledAtomsForRandomParticle.empty() == false)
                GetParticleFromIndexForGenerator(LocalNewParticleIndex).ListOfAtoms = FilledAtomsForRandomParticle;

            CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<float, CellEngineAtom>(GetParticleFromIndexForGenerator(LocalNewParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, true);

            return true;
        }
    }
    CATCH_AND_THROW("generate particle atoms when selected space is free")

    return false;
}
