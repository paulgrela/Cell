
#include "CellEngineParticlesFullAtomShapesGenerator.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, const RealType PosXStart, const RealType PosYStart, const RealType PosZStart, const RealType StepX, const RealType StepY, const RealType StepZ, const RealType SizeOfParticleX, const RealType SizeOfParticleY, const RealType SizeOfParticleZ, const UniqueIdInt ValueToCheck)
{
    try
    {
        vector3_Real32 Center{};
        const RealType Radius = (SizeOfParticleX + SizeOfParticleY + SizeOfParticleZ) / 3.0f;
        ListOfAtomsType ListOfAtoms;

        if (CellEngineConfigDataObject.CheckOnlyParticlesCenters == true)
            Center = { PosXStart + SizeOfParticleX / 2, PosYStart + SizeOfParticleY / 2, PosZStart  + SizeOfParticleZ / 2 };
        else
            SetValueToAtomsForCuboidSelectedSpace(ListOfAtoms, 0, PosXStart, PosYStart, PosZStart, StepX, StepY, StepZ, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, Center, ParticlesParam, CellEngineUseful::GetSectorPos(Center.X, Center.Y, Center.Z), 0, 0, 0, SimulationSpaceSectorBounds{}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
    CATCH("checking free space in cuboid selected space")

    return false;
}

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForCuboidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, const EntityIdInt EntityId, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType StepXParam, const RealType StepYParam, const RealType StepZParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam)
{
    try
    {
        for (RealType PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += 1.5f * StepXParam)
            for (RealType PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += 1.5f * StepYParam)
                for (RealType PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += 1.5f * StepZParam)
                    FilledSpaceAtoms.emplace_back(EntityId, PosX, PosY, PosZ, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
    }
    CATCH("setting value to full atoms for cuboid selected space")
}

bool CellEngineParticlesFullAtomShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace(const ParticlesContainer<Particle>& ParticlesParam, const RealType PosXStart, const RealType PosYStart, const RealType PosZStart, const RealType StepX, const RealType StepY, const RealType StepZ, const RealType RadiusXParam, const RealType RadiusYParam, const RealType RadiusZParam, const UniqueIdInt ValueToCheck)
{
    try
    {
        vector3_Real32 Center{};
        const RealType Radius = (RadiusXParam + RadiusYParam + RadiusZParam) / 3.0f;
        ListOfAtomsType ListOfAtoms;

        if (CellEngineConfigDataObject.CheckOnlyParticlesCenters == true)
            Center = { PosXStart, PosYStart, PosZStart };
        else
            SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtoms, 0, PosXStart, PosYStart, PosZStart, StepX, StepY, StepZ, RadiusXParam, RadiusYParam, RadiusZParam);

        return CheckFreeSpaceAndBoundsForParticleMovedByVector(ListOfAtoms, Radius, 0, Center, ParticlesParam, CellEngineUseful::GetSectorPos(Center.X, Center.Y, Center.Z), 0, 0, 0, SimulationSpaceSectorBounds{}, CellEngineConfigDataObject.CheckOnlyParticlesCenters, false, false, false);
    }
    CATCH("checking free space in ellipsoid selected space")

    return false;
};

void CellEngineParticlesFullAtomShapesGenerator::SetValueToAtomsForEllipsoidSelectedSpace(ListOfAtomsType& FilledSpaceAtoms, const EntityIdInt EntityId, const RealType PosXStart, const RealType PosYStart, const RealType PosZStart, const RealType StepX, const RealType StepY, const RealType StepZ, const RealType RadiusXParam, const RealType RadiusYParam, const RealType RadiusZParam)
{
    try
    {
        for (RealType x = 0; x < RadiusXParam * 2; x += 1.5f * StepX)
            for (RealType y = 0; y < RadiusYParam * 2; y += 1.5f * StepY)
                for (RealType z = 0; z < RadiusZParam * 2; z += 1.5f * StepZ)
                {
                    const RealType dx = RadiusXParam - x;
                    const RealType dy = RadiusYParam - y;
                    const RealType dz = RadiusZParam - z;

                    if ((dx * dx * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam + dy * dy * RadiusXParam * RadiusXParam * RadiusZParam * RadiusZParam + dz * dz * RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam) <= (RadiusXParam * RadiusXParam * RadiusYParam * RadiusYParam * RadiusZParam * RadiusZParam))
                        FilledSpaceAtoms.emplace_back(EntityId, PosXStart + dx, PosYStart + dy, PosZStart + dz, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                }
    }
    CATCH("setting value to voxels for ellipsoid selected space")
};

bool CellEngineParticlesFullAtomShapesGenerator::GenerateParticleAtomsWhenSelectedSpaceIsFree(const ParticlesContainer<Particle>& ParticlesParam, const UnsignedInt LocalNewParticleIndex, const RealType PosXStart, const RealType PosYStart, const RealType PosZStart, const RealType SizeOfParticleX, const RealType SizeOfParticleY, const RealType SizeOfParticleZ, const RealType StartXPosParam, const RealType StartYPosParam, const RealType StartZPosParam, const RealType SizeXParam, const RealType SizeYParam, const RealType SizeZParam, const CheckFreeSpaceForSelectedSpaceType CheckFreeSpaceForSelectedSpace, const SetValueToAtomsForSelectedSpaceType SetValueToAtomsForSelectedSpace)
{
    try
    {
        ListOfAtomsType FilledAtomsForRandomParticle;

        if ((this->*CheckFreeSpaceForSelectedSpace)(ParticlesParam, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ, 0) == true)
        {
            (this->*SetValueToAtomsForSelectedSpace)(GetParticleFromIndexForGenerator(LocalNewParticleIndex).ListOfAtoms, GetParticleFromIndexForGenerator(LocalNewParticleIndex).EntityId, PosXStart, PosYStart, PosZStart, 1, 1, 1, SizeOfParticleX, SizeOfParticleY, SizeOfParticleZ);

            CellEngineBasicParticlesOperations::GetMinMaxCoordinatesForParticle<RealType, CellEngineAtom>(GetParticleFromIndexForGenerator(LocalNewParticleIndex), &Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, true);

            return true;
        }
    }
    CATCH_AND_THROW("generate particle atoms when selected space is free")

    return false;
}
