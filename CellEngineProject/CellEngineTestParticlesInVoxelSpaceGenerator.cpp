
#include "CellEngineTestParticlesInVoxelSpaceGenerator.h"

using namespace std;

void CellEngineTestParticlesInVoxelSpaceGenerator::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << MaxParticleIndex << " " << (MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(0, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);
        uniform_int_distribution<ElectricChargeType> UniformDistributionObjectElectricChargeParticle_int64t(-5, 5);

        vector<UnsignedInt> LocalNewParticlesIndexes;

        for (UniqueIdInt ParticleNumber = 1; ParticleNumber <= NumberOfRandomParticles; ParticleNumber++)
            LocalNewParticlesIndexes.emplace_back(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), UniformDistributionObjectObjectOfParticle_Uint64t(mt64R), 0, -1, 1, UniformDistributionObjectElectricChargeParticle_int64t(mt64R), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));

        vector<vector3_16> FilledVoxelsForRandomParticle;

        for (const auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
        {
            UnsignedInt RandomSizeOfParticles = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);
            GenerateParticleVoxelsWhenSelectedSpaceIsFree(LocalNewParticleIndex, UniformDistributionObjectX_Uint64t(mt64R), UniformDistributionObjectY_Uint64t(mt64R), UniformDistributionObjectZ_Uint64t(mt64R), RandomSizeOfParticles, RandomSizeOfParticles, RandomSizeOfParticles, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        }
    }
    CATCH("generating random particles in selected space")
}

void CellEngineTestParticlesInVoxelSpaceGenerator::GeneratePlanedCuboidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 0, 1, -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 1, StartYPosParam + 3, StartZPosParam + 7, 1, 1, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 1, 1, -1, 1, 2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 3, StartYPosParam + 13, StartZPosParam + 6, 1, 1, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 2, 1, -1, 1, 3, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 6, StartYPosParam + 12, StartZPosParam + 5, 1, 1, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 3, 1, -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 9, StartYPosParam + 9, StartZPosParam + 9, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 4, 1, -1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 17, StartYPosParam + 15, StartZPosParam + 15, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 5, 1, -1, 1, -3, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 17, StartYPosParam + 19, StartZPosParam + 4, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, -2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 21, StartYPosParam + 20, StartZPosParam + 19, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, 2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 3, StartYPosParam + 20, StartZPosParam + 7, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, -2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 5, StartYPosParam + 15, StartZPosParam + 20, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 10, 1, -1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 21, StartYPosParam + 5, StartZPosParam + 3, 2, 2, 2, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);

        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('T'), -1, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 9, StartYPosParam + 12, StartZPosParam + 3, 2, 2, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('A'), -1, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 9, StartYPosParam + 16, StartZPosParam + 3, 2, 2, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('C'), -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 9, StartYPosParam + 19, StartZPosParam + 3, 2, 2, 1, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);

        GenerateOneStrand(CellEngineConfigDataObject.RNAIdentifier, "UCGAGAA", StartXPosParam + 5, StartYPosParam + 5, StartZPosParam + 3, 2, 1, 2, 2, 0, 0);
    }
    CATCH("generating planed cuboid particles in selected space")
}

void CellEngineTestParticlesInVoxelSpaceGenerator::GeneratePlanedEllipsoidParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt StepXParam, UnsignedInt StepYParam, UnsignedInt StepZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 2, 1, -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 9, StartYPosParam + 9, StartZPosParam + 10, 3, 4, 7, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForEllipsoidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 3, 1, -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 18, StartYPosParam + 18, StartZPosParam + 18, 3, 3, 3, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForEllipsoidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForEllipsoidSelectedSpace);
        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 4, 1, -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), StartXPosParam + 13, StartYPosParam + 13, StartZPosParam + 18, 3, 3, 3, StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
    }
    CATCH("generating planed ellipsoid particles in selected space")
}
