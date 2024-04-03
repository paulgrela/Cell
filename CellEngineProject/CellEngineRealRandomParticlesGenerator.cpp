
#include "CellEngineRealRandomParticlesGenerator.h"

void CellEngineRealRandomParticlesGenerator::GenerateAllRealRandomParticles()
{
    try
    {
        ClearVoxelSpaceAndParticles();
        GenerateRealRandomMembraneParticles();
        GenerateRealRandomRibosomesParticles();
    }
    CATCH("generating all random particles")
}

void CellEngineRealRandomParticlesGenerator::GenerateRealRandomMembraneParticles()
{
    UnsignedInt Radius = 400;
    UnsignedInt RadiusSize = 30;
    UnsignedInt PosXStart = 512;
    UnsignedInt PosYStart = 512;
    UnsignedInt PosZStart = 512;
    UnsignedInt StepX = 10;
    UnsignedInt StepY = 10;
    UnsignedInt StepZ = 10;

    try
    {
        for (SignedInt x = 0; x < Radius * 2; x += static_cast<SignedInt>(StepX))
            for (SignedInt y = 0; y < Radius * 2; y += static_cast<SignedInt>(StepY))
                for (SignedInt z = 0; z < Radius * 2; z += static_cast<SignedInt>(StepZ))
                {
                    SignedInt dx = static_cast<SignedInt>(Radius) - x;
                    SignedInt dy = static_cast<SignedInt>(Radius) - y;
                    SignedInt dz = static_cast<SignedInt>(Radius) - z;
                    if ((dx * dx + dy * dy + dz * dz >= (Radius - RadiusSize) * (Radius - RadiusSize)) && (dx * dx + dy * dy + dz * dz) <= (Radius * Radius))
                    {
                        GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 4, 1, -1, 1, -1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosXStart + dx, PosYStart + dy, PosZStart + dz, 3, 3, 3, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
                    }
                }
    }
    CATCH("generating random membrane particles")
}

void CellEngineRealRandomParticlesGenerator::GenerateRealRandomRibosomesParticles()
{
    try
    {

    }
    CATCH("generating random ribosomes particles")
}
