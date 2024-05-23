
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateAllRealRandomParticles()
{
    try
    {
        ClearVoxelSpaceAndParticles();

        GenerateRealRandomOtherParticles();
        GenerateRealRandomBasicParticles();
        GenerateRealRandomLipidParticles();
        GenerateRealRandom_tRNAParticles();
        GenerateRealRandom_mRNAParticles();
        GenerateRealRandom_rRNAParticles();
        GenerateRealRandomMembraneParticles();
        GenerateRealRandomRibosomesParticles();
        GenerateRealRandomPolymeraseParticles();
        GenerateRealRandomRNAPolymeraseParticles();
    }
    CATCH("generating all random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomOtherParticles()
{
    try
    {
        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == true)
            {
                //DODAJ 10 CZASTEK
            }
    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomBasicParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomLipidParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandom_tRNAParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandom_mRNAParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandom_rRNAParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}


void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomMembraneParticles()
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
        UnsignedInt MembraneProteinCounter = 0;
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
            {
                for (const auto& ParticleKindSpecialDataObject:  ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticlesTypes::MembraneProtein)
                    {
                        MembraneProteinCounter++;
                    }
            }
        }
        LoggersManagerObject.Log(STREAM("Number Of Membrane Proteins = " << MembraneProteinCounter));


        UnsignedInt RibosomesProteinCounter = 0;
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
            {
                for (const auto& ParticleKindSpecialDataObject:  ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticlesTypes::RibosomesProtein)
                    {
                        RibosomesProteinCounter++;
                    }
            }
        }
        LoggersManagerObject.Log(STREAM("Number Of Ribosomes Proteins = " << RibosomesProteinCounter));

        UnsignedInt RNAPolymeraseProteinCounter = 0;
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
            {
                for (const auto& ParticleKindSpecialDataObject:  ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticlesTypes::RNAPolymeraseProtein)
                    {
                        RNAPolymeraseProteinCounter++;
                    }
            }
        }
        LoggersManagerObject.Log(STREAM("Number Of RNA Polymerase Proteins = " << RNAPolymeraseProteinCounter));

        UnsignedInt ProteinFracCounter = 0;
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
        {
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
            {
                for (const auto& ParticleKindSpecialDataObject:  ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticlesTypes::ProteinFrac)
                    {
                        ProteinFracCounter++;
                    }
            }
        }
        LoggersManagerObject.Log(STREAM("Number Of Proteins Frac = " << ProteinFracCounter));

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

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomRibosomesParticles()
{
    try
    {

    }
    CATCH("generating random ribosomes particles")
}


void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomPolymeraseParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateRealRandomRNAPolymeraseParticles()
{
    try
    {

    }
    CATCH("generating real random particles")
}

