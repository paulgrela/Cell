
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"

using namespace std;

tuple<UnsignedInt, UnsignedInt> CellEngineRealRandomParticlesInVoxelSpaceGenerator::GetNumberOfParticlesKind(ParticlesTypes ParticleTypeParam)
{
    UnsignedInt ParticlesCounter = 0;
    UnsignedInt ParticlesKindsCounter = 0;

    try
    {
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                for (const auto& ParticleKindSpecialDataObject: ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticleTypeParam)
                    {
                        TotalNumberOfAllParticles += ParticleKindSpecialDataObject.CounterAtStartOfSimulation;
                        ParticlesCounter += ParticleKindSpecialDataObject.CounterAtStartOfSimulation;
                        ParticlesKindsCounter++;
                    }
    }
    CATCH("getting number of particles of type")

    return { ParticlesKindsCounter, ParticlesCounter };
}

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
    try
    {
        LoggersManagerObject.Log(STREAM("Number Of Other Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::OtherProtein)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::OtherProtein))));
        LoggersManagerObject.Log(STREAM("Number Of Membrane Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::MembraneProtein)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::MembraneProtein))));
        LoggersManagerObject.Log(STREAM("Number Of Ribosomes Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::RibosomesProtein)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::RibosomesProtein))));
        LoggersManagerObject.Log(STREAM("Number Of RNA Polymerase Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymeraseProtein)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymeraseProtein))));
        LoggersManagerObject.Log(STREAM("Number Of RNA Polymerase Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::PolymeraseProtein)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::PolymeraseProtein))));
        LoggersManagerObject.Log(STREAM("Number Of Proteins Frac = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::ProteinFrac)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::ProteinFrac))));

        LoggersManagerObject.Log(STREAM("Number Of TRNA = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::tRNA)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::tRNA))));
        LoggersManagerObject.Log(STREAM("Number Of MRNA = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::mRNA)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::mRNA))));
        LoggersManagerObject.Log(STREAM("Number Of RRNA = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::rRNA)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::rRNA))));

        LoggersManagerObject.Log(STREAM("Number Of Basic = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::Basic)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::Basic))));
        LoggersManagerObject.Log(STREAM("Total Number Of All Particles = " << TotalNumberOfAllParticles));

        UnsignedInt Radius = 400;
        UnsignedInt RadiusSize = 30;
        UnsignedInt PosXStart = 512;
        UnsignedInt PosYStart = 512;
        UnsignedInt PosZStart = 512;
        UnsignedInt StepX = 10;
        UnsignedInt StepY = 10;
        UnsignedInt StepZ = 10;

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

