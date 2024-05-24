
#include <cmath>

#include "CellEngineConstants.h"

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

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineRealRandomParticlesInVoxelSpaceGenerator::GetRandomPositionInsideSphere()
{
    try
    {
        UnsignedInt Radius = 400;
        UnsignedInt RadiusSize = 30;
        UnsignedInt PosXStart = 512;
        UnsignedInt PosYStart = 512;
        UnsignedInt PosZStart = 512;

        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(0, 1024);

        UnsignedInt TryFindNewPositionCounter = 0;
        while (TryFindNewPositionCounter < 1000)
        {
            TryFindNewPositionCounter++;

            SignedInt dx = static_cast<SignedInt>(Radius) - UniformDistributionObjectMoveParticleDirection_int64t(mt64R);
            SignedInt dy = static_cast<SignedInt>(Radius) - UniformDistributionObjectMoveParticleDirection_int64t(mt64R);
            SignedInt dz = static_cast<SignedInt>(Radius) - UniformDistributionObjectMoveParticleDirection_int64t(mt64R);

            if ((dx * dx + dy * dy + dz * dz >= (Radius - RadiusSize) * (Radius - RadiusSize)) && (dx * dx + dy * dy + dz * dz) <= (Radius * Radius))
                return { PosXStart + dx, PosYStart + dy, PosZStart + dz };
        }
    }
    CATCH("getting random position inside cell")

    return { 0, 0, 0 };
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

        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.EntityId >= StartParticleKindId)
                if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                    if (ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType == ParticlesTypes::MembraneProtein)
                        for (const auto& ParticleKindSpecialDataObject: ParticleKindObject.second.ParticleKindSpecialDataSector)
                            for (UnsignedInt ParticleCounter = 1; ParticleCounter <= ParticleKindSpecialDataObject.CounterAtStartOfSimulation; ParticleCounter++)
                            {
                                LoggersManagerObject.Log(STREAM("Particle Counter = " << ParticleCounter));

                                auto GeneIter = ParticlesKindsManagerObject.Genes.find(ParticleKindSpecialDataObject.GeneId);
                                if (GeneIter != ParticlesKindsManagerObject.Genes.end())
                                {
                                    LoggersManagerObject.Log(STREAM("Gene Number = " << GeneIter->second.NumId));

                                    auto Size = static_cast<UnsignedInt>(pow(GeneIter->second.Sequence.length(), (1.0 / 3.0)));
                                    UnsignedInt PosX, PosY, PosZ;
                                    bool TryResult = false;

                                    UnsignedInt TryInsertNewParticleCounter = 0;

                                    while (TryInsertNewParticleCounter < 100 && TryResult == false)
                                    {
                                        //auto[ PosX, PosY, PosZ] = GetRandomPositionInsideSphere();
                                        tie(PosX, PosY, PosZ) = GetRandomPositionInsideSphere();
                                        TryResult = GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ParticleKindObject.second.EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosX, PosY, PosZ, Size, Size, Size, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
                                        //If (TryResult == true) DODAJ Sequence z Genu
                                        TryInsertNewParticleCounter++;

                                        //LoggersManagerObject.Log(STREAM("Try Insert = " << TryInsertNewParticleCounter << " Gene Length = " << GeneIter->second.Sequence.length() << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size));
                                    }
                                    LoggersManagerObject.Log(STREAM("Try Insert = " << TryInsertNewParticleCounter << " Gene Length = " << GeneIter->second.Sequence.length() << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size));
                                    LoggersManagerObject.Log(STREAM(""));
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

