
#include <cmath>

#include "DateTimeUtils.h"

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

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::PrintNumberOfParticlesForAllMainTypesOfParticles()
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
    }
    CATCH("printing number of particles for all main types of particles")
}

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineRealRandomParticlesInVoxelSpaceGenerator::GetRandomPositionInsideSphere(const UnsignedInt Radius, const UnsignedInt RadiusSize)
{
    try
    {
        UnsignedInt PosXStart = 512;
        UnsignedInt PosYStart = 512;
        UnsignedInt PosZStart = 512;

        uniform_int_distribution<SignedInt> UniformDistributionObjectMoveParticleDirection_int64t(0, 1024);

        UnsignedInt TryFindNewRandomPositionCounter = 0;

        while (TryFindNewRandomPositionCounter < MaxNumberOfTriesToFindARandomPositionInCell)
        {
            TryFindNewRandomPositionCounter++;

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

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::InsertNewRandomParticlesForType(ParticlesTypes ParticleTypeParam, const UnsignedInt Radius, const UnsignedInt RadiusSize)
{
    try
    {
        UnsignedInt NumberOfErrors = 0;

        CellEngineUseful::SwitchOffLogs();

        const auto start_time = chrono::high_resolution_clock::now();

        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.EntityId >= StartParticleKindId)
                if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                    if (ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType == ParticleTypeParam)
                        for (const auto& ParticleKindSpecialDataObject: ParticleKindObject.second.ParticleKindSpecialDataSector)
                            for (UnsignedInt ParticleCounter = 1; ParticleCounter <= ParticleKindSpecialDataObject.CounterAtStartOfSimulation; ParticleCounter++)
                            {
                                LoggersManagerObject.Log(STREAM("Particle Type = " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleTypeParam) << " Counter = " << ParticleCounter));

                                auto GeneIter = ParticlesKindsManagerObject.Genes.find(ParticleKindSpecialDataObject.GeneId);
                                if (GeneIter != ParticlesKindsManagerObject.Genes.end())
                                {
                                    LoggersManagerObject.Log(STREAM("Gene Number = " << GeneIter->second.NumId));

                                    auto GeneralSize = static_cast<UnsignedInt>(pow(GeneIter->second.Sequence.length(), (1.0 / 3.0)));
                                    UnsignedInt Size = GeneralSize > 10 ? GeneralSize - 3 : GeneralSize;

                                    UnsignedInt PosX, PosY, PosZ;
                                    bool TryResult = false;

                                    UnsignedInt TryInsertNewParticleCounter = 0;

                                    while (TryInsertNewParticleCounter < MaxNumberOfTriesToInsertNewParticle && TryResult == false)
                                    {
                                        tie(PosX, PosY, PosZ) = GetRandomPositionInsideSphere(Radius, RadiusSize);
                                        TryResult = GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ParticleKindObject.second.EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosX, PosY, PosZ, Size, Size, Size, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
                                        //If (TryResult == true) DODAJ Sequence z Genu
                                        TryInsertNewParticleCounter++;
                                    }

                                    if (TryInsertNewParticleCounter == MaxNumberOfTriesToInsertNewParticle)
                                    {
                                        NumberOfErrors++;
                                        LoggersManagerObject.Log(STREAM("ERROR Tried insert too many times - Gene Length = " << GeneIter->second.Sequence.length() << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size << endl));
                                    }
                                    else
                                        LoggersManagerObject.Log(STREAM("Try Insert = " << TryInsertNewParticleCounter << " Gene Length = " << GeneIter->second.Sequence.length() << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size << endl));
                                }
                            }

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, string("Execution of generating all " + ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleTypeParam) + " particles from file has taken time: ").c_str(), "executing printing duration_time")));

        LoggersManagerObject.Log(STREAM("NUMBER OF ERRORS = " << NumberOfErrors));
    }
    CATCH("generating random membrane particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateAllRealRandomParticles()
{
    try
    {
        ClearVoxelSpaceAndParticles();

        CellEngineUseful::SwitchOffLogs();

        const auto start_time = chrono::high_resolution_clock::now();

        InsertNewRandomParticlesForType(ParticlesTypes::MembraneProtein, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::RibosomesProtein, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::PolymeraseProtein, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymeraseProtein, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::OtherProtein, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::ProteinFrac, 400, 30);

        InsertNewRandomParticlesForType(ParticlesTypes::tRNA, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::mRNA, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::rRNA, 400, 30);

        InsertNewRandomParticlesForType(ParticlesTypes::Basic, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::Lipid, 400, 30);
        InsertNewRandomParticlesForType(ParticlesTypes::Other, 400, 30);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating all random particles from file has taken time: ", "executing printing duration_time")));

        PrintNumberOfParticlesForAllMainTypesOfParticles();
    }
    CATCH("generating all random particles")
}