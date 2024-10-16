
#include <cmath>

#include "DateTimeUtils.h"

#include "CellEngineConstants.h"
#include "CellEngineDataFile.h"
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"

using namespace std;

tuple<UnsignedInt, UnsignedInt> CellEngineRealRandomParticlesInVoxelSpaceGenerator::GetNumberOfParticlesKind(ParticlesTypes ParticleTypeParam, const bool AddToTotalNumberOfAllParticles)
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
                        if (AddToTotalNumberOfAllParticles == true)
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
        LoggersManagerObject.Log(STREAM("Number Of Ribosomes = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::Ribosome, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::Ribosome, false))));
        LoggersManagerObject.Log(STREAM("Number Of RNAPolymerases = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymerase, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymerase, false))));
        LoggersManagerObject.Log(STREAM("Number Of DNAPolymerases = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::DNAPolymerase, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::DNAPolymerase, false))));

        LoggersManagerObject.Log(STREAM("Number Of Membrane Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::MembraneProtein, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::MembraneProtein, false))));
        LoggersManagerObject.Log(STREAM("Number Of Ribosomes Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::RibosomeProtein, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::RibosomeProtein, false))));
        LoggersManagerObject.Log(STREAM("Number Of RNA Polymerase Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymeraseProtein, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::RNAPolymeraseProtein, false))));
        LoggersManagerObject.Log(STREAM("Number Of Polymerase Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::PolymeraseProtein, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::PolymeraseProtein, false))));
        LoggersManagerObject.Log(STREAM("Number Of Proteins Frac = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::ProteinFrac, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::ProteinFrac, false))));
        LoggersManagerObject.Log(STREAM("Number Of Other Proteins = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::OtherProtein, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::OtherProtein, false))));

        LoggersManagerObject.Log(STREAM("Number Of TRNA_uncharged = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::tRNA_uncharged, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::tRNA_uncharged, false))));
        LoggersManagerObject.Log(STREAM("Number Of TRNA_charged = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::tRNA_charged, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::tRNA_charged, false))));
        LoggersManagerObject.Log(STREAM("Number Of MRNA = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::mRNA, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::mRNA, false))));
        LoggersManagerObject.Log(STREAM("Number Of RRNA = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::rRNA, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::rRNA, false))));

        LoggersManagerObject.Log(STREAM("Number Of Basic = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::Basic, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::Basic, false))));
        LoggersManagerObject.Log(STREAM("Number Of Lipids = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::Lipid, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::Lipid, false))));
        LoggersManagerObject.Log(STREAM("Number Of Other = " << get<0>(GetNumberOfParticlesKind(ParticlesTypes::Other, true)) << " Total Number = " << get<1>(GetNumberOfParticlesKind(ParticlesTypes::Other, false))));

        LoggersManagerObject.Log(STREAM("Total Number Of All Particles Kinds = " << TotalNumberOfAllParticles));
        LoggersManagerObject.Log(STREAM("Total Number Of All Particles = " << Particles.size()));

        PrintInformationAboutRibosomesProteins();
    }
    CATCH("printing number of particles for all main types of particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::PrintInformationAboutRibosomesProteins()
{
    try
    {
        double Average30SLength = 0;
        for (const auto& GeneId : ParticlesKindsManagerObject.Ribosomes30SProteinsList)
            Average30SLength += (double) ParticlesKindsManagerObject.Genes[GeneId].Sequence.length();
        Average30SLength /= (double)ParticlesKindsManagerObject.Ribosomes30SProteinsList.size();

        double Average50SLength = 0;
        for (const auto& GeneId : ParticlesKindsManagerObject.Ribosomes50SProteinsList)
            Average50SLength += (double) ParticlesKindsManagerObject.Genes[GeneId].Sequence.length();
        Average50SLength /= (double)ParticlesKindsManagerObject.Ribosomes50SProteinsList.size();

        LoggersManagerObject.Log(STREAM("Number of 30S proteins = " << ParticlesKindsManagerObject.Ribosomes30SProteinsList.size() << " AVERAGE SIZE = " << (UnsignedInt)Average30SLength << " " << static_cast<UnsignedInt>(pow(Average30SLength, (1.0 / 3.0)))));
        LoggersManagerObject.Log(STREAM("Number of 50S proteins = " << ParticlesKindsManagerObject.Ribosomes50SProteinsList.size() << " AVERAGE SIZE = " << (UnsignedInt)Average50SLength << " " << static_cast<UnsignedInt>(pow(Average50SLength, (1.0 / 3.0)))));
    }
    CATCH("printing information about ribosomes proteins")
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

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::TryToGenerateRandomParticlesForType(const pair<const EntityIdInt, ParticleKind>& ParticleKindObject, const UnsignedInt Radius, const UnsignedInt RadiusSize, UnsignedInt& NumberOfErrors, const GeneIdInt GeneNumId, const string& GeneSequence, const UnsignedInt GeneSequenceLength)
{
    try
    {
        LoggersManagerObject.Log(STREAM("Gene Number = " << GeneNumId));

        auto GeneralSize = static_cast<UnsignedInt>(pow(GeneSequenceLength, (1.0 / 3.0)));
        UnsignedInt Size = GeneralSize > 10 ? GeneralSize - 3 : GeneralSize;

        UnsignedInt PosX, PosY, PosZ;
        bool TryResult = false;

        UnsignedInt TryInsertNewParticleCounter = 0;

        AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ParticleKindObject.second.EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, GeneSequence, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        while (TryInsertNewParticleCounter < MaxNumberOfTriesToInsertNewParticle && TryResult == false)
        {
            tie(PosX, PosY, PosZ) = GetRandomPositionInsideSphere(Radius, RadiusSize);
            TryResult = GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ParticleKindObject.second.EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, GeneSequence, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosX, PosY, PosZ, Size, Size, Size, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
            if (TryResult == false)
                RemoveParticle(MaxParticleIndex, true);
            TryInsertNewParticleCounter++;
        }

        if (TryInsertNewParticleCounter == MaxNumberOfTriesToInsertNewParticle)
        {
            NumberOfErrors++;
            LoggersManagerObject.Log(STREAM("ERROR Tried insert too many times - Gene Length = " << GeneSequenceLength << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size << endl));
        }
        else
            LoggersManagerObject.Log(STREAM("Try Insert = " << TryInsertNewParticleCounter << " Gene Length = " << GeneSequenceLength << " PosX = " << PosX << " PosY = " << PosY << " PosZ = " << PosZ << " Size " << Size << endl));
    }
    CATCH("trying to generate random particle for type")
}

UnsignedInt GetSizeOfGeneratedParticle(const ParticlesTypes ParticlesTypesObject)
{
    switch (ParticlesTypesObject)
    {
        case ParticlesTypes::Ribosome : return static_cast<UnsignedInt>(pow(40, 3));
        case ParticlesTypes::RNAPolymerase : return static_cast<UnsignedInt>(pow(13, 3));
        case ParticlesTypes::DNAPolymerase : return static_cast<UnsignedInt>(pow(15, 3));
        default : return 8;
    }
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
                    for (const auto& ParticleKindSpecialDataObject: ParticleKindObject.second.ParticleKindSpecialDataSector)
                        if (ParticleKindSpecialDataObject.ParticleType == ParticleTypeParam)
                            for (UnsignedInt ParticleCounter = 1; ParticleCounter <= ParticleKindSpecialDataObject.CounterAtStartOfSimulation; ParticleCounter++)
                            {
                                LoggersManagerObject.Log(STREAM("Particle Type = " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleTypeParam) << " Counter = " << ParticleCounter));

                                auto GeneIter = ParticlesKindsManagerObject.Genes.find(ParticleKindSpecialDataObject.GeneId);

                                if (GeneIter != ParticlesKindsManagerObject.Genes.end())
                                    TryToGenerateRandomParticlesForType(ParticleKindObject, Radius, RadiusSize, NumberOfErrors, GeneIter->second.NumId, GeneIter->second.Sequence, GeneIter->second.Sequence.length());
                                else
                                    TryToGenerateRandomParticlesForType(ParticleKindObject, Radius, RadiusSize, NumberOfErrors, 0, "NO GENE", GetSizeOfGeneratedParticle(ParticleKindSpecialDataObject.ParticleType));
                            }

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, string("Execution of generating all " + ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleTypeParam) + " particles from file has taken time: ").c_str(), "executing printing duration_time")));

        LoggersManagerObject.Log(STREAM("NUMBER OF ERRORS = " << NumberOfErrors));
    }
    CATCH("inserting new random particle for type")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::GenerateAllRealRandomParticles()
{
    try
    {
        ClearVoxelSpaceAndParticles();

        CellEngineUseful::SwitchOffLogs();

        const auto start_time = chrono::high_resolution_clock::now();

        InsertNewRandomParticlesForType(ParticlesTypes::Ribosome, 400, 400);

        InsertNewRandomParticlesForType(ParticlesTypes::MembraneProtein, 420, 45);
        InsertNewRandomParticlesForType(ParticlesTypes::RibosomeProtein, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::PolymeraseProtein, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymeraseProtein, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::ProteinFrac, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::OtherProtein, 400, 400);

        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell(579990, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartZPos, 2, 2, 2, 2, 2, 2, 2, 2);

        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::mRNA, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::rRNA, 400, 400);

        InsertNewRandomParticlesForType(ParticlesTypes::Basic, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::Lipid, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::Other, 400, 400);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating all random particles from file has taken time: ", "executing printing duration_time")));

        PrintNumberOfParticlesForAllMainTypesOfParticles();
    }
    CATCH("generating all random particles")
}