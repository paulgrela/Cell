
#include <cmath>

#include "DateTimeUtils.h"

#include "CellEngineConstants.h"
#include "CellEngineDataFile.h"
#include "CellEngineParticlesKindsManager.h"
#include "CellEngineRealRandomParticlesInVoxelSpaceGenerator.h"

using namespace std;

UnsignedInt CellEngineRealRandomParticlesInVoxelSpaceGenerator::GetNumberOfRealParticlesOfKind(ParticlesTypes ParticleTypeParam)
{
    UnsignedInt ParticlesCounter = 0;

    try
    {
        for (const auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId != 0)
            {
                auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.second.EntityId);
                if (ParticleKindObject.ParticleKindSpecialDataSector.empty() == false)
                    for (const auto& ParticleKindSpecialDataObject : ParticleKindObject.ParticleKindSpecialDataSector)
                        if (ParticleKindSpecialDataObject.ParticleType == ParticleTypeParam)
                            ParticlesCounter++;
            }
    }
    CATCH("getting number of particles of type")

    return ParticlesCounter;
}

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
    CATCH("getting number of particles kinds")

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

        PrintInformationAboutRibosomesProteins();

        LoggersManagerObject.Log(STREAM("Total Number Of All Particles Kinds = " << TotalNumberOfAllParticles));
        LoggersManagerObject.Log(STREAM("Total Number Of All Particles = " << Particles.size()));

        UnsignedInt OldParticlesKindsCounter = 0;
        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.EntityId < StartParticleKindId)
                OldParticlesKindsCounter++;
        LoggersManagerObject.Log(STREAM("Old Particles Kinds Counter = " << OldParticlesKindsCounter));

        UniqueIdInt OldParticlesCounter = 0;
        set<UniqueIdInt> TestSet;
        for (const auto& ParticleObject : Particles)
        {
            if (ParticleObject.second.EntityId < StartParticleKindId)
            {
                TestSet.insert(ParticleObject.second.EntityId);
                OldParticlesCounter++;
            }
        }
        LoggersManagerObject.Log(STREAM("Old Particles Counter = " << OldParticlesCounter));

        for (const auto& TestSetElementObject : TestSet)
            LoggersManagerObject.Log(STREAM("Test set element id = " << TestSetElementObject));

        UnsignedInt CountMultiParticleKind = 0;
        for (auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.size() > 1)
                if (ParticleKindObject.second.ParticleKindSpecialDataSector[0].ParticleType == ParticlesTypes::Basic && ParticleKindObject.second.ParticleKindSpecialDataSector[1].ParticleType != ParticlesTypes::Basic)
                {
                    CountMultiParticleKind++;
                    ParticleKindObject.second.ParticleKindSpecialDataSector.erase(ParticleKindObject.second.ParticleKindSpecialDataSector.begin());
                }

        LoggersManagerObject.Log(STREAM("CountMultiParticleKind = " << CountMultiParticleKind));

        LoggersManagerObject.Log(STREAM("Number Of Ribosomes P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::Ribosome)));
        LoggersManagerObject.Log(STREAM("Number Of RNAPolymerases P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::RNAPolymerase)));
        LoggersManagerObject.Log(STREAM("Number Of DNAPolymerases P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::DNAPolymerase)));

        LoggersManagerObject.Log(STREAM("Number Of Membrane Proteins P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::MembraneProtein)));
        LoggersManagerObject.Log(STREAM("Number Of Ribosomes Proteins P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::RibosomeProtein)));
        LoggersManagerObject.Log(STREAM("Number Of RNA Polymerase Proteins P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::RNAPolymeraseProtein)));
        LoggersManagerObject.Log(STREAM("Number Of Polymerase Proteins P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::PolymeraseProtein)));
        LoggersManagerObject.Log(STREAM("Number Of Proteins Frac P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::ProteinFrac)));
        LoggersManagerObject.Log(STREAM("Number Of Other Proteins P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::OtherProtein)));

        LoggersManagerObject.Log(STREAM("Number Of TRNA_uncharged P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::tRNA_uncharged)));
        LoggersManagerObject.Log(STREAM("Number Of TRNA_charged P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::tRNA_charged)));
        LoggersManagerObject.Log(STREAM("Number Of MRNA P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::mRNA)));
        LoggersManagerObject.Log(STREAM("Number Of RRNA P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::rRNA)));

        LoggersManagerObject.Log(STREAM("Number Of Basic P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::Basic)));
        LoggersManagerObject.Log(STREAM("Number Of Lipids P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::Lipid)));
        LoggersManagerObject.Log(STREAM("Number Of Other P = " << GetNumberOfRealParticlesOfKind(ParticlesTypes::Other)));
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

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::TryToGenerateRandomParticlesForType(const pair<const EntityIdInt, ParticleKind>& ParticleKindObject, const EntityIdInt EntityId, const UnsignedInt Radius, const UnsignedInt RadiusSize, UnsignedInt& NumberOfErrors, const GeneIdInt GeneNumId, const string& GeneSequence, const UnsignedInt GeneSequenceLength)
{
    try
    {
        LoggersManagerObject.Log(STREAM("Gene Number = " << GeneNumId));

        auto GeneralSize = static_cast<UnsignedInt>(pow(GeneSequenceLength, (1.0 / 3.0)));
        UnsignedInt Size = GeneralSize > 10 ? GeneralSize - 3 : GeneralSize;

        UnsignedInt PosX, PosY, PosZ;
        bool TryResult = false;

        UnsignedInt TryInsertNewParticleCounter = 0;

        while (TryInsertNewParticleCounter < MaxNumberOfTriesToInsertNewParticle && TryResult == false)
        {
            tie(PosX, PosY, PosZ) = GetRandomPositionInsideSphere(Radius, RadiusSize);
            //TryResult = GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ParticleKindObject.second.EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, GeneSequence, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosX, PosY, PosZ, Size, Size, Size, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
            TryResult = GenerateParticleVoxelsWhenSelectedSpaceIsFree(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), EntityId, 1, -1, 1, ParticleKindObject.second.ElectricCharge, GeneSequence, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))), PosX, PosY, PosZ, Size, Size, Size, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceForSphereSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForSphereSelectedSpace);
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

EntityIdInt GetParticleKindIdForRNA(const EntityIdInt EntityId, const ParticlesTypes ParticleTypeParam, const bool ModifyRNAParticleKindId)
{
    if (ModifyRNAParticleKindId == true)
    {
        if (ParticleTypeParam == ParticlesTypes::mRNA || ParticleTypeParam == ParticlesTypes::rRNA || ParticleTypeParam == ParticlesTypes::tRNA_uncharged || ParticleTypeParam == ParticlesTypes::tRNA_charged)
            return CellEngineConfigDataObject.RNAIdentifier;
        else
            return EntityId;
    }
    else
        return EntityId;
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::InsertNewRandomParticlesForType(const ParticlesTypes ParticleTypeParam, bool ModifyParticleKindIdForRNA, const UnsignedInt Radius, const UnsignedInt RadiusSize)
{
    try
    {
        UnsignedInt NumberOfErrors = 0;

        CellEngineUseful::SwitchOffLogs();

        const auto start_time = chrono::high_resolution_clock::now();

        for (const auto& ParticleKindObject : ParticlesKindsManagerObject.ParticlesKinds)
            if (ParticleKindObject.second.ParticleKindSpecialDataSector.empty() == false)
                for (const auto& ParticleKindSpecialDataObject: ParticleKindObject.second.ParticleKindSpecialDataSector)
                    if (ParticleKindSpecialDataObject.ParticleType == ParticleTypeParam)
                        for (UnsignedInt ParticleCounter = 1; ParticleCounter <= ParticleKindSpecialDataObject.CounterAtStartOfSimulation; ParticleCounter++)
                        {
                            LoggersManagerObject.Log(STREAM("Particle Type = " << ParticlesKindsManagerObject.ConvertParticleTypeToString(ParticleKindSpecialDataObject.ParticleType) << " Counter = " << ParticleCounter));

                            EntityIdInt EntityId = GetParticleKindIdForRNA(ParticleKindObject.second.EntityId, ParticleKindSpecialDataObject.ParticleType, ModifyParticleKindIdForRNA);

                            auto GeneIter = ParticlesKindsManagerObject.Genes.find(ParticleKindSpecialDataObject.GeneId);
                            if (GeneIter != ParticlesKindsManagerObject.Genes.end())
                                TryToGenerateRandomParticlesForType(ParticleKindObject, EntityId, Radius, RadiusSize, NumberOfErrors, GeneIter->second.NumId, GeneIter->second.Sequence, GeneIter->second.Sequence.length());
                            else
                                TryToGenerateRandomParticlesForType(ParticleKindObject, EntityId, Radius, RadiusSize, NumberOfErrors, 0, "NO GENE", GetSizeOfGeneratedParticle(ParticleKindSpecialDataObject.ParticleType));
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

        InsertNewRandomParticlesForType(ParticlesTypes::Ribosome, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymerase, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::DNAPolymerase, false, 400, 400);

        InsertNewRandomParticlesForType(ParticlesTypes::MembraneProtein, false, 420, 45);
        InsertNewRandomParticlesForType(ParticlesTypes::RibosomeProtein, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::PolymeraseProtein, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::RNAPolymeraseProtein, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::ProteinFrac, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::OtherProtein, false, 400, 400);

        CellEngineDataFileObjectPointer->CellEngineVoxelSimulationSpaceObjectPointer->GenerateRandomDNAInWholeCell(579990, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartXPos + 3, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartYPos, CellEngineConfigDataObject.VoxelSimulationSpaceSelectionStartZPos, 2, 2, 2, 2, 2, 2, 2, 2);

        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::mRNA, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::rRNA, false, 400, 400);

        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_uncharged, false, true, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::tRNA_charged, false, true, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::mRNA, false, true, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::rRNA, false, true, 400);

        InsertNewRandomParticlesForType(ParticlesTypes::Basic, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::Lipid, false, 400, 400);
        InsertNewRandomParticlesForType(ParticlesTypes::Other, false, 400, 400);

        const auto stop_time = chrono::high_resolution_clock::now();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM(GetDurationTimeInOneLineStr(start_time, stop_time, "Execution of generating all random particles from file has taken time: ", "executing printing duration_time")));

        PrintNumberOfParticlesForAllMainTypesOfParticles();
    }
    CATCH("generating all random particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveParticlesWithChosenEntityId(const EntityIdInt EntityId, const UnsignedInt NumberOfParticlesToBeRemoved)
{
    try
    {
        UnsignedInt ParticleIndexCounter = 0;
        for (const auto& ParticleIndex : GetAllParticlesWithChosenEntityId(EntityId))
        {
            if (ParticleIndexCounter + 1 > NumberOfParticlesToBeRemoved)
                break;

            RemoveParticle(ParticleIndex, true);
            ParticleIndexCounter++;
        }
    }
    CATCH("removing number of particles with chosen entity id")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAllParticlesWithChosenEntityId(const EntityIdInt EntityId)
{
    try
    {
        for (const auto& ParticleIndex : GetAllParticlesWithChosenEntityId(EntityId))
            RemoveParticle(ParticleIndex, true);
    }
    CATCH("removing all of particles with chosen entity id")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveParticlesWithChosenParticleType(const ParticlesTypes ParticleTypeParam, const UnsignedInt NumberOfParticlesToBeRemoved)
{
    try
    {
        UnsignedInt ParticleIndexCounter = 0;

        for (const auto& ParticleIndex : GetAllParticlesWithChosenParticleType(ParticleTypeParam))
        {
            if (ParticleIndexCounter + 1 > NumberOfParticlesToBeRemoved)
                break;

            RemoveParticle(ParticleIndex, true);
            ParticleIndexCounter++;
        }
    }
    CATCH("removing number of particles with chosen entity id")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAllParticlesWithChosenParticleType(const ParticlesTypes ParticleTypeParam)
{
    try
    {
        for (const auto& ParticleIndex : GetAllParticlesWithChosenParticleType(ParticleTypeParam))
            RemoveParticle(ParticleIndex, true);
    }
    CATCH("removing all of particles with chosen particle type")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAllRNAParticles()
{
    try
    {
        RemoveAllmRNAParticles();
        RemoveAlltRNAParticles();
        RemoveAllrRNAParticles();
    }
    CATCH("removing all rna particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAllmRNAParticles()
{
    try
    {
        RemoveAllParticlesWithChosenParticleType(ParticlesTypes::mRNA);
    }
    CATCH("removing all mrna particles")
}
void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAlltRNAParticles()
{
    try
    {
        RemoveAllParticlesWithChosenParticleType(ParticlesTypes::tRNA_uncharged);
        RemoveAllParticlesWithChosenParticleType(ParticlesTypes::tRNA_charged);
    }
    CATCH("removing all trna particles")
}

void CellEngineRealRandomParticlesInVoxelSpaceGenerator::RemoveAllrRNAParticles()
{
    try
    {
        RemoveAllParticlesWithChosenParticleType(ParticlesTypes::rRNA);
    }
    CATCH("removing all trna particles")
}