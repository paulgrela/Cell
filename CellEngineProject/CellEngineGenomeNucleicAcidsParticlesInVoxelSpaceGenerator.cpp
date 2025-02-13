
#include "CellEngineParticlesKindsManager.h"

#include "CellEngineGenesPromotersAndGenesStartCodonsFinder.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

using namespace std;

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateParticleVoxelsWhenSelectedSpaceIsFreeBasic(const UnsignedInt ParticleIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleStepX, const UnsignedInt ParticleStepY, const UnsignedInt ParticleStepZ)
{
    GenerateParticleVoxelsWhenSelectedSpaceIsFree(ParticleIndex, StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, 0, 0, 0, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::UpdateRandomPositions(const UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, const UnsignedInt Size)
{
    switch (RandomMoveDirection)
    {
        case 1 : RandomPosX += Size; break;
        case 2 : RandomPosX -= Size; break;
        case 3 : RandomPosY += Size; break;
        case 4 : RandomPosY -= Size; break;
        case 5 : RandomPosZ += Size; break;
        case 6 : RandomPosZ -= Size; break;
        default : break;
    }
}

bool CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::TestFormerForbiddenPositions(unordered_set<string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt Size)
{
    UnsignedInt PosX = RandomPosX;
    UnsignedInt PosY = RandomPosY;
    UnsignedInt PosZ = RandomPosZ;

    UpdateRandomPositions(RandomMoveDirection, PosX, PosY, PosZ, Size);

    return TestedFormerForbiddenPositions.contains(to_string(PosX) + "|" + to_string(PosY) + "|" + to_string(PosZ));
}

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::EraseLastRandomDNAParticle(vector<UniqueIdInt>& Genome)
{
    vector3_16 LastLocalRandomPos;

    try
    {
        UnsignedInt PreviousParticleIndex = Genome.back();
        Genome.pop_back();

        LastLocalRandomPos = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0];

        RemoveParticle(PreviousParticleIndex, true);

        LoggersManagerObject.Log(STREAM("ERASED PARTICLE PreviousParticleIndex = " << PreviousParticleIndex << " RandomPosX = " << LastLocalRandomPos.X << " RandomPosY = " << LastLocalRandomPos.Y << " RandomPosZ = " << LastLocalRandomPos.Z));
    }
    CATCH("erasing last random particle")

    return { LastLocalRandomPos.X, LastLocalRandomPos.Y, LastLocalRandomPos.Z };
}

UnsignedInt Sqr(const UnsignedInt Value)
{
    return Value * Value;
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::AddPairedNucleotidesToExistingNucleotidesFromGenome(const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ)
{
    try
    {
        Particle* ParticlePrev2 = nullptr;
        UnsignedInt NotFreeSpaceFoundCounter = 0;
        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
        {
            bool FoundFreeSpace = false;
            auto ParticleObject = GetParticleFromIndex(Genomes[0][GenomeIndex]);
            for (auto PosX = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].X - ParticleSizeX); PosX <= ParticleObject.ListOfVoxels[0].X + ParticleSizeX; PosX++)
                for (auto PosY = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].Y - ParticleSizeY); PosY <= ParticleObject.ListOfVoxels[0].Y + ParticleSizeY; PosY++)
                    for (auto PosZ = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].Z - ParticleSizeZ); PosZ <= ParticleObject.ListOfVoxels[0].Z + ParticleSizeZ; PosZ++)

                        if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
                        {
                            ParticlePrev2 = GenerateNucleotideParticle(ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainIdForDNAorRNA(ParticleObject.ChainId), 1, GenomeIndex, PosX, PosY, PosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[1], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
                            FoundFreeSpace = true;
                            goto Outside;
                        }
            Outside:
            if (FoundFreeSpace == false)
                NotFreeSpaceFoundCounter++;
        }

        LoggersManagerObject.Log(STREAM("NOT FREE SPACE COUNTER = " << NotFreeSpaceFoundCounter));
    }
    CATCH("adding paired nucleotides to existing nucleotides from genome")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell1or2RandomTurn(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleSize1, const UnsignedInt ParticleSize2, const UnsignedInt ParticleSize3, const UnsignedInt ParticleSize4, const UnsignedInt ParticleSize5, const bool Paired)
{
    try
    {
        NumberOfNucleotidesToBeGenerated = (NumberOfNucleotidesToBeGenerated == 0 ? GenomeLength : NumberOfNucleotidesToBeGenerated);

        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        CellEngineUseful::SwitchOffLogs();

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveAdvanceOfParticle_Uint64t(1, CellEngineConfigDataObject.SizeOfAdvanceDuringGenerationOfRandomDNA);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [&RandomMovesDirections]() { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

        auto timeStart = clock();

        Particle* ParticlePrev1 = nullptr;

        while (NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
        {
            if ((clock() - timeStart) / CLOCKS_PER_SEC >= 15)
                break;

            UnsignedInt RandomMoveDirection = 0;

            do
            {
                RandomMoveDirection = UniformDistributionObjectMoveOfParticle_Uint64t(mt64R);

                if (TestFormerForbiddenPositions(TestedFormerForbiddenPositions, RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize1) == true)
                    RandomMovesDirections[RandomMoveDirection - 1] = 1;

                LoggersManagerObject.Log(STREAM("RandomMoveDirection = " << RandomMoveDirection << " " << RandomMovesDirections[RandomMoveDirection - 1] << " " << CheckIfAllRandomMovesDirectionsWereChecked()));
            }
            while (RandomMovesDirections[RandomMoveDirection - 1] == 1 && CheckIfAllRandomMovesDirectionsWereChecked() == false);

            if (RandomMovesDirections[RandomMoveDirection - 1] == 0)
            {
                UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize2);

                RandomMovesDirections[RandomMoveDirection - 1] = 1;
            }

            UnsignedInt MaxNumberOfMoveAdvanceSteps = UniformDistributionObjectMoveAdvanceOfParticle_Uint64t(mt64R);
            UnsignedInt NumberOfMoveAdvanceSteps = 0;

            bool EmptyVoxelSpaceForNewNucleotideBool = true;

            while (EmptyVoxelSpaceForNewNucleotideBool == true && NumberOfMoveAdvanceSteps < MaxNumberOfMoveAdvanceSteps && NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
            {
                bool EmptyVoxelSpaceForNewNucleotideBoolX = CheckFreeSpaceInCuboidSelectedSpace(RandomPosX, RandomPosY, RandomPosZ, 1, 1, 1, 1 * ParticleSizeX, ParticleSizeY, ParticleSizeZ, GetZeroSimulationSpaceVoxel());
                bool EmptyVoxelSpaceForNewNucleotideBoolY = CheckFreeSpaceInCuboidSelectedSpace(RandomPosX, RandomPosY, RandomPosZ, 1, 1, 1, ParticleSizeX, 1 * ParticleSizeY, ParticleSizeZ, GetZeroSimulationSpaceVoxel());
                bool EmptyVoxelSpaceForNewNucleotideBoolZ = CheckFreeSpaceInCuboidSelectedSpace(RandomPosX, RandomPosY, RandomPosZ, 1, 1, 1, ParticleSizeX, ParticleSizeY, 1 * ParticleSizeZ, GetZeroSimulationSpaceVoxel());

                EmptyVoxelSpaceForNewNucleotideBool = EmptyVoxelSpaceForNewNucleotideBoolX || EmptyVoxelSpaceForNewNucleotideBoolY || EmptyVoxelSpaceForNewNucleotideBoolZ;

                if (EmptyVoxelSpaceForNewNucleotideBool == false)
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize3);

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked()));

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    NumberOfGeneratedNucleotides++;
                    NumberOfMoveAdvanceSteps++;

                    ChainIdInt ChainId = UniformDistributionObjectChainOfParticle_Uint64t(mt64R);

                    ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, CellEngineConfigDataObject.DNAIdentifier, ChainId, 0, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false);

                    fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);

                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize5);

                    LoggersManagerObject.Log(STREAM("ADDED PARTICLE NumberOfGeneratedNucleotides = " << NumberOfGeneratedNucleotides << " RX = " << RandomPosX << " RY = " << RandomPosY << " RZ = " << RandomPosZ));
                }
            }

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked() == true)
            {
                NumberOfGeneratedNucleotides--;

                tie(RandomPosX, RandomPosY, RandomPosZ) = EraseLastRandomDNAParticle(Genomes[0]);

                TestedFormerForbiddenPositions.insert(to_string(RandomPosX) + "|" + to_string(RandomPosY) + "|" + to_string(RandomPosZ));

                fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);
            }

            LoggersManagerObject.Log(STREAM("END OF GOING IN ONE DIRECTION"));
        }

        LoggersManagerObject.Log(STREAM("END OF FINDING SINGLE STRAND DNA"));

        CellEngineUseful::SwitchOnLogs();

        if (Paired == true)
            AddPairedNucleotidesToExistingNucleotidesFromGenome(ParticleSizeX, ParticleSizeY, ParticleSizeZ);

        GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell1or2Vertical(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleSize1, const UnsignedInt ParticleSize2, const UnsignedInt ParticleSize3, const UnsignedInt ParticleSize4, const UnsignedInt ParticleSize5, const bool Paired)
{
    try
    {
        NumberOfNucleotidesToBeGenerated = (NumberOfNucleotidesToBeGenerated == 0 ? GenomeLength : NumberOfNucleotidesToBeGenerated);

        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        CellEngineUseful::SwitchOffLogs();

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [&RandomMovesDirections]() { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

        auto timeStart = clock();

        Particle* ParticlePrev1 = nullptr;
        Particle* ParticlePrev2 = nullptr;

        while (NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
        {
            if ((clock() - timeStart) / CLOCKS_PER_SEC >= 15)
                break;

            UnsignedInt RandomMoveDirection = 0;

            do
            {
                RandomMoveDirection = UniformDistributionObjectMoveOfParticle_Uint64t(mt64R);

                if (TestFormerForbiddenPositions(TestedFormerForbiddenPositions, RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize1) == true)
                    RandomMovesDirections[RandomMoveDirection - 1] = 1;

                LoggersManagerObject.Log(STREAM("RandomMoveDirection = " << RandomMoveDirection << " " << RandomMovesDirections[RandomMoveDirection - 1] << " " << CheckIfAllRandomMovesDirectionsWereChecked()));
            }
            while (RandomMovesDirections[RandomMoveDirection - 1] == 1 && CheckIfAllRandomMovesDirectionsWereChecked() == false);

            if (RandomMovesDirections[RandomMoveDirection - 1] == 0)
            {
                UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize2);

                RandomMovesDirections[RandomMoveDirection - 1] = 1;
            }

            bool EmptyVoxelSpaceForNewNucleotideBool = true;
            while (EmptyVoxelSpaceForNewNucleotideBool == true && NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
            {
                EmptyVoxelSpaceForNewNucleotideBool = CheckFreeSpaceInCuboidSelectedSpace(RandomPosX, RandomPosY, RandomPosZ, 1, 1, 1, ParticleSizeX, ParticleSizeY, ParticleSizeZ, GetZeroSimulationSpaceVoxel());

                if (EmptyVoxelSpaceForNewNucleotideBool == false)
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize3);

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.SizeOfSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked()));

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    NumberOfGeneratedNucleotides++;

                    ChainIdInt ChainId = UniformDistributionObjectChainOfParticle_Uint64t(mt64R);

                    if (Paired == true)
                    {
                        if (RandomMoveDirection == 1 || RandomMoveDirection == 2)
                            tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, 1, ParticleSizeZ, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false, false);
                        else
                        if (RandomMoveDirection == 3 || RandomMoveDirection == 4 || RandomMoveDirection == 5 || RandomMoveDirection == 6)
                            tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, 1, ParticleSizeY, ParticleSizeZ, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false, false);
                    }
                    else
                        ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, CellEngineConfigDataObject.DNAIdentifier, ChainId, 0, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false);

                    fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);

                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize5);

                    LoggersManagerObject.Log(STREAM("ADDED PARTICLE NumberOfGeneratedNucleotides = " << NumberOfGeneratedNucleotides << " RX = " << RandomPosX << " RY = " << RandomPosY << " RZ = " << RandomPosZ));
                }
            }

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked() == true)
            {
                NumberOfGeneratedNucleotides--;

                tie(RandomPosX, RandomPosY, RandomPosZ) = EraseLastRandomDNAParticle(Genomes[0]);

                TestedFormerForbiddenPositions.insert(to_string(RandomPosX) + "|" + to_string(RandomPosY) + "|" + to_string(RandomPosZ));

                if (Paired == true)
                    EraseLastRandomDNAParticle(Genomes[1]);

                fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);
            }

            LoggersManagerObject.Log(STREAM("END OF GOING IN ONE DIRECTION"));
        }

        GetMinMaxCoordinatesForDNA();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}
