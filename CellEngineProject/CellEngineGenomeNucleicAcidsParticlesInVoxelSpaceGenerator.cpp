
#include "CellEngineParticlesKindsManager.h"

#include "CellEngineGenesPromotersAndGenesStartCodonsFinder.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

using namespace std;

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::EraseAllDNAParticles()
{
    try
    {
        LoggersManagerObject.Log(STREAM("Number of DNA particles to be removed = " << GetNumberOfParticlesWithChosenEntityId(CellEngineConfigDataObject.DNAIdentifier)));

        for (const auto& ParticleIndex : GetAllParticlesWithChosenEntityId(CellEngineConfigDataObject.DNAIdentifier))
            RemoveParticle(ParticleIndex, true);

        LoggersManagerObject.Log(STREAM("Number of DNA particles after remove = " << GetNumberOfParticlesWithChosenEntityId(CellEngineConfigDataObject.DNAIdentifier)));

        Genomes[0].clear();
        Genomes[1].clear();
    }
    CATCH("erasing all dna particles")
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

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateOneStrand(const EntityIdInt EntityId, string_view Sequence, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleStepX, const UnsignedInt ParticleStepY, const UnsignedInt ParticleStepZ)
{
    try
    {
        vector<UniqueIdInt> MockVector;
        Particle* ParticlePrev1 = nullptr;

        for (UnsignedInt SequenceIndex = 0; SequenceIndex < Sequence.size(); SequenceIndex++)
            ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, EntityId, CellEngineUseful::GetChainIdFromLetterForDNAorRNA(Sequence[SequenceIndex]), 0, 1, StartPosX + SequenceIndex * ParticleStepX, StartPosY + SequenceIndex * ParticleStepY, StartPosZ + SequenceIndex * ParticleStepZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, false, MockVector, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
    }
    CATCH("generating one strand")
}

Particle* CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateNucleotideParticle(Particle* ParticlePrev, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeThread, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, bool AddToGenome, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam, const bool LinkWithPreviousNucleotide)
{
    UnsignedInt ParticleIndex;

    try
    {
        ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), EntityId, ChainId, GenomeThread, GenomeIndex, 0, UniqueColorParam));

        if (LinkWithPreviousNucleotide == true)
            Particle::InsertAfterGivenNode(ParticlePrev, &GetParticleFromIndex(ParticleIndex));

        GenerateParticleVoxelsWhenSelectedSpaceIsFree(ParticleIndex, StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, 0, 0, 0, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, &CellEngineParticlesVoxelsShapesGenerator::CheckFreeSpaceInCuboidSelectedSpace, &CellEngineParticlesVoxelsShapesGenerator::SetValueToVoxelsForCuboidSelectedSpace);

        if (AddToGenome == true)
            Genome.emplace_back(ParticleIndex);
    }
    CATCH_AND_THROW("generating nucleotide particle")

    return &GetParticleFromIndex(ParticleIndex);
}

tuple<Particle*, Particle*> CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam, const bool Linked, const bool LinkWithPreviousNucleotide)
{
    try
    {
        ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, CellEngineConfigDataObject.DNAIdentifier, ChainId, 0, Genomes[0].size(), StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), LinkWithPreviousNucleotide);
        ParticlePrev2 = GenerateNucleotideParticle(ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainIdForDNAorRNA(ChainId), 1, Genomes[1].size(), StartPosX + AddSizeX, StartPosY + AddSizeY, StartPosZ + AddSizeZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[1], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), LinkWithPreviousNucleotide);

        if (Linked == true)
            PairedNucleotide<Particle>::LinkPairedNucleotides(ParticlePrev1, ParticlePrev2);
    }
    CATCH("generating two paired nucleotides")

    return { ParticlePrev1, ParticlePrev2 };
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

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell1RandomTurn(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
{
    try
    {
        NumberOfNucleotidesToBeGenerated = (NumberOfNucleotidesToBeGenerated == 0 ? GenomeLength : NumberOfNucleotidesToBeGenerated);

        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        CellEngineUseful::SwitchOffLogs();

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveAdvanceOfParticle_Uint64t(1, CellEngineConfigDataObject.NumberOfAdvanceVoxelsDuringGenerationOfRandomDNA);

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

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
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

        // Particle* ParticlePrev2 = nullptr;
        // UnsignedInt NotFreeSpaceFoundCounter = 0;
        // for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
        // {
        //     bool FoundFreeSpace = false;
        //     auto ParticleObject = GetParticleFromIndex(Genomes[0][GenomeIndex]);
        //     for (SignedInt PosX = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].X) - ParticleSizeX; PosX <= ParticleObject.ListOfVoxels[0].X + ParticleSizeX; PosX++)
        //         for (SignedInt PosY = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].Y) - ParticleSizeY; PosY <= ParticleObject.ListOfVoxels[0].Y + ParticleSizeY; PosY++)
        //             for (SignedInt PosZ = static_cast<SignedInt>(ParticleObject.ListOfVoxels[0].Z) - ParticleSizeZ; PosZ <= ParticleObject.ListOfVoxels[0].Z + ParticleSizeZ; PosZ++)
        //                 if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
        //                 {
        //                     ParticlePrev2 = GenerateNucleotideParticle(ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainIdForDNAorRNA(ParticleObject.ChainId), 1, GenomeIndex, PosX, PosY, PosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[1], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
        //                     FoundFreeSpace = true;
        //                     goto Outside;
        //                 }
        //     Outside:
        //     if (FoundFreeSpace == false)
        //         NotFreeSpaceFoundCounter++;
        //
        // }
        //
        // LoggersManagerObject.Log(STREAM("NUMBER OF NO FREE SPACE = " << NotFreeSpaceFoundCounter));

        GetMinMaxCoordinatesForDNA();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell1Vertical(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
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

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked()));

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    NumberOfGeneratedNucleotides++;

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

        GetMinMaxCoordinatesForDNA();

        CellEngineUseful::SwitchOnLogs();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell2Vertical(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
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

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked()));

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    NumberOfGeneratedNucleotides++;

                    ChainIdInt ChainId = UniformDistributionObjectChainOfParticle_Uint64t(mt64R);

                    if (RandomMoveDirection == 1 || RandomMoveDirection == 2)
                        tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, 1, ParticleSizeZ, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false, false);
                    else
                    if (RandomMoveDirection == 3 || RandomMoveDirection == 4 || RandomMoveDirection == 5 || RandomMoveDirection == 6)
                        tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, 1, ParticleSizeY, ParticleSizeZ, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), false, false);

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

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GetMinMaxCoordinatesForDNA()
{
    try
    {
        for (auto& LocalNewParticleIndex : Genomes[0])
            GetMinMaxCoordinatesForParticle(GetParticleFromIndex(LocalNewParticleIndex), true);
        for (auto& LocalNewParticleIndex : Genomes[1])
            GetMinMaxCoordinatesForParticle(GetParticleFromIndex(LocalNewParticleIndex), true);
    }
    CATCH("getting min max of coordinates for DNA")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::SaveGenomeDataToFile(UnsignedInt ParticleSize)
{
    try
    {
        ofstream FileToWriteGenome;
        FileToWriteGenome.open(CellEngineConfigDataObject.CellGenomePositionsFileName, ios_base::out | ios_base::trunc);
        FileToWriteGenome << to_string(ParticleSize) << endl;
        FileToWriteGenome << to_string(Genomes[0].size()) << endl;
        for (const auto& Nucleotide : Genomes[0])
        {
            EntityIdInt EntityId = GetParticleFromIndex(Nucleotide).EntityId;
            ChainIdInt ChainId = GetParticleFromIndex(Nucleotide).ChainId;
            UniqueIdInt GenomeIndex = GetParticleFromIndex(Nucleotide).GenomeIndex;
            UnsignedInt PosX = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].X;
            UnsignedInt PosY = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Y;
            UnsignedInt PosZ = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Z;
            FileToWriteGenome << to_string(EntityId) << "," << to_string(ChainId)  << "," << to_string(GenomeIndex) << "," << to_string(PosX) << "," << to_string(PosY) << "," << to_string(PosZ) << endl;
        }
        FileToWriteGenome.close();
    }
    CATCH("saving genome data to file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::ReadGenomeDataFromFile(bool Paired)
{
    try
    {
        CellEngineConfigDataObject.GenomeReadFromFile = true;

        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        string Line, Word;
        vector<string> Row;

        fstream FileToReadGenome(CellEngineConfigDataObject.CellGenomePositionsFileName, ios::in);

        getline(FileToReadGenome, Line);
        UnsignedInt ParticleSize = stoi(Line);
        getline(FileToReadGenome, Line);
        UnsignedInt GenomeSize = stoi(Line);

        UnsignedInt PrevStartPosX = 0;
        UnsignedInt GenomeIndex = 0;

        Particle* ParticlePrev1 = nullptr;
        Particle* ParticlePrev2 = nullptr;

        while(getline(FileToReadGenome, Line))
        {
            Row.clear();
            stringstream Str(Line);

            while(getline(Str, Word, ','))
                Row.push_back(Word);

            auto StartPosX = static_cast<UnsignedInt>(static_cast<float>(stoi(Row[3])) / CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles);
            auto StartPosY = static_cast<UnsignedInt>(static_cast<float>(stoi(Row[4])) / CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles);
            auto StartPosZ = static_cast<UnsignedInt>(static_cast<float>(stoi(Row[5])) / CellEngineConfigDataObject.DivisionFactorForReadingPositionsOfParticles);

            if (Paired == true)
            {
                if (abs(static_cast<long>(PrevStartPosX - StartPosX)) > 0)
                    tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, ParticleSize, 1, ParticleSize, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true, true);
                else
                    tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, 1, ParticleSize, ParticleSize, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true, true);
            }
            else
                ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, stoi(Row[0]), stoi(Row[1]), 0, stoi(Row[2]), stoi(Row[3]), stoi(Row[4]), stoi(Row[5]), ParticleSize, ParticleSize, ParticleSize, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);

            PrevStartPosX = StartPosX;

            GenomeIndex++;
        }

        FileToReadGenome.close();

        //GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("reading genome data from file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::ReadGenomeSequenceFromFile(const bool Paired)
{
    try
    {
        fstream FileToReadGenome(CellEngineConfigDataObject.CellGenomeSequenceFileName, ios::in);

        getline(FileToReadGenome, GenomesLines[0]);
        if (Paired == true)
        {
            GenomesLines[1] = "";
            for (const auto& GenomeLetter : GenomesLines[0])
                GenomesLines[1] += CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomeLetter)));
        }

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
        {
            GetParticleFromIndex(Genomes[0][GenomeIndex + 1]).ChainId = CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]);
            if (Paired == true)
                GetParticleFromIndex(Genomes[1][GenomeIndex + 1]).ChainId = CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]));
        }
    }
    CATCH("reading real genome data from file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::ReadGenomeSequenceFromFastaFile()
{
    try
    {
        fstream FileToReadGenome(CellEngineConfigDataObject.CellGenomeSequenceFastaFileName, ios::in);

        string Line;
        while(getline(FileToReadGenome, Line))
            GenomesLines[0] += Line.substr(0, Line.size() - 1);

        ofstream FileToSaveGenome(CellEngineConfigDataObject.CellGenomeSequenceFastaOneFileName, ios::out);
        FileToSaveGenome << GenomesLines[0];
        FileToSaveGenome.close();
    }
    CATCH("reading real genome data from file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::TestGeneratedGenomeCorrectness(const UnsignedInt ParticleSize)
{
    try
    {
        bool FoundBreak = false;

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
            if (GenomeIndex > 3)
            {
                if ((abs(GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].X - GetParticleFromIndex(Genomes[0][GenomeIndex - 1]).ListOfVoxels[0].X) != ParticleSize) &&
                    (abs(GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].Y - GetParticleFromIndex(Genomes[0][GenomeIndex - 1]).ListOfVoxels[0].Y) != ParticleSize) &&
                    (abs(GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].Z - GetParticleFromIndex(Genomes[0][GenomeIndex - 1]).ListOfVoxels[0].Z) != ParticleSize))
                {
                    LoggersManagerObject.Log(STREAM("GenomeIndex = " << GenomeIndex << " ChainId = " << GetParticleFromIndex(Genomes[0][GenomeIndex]).ChainId << " Letter = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(Genomes[0][GenomeIndex]).ChainId)));
                    LoggersManagerObject.Log(STREAM("DIFF X = " << GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].X << " " << GetParticleFromIndex(Genomes[0][GenomeIndex - ParticleSize]).ListOfVoxels[0].X));
                    LoggersManagerObject.Log(STREAM("DIFF Y = " << GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].Y << " " << GetParticleFromIndex(Genomes[0][GenomeIndex - ParticleSize]).ListOfVoxels[0].Y));
                    LoggersManagerObject.Log(STREAM("DIFF Z = " << GetParticleFromIndex(Genomes[0][GenomeIndex]).ListOfVoxels[0].Z << " " << GetParticleFromIndex(Genomes[0][GenomeIndex - ParticleSize]).ListOfVoxels[0].Z));
                    FoundBreak = true;
                    break;
                }
            }

        if (FoundBreak == false)
            LoggersManagerObject.Log(STREAM("Genome is continuous and correctly generated - OK!"));
    }
    CATCH("testing generated genome correctness")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms() const
{
    TestSeveralDifferentKindsOfPromotersFindingsAndTerminatorFindingsAlgorithms(GenomesLines, Genomes);
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::FindInterGenesSequences()
{
    PrintInterGenesSequencesFromGenesData();
}