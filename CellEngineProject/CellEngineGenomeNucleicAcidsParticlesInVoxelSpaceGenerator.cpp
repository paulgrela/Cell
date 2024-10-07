
#include "CellEngineConstants.h"
#include "CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator.h"

#include "CellEngineParticlesKindsManager.h"

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
    CATCH("generating particle")

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

    return (TestedFormerForbiddenPositions.find(to_string(PosX) + "|" + to_string(PosY) + "|" + to_string(PosZ)) != TestedFormerForbiddenPositions.end());
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

UnsignedInt Sqr(UnsignedInt Value)
{
    return Value * Value;
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
{
    try
    {
        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        CellEngineUseful::SwitchOffLogs();

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [](const vector<UnsignedInt>& RandomMovesDirections) { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

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

                LoggersManagerObject.Log(STREAM("RandomMoveDirection = " << RandomMoveDirection << " " << RandomMovesDirections[RandomMoveDirection - 1] << " " << CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections)));
            }
            while (RandomMovesDirections[RandomMoveDirection - 1] == 1 && CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections) == false);

            if (RandomMovesDirections[RandomMoveDirection - 1] == 0)
            {
                UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize2);

                RandomMovesDirections[RandomMoveDirection - 1] = 1;
            }

            bool EmptyVoxelSpaceForNewNucleotideBool = true;
            while (EmptyVoxelSpaceForNewNucleotideBool == true && NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
            {
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + ParticleSizeX; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + ParticleSizeY; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + ParticleSizeZ; PosZ++)
                            if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                            {
                                LoggersManagerObject.Log(STREAM("NOT EMPTY A POS = " << PosX << " " << PosY << " " << PosZ << " " << GetSpaceVoxel(PosX, PosY, PosZ)));

                                EmptyVoxelSpaceForNewNucleotideBool = false;

                                goto BreakOutOfLoop;
                            }

                BreakOutOfLoop:;

                if (EmptyVoxelSpaceForNewNucleotideBool == false)
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize3);

                if (EmptyVoxelSpaceForNewNucleotideBool == true && sqrt(Sqr(RandomPosX - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosY - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2)) + Sqr(RandomPosZ - (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) >= CellEngineConfigDataObject.RadiusOfCellForDNA)
                {
                    EmptyVoxelSpaceForNewNucleotideBool = false;
                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, -ParticleSize4);
                }

                LoggersManagerObject.Log(STREAM("EmptyVoxelSpaceForNewNucleotideBool = " << EmptyVoxelSpaceForNewNucleotideBool << " " << CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections)));

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

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections) == true)
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
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        ofstream FileToWriteGenome;
        FileToWriteGenome.open(FileName, ios_base::out | ios_base::trunc);
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

        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        fstream FileToReadGenome(FileName, ios::in);

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

            UnsignedInt StartPosX = stoi(Row[3]);
            UnsignedInt StartPosY = stoi(Row[4]);
            UnsignedInt StartPosZ = stoi(Row[5]);

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

        GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("reading genome data from file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::ReadGenomeSequenceFromFile()
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_SEQUENCE.DAT");

        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome, GenomesLines[0]);
        GenomesLines[1] = "";
        for (const auto& GenomeLetter : GenomesLines[0])
            GenomesLines[1] += CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomeLetter)));

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
        {
            GetParticleFromIndex(Genomes[0][GenomeIndex + 1]).ChainId = CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]);
            GetParticleFromIndex(Genomes[1][GenomeIndex + 1]).ChainId = CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]));
        }
    }
    CATCH("reading real genome data from file")
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::ReadGenomeSequenceFromFastaFile()
{
    try
    {
        string FileNameIn = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("sequence.fasta");

        fstream FileToReadGenome(FileNameIn, ios::in);

        string Line;
        while(getline(FileToReadGenome, Line))
            GenomesLines[0] += Line.substr(0, Line.size() - 1);

        string FileNameOut = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("sequence.fasta.one");
        ofstream FileToSaveGenome(FileNameOut, ios::out);
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













inline void CompareSequences1(const std::vector<ChainIdInt>& TemplateSequence, const std::vector<ChainIdInt>& NucleotidesSequenceToCompareVector, bool& FoundSequenceNotFit)
{
    if (NucleotidesSequenceToCompareVector.size() >= TemplateSequence.size())
    {
        //LoggersManagerObject.Log(STREAM("LOOP COMPARISON SIZE = " << std::to_string(NucleotidesSequenceToCompareVector.size()) << " " << std::to_string(TemplateSequence.size())));

        for (UnsignedInt NucleotideNum = 0; NucleotideNum < TemplateSequence.size(); NucleotideNum++)
            if (CellEngineUseful::CompareIUPACNucleotideCode(TemplateSequence[NucleotideNum], NucleotidesSequenceToCompareVector[NucleotideNum]) == false)
            {
                //LoggersManagerObject.Log(STREAM("LOOP COMPARISON BREAK = " << std::to_string(NucleotideNum) << "#"));

                FoundSequenceNotFit = true;
                break;
            }
    }
    else
        FoundSequenceNotFit = true;
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::FindInterGenesSequences() const
{
    try
    {
        vector<Gene> LocalGenesVector;
        transform(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), back_inserter(LocalGenesVector), [](const auto& Gene){ return Gene.second; });
        sort(LocalGenesVector.begin(), LocalGenesVector.end(), [](const Gene& G1, const Gene& G2){ return G1.StartPosInGenome < G2.StartPosInGenome; });

        for (const auto& Gene : LocalGenesVector)
        {
            UnsignedInt StartPos = Gene.StartPosInGenome;
            UnsignedInt EndPos = Gene.EndPosInGenome;
            LoggersManagerObject.Log(STREAM("Gene = " << Gene.NumId << " " << StartPos << " " << EndPos));
        }
    }
    CATCH("finding inter genes sequences")
}

void FindPromotersAndStartCodons1(const std::string& genome)
{
    const std::string promoter_35_box = "TTGACA";  // Consensus -35 region
    //const std::string promoter_10_box = "TATAAT";  // Consensus -10 region
    const std::string promoter_10_box = "TATAA";  // Consensus -10 region
    const std::vector<std::string> start_codons = {"ATG", "GTG", "TTG"};  // Common bacterial start codons

    const int promoter_search_distance = 50;  // Max distance between -35 and -10 boxes
    const int max_distance_to_start = 50;     // Max distance from -10 box to start codon

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t i = 0; i < genome.size() - promoter_35_box.size(); ++i)
    {
        // Look for -35 box
        if (genome.substr(i, promoter_35_box.size()) == promoter_35_box)
        {
            // Search for -10 box downstream
            for (size_t j = i + promoter_35_box.size(); j < i + promoter_search_distance && j < genome.size(); ++j)
            {
                if (genome.substr(j, promoter_10_box.size()) == promoter_10_box)
                {
                    // Search for start codon downstream of -10 box
                    for (size_t k = j + promoter_10_box.size(); k < j + promoter_10_box.size() + max_distance_to_start && k < genome.size(); ++k)
                    {
                        for (const std::string& start_codon : start_codons)
                        {
                            if (genome.substr(k, start_codon.size()) == start_codon)
                            {
                                // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                                // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                                auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k, &ScopeForStartGeneIndex](const auto& G1) { return k - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= k + ScopeForStartGeneIndex; });
                                //auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k](const auto& G1) { return G1.second.StartPosInGenome == k - 1 || G1.second.StartPosInGenome == k || G1.second.StartPosInGenome == k + 1; });
                                if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                {
                                    //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                    NumberOfFoundPromotersConfirmed++;
                                }
                                NumberOfFoundPromoters++;

                                break;
                            }
                        }
                    }
                    break;  // Found the -10 box, no need to search further in this window
                }
            }
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found ConfirmedPromoters = " << NumberOfFoundPromotersConfirmed << endl));
}

int CalculateSimilarity2(const std::string& sequence1, const std::string& sequence2)
{
    int mismatches = 0;
    for (size_t i = 0; i < sequence1.size(); ++i)
    {
        if (sequence1[i] != sequence2[i])
        {
            mismatches++;
        }
    }
    return mismatches;
}

// Function to find promoters and start codons with fuzzy matching
void FindPromotersAndStartCodons2(const std::string& genome)
{
    const std::vector<std::string> promoter_35_boxes = {"TTGACA", "TTTACA", "CTGACA", "TTGATA"};  // Common -35 regions with variability
    //const std::vector<std::string> promoter_35_boxes = {"TTGACA"};  // Common -35 regions with variability
    //const std::vector<std::string> promoter_35_boxes = {"TTGACA", "TTTACA"};  // Common -35 regions with variability

    const std::vector<std::string> promoter_10_boxes = {"TATAAT", "TATAAA", "TATGAT", "TATATT"};  // Common -10 regions with variability
    //const std::vector<std::string> promoter_10_boxes = {"TATAAT", "TATAAA"};  // Common -10 regions with variability
    //const std::vector<std::string> promoter_10_boxes = {"TATAA"};  // Common -10 regions with variability

    const std::vector<std::string> start_codons = {"ATG", "GTG", "TTG"};  // Common bacterial start codons

    const int promoter_search_distance = 50;  // Max distance between -35 and -10 boxes 50
    const int max_distance_to_start = 100;    // Max distance from -10 box to start codon 100
    const int max_mismatches_35 = 1;          // Allowable mismatches for -35 region
    const int max_mismatches_10 = 1;          // Allowable mismatches for -10 region

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 5;

    for (size_t i = 0; i < genome.size() - 6; ++i)
    {
        // Look for -35 box with some flexibility (fuzzy match)
        for (const std::string& promoter_35_box : promoter_35_boxes)
        {
            if (CalculateSimilarity2(genome.substr(i, 6), promoter_35_box) <= max_mismatches_35)
            {
                // Search for -10 box downstream with flexibility
                for (size_t j = i + 6; j < i + promoter_search_distance && j < genome.size() - 6; ++j)
                {
                    for (const std::string& promoter_10_box : promoter_10_boxes)
                    {
                        if (CalculateSimilarity2(genome.substr(j, 6), promoter_10_box) <= max_mismatches_10)
                        {
                            // Search for start codon downstream of -10 box
                            for (size_t k = j + 6; k < j + 6 + max_distance_to_start && k < genome.size() - 3; ++k)
                            {
                                for (const std::string& start_codon : start_codons)
                                {
                                    if (genome.substr(k, 3) == start_codon)
                                    //if (genome.substr(k, start_codon.size()) == start_codon)
                                    {
                                        // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                                        // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                                        //auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k](const auto& G1) { return G1.second.StartPosInGenome == k - 1 || G1.second.StartPosInGenome == k || G1.second.StartPosInGenome == k + 1; });
                                        auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k, &ScopeForStartGeneIndex](const auto& G1) { return k - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= k + ScopeForStartGeneIndex; });
                                        if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                        {
                                            //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                            FoundGenesCounter[FindResult->second.NumId]++;
                                            NumberOfFoundPromotersConfirmed++;
                                        }
                                        NumberOfFoundPromoters++;

                                        break;
                                    }
                                }
                            }
                            break;  // Found a match for -10 box, move on
                        }
                    }
                }
            }
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

int CalculateSimilarity3(const std::string& sequence1, const std::string& sequence2)
{
    int matches = 0;
    for (size_t i = 0; i < sequence1.size(); ++i)
    {
        if (sequence1[i] == sequence2[i])
        {
            matches++;
        }
    }
    return matches;
}

void FindPromotersAndStartCodons3(const std::string& genome)
{
    // Consensus sequences for -35 and -10 regions
    const std::vector<std::string> promoter_35_boxes = {"TTGACA", "TTTACA", "CTGACA", "TTGATA"};
    const std::vector<std::string> promoter_10_boxes = {"TATAAT", "TATAAA", "TATGAT", "TATATT"};
    const std::vector<std::string> start_codons = {"ATG", "GTG", "TTG"};  // Common bacterial start codons

    const int promoter_search_distance = 19;  // Distance between -35 and -10 should be 16-18bp
    const int max_distance_to_start = 50;    // Max distance from -10 box to start codon
    const int min_distance_to_start = 10;    // Min distance from -10 box to start codon
    const int threshold_similarity_35 = 5;   // Minimum similarity score for -35 region (out of 6)
    const int threshold_similarity_10 = 5;   // Minimum similarity score for -10 region (out of 6)

    UniqueIdInt NumberOfFoundPromoters = 0;
    UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
    map<GeneIdInt, UnsignedInt> FoundGenesCounter;
    UnsignedInt ScopeForStartGeneIndex = 3;

    for (size_t i = 0; i < genome.size() - 6; ++i)
    {
        // Look for -35 box with fuzzy matching
        for (const std::string& promoter_35_box : promoter_35_boxes)
        {
            if (CalculateSimilarity3(genome.substr(i, 6), promoter_35_box) >= threshold_similarity_35)
            {
                // Search for -10 box downstream within the correct distance window
                for (size_t j = i + 16; j <= i + promoter_search_distance && j < genome.size() - 6; ++j)
                {
                    for (const std::string& promoter_10_box : promoter_10_boxes)
                    {
                        if (CalculateSimilarity3(genome.substr(j, 6), promoter_10_box) >= threshold_similarity_10)
                        {
                            // Search for start codon downstream of -10 box within valid distance
                            for (size_t k = j + min_distance_to_start; k < j + max_distance_to_start && k < genome.size() - 3; ++k)
                            {
                                for (const std::string& start_codon : start_codons)
                                {
                                    if (genome.substr(k, 3) == start_codon)
                                    {
                                        // std::cout << "Promoter found: -35 at " << i << " and -10 at " << j << std::endl;
                                        // std::cout << "Start codon (" << start_codon << ") found at: " << k << std::endl;
                                        // std::cout << "----\n";

                                        auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k, &ScopeForStartGeneIndex](const auto& G1) { return k - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= k + ScopeForStartGeneIndex; });
                                        if (FindResult != ParticlesKindsManagerObject.Genes.end())
                                        {
                                            //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                            FoundGenesCounter[FindResult->second.NumId]++;
                                            NumberOfFoundPromotersConfirmed++;
                                        }
                                        NumberOfFoundPromoters++;

                                        break;
                                    }
                                }
                            }
                            break;  // Found valid -10 box
                        }
                    }
                }
            }
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));
}

int CalculateSimilarity4(const std::string& sequence1, const std::string& sequence2)
{
    int matches = 0;
    for (size_t i = 0; i < sequence1.size(); ++i)
    {
        if (sequence1[i] == sequence2[i])
        {
            matches++;
        }
    }
    return matches;
}

// Function to search for promoters for each gene start
void FindPromotersForGenes1(const std::string& genome, const std::vector<size_t>& GenesStarts)
{
    // Promoter sequences for -10 and -35 boxes (you can add more variants as needed)
    const std::vector<std::string> promoter_35_boxes = {"TTGACA", "TTTACA", "CTGACA", "TTGATA"};
    const std::vector<std::string> promoter_10_boxes = {"TATAAT", "TATAAA", "TATGAT", "TATATT"};

    const int min_distance_to_start = 10;    // Min distance from -10 box to gene start
    const int max_distance_to_start = 50;    // Max distance from -10 box to gene start
    const int promoter_search_distance = 19; // Max distance between -35 and -10 boxes
    const int threshold_similarity_35 = 5;   // Minimum similarity score for -35 region (out of 6)
    const int threshold_similarity_10 = 5;   // Minimum similarity score for -10 region (out of 6)

    UniqueIdInt NumberOfFoundPromoters = 0;
    map<GeneIdInt, UnsignedInt> FoundPromotersCounter;

    // Iterate over each gene start
    for (size_t geneStart : GenesStarts)
    {
        bool promoter_found = false;

        // Search for the -10 box (reverse search from gene start)
        for (int i = geneStart - min_distance_to_start; i >= std::max(0, static_cast<int>(geneStart - max_distance_to_start)); --i)
        {
            for (const std::string& promoter_10_box : promoter_10_boxes)
            {
                if (CalculateSimilarity4(genome.substr(i, 6), promoter_10_box) >= threshold_similarity_10)
                {
                    // We found a -10 box, now look upstream for the -35 box
                    for (int j = i - promoter_search_distance; j >= std::max(0, i - promoter_search_distance - 5); --j)
                    {
                        for (const std::string& promoter_35_box : promoter_35_boxes)
                        {
                            if (CalculateSimilarity4(genome.substr(j, 6), promoter_35_box) >= threshold_similarity_35)
                            {
                                // Found both -35 and -10 boxes, print the promoter info
                                LoggersManagerObject.Log(STREAM("Promoter found for gene starting at " << geneStart << std::endl));
                                LoggersManagerObject.Log(STREAM("-35 box at position: " << j << std::endl));
                                LoggersManagerObject.Log(STREAM("-10 box at position: " << i << std::endl));
                                FoundPromotersCounter[i]++;
                                NumberOfFoundPromoters++;

                                promoter_found = true;
                                break;
                            }
                        }
                        if (promoter_found)
                            break;
                    }
                }
                if (promoter_found)
                    break;
            }
            if (promoter_found)
                break;
        }

        if (!promoter_found)
        {
            LoggersManagerObject.Log(STREAM("No promoter found for gene starting at " << geneStart << std::endl));
        }
    }

    LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
    LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundPromotersCounter.size() << endl));
}

void CellEngineGenomeNucleicAcidsParticlesInVoxelSpaceGenerator::CheckGenomePromoters() const
{
    try
    {
        vector<Gene> LocalGenesVector;
        transform(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), back_inserter(LocalGenesVector), [](const auto& Gene){ return Gene.second; });
        sort(LocalGenesVector.begin(), LocalGenesVector.end(), [](const Gene& G1, const Gene& G2){ return G1.StartPosInGenome < G2.StartPosInGenome; });

        for (const auto& Gene : LocalGenesVector)
        {
            UnsignedInt StartPos = Gene.StartPosInGenome;
            UnsignedInt EndPos = Gene.EndPosInGenome;
            LoggersManagerObject.Log(STREAM("Gene = " << Gene.NumId << " " << StartPos << " " << EndPos));
        }

        UnsignedInt NumberOfFoundPromoterSequences = 0;
        string AttachPolymeraseToDNAStartSequenceStr = "TATAAT";

        UniqueIdInt NumberOfFoundPromoters = 0;
        UniqueIdInt NumberOfFoundPromotersConfirmed = 0;
        map<GeneIdInt, UnsignedInt> FoundGenesCounter;
        UnsignedInt ScopeForStartGeneIndex = 5;
        const std::vector<std::string> start_codons = {"ATG", "GTG", "TTG"};
        const int max_distance_to_start = 100;

        auto AttachPolymeraseToDNAStartSequence = CellEngineUseful::ConvertStringSequenceToChainIdSequence(AttachPolymeraseToDNAStartSequenceStr);

        for (UnsignedInt StartGenomeIndex = 0; StartGenomeIndex < Genomes[0].size() - AttachPolymeraseToDNAStartSequenceStr.size(); StartGenomeIndex++)
        {
            auto NucleotidesSequenceToCompareVector = CellEngineUseful::ConvertStringSequenceToChainIdSequence(GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()));

            bool FoundSequenceNotFit = false;
            CompareSequences1(AttachPolymeraseToDNAStartSequence, NucleotidesSequenceToCompareVector, FoundSequenceNotFit);
            if (FoundSequenceNotFit == false)
            {
                //LoggersManagerObject.Log(STREAM("Promoter sequence found = " << GenomesLines[0].substr(StartGenomeIndex, AttachPolymeraseToDNAStartSequenceStr.size()) << " StartGenomeIndex = " << StartGenomeIndex));

                // for (UnsignedInt GeneIndex = 1; GeneIndex < LocalGenesVector.size(); GeneIndex++)
                //     if (StartGenomeIndex < LocalGenesVector[GeneIndex].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length() && StartGenomeIndex > LocalGenesVector[GeneIndex - 1].StartPosInGenome - AttachPolymeraseToDNAStartSequenceStr.length())
                //     {
                //         LoggersManagerObject.Log(STREAM("Promoter sequence for gene = " << GeneIndex << " StartGenomeIndex for sequence " << AttachPolymeraseToDNAStartSequenceStr << " = " << StartGenomeIndex << " Gene.StartPosInGenome " << LocalGenesVector[GeneIndex].StartPosInGenome << " Diff = " << LocalGenesVector[GeneIndex].StartPosInGenome - StartGenomeIndex));
                //         break;
                //     }

                for (size_t k = StartGenomeIndex + 6; k < StartGenomeIndex + 6 + max_distance_to_start && k < Genomes[0].size() - 3; ++k)
                {
                    for (const std::string& start_codon : start_codons)
                    {
                        if (GenomesLines[0].substr(k, 3) == start_codon)
                        //if (genome.substr(k, start_codon.size()) == start_codon)
                        {
                            // LoggersManagerObject.Log(STREAM("Promoter found at -35: " << i << " and -10: " << j));
                            // LoggersManagerObject.Log(STREAM("Start codon (" << start_codon << ") found at: " << k << std::endl));

                            //auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k](const auto& G1) { return G1.second.StartPosInGenome == k - 1 || G1.second.StartPosInGenome == k || G1.second.StartPosInGenome == k + 1; });
                            auto FindResult = find_if(ParticlesKindsManagerObject.Genes.begin(), ParticlesKindsManagerObject.Genes.end(), [&k, &ScopeForStartGeneIndex](const auto& G1) { return k - ScopeForStartGeneIndex <= G1.second.StartPosInGenome && G1.second.StartPosInGenome <= k + ScopeForStartGeneIndex; });
                            if (FindResult != ParticlesKindsManagerObject.Genes.end())
                            {
                                //LoggersManagerObject.Log(STREAM("GENE NR = " << FindResult->second.NumId << " " << FindResult->second.StrId << std::endl));
                                FoundGenesCounter[FindResult->second.NumId]++;
                                NumberOfFoundPromotersConfirmed++;
                            }
                            NumberOfFoundPromoters++;

                            break;
                        }
                    }
                }

                NumberOfFoundPromoterSequences++;
            }
        }

        LoggersManagerObject.Log(STREAM("Number Of Found Promoter Sequences = " << NumberOfFoundPromoterSequences << endl));

        LoggersManagerObject.Log(STREAM("Number Of Found Promoters = " << NumberOfFoundPromoters));
        LoggersManagerObject.Log(STREAM("Number Of Found Confirmed Promoters = " <<  NumberOfFoundPromotersConfirmed));
        LoggersManagerObject.Log(STREAM("Number Of Different Confirmed Promoters = " << FoundGenesCounter.size() << endl));

        FindPromotersAndStartCodons1(GenomesLines[0]);
        FindPromotersAndStartCodons2(GenomesLines[0]);
        FindPromotersAndStartCodons3(GenomesLines[0]);

        std::vector<size_t> GenesStarts;
        for (const auto& Gene : ParticlesKindsManagerObject.Genes)
            GenesStarts.emplace_back(Gene.second.StartPosInGenome);
        FindPromotersForGenes1(GenomesLines[0], GenesStarts);
    }
    CATCH("checking genome promoters")
}
