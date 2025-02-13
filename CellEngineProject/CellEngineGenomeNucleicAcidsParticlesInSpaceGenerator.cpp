
#include "CellEngineParticlesKindsManager.h"

#include "CellEngineGenesPromotersAndGenesStartCodonsFinder.h"
#include "CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator.h"

#ifdef USING_MODULES
import CellEngineColors;
#else
#include "CellEngineColors.h"
#endif

using namespace std;

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::EraseAllDNAParticles()
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::GenerateOneStrand(const EntityIdInt EntityId, const string_view Sequence, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleStepX, const UnsignedInt ParticleStepY, const UnsignedInt ParticleStepZ)
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

Particle* CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::GenerateNucleotideParticle(Particle* ParticlePrev, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeThread, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const bool AddToGenome, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam, const bool LinkWithPreviousNucleotide)
{
    UnsignedInt ParticleIndex;

    try
    {
        ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), EntityId, ChainId, GenomeThread, GenomeIndex, 0, UniqueColorParam));

        if (LinkWithPreviousNucleotide == true)
            Particle::InsertAfterGivenNode(ParticlePrev, &GetParticleFromIndex(ParticleIndex));

        GenerateParticleVoxelsWhenSelectedSpaceIsFreeBasic(ParticleIndex, StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, 0, 0, 0);

        if (AddToGenome == true)
            Genome.emplace_back(ParticleIndex);
    }
    CATCH_AND_THROW("generating nucleotide particle")

    return &GetParticleFromIndex(ParticleIndex);
}

tuple<Particle*, Particle*> CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam, const bool Linked, const bool LinkWithPreviousNucleotide)
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::GetMinMaxCoordinatesForDNA()
{
    try
    {
        for (const auto& LocalNewParticleIndex : Genomes[0])
            GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(GetParticleFromIndex(LocalNewParticleIndex), &Particle::ListOfVoxels, &ParticleKind::ListOfVoxels, true);
        for (const auto& LocalNewParticleIndex : Genomes[1])
            GetMinMaxCoordinatesForParticle<UnsignedInt, vector3_16>(GetParticleFromIndex(LocalNewParticleIndex), &Particle::ListOfVoxels, &ParticleKind::ListOfVoxels, true);
    }
    CATCH("getting min max of coordinates for DNA")
}

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::SaveGenomeDataToFile(UnsignedInt ParticleSize)
{
    try
    {
        ofstream FileToWriteGenome;
        FileToWriteGenome.open(CellEngineConfigDataObject.CellGenomePositionsFileName, ios_base::out | ios_base::trunc);
        FileToWriteGenome << to_string(ParticleSize) << endl;
        FileToWriteGenome << to_string(Genomes[0].size()) << endl;
        for (const auto& Nucleotide : Genomes[0])
        {
            const EntityIdInt EntityId = GetParticleFromIndex(Nucleotide).EntityId;
            const ChainIdInt ChainId = GetParticleFromIndex(Nucleotide).ChainId;
            const UniqueIdInt GenomeIndex = GetParticleFromIndex(Nucleotide).GenomeIndex;
            const UnsignedInt PosX = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].X;
            const UnsignedInt PosY = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Y;
            const UnsignedInt PosZ = GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Z;
            FileToWriteGenome << to_string(EntityId) << "," << to_string(ChainId)  << "," << to_string(GenomeIndex) << "," << to_string(PosX) << "," << to_string(PosY) << "," << to_string(PosZ) << endl;
        }
        FileToWriteGenome.close();
    }
    CATCH("saving genome data to file")
}

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::ReadGenomeDataFromFile(const bool PairedBool)
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

            if (PairedBool == true)
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::ReadGenomeSequenceFromFile(const bool Paired)
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::ReadGenomeSequenceFromFastaFile()
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::TestGeneratedGenomeCorrectness(const UnsignedInt ParticleSize)
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

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::TestDifferentKindsOfPromotersFindingsAndTerminatorsFindingsAlgorithms() const
{
    TestSeveralDifferentKindsOfPromotersFindingsAndTerminatorFindingsAlgorithms(GenomesLines);
}

void CellEngineGenomeNucleicAcidsParticlesInSpaceGenerator::FindInterGenesSequences()
{
    PrintInterGenesSequencesFromGenesData();
}