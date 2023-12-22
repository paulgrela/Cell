
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "FileUtils.h"
#include "DoublyLinkedList.h"
#include "DestinationPlatform.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"
#include "CellEngineOpenGLVisualiserOfVoxelSimulationSpace.h"

using namespace std;

[[nodiscard]] float CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace;
};

[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace) + (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension / 2);
};

inline SimulationSpaceVoxel& CellEngineVoxelSimulationSpace::GetSpaceVoxel(UnsignedInt x, UnsignedInt y, UnsignedInt z)
{
    return CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? (*static_cast<Space_2048_2048_2048*>(SpacePointer))[x][y][z] : (*static_cast<Space_1024_1024_1024*>(SpacePointer))[x][y][z];
}

inline Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndex(const UniqueIdInt ParticleIndex)
{
    return Particles[ParticleIndex];
}

inline SimulationSpaceVoxel GetZeroSimulationSpaceVoxel()
{
    return 0;
}

SimulationSpaceVoxel CellEngineVoxelSimulationSpace::GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndexInSimulationSpaceVoxel(const UniqueIdInt ParticleIndex)
{
    return GetParticleFromIndex(ParticleIndex);
}

void CellEngineVoxelSimulationSpace::MakeZeroVoxelsForSelectedSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += 1)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += 1)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += 1)
                    GetSpaceVoxel(PosX, PosY, PosZ) = GetZeroSimulationSpaceVoxel();
    }
    CATCH("making zero voxels for selected space")
}

void CellEngineVoxelSimulationSpace::SetAllVoxelsInListOfVoxelsToValue(std::vector<vector3_16>& ListOfVoxels, SimulationSpaceVoxel SimulationSpaceVoxelValue)
{
    try
    {
        for (auto& Voxel : ListOfVoxels)
            GetSpaceVoxel(Voxel.X, Voxel.Y, Voxel.Z) = SimulationSpaceVoxelValue;
    }
    CATCH("making all zero voxels in list of voxels")
}

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        MakeZeroVoxelsForSelectedSpace(0, 0, 0, 1, 1, 1, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension, CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension);

        Genomes.resize(2);
        GenomesLines.resize(1);
    }
    CATCH("execution of constructor of voxel simulation space")
}

CellEngineVoxelSimulationSpace::~CellEngineVoxelSimulationSpace()
{
    try
    {
        free(SpacePointer);
    }
    CATCH("execution of destructor of voxel simulation space")
}

void CellEngineVoxelSimulationSpace::InitiateFreeParticleIndexes()
{
    try
    {
        for (UnsignedInt FreeIndex = MaxParticleIndex + 100'000'000; FreeIndex >= MaxParticleIndex + 1; FreeIndex--)
            FreeIndexesOfParticles.push(FreeIndex);
    }
    CATCH("initiating free particle indexes")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForAllParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            GetMinMaxCoordinatesForParticle(ParticleObject.second);
    }
    CATCH("getting min max coordinates for all particles")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForParticle(Particle& ParticleObject)
{
    try
    {
        UnsignedInt ParticleXMin{}, ParticleXMax{}, ParticleYMin{}, ParticleYMax{}, ParticleZMin{}, ParticleZMax{};
        ParticleXMin = ParticleYMin = ParticleZMin = 10000;
        ParticleXMax = ParticleYMax = ParticleZMax = 0;

        for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
            GetMinMaxOfCoordinates(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z, ParticleXMin, ParticleXMax, ParticleYMin, ParticleYMax, ParticleZMin, ParticleZMax);

        ParticleObject.SetCenterCoordinates(ParticleXMin + (ParticleXMax - ParticleXMin) / 2, ParticleYMin + (ParticleYMax - ParticleYMin) / 2, ParticleZMin + (ParticleZMax - ParticleZMin) / 2);

        auto& ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);
        if (ParticleKindObject.ListOfVoxels.empty() == true)
        {
            for (auto& VoxelCoordinates : ParticleObject.ListOfVoxels)
                ParticleKindObject.ListOfVoxels.emplace_back(VoxelCoordinates.X - ParticleXMin, VoxelCoordinates.Y - ParticleYMin, VoxelCoordinates.Z - ParticleZMin);
            ParticleKindObject.XSizeDiv2 = (ParticleXMax - ParticleXMin) / 2;
            ParticleKindObject.YSizeDiv2 = (ParticleYMax - ParticleYMin) / 2;
            ParticleKindObject.ZSizeDiv2 = (ParticleZMax - ParticleZMin) / 2;
        }
    }
    CATCH("getting min max coordinates for one particle")
}

void CellEngineVoxelSimulationSpace::PreprocessData()
{
    try
    {
        InitiateFreeParticleIndexes();
        GetMinMaxCoordinatesForAllParticles();
    }
    CATCH("preprocessing data for voxel simulation space")
}

void CellEngineVoxelSimulationSpace::SetStartValuesForSpaceMinMax()
{
    XMin = YMin = ZMin = 10000;
    XMax = YMax = ZMax = 0;
}

inline void CellEngineVoxelSimulationSpace::GetMinMaxOfCoordinates(const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ, UnsignedInt& XMinParam, UnsignedInt& XMaxParam, UnsignedInt& YMinParam, UnsignedInt& YMaxParam, UnsignedInt& ZMinParam, UnsignedInt& ZMaxParam)
{
    try
    {
        XMinParam = min(PosX, XMinParam);
        XMaxParam = max(PosX, XMaxParam);
        YMinParam = min(PosY, YMinParam);
        YMaxParam = max(PosY, YMaxParam);
        ZMinParam = min(PosZ, ZMinParam);
        ZMaxParam = max(PosZ, ZMaxParam);
    }
    CATCH("getting min max of coordinates")
}

[[nodiscard]] stringstream CellEngineVoxelSimulationSpace::PrintSpaceMinMaxValues() const
{
    stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << to_string(XMin) << " ][ Xmax = " << to_string(XMax) << " ][ Ymin = " << to_string(YMin) << " ][ Ymax = " << to_string(YMax) << " ][ Zmin = " << to_string(ZMin) << " ][ Zmax = " << to_string(XMax) << " ] " << endl;
    return ss;
}

void CellEngineVoxelSimulationSpace::CountStatisticsOfVoxelSimulationSpace()
{
    try
    {
        for(UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for(UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for(UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                {
                    if (GetParticleFromIndex(GetSpaceVoxel(PosX, PosY, PosZ)).EntityId != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
    CATCH("counting statistics of voxel simulation space")
}

void CellEngineVoxelSimulationSpace::SetAtomInVoxelSimulationSpace(const UniqueIdInt ParticleIndex, const CellEngineAtom& AppliedAtom)
{
    try
    {
        UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(PosX, PosY, PosZ, XMin, XMax, YMin, YMax, ZMin, ZMax);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
        }
    }
    CATCH("setting atom in voxel simulation space")
}

UniqueIdInt CellEngineVoxelSimulationSpace::AddNewParticle(const Particle& ParticleParam)
{
    Particles[ParticleParam.Index] = ParticleParam;
    return MaxParticleIndex = ParticleParam.Index;
}

UniqueIdInt CellEngineVoxelSimulationSpace::GetFreeIndexesOfParticleSize()
{
    return FreeIndexesOfParticles.size();
}

UniqueIdInt CellEngineVoxelSimulationSpace::GetNewFreeIndexOfParticle()
{
    if (FreeIndexesOfParticles.empty() == false)
    {
        UniqueIdInt FreeIndexOfParticle = FreeIndexesOfParticles.top();
        FreeIndexesOfParticles.pop();
        return FreeIndexOfParticle;
    }
    else
    {
        LoggersManagerObject.Log(STREAM("Lack of new free indexes of particles"));
        return MaxParticleIndex + 1;
    }
}

void CellEngineVoxelSimulationSpace::AddBasicParticlesKindsAndReactions()
{
    try
    {
        ParticlesKindsManagerObject.AddParticleKind({ 0, "Water", "H2O", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 1, "Glucose", "C6H12O6", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 2, "Oxygen2", "02", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 3, "Carbon dioxide", "CO2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 4, "Eten", "CH2CH2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 5, "Ethanol", "CH3CH2(OH)", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 6, "Propen", "CH3CHCH2", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 7, "HX", "HX", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 8, "2Halogenopropan", "CH3CHXCH3", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 9, "Test", "TEST", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 10, "Ethylene", "CH2CH2O", 0 });
        ParticlesKindsManagerObject.AddParticleKind({ 11, "Oxygen", "0", 0 });

        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", 0 });

        const string DNASequenceForTestFindingDNA = "TACAAAAAAAGAGGTGTTAGC";
        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", DNASequenceForTestFindingDNA, 0 });

        const string DNASequence1ForTestCutLink1 = "TACAAAAAAAGAGGTGTT";
        ParticlesKindsManagerObject.AddParticleKind({ 10002, "DNA", DNASequence1ForTestCutLink1, 0 });
        //const string DNASequence2ForTestCutLink1 = "AGCTCTTATTA";
        const string DNASequence2ForTestCutLink1 = "GCTCTTATTAT";
        ParticlesKindsManagerObject.AddParticleKind({ 10003, "DNA", DNASequence2ForTestCutLink1, 0 });

        const string RNASequence1ForTestCutLink1 = "TACAAAAAAAGAGGTGTT";
        ParticlesKindsManagerObject.AddParticleKind({ 10004, "RNA", RNASequence1ForTestCutLink1, 0 });
        const string RNASequence2ForTestCutLink1 = "AGCTCTTATTA";
        ParticlesKindsManagerObject.AddParticleKind({ 10005, "RNA", RNASequence2ForTestCutLink1, 0 });

        const string DNASequence1ForTestCutLink2 = "TACAAAAAAAGAGGTGTT";
        ParticlesKindsManagerObject.AddParticleKind({ 10006, "DNA", DNASequence1ForTestCutLink2, 0 });
        const string DNASequence2ForTestCutLink2 = "AGC";
        ParticlesKindsManagerObject.AddParticleKind({ 10007, "DNA", DNASequence2ForTestCutLink2, 0 });

        const string RNASequenceForTestCutCrisper = "AGC";
        ParticlesKindsManagerObject.AddParticleKind({ 10008, "RNA", RNASequenceForTestCutCrisper, 0 });

        AddChemicalReaction(Reaction(101, "", "CH3CH2(OH) + DNA + ", { { 5, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 10, 1, "", true } }));


        AddChemicalReaction(Reaction(1, "CUT 1", "CH3CH2(OH) + DNA + ", { { 5, 1, "", false }, { 10002, 1, DNASequence1ForTestCutLink1, false } }, {}));

        AddChemicalReaction(Reaction(2, "LINK 1", "CH3CH2(OH) + DNA + DNA + ", { { 5, 1, "", false }, { 10002, 1, DNASequence1ForTestCutLink1, false }, { 10003, 1, DNASequence2ForTestCutLink1, false } }, {}));
//
//
//        AddChemicalReaction(Reaction(3, "CUT 1", "CH3CHCH2 + RNA + ", { { 6, 1, "", false }, { 10004, 1, RNASequence1ForTestCutLink1, false } }, {}));
//
//        AddChemicalReaction(Reaction(4, "LINK 1", "CH3CHCH2 + RNA + RNA + ", { { 6, 1, "", false }, { 10004, 1, RNASequence1ForTestCutLink1, false }, { 10005, 1, DNASequence2ForTestCutLink1, false } }, {}));
//
//
//        AddChemicalReaction(Reaction(5, "CUT 2", "CH3CHCH2 + DNA + ", 5, { { 6, 1, "", false }, { 10006, 1, DNASequence1ForTestCutLink2, false }, { 10007, 1, DNASequence2ForTestCutLink1, false } }, { { 86, 1, "", true } }));
//
//        AddChemicalReaction(Reaction(6, "LINK 2", "CH3CHCH2 + DNA + DNA + ", { { 6, 1, "", false }, { 10006, 1, DNASequence1ForTestCutLink2, false }, { 10007, 1, DNASequence2ForTestCutLink1, false } }, {}));
//
//        AddChemicalReaction(Reaction(5, "CUT 2", "CH3CHCH2 + DNA + ", { { 6, 1, "", false }, { 10006, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }));
//
//        AddChemicalReaction(Reaction(6, "LINK 2", "CH3CHCH2 + DNA + DNA + ", { { 6, 1, "", false }, { 10006, 1, DNASequence1ForTestCutLink2, false }, { 10007, 1, DNASequence2ForTestCutLink1, false } }, {}));
//
//
//        AddChemicalReaction(Reaction(7, "CUT CRISPER 2", "CH3CHCH2 + RNA + DNA + ", { { 6, 1, "", false }, { 10008, 1, RNASequenceForTestCutCrisper, false } }, { { 90, 1, "", true } }));
//
//
//        AddChemicalReaction(Reaction(8, "POLYMERASE DNA RUN", "POLYMERASE + DNA + NUCLEOTIDE_DNA + ", { { 100, 1, "", false }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, {}));
//        AddChemicalReaction(Reaction(9, "POLYMERASE DNA ADD", "POLYMERASE + DNA + NUCLEOTIDE_DNA + ", { { 100, 1, "", false }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, {}));
//        AddChemicalReaction(Reaction(10, "POLYMERASE DNA END", "POLYMERASE + DNA + NUCLEOTIDE_DNA + ", { { 100, 1, "", false }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, {}));
//
//        AddChemicalReaction(Reaction(11, "POLYMERASE RNA RUN", "POLYMERASE + DNA + NUCLEOTIDE_RNA + ", { { 101, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 101, 1, "", true } }));
//        AddChemicalReaction(Reaction(12, "POLYMERASE RNA ADD", "POLYMERASE + DNA + NUCLEOTIDE_RNA + ", { { 101, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 101, 1, "", true } }));
//        AddChemicalReaction(Reaction(13, "POLYMERASE RNA END", "POLYMERASE + DNA + NUCLEOTIDE_RNA + ", { { 101, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 101, 1, "", true } }));
//
//        AddChemicalReaction(Reaction(14, "SEPARATE DNA STRAND", "SEPARATOR_PROTEIN + DNA + ", { { 102, 1, "", false }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, {}));
//        AddChemicalReaction(Reaction(15, "LINK DNA STRAND", "LINKER_PROTEIN + DNA + DNA + ", { { 103, 1, "", false }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, {}));
//
//        AddChemicalReaction(Reaction(1001, "", "C6H12O6 + O2 + ", { { 1, 1, "", true }, { 2, 6, "", true } }, { { 3, 6, "", true }, { 0, 6, "", true } }));
//        AddChemicalReaction(Reaction(1002, "", "CH2CH2 + H2O + ", { { 4, 1, "", true }, { 0, 1, "", true } }, { { 5, 1, "", true } }));
//        AddChemicalReaction(Reaction(1003, "", "CH3CHCH2 + HX + ", { { 6, 1, "", true }, { 7, 1, "", true } }, { { 8, 1, "", true } }));
//        AddChemicalReaction(Reaction(1004, "", "CH2CH2 + O + ", { { 4,  1, "", true }, { 11, 1, "", true } }, { { 10, 1, "", true } }));

        PreprocessChemicalReactions();
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << MaxParticleIndex << " " << (MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(0, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        AddBasicParticlesKindsAndReactions();

        vector<UnsignedInt> LocalNewParticlesIndexes;

        for (UniqueIdInt ParticleNumber = 1; ParticleNumber <= NumberOfRandomParticles; ParticleNumber++)
            LocalNewParticlesIndexes.emplace_back(AddNewParticle(Particle(GetNewFreeIndexOfParticle(), UniformDistributionObjectObjectOfParticle_Uint64t(mt64R), 0, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));

        vector<vector3_16> FilledVoxelsForRandomParticle;

        MakeZeroVoxelsForSelectedSpace(StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);

        for (const auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
        {
            UnsignedInt RandomPosX = UniformDistributionObjectX_Uint64t(mt64R);
            UnsignedInt RandomPosY = UniformDistributionObjectY_Uint64t(mt64R);
            UnsignedInt RandomPosZ = UniformDistributionObjectZ_Uint64t(mt64R);

            UnsignedInt RandomSizeOfParticle = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

            FilledVoxelsForRandomParticle.clear();

            if (RandomPosX + RandomSizeOfParticle < StartXPosParam + SizeXParam && RandomPosY + RandomSizeOfParticle < StartYPosParam + SizeYParam && RandomPosZ + RandomSizeOfParticle < StartZPosParam + SizeZParam)
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + RandomSizeOfParticle; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + RandomSizeOfParticle; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + RandomSizeOfParticle; PosZ++)
                        {
                            if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
                            {
                                FilledVoxelsForRandomParticle.emplace_back(PosX, PosY, PosZ);
                                GetSpaceVoxel(PosX, PosY, PosZ) = LocalNewParticleIndex;
                            }
                            else
                            {
                                SetAllVoxelsInListOfVoxelsToValue(FilledVoxelsForRandomParticle,GetZeroSimulationSpaceVoxel());

                                FilledVoxelsForRandomParticle.clear();

                                goto NextRandomParticleOutsideLoopLabel;
                            }
                        }
            NextRandomParticleOutsideLoopLabel:;

            if (FilledVoxelsForRandomParticle.empty() == false)
            {
                LoggersManagerObject.Log(STREAM("LocalNewParticleIndex = " << LocalNewParticleIndex));
                Particles[LocalNewParticleIndex].ListOfVoxels = FilledVoxelsForRandomParticle;
            }
        }

        for (auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
    }
    CATCH("generating random particles in selected space")
}

inline void GetRangeOfParticlesForRandomParticles(UniqueIdInt& StartParticleIndexParam, UniqueIdInt& EndParticleIndexParam, UniqueIdInt MaxParticleIndex)
{
    if (EndParticleIndexParam == 0)
    {
        StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
        EndParticleIndexParam = MaxParticleIndex;
    }
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusionForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        uniform_int_distribution<SignedInt> UniformDistributionObjectSizeOfParticle_int64t(-1, 1);

        vector<vector3_16> NewVoxelsForParticle;

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto &LocalParticleObject = Particles[ParticleIndex];

            SetAllVoxelsInListOfVoxelsToValue(LocalParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());

            SignedInt ShiftX = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftY = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftZ = UniformDistributionObjectSizeOfParticle_int64t(mt64R);

            NewVoxelsForParticle.clear();
            bool Collision = false;

            for (auto &VoxelForParticle: LocalParticleObject.ListOfVoxels)
                if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY,VoxelForParticle.Z + ShiftZ) == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                {
                    GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) = ParticleIndex;
                    NewVoxelsForParticle.emplace_back(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ);
                }
                else
                {
                    SetAllVoxelsInListOfVoxelsToValue(NewVoxelsForParticle, GetZeroSimulationSpaceVoxel());
                    SetAllVoxelsInListOfVoxelsToValue(LocalParticleObject.ListOfVoxels, ParticleIndex);
                    Collision = true;
                    break;
                }

            if (Collision == false)
                LocalParticleObject.ListOfVoxels = NewVoxelsForParticle;
        }
    }
    CATCH("generating one step of diffusion")
}

tuple<Particle*, UnsignedInt, string, vector<ChainIdInt>> GetNucleotidesSequenceForward(const UnsignedInt LengthOfSequence, Particle& ParticleObjectForReaction, const bool ToString, const bool ToVector)
{
    string NucleotidesSequenceToCompareString;
    vector<ChainIdInt> NucleotidesSequenceToCompareVector;
    Particle* ParticlePtr = &ParticleObjectForReaction;
    UnsignedInt NucleotidesCounter = 1;
    while (NucleotidesCounter < LengthOfSequence + 1 && ParticleObjectForReaction.Next != nullptr && ParticlePtr != nullptr)
    {
        if (ToString)
            NucleotidesSequenceToCompareString += CellEngineUseful::GetLetterForDNAChainId(ParticlePtr->ChainId);
        if (ToVector)
            NucleotidesSequenceToCompareVector.emplace_back(ParticlePtr->ChainId);
        ParticlePtr = ParticlePtr->Next;
        NucleotidesCounter++;
    }
    return { ParticlePtr, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector };
}

void CutDNA(Particle* NucleotideObjectForReactionPtr)
{
    NucleotideObjectForReactionPtr->Next->Prev = nullptr;
    NucleotideObjectForReactionPtr->Next = nullptr;
}

void LinkDNA(Particle* NucleotideObjectForReactionPtr1, Particle* NucleotideObjectForReactionPtr2)
{
    NucleotideObjectForReactionPtr1->Prev = NucleotideObjectForReactionPtr2;
    NucleotideObjectForReactionPtr2->Next = NucleotideObjectForReactionPtr1;
}

bool CellEngineVoxelSimulationSpace::CompareFitnessOfDNASequenceByString(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    return get<2>(GetNucleotidesSequenceForward(ParticleKindForReactionObject.SequenceStr.length(), ParticleObjectForReaction, true, false)) == ParticleKindForReactionObject.SequenceStr;
}

bool CellEngineVoxelSimulationSpace::CompareFitnessOfDNASequenceByNucleotidesLoop(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    bool FoundSequenceNotFit = false;

    try
    {
        auto [ParticlePtr, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector] = GetNucleotidesSequenceForward(ParticleKindForReactionObject.SequenceStr.length(), ParticleObjectForReaction, true, true);

        LoggersManagerObject.Log(STREAM("DNA SEQUENCE COMPARE = #" << NucleotidesSequenceToCompareString << "#" << ParticleKindForReactionObject.SequenceStr << "#" << to_string(ParticleKindForReactionObject.EntityId)));

                                                                                                                        if (NucleotidesSequenceToCompareVector.size() >= ParticleKindForReactionObject.Sequence.size())
                                                                                                                        {
        for (UnsignedInt NucleotideNum = 0; NucleotideNum < ParticleKindForReactionObject.Sequence.size(); NucleotideNum++)
            if (ParticleKindForReactionObject.Sequence[NucleotideNum] != NucleotidesSequenceToCompareVector[NucleotideNum])
            {
                if (NucleotidesSequenceToCompareString == ParticleKindForReactionObject.SequenceStr)
                    LoggersManagerObject.Log(STREAM("SEQUENCE DIFF = " << to_string(NucleotideNum) << " " << NucleotidesSequenceToCompareVector[NucleotideNum] << " " << ParticleKindForReactionObject.Sequence[NucleotideNum]));

                FoundSequenceNotFit = true;
                break;
            }
                                                                                                                        }
                                                                                                                        else
                                                                                                                            FoundSequenceNotFit = true;
        if (FoundSequenceNotFit == false)
            LoggersManagerObject.Log(STREAM("DNA SEQUENCE FOUND"));
    }
    CATCH("comparing fitness of dna sequence by nucleotides loop")

    return !FoundSequenceNotFit;
}

tuple<vector<UniqueIdInt>, bool> CellEngineVoxelSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(Reaction& ReactionObject)
{
    vector<pair<UniqueIdInt, UnsignedInt>> NucleotidesIndexesChosenForReaction;

    vector<UniqueIdInt> ParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : ParticlesSortedByCapacityFoundInParticlesProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.XCenter) << " Y = " << to_string(ParticleObjectTestedForReaction.YCenter) << " Z = " << to_string(ParticleObjectTestedForReaction.ZCenter)));

            vector<ParticleKindForReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [&ParticleObjectTestedForReaction](ParticleKindForReaction& ParticleKindObject){ return ParticleKindObject.EntityId == ParticleObjectTestedForReaction.EntityId; });
            else
            {
                ReactantIterator = find_if(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [&ParticleObjectTestedForReaction](ParticleKindForReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsSpecialDNA(ParticleKindForReactionObjectParam.EntityId) && CompareFitnessOfDNASequenceByNucleotidesLoop(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
                if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[ReactantIterator - ReactionObject.Reactants.begin()] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, ReactantIterator - ReactionObject.Reactants.begin());
            }

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[ReactantIterator - ReactionObject.Reactants.begin()] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[ReactantIterator - ReactionObject.Reactants.begin()] > 0)
            {
                LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.XCenter) << " Y = " << to_string(ParticleObjectTestedForReaction.YCenter) << " Z = " << to_string(ParticleObjectTestedForReaction.ZCenter) << endl));
                ReactantsCounters[ReactantIterator - ReactionObject.Reactants.begin()]--;
            }

            if (all_of(ReactantsCounters.begin(), ReactantsCounters.end(), [this](const UnsignedInt& Counter){ return Counter == 0; }) == true)
            {
                LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));
                break;
            }
            LoggersManagerObject.Log(STREAM(""));
        }

        if (ReactionObject.Id == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT 1"));

            if (NucleotidesIndexesChosenForReaction.size() == 1)
            {
                LoggersManagerObject.Log(STREAM("CUT 1 inside 1"));

                CutDNA(get<0>(GetNucleotidesSequenceForward(ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false)));
            }
            else
                return {};
        }
        else
        if (ReactionObject.Id == 2)
        {
            if (NucleotidesIndexesChosenForReaction.size() == 2)
            {
                LoggersManagerObject.Log(STREAM("LINK 1 inside 1"));
                                                                                                                        Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first);
                                                                                                                        Particle* NucleotideObjectForReactionPtr2 = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first);
                                                                                                                        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
                                                                                                                        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));

                LinkDNA(&GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), get<0>(GetNucleotidesSequenceForward(ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[1].second].SequenceStr.length(), GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first), false, false)));
            }
            else
                return {};
        }
//        else
//        if (ReactionObject.Id == 5)
//        {
//            if (NucleotidesIndexesChosenForReaction.size() == 1)
//            {
//                //PRZECINA OBIE NICI w JEDNYM MIEJSCU JAK ReactionObject.Id == 1 + AdditionalParameter
//            }
//            else
//            if (NucleotidesIndexesChosenForReaction.size() == 2)
//            {
//                //PRZECINA OBIE NICI:
//                //GORNA w MIEJSCU PO WYZNACZONYM NucleotidesIndexesChosenForReaction[0] + SequenceStr1 + AdditionalParameter
//                //A DOLNA w MIEJSCU NucleotidesIndexesChosenForReaction[0] + SequenceStr1 + SequenceStr2 + AdditionalParameter
//            }
//            else
//                return {};
//        }
    }
    CATCH("choosing particles for reaction from all particles in proximity")

    if (all_of(ReactantsCounters.begin(), ReactantsCounters.end(), [this](const UnsignedInt& Counter){ return Counter == 0; }) == true)
    {
        LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size())));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { vector<UniqueIdInt>(), false };
}

void CellEngineVoxelSimulationSpace::EraseParticlesChosenForReactionAndGetCentersForNewProductsOfReaction(const UnsignedInt ParticleIndexChosenForReaction, vector<vector3_16>& Centers)
{
    try
    {
        auto& ParticleObjectToBeErased = GetParticleFromIndex(ParticleIndexChosenForReaction);

        Centers.emplace_back(ParticleObjectToBeErased.XCenter, ParticleObjectToBeErased.YCenter, ParticleObjectToBeErased.ZCenter);
        LoggersManagerObject.Log(STREAM("Centers - X = " << to_string(ParticleObjectToBeErased.XCenter) << " Y = " << to_string(ParticleObjectToBeErased.YCenter) << " Z = " << to_string(ParticleObjectToBeErased.ZCenter) << endl));

        SetAllVoxelsInListOfVoxelsToValue(ParticleObjectToBeErased.ListOfVoxels, GetZeroSimulationSpaceVoxel());

        FreeIndexesOfParticles.push(ParticleIndexChosenForReaction);
        Particles.erase(ParticleIndexChosenForReaction);
    }
    CATCH("erasing particles chosen for reaction and get centers for new products of reaction")
}

bool CellEngineVoxelSimulationSpace::MakeChemicalReaction(Reaction& ReactionObject)
{
    try
    {
        auto [ParticlesIndexesChosenForReaction, FoundInProximity] = ChooseParticlesForReactionFromAllParticlesInProximity(ReactionObject);

        if (FoundInProximity == false)
            return false;

        LoggersManagerObject.Log(STREAM("Reaction Step 1 - chosen particles for reaction from all particles in proximity" << endl));

        vector<vector3_16> Centers;
        for (const auto& ParticleIndexChosenForReaction : ParticlesIndexesChosenForReaction)
            EraseParticlesChosenForReactionAndGetCentersForNewProductsOfReaction(ParticleIndexChosenForReaction, Centers);

        LoggersManagerObject.Log(STREAM("Reaction Step 2 - erasing particles chosen for reaction" << endl));

        LoggersManagerObject.Log(STREAM("Centers size = " << to_string(Centers.size()) << endl));

        UnsignedInt CenterIndex = 0;
        for (const auto& ReactionProduct : ReactionObject.Products)
        {
            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            if (CenterIndex < Centers.size())
                for (const auto& ParticleKindVoxel : ParticleKindObjectForProduct.ListOfVoxels)
                {
                    vector3_64 NewVoxel(Centers[CenterIndex].X - ParticleKindObjectForProduct.XSizeDiv2 + ParticleKindVoxel.X, Centers[CenterIndex].Y - ParticleKindObjectForProduct.YSizeDiv2 + ParticleKindVoxel.Y, Centers[CenterIndex].Z - ParticleKindObjectForProduct.ZSizeDiv2 + ParticleKindVoxel.Z);
                    GetSpaceVoxel(NewVoxel.X, NewVoxel.Y, NewVoxel.Z) = ParticleIndex;
                    GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(NewVoxel.X, NewVoxel.Y, NewVoxel.Z);
                    GetMinMaxCoordinatesForParticle(GetParticleFromIndex(ParticleIndex));
                    LoggersManagerObject.Log(STREAM("New Centers From Product Added X = " << to_string(NewVoxel.X) << " Y = " << to_string(NewVoxel.Y) << " Z = " << to_string(NewVoxel.Z) << endl));
                }
            else
            {
            }

            CenterIndex++;
        }

        LoggersManagerObject.Log(STREAM("Reaction Step 3 - Reaction finished" << endl));
    }
    CATCH("making reaction")

    return true;
};

std::vector<UnsignedInt> CellEngineVoxelSimulationSpace::GetRandomParticles(const UnsignedInt NumberOfReactants)
{
    vector<UnsignedInt> RandomParticlesTypes;

    try
    {
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKindsFoundInParticlesProximity.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKindsFoundInParticlesProximity), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);

            if (CellEngineUseful::IsSpecialDNA(RandomParticlesTypes.back()) == true)
                RandomParticlesTypes.back() = CellEngineConfigDataObject.DNAIdentifier;

            LoggersManagerObject.Log(STREAM("ParticleKind Reactant " << to_string(ReactantNumber) << " (" << to_string(RandomParticlesTypes.back()) << ")"));
        }
    }
    CATCH("getting random particles kind")

    return RandomParticlesTypes;
}

bool CellEngineVoxelSimulationSpace::IsChemicalReactionPossible(const Reaction& ReactionObject)
{
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInParticlesProximity[CellEngineUseful::IfSpecialDNAThenReturnNormalDNACode(ReactionReactant.EntityId)]; });
}

void CellEngineVoxelSimulationSpace::PrepareRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMainRandomCondition_Uint64t(0, 1);
        if (UniformDistributionObjectMainRandomCondition_Uint64t(mt64R) == 0)
            return;

        LoggersManagerObject.Log(STREAM("Looking for particles in proximity"));
    }
    CATCH("preparing random reaction")
}

void CellEngineVoxelSimulationSpace::FindAndExecuteRandomReaction()
{
    try
    {
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectNumberOfReactants_Uint64t(2, 2);

        UnsignedInt NumberOfTries = 0;
        while (NumberOfTries <= 100)
        {
            NumberOfTries++;

            UnsignedInt NumberOfReactants = UniformDistributionObjectNumberOfReactants_Uint64t(mt64R);

            if (TryToMakeRandomChemicalReaction(NumberOfReactants) == true)
                break;
        }
    }
    CATCH("finding and executing random reaction")
}

bool CellEngineVoxelSimulationSpace::FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        set<UnsignedInt> FoundParticleIndexes;

        ParticlesKindsFoundInParticlesProximity.clear();
        ParticlesSortedByCapacityFoundInParticlesProximity.clear();
        NucleotidesWithFreeEndingsFoundInParticlesProximity.clear();
        NucleotidesFreeFoundInParticlesProximity.clear();

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX >= 0 && PosY >= 0 && PosZ >= 0 && PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                        {
                            UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosY, PosZ);
                            if (FoundParticleIndexes.find(ParticleIndex) == FoundParticleIndexes.end())
                            {
                                ParticlesSortedByCapacityFoundInParticlesProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(ParticleIndex).EntityId) && (GetParticleFromIndex(ParticleIndex).Next == nullptr || GetParticleFromIndex(ParticleIndex).Prev == nullptr))
                                    NucleotidesWithFreeEndingsFoundInParticlesProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(ParticleIndex).EntityId) && GetParticleFromIndex(ParticleIndex).Next == nullptr && GetParticleFromIndex(ParticleIndex).Prev == nullptr)
                                    NucleotidesFreeFoundInParticlesProximity.emplace_back(ParticleIndex);
                                ParticlesKindsFoundInParticlesProximity[GetParticleFromIndex(ParticleIndex).EntityId]++;
                                FoundParticleIndexes.insert(ParticleIndex);
                            }
                        }

        if (ParticlesSortedByCapacityFoundInParticlesProximity.empty() == false)
        {
            sort(ParticlesSortedByCapacityFoundInParticlesProximity.begin(), ParticlesSortedByCapacityFoundInParticlesProximity.end(), [this](UnsignedInt PK1, UnsignedInt PK2) { return GetParticleFromIndex(PK1).ListOfVoxels.size() > GetParticleFromIndex(PK2).ListOfVoxels.size(); });

            LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
            for (const auto& LocalParticleIndexObjectToWrite : ParticlesSortedByCapacityFoundInParticlesProximity)
                LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId)));
            LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInParticlesProximity List"));
            for (const auto& LocalParticleKindObjectToWrite : ParticlesKindsFoundInParticlesProximity)
                LoggersManagerObject.Log(STREAM("ParticleKind EntityId = " << to_string(LocalParticleKindObjectToWrite.first) << " in quantity = " << to_string(LocalParticleKindObjectToWrite.second)));
        }
        else
        {
            LoggersManagerObject.Log(STREAM(endl << "No particle found in proximity "));
            return false;
        }
    }
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

bool CellEngineVoxelSimulationSpace::FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(const Particle& ParticleObject, const UnsignedInt AdditionalBoundFactor)
{
    try
    {
        LoggersManagerObject.Log(STREAM("EntityId = " << to_string(ParticleObject.EntityId)));

        auto ParticleKindObject = ParticlesKindsManagerObject.GetParticleKind(ParticleObject.EntityId);

        FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(ParticleObject.XCenter - ParticleKindObject.XSizeDiv2 - AdditionalBoundFactor, ParticleObject.YCenter - ParticleKindObject.YSizeDiv2 - AdditionalBoundFactor, ParticleObject.ZCenter - ParticleKindObject.ZSizeDiv2 - AdditionalBoundFactor, 2 * ParticleKindObject.XSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.YSizeDiv2 + 2 * AdditionalBoundFactor, 2 * ParticleKindObject.ZSizeDiv2 + 2 * AdditionalBoundFactor);
    }
    CATCH("finding particles in proximity of voxel simulation space for chosen particle")

    return true;
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForSelectedVoxelSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForParticle(Particle& ParticleObject)
{
    try
    {
        PrepareRandomReaction();

        if (FindParticlesInProximityOfVoxelSimulationSpaceForChosenParticle(ParticleObject, 20) == true)
            FindAndExecuteRandomReaction();
    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::FillSquareParticle(UnsignedInt RPosX, UnsignedInt RPosY, UnsignedInt RPosZ, UniqueIdInt LocalNewParticleIndex, UniqueIdInt Size)
{
    try
    {
        vector<vector3_16> FilledVoxelsForRandomParticle;

        for (UnsignedInt PosX = RPosX; PosX < RPosX + Size; PosX++)
            for (UnsignedInt PosY = RPosY; PosY < RPosY + Size; PosY++)
                for (UnsignedInt PosZ = RPosZ; PosZ < RPosZ + Size; PosZ++)
                {
                    GetSpaceVoxel(PosX, PosY, PosZ) = LocalNewParticleIndex;
                    FilledVoxelsForRandomParticle.emplace_back(PosX, PosY, PosZ);
                }

        Particles[LocalNewParticleIndex].ListOfVoxels = FilledVoxelsForRandomParticle;

        GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
    }
    CATCH("filling square particles")
}

void CellEngineVoxelSimulationSpace::GeneratePlanedParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        AddBasicParticlesKindsAndReactions();

        MakeZeroVoxelsForSelectedSpace(StartXPosParam, StartYPosParam, StartZPosParam, StepXParam, StepYParam, StepZParam, SizeXParam, SizeYParam, SizeZParam);

        auto P0 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 0, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 1, StartYPosParam + 3, StartZPosParam + 7, P0, 1);

        auto P1 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 1, 1,  -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 3, StartYPosParam + 13, StartZPosParam + 6, P1, 1);

        auto P2 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 2, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 6, StartYPosParam + 12, StartZPosParam + 5, P2, 1);

        auto P3 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 3, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 9, StartYPosParam + 9, StartZPosParam + 9, P3, 2);

        auto P4 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 4, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 17, StartYPosParam + 15, StartZPosParam + 15, P4, 2);

        auto P5 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 5, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 17, StartYPosParam + 19, StartZPosParam + 4, P5, 2);

        auto P6 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 21, StartYPosParam + 20, StartZPosParam + 19, P6, 2);

        auto P7 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 3, StartYPosParam + 20, StartZPosParam + 7, P7, 2);

        auto P8 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 5, StartYPosParam + 15, StartZPosParam + 20, P8, 2);

        auto P9 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 10, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 21, StartYPosParam + 5, StartZPosParam + 3, P9, 2);
    }
    CATCH("generating planed particles in selected space")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionsForAllParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.SelectedForReaction == false)
                GenerateRandomReactionForParticle(ParticleObject.second);
    }
    CATCH("generating random reactions for all particles")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto ParticlesIterator = Particles.find(ParticleIndex);
            if (ParticlesIterator != Particles.end())
                GenerateRandomReactionForParticle(ParticlesIterator->second);
        }
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForOneParticle(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, MaxParticleIndex);
        GenerateRandomReactionForParticle(GetParticleFromIndex(StartParticleIndexParam + 4));
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineVoxelSimulationSpace::EraseAllDNAParticles()
{
    try
    {
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
            {
                SetAllVoxelsInListOfVoxelsToValue(ParticleObject.second.ListOfVoxels, GetZeroSimulationSpaceVoxel());

                FreeIndexesOfParticles.push(ParticleObject.first);
            }

        const auto RemovedDNAParticlesCounter = erase_if(Particles, [](const pair<UniqueIdInt, Particle>& item) { auto const& [key, value] = item; return (value.EntityId == CellEngineConfigDataObject.DNAIdentifier); });
        LoggersManagerObject.Log(STREAM("RemovedDNAParticlesCounter = " << RemovedDNAParticlesCounter));

        UnsignedInt DNAParticleCounter = 0;
        for (auto& ParticleObject : Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                DNAParticleCounter++;
        LoggersManagerObject.Log(STREAM("DNAParticleCounter = " << DNAParticleCounter));

        Genomes[0].clear();
        Genomes[1].clear();
    }
    CATCH("erasing all dna particles")
}

void CellEngineVoxelSimulationSpace::UpdateRandomPositions(const UnsignedInt RandomMoveDirection, UnsignedInt& RandomPosX, UnsignedInt& RandomPosY, UnsignedInt& RandomPosZ, const UnsignedInt Size)
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

Particle* CellEngineVoxelSimulationSpace::GenerateDNAParticle(Particle* ParticlePrev, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeThread, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam)
{
    UnsignedInt ParticleIndex;

    try
    {
        ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), EntityId, ChainId, GenomeThread, GenomeIndex, UniqueColorParam));

        Particle::InsertAfterGivenNode(ParticlePrev, &GetParticleFromIndex(ParticleIndex));

        for (UnsignedInt PosX = StartPosX; PosX < StartPosX + ParticleSizeX; PosX++)
            for (UnsignedInt PosY = StartPosY; PosY < StartPosY + ParticleSizeY; PosY++)
                for (UnsignedInt PosZ = StartPosZ; PosZ < StartPosZ + ParticleSizeZ; PosZ++)
                {
                    GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
                    GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
                }

        Genome.emplace_back(ParticleIndex);
    }
    CATCH("generating particle")

    return &GetParticleFromIndex(ParticleIndex);
}

tuple<Particle*, Particle*> CellEngineVoxelSimulationSpace::GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam, const bool Linked)
{
    try
    {
        ParticlePrev1 = GenerateDNAParticle(ParticlePrev1, CellEngineConfigDataObject.DNAIdentifier, ChainId, 0, Genomes[0].size(), StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
        ParticlePrev2 = GenerateDNAParticle(ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainId(ChainId), 1, Genomes[1].size(), StartPosX + AddSizeX, StartPosY + AddSizeY, StartPosZ + AddSizeZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, Genomes[1], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
        if (Linked == true)
            PairedNucleotide<Particle>::LinkPairedNucleotides(ParticlePrev1, ParticlePrev2);
    }
    CATCH("generating two paired nucleotides")

    return { ParticlePrev1, ParticlePrev2 };
}

bool CellEngineVoxelSimulationSpace::TestFormerForbiddenPositions(unordered_set<string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt Size)
{
    UnsignedInt PosX = RandomPosX;
    UnsignedInt PosY = RandomPosY;
    UnsignedInt PosZ = RandomPosZ;

    UpdateRandomPositions(RandomMoveDirection, PosX, PosY, PosZ, Size);

    return (TestedFormerForbiddenPositions.find(to_string(PosX) + "|" + to_string(PosY) + "|" + to_string(PosZ)) != TestedFormerForbiddenPositions.end());
}

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineVoxelSimulationSpace::EraseLastRandomDNAParticle(vector<UniqueIdInt>& Genome)
{
    UnsignedInt LocalRandomPosX, LocalRandomPosY, LocalRandomPosZ;

    try
    {
        UnsignedInt PreviousParticleIndex = Genome.back();
        Genome.pop_back();

        LocalRandomPosX = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].X;
        LocalRandomPosY = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Y;
        LocalRandomPosZ = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Z;

        SetAllVoxelsInListOfVoxelsToValue(GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels, GetZeroSimulationSpaceVoxel());

        Particle::DeleteNode(nullptr, &GetParticleFromIndex(PreviousParticleIndex));

        Particles.erase(PreviousParticleIndex);

        FreeIndexesOfParticles.push(PreviousParticleIndex);

        LoggersManagerObject.Log(STREAM("ERASED PARTICLE PreviousParticleIndex = " << PreviousParticleIndex << " RandomPosX = " << LocalRandomPosX << " RandomPosY = " << LocalRandomPosY << " RandomPosZ = " << LocalRandomPosZ));
    }
    CATCH("erasing last random particle")

    return { LocalRandomPosX, LocalRandomPosY, LocalRandomPosZ };
}

UnsignedInt Sqr(UnsignedInt Value)
{
    return Value * Value;
}

void CellEngineVoxelSimulationSpace::GenerateRandomDNAInWholeCell(UnsignedInt NumberOfNucleotidesToBeGenerated, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, UnsignedInt ParticleSizeX, UnsignedInt ParticleSizeY, UnsignedInt ParticleSizeZ, UnsignedInt ParticleSize1, UnsignedInt ParticleSize2, UnsignedInt ParticleSize3, UnsignedInt ParticleSize4, UnsignedInt ParticleSize5)
{
    try
    {
        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = Particles.size();

        LoggersManagerObject.InitializePrintingParameters(false, false, false, false, false, false, false, false, false, false, false, false, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [](vector<UnsignedInt>& RandomMovesDirections) { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

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
                        tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, 1, ParticleSizeZ, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
                    else
                    if (RandomMoveDirection == 3 || RandomMoveDirection == 4 || RandomMoveDirection == 5 || RandomMoveDirection == 6)
                        tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, ChainId, Genomes[0].size(), RandomPosX, RandomPosY, RandomPosZ, 1, ParticleSizeY, ParticleSizeZ, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);

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

        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineVoxelSimulationSpace::GetMinMaxCoordinatesForDNA()
{
    try
    {
        for (auto& LocalNewParticleIndex : Genomes[0])
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
        for (auto& LocalNewParticleIndex : Genomes[1])
            GetMinMaxCoordinatesForParticle(Particles[LocalNewParticleIndex]);
    }
    CATCH("getting min max of coordinates for DNA")
}

void CellEngineVoxelSimulationSpace::SaveGenomeDataToFile(UnsignedInt ParticleSize)
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

void CellEngineVoxelSimulationSpace::ReadGenomeDataFromFile(bool Paired)
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
                    tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, ParticleSize, 1, ParticleSize, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
                else
                    tie(ParticlePrev1, ParticlePrev2) = GenerateTwoPairedNucleotides(ParticlePrev1, ParticlePrev2, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, 1, ParticleSize, ParticleSize, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()), true);
            }
            else
                ParticlePrev1 = GenerateDNAParticle(ParticlePrev1, stoi(Row[0]), stoi(Row[1]), 0, stoi(Row[2]), stoi(Row[3]), stoi(Row[4]), stoi(Row[5]), ParticleSize, ParticleSize, ParticleSize, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

            PrevStartPosX = StartPosX;

            GenomeIndex++;
        }

        FileToReadGenome.close();

        GetMinMaxCoordinatesForDNA();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("reading genome data from file")
}

void CellEngineVoxelSimulationSpace::ReadGenomeSequenceFromFile()
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_SEQUENCE.DAT");

        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome, GenomesLines[0]);

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < Genomes[0].size(); GenomeIndex++)
        {
            GetParticleFromIndex(Genomes[0][GenomeIndex]).ChainId = CellEngineUseful::GetChainIdFromLetter(GenomesLines[0][GenomeIndex]);
            GetParticleFromIndex(Genomes[1][GenomeIndex]).ChainId = CellEngineUseful::GetPairedChainId(CellEngineUseful::GetChainIdFromLetter(GenomesLines[0][GenomeIndex]));
        }
    }
    CATCH("reading real genome data from file")
}

void CellEngineVoxelSimulationSpace::TestGeneratedGenomeCorrectness(const UnsignedInt ParticleSize)
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
                    LoggersManagerObject.Log(STREAM("GenomeIndex = " << GenomeIndex << " ChainId = " << GetParticleFromIndex(Genomes[0][GenomeIndex]).ChainId << " Letter = " << CellEngineUseful::GetLetterForDNAChainId(GetParticleFromIndex(Genomes[0][GenomeIndex]).ChainId)));
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