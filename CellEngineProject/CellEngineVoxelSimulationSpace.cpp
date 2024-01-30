
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

#include "TerminalColorsUtils.h"

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
        ParticlesKindsManagerObject.AddParticleKind({ 12, "Polymerase", "POL", 0 });

        ParticlesKindsManagerObject.AddParticleKind({ CellEngineConfigDataObject.DNAIdentifier, "DNA", "DNA", 0 });

        ParticlesKindsManagerObject.AddParticleKind({ 10001, "DNA", "?", 0 });

        const string DNASequenceForTestFindingDNA = "TACAAAAAAAGAGGTGTTAGC";

        const string DNASequence1ForTestCutLink1 = "TACAAAAAAAGAGGTGTT";
        const string DNASequence2ForTestCutLink1 = "AGCTCTTATTA";

        const string DNASequence1ForTestCutLink1Any = "ANY";

        const string DNASequence1ForTestCutLink2 = "TACAAAAAAAGAGGTGTT";

        const string DNASequenceForTestCutCrisper = "RNA";
        const string RNASequenceForTestCutCrisper = "ANY";

        AddChemicalReaction(Reaction(1101, "STD ONLY WITH SEQ", "CH3CH2(OH) + DNA + ", { { 5, 1, "", true }, { 10001, 1, DNASequenceForTestFindingDNA, false } }, { { 10, 1, "", true } }, nullptr, "standard normal reaction example"));

        AddChemicalReaction(Reaction(10, "CUT 1 SEQ", "CH3CH2(OH) + DNA + ", { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink1, false } }, {}, &CellEngineVoxelSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "only in presence of chosen sequence of DNA"));

        AddChemicalReaction(Reaction(20, "LINK 1 SEQ", "CH3CH2(OH) + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink1, false }, { 10001, 1, DNASequence2ForTestCutLink1, false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNAInChosenPlaceSpecialReactionFunction, "links one strand of DNA with 2 endings with chosen sequences"));
        AddChemicalReaction(Reaction(30, "LINK 1 ANY", "CH3CH2(OH) + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 2, DNASequence1ForTestCutLink1Any, false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "link one strand of DNA with 2 endings with any sequence"));

        AddChemicalReaction(Reaction(40, "CUT 2 SEQ SHIFT 3 10", "CH3CHCH2 + DNA + ", 3, 10, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineVoxelSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=3 shift=10"));
        AddChemicalReaction(Reaction(41, "CUT 2 SEQ SHIFT 7 3", "CH3CHCH2 + DNA + ", 7, 3, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineVoxelSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=7 shift2=3"));
        AddChemicalReaction(Reaction(42, "CUT 2 SEQ SHIFT 10 3", "CH3CHCH2 + DNA + ", 10, 3, { { 5, 1, "", false }, { 10001, 1, DNASequence1ForTestCutLink2, false } }, { { 86, 1, "", true } }, &CellEngineVoxelSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction, "cuts both strands when found first sequence plus additional shift in both strands defined by additional parameters shift1=10 shift2=3"));

        AddChemicalReaction(Reaction(60, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGC", false }, { 10001, 1, "TCTTATT", false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));
        AddChemicalReaction(Reaction(61, "LINK 2 SEQ COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGCTCTT", false }, { 10001, 1, "ATTATGA", false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction, "Links both DNA strands if first strand has first sequence and second strand has second sequence and opposite joining strands and complementary"));

        AddChemicalReaction(Reaction(70, "LINK 2 ANY COMPLEMENT", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "ANY", false }, { 10001, 1, "ANY", false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction, "Links both DNA strands with any sequence if opposite joining strands when complementary"));

        AddChemicalReaction(Reaction(80, "LINK 2 ANY EQU SAME", "CH3CHCH2 + DNA + DNA + ", { { 5, 1, "", false }, { 10001, 1, "ANY", false }, { 10001, 1, "ANY", false } }, {}, &CellEngineVoxelSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction, "links both strands of DNA if they are cut equally in the same place"));

        AddChemicalReaction(Reaction(100, "CUT CRISPER 1", "CH3CHCH2 + RNA + DNA + ", 3, 7, { { 5, 1, "", false }, { 10001, 1, DNASequenceForTestCutCrisper, false }, { 10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineVoxelSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of one strand of DNA with RNA as template"));
        AddChemicalReaction(Reaction(110, "CUT CRISPER 2", "CH3CHCH2 + RNA + DNA + ", 3, 7, { { 5, 1, "", false }, { 10001, 1, DNASequenceForTestCutCrisper, false }, { 10001, 0, RNASequenceForTestCutCrisper, false } }, {}, &CellEngineVoxelSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction, "cut crisper of two strands of DNA with RNA as template"));

        AddChemicalReaction(Reaction(150, "POLYMERASE DNA START SEQ", "CH3CHCH2 + DNA + ", { { 5, 1, "", false }, { 10001, 1, "TACAAAAAAAGAGGTGTTAGC", false } }, {}, &CellEngineVoxelSimulationSpace::PolymeraseDNAStartSpecialReactionFunction, "links particle with DNA when found sequence and joins first nucleotide "));
        AddChemicalReaction(Reaction(160, "POLYMERASE DNA CONTINUE", "CH3CHCH2 + DNA + ", { { 5, 1, "", false, { CellEngineConfigDataObject.RNAIdentifier, CellEngineConfigDataObject.DNAIdentifier } } }, {}, &CellEngineVoxelSimulationSpace::PolymeraseDNAContinueSpecialReactionFunction, "links new nucleotide that fits nucleotide in DNA if free nucleotides are found in proximity"));

        AddChemicalReaction(Reaction(1001, "STD", "C6H12O6 + O2 + ", { { 1, 1, "", true }, { 2, 6, "", true } }, { { 3, 6, "", true }, { 0, 6, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1002, "STD", "CH2CH2 + H2O + ", { { 4, 1, "", true }, { 0, 1, "", true } }, { { 5, 1, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1003, "STD", "CH3CHCH2 + HX + ", { { 6, 1, "", true }, { 7, 1, "", true } }, { { 8, 1, "", true } }, nullptr));
        AddChemicalReaction(Reaction(1004, "STD", "CH2CH2 + O + ", { { 4,  1, "", true }, { 11, 1, "", true } }, { { 10, 1, "", true } }, nullptr));

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
                            //REDRAW MOVE
                            if (GetSpaceVoxel(PosX, PosY, PosZ) == GetZeroSimulationSpaceVoxel())
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

            //MOVE OBJECT - najpierw sprawdz czy wolny obiekt osobna funkcja!
            for (auto &VoxelForParticle: LocalParticleObject.ListOfVoxels)
                if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY,VoxelForParticle.Z + ShiftZ) == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                {
                    //REDRAW
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

double DistanceOfParticles(Particle& Particle1, Particle& Particle2)
{
    return sqrt(pow(Particle1.XCenter - Particle2.XCenter, 2.0) + pow(Particle1.YCenter - Particle2.YCenter, 2.0) + pow(Particle1.ZCenter - Particle2.ZCenter, 2.0));
}

string GetPairedSequenceStr(string SequenceStr)
{
    try
    {
        for (auto& SequenceNucleotideChar : SequenceStr)
            SequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(SequenceNucleotideChar)));
    }
    CATCH("getting paired sequence")

    return SequenceStr;
}

tuple<Particle*, Particle*, UnsignedInt, string, vector<ChainIdInt>> GetNucleotidesSequence(Particle* Particle::*Direction, const UnsignedInt LengthOfSequence, Particle& ParticleObjectForReaction, const bool ToString, const bool ToVector, bool (*Predicate)(const Particle*))
{
    string NucleotidesSequenceToCompareString;
    vector<ChainIdInt> NucleotidesSequenceToCompareVector;

    Particle* ParticlePtr = &ParticleObjectForReaction;
    Particle* ParticlePtrPrev = &ParticleObjectForReaction;
    UnsignedInt NucleotidesCounter = 1;

    try
    {
        while (NucleotidesCounter < LengthOfSequence + 1 && ParticleObjectForReaction.*Direction != nullptr && ParticlePtr != nullptr && Predicate(ParticlePtr) == true)
        {
            if (ToString == true)
                NucleotidesSequenceToCompareString += CellEngineUseful::GetLetterFromChainIdForDNAorRNA(ParticlePtr->ChainId);
            if (ToVector == true)
                NucleotidesSequenceToCompareVector.emplace_back(ParticlePtr->ChainId);
            ParticlePtrPrev = ParticlePtr;
            ParticlePtr = ParticlePtr->*Direction;
            NucleotidesCounter++;
        }
    }
    CATCH("getting nucleotides sequence forward")

    return { ParticlePtr, ParticlePtrPrev, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector };
}

void CutDNA(Particle* NucleotideObjectForReactionPtr)
{
    try
    {
        if (NucleotideObjectForReactionPtr->Next != nullptr)
        {
            NucleotideObjectForReactionPtr->Next->Prev = nullptr;
            NucleotideObjectForReactionPtr->Next = nullptr;
        }
    }
    CATCH("cutting DNA")
}

void LinkDNA(Particle* NucleotideObjectForReactionPtr1, Particle* NucleotideObjectForReactionPtr2)
{
    try
    {
        if (NucleotideObjectForReactionPtr1 != nullptr && NucleotideObjectForReactionPtr2 != nullptr)
        {
            NucleotideObjectForReactionPtr1->Prev = NucleotideObjectForReactionPtr2;
            NucleotideObjectForReactionPtr2->Next = NucleotideObjectForReactionPtr1;
        }
    }
    CATCH("linking DNA")
}

void SeparateTwoPairedDNANucleotides(Particle* ParticleNucleotide)
{
    try
    {
        LoggersManagerObject.Log(STREAM("SEPARATE PAIRED NUCLEOTIDE"));

        ParticleNucleotide->PairedNucleotide->PairedNucleotide = nullptr;
        ParticleNucleotide->PairedNucleotide = nullptr;
    }
    CATCH("separating dna strands")
}

void SeparateDNAStrands(Particle* Particle::*Direction, Particle* NucleotidePtr, const UnsignedInt LengthOfStrand)
{
    try
    {
        for (UnsignedInt DiffPos = 0; DiffPos < LengthOfStrand; DiffPos++)
        {
            SeparateTwoPairedDNANucleotides(NucleotidePtr);
            NucleotidePtr = NucleotidePtr->*Direction;
        }
    }
    CATCH("separating dna strands")
}

void JoinDNAStrands(Particle* Particle::*Direction, Particle* Strand1, Particle* Strand2)
{
    try
    {
        while (Strand1->PairedNucleotide == nullptr)
        {
            Strand1->PairedNucleotide = Strand2;
            Strand2->PairedNucleotide = Strand1;
            Strand1 = Strand1->*Direction;
            Strand2 = Strand2->*Direction;
        }
    }
    CATCH("joining dna strands")
}

bool CellEngineVoxelSimulationSpace::CutDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT 10 or 40 inside 1"));

            auto NucleotidePtr1 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }));

            if (NucleotidePtr1 != nullptr && NucleotidePtr1->Prev != nullptr)
            {
                CutDNA(NucleotidePtr1->Prev);

                if (ReactionObject.Id == 40 || ReactionObject.Id == 41 || ReactionObject.Id == 42)
                {
                    auto NucleotidePtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() + ReactionObject.AdditionalParameter2, *GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).PairedNucleotide, false, false, [](const Particle*){ return true; }));
                    if (NucleotidePtr2 != nullptr && NucleotidePtr2->Prev != nullptr)
                    {
                        CutDNA(NucleotidePtr2->Prev);

                        SeparateDNAStrands(&Particle::Next, ReactionObject.AdditionalParameter1 < ReactionObject.AdditionalParameter2 ? NucleotidePtr1 : NucleotidePtr2, std::abs(static_cast<SignedInt>(ReactionObject.AdditionalParameter2) - static_cast<SignedInt>(ReactionObject.AdditionalParameter1)));
                    }
                }
            }

            return true;
        }
    }
    CATCH("cutting dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::LinkDNAInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 2)
        {
            LoggersManagerObject.Log(STREAM("LINK 1 inside 1"));

            Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first);
            Particle* NucleotideObjectForReactionPtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[1].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first), false, false, [](const Particle*){ return true; }));

            Particle* NucleotideObjectForReactionPtr1Inv = &GetParticleFromIndex(NucleotidesIndexesChosenForReaction[1].first);
            Particle* NucleotideObjectForReactionPtr2Inv = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReaction[0].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }));

            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2INV GENOME INDEX = " << NucleotideObjectForReactionPtr2Inv->GenomeIndex));

            if (NucleotideObjectForReactionPtr1->Prev == nullptr && NucleotideObjectForReactionPtr2->Next == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
            else
            if (NucleotideObjectForReactionPtr2Inv->Prev == nullptr && NucleotideObjectForReactionPtr1Inv->Next == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr2, NucleotideObjectForReactionPtr1);
            else
            if (NucleotideObjectForReactionPtr2Inv->Next == nullptr && NucleotideObjectForReactionPtr1Inv->Prev == nullptr)
                LinkDNA(NucleotideObjectForReactionPtr1Inv, NucleotideObjectForReactionPtr2Inv);

            return true;
        }
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::LinkDNAInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesWithFreePrevEndingsFoundInProximity.empty() == false && NucleotidesWithFreeNextEndingsFoundInProximity.empty() == false)
        {
            LoggersManagerObject.Log(STREAM("LINK 1 or 2 ANY inside 1"));

            if (DistanceOfParticles(GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]), GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0])) <= 2.0)
            {
                LoggersManagerObject.Log(STREAM("LINK 1 or 2 ANY CLOSE ENOUGH"));

                LinkDNA(&GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]), &GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0]));

                if (ReactionObject.Id == 80)
                    LinkDNA(GetParticleFromIndex(NucleotidesWithFreePrevEndingsFoundInProximity[0]).PairedNucleotide, GetParticleFromIndex(NucleotidesWithFreeNextEndingsFoundInProximity[0]).PairedNucleotide);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::CutDNACrisperInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesIndexesChosenForReaction.size() == 1)
        {
            LoggersManagerObject.Log(STREAM("CUT CRISPER inside 1"));

            Particle* FirstStrand = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.AdditionalParameter1, GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first), false, false, [](const Particle*){ return true; }))->Prev;
            CutDNA(FirstStrand);
            if (ReactionObject.Id == 110)
                CutDNA(FirstStrand->PairedNucleotide);
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::LinkDNALigaseInAnyPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        if (NucleotidesWithFreePrevEndingsFoundInProximity.empty() == false && NucleotidesWithFreeNextEndingsFoundInProximity.empty() == false)
        {
            LoggersManagerObject.Log(STREAM("LINK 2 ANY COMPATIBLE inside 1"));

            if (DistanceOfParticles(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]), GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0])) <= 2.0 && (DistanceOfParticles(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]), GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1])) <= 2.0))
            {
                LoggersManagerObject.Log(STREAM("LINK 2 ANY COMPATIBLE CLOSE ENOUGH - " << to_string(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).GenomeIndex) << " " << to_string(GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]).GenomeIndex) << " "));

                LinkDNA(&GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[0]), &GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[0]));
                LinkDNA(&GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]), &GetParticleFromIndex(DNANucleotidesWithFreeNextEndingsFoundInProximity[1]));

                JoinDNAStrands(&Particle::Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).Prev, GetParticleFromIndex(DNANucleotidesWithFreePrevEndingsFoundInProximity[1]).PairedNucleotide->Prev);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::LinkDNALigaseInChosenPlaceSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("LINK 2 REACTION 60"));

        auto NucleotidesIndexesChosenForReactionCopy(NucleotidesIndexesChosenForReaction);

        Particle* NucleotideObjectForReactionPtr1 = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first);
        Particle* NucleotideObjectForReactionPtr1Paired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[0].first).PairedNucleotide;
        Particle* NucleotideObjectForReactionPtr1Inv = &GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first);
        Particle* NucleotideObjectForReactionPtr1InvPaired = GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotide;

        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));

        if ((NucleotideObjectForReactionPtr1Inv->Next == nullptr || NucleotideObjectForReactionPtr1Inv->Prev == nullptr) && (NucleotideObjectForReactionPtr1->Next != nullptr && NucleotideObjectForReactionPtr1->Prev != nullptr))
        {
            swap(NucleotidesIndexesChosenForReactionCopy[0], NucleotidesIndexesChosenForReactionCopy[1]);

            swap(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr1Inv);
            swap(NucleotideObjectForReactionPtr1Paired, NucleotideObjectForReactionPtr1InvPaired);

            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 AFTER SWAP GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
            LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1INV AFTER SWAP GENOME INDEX = " << NucleotideObjectForReactionPtr1Inv->GenomeIndex));
        }

        Particle* NucleotideObjectForReactionPtr2 = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReactionCopy[1].second].SequenceStr.length() - 1, GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first), false, false, [](const Particle*){ return true; }));
        Particle* NucleotideObjectForReactionPtr2Paired = get<0>(GetNucleotidesSequence(&Particle::Next, ReactionObject.Reactants[NucleotidesIndexesChosenForReactionCopy[1].second].SequenceStr.length() - 1, *GetParticleFromIndex(NucleotidesIndexesChosenForReactionCopy[1].first).PairedNucleotide, false, false, [](const Particle*){ return true; }));

        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 1 GENOME INDEX = " << NucleotideObjectForReactionPtr1->GenomeIndex));
        LoggersManagerObject.Log(STREAM("NUCLEOTIDE 2 GENOME INDEX = " << NucleotideObjectForReactionPtr2->GenomeIndex));

        if (NucleotideObjectForReactionPtr1Paired == nullptr && NucleotideObjectForReactionPtr2Paired != nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X1"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Next, 32, *NucleotideObjectForReactionPtr1, true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Next, 32, *(NucleotideObjectForReactionPtr2Paired->Next), true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtr1->GenomeIndex << " " << ParticlePtr1->PairedNucleotide->GenomeIndex << " " << ParticlePtrPrev2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtr1->PairedNucleotide, ParticlePtrPrev2);
                JoinDNAStrands(&Particle::Next, NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2Paired->Next);
            }
        }
        else
        if (NucleotideObjectForReactionPtr1Paired != nullptr && NucleotideObjectForReactionPtr2Paired == nullptr)
        {
            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY C2X2"));

            auto [ParticlePtr1, ParticlePtrPrev1, Counter1, SequenceStr1, SequenceVector1] = GetNucleotidesSequence(&Particle::Prev, 32, *NucleotideObjectForReactionPtr1Paired->Prev, true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });
            auto [ParticlePtr2, ParticlePtrPrev2, Counter2, SequenceStr2, SequenceVector2] = GetNucleotidesSequence(&Particle::Prev, 32, *(NucleotideObjectForReactionPtr2), true, false, [](const Particle* P){ return P->PairedNucleotide == nullptr; });

            LoggersManagerObject.Log(STREAM("CHECKING COMPLEMENTARY Sequence1 = [" << SequenceStr1 << "] Sequence2 = [" << SequenceStr2 << "]"));

            if (SequenceStr1 == GetPairedSequenceStr(SequenceStr2))
            {
                LoggersManagerObject.Log(STREAM("LINKING PAIRED " << ParticlePtrPrev1->GenomeIndex << " " << ParticlePtr2->GenomeIndex));

                LinkDNA(NucleotideObjectForReactionPtr1, NucleotideObjectForReactionPtr2);
                LinkDNA(ParticlePtrPrev1, ParticlePtr2->PairedNucleotide);
                JoinDNAStrands(&Particle::Prev, NucleotideObjectForReactionPtr1Paired->Prev, NucleotideObjectForReactionPtr2);
            }
        }
        else
            return true;
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

void CellEngineVoxelSimulationSpace::MoveParticleByVector(Particle& ParticleObject, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ)
{
    try
    {
        SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, GetZeroSimulationSpaceVoxel());
        for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
        {
            VoxelOfParticle.X += VectorX;
            VoxelOfParticle.Y += VectorY;
            VoxelOfParticle.Z += VectorZ;
        }
        SetAllVoxelsInListOfVoxelsToValue(ParticleObject.ListOfVoxels, ParticleObject.Index);
    }
    CATCH("moving particle by vector")
}

void CellEngineVoxelSimulationSpace::MoveParticleNearOtherParticle(Particle& ParticleObject, const Particle& NewPositionParticleObject, const SignedInt AddX,  const SignedInt AddY,  const SignedInt AddZ)
{
    try
    {
        MoveParticleByVector(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
    }
    CATCH("moving particle near other particles")
}

bool CellEngineVoxelSimulationSpace::CheckFreeSpaceForParticleMovedByVector(Particle& ParticleObject, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ)
{
    try
    {
        for (auto &VoxelOfParticle: ParticleObject.ListOfVoxels)
            if (GetSpaceVoxel(VoxelOfParticle.X + VectorX, VoxelOfParticle.Y + VectorY, VoxelOfParticle.Z + VectorZ) != GetZeroSimulationSpaceVoxel())
                return false;

        return true;
    }
    CATCH("checking free space for particle moved by vector")
}

void CellEngineVoxelSimulationSpace::MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(Particle& ParticleObject, const Particle& NewPositionParticleObject, const SignedInt AddX,  const SignedInt AddY,  const SignedInt AddZ)
{
    try
    {
        MoveParticleByVectorIfVoxelSpaceIsEmpty(ParticleObject, NewPositionParticleObject.ListOfVoxels[0].X - ParticleObject.ListOfVoxels[0].X + AddX, NewPositionParticleObject.ListOfVoxels[0].Y - ParticleObject.ListOfVoxels[0].Y + AddY, NewPositionParticleObject.ListOfVoxels[0].Z - ParticleObject.ListOfVoxels[0].Z + AddZ);
    }
    CATCH("moving particle near other particles")
}

void CellEngineVoxelSimulationSpace::MoveParticleByVectorIfVoxelSpaceIsEmpty(Particle& ParticleObject, const UnsignedInt VectorX, const UnsignedInt VectorY, const UnsignedInt VectorZ)
{
    try
    {
        if (CheckFreeSpaceForParticleMovedByVector(ParticleObject, VectorX, VectorY, VectorZ) == true)
            MoveParticleByVector(ParticleObject, VectorX, VectorY, VectorZ);
    }
    CATCH("moving particle by vector if voxel space is empty")
}

bool CellEngineVoxelSimulationSpace::PolymeraseDNAStartSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("POLYMERASE DNA START REACTION"));

        if (NucleotidesFreeFoundInProximity.empty() == false && GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList.empty() == true)
        {
            LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(ParticlesIndexesChosenForReaction[0].first) << " Nucleotide = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).ChainId) << " Nucleotide Index = " << GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first).GenomeIndex));

            GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).AddNewLinkToParticle(&GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]));
            GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).AddNewLinkToParticle(&GetParticleFromIndex(NucleotidesIndexesChosenForReaction[0].first));

            MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first), *GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1], 4, 0, 0);
            MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]), GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first), 2, 0, 0);
        }
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::PolymeraseDNAContinueSpecialReactionFunction(const std::vector<std::pair<UniqueIdInt, UnsignedInt>>& ParticlesIndexesChosenForReaction, const vector<pair<UniqueIdInt, UnsignedInt>>& NucleotidesIndexesChosenForReaction, const Reaction& ReactionObject)
{
    try
    {
        LoggersManagerObject.Log(STREAM("POLYMERASE DNA CONTINUE REACTION"));

        LoggersManagerObject.Log(STREAM("Letter = " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]).ChainId) << " " << CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1]->ChainId))));

        if (GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]).ChainId == CellEngineUseful::GetPairedChainIdForDNAorRNA(GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1]->ChainId))
        {
            GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[0]->Next = &GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]);
            GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1] = GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1]->Next;

            //ZA DRUGIM CZYLI PIERWSZYM TUTAJ ZNAJDUJE TEN SAM NUKLEOTYD WIEC GO PRZESUWA A PRZECIEZ ON JUZ NIE WOLNY BO PRZYSPAWANY DO BIALKA
            MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(GetParticleFromIndex(NucleotidesFreeFoundInProximity[0]), GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first), 2, 0, 0);
            MoveParticleNearOtherParticleIfVoxelSpaceIsEmpty(GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first), *GetParticleFromIndex(ParticlesIndexesChosenForReaction[0].first).LinkedParticlesPointersList[1], 4, 0, 0);
        }
    }
    CATCH("linking dna in chosen place in special reaction function")

    return false;
}

bool CellEngineVoxelSimulationSpace::CompareFitnessOfParticle(const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    return (ParticleKindForReactionObject.LinkedParticleTypes.empty() == true || (ParticleKindForReactionObject.LinkedParticleTypes.empty() == false && ParticleObjectForReaction.LinkedParticlesPointersList.size() == ParticleKindForReactionObject.LinkedParticleTypes.size() && all_of(ParticleObjectForReaction.LinkedParticlesPointersList.begin(), ParticleObjectForReaction.LinkedParticlesPointersList.end(), [](const Particle* PointerToParticle){ return PointerToParticle != nullptr; }) && equal(ParticleObjectForReaction.LinkedParticlesPointersList.begin(), ParticleObjectForReaction.LinkedParticlesPointersList.end(), ParticleKindForReactionObject.LinkedParticleTypes.begin(), [](const Particle* PointerToParticle, UniqueIdInt ParticleType){ return PointerToParticle->EntityId == ParticleType; })));
}

bool CellEngineVoxelSimulationSpace::CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType TypeOfComparison, const ParticleKindForReaction& ParticleKindForReactionObject, Particle& ParticleObjectForReaction)
{
    bool FoundSequenceNotFit = false;

    try
    {
        string TemplateSequenceStr = ParticleKindForReactionObject.SequenceStr;

        if (TemplateSequenceStr == "ANY")
        {
            LoggersManagerObject.Log(STREAM("DNA SEQUENCE FOUND = ANY => ParticleKindForReactionObject.EntityId = " << to_string(ParticleKindForReactionObject.EntityId)));
            return true;
        }

        vector<ChainIdInt> TemplateSequence = ParticleKindForReactionObject.Sequence;

        string OriginalTemplateRNASequenceStr;

        if (TemplateSequenceStr == "RNA")
            if (RNANucleotidesFoundInProximity.empty() == false)
            {
                UnsignedInt LengthOfTemplateForRNA =  32;

                auto [ParticlePtrNext, ParticlePtrPrevNext, NucleotidesCounterNext, NucleotidesSequenceToCompareStringNext, NucleotidesSequenceToCompareVectorNext] = GetNucleotidesSequence(&Particle::Next, LengthOfTemplateForRNA, GetParticleFromIndex(RNANucleotidesFoundInProximity[0]), true, true, [](const Particle*){ return true; });
                auto [ParticlePtrBack, ParticlePtrBackNext, NucleotidesCounterBack, NucleotidesSequenceToCompareStringBack, NucleotidesSequenceToCompareVectorBack] = GetNucleotidesSequence(&Particle::Prev, LengthOfTemplateForRNA, GetParticleFromIndex(RNANucleotidesFoundInProximity[0]), true, true, [](const Particle*){ return true; });

                reverse(NucleotidesSequenceToCompareVectorBack.begin(), NucleotidesSequenceToCompareVectorBack.end());
                TemplateSequence = NucleotidesSequenceToCompareVectorBack;
                TemplateSequence.insert(end(TemplateSequence), begin(NucleotidesSequenceToCompareVectorNext), end(NucleotidesSequenceToCompareVectorNext));
                for (auto& Nucleotide : TemplateSequence)
                    Nucleotide = CellEngineUseful::GetPairedChainIdForDNAorRNA(Nucleotide);

                reverse(NucleotidesSequenceToCompareStringBack.begin(), NucleotidesSequenceToCompareStringBack.end());
                TemplateSequenceStr = NucleotidesSequenceToCompareStringBack + NucleotidesSequenceToCompareStringNext;
                OriginalTemplateRNASequenceStr = TemplateSequenceStr;
                for (auto& TemplateSequenceNucleotideChar : TemplateSequenceStr)
                    TemplateSequenceNucleotideChar = CellEngineUseful::GetLetterFromChainIdForDNAorRNA(CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(TemplateSequenceNucleotideChar)));
            }

        auto [ParticlePtr, ParticlePtrPrev, NucleotidesCounter, NucleotidesSequenceToCompareString, NucleotidesSequenceToCompareVector] = GetNucleotidesSequence(&Particle::Next, TemplateSequenceStr.length(), ParticleObjectForReaction, true, true, [](const Particle*){ return true; });

        LoggersManagerObject.Log(STREAM("DNA SEQUENCE COMPARE = #" << NucleotidesSequenceToCompareString << "#" << TemplateSequenceStr << "#" << OriginalTemplateRNASequenceStr << "#" << to_string(ParticleKindForReactionObject.EntityId)));

        if (TypeOfComparison == ComparisonType::ByVectorLoop)
        {
            if (NucleotidesSequenceToCompareVector.size() >= TemplateSequence.size())
            {
                LoggersManagerObject.Log(STREAM("LOOP COMPARISON SIZE = " << to_string(NucleotidesSequenceToCompareVector.size()) << " " << to_string(TemplateSequence.size())));

                for (UnsignedInt NucleotideNum = 0; NucleotideNum < TemplateSequence.size(); NucleotideNum++)
                    if (TemplateSequence[NucleotideNum] != NucleotidesSequenceToCompareVector[NucleotideNum])
                    {
                        LoggersManagerObject.Log(STREAM("LOOP COMPARISON BREAK = " << to_string(NucleotideNum) << "#"));

                        FoundSequenceNotFit = true;
                        break;
                    }
            }
            else
                FoundSequenceNotFit = true;
        }
        else
        if (TypeOfComparison == ComparisonType::ByString)
            FoundSequenceNotFit = !(NucleotidesSequenceToCompareString == TemplateSequenceStr);

        if (FoundSequenceNotFit == false)
            LoggersManagerObject.Log(STREAM(terminal_colors_utils::green << "DNA SEQUENCE FOUND FIT" << terminal_colors_utils::white));
    }
    CATCH("comparing fitness of dna sequence by nucleotides loop")

    return !FoundSequenceNotFit;
}

tuple<vector<pair<UniqueIdInt, UnsignedInt>>, bool> CellEngineVoxelSimulationSpace::ChooseParticlesForReactionFromAllParticlesInProximity(const Reaction& ReactionObject)
{
    bool AllAreZero = false;

    vector<pair<UniqueIdInt, UnsignedInt>> NucleotidesIndexesChosenForReaction, ParticlesIndexesChosenForReaction, AllParticlesIndexesChosenForReaction;

    vector<UnsignedInt> ReactantsCounters(ReactionObject.Reactants.size());

    try
    {
        for (UnsignedInt ReactantIndex = 0; ReactantIndex < ReactionObject.Reactants.size(); ReactantIndex++)
            ReactantsCounters[ReactantIndex] = ReactionObject.Reactants[ReactantIndex].Counter;

        for (const auto& ParticleObjectIndex : ParticlesSortedByCapacityFoundInProximity)
        {
            auto& ParticleObjectTestedForReaction = GetParticleFromIndex(ParticleObjectIndex);

            LoggersManagerObject.Log(STREAM("ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.XCenter) << " Y = " << to_string(ParticleObjectTestedForReaction.YCenter) << " Z = " << to_string(ParticleObjectTestedForReaction.ZCenter)));

            vector<ParticleKindForReaction>::const_iterator ReactantIterator;
            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == false)
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return ParticleKindForReactionObjectParam.EntityId == ParticleObjectTestedForReaction.EntityId && CompareFitnessOfParticle(ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });
            else
                ReactantIterator = find_if(ReactionObject.Reactants.cbegin(), ReactionObject.Reactants.cend(), [&ParticleObjectTestedForReaction, this](const ParticleKindForReaction& ParticleKindForReactionObjectParam){ return CellEngineUseful::IsSpecialDNA(ParticleKindForReactionObjectParam.EntityId) && CompareFitnessOfDNASequenceByNucleotidesLoop(ComparisonType::ByVectorLoop, ParticleKindForReactionObjectParam, ParticleObjectTestedForReaction) == true; });

            auto PositionInReactants = ReactantIterator - ReactionObject.Reactants.begin();

            if (CellEngineUseful::IsDNAorRNA(ParticleObjectTestedForReaction.EntityId) == true)
                if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == false)
                    NucleotidesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0 && ReactantIterator->ToRemoveInReaction == true)
                ParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);

            if (ReactantIterator != ReactionObject.Reactants.end() && ReactantsCounters[PositionInReactants] > 0)
            {
                AllParticlesIndexesChosenForReaction.emplace_back(ParticleObjectIndex, PositionInReactants);
                LoggersManagerObject.Log(STREAM("CHOSEN ParticleObjectIndex = " << to_string(ParticleObjectIndex) <<" EntityId = " << to_string(ParticleObjectTestedForReaction.EntityId) << " X = " << to_string(ParticleObjectTestedForReaction.XCenter) << " Y = " << to_string(ParticleObjectTestedForReaction.YCenter) << " Z = " << to_string(ParticleObjectTestedForReaction.ZCenter) << endl));
                ReactantsCounters[PositionInReactants]--;
            }

            AllAreZero = all_of(ReactantsCounters.begin(), ReactantsCounters.end(), [this](const UnsignedInt& Counter){ return Counter == 0; });
            if (AllAreZero == true)
            {
                LoggersManagerObject.Log(STREAM("ALL ARE ZERO"));
                break;
            }
            LoggersManagerObject.Log(STREAM(""));
        }

        if (AllAreZero == true || (AllAreZero == false && (ReactionObject.Id == 30 || ReactionObject.Id == 80 || ReactionObject.Id == 70)))
            ReactionObject.SpecialReactionFunction(this, AllParticlesIndexesChosenForReaction, NucleotidesIndexesChosenForReaction, ReactionObject);
    }
    CATCH("choosing particles for reaction from all particles in proximity")

    if (AllAreZero == true)
    {
        LoggersManagerObject.Log(STREAM("ALL ARE ZERO AT END = " << to_string(ParticlesIndexesChosenForReaction.size())));
        return { ParticlesIndexesChosenForReaction, true };
    }
    else
        return { vector<pair<UniqueIdInt, UnsignedInt>>(), false };
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
            EraseParticlesChosenForReactionAndGetCentersForNewProductsOfReaction(ParticleIndexChosenForReaction.first, Centers);

        LoggersManagerObject.Log(STREAM("Reaction Step 2 - erasing particles chosen for reaction" << endl));

        LoggersManagerObject.Log(STREAM("Centers size = " << to_string(Centers.size()) << endl));

        UnsignedInt CenterIndex = 0;
        for (const auto& ReactionProduct : ReactionObject.Products)
        {
            UnsignedInt ParticleIndex = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), ReactionProduct.EntityId, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.clear();

            auto& ParticleKindObjectForProduct = ParticlesKindsManagerObject.GetParticleKind(ReactionProduct.EntityId);

            //REDRAW
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
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectUint64t(0, ParticlesKindsFoundInProximity.size() - 1);

        for (UnsignedInt ReactantNumber = 1; ReactantNumber <= NumberOfReactants; ReactantNumber++)
        {
            RandomParticlesTypes.emplace_back(std::next(std::begin(ParticlesKindsFoundInProximity), static_cast<int>(UniformDistributionObjectUint64t(mt64R)))->first);

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
    return all_of(ReactionObject.Reactants.begin(), ReactionObject.Reactants.end(), [this](const ParticleKindForReaction& ReactionReactant){ return ReactionReactant.Counter <= ParticlesKindsFoundInProximity[CellEngineUseful::IfSpecialDNAThenReturnNormalDNACode(ReactionReactant.EntityId)]; });
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

void CellEngineVoxelSimulationSpace::FindAndExecuteChosenReaction(const UnsignedInt ReactionId)
{
    try
    {
        auto ReactionIter = ReactionsPosFromId.find(ReactionId);
        if (ReactionIter != ReactionsPosFromId.end())
        {
            auto& ReactionObject = Reactions[ReactionIter->second];

            bool IsPossible = IsChemicalReactionPossible(ReactionObject);
            if (IsPossible == true)
            {
                LoggersManagerObject.Log(STREAM("CHOSEN REACTION POSSIBLE" << endl));

                if (MakeChemicalReaction(ReactionObject) == false)
                    LoggersManagerObject.Log(STREAM("Chosen reaction not executed!"));
            }
            else
                LoggersManagerObject.Log(STREAM("Chosen reaction impossible!"));
        }
        else
            LoggersManagerObject.Log(STREAM("Chosen reaction Id not found!"));
    }
    CATCH("finding and executing chosen reaction")
}

bool CellEngineVoxelSimulationSpace::FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        set<UnsignedInt> FoundParticleIndexes;

        ParticlesKindsFoundInProximity.clear();
        ParticlesSortedByCapacityFoundInProximity.clear();
        NucleotidesWithFreeNextEndingsFoundInProximity.clear();
        NucleotidesWithFreePrevEndingsFoundInProximity.clear();
        DNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        DNANucleotidesWithFreePrevEndingsFoundInProximity.clear();
        RNANucleotidesWithFreeNextEndingsFoundInProximity.clear();
        RNANucleotidesWithFreePrevEndingsFoundInProximity.clear();
        NucleotidesFreeFoundInProximity.clear();
        RNANucleotidesFoundInProximity.clear();

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX++)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY++)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ++)
                    if (PosX >= 0 && PosY >= 0 && PosZ >= 0 && PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension && PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension)
                        if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                        {
                            UniqueIdInt ParticleIndex = GetSpaceVoxel(PosX, PosY, PosZ);
                            if (FoundParticleIndexes.find(ParticleIndex) == FoundParticleIndexes.end())
                            {
                                ParticlesSortedByCapacityFoundInProximity.emplace_back(ParticleIndex);
                                Particle& ParticleFromIndex = GetParticleFromIndex(ParticleIndex);

                                if (CellEngineUseful::IsDNAorRNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Next == nullptr && ParticleFromIndex.Prev != nullptr)
                                    NucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsDNAorRNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Prev == nullptr && ParticleFromIndex.Next != nullptr)
                                    NucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

                                if (CellEngineUseful::IsDNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Next == nullptr && ParticleFromIndex.Prev != nullptr)
                                    DNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsDNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Prev == nullptr && ParticleFromIndex.Next != nullptr)
                                    DNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

                                if (CellEngineUseful::IsRNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Next == nullptr && ParticleFromIndex.Prev != nullptr)
                                    RNANucleotidesWithFreeNextEndingsFoundInProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsRNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Prev == nullptr && ParticleFromIndex.Next != nullptr)
                                    RNANucleotidesWithFreePrevEndingsFoundInProximity.emplace_back(ParticleIndex);

                                if (CellEngineUseful::IsDNAorRNA(ParticleFromIndex.EntityId) && ParticleFromIndex.Next == nullptr && ParticleFromIndex.Prev == nullptr)
                                    NucleotidesFreeFoundInProximity.emplace_back(ParticleIndex);
                                if (CellEngineUseful::IsRNA(ParticleFromIndex.EntityId))
                                    RNANucleotidesFoundInProximity.emplace_back(ParticleIndex);

                                ParticlesKindsFoundInProximity[ParticleFromIndex.EntityId]++;
                                FoundParticleIndexes.insert(ParticleIndex);
                            }
                        }

        if (ParticlesSortedByCapacityFoundInProximity.empty() == false)
        {
            sort(ParticlesSortedByCapacityFoundInProximity.begin(), ParticlesSortedByCapacityFoundInProximity.end(), [this](UnsignedInt PK1, UnsignedInt PK2) { return GetParticleFromIndex(PK1).ListOfVoxels.size() > GetParticleFromIndex(PK2).ListOfVoxels.size(); });

            LoggersManagerObject.Log(STREAM(endl << "ParticlesSortedByCapacityFoundInParticlesProximity List"));
            for (const auto& LocalParticleIndexObjectToWrite : ParticlesSortedByCapacityFoundInProximity)
                LoggersManagerObject.Log(STREAM("ParticleIndex = " << to_string(LocalParticleIndexObjectToWrite) << " EntityId = " << to_string(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) << " NUCLEOTIDE = " << ((CellEngineUseful::IsDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).EntityId) == true) ? CellEngineUseful::GetLetterFromChainIdForDNAorRNA(GetParticleFromIndex(LocalParticleIndexObjectToWrite).ChainId) : '0') << " GENOME INDEX = " << GetParticleFromIndex(LocalParticleIndexObjectToWrite).GenomeIndex));
            LoggersManagerObject.Log(STREAM(endl << "ParticlesKindsFoundInProximity List"));
            for (const auto& LocalParticleKindObjectToWrite : ParticlesKindsFoundInProximity)
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
    CATCH("generating random reaction for selected voxel space")
}

void CellEngineVoxelSimulationSpace::GenerateChosenReactionForSelectedVoxelSpace(const UnsignedInt ReactionId, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, const UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        if (FindParticlesInProximityOfVoxelSimulationSpaceForSelectedVoxelSpace(StartXPosParam, StartYPosParam, StartZPosParam, SizeXParam, SizeYParam, SizeZParam) == true)
            FindAndExecuteChosenReaction(ReactionId);
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

void CellEngineVoxelSimulationSpace::FillSquareParticle(UnsignedInt RPosX, UnsignedInt RPosY, UnsignedInt RPosZ, UniqueIdInt LocalNewParticleIndex, UniqueIdInt SizeX, UniqueIdInt SizeY, UniqueIdInt SizeZ)
{
    try
    {
        vector<vector3_16> FilledVoxelsForRandomParticle;

        //REDRAW
        for (UnsignedInt PosX = RPosX; PosX < RPosX + SizeX; PosX++)
            for (UnsignedInt PosY = RPosY; PosY < RPosY + SizeY; PosY++)
                for (UnsignedInt PosZ = RPosZ; PosZ < RPosZ + SizeZ; PosZ++)
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
        FillSquareParticle(StartXPosParam + 1, StartYPosParam + 3, StartZPosParam + 7, P0, 1, 1, 1);

        auto P1 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 1, 1,  -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 3, StartYPosParam + 13, StartZPosParam + 6, P1, 1, 1, 1);

        auto P2 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 2, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 6, StartYPosParam + 12, StartZPosParam + 5, P2, 1, 1, 1);

        auto P3 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 3, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 9, StartYPosParam + 9, StartZPosParam + 9, P3, 2, 2, 2);

        auto P4 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 4, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 17, StartYPosParam + 15, StartZPosParam + 15, P4, 2, 2, 2);

        auto P5 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 5, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 17, StartYPosParam + 19, StartZPosParam + 4, P5, 2, 2, 2);

        auto P6 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 21, StartYPosParam + 20, StartZPosParam + 19, P6, 2, 2, 2);

        auto P7 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 3, StartYPosParam + 20, StartZPosParam + 7, P7, 2, 2, 2);

        auto P8 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 6, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 5, StartYPosParam + 15, StartZPosParam + 20, P8, 2, 2, 2);

        auto P9 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), 10, 1, -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 21, StartYPosParam + 5, StartZPosParam + 3, P9, 2, 2, 2);

        auto P10 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('T'), -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 9, StartYPosParam + 12, StartZPosParam + 3, P10, 2, 2, 1);

        auto P11 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('A'), -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 9, StartYPosParam + 16, StartZPosParam + 3, P11, 2, 2, 1);

        auto P12 = AddNewParticle(Particle(GetNewFreeIndexOfParticle(), CellEngineConfigDataObject.RNAIdentifier, CellEngineUseful::GetChainIdFromLetterForDNAorRNA('C'), -1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
        FillSquareParticle(StartXPosParam + 9, StartYPosParam + 19, StartZPosParam + 3, P12, 2, 2, 1);

        GenerateOneStrand(CellEngineConfigDataObject.RNAIdentifier, "UCGAGAA", StartXPosParam + 5, StartYPosParam + 5, StartZPosParam + 3, 2, 1, 2, 2, 0, 0);
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
        //REDRAW
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

void CellEngineVoxelSimulationSpace::GenerateOneStrand(const EntityIdInt EntityId, string_view Sequence, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt ParticleStepX, const UnsignedInt ParticleStepY, const UnsignedInt ParticleStepZ)
{
    try
    {
        vector<UniqueIdInt> MockVector;
        Particle* ParticlePrev1 = nullptr;

        for (UnsignedInt SequenceIndex = 0; SequenceIndex < Sequence.size(); SequenceIndex++)
            ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, EntityId, CellEngineUseful::GetChainIdFromLetterForDNAorRNA(Sequence[SequenceIndex]), 0, 1, StartPosX + SequenceIndex * ParticleStepX, StartPosY + SequenceIndex * ParticleStepY, StartPosZ + SequenceIndex * ParticleStepZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, false, MockVector, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
    }
    CATCH("generating one strand")
}

Particle* CellEngineVoxelSimulationSpace::GenerateNucleotideParticle(Particle* ParticlePrev, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeThread, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, bool AddToGenome, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam)
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
                    //REDRAW
                    GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
                    GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
                }

        if (AddToGenome == true)
            Genome.emplace_back(ParticleIndex);
    }
    CATCH("generating particle")

    return &GetParticleFromIndex(ParticleIndex);
}

tuple<Particle*, Particle*> CellEngineVoxelSimulationSpace::GenerateTwoPairedNucleotides(Particle* ParticlePrev1, Particle* ParticlePrev2, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam, const bool Linked)
{
    try
    {
        ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, CellEngineConfigDataObject.DNAIdentifier, ChainId, 0, Genomes[0].size(), StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
        ParticlePrev2 = GenerateNucleotideParticle(ParticlePrev2, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainIdForDNAorRNA(ChainId), 1, Genomes[1].size(), StartPosX + AddSizeX, StartPosY + AddSizeY, StartPosZ + AddSizeZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, true, Genomes[1], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

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
                ParticlePrev1 = GenerateNucleotideParticle(ParticlePrev1, stoi(Row[0]), stoi(Row[1]), 0, stoi(Row[2]), stoi(Row[3]), stoi(Row[4]), stoi(Row[5]), ParticleSize, ParticleSize, ParticleSize, true, Genomes[0], CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

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
            GetParticleFromIndex(Genomes[0][GenomeIndex + 1]).ChainId = CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]);
            GetParticleFromIndex(Genomes[1][GenomeIndex + 1]).ChainId = CellEngineUseful::GetPairedChainIdForDNAorRNA(CellEngineUseful::GetChainIdFromLetterForDNAorRNA(GenomesLines[0][GenomeIndex]));
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