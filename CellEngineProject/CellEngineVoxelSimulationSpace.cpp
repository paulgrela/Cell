
#include "FileUtils.h"
#include "DestinationPlatform.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
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

//inline Particle& CellEngineVoxelSimulationSpace::GetParticleFromIndex(const UniqueIdInt ParticleIndex)
//{
//    return Particles[ParticleIndex];
//}

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
    return CellEngineSimulationManagerObject.GetParticleFromIndex(ParticleIndex);
}

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    GetSpaceVoxel(PosX, PosZ, PosY) = GetZeroSimulationSpaceVoxel();
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

void CellEngineVoxelSimulationSpace::SetStartValuesForSpaceMinMax()
{
    XMin = YMin = ZMin = 10000;
    XMax = YMax = ZMax = 0;
}

void CellEngineVoxelSimulationSpace::GetMinMaxOfCoordinates(const UnsignedInt PosX, const UnsignedInt PosY, const UnsignedInt PosZ)
{
    try
    {
        XMin = min(PosX, XMin);
        XMax = max(PosX, XMax);
        YMin = min(PosY, YMin);
        YMax = max(PosY, YMax);
        ZMin = min(PosZ, ZMin);
        ZMax = max(PosZ, ZMax);
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
                    if (CellEngineSimulationManagerObject.GetParticleFromIndex(GetSpaceVoxel(PosX, PosY, PosZ)).EntityId != 0)
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

        GetMinMaxOfCoordinates(PosX, PosY, PosZ);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
            CellEngineSimulationManagerObject.GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
        }
    }
    CATCH("setting atom in voxel simulation space")
}

//void CellEngineVoxelSimulationSpace::AddParticleKind(const ParticleKind& ParticleParam)
//{
//    ParticlesKinds.emplace_back(ParticleParam);
//}
//
//UniqueIdInt CellEngineVoxelSimulationSpace::AddNewParticle(UniqueIdInt ParticleIndex, const Particle& ParticleParam)
//{
//    #ifdef PARTICLES_IN_VECTOR
//    Particles.emplace_back(ParticleParam);
//    return MaxParticleIndex = Particles.size() - 1;
//    #else
//    Particles[ParticleIndex] = ParticleParam;
//    return MaxParticleIndex = ParticleIndex;
//    #endif
//}

//void CellEngineVoxelSimulationSpace::AddReaction(const Reaction& ReactionParam)
//{
//    try
//    {
//        Reactions.emplace_back(ReactionParam);
//        ReactionsIdByString.insert(make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
//    }
//    CATCH("adding reaction")
//}

void CellEngineVoxelSimulationSpace::AddBasicParticlesKindsAndReactions()
{
    try
    {
        CellEngineSimulationManagerObject.AddParticleKind({0, "Water", "H2O", 100 });
        CellEngineSimulationManagerObject.AddParticleKind({1, "Glucose", "C6H12O6", 50 });
        CellEngineSimulationManagerObject.AddParticleKind({2, "Oxygen6", "06", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({3, "Carbon dioxide", "CO2", 5 });
        CellEngineSimulationManagerObject.AddParticleKind({4, "Eten", "CH2CH2", 15 });
        CellEngineSimulationManagerObject.AddParticleKind({5, "Ethanol", "CH3CH2(OH)", 25 });
        CellEngineSimulationManagerObject.AddParticleKind({6, "Propen", "CH3CHCH2", 5 });
        CellEngineSimulationManagerObject.AddParticleKind({7, "HX", "HX", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({8, "2Halogenopropan", "CH3CHXCH3", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({9, "Eten", "CH2CH2", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({10, "Ethylene", "CH2CH2O", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({11, "Oxygen", "0", 10 });
        CellEngineSimulationManagerObject.AddParticleKind({12, "DNA", "CGATATTAAATAGGGCCT", 10 });

        CellEngineSimulationManagerObject.AddReaction(Reaction("C6H12O6 + O6 + ", { { 1, 1 }, { 2, 2 } }, { { 0, 6 }, { 0, 6 } }));
        CellEngineSimulationManagerObject.AddReaction(Reaction("CH2CH2 + H2O + ", { { 4, 1 }, { 0, 1 } }, { { 5, 1 } }));
        CellEngineSimulationManagerObject.AddReaction(Reaction("CH3CHCH2 + HX + ", { { 6, 1 }, { 7, 1 } }, { { 8, 1 } }));
        CellEngineSimulationManagerObject.AddReaction(Reaction("CH2CH2 + O + ", { { 9, 1 }, { 11, 1 } }, { { 10, 1 } }));
    }
    CATCH("adding particles kinds and reactions")
};

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << CellEngineSimulationManagerObject.MaxParticleIndex << " " << (CellEngineSimulationManagerObject.MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(1, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        AddBasicParticlesKindsAndReactions();

        vector<UnsignedInt> LocalNewParticlesIndexes;

        UnsignedInt LocalMaxParticleIndex = CellEngineSimulationManagerObject.MaxParticleIndex;
        for (UniqueIdInt ParticleIndex = LocalMaxParticleIndex + 1; ParticleIndex <= LocalMaxParticleIndex + NumberOfRandomParticles; ParticleIndex++)
            //LocalNewParticlesIndexes.emplace_back(CellEngineSimulationManagerObject.AddNewParticle(ParticleIndex, Particle(ParticleIndex, UniformDistributionObjectObjectOfParticle_Uint64t(mt64R), 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));
            LocalNewParticlesIndexes.emplace_back(CellEngineSimulationManagerObject.AddNewParticle(ParticleIndex, Particle(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, 1, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));
            //blad segmentation fault bo nie ma takiego ParticleKind podczas rysowania

        vector<vector3_16> FilledVoxelsForRandomParticle;

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    GetSpaceVoxel(PosX, PosY, PosZ) = GetZeroSimulationSpaceVoxel();

        for (auto& LocalNewParticleIndex : LocalNewParticlesIndexes)
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
                                for (auto& VoxelForRandomParticle : FilledVoxelsForRandomParticle)
                                    GetSpaceVoxel(VoxelForRandomParticle.X, VoxelForRandomParticle.Y, VoxelForRandomParticle.Z) = GetZeroSimulationSpaceVoxel();

                                FilledVoxelsForRandomParticle.clear();

                                goto NextRandomParticleOutsideLoopLabel;
                            }
                        }
            NextRandomParticleOutsideLoopLabel:;

            if (FilledVoxelsForRandomParticle.empty() == false)
            {
                LoggersManagerObject.Log(STREAM("LocalNewParticleIndex = " << LocalNewParticleIndex));
                CellEngineSimulationManagerObject.Particles[LocalNewParticleIndex].ListOfVoxels = FilledVoxelsForRandomParticle;
            }
        }
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

        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, CellEngineSimulationManagerObject.MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto& LocalParticleObject = CellEngineSimulationManagerObject.Particles[ParticleIndex];

            for (auto& VoxelForParticle : LocalParticleObject.ListOfVoxels)
                GetSpaceVoxel(VoxelForParticle.X, VoxelForParticle.Y, VoxelForParticle.Z) = GetZeroSimulationSpaceVoxel();

            SignedInt ShiftX = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftY = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
            SignedInt ShiftZ = UniformDistributionObjectSizeOfParticle_int64t(mt64R);

            NewVoxelsForParticle.clear();
            bool Collision = false;

            for (auto& VoxelForParticle : LocalParticleObject.ListOfVoxels)
                if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                {
                    GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) = ParticleIndex;
                    NewVoxelsForParticle.emplace_back(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ);
                }
                else
                {
                    for (auto& NewVoxelForParticle : NewVoxelsForParticle)
                        GetSpaceVoxel(NewVoxelForParticle.X, NewVoxelForParticle.Y, NewVoxelForParticle.Z) = GetZeroSimulationSpaceVoxel();

                    for (auto& OldVoxelForParticle : LocalParticleObject.ListOfVoxels)
                        GetSpaceVoxel(OldVoxelForParticle.X, OldVoxelForParticle.Y, OldVoxelForParticle.Z) = ParticleIndex;

                    Collision = true;
                    break;
                }

            if (Collision == false)
                LocalParticleObject.ListOfVoxels = NewVoxelsForParticle;
        }
    }
    CATCH("generating one step of diffusion")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionForParticle(Particle& ParticleObject)
{
    try
    {

    }
    CATCH("generating random reaction for particle")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfRandomReactionsForSelectedRangeOfParticles(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, UnsignedInt StartXPosParam, UnsignedInt StartYPosParam, UnsignedInt StartZPosParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam)
{
    try
    {
        GetRangeOfParticlesForRandomParticles(StartParticleIndexParam, EndParticleIndexParam, CellEngineSimulationManagerObject.MaxParticleIndex);

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto ParticlesIterator = CellEngineSimulationManagerObject.Particles.find(ParticleIndex);
            if (ParticlesIterator != CellEngineSimulationManagerObject.Particles.end())
                GenerateRandomReactionForParticle(ParticlesIterator->second);
        }
    }
    CATCH("generating one step of random reactions for selected range of particles")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionsForAllParticles()
{
    try
    {
        for (auto& ParticleObject : CellEngineSimulationManagerObject.Particles)
            if (ParticleObject.second.SelectedForReaction == false)
                GenerateRandomReactionForParticle(ParticleObject.second);
    }
    CATCH("generating random reactions for all particles")
}

void CellEngineVoxelSimulationSpace::GenerateRandomReactionsInWholeVoxelSimulationSpace()
{
    try
    {
        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    if (GetSpaceVoxel(PosX, PosZ, PosY) != 0)
                    {
                        auto& ParticleObject = CellEngineSimulationManagerObject.GetParticleFromIndex(GetSpaceVoxel(PosX, PosZ, PosY));
                        if (ParticleObject.SelectedForReaction == false)
                            GenerateRandomReactionForParticle(ParticleObject);
                    }
    }
    CATCH("generation random reactions in whole voxel simulation space")
}

void CellEngineVoxelSimulationSpace::EraseAllDNAParticles()
{
    try
    {
        for (auto& ParticleObject : CellEngineSimulationManagerObject.Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                for (auto& VoxelCoordinates : ParticleObject.second.ListOfVoxels)
                    GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZeroSimulationSpaceVoxel();

        const auto RemovedDNAParticlesCounter = erase_if(CellEngineSimulationManagerObject.Particles, [](const pair<UniqueIdInt, Particle>& item) { auto const& [key, value] = item; return (value.EntityId == CellEngineConfigDataObject.DNAIdentifier); });
        LoggersManagerObject.Log(STREAM("RemovedDNAParticlesCounter = " << RemovedDNAParticlesCounter));

        UnsignedInt DNAParticleCounter = 0;
        for (auto& ParticleObject : CellEngineSimulationManagerObject.Particles)
            if (ParticleObject.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                DNAParticleCounter++;
        LoggersManagerObject.Log(STREAM("DNAParticleCounter = " << DNAParticleCounter));

        CellEngineSimulationManagerObject.Genome1.clear();
        CellEngineSimulationManagerObject.Genome2.clear();
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

void CellEngineVoxelSimulationSpace::GenerateParticle(UnsignedInt& ParticleIndex, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, vector<UniqueIdInt>& Genome, const vector3_16 UniqueColorParam)
{
    try
    {
        ParticleIndex++;
        CellEngineSimulationManagerObject.AddNewParticle(ParticleIndex, Particle(ParticleIndex, EntityId, ChainId, GenomeIndex, UniqueColorParam));

        for (UnsignedInt PosX = StartPosX; PosX < StartPosX + ParticleSizeX; PosX++)
            for (UnsignedInt PosY = StartPosY; PosY < StartPosY + ParticleSizeY; PosY++)
                for (UnsignedInt PosZ = StartPosZ; PosZ < StartPosZ + ParticleSizeZ; PosZ++)
                {
                    GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
                    CellEngineSimulationManagerObject.GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
                }

        Genome.emplace_back(ParticleIndex);
    }
    CATCH("generating particle")
}

void CellEngineVoxelSimulationSpace::GenerateTwoPairedNucleotides(UnsignedInt& ParticleIndex, const EntityIdInt EntityId, const ChainIdInt ChainId, const UnsignedInt GenomeIndex, const UnsignedInt StartPosX, const UnsignedInt StartPosY, const UnsignedInt StartPosZ, const UnsignedInt ParticleSizeX, const UnsignedInt ParticleSizeY, const UnsignedInt ParticleSizeZ, const UnsignedInt AddSizeX, const UnsignedInt AddSizeY, const UnsignedInt AddSizeZ, const vector3_16 UniqueColorParam)
{
    try
    {
        GenerateParticle(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, ChainId, CellEngineSimulationManagerObject.Genome1.size(), StartPosX, StartPosY, StartPosZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, CellEngineSimulationManagerObject.Genome1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
        GenerateParticle(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, CellEngineUseful::GetPairedChainId(ChainId), CellEngineSimulationManagerObject.Genome2.size(), StartPosX + AddSizeX, StartPosY + AddSizeY, StartPosZ + AddSizeZ, ParticleSizeX, ParticleSizeY, ParticleSizeZ, CellEngineSimulationManagerObject.Genome2, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
    }
    CATCH("generating two paired nucleotides")
}

bool CellEngineVoxelSimulationSpace::TestFormerForbiddenPositions(unordered_set<string>& TestedFormerForbiddenPositions, UnsignedInt RandomMoveDirection, UnsignedInt RandomPosX, UnsignedInt RandomPosY, UnsignedInt RandomPosZ, const UnsignedInt Size)
{
    UnsignedInt PosX = RandomPosX;
    UnsignedInt PosY = RandomPosY;
    UnsignedInt PosZ = RandomPosZ;

    UpdateRandomPositions(RandomMoveDirection, PosX, PosY, PosZ, Size);

    return (TestedFormerForbiddenPositions.find(to_string(PosX) + "|" + to_string(PosY) + "|" + to_string(PosZ)) != TestedFormerForbiddenPositions.end());
}

tuple<UnsignedInt, UnsignedInt, UnsignedInt> CellEngineVoxelSimulationSpace::EraseLastRandomParticle(vector<UniqueIdInt>& Genome)
{
    UnsignedInt LocalRandomPosX, LocalRandomPosY, LocalRandomPosZ;

    try
    {
        UnsignedInt PreviousParticleIndex = Genome.back();
        Genome.pop_back();

        LocalRandomPosX = CellEngineSimulationManagerObject.GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].X;
        LocalRandomPosY = CellEngineSimulationManagerObject.GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Y;
        LocalRandomPosZ = CellEngineSimulationManagerObject.GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Z;

        for (auto& VoxelCoordinates : CellEngineSimulationManagerObject.GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels)
            GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZeroSimulationSpaceVoxel();

        CellEngineSimulationManagerObject.Particles.erase(PreviousParticleIndex);

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

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = CellEngineSimulationManagerObject.Particles.size();

        LoggersManagerObject.InitializePrintingParameters(false, false, false, false, false, false, false, false, false, false, false, false, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);

        UnsignedInt NumberOfGeneratedNucleotides = 0;

        unordered_set<string> TestedFormerForbiddenPositions;

        vector<UnsignedInt> RandomMovesDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = [](vector<UnsignedInt>& RandomMovesDirections) { return all_of(RandomMovesDirections.cbegin(), RandomMovesDirections.cend(), [] (const UnsignedInt Element) { return Element == 1; }); };

        auto timeStart = clock();

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

                    UnsignedInt ParticleIndex = CellEngineSimulationManagerObject.MaxParticleIndex;

                    if (RandomMoveDirection == 1 || RandomMoveDirection == 2)
                        GenerateTwoPairedNucleotides(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, ChainId, CellEngineSimulationManagerObject.Genome1.size(), RandomPosX, RandomPosY, RandomPosZ, ParticleSizeX, 1, ParticleSizeZ, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                    else
                    if (RandomMoveDirection == 3 || RandomMoveDirection == 4 || RandomMoveDirection == 5 || RandomMoveDirection == 6)
                        GenerateTwoPairedNucleotides(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, ChainId, CellEngineSimulationManagerObject.Genome1.size(), RandomPosX, RandomPosY, RandomPosZ, 1, ParticleSizeY, ParticleSizeZ, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

                    fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);

                    UpdateRandomPositions(RandomMoveDirection, RandomPosX, RandomPosY, RandomPosZ, ParticleSize5);

                    LoggersManagerObject.Log(STREAM("ADDED PARTICLE NumberOfGeneratedNucleotides = " << NumberOfGeneratedNucleotides << " ParticleIndex = " << ParticleIndex << " RX = " << RandomPosX << " RY = " << RandomPosY << " RZ = " << RandomPosZ));
                }
            }

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked(RandomMovesDirections) == true)
            {
                NumberOfGeneratedNucleotides--;

                tie(RandomPosX, RandomPosY, RandomPosZ) = EraseLastRandomParticle(CellEngineSimulationManagerObject.Genome1);

                TestedFormerForbiddenPositions.insert(to_string(RandomPosX) + "|" + to_string(RandomPosY) + "|" + to_string(RandomPosZ));

                EraseLastRandomParticle(CellEngineSimulationManagerObject.Genome2);

                fill(RandomMovesDirections.begin(), RandomMovesDirections.end(), 0);
            }

            LoggersManagerObject.Log(STREAM("END OF GOING IN ONE DIRECTION"));
        }

        LoggersManagerObject.InitializePrintingParameters(CellEngineConfigDataObject.PrintLogToConsole, CellEngineConfigDataObject.PrintLogToFiles, CellEngineConfigDataObject.PrintLogLineNumberToConsole, CellEngineConfigDataObject.PrintLogDateTimeToConsole, CellEngineConfigDataObject.PrintLogProcessIdToConsole, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToConsole, CellEngineConfigDataObject.PrintLogThreadIdToConsole, CellEngineConfigDataObject.PrintLogLineNumberToFile, CellEngineConfigDataObject.PrintLogDateTimeToFile, CellEngineConfigDataObject.PrintLogProcessIdToFile, CellEngineConfigDataObject.PrintLogProcessPriorityLevelToFile, CellEngineConfigDataObject.PrintLogThreadIdToFile, CellEngineConfigDataObject.MaximalNumberOfLinesInOneFile);

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << CellEngineSimulationManagerObject.Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("generating random dna in whole cell")
}

void CellEngineVoxelSimulationSpace::SaveGenomeDataToFile(UnsignedInt ParticleSize)
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        ofstream FileToWriteGenome;
        FileToWriteGenome.open(FileName, ios_base::out | ios_base::trunc);
        FileToWriteGenome << to_string(ParticleSize) << endl;
        FileToWriteGenome << to_string(CellEngineSimulationManagerObject.Genome1.size()) << endl;
        for (const auto& Nucleotide : CellEngineSimulationManagerObject.Genome1)
        {
            EntityIdInt EntityId = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).EntityId;
            ChainIdInt ChainId = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).ChainId;
            UniqueIdInt GenomeIndex = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).GenomeIndex;
            UnsignedInt PosX = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).ListOfVoxels[0].X;
            UnsignedInt PosY = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Y;
            UnsignedInt PosZ = CellEngineSimulationManagerObject.GetParticleFromIndex(Nucleotide).ListOfVoxels[0].Z;
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
        EraseAllDNAParticles();

        UnsignedInt ParticlesSizeBeforeAddingRandomDNA = CellEngineSimulationManagerObject.Particles.size();

        string Line, Word;
        vector<string> Row;

        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_POSITIONS.DAT");
        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome, Line);
        UnsignedInt ParticleSize = stoi(Line);
        getline(FileToReadGenome, Line);
        UnsignedInt GenomeSize = stoi(Line);

        UnsignedInt ParticleIndex = CellEngineSimulationManagerObject.MaxParticleIndex;
        UnsignedInt PrevStartPosX = 0;
        UnsignedInt GenomeIndex = 0;

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
                    GenerateTwoPairedNucleotides(ParticleIndex, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, ParticleSize, 1, ParticleSize, 0, 1, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
                else
                    GenerateTwoPairedNucleotides(ParticleIndex, stoi(Row[0]), stoi(Row[1]), GenomeIndex, StartPosX, StartPosY, StartPosZ, 1, ParticleSize, ParticleSize, 1, 0, 0, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));
            }
            else
                GenerateParticle(ParticleIndex, stoi(Row[0]), stoi(Row[1]), stoi(Row[2]), stoi(Row[3]), stoi(Row[4]), stoi(Row[5]), ParticleSize, ParticleSize, ParticleSize, CellEngineSimulationManagerObject.Genome1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()));

            PrevStartPosX = StartPosX;

            GenomeIndex++;
        }

        FileToReadGenome.close();

        LoggersManagerObject.Log(STREAM("NUMBER OF ADDED PARTICLES = " << CellEngineSimulationManagerObject.Particles.size() - ParticlesSizeBeforeAddingRandomDNA));
    }
    CATCH("reading genome data from file")
}

void CellEngineVoxelSimulationSpace::ReadGenomeSequenceFromFile()
{
    try
    {
        string FileName = string(".") + OS_DIR_SEP + string("data") + OS_DIR_SEP + string("genome") + OS_DIR_SEP + string("GENOME_SEQUENCE.DAT");

        fstream FileToReadGenome(FileName, ios::in);

        getline(FileToReadGenome,  CellEngineSimulationManagerObject.GenomeLine);

        for (UnsignedInt GenomeIndex = 0; GenomeIndex <  CellEngineSimulationManagerObject.Genome1.size(); GenomeIndex++)
        {
            CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ChainId = CellEngineUseful::GetChainIdFromLetter(CellEngineSimulationManagerObject.GenomeLine[GenomeIndex]);
            CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome2[GenomeIndex]).ChainId = CellEngineUseful::GetPairedChainId(CellEngineUseful::GetChainIdFromLetter(CellEngineSimulationManagerObject.GenomeLine[GenomeIndex]));
        }
    }
    CATCH("reading real genome data from file")
}

void CellEngineVoxelSimulationSpace::TestGeneratedGenomeCorrectness(const UnsignedInt ParticleSize)
{
    try
    {
        bool FoundBreak = false;

        for (UnsignedInt GenomeIndex = 0; GenomeIndex < CellEngineSimulationManagerObject.Genome1.size(); GenomeIndex++)
            if (GenomeIndex > 3)
            {
                if ((abs(CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].X - CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - 1]).ListOfVoxels[0].X) != ParticleSize) &&
                    (abs(CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].Y - CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - 1]).ListOfVoxels[0].Y) != ParticleSize) &&
                    (abs(CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].Z - CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - 1]).ListOfVoxels[0].Z) != ParticleSize))
                {
                    LoggersManagerObject.Log(STREAM("GenomeIndex = " << GenomeIndex << " ChainId = " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ChainId << " Letter = " << CellEngineUseful::GetLetterForDNAChainId(CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ChainId)));
                    LoggersManagerObject.Log(STREAM("DIFF X = " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].X << " " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].X));
                    LoggersManagerObject.Log(STREAM("DIFF Y = " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].Y << " " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].Y));
                    LoggersManagerObject.Log(STREAM("DIFF Z = " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex]).ListOfVoxels[0].Z << " " << CellEngineSimulationManagerObject.GetParticleFromIndex(CellEngineSimulationManagerObject.Genome1[GenomeIndex - ParticleSize]).ListOfVoxels[0].Z));
                    FoundBreak = true;
                    break;
                }
            }

        if (FoundBreak == false)
            LoggersManagerObject.Log(STREAM("Genome is continuous and correctly generated - OK!"));
    }
    CATCH("testing generated genome correctness")
}