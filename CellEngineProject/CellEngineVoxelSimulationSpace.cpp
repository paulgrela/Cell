
#include "FileUtils.h"

#include "CellEngineAtom.h"
#include "CellEngineColors.h"
#include "CellEngineUseful.h"
#include "CellEngineVoxelSimulationSpace.h"

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

inline SimulationSpaceVoxel GetZerSimulationSpaceVoxel()
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

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    GetSpaceVoxel(PosX, PosZ, PosY) = GetZerSimulationSpaceVoxel();
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

[[nodiscard]] std::stringstream CellEngineVoxelSimulationSpace::PrintSpaceMinMaxValues() const
{
    std::stringstream ss;
    ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << std::to_string(XMin) << " ][ Xmax = " << std::to_string(XMax) << " ][ Ymin = " << std::to_string(YMin) << " ][ Ymax = " << std::to_string(YMax) << " ][ Zmin = " << std::to_string(ZMin) << " ][ Zmax = " << std::to_string(XMax) << " ] " << std::endl;
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

        GetMinMaxOfCoordinates(PosX, PosY, PosZ);

        if (GetSpaceVoxel(PosX, PosY, PosZ) == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
            GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
        }
    }
    CATCH("setting atom in voxel simulation space")
}

void CellEngineVoxelSimulationSpace::AddParticleKind(const ParticleKind& ParticleParam)
{
    ParticlesKinds.emplace_back(ParticleParam);
}

UniqueIdInt CellEngineVoxelSimulationSpace::AddNewParticle(UniqueIdInt ParticleIndex, const Particle& ParticleParam)
{
    #ifdef PARTICLES_IN_VECTOR
    Particles.emplace_back(ParticleParam);
    return MaxParticleIndex = Particles.size() - 1;
    #else
    Particles[ParticleIndex] = ParticleParam;
    return MaxParticleIndex = ParticleIndex;
    #endif
}

void CellEngineVoxelSimulationSpace::AddReaction(const Reaction& ReactionParam)
{
    try
    {
        Reactions.emplace_back(ReactionParam);
        ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
    }
    CATCH("adding reaction")
}

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        vector<UnsignedInt> LocalNewParticlesIndexes;

        LoggersManagerObject.Log(STREAM("PARTICLES RANGE = " << MaxParticleIndex << " " << (MaxParticleIndex + NumberOfRandomParticles) << " " << NumberOfRandomParticles));

        UnsignedInt LocalMaxParticleIndex = MaxParticleIndex;
        for (UniqueIdInt ParticleIndex = LocalMaxParticleIndex + 1; ParticleIndex <= LocalMaxParticleIndex + NumberOfRandomParticles; ParticleIndex++)
            LocalNewParticlesIndexes.emplace_back(AddNewParticle(ParticleIndex, Particle(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, 1, CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor()))));

        vector<vector3_16> FilledVoxelsForRandomParticle;

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectObjectOfParticle_Uint64t(1, NumberOfRandomParticles);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    GetSpaceVoxel(PosX, PosY, PosZ) = GetZerSimulationSpaceVoxel();

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
                                    GetSpaceVoxel(VoxelForRandomParticle.X, VoxelForRandomParticle.Y, VoxelForRandomParticle.Z) = GetZerSimulationSpaceVoxel();

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
    }
    CATCH("generating random particles in selected space")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusion(UniqueIdInt StartParticleIndexParam, UniqueIdInt EndParticleIndexParam, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        std::uniform_int_distribution<SignedInt> UniformDistributionObjectSizeOfParticle_int64t(-1, 1);

        vector<vector3_16> NewVoxelsForParticle;

        if (EndParticleIndexParam == 0)
        {
            StartParticleIndexParam = MaxParticleIndex - StartParticleIndexParam;
            EndParticleIndexParam = MaxParticleIndex;
        }

        for (UniqueIdInt ParticleIndex = StartParticleIndexParam; ParticleIndex <= EndParticleIndexParam; ParticleIndex++)
        {
            auto& LocalParticleObject = Particles[ParticleIndex];

            for (auto& VoxelForParticle : LocalParticleObject.ListOfVoxels)
                GetSpaceVoxel(VoxelForParticle.X, VoxelForParticle.Y, VoxelForParticle.Z) = GetZerSimulationSpaceVoxel();

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
                        GetSpaceVoxel(NewVoxelForParticle.X, NewVoxelForParticle.Y, NewVoxelForParticle.Z) = GetZerSimulationSpaceVoxel();

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

std::mt19937_64 mt64RR;
vector<UniqueIdInt> Genome;

void CellEngineVoxelSimulationSpace::GenerateRandomDNAInWholeCell()
{
    try
    {
        for (auto& ParticleObjectLoop : Particles)
            if (ParticleObjectLoop.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                for (auto& VoxelCoordinates : ParticleObjectLoop.second.ListOfVoxels)
                    GetSpaceVoxel(ParticleObjectLoop.second.ListOfVoxels[0].X, ParticleObjectLoop.second.ListOfVoxels[0].Y, ParticleObjectLoop.second.ListOfVoxels[0].Z) = 0;

        const auto RemovedDNAParticlesCounter = std::erase_if(Particles, [](const pair<UniqueIdInt, Particle>& item) { auto const& [key, value] = item; return (value.EntityId == CellEngineConfigDataObject.DNAIdentifier); });
        LoggersManagerObject.Log(STREAM("RemovedDNAParticlesCounter = " << RemovedDNAParticlesCounter));

        UnsignedInt DNAParticleCounter = 0;
        for (auto& ParticleObjectLoop : Particles)
            if (ParticleObjectLoop.second.EntityId == CellEngineConfigDataObject.DNAIdentifier)
                DNAParticleCounter++;
        LoggersManagerObject.Log(STREAM("DNAParticleCounter = " << DNAParticleCounter));

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectChainOfParticle_Uint64t(1, 4);

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectMoveOfParticle_Uint64t(1, 6);
        UnsignedInt RandomPosX = 539, RandomPosY = 727, RandomPosZ = 587;
        UnsignedInt SizeX = 2, SizeY = 2, SizeZ = 2;

        UnsignedInt NumberOfGeneratedNucleotides = 0;
        vector<UnsignedInt> RandomMoveDirections = { 0, 0, 0, 0, 0, 0 };

        auto CheckIfAllRandomMovesDirectionsWereChecked = std::all_of(RandomMoveDirections.begin(), RandomMoveDirections.end(), [&] (const UnsignedInt &Element) { return Element == 1; });

        UnsignedInt NumberOfNucleotidesToBeGenerated = 10;

        while (NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
        {
            UnsignedInt RandomMoveDirection = 0;
            //losuj tak dlugo az trafisz w 0
            //osobno sprawdzic ta petle - niech wypisze Log bez pozostalych
            do
            {
                RandomMoveDirection = UniformDistributionObjectMoveOfParticle_Uint64t(mt64RR);
            }
            while (RandomMoveDirections[RandomMoveDirection] == 1 && CheckIfAllRandomMovesDirectionsWereChecked == false);

            if (RandomMoveDirections[RandomMoveDirection] == 0)
            {
                switch (RandomMoveDirection)
                {
                    case 1 : RandomPosX++; break;
                    case 2 : RandomPosX--; break;
                    case 3 : RandomPosY++; break;
                    case 4 : RandomPosY--; break;
                    case 5 : RandomPosZ++; break;
                    case 6 : RandomPosZ--; break;
                    default : break;
                }
                RandomMoveDirections[RandomMoveDirection] = 1;
            }

            bool EmptyVoxelSpaceForNewNucleotideBool = true;
            while (EmptyVoxelSpaceForNewNucleotideBool == true && NumberOfGeneratedNucleotides < NumberOfNucleotidesToBeGenerated)
            {
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + SizeX; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + SizeY; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + SizeZ; PosZ++)
                            if (GetSpaceVoxel(PosX, PosY, PosZ) != 0)
                            {
                                EmptyVoxelSpaceForNewNucleotideBool = false;
                                break;
                            }

                if (EmptyVoxelSpaceForNewNucleotideBool == true)
                {
                    UnsignedInt ParticleIndex = Particles.size();

                    NumberOfGeneratedNucleotides++;

                    ParticleIndex++;
                    AddNewParticle(ParticleIndex, Particle(ParticleIndex, CellEngineConfigDataObject.DNAIdentifier, UniformDistributionObjectChainOfParticle_Uint64t(mt64RR), CellEngineUseful::GetVector3FormVMathVec3ForColor(CellEngineColorsObject.GetRandomColor())));
                    for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + SizeX; PosX++)
                        for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + SizeY; PosY++)
                            for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + SizeZ; PosZ++)
                            {
                                GetSpaceVoxel(PosX, PosY, PosZ) = ParticleIndex;
                                GetParticleFromIndex(ParticleIndex).ListOfVoxels.emplace_back(PosX, PosY, PosZ);
                            }
                    Genome.emplace_back(ParticleIndex);
                }
            }

            if (EmptyVoxelSpaceForNewNucleotideBool == false && CheckIfAllRandomMovesDirectionsWereChecked == true)
            {
                UnsignedInt PreviousParticleIndex = Particles.size() - 1;
                RandomPosX = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].X;
                RandomPosY = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Y;
                RandomPosZ = GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels[0].Z;
                for (auto& VoxelCoordinates : GetParticleFromIndex(PreviousParticleIndex).ListOfVoxels)
                    GetSpaceVoxel(VoxelCoordinates.X, VoxelCoordinates.Y, VoxelCoordinates.Z) = GetZerSimulationSpaceVoxel();
                Particles.erase(Particles.size() - 1);

                fill(RandomMoveDirections.begin(), RandomMoveDirections.end(), 0);
            }
        }
    }
    CATCH("generating random dna in whole cell")
}