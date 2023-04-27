
#include "CellEngineAtom.h"
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

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelsInVoxelSimulationSpaceInEachDimension; PosZ++)
                    GetSpaceVoxel(PosX, PosZ, PosY) = {0, 0 };
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
                    if (GetSpaceVoxel(PosX, PosY, PosZ).EntityId != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
    CATCH("counting statistics of voxel simulation space")
}

void CellEngineVoxelSimulationSpace::SetAtomInVoxelSimulationSpace(const CellEngineAtom& AppliedAtom)
{
    try
    {
        UnsignedInt PosX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt PosY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt PosZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(PosX, PosY, PosZ);

        if (GetSpaceVoxel(PosX, PosY, PosZ).EntityId == 0)
        {
            GetSpaceVoxel(PosX, PosY, PosZ).EntityId = AppliedAtom.EntityId;
            if (CellEngineConfigDataObject.IsDNAorRNA(AppliedAtom.EntityId) == true)
                GetSpaceVoxel(PosX, PosY, PosZ).ChainId = stoi(std::string(AppliedAtom.Chain).substr(2, 2));
        }
    }
    CATCH("setting atom in voxel simulation space")
}

SimulationSpaceVoxel CellEngineVoxelSimulationSpace::GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z)
{
    return GetSpaceVoxel(X, Y, Z);
}

void CellEngineVoxelSimulationSpace::AddParticleKind(const ParticleKind& ParticleParam)
{
    Particles.emplace_back(ParticleParam);
}

void CellEngineVoxelSimulationSpace::AddReaction(const Reaction& ReactionParam)
{
    Reactions.emplace_back(ReactionParam);
    ReactionsIdByString.insert(std::make_pair(ReactionParam.ReactantsStr, Reactions.size() - 1));
}

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt StepXParam, const UnsignedInt StepYParam, const UnsignedInt StepZParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        Particles.clear();

        AddParticleKind({ 1, "A", "A", 0, { { "A + B + ", 0 }, { "D + A + ", 3 } } });
        AddParticleKind({ 2, "B", "B", 0, { { "A + B + ", 0 }, { "B + C + ", 1 } } });
        AddParticleKind({ 3, "C", "C", 0, { { "B + C + ", 1 }, { "C + D + ", 2 } } });
        AddParticleKind({ 4, "D", "D", 0, { { "C + D + ", 2 }, { "D + A + ", 3 } } });

        AddReaction(Reaction("A + B + ", { {1, 1}, {2, 1} }, { {3, 1}, {4, 1} }));
        AddReaction(Reaction("B + C + ", { {2, 1}, {3, 1} }, { {4, 1} }));
        AddReaction(Reaction("C + D + ", { {3, 1}, {4, 1} }, { {1, 1} }));
        AddReaction(Reaction("D + A + ", { {4, 1}, {1, 1} }, { {2, 1} }));

        vector<vector3_64> FilledVoxelsForRandomParticle;

        uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectTypeOfParticle_Uint64t(1, 4);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(StartXPosParam, StartXPosParam + SizeXParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(StartYPosParam, StartYPosParam + SizeYParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(StartZPosParam, StartZPosParam + SizeZParam);

        for (UnsignedInt PosX = StartXPosParam; PosX < StartXPosParam + SizeXParam; PosX += StepXParam)
            for (UnsignedInt PosY = StartYPosParam; PosY < StartYPosParam + SizeYParam; PosY += StepYParam)
                for (UnsignedInt PosZ = StartZPosParam; PosZ < StartZPosParam + SizeZParam; PosZ += StepZParam)
                    GetSpaceVoxel(PosX, PosY, PosZ) = {0, 0 };

        for (UnsignedInt ParticleNum = 1;  ParticleNum <= NumberOfRandomParticles; ParticleNum++)
        {
            UnsignedInt RandomPosX = UniformDistributionObjectX_Uint64t(mt64R);
            UnsignedInt RandomPosY = UniformDistributionObjectY_Uint64t(mt64R);
            UnsignedInt RandomPosZ = UniformDistributionObjectZ_Uint64t(mt64R);

            UnsignedInt RandomSizeOfParticle = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

            ChainIdInt RandomChainId = UniformDistributionObjectTypeOfParticle_Uint64t(mt64R);

            FilledVoxelsForRandomParticle.clear();

            if (RandomPosX + RandomSizeOfParticle < StartXPosParam + SizeXParam && RandomPosY + RandomSizeOfParticle < StartYPosParam + SizeYParam && RandomPosZ + RandomSizeOfParticle < StartZPosParam + SizeZParam)
                for (UnsignedInt PosX = RandomPosX; PosX < RandomPosX + RandomSizeOfParticle; PosX++)
                    for (UnsignedInt PosY = RandomPosY; PosY < RandomPosY + RandomSizeOfParticle; PosY++)
                        for (UnsignedInt PosZ = RandomPosZ; PosZ < RandomPosZ + RandomSizeOfParticle; PosZ++)
                        {
                            if (GetSpaceVoxel(PosX, PosY, PosZ).EntityId == 0)
                            {
                                FilledVoxelsForRandomParticle.emplace_back(PosX, PosY, PosZ);
                                GetSpaceVoxel(PosX, PosY, PosZ) = {static_cast<EntityIdInt>(CellEngineConfigDataObject.DNAIdentifier), RandomChainId };
                            }
                            else
                            {
                                for (auto& VoxelForRandomParticle : FilledVoxelsForRandomParticle)
                                    GetSpaceVoxel(VoxelForRandomParticle.X, VoxelForRandomParticle.Y, VoxelForRandomParticle.Z) = {0, 0 };


                                goto NextRandomParticleOutsideLoopLabel;
                            }
                        }
            NextRandomParticleOutsideLoopLabel:;

            if (FilledVoxelsForRandomParticle.empty() == false)
                Particles[RandomChainId - 1].ParticlesObjects.emplace_back(FilledVoxelsForRandomParticle);
        }
    }
    CATCH("generating random particles in selected space")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusion(const UnsignedInt StartXPosParam, const UnsignedInt StartYPosParam, const UnsignedInt StartZPosParam, const UnsignedInt SizeXParam, UnsignedInt SizeYParam, const UnsignedInt SizeZParam)
{
    try
    {
        std::uniform_int_distribution<SignedInt> UniformDistributionObjectSizeOfParticle_int64t(-1, 1);

        vector<vector3_64> NewVoxelsForParticle;

        for (auto& ParticleKindObject : Particles)
            for (auto& ParticleObject : ParticleKindObject.ParticlesObjects)
            {
                ChainIdInt RememberChainId = GetSpaceVoxel(ParticleObject.ListOfVoxels[0].X, ParticleObject.ListOfVoxels[0].Y, ParticleObject.ListOfVoxels[0].Z).ChainId;
                for (auto& VoxelForParticle : ParticleObject.ListOfVoxels)
                    GetSpaceVoxel(VoxelForParticle.X, VoxelForParticle.Y, VoxelForParticle.Z) = {0, 0 };

                SignedInt ShiftX = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
                SignedInt ShiftY = UniformDistributionObjectSizeOfParticle_int64t(mt64R);
                SignedInt ShiftZ = UniformDistributionObjectSizeOfParticle_int64t(mt64R);

                NewVoxelsForParticle.clear();
                bool Collision = false;

                for (auto& VoxelForParticle : ParticleObject.ListOfVoxels)
                    if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ).EntityId == 0 && VoxelForParticle.X + ShiftX >= StartXPosParam && VoxelForParticle.X + ShiftX < StartXPosParam + SizeXParam && VoxelForParticle.Y + ShiftY >= StartYPosParam && VoxelForParticle.Y + ShiftY < StartYPosParam + SizeYParam && VoxelForParticle.Z + ShiftZ >= StartZPosParam && VoxelForParticle.Z + ShiftZ < StartZPosParam + SizeZParam)
                    {
                        GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ) = {static_cast<EntityIdInt>(CellEngineConfigDataObject.DNAIdentifier), RememberChainId };
                        NewVoxelsForParticle.emplace_back(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ);
                    }
                    else
                    {
                        for (auto& NewVoxelForParticle : NewVoxelsForParticle)
                            GetSpaceVoxel(NewVoxelForParticle.X, NewVoxelForParticle.Y, NewVoxelForParticle.Z) = {0, 0 };

                        for (auto& OldVoxelForParticle : ParticleObject.ListOfVoxels)
                            GetSpaceVoxel(OldVoxelForParticle.X, OldVoxelForParticle.Y, OldVoxelForParticle.Z) = {static_cast<EntityIdInt>(CellEngineConfigDataObject.DNAIdentifier), RememberChainId };

                        Collision = true;
                        break;
                    }

                if (Collision == false)
                    ParticleObject.ListOfVoxels = NewVoxelsForParticle;
            }
    }
    CATCH("generating one step of diffusion")
}
