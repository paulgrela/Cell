
#include "CellEngineAtom.h"
#include "CellEngineVoxelSimulationSpace.h"

using namespace std;

[[nodiscard]] float CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<SignedInt>(CoordinateParam) - (static_cast<SignedInt>(CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2))) * CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace;
};
[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam) / CellEngineConfigDataObject.DivisionFactorForVoxelSimulationSpace) + (CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2);
};

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SpacePointer = (CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension == 2048 ? malloc(sizeof(Space_2048_2048_2048)) : malloc(sizeof(Space_1024_1024_1024)));

        SetStartValuesForSpaceMinMax();

        for (UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosX++)
            for (UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosY++)
                for (UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosZ++)
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

void CellEngineVoxelSimulationSpace::GetMinMaxOfCoordinates(const UnsignedInt SpaceX, const UnsignedInt SpaceY, const UnsignedInt SpaceZ)
{
    try
    {
        XMin = min(SpaceX, XMin);
        XMax = max(SpaceX, XMax);
        YMin = min(SpaceY, YMin);
        YMax = max(SpaceY, YMax);
        ZMin = min(SpaceZ, ZMin);
        ZMax = max(SpaceZ, ZMax);
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
        for(UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosX++)
            for(UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosY++)
                for(UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosZ++)
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
        UnsignedInt SpaceX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt SpaceY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt SpaceZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        GetMinMaxOfCoordinates(SpaceX, SpaceY, SpaceZ);

        if (GetSpaceVoxel(SpaceX, SpaceY, SpaceZ).EntityId == 0)
        {
            GetSpaceVoxel(SpaceX, SpaceY, SpaceZ).EntityId = AppliedAtom.EntityId;
            if (CellEngineConfigDataObject.IsDNAorRNA(AppliedAtom.EntityId) == true)
                GetSpaceVoxel(SpaceX, SpaceY, SpaceZ).ChainId = stoi(std::string(AppliedAtom.Chain).substr(2, 2));
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

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
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
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(XStartParam, XStartParam + XSizeParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(YStartParam, YStartParam + YSizeParam);
        uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(ZStartParam, ZStartParam + ZSizeParam);

        for (UnsignedInt SpaceXP = XStartParam; SpaceXP < XStartParam + XSizeParam; SpaceXP += XStepParam)
            for (UnsignedInt SpaceYP = YStartParam; SpaceYP < YStartParam + YSizeParam; SpaceYP += YStepParam)
                for (UnsignedInt SpaceZP = ZStartParam; SpaceZP < ZStartParam + ZSizeParam; SpaceZP += ZStepParam)
                    GetSpaceVoxel(SpaceXP, SpaceYP, SpaceZP) = { 0, 0 };

        for (UnsignedInt ParticleNum = 1;  ParticleNum <= NumberOfRandomParticles; ParticleNum++)
        {
            UnsignedInt RandomPosX = UniformDistributionObjectX_Uint64t(mt64R);
            UnsignedInt RandomPosY = UniformDistributionObjectY_Uint64t(mt64R);
            UnsignedInt RandomPosZ = UniformDistributionObjectZ_Uint64t(mt64R);

            UnsignedInt RandomSizeOfParticle = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

            ChainIdInt RandomChainId = UniformDistributionObjectTypeOfParticle_Uint64t(mt64R);

            FilledVoxelsForRandomParticle.clear();

            if (RandomPosX + RandomSizeOfParticle < XStartParam + XSizeParam && RandomPosY + RandomSizeOfParticle < YStartParam + YSizeParam && RandomPosZ + RandomSizeOfParticle < ZStartParam + ZSizeParam)
                for (UnsignedInt SpaceXP = RandomPosX; SpaceXP < RandomPosX + RandomSizeOfParticle; SpaceXP++)
                    for (UnsignedInt SpaceYP = RandomPosY; SpaceYP < RandomPosY + RandomSizeOfParticle; SpaceYP++)
                        for (UnsignedInt SpaceZP = RandomPosZ; SpaceZP < RandomPosZ + RandomSizeOfParticle; SpaceZP++)
                        {
                            if (GetSpaceVoxel(SpaceXP, SpaceYP, SpaceZP).EntityId == 0)
                            {
                                FilledVoxelsForRandomParticle.emplace_back(SpaceXP, SpaceYP, SpaceZP);
                                GetSpaceVoxel(SpaceXP, SpaceYP, SpaceZP) = {static_cast<EntityIdInt>(CellEngineConfigDataObject.DNAIdentifier), RandomChainId };
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

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusion(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
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
                    if (GetSpaceVoxel(VoxelForParticle.X + ShiftX, VoxelForParticle.Y + ShiftY, VoxelForParticle.Z + ShiftZ).EntityId == 0 && VoxelForParticle.X + ShiftX >= XStartParam && VoxelForParticle.X + ShiftX < XStartParam + XSizeParam && VoxelForParticle.Y + ShiftY >= YStartParam && VoxelForParticle.Y + ShiftY < YStartParam + YSizeParam && VoxelForParticle.Z + ShiftZ >= ZStartParam && VoxelForParticle.Z + ShiftZ < ZStartParam + ZSizeParam)
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
