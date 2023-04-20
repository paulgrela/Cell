
#include "CellEngineAtom.h"
#include "CellEngineVoxelSimulationSpace.h"

using namespace std;

[[nodiscard]] float CellEngineVoxelSimulationSpace::ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
{
    return static_cast<float>(static_cast<Int>(CoordinateParam) - (static_cast<Int>(CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2))) * 4;
};
[[nodiscard]] UnsignedInt CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(double CoordinateParam)
{
    return static_cast<UnsignedInt>(round(CoordinateParam / 4.0)) + (CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2);
};

CellEngineVoxelSimulationSpace::CellEngineVoxelSimulationSpace()
{
    try
    {
        SetStartValuesForSpaceMinMax();

        for (auto& SelectedX : Space)
            for (auto& SelectedXY : SelectedX)
                for (SimulationSpaceVoxel& SelectedXYZ : SelectedXY)
                    SelectedXYZ.EntityId = SelectedXYZ.ChainId = 0;
    }
    CATCH("execution of constructor of voxel simulation space")
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
        XMin = std::min(SpaceX, XMin);
        XMax = std::max(SpaceX, XMax);
        YMin = std::min(SpaceY, YMin);
        YMax = std::max(SpaceY, YMax);
        ZMin = std::min(SpaceZ, ZMin);
        ZMax = std::max(SpaceZ, ZMax);
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
                    if (Space[PosX][PosY][PosZ].EntityId != 0)
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

        if (Space[SpaceX][SpaceY][SpaceZ].EntityId == 0)
        {
            Space[SpaceX][SpaceY][SpaceZ].EntityId = AppliedAtom.EntityId;
            if (CellEngineConfigDataObject.IsDNAorRNA(AppliedAtom.EntityId) == true)
                Space[SpaceX][SpaceY][SpaceZ].ChainId = stoi(std::string(AppliedAtom.Chain).substr(2,2));
        }
    }
    CATCH("setting atom in voxel simulation space")
}

void CellEngineVoxelSimulationSpace::AddParticleKind(const ParticleKind& ParticleParam)
{
    Particles.emplace_back(ParticleParam);
}

void CellEngineVoxelSimulationSpace::GenerateRandomParticlesInSelectedSpace(const UnsignedInt NumberOfRandomParticles, const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {
        AddParticleKind({1, "A", "A1", 0});
        AddParticleKind({2, "B", "B1", 0});
        AddParticleKind({3, "C", "C1", 0});
        AddParticleKind({4, "D", "D1", 0});

        std::vector<vector3_64> FilledVoxelsForRandomParticle;

        std::mt19937_64 mt64R{ std::random_device{}() };

        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectSizeOfParticle_Uint64t(1, 2);
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectTypeOfParticle_Uint64t(1, 4);
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectX_Uint64t(XStartParam, XStartParam + XSizeParam);
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectY_Uint64t(YStartParam, YStartParam + YSizeParam);
        std::uniform_int_distribution<UnsignedInt> UniformDistributionObjectZ_Uint64t(ZStartParam, ZStartParam + ZSizeParam);

        for (UnsignedInt SpaceXP = XStartParam; SpaceXP < XStartParam + XSizeParam; SpaceXP += XStepParam)
            for (UnsignedInt SpaceYP = YStartParam; SpaceYP < YStartParam + YSizeParam; SpaceYP += YStepParam)
                for (UnsignedInt SpaceZP = ZStartParam; SpaceZP < ZStartParam + ZSizeParam; SpaceZP += ZStepParam)
                    Space[SpaceXP][SpaceYP][SpaceZP].EntityId = Space[SpaceXP][SpaceYP][SpaceZP].ChainId = 0;

        for (UnsignedInt ParticleNum = 1;  ParticleNum <= NumberOfRandomParticles; ParticleNum++)
        {
            UnsignedInt RandomPosX = UniformDistributionObjectX_Uint64t(mt64R);
            UnsignedInt RandomPosY = UniformDistributionObjectY_Uint64t(mt64R);
            UnsignedInt RandomPosZ = UniformDistributionObjectZ_Uint64t(mt64R);

            UnsignedInt RandomSizeOfParticle = UniformDistributionObjectSizeOfParticle_Uint64t(mt64R);

            UnsignedInt RandomChainId = UniformDistributionObjectTypeOfParticle_Uint64t(mt64R);

            FilledVoxelsForRandomParticle.clear();

            for (UnsignedInt SpaceXP = RandomPosX; SpaceXP < RandomPosX + RandomSizeOfParticle; SpaceXP++)
                for (UnsignedInt SpaceYP = RandomPosY; SpaceYP < RandomPosY + RandomSizeOfParticle; SpaceYP++)
                    for (UnsignedInt SpaceZP = RandomPosZ; SpaceZP < RandomPosZ + RandomSizeOfParticle; SpaceZP++)
                    {
                        if (Space[SpaceXP][SpaceXP][SpaceXP].EntityId == 0)
                        {
                            FilledVoxelsForRandomParticle.emplace_back(SpaceXP, SpaceYP, SpaceZP);
                            Space[SpaceXP][SpaceYP][SpaceZP].EntityId = CellEngineConfigDataObject.DNAIdentifier;
                            Space[SpaceXP][SpaceYP][SpaceZP].ChainId = RandomChainId;
                        }
                        else
                        {
                            for (auto& VoxelForRandomParticle : FilledVoxelsForRandomParticle)
                                Space[VoxelForRandomParticle.X][VoxelForRandomParticle.Y][VoxelForRandomParticle.Z].EntityId = Space[VoxelForRandomParticle.X][VoxelForRandomParticle.Y][VoxelForRandomParticle.Z].ChainId = 0;
                            goto NextRandomParticleOutsideLoopLabel;
                        }
                    }
            NextRandomParticleOutsideLoopLabel:;

            if (FilledVoxelsForRandomParticle.empty() == false)
                Particles[RandomChainId - 1].ParticlesObjects.emplace_back(RandomPosX, RandomPosY, RandomPosZ, FilledVoxelsForRandomParticle);
        }
    }
    CATCH("generating random particles in selected space")
}

void CellEngineVoxelSimulationSpace::GenerateOneStepOfDiffusion(const UnsignedInt XStartParam, const UnsignedInt YStartParam, const UnsignedInt ZStartParam, const UnsignedInt XStepParam, const UnsignedInt YStepParam, const UnsignedInt ZStepParam, const UnsignedInt XSizeParam, UnsignedInt YSizeParam, const UnsignedInt ZSizeParam)
{
    try
    {

    }
    CATCH("generating one step of diffusion")
}
