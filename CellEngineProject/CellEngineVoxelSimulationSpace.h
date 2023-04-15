#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineConfigData.h"

constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst = 1024;

struct SimulationSpaceVoxel
{
    EntityIdInt EntityId;
    ChainIdInt ChainId;
};

class CellEngineVoxelSimulationSpace
{
public:
    SimulationSpaceVoxel Space[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst]{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
    {
        return static_cast<float>(static_cast<Int>(CoordinateParam) - (static_cast<Int>(CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2))) * 4;
    };
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam)
    {
        return static_cast<UnsignedInt>(round(CoordinateParam / 4.0)) + (CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension / 2);
    };
public:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    void SetStartValuesForSpaceMinMax()
    {
        XMin = YMin = ZMin = 10000;
        XMax = YMax = ZMax = 0;
    }
public:
    void CompareAndGetSpaceMinMax(const UnsignedInt SpaceX, const UnsignedInt SpaceY, const UnsignedInt SpaceZ)
    {
        XMin = std::min(SpaceX, XMin);
        XMax = std::max(SpaceX, XMax);
        YMin = std::min(SpaceY, YMin);
        YMax = std::max(SpaceY, YMax);
        ZMin = std::min(SpaceZ, YMin);
        ZMax = std::max(SpaceZ, YMax);
    }
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const
    {
        std::stringstream ss;
        ss << "CELL SPACE LIMITS PARAMETERS [ Xmin = " << std::to_string(XMin) << " ][ Xmax = " << std::to_string(XMax) << " ][ Ymin = " << std::to_string(YMin) << " ][ Ymax = " << std::to_string(YMax) << " ][ Zmin = " << std::to_string(ZMin) << " ][ Zmax = " << std::to_string(XMax) << " ] " << std::endl;
        return ss;
    }
public:
    UnsignedInt SumOfNotEmptyVoxels{};
public:
    void CountStatisticsOfSpace()
    {
        for(UnsignedInt PosX = 0; PosX < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosX++)
            for(UnsignedInt PosY = 0; PosY < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosY++)
                for(UnsignedInt PosZ = 0; PosZ < CellEngineConfigDataObject.NumberOfVoxelSimulationSpaceInEachDimension; PosZ++)
                {
                    if (Space[PosX][PosY][PosZ].EntityId != 0)
                        SumOfNotEmptyVoxels++;
                }

    }
public:
    void SetAtomInVoxelSpace(const CellEngineAtom& AppliedAtom)
    {
        UnsignedInt SpaceX = ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedInt SpaceY = ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedInt SpaceZ = ConvertToSpaceCoordinate(AppliedAtom.Z);

        CompareAndGetSpaceMinMax(SpaceX, SpaceY, SpaceZ);

        if (Space[SpaceX][SpaceY][SpaceZ].EntityId == 0)
        {
            Space[SpaceX][SpaceY][SpaceZ].EntityId = AppliedAtom.EntityId;
            if (CellEngineConfigDataObject.IsDNAorRNA(AppliedAtom.EntityId) == true)
                Space[SpaceX][SpaceY][SpaceZ].ChainId = stoi(std::string(AppliedAtom.Chain).substr(2,2));
        }
    }

public:
    CellEngineVoxelSimulationSpace()
    {
        SetStartValuesForSpaceMinMax();

        for (auto& SelectedX : Space)
            for (auto& SelectedXY : SelectedX)
                for (SimulationSpaceVoxel& SelectedXYZ : SelectedXY)
                    SelectedXYZ.EntityId = SelectedXYZ.ChainId = 0;
    }
};

#endif
