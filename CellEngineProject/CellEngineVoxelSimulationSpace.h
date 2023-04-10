#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"

class CellEngineVoxelSimulationSpace
{
public:
    uint64_t Space[1024][1024][1024]{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(IntType CoordinateParam)
    {
        return static_cast<float>(CoordinateParam - 512) * 4;
    };
    [[nodiscard]] static IntType ConvertToSpaceCoordinate(double CoordinateParam)
    {
        return static_cast<IntType>(round(CoordinateParam / 4.0)) + 512;
    };
public:
    IntType XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    void SetStartValuesForSpaceMinMax()
    {
        XMin = 10000;
        XMax = -10000;
        YMin = 10000;
        YMax = -10000;
        ZMin = 10000;
        ZMax = -10000;
    }
public:
    void CompareAndGetSpaceMinMax(const IntType SpaceX, const IntType SpaceY, const IntType SpaceZ)
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
    UnsignedIntType SumOfNotEmptyVoxels{};
public:
    void CountStatisticsOfSpace()
    {
        for (auto& SelectedX : Space)
            for (auto& SelectedXY : SelectedX)
                for (uint64_t& SelectedXYZ : SelectedXY)
                {
                    if (SelectedXYZ != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
public:
    void SetAtomInVoxelSpace(const CellEngineAtom& AppliedAtom)
    {
        IntType SpaceX = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.X);
        IntType SpaceY = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.Y);
        IntType SpaceZ = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.Z);

        CompareAndGetSpaceMinMax(SpaceX, SpaceY, SpaceZ);

        if (Space[SpaceX][SpaceY][SpaceZ] == 0)
            Space[SpaceX][SpaceY][SpaceZ] = AppliedAtom.EntityId;
    }

public:
    CellEngineVoxelSimulationSpace()
    {
        SetStartValuesForSpaceMinMax();

        for (auto& SelectedX : Space)
            for (auto& SelectedXY : SelectedX)
                for (uint64_t& SelectedXYZ : SelectedXY)
                    SelectedXYZ = 0;
    }
};

#endif
