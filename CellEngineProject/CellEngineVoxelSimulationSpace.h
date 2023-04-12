#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"

constexpr UnsignedIntType NumberOfVoxelsInOneDimension = 1024;
//TYLKO 2 PONIZSZE BEZ KUBELKOWZ Z PLIKU KONFIGURACYJNEGO
using T1024 = uint64_t[1024][1024][1024]; //voxel 0.5 nm
using T512 = uint64_t[512][512][512]; //voxel 1 nm
using T256 = uint64_t[256][256][256]; // voxel 2 nm

//TE 2 z kubelkami
using T128 = uint64_t[128][128][128]; //voxel 4 nm
using T64 = uint64_t[64][64][64]; // voxel 8 nm

inline void Test1()
{
    auto a = new T256;
    a[1][2][3] = 2;
}


class CellEngineVoxelSimulationSpace
{
public:
    uint64_t Space[1024][1024][1024]{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedIntType CoordinateParam)
    {
        return static_cast<float>(static_cast<IntType>(CoordinateParam) - 512) * 4;
    };
    [[nodiscard]] static UnsignedIntType ConvertToSpaceCoordinate(double CoordinateParam)
    {
        return static_cast<UnsignedIntType>(round(CoordinateParam / 4.0)) + 512;
    };
public:
    UnsignedIntType XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    void SetStartValuesForSpaceMinMax()
    {
        XMin = 10000;
        XMax = 0;
        YMin = 10000;
        YMax = 0;
        ZMin = 10000;
        ZMax = 0;
    }
public:
    void CompareAndGetSpaceMinMax(const UnsignedIntType SpaceX, const UnsignedIntType SpaceY, const UnsignedIntType SpaceZ)
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
        UnsignedIntType SpaceX = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.X);
        UnsignedIntType SpaceY = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.Y);
        UnsignedIntType SpaceZ = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinate(AppliedAtom.Z);

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
