#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"

constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimension = 1024;
constexpr UnsignedInt NumberOfVoxelsInEachDimension1 = 512;
constexpr UnsignedInt NumberOfVoxelsInEachDimension2 = 256;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInDimensionX = NumberOfVoxelSimulationSpaceInEachDimension;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInDimensionY = NumberOfVoxelSimulationSpaceInEachDimension;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInDimensionZ = NumberOfVoxelSimulationSpaceInEachDimension;

//TYLKO 2 PONIZSZE BEZ KUBELKOW Z PLIKU KONFIGURACYJNEGO

//using VoxelSimulationSpaceType1024 = UnsignedInt[NumberOfVoxelSimulationSpaceInDimensionX][NumberOfVoxelSimulationSpaceInDimensionY][NumberOfVoxelSimulationSpaceInDimensionZ];
//using VoxelSimulationSpaceType512 = UnsignedInt[512][512][512];
//using VoxelSimulationSpaceType256 = UnsignedInt[256][256][256];
//
////TE 2 z kubelkami
//using T128 = uint64_t[128][128][128];
//using T64 = uint64_t[64][64][64];
//
//inline void Test1()
//{
//    auto a = new T256;
//    a[1][2][3] = 2;
//}


class CellEngineVoxelSimulationSpace
{
public:
    UnsignedInt Space[NumberOfVoxelSimulationSpaceInDimensionX][NumberOfVoxelSimulationSpaceInDimensionY][NumberOfVoxelSimulationSpaceInDimensionZ]{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam)
    {
        return static_cast<float>(static_cast<Int>(CoordinateParam) - (static_cast<Int>(NumberOfVoxelSimulationSpaceInEachDimension / 2))) * 4;
    };
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam)
    {
        return static_cast<UnsignedInt>(round(CoordinateParam / 4.0)) + (NumberOfVoxelSimulationSpaceInEachDimension / 2);
    };
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinateX(UnsignedInt CoordinateParam)
    {
        return ConvertToGraphicsCoordinate(CoordinateParam);
    };
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinateX(double CoordinateParam)
    {
        return ConvertToSpaceCoordinate(CoordinateParam);
    };
    [[nodiscard]] static float ConvertToGraphicsCoordinateY(UnsignedInt CoordinateParam)
    {
        return ConvertToGraphicsCoordinate(CoordinateParam);
    };
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinateY(double CoordinateParam)
    {
        return ConvertToSpaceCoordinate(CoordinateParam);
    };
    [[nodiscard]] static float ConvertToGraphicsCoordinateZ(UnsignedInt CoordinateParam)
    {
        return ConvertToGraphicsCoordinate(CoordinateParam);
    };
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinateZ(double CoordinateParam)
    {
        return ConvertToSpaceCoordinate(CoordinateParam);
    };
public:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
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
        for (auto& SelectedX : Space)
            for (auto& SelectedXY : SelectedX)
                for (UnsignedInt& SelectedXYZ : SelectedXY)
                {
                    if (SelectedXYZ != 0)
                        SumOfNotEmptyVoxels++;
                }
    }
public:
    void SetAtomInVoxelSpace(const CellEngineAtom& AppliedAtom)
    {
        UnsignedInt SpaceX = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinateX(AppliedAtom.X);
        UnsignedInt SpaceY = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinateY(AppliedAtom.Y);
        UnsignedInt SpaceZ = CellEngineVoxelSimulationSpace::ConvertToSpaceCoordinateZ(AppliedAtom.Z);

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
                for (UnsignedInt& SelectedXYZ : SelectedXY)
                    SelectedXYZ = 0;
    }
};

#endif
