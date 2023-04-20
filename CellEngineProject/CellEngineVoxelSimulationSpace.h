#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
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
    std::vector<ParticleKind> Particles;
public:
    SimulationSpaceVoxel Space[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst]{};
public:
    UnsignedInt XMin{}, XMax{}, YMin{}, YMax{}, ZMin{}, ZMax{};
public:
    [[nodiscard]] static float ConvertToGraphicsCoordinate(UnsignedInt CoordinateParam);
    [[nodiscard]] static UnsignedInt ConvertToSpaceCoordinate(double CoordinateParam);
public:
    void SetStartValuesForSpaceMinMax();
public:
    void GetMinMaxOfCoordinates(UnsignedInt SpaceX, UnsignedInt SpaceY, UnsignedInt SpaceZ);
public:
    [[nodiscard]] std::stringstream PrintSpaceMinMaxValues() const;
public:
    UnsignedInt SumOfNotEmptyVoxels{};
public:
    void CountStatisticsOfVoxelSimulationSpace();
    void SetAtomInVoxelSimulationSpace(const CellEngineAtom& AppliedAtom);
public:
    void AddParticleKind(const ParticleKind& ParticleParam);
public:
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
    void GenerateOneStepOfDiffusion(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
public:
    CellEngineVoxelSimulationSpace();
};

#endif
