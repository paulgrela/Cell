#pragma once

#ifndef CELL_ENGINE_SIMULATION_SPACE_H
#define CELL_ENGINE_SIMULATION_SPACE_H

#include "CellEngineTypes.h"
#include "CellEngineParticle.h"
#include "CellEngineReaction.h"
#include "CellEngineConfigData.h"

//constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst = 1024;
constexpr UnsignedInt NumberOfVoxelSimulationSpaceInEachDimensionMaxConst = 2048;

struct SimulationSpaceVoxel
{
    EntityIdInt EntityId;
    ChainIdInt ChainId;
};

class CellEngineVoxelSimulationSpace
{
private:
    std::mt19937_64 mt64R{ std::random_device{}() };
private:
    std::vector<ParticleKind> Particles;
    std::vector<Reaction> Reactions;
    std::unordered_map<std::string, UnsignedInt> ReactionsIdByString;
private:
    SimulationSpaceVoxel Space[NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst][NumberOfVoxelSimulationSpaceInEachDimensionMaxConst]{};
private:
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
    SimulationSpaceVoxel GetSimulationSpaceVoxel(UnsignedInt X, UnsignedInt Y, UnsignedInt Z);
public:
    void AddParticleKind(const ParticleKind& ParticleParam);
    void AddReaction(const Reaction& ReactionParam);
public:
    void GenerateRandomParticlesInSelectedSpace(UnsignedInt NumberOfRandomParticles, UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
    void GenerateOneStepOfDiffusion(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam);
public:
    CellEngineVoxelSimulationSpace();
};

#endif
