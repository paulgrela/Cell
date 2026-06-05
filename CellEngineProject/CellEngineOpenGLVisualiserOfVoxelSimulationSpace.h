#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include "../Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineSimulationSpace.h"

constexpr UnsignedInt MaxNumberOfSectors = 16;

class CellEngineOpenGLVisualiserOfVoxelSimulationSpace : public CellEngineOpenGLVisualiser
{
private:
    CellEngineAtom TempAtomObjectInSectors[MaxNumberOfSectors][MaxNumberOfSectors][MaxNumberOfSectors] = {};
    std::uint32_t AtomOffsetInSectors[MaxNumberOfSectors][MaxNumberOfSectors][MaxNumberOfSectors] = {};
    std::vector<GPUParticle> GPUParticlesInSectors[MaxNumberOfSectors][MaxNumberOfSectors][MaxNumberOfSectors];
    std::vector<GPUAtom> GPUAtomsInSectors[MaxNumberOfSectors][MaxNumberOfSectors][MaxNumberOfSectors];
    std::vector<GPUAtomLocal> GPUAtomsLocalInSectors[MaxNumberOfSectors][MaxNumberOfSectors][MaxNumberOfSectors];

public:
    enum class VoxelSpaceDrawingTypes : UnsignedInt
    {
        DrawVoxelSpaceFull = 1,
        DrawVoxelSpaceSelected = 2
    };
public:
    VoxelSpaceDrawingTypes SpaceDrawingType = VoxelSpaceDrawingTypes::DrawVoxelSpaceFull;
    bool DrawEmptyVoxels = false;
private:
    UnsignedInt SelectionStartXPos{}, SelectionStartYPos{}, SelectionStartZPos{};
    UnsignedInt SelectionStepX{}, SelectionStepY{}, SelectionStepZ{};
    UnsignedInt SelectionSizeX{}, SelectionSizeY{}, SelectionSizeZ{};
public:
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetStartPositions();
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetSteps();
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetSizes();
    void SetVoxelSpaceSelection(UnsignedInt SelectionStartXParam, UnsignedInt SelectionStartYParam, UnsignedInt SelectionStartZParam, UnsignedInt SelectionStepXParam, UnsignedInt SelectionStepYParam, UnsignedInt SelectionStepZParam, UnsignedInt SelectionSizeXParam, UnsignedInt SelectionSizeYParam, UnsignedInt SelectionSizeZParam);
public:
    CellEngineOpenGLVisualiserOfVoxelSimulationSpace();
private:
    struct TemporaryRenderedVoxel
    {
        CellEngineAtom CellEngineAtomObject;
        ParticlesDetailedContainer<Particle>::iterator CellEngineParticleObject;
        UnsignedInt X, Y, Z;
    };
    UnsignedInt SaveXMousePosition{}, SaveYMousePosition{}, SaveZMousePosition{};
public:
    void SaveVoxelPositionChosenByMouse();
private:
    void SetSaveXYZPositions(UnsignedInt SaveXParam, UnsignedInt SaveYParam, UnsignedInt SaveZParam);
private:
    inline void ConvertAtomPosToGraphicCoordinate(CellEngineAtom& CellEngineAtomObjectParam, UnsignedInt StartXParam, UnsignedInt StartYParam, UnsignedInt StartZParam, UnsignedInt SpaceXParam, UnsignedInt SpaceYParam, UnsignedInt SpaceZParam, UnsignedInt SizeXParam, UnsignedInt SizeYParam, UnsignedInt SizeZParam) const;
    static inline void SetParticleParametersToDraw(CellEngineAtom& TempAtomObject, Particle& ParticleObject);
private:
    void GenerateVoxelsForGPUParallel(UnsignedInt PosXParam, UnsignedInt XEndParam, UnsignedInt PosYParam, UnsignedInt YEndParam, UnsignedInt PosZParam, UnsignedInt ZEndParam, UnsignedInt MainPosX, UnsignedInt MainPosY, UnsignedInt MainPosZ, SimulationSpaceVoxel SimulationSpaceVoxelObject, SimulationSpaceVoxel LastSimulationSpaceVoxel, UnsignedInt& AtomsCounter, const CellEngineAtom& TempAtomObject, const Particle& ParticleObject);
    void GenerateVoxelsForGPU(SimulationSpaceVoxel SimulationSpaceVoxelObject, SimulationSpaceVoxel LastSimulationSpaceVoxel, UnsignedInt& AtomsCounter, const CellEngineAtom& TempAtomObject, const Particle& ParticleObject);
    void RenderSelectedSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, CellEngineAtom& TempAtomObject);
protected:
    void RenderSpace(const vmath::mat4& ViewMatrix) override;
protected:
    void FindAndDrawAllBondsBetweenAtoms(const vmath::mat4& ViewMatrix) override
    {
    }
public:
    inline void DrawChosenAtomUsingStencilBuffer1(GLuint ChosenParticleCenterIndex) override;
protected:
    void GetStartCenterPoint() override;
};

#endif

#endif