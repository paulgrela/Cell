#pragma once

#ifndef CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_VOXEL_SIMULATION_SPACE_VISUALISER_H

#include "CellEngineConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include <string>

#include <sb7.h>
#include <tuple>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"
#include "CellEngineSimulationSpace.h"

class CellEngineOpenGLVisualiserOfVoxelSimulationSpace : public CellEngineOpenGLVisualiser
{
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
    void RenderSelectedSpace(UnsignedInt XStartParam, UnsignedInt YStartParam, UnsignedInt ZStartParam, UnsignedInt XStepParam, UnsignedInt YStepParam, UnsignedInt ZStepParam, UnsignedInt XSizeParam, UnsignedInt YSizeParam, UnsignedInt ZSizeParam, UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, CellEngineAtom& TempAtomObject, std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList, UnsignedInt StencilBufferLoopCounter);
    void RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, vmath::mat4& ViewMatrix) override;
private:
    inline void DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<TemporaryRenderedVoxel>& TemporaryRenderedVoxelsList);
public:
    void GetStartCenterPoint() override;
public:
    static inline std::recursive_mutex RenderMenuAndVoxelSimulationSpaceMutexObject;
};

#endif

#endif