#pragma once

#ifndef CELL_ENGINE_OPENGL_FULL_ATOM_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_FULL_ATOM_SIMULATION_SPACE_VISUALISER_H

#include "../Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include <string>

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"

class CellEngineOpenGLVisualiserOfFullAtomSimulationSpace : public CellEngineOpenGLVisualiser
{
public:
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetStartPositions();
    std::tuple<UnsignedInt, UnsignedInt, UnsignedInt> GetSizes();
private:
    UnsignedInt SelectionStartXPos{}, SelectionStartYPos{}, SelectionStartZPos{};
    UnsignedInt SelectionSizeX{}, SelectionSizeY{}, SelectionSizeZ{};
private:
    void RenderSpace0(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix);
    void RenderSpace1(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, uint32_t AtomOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms);
    void RenderSpace2(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, uint32_t AtomOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms);
    void RenderSpace3(UnsignedInt& NumberOfAllRenderedAtoms, const vmath::mat4& ViewMatrix, uint32_t AtomOffsetTotal, uint32_t ParticlesOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms);

    void RenderSpace(UnsignedInt& NumberOfAllRenderedAtoms, UnsignedInt& NumberOfFoundParticlesCenterToBeRenderedInAtomDetails, const vmath::mat4& ViewMatrix) override;
    inline void DrawChosenAtomUsingStencilBuffer(const vmath::mat4& ViewMatrix, const GLuint* PartOfStencilBufferIndex, UnsignedInt& NumberOfAllRenderedAtoms, const std::vector<std::tuple<UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt, UnsignedInt>>& TemporaryRenderedAtomsList);
public:
    void GetStartCenterPoint() override;
public:
    static inline std::mutex RenderMenuAndFullAtomSimulationSpaceMutexObject;

    static inline std::mutex PrepareParticlesFullAtomSimulationSpaceMutexObject;
};

#endif

#endif