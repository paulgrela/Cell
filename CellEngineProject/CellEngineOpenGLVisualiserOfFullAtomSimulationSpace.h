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
protected:
    void RenderSpace1(const vmath::mat4& ViewMatrix, uint32_t& ParticlesOffsetTotal, uint32_t& AtomOffsetTotal, uint32_t& AtomLocalOffsetTotal, vector<GPUParticle>& GPUParticles, vector<GPUAtom>& GPUAtoms, vector<GPUAtomLocal>& GPUAtomsLocal);
    void RenderSpace2(const vmath::mat4& ViewMatrix, uint32_t& ParticlesOffsetTotal, uint32_t& AtomOffsetTotal, uint32_t& AtomLocalOffsetTotal);
    void RenderSpace(const vmath::mat4& ViewMatrix) override;
protected:
    void FindAndDrawAllBondsBetweenAtoms(const vmath::mat4& ViewMatrix) override;
public:
    inline void DrawChosenAtomUsingStencilBuffer1(GLuint ChosenParticleCenterIndex) override;
protected:
    void GetStartCenterPoint() override;
};

#endif

#endif