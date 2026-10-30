#pragma once

#ifndef CELL_ENGINE_OPENGL_FULL_ATOM_SIMULATION_SPACE_VISUALISER_H
#define CELL_ENGINE_OPENGL_FULL_ATOM_SIMULATION_SPACE_VISUALISER_H

#include "../Compilation/ConditionalCompilationConstants.h"

#ifdef USE_OPENGL

#include "CellEngineDataFile.h"
#include "CellEngineConfigData.h"
#include "CellEngineOpenGLVisualiser.h"

class CellEngineOpenGLVisualiserOfFullAtomSimulationSpace : public CellEngineOpenGLVisualiser
{
private:
    std::vector<std::vector<std::vector<uint32_t>>> AtomOffsetInSectors = {};
    std::vector<std::vector<std::vector<vector<GPUParticle>>>> GPUParticlesInSectors;
    std::vector<std::vector<std::vector<vector<GPUAtom>>>> GPUAtomsInSectors;
    std::vector<std::vector<std::vector<vector<GPUAtomLocal>>>> GPUAtomsLocalInSectors;
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
public:
    CellEngineOpenGLVisualiserOfFullAtomSimulationSpace()
    {
        AtomOffsetInSectors.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
        GPUParticlesInSectors.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
        GPUAtomsInSectors.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);
        GPUAtomsLocalInSectors.resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInX);

        for (UnsignedInt ParticleSectorXIndex = 0; ParticleSectorXIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInX; ParticleSectorXIndex++)
        {
            AtomOffsetInSectors[ParticleSectorXIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInY);
            GPUParticlesInSectors[ParticleSectorXIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInY);
            GPUAtomsInSectors[ParticleSectorXIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInY);
            GPUAtomsLocalInSectors[ParticleSectorXIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);

            for (UnsignedInt ParticleSectorYIndex = 0; ParticleSectorYIndex < CellEngineConfigDataObject.NumberOfParticlesSectorsInY; ParticleSectorYIndex++)
            {
                AtomOffsetInSectors[ParticleSectorXIndex][ParticleSectorYIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
                GPUParticlesInSectors[ParticleSectorXIndex][ParticleSectorYIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
                GPUAtomsInSectors[ParticleSectorXIndex][ParticleSectorYIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
                GPUAtomsLocalInSectors[ParticleSectorXIndex][ParticleSectorYIndex].resize(CellEngineConfigDataObject.NumberOfParticlesSectorsInZ);
            }
        }
    }
};

#endif

#endif