#pragma once

#ifndef CELL_ENGINE_PROJECT_BUILD_PARTICLES_DATA_OPERATIONS_H
#define CELL_ENGINE_PROJECT_BUILD_PARTICLES_DATA_OPERATIONS_H

#include "CellEngineDataFile.h"

class CellEngineBuildParticlesDataOperations
{
protected:
    virtual void SetStartValues() = 0;
    virtual UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) = 0;
    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, UniqueIdInt ParticleIndex) = 0;
    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
    virtual void PreprocessData() = 0;
    virtual void PrintStatistics() = 0;
};

#endif
