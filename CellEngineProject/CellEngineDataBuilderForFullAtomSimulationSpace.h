#pragma once

#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineParticlesDataFile.h"
#include "CellEngineSimulationParallelExecutionManager.h"

class CellEngineDataBuilderForFullAtomSimulationSpace : public CellEngineParticlesDataFile
{
protected:
    void SetStartValues() override
    {
        GetMemoryForParticlesInSectors();

        CellEngineFullAtomSimulationSpaceObjectPointer = std::make_unique<CellEngineFullAtomSimulationSpace>(Particles, true, 0, CurrentThreadPosType{ 0, 0, 0 });

        CellEngineSimulationParallelExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineFullAtomSimulationSpace>(CellEngineSimulationSpaceForThreadsObjectsPointer, Particles);
    }
protected:
    UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) override
    {
        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false)
            return 0;
        else
            return CellEngineFullAtomSimulationSpaceObjectPointer->AddNewParticle(ParticleObjectParam);
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, const UniqueIdInt ParticleIndex) override
    {
        LocalCellEngineAllAtomsObject.emplace_back(AppliedAtom);
    }
protected:
    void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) override
    {
        LocalCellEngineParticlesCentersObject.emplace_back(LocalCellEngineAllAtomsObject.front().X, LocalCellEngineAllAtomsObject.front().Y, LocalCellEngineAllAtomsObject.front().Z, AllAtoms.size(), LocalCellEngineParticlesCentersObject.size(), LocalCellEngineAllAtomsObject.front().EntityId, LocalCellEngineAllAtomsObject.front().Name, LocalCellEngineAllAtomsObject.front().ResName, LocalCellEngineAllAtomsObject.front().Chain, LocalCellEngineAllAtomsObject.front().ParticleColor);
        AllAtoms.emplace_back(LocalCellEngineAllAtomsObject);
    }
protected:
    void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) override
    {
        ParticlesCenters.emplace_back(LocalCellEngineParticlesCentersObject);
    }
protected:
    void PreprocessData(bool Update) override
    {
    }
protected:
    void PrintStatistics() override
    {
    }
};

#endif
