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
    void SetCurrentSectorPos(const SectorPosType& CurrentSectorPos) override
    {
        CellEngineFullAtomSimulationSpaceObjectPointer->SetCurrentSectorPos(CurrentSectorPos);
    }
protected:
    UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) override
    {
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
    }
protected:
    void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) override
    {
    }
protected:
    void PreprocessData(bool Update) override
    {
        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false)
            CellEngineFullAtomSimulationSpaceObjectPointer->PreprocessData<RealType, CellEngineAtom>(&Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, Update);
    }
protected:
    void PrintStatistics() override
    {
    }
};

#endif
