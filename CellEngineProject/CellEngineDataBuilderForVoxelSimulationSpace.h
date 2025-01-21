
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineParticlesDataFile.h"
#include "CellEngineSimulationParallelExecutionManager.h"

class CellEngineDataBuilderForVoxelSimulationSpace : public CellEngineParticlesDataFile
{
protected:
    void SetStartValues() override
    {
        GetMemoryForParticlesInSectors();

        CellEngineVoxelSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpace>(Particles, true, 0, CurrentThreadPosType{ 0, 0, 0 });

        CellEngineSimulationParallelExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineVoxelSimulationSpace>(CellEngineSimulationSpaceForThreadsObjectsPointer, Particles);
    }
protected:
    UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) override
    {
        return CellEngineVoxelSimulationSpaceObjectPointer->AddNewParticle(ParticleObjectParam);
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, const UniqueIdInt ParticleIndex) override
    {
        CellEngineVoxelSimulationSpaceObjectPointer->SetAtomInVoxelSimulationSpace(ParticleIndex, AppliedAtom);
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
        CellEngineVoxelSimulationSpaceObjectPointer->PreprocessData<UnsignedInt>(Update);
    }
protected:
    void PrintStatistics() override
    {
        LoggersManagerObject.Log(CellEngineVoxelSimulationSpaceObjectPointer->PrintSpaceMinMaxValues());
        CellEngineVoxelSimulationSpaceObjectPointer->CountStatisticsOfVoxelSimulationSpace();
        LoggersManagerObject.Log(STREAM("Sum Of Not Empty Voxels = " << CellEngineVoxelSimulationSpaceObjectPointer->GetSumOfNotEmptyVoxels()));
    }
};

#endif
