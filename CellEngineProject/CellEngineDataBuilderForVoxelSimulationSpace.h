
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineParticlesDataFile.h"

class CellEngineDataBuilderForVoxelSimulationSpace : public CellEngineParticlesDataFile
{
protected:
    void SetStartValues() override
    {
        CellEngineVoxelSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpaceD1>(Particles);
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
        CellEngineVoxelSimulationSpaceObjectPointer->PreprocessData(Update);
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
