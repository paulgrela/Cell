
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFileReader.h"

class CellEngineDataBuilderForVoxelSimulationSpace : virtual public CellEngineDataFile
{
protected:
    void SetStartValues() override
    {
        CellEngineVoxelSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpace>();
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
    void PreprocessData() override
    {
        CellEngineVoxelSimulationSpaceObjectPointer->PreprocessData(true);
    }
protected:
    void PrintStatistics() override
    {
        LoggersManagerObject.Log(CellEngineVoxelSimulationSpaceObjectPointer->PrintSpaceMinMaxValues());
        CellEngineVoxelSimulationSpaceObjectPointer->CountStatisticsOfVoxelSimulationSpace();
        LoggersManagerObject.Log(STREAM("Sum Of Not Empty Voxels = " << CellEngineVoxelSimulationSpaceObjectPointer->SumOfNotEmptyVoxels));
    }
};

class CellEngineCIFDataFileReaderOfVoxelSimulationSpace : public CellEngineCIFDataFileReader, public CellEngineDataBuilderForVoxelSimulationSpace
{
};

class CellEngineParticlesDataFileReaderForVoxelSimulationSpace : public CellEngineParticlesDataFileReader, public CellEngineDataBuilderForVoxelSimulationSpace
{
};

#endif
