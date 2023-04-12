
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFile.h"

class CellEngineVoxelSimulationSpaceCIFDataFileReader : public CellEngineCIFDataFile
{
protected:
    void SetStartValues() override
    {
        CellEngineVoxelSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpace>();

        CellEngineVoxelSimulationSpaceObjectPointer->SetStartValuesForSpaceMinMax();
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom) override
    {
        CellEngineVoxelSimulationSpaceObjectPointer->SetAtomInVoxelSpace(AppliedAtom);
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
    void PrintStatistics() override
    {
        LoggersManagerObject.Log(CellEngineVoxelSimulationSpaceObjectPointer->PrintSpaceMinMaxValues());
        CellEngineVoxelSimulationSpaceObjectPointer->CountStatisticsOfSpace();
        LoggersManagerObject.Log(STREAM("Sum Of Not Empty Voxels = " << CellEngineVoxelSimulationSpaceObjectPointer->SumOfNotEmptyVoxels));
    }
};


#endif
