
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFile.h"

class CellEngineVoxelSimulationSpaceCIFDataFileReader : public CellEngineCIFDataFile
{
protected:
    void SetStartValues() override
    {
        CellEngineSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpace>();

        CellEngineSimulationSpaceObjectPointer->SetStartValuesForSpaceMinMax();
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom) override
    {
        CellEngineSimulationSpaceObjectPointer->SetAtomInVoxelSpace(AppliedAtom);
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
        LoggersManagerObject.Log(CellEngineSimulationSpaceObjectPointer->PrintSpaceMinMaxValues());
        CellEngineSimulationSpaceObjectPointer->CountStatisticsOfSpace();
        LoggersManagerObject.Log(STREAM("Sum Of Not Empty Voxels = " << CellEngineSimulationSpaceObjectPointer->SumOfNotEmptyVoxels));
    }
};


#endif
