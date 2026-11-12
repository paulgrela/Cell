
#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineParticlesDataFile.h"
#include "CellEngineSimulationParallelTasksExecutionManager.h"

class CellEngineDataBuilderForFullAtomSimulationSpace : public CellEngineParticlesDataFile
{
protected:
    void SetStartValues() override
    {
        GetMemoryForParticlesInSectors();

        CellEngineFullAtomSimulationSpaceObjectPointer = std::make_unique<CellEngineFullAtomSimulationSpace>(Particles, true, 0, ParallelTaskPosType{ .PosX = 0, .PosY = 0, .PosZ = 0 });

        CellEngineSimulationParallelTasksExecutionManager::CreateSimulationSpaceForParallelExecution<CellEngineFullAtomSimulationSpace>(CellEngineSimulationSpaceForThreadsObjectsPointer, Particles);

        CellEngineSimulationParallelTasksExecutionManager::SetProperParallelTasksPositionsAndProperMPIProcessIndexesForEveryParticlesSector(Particles);

        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == true)
            CellEngineFullAtomSimulationSpaceObjectPointer->CreateDataForEveryMPIProcessForParallelExecution();
        else
            CellEngineFullAtomSimulationSpaceObjectPointer->CreateDataForEveryThreadForParallelExecution();
    }
protected:
    void SetCurrentSectorPos(const SectorPosType& CurrentSectorPos) override
    {
        CellEngineFullAtomSimulationSpaceObjectPointer->SetCurrentSectorPos(CurrentSectorPos);
    }
protected:
    UniqueIdUInt64Type AddNewParticle(const Particle& ParticleObjectParam) override
    {
        return CellEngineFullAtomSimulationSpaceObjectPointer->AddNewParticle(ParticleObjectParam);
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, const UniqueIdUInt64Type ParticleUniqueId) override
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
    void PreprocessData(const bool Update) override
    {
        if (CellEngineConfigDataObject.MixedFullAtomWithVoxelSpace == false)
            CellEngineFullAtomSimulationSpaceObjectPointer->PreprocessData<RealType, CellEngineAtom>(&Particle::ListOfAtoms, &ParticleKind::ListOfAtoms, Update);

        if (CellEngineConfigDataObject.FullAtomMPIParallelProcessesExecution == false)
            CellEngineSimulationParallelTasksExecutionManager::MoveParticlesOfEveryThreadToMemoryOfItsNUMANode(Particles);
    }
protected:
    void PrintStatistics() override
    {
        CellEngineFullAtomSimulationSpaceObjectPointer->CountMinMaxValuesOfBordersAtomsPositionsInTheCell();
        LoggersManagerObject.Log(CellEngineFullAtomSimulationSpaceObjectPointer->PrintSpaceMinMaxValues());
    }
};

#endif
