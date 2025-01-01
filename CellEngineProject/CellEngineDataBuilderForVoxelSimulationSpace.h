
#ifndef CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_VOXEL_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineParticlesDataFile.h"

class CellEngineDataBuilderForVoxelSimulationSpace : public CellEngineParticlesDataFile
{
protected:
    void SetStartValues() override
    {
        CellEngineVoxelSimulationSpaceObjectPointer = std::make_unique<CellEngineVoxelSimulationSpace>(Particles, true, 0, CurrentThreadPosType{ 0, 0, 0 });

        UnsignedInt ThreadIndexPos = 1;

        CellEngineSimulationSpaceForThreadsObjectsPointer.clear();
        CellEngineSimulationSpaceForThreadsObjectsPointer.resize(CellEngineConfigDataObject.NumberOfXThreadsInSimulation);
        UnsignedInt ThreadXPos = 1;
        for (auto& ThreadLocalParticlesInProximityXPos : CellEngineSimulationSpaceForThreadsObjectsPointer)
        {
            ThreadLocalParticlesInProximityXPos.resize(CellEngineConfigDataObject.NumberOfYThreadsInSimulation);
            UnsignedInt ThreadYPos = 1;
            for (auto& ThreadLocalParticlesInProximityYPos : ThreadLocalParticlesInProximityXPos)
            {
                ThreadLocalParticlesInProximityYPos.resize(CellEngineConfigDataObject.NumberOfZThreadsInSimulation);
                UnsignedInt ThreadZPos = 1;
                for (auto& ThreadLocalParticlesInProximityZPos : ThreadLocalParticlesInProximityYPos)
                {
                    LoggersManagerObject.Log(STREAM("THREAD INDEXES = " << ThreadIndexPos << " (" << ThreadXPos << ", " << ThreadYPos << ", " << ThreadZPos << ")"));

                    ThreadLocalParticlesInProximityZPos = std::make_unique<CellEngineVoxelSimulationSpace>(Particles, false, ThreadIndexPos, CurrentThreadPosType{ ThreadXPos, ThreadYPos, ThreadZPos });
                    ThreadIndexPos++;
                    ThreadZPos++;
                }
                ThreadYPos++;
            }
            ThreadXPos++;
        }
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
