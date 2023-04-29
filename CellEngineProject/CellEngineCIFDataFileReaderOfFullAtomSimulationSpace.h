#pragma once

#ifndef CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H
#define CELL_ENGINE_FULL_ATOM_SIMULATION_SPACE_CIF_DATA_FILE_READER_H

#include "CellEngineConfigData.h"
#include "CellEngineCIFDataFile.h"

class CellEngineCIFDataFileReaderOfFullAtomSimulationSpace : public CellEngineCIFDataFile
{
protected:
    void SetStartValues() override
    {
    }
protected:
    void SetParticleKindData(const EntityIdInt EntityId, const ChainIdInt ChainId) override
    {
    }
protected:
    void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom) override
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
    void PrintStatistics() override
    {
    }
};

#endif
