
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"
#include "CellEngineSimulationSpace.h"

class CellEngineDataFile
{
public:
    CellEngineSimulationSpace CellEngineSimulationSpaceObject;
public:
    bool FilmOfStructuresActive = false;
protected:
    std::vector<std::vector<CellEngineAtom>> ParticlesCenters;
    std::vector<std::vector<CellEngineAtom>> AllAtoms;
public:
    virtual void ReadDataFromFile() = 0;
public:
    static vmath::vec3 GetCenter(const std::vector<CellEngineAtom>& AtomsParam)
    {
        vmath::vec3 Center(0.0, 0.0, 0.0);

        try
        {
            for (const CellEngineAtom& AtomObject : AtomsParam)
                Center += AtomObject.Position();
            Center /= static_cast<float>(AtomsParam.size());
        }
        CATCH_AND_THROW("counting mass center")

        return Center;
    }
public:
    [[nodiscard]] std::vector<std::vector<CellEngineAtom>>& GetAllAtoms()
    {
        return AllAtoms;
    }
    [[nodiscard]] std::vector<CellEngineAtom>& GetParticlesCenters()
    {
        return ParticlesCenters[CellEngineConfigDataObject.ChosenStructureIndex];
    }
    [[nodiscard]] UnsignedIntType GetNumberOfStructures()
    {
        return ParticlesCenters.size();
    }
public:
    void ShowNextStructureFromActiveFilm()
    {
        try
        {
            if (FilmOfStructuresActive == true)
                ShowNextStructure();
        }
        CATCH("showing next structure from film")
    }
public:
    void StartFilmOfStructures()
    {
        try
        {
            std::lock_guard<std::mutex> LockGuardObject{ChosenStructureMutexObject};

            CellEngineConfigDataObject.ChosenStructureIndex = 0;
            FilmOfStructuresActive = true;
        }
        CATCH("starting film of structures")
    }
public:
    void StopFilmOfStructures()
    {
        try
        {
            FilmOfStructuresActive = false;
        }
        CATCH("starting film of structures")
    }
public:
    void ShowNextStructure()
    {
        try
        {
            std::lock_guard<std::mutex> LockGuardObject{ChosenStructureMutexObject};

            if (CellEngineConfigDataObject.ChosenStructureIndex < GetNumberOfStructures() - 1)
                CellEngineConfigDataObject.ChosenStructureIndex++;
            else
                FilmOfStructuresActive = false;
        }
        CATCH("showing next structure")
    }
public:
    static void ShowPrevStructure()
    {
        try
        {
            std::lock_guard<std::mutex> LockGuardObject{ChosenStructureMutexObject};

            if (CellEngineConfigDataObject.ChosenStructureIndex > 0)
                CellEngineConfigDataObject.ChosenStructureIndex--;
        }
        CATCH("showing previous structure")
    }
public:
    static inline std::mutex ChosenStructureMutexObject;
};

inline std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;

#endif
