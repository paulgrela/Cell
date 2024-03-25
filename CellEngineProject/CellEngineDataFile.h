
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include <memory>


#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"
#include "CellEngineVoxelSimulationSpace.h"

class CellEngineDataFile
{
public:
    CellEngineDataFile() = default;
public:
    std::unique_ptr<CellEngineVoxelSimulationSpace> CellEngineVoxelSimulationSpaceObjectPointer;
public:
    bool FilmOfStructuresActive = false;
protected:
    std::unordered_map<UniqueIdInt, Particle>* Particles = nullptr;
protected:
    std::vector<std::vector<CellEngineAtom>> ParticlesCenters;
    std::vector<std::vector<CellEngineAtom>> AllAtoms;
public:
    virtual void ReadDataFromFile(bool StartValuesBool) = 0;
protected:
    virtual UniqueIdInt AddNewParticle(const Particle& ParticleObjectParam) = 0;
    virtual void SetStartValues() = 0;
    virtual void InsertAtom(std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject, const CellEngineAtom& AppliedAtom, UniqueIdInt ParticleIndex) = 0;
    virtual void InsertGroupOfAtoms(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject, std::vector<CellEngineAtom>& LocalCellEngineAllAtomsObject) = 0;
    virtual void InsertParticlesCenters(std::vector<CellEngineAtom>& LocalCellEngineParticlesCentersObject) = 0;
    virtual void PreprocessData() = 0;
    virtual void PrintStatistics() = 0;
public:
    void SetMainDataPointers(std::unordered_map<UniqueIdInt, Particle>* ParticlesParam)
    {
        Particles = ParticlesParam;
    }
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
    [[nodiscard]] UnsignedInt GetNumberOfStructures()
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
