
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"

class CellEngineDataFile
{
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
            Center /= AtomsParam.size();
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
            if (CellEngineConfigDataObject.ChosenStructureIndex < GetNumberOfStructures() - 1)
                CellEngineConfigDataObject.ChosenStructureIndex++;
        }
        CATCH("showing next structure")
    }
public:
    void ShowPrevStructure()
    {
        try
        {
            if (CellEngineConfigDataObject.ChosenStructureIndex > 0)
                CellEngineConfigDataObject.ChosenStructureIndex--;
        }
        CATCH("showing previous structure")
    }
};

inline std::unique_ptr<CellEngineDataFile> CellEngineDataFileObjectPointer;

#endif
