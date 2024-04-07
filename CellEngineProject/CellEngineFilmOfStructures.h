
#ifndef CELL_ENGINE_FILM_OF_STRUCTURES_H
#define CELL_ENGINE_FILM_OF_STRUCTURES_H

#include <memory>

#include "CellEngineAtom.h"
#include "CellEngineConfigData.h"
#include "CellEngineVoxelSimulationSpace.h"

class CellEngineFilmOfStructures
{
public:
    bool FilmOfStructuresActive = false;
public:
    virtual UnsignedInt GetNumberOfStructures() = 0;
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

#endif
