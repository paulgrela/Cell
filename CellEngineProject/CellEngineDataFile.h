
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include "CellEngineAtom.h"

class CellEngineDataFile
{
public:
    UnsignedIntType ChosenStructureIndex = 0;
    bool FilmOfStructuresActive = false;
public:
    virtual std::vector<CellEngineAtom>& GetAtoms() = 0;
    [[nodiscard]] virtual FloatVectorType MassCenter() = 0;
public:
    [[nodiscard]] virtual IntType GetNumberOfStructures() = 0;
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

    void StartFilmOfStructures()
    {
        try
        {
            ChosenStructureIndex = 0;
            FilmOfStructuresActive = true;
        }
        CATCH("starting film of structures")
    }

    void StopFilmOfStructures()
    {
        try
        {
            FilmOfStructuresActive = false;
        }
        CATCH("starting film of structures")
    }

    void ShowNextStructure()
    {
        try
        {
            if (ChosenStructureIndex < GetNumberOfStructures() - 1)
                ChosenStructureIndex++;
        }
        CATCH("showing next structure")
    }

    void ShowPrevStructure()
    {
        try
        {
            if (ChosenStructureIndex > 0)
                ChosenStructureIndex--;
        }
        CATCH("showing previous structure")
    }
};

#endif
