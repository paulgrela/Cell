
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include <cinttypes>
#include "CellEngineAtom.h"

class CellEngineDataFile
{
public:
    float SizeX;
    float SizeY;
    float SizeZ;
    float SizeStep;
    uint64_t LoadOfAtomsStep;
    float ViewStep;
    float CameraXMoveStep;
    float CameraYMoveStep;
    float CameraZMoveStep;
    float ViewXMoveStep;
    float ViewYMoveStep;
    float ViewZMoveStep;
    bool ShowDetailsInAtomScale = false;
public:
    UnsignedIntType ChosenStructureIndex = 0;
    bool FilmOfStructuresActive = false;
protected:
    std::vector<std::vector<CellEngineAtom>> ParticlesCenters;
    std::vector<std::vector<CellEngineAtom>> AllAtoms;
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
        return ParticlesCenters[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures()
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
            ChosenStructureIndex = 0;
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
            if (ChosenStructureIndex < GetNumberOfStructures() - 1)
                ChosenStructureIndex++;
        }
        CATCH("showing next structure")
    }
public:
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
