
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include <cinttypes>
#include "CellEngineAtom.h"

struct ParticleKind
{
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vmath::vec3 Color;
    std::string NameFromXML;
    std::string NameFromDataFile;
};

struct AtomKind
{
    float SizeX;
    float SizeY;
    float SizeZ;
    vmath::vec3 Color;
};

class CellEngineDataFile
{
public:
    std::string CellStateFileName;
public:
    float SpecularPower = 8.0f;
    float SpecularAlbedo = 0.2222f;
public:
    enum class MakeColorsType
    {
        DrawColorForEveryAtom = 1,
        DrawColorForEveryParticle = 2,
        DrawRandomColorForEveryParticle = 3
    };
    MakeColorsType MakeColorsTypeObject;
public:
    enum class StencilForDrawingObjectsTypes
    {
        StencilForDrawingOnlyParticlesCenters = 1,
        StencilForDrawingOnlyInAtomScale = 2
    };
    StencilForDrawingObjectsTypes StencilForDrawingObjectsTypesObject;
    std::uint64_t NumberOfStencilBufferLoops;
public:
    bool DrawBondsBetweenParticlesCenters;
    bool DrawBondsBetweenAtoms;
public:
    bool ShowDetailsInAtomScale = false;
    bool CheckAtomVisibility;
    float CutZ;
    float Distance;
    std::uint64_t LoadOfAtomsStep;
public:
    float XLowToDrawInAtomScale;
    float XHighToDrawInAtomScale;
    float YLowToDrawInAtomScale;
    float YHighToDrawInAtomScale;
    float ZLowToDrawInAtomScale;
    float ZHighToDrawInAtomScale;
public:
    enum class SizeOfAtomsDrawingTypes
    {
        AtomSize = 1,
        ParticleSize = 2,
        AutomaticChangeSize = 3
    };
    SizeOfAtomsDrawingTypes SizeOfAtomsDrawingTypesObject;
public:
    float SizeStep;
    float SizeOfAtomX;
    float SizeOfAtomY;
    float SizeOfAtomZ;
    float CameraXPosition;
    float CameraYPosition;
    float CameraZPosition;
    float CameraXMoveShortStep;
    float CameraYMoveShortStep;
    float CameraZMoveShortStep;
    float CameraXMoveLongStep;
    float CameraYMoveLongStep;
    float CameraZMoveLongStep;
    float ViewXMoveShortStep;
    float ViewYMoveShortStep;
    float ViewZMoveShortStep;
    float ViewXMoveLongStep;
    float ViewYMoveLongStep;
    float ViewZMoveLongStep;
public:
    bool ViewChangeUsingLongStep;
    bool AutomaticChangeOfSizeOfAtom;
    bool AutomaticChangeOfLoadAtomsStep;
public:
    UnsignedIntType ChosenStructureIndex = 0;
    bool FilmOfStructuresActive = false;
public:
    vmath::vec3 BackgroundColors[4];
    std::uint64_t ChosenBackgroundColor;
public:
    std::unordered_map<std::string, AtomKind> AtomsKinds;
    std::unordered_map<UnsignedIntType, ParticleKind> ParticlesKinds;
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
