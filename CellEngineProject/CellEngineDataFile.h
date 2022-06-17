
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include <cinttypes>
#include "CellEngineAtom.h"

//using namespace sb7;

////czytac z pliku XML
//const uint64_t DNACode = 694;
//const uint64_t RNACode = 695;
//const uint64_t RIBOSOME_70SCode = 682;
//const uint64_t RIBOSOME_50SCode = 681;
//const uint64_t RIBOSOME_30SCode = 679;
//const uint64_t DNA_POLYMERASE_GAMMA_COMPLEXCode = 516;
//const uint64_t DNA_POLYMERASE_CORECode = 513;
//const uint64_t RNA_POLYMERASECode = 683;
//const uint64_t tRNACode = 688;
//const uint64_t tRNA_aminoacylatedCode = 689;

//inline vmath::vec3 FromVec4ToVec3(const vmath::vec4& Vector4)
//{
//    return vmath::vec3(Vector4.data[0], Vector4.data[1], Vector4.data[2]);
//}
/*
inline vmath::vec3 ChooseColorForParticle(const CellEngineAtom& AtomObject)
{
    vmath::vec3 ChosenColor;

    try
    {
        if (AtomObject.EntityId == RNACode)
            ChosenColor = FromVec4ToVec3(sb7::color::Blue);
        else
        if (AtomObject.EntityId == tRNA_aminoacylatedCode)
            ChosenColor = FromVec4ToVec3(sb7::color::Blue);
        else
        if (AtomObject.EntityId == tRNACode)
            ChosenColor = FromVec4ToVec3(sb7::color::Blue);
        else
        if(AtomObject.EntityId == DNACode)
            ChosenColor = FromVec4ToVec3(sb7::color::Red);
        else
        if(AtomObject.EntityId == RIBOSOME_70SCode)
            ChosenColor = FromVec4ToVec3(sb7::color::Lime);
        else
        if(AtomObject.EntityId == RIBOSOME_50SCode)
            ChosenColor = FromVec4ToVec3(sb7::color::Pink);
        else
        if(AtomObject.EntityId == RIBOSOME_30SCode)
            ChosenColor = FromVec4ToVec3(sb7::color::Orange);
        else
        if(AtomObject.EntityId == DNA_POLYMERASE_GAMMA_COMPLEXCode)
            ChosenColor = FromVec4ToVec3(sb7::color::Yellow);
        else
        if(AtomObject.EntityId == DNA_POLYMERASE_CORECode)
            ChosenColor = FromVec4ToVec3(sb7::color::Cyan);
        else
                            //if(AtomObject.EntityId == RNA_POLYMERASECode)
                            //    ChosenColor = FromVec4ToVec3(sb7::color::Purple);
                            //else
                            //if(AtomObject.EntityId == DNA_POLYMERASE_CORECode)
                            //    ChosenColor = FromVec4ToVec3(sb7::color::Cyan);
                            //else
        if(AtomObject.EntityId == RNA_POLYMERASECode)
            ChosenColor = FromVec4ToVec3(sb7::color::Purple);
        else
            ChosenColor = FromVec4ToVec3(sb7::color::Green);
    }
    CATCH("chosing color for atom for cell visualization")

    return ChosenColor;
}

inline vmath::vec3 ChooseColorForAtom(const CellEngineAtom& AtomObject)
{
    vmath::vec3 ChosenColor;

    try
    {
        switch(AtomObject.Name[0])
        {
            case 'C': ChosenColor = FromVec4ToVec3(sb7::color::Turquoise);  break;
            case 'O': ChosenColor = FromVec4ToVec3(sb7::color::Red); break;
            case 'H': ChosenColor = FromVec4ToVec3(sb7::color::White); break;
            case 'N': ChosenColor = FromVec4ToVec3(sb7::color::Blue); break;
            case 'P': ChosenColor = FromVec4ToVec3(sb7::color::Green); break;
            default: ChosenColor = FromVec4ToVec3(sb7::color::Gray); break;
        }
    }
    CATCH("chosing color for atom for cell visualization")

    return ChosenColor;
}
*/
struct ParticleKind
{
    bool Visible;
    float SizeX;
    float SizeY;
    float SizeZ;
    vmath::vec3 Color;
    std::string Name;
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
