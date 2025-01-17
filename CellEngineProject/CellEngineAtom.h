#pragma once

#ifndef CELL_ENGINE_PROJECT_ATOM_H
#define CELL_ENGINE_PROJECT_ATOM_H

#include <utility>
#include <cstring>
#include <string>

#include "vmath.h"
#include "CellEngineTypes.h"
#include "CellEngineConstants.h"

class CellEngineAtom
{
public:
    UnsignedInt AtomIndex{};
    UnsignedInt Serial{};
    EntityIdInt EntityId{};
    char Name[5 + 1]{};
    char ResName[4 + 1]{};
    char Chain[6 + 1]{};
    float X{};
    float Y{};
    float Z{};
    // std::uint16_t X{};
    // std::uint16_t Y{};
    // std::uint16_t Z{};
    // float XR{};
    // float YR{};
    // float ZR{};
    vector3_16 AtomColor{};
    vector3_16 ParticleColor{};
    vector3_16 UniqueParticleColor{};
    vector3_16 RandomParticleKindColor{};
    bool Visible{};
    UniqueIdInt GenomeIndex = 0;
    UniqueIdInt GenomeIndexPrev = 0;
    UniqueIdInt GenomeIndexNext = 0;
    char Nucleotide = '?';
    #ifdef RNA_IN_ONE_PARTICLE
    std::string SequenceStr;
    #endif
    #ifdef EXTENDED_RAM_MEMORY
    float SizeXAtom{};
    float SizeYAtom{};
    float SizeZAtom{};
    float SizeXParticle{};
    float SizeYParticle{};
    float SizeZParticle{};
    #endif
public:
    CellEngineAtom(const float XParam, const float YParam, const float ZParam, const UnsignedInt AtomIndexParam, const UnsignedInt SerialParam, const UnsignedInt EntityIdParam, char NameParam[2], char ResNameParam[4], char ChainParam[6], const vector3_16 ColorParam) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam), EntityId(EntityIdParam), ParticleColor(ColorParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    CellEngineAtom(const std::uint16_t XParam, const std::uint16_t YParam, const std::uint16_t ZParam) : X(XParam), Y(YParam), Z(ZParam)
    {
    }
    CellEngineAtom() = default;
    ~CellEngineAtom() = default;
public:
    void SetAtomPositionsData(float XParam, float YParam, float ZParam);
public:
    [[nodiscard]] vmath::vec3 Position() const
    {
        //return { XR, YR, ZR };
        return { X, Y, Z };
    }
};

#endif