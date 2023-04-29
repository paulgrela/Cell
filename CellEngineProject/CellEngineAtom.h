#pragma once

#ifndef CELL_ENGINE_PROJECT_ATOM_H
#define CELL_ENGINE_PROJECT_ATOM_H

#include <utility>
#include <cstring>

#include "vmath.h"
#include "CellEngineTypes.h"

#define EXTENDED_RAM_MEMORY_

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
    vector3_16 AtomColor{};
    vector3_16 ParticleColor{};
    vector3_16 UniqueParticleColor{};
    vector3_16 RandomParticleKindColor{};
    bool Visible{};
    #ifdef EXTENDED_RAM_MEMORY
    float SizeXAtom{};
    float SizeYAtom{};
    float SizeZAtom{};
    float SizeXParticle{};
    float SizeYParticle{};
    float SizeZParticle{};
    #endif
public:
    CellEngineAtom(float XParam, float YParam, float ZParam, UnsignedInt AtomIndexParam, UnsignedInt SerialParam, UnsignedInt EntityIdParam, char NameParam[2], char ResNameParam[4], char ChainParam[6], vector3_16 ColorParam) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam), EntityId(EntityIdParam), ParticleColor(ColorParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    CellEngineAtom() = default;
    ~CellEngineAtom() = default;
public:
    void SetAtomPositionsData(float XParam, float YParam, float ZParam);
public:
    [[nodiscard]] vmath::vec3 Position() const
    {
        return { X, Y, Z };
    }
};

#endif