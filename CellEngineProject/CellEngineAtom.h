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
    EntityIdInt EntityId{};
    char Name[5 + 1]{};
    char ResName[4 + 1]{};
    char Chain[6 + 1]{};
    float X{};
    float Y{};
    float Z{};
    vector3_16 AtomColor{};
public:
    #ifdef EXTENDED_RAM_MEMORY
    float SizeXAtom{};
    float SizeYAtom{};
    float SizeZAtom{};
    float SizeXParticle{};
    float SizeYParticle{};
    float SizeZParticle{};
    #endif
public:
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