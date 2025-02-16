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
    RealType X{};
    RealType Y{};
    RealType Z{};
    vector3_16 AtomColor{};
public:
    #ifdef EXTENDED_RAM_MEMORY
    RealType SizeXAtom{};
    RealType SizeYAtom{};
    RealType SizeZAtom{};
    RealType SizeXParticle{};
    RealType SizeYParticle{};
    RealType SizeZParticle{};
    #endif
public:
    CellEngineAtom() = default;
    ~CellEngineAtom() = default;
public:
    CellEngineAtom(const EntityIdInt EntityIdParam, const RealType XParam, const RealType YParam, const RealType ZParam, const vector3_16 AtomColorParam) : EntityId{ EntityIdParam }, X { XParam }, Y { YParam }, Z { ZParam }, AtomColor { AtomColorParam }
    {
    }
    CellEngineAtom(const RealType XParam, const RealType YParam, const RealType ZParam) : X { XParam }, Y { YParam }, Z { ZParam }
    {
    }
public:
    void SetAtomPositionsData(RealType XParam, RealType YParam, RealType ZParam);
public:
    [[nodiscard]] vmath::vec3 Position() const
    {
        return { X, Y, Z };
    }
};

#endif