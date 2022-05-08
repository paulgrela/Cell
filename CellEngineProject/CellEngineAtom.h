#pragma once

#ifndef CELL_ENGINE_PROJECT_ATOM_H
#define CELL_ENGINE_PROJECT_ATOM_H

#include <utility>

#include "vmath.h"
#include "VectorType.h"
#include "CellEngineTypes.h"

struct CellEngineAtom
{
public:
    UnsignedIntType AtomIndex;
    IntType Serial;
    char Name[2];
    char ResName[4];
    float X;
    float Y;
    float Z;
    char Chain[6];
    UnsignedIntType EntityId;
    vmath::vec3 Color;
public:
    CellEngineAtom(float XParam, float YParam, float ZParam, UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6], vmath::vec3 ColorParam) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam), Color(std::move(ColorParam))
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    CellEngineAtom() = default;
    ~CellEngineAtom() = default;
public:
    [[nodiscard]] FloatVectorType Position() const
    {
        return FloatVectorType(X, Y, Z);
    }
};

#endif
