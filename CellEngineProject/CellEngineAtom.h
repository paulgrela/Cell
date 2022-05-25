#pragma once

#ifndef CELL_ENGINE_PROJECT_ATOM_H
#define CELL_ENGINE_PROJECT_ATOM_H

#include <utility>

#include "vmath.h"
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
    vmath::vec3 AtomColor;
    vmath::vec3 ParticleColor;
    vmath::vec3 RandomParticleColor;
public:
    CellEngineAtom(float XParam, float YParam, float ZParam, UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6], vmath::vec3 ColorParam) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam), ParticleColor(std::move(ColorParam))
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    CellEngineAtom() = default;
    ~CellEngineAtom() = default;
public:
    [[nodiscard]] vmath::vec3 Position() const
    {
        return vmath::vec3(X, Y, Z);
    }
};

#endif
