#pragma once

#ifndef CELL_ENGINE_PROJECT_ATOM_H
#define CELL_ENGINE_PROJECT_ATOM_H

struct AtomBase
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
public:
    AtomBase(float XParam, float YParam, float ZParam,  UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6]) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    AtomBase() = default;
    ~AtomBase() = default;
public:
    [[nodiscard]] FloatVectorType Position() const
    {
        return FloatVectorType(X, Y, Z);
    }
};

#endif
