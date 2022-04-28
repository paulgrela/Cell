#pragma once

#ifndef CELL_ENGINE_CIF_DATA_FILE_H
#define CELL_ENGINE_CIF_DATA_FILE_H

#include <string>
#include <vector>

#include <fstream>
#include <stdexcept>

#include <memory>
#include <unordered_map>

#include "VectorType.h"

#include "CellEngineTypes.h"
#include "CellEngineAtom.h"

/*
struct Atom
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
    Atom(const char* CIFRecord, UnsignedIntType AtomIndex);
    Atom(float XParam, float YParam, float ZParam,  UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6]) : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam)
    {
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    ~Atom() = default;
public:
    [[nodiscard]] FloatVectorType Position() const;
private:
    void ParseRecord(const char* LocalPDBRecord);
};
*/
/*
class AtomCIF: public AtomBase
{
public:
    AtomCIF() = default;
    AtomCIF(float XParam, float YParam, float ZParam,  UnsignedIntType AtomIndexParam, IntType SerialParam, char NameParam[2], char ResNameParam[4], char ChainParam[6])// : X(XParam), Y(YParam), Z(ZParam), AtomIndex(AtomIndexParam), Serial(SerialParam)
    {
        X = XParam;
        Y = YParam;
        Z = ZParam;
        AtomIndex = AtomIndexParam;
        Serial = SerialParam;
        strncpy(Name, NameParam, 2);
        strncpy(ResName, ResNameParam, 4);
        strncpy(Chain, ChainParam, 6);
    }
    AtomCIF(const char* CIFRecord, UnsignedIntType AtomIndex);
public:
    //[[nodiscard]] FloatVectorType Position() const;
    //void ParseRecord(const char* LocalPDBRecord);
};
*/
struct Matrix
{
public:
    UnsignedIntType MatrixId;
    float Matrix[3][4];
};

class CIFDataFile
{
private:
    std::vector<Matrix> Matrixes;
    std::unordered_map<std::string, std::vector<AtomBase>> ChainsNames;
private:
    std::vector<std::vector<AtomBase>> Atoms;
private:
    void ReadDataFromFile(const std::string_view LocalCIFRecord);
public:
    UnsignedIntType ChosenStructureIndex;
public:
    explicit CIFDataFile(const std::string_view FileName);
    ~CIFDataFile() = default;
public:
    static AtomBase ParseRecord(const char* LocalPDBRecord);
    [[nodiscard]] std::vector<AtomBase>& GetAtoms()
    {
        return Atoms[ChosenStructureIndex];
    }
    [[nodiscard]] IntType GetNumberOfStructures() const
    {
        return Atoms.size();
    }
    [[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const AtomBase& GetAtom(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() const;
};

#endif