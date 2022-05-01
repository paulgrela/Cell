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

#include "CellEngineDataFile.h"

struct AtomsPositionMatrix3x4
{
public:
    UnsignedIntType MatrixId;
    float Matrix[3][4];
};

class CellEngineCIFDataFile: public CellEngineDataFile
{
private:
    std::vector<AtomsPositionMatrix3x4> AtomsPositiosnMatrixes;
    std::unordered_map<std::string, std::vector<CellEngineAtom>> ChainsNames;
//private:
//    std::vector<std::vector<CellEngineAtom>> Atoms;
//    std::vector<std::vector<CellEngineAtom>> AllAtoms;
private:
    void ReadDataFromFile(std::string_view LocalCIFRecord);
public:
    explicit CellEngineCIFDataFile(std::string_view FileName);
    ~CellEngineCIFDataFile() = default;
public:
    static CellEngineAtom ParseRecord(const char* LocalPDBRecord);
//    [[nodiscard]] std::vector<CellEngineAtom>& GetAtoms() override
//    {
//        return Atoms[ChosenStructureIndex];
//    }
//    [[nodiscard]] IntType GetNumberOfStructures() override
//    {
//        return Atoms.size();
//    }
    [[nodiscard]] IntType GetNumberOfAtoms() const;
    [[nodiscard]] const CellEngineAtom& GetAtom(IntType Index) const;
    [[nodiscard]] FloatVectorType MassCenter() override;
};

#endif