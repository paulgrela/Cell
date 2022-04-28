
#ifndef CELL_ENGINE_PROJECT_DATA_FILE_H
#define CELL_ENGINE_PROJECT_DATA_FILE_H

#include "CellEngineAtom.h"

class CellEngineDataFile
{
public:
    virtual std::vector<CellEngineAtom>& GetAtoms() = 0;
    [[nodiscard]] virtual FloatVectorType MassCenter() = 0;
};

#endif
