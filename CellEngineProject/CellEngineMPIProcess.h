
#ifndef CELL_ENGINE_MPI_PROCESS_H
#define CELL_ENGINE_MPI_PROCESS_H

#include "CellEngineTypes.h"

class MPIProcessData
{
public:
    UnsignedInt CurrentMPIProcessIndex{ 0 };
    UnsignedInt NumberOfMPIProcesses{ 0 };
};

inline MPIProcessData MPIProcessDataObject;

#endif
