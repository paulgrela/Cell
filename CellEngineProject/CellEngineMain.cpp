
#include "mpi/CellEngineMPITests.h"

#include "../Compilation/ConditionalCompilationConstants.h"

#include "CellEngineImGuiMenu.h"

#include "./mds/CellEngineMolecularDynamicsSimulationForceField2.h"

int main(const int argc, const char** argv)
{
    CellEngineImGuiMenu CellEngineImGuiMenuObject(argc, argv);

    #ifdef COMPUTE_MOLLECULAR_DYNAMICS
    ComputeMolecularDynamicsSimulationForceField2();
    #endif

    return 0;
}
