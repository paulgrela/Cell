
#define USE_MPI

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "CellEngineImGuiMenu.h"

#include "./mds/CellEngineMolecularDynamicsSimulationForceField1.h"

#ifdef USE_MPI
int NumberOfMPIProcesses;
int MPIProcessIdentifier;
#endif

void StartMPI()
{
    MPI_Init(0, {});
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIProcessIdentifier);
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfMPIProcesses);
    cout << "MPI Process Identifier: " << NumberOfMPIProcesses << endl;
    cout << "Number of MPI processes: " << NumberOfMPIProcesses << endl;
}

int main(const int argc, const char ** argv)
{
    StartMPI();

    if(MPIProcessIdentifier == 0)
    {
        CellEngineImGuiMenu CellEngineImGuiMenuObject(argc, argv);
    }

    MPI_Finalize();

    return 0;
}