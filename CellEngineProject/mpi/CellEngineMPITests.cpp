
#include <mpi.h>

#include <cstring>
#include <iostream>

#include "StringUtils.h"
#include "DateTimeUtils.h"

using namespace std;

void StartMPI()
{
    int NumberOfMPIProcesses;
    int MPIProcessIdentifier;

    MPI_Init(0, {});
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIProcessIdentifier);
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfMPIProcesses);
    cout << "MPI Process Identifier: " << MPIProcessIdentifier << endl;
    cout << "Number of MPI processes: " << NumberOfMPIProcesses << endl;
}

int GenerateRandom(int size, int rank)
{
    int rdest = 0;
    if (rank < size - 1)
        rdest = rank + 1;
    else
        rdest = 0;
    return rdest;
}

inline void MPIMessagesTest()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    srand(time(NULL));
    int r = 0;
    if (rank == 0)
        r = rand() % size;

    int IndexLoop = 0;

    r = 3;

    MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Request request = MPI_REQUEST_NULL;
    if (rank == r)
        MPI_Isend(&r,1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

    if (rank == 0)
    {
        int flag = 0;

        while(!flag)
        {
            int i;
            int coffee = 42;
            for (i = 0; i < 1000000; i++)
                coffee += 1;

            MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
        }

        MPI_Recv(&r, 1, MPI_INT,MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &status);

        printf("process 0 got %d\n", r);
    }

    MPI_Wait(&request, &status);
}

void MPIMessagesTest1()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    srand(time(NULL));
    int r = 0;
    if (rank == 0)
        r = rand() % size;

    int IndexLoop = 0;

    if (rank > 0)
        while (IndexLoop < 10)
        {
            r = rank;

            MPI_Request request = MPI_REQUEST_NULL;
            if (rank == r)
                MPI_Isend(&r,1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

            IndexLoop++;
        }

    if (rank == 0)
        while (IndexLoop < 10)
        {
            int flag = 0;

            while (!flag)
            {
                int i;
                int coffee = 42;
                for (i = 0; i < 1000000; i++)
                    coffee += 1;

                MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
            }

            MPI_Recv(&r, 1, MPI_INT,MPI_ANY_SOURCE , 0, MPI_COMM_WORLD, &status);

            printf("process 0 got %d\n", r);

            IndexLoop++;
        }
}

void MPIMessagesTest2()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    srand(time(NULL));
    int r = 0;
    if (rank == 0)
        r = rand() % size;

    MPI_Request request = MPI_REQUEST_NULL;

    int IndexLoop = 0;

    if (rank > 0)
        while (IndexLoop < 10)
        {
            r = rank;

            MPI_Isend(&r,1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

            IndexLoop++;
        }

    if (rank == 0)
        while (IndexLoop < 10)
        {
            int flag = 0;

            MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

            MPI_Recv(&r, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            printf("process 0 got %d\n", r);

            IndexLoop++;
        }
}

void MPIMessagesTest3()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int rdest = 0;
    if (rank < size - 1)
        rdest = rank + 1;
    else
        rdest = 0;

    int r = rank;

    MPI_Isend(&rank, 1, MPI_INT, rdest, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, rank, rdest);

    int flag = 0;
    while (flag == 0)
    {
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag == true)
        {
            MPI_Recv(&r, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            printf("process %d got %d from process %d\n", rank, r, status.MPI_SOURCE);
        }
    }
}

void MPIMessagesTest4()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int rdest = 0;
    if (rank < size - 1)
        rdest = rank + 1;
    else
        rdest = 0;

    int r = rank;

    MPI_Bsend(&rank, 1, MPI_INT, rdest, 0, MPI_COMM_WORLD);
    printf("process %d sent %d to dest %d\n", rank, rank, rdest);

    int flag = 0;
    while (flag == 0)
    {
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag == true)
        {
            MPI_Recv(&r, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            printf("process %d got %d from process %d\n", rank, r, status.MPI_SOURCE);
        }
    }
}

void MPIMessagesTest5()
{
    int size, rank;
    MPI_Init(0, {});
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r1 = 0;
    MPI_Isend(&r1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r1, 1);

    int r2 = 1;
    MPI_Isend(&r2, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r2, 2);

    int r3 = 2;
    MPI_Isend(&r3, 1, MPI_INT, 3, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r3, 3);

    int r4 = 3;
    MPI_Isend(&r4, 1, MPI_INT, 4, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r4, 4);

    int r5 = 4;
    MPI_Isend(&r5, 1, MPI_INT, 5, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r5, 5);

    int r6 = 5;
    MPI_Isend(&r6, 1, MPI_INT, 6, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r6, 6);

    int r7 = 6;
    MPI_Isend(&r7, 1, MPI_INT, 7, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r7, 7);

    int r8 = 7;
    MPI_Isend(&r8, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    printf("process %d sent %d to dest %d\n", rank, r8, 0);

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag == true)
        {
            int r;
            MPI_Recv(&r, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            printf("process %d got %d from process %d\n", rank, r, status.MPI_SOURCE);

            Counter++;
        }
    }
}

void MPIMessagesTest6(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            MPI_Isend(&r[p], 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], p + 1);
        }
        else
        {
            MPI_Isend(&r[p], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag == true)
        {
            int rr;
            MPI_Recv(&rr, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (PrintBool == true)
                printf("process %d got %d from process %d\n", rank, rr, status.MPI_SOURCE);

            Counter++;
        }
    }
}

void MPIMessagesTest6_1(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            MPI_Isend(&r[p], 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], p + 1);
        }
        else
        {
            MPI_Isend(&r[p], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], 0);
        }
    }

    MPI_Request request1 = MPI_REQUEST_NULL;

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        if (flag == true)
        {
            int rr;

            MPI_Irecv(&rr, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request1);
            if (PrintBool == true)
                printf("process %d got %d from process %d\n", rank, rr, status.MPI_SOURCE);

            Counter++;
        }
    }
}

void MPIMessagesTest7(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            MPI_Send(&r[p], 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], p + 1);
        }
        else
        {
            MPI_Send(&r[p], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, r[p], 0);
        }
    }

    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    int rr;
    MPI_Recv(&rr, 1, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

    if (PrintBool == true)
        printf("process %d got %d from process %d\n", rank, rr, status.MPI_SOURCE);
}

void MPIMessagesTest6_2_1(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8][9] = { { 8,8,8, 8,8,8, 8,8,8 }, { 1,1,1, 1,1,1, 1,1,1 }, { 2,2,2, 2,2,2, 2,2,2 }, { 3,3,3, 3,3,3, 3,3,3 }, { 4,4,4, 4,4,4, 4,4,4 }, { 5,5,5, 5,5,5, 5,5,5 }, { 6,6,6, 6,6,6, 6,6,6 }, { 7,7,7, 7,7,7, 7,7,7 } };

    int rr[8][9];

    for (auto& r1 : rr)
        for (auto& r2 : r1)
            r2 = 0;

    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            r[p][0] = rank;
            MPI_Isend(&r[p][0], p + 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, p + 1);
        }
        else
        {
            r[0][0] = rank;
            MPI_Isend(&r[0][0], p + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        int count;
        MPI_Get_count(&status, MPI_INT, &count);

        if (flag == true)
        {
            int rrr[9];
            MPI_Recv(&rrr, count, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            int ProcessSenderOfMessage = rrr[0];
            for (int k = 0; k < count; k++)
                rr[ProcessSenderOfMessage][k] = rrr[k];

            if (PrintBool == true)
            {
                printf("process %d got %d bytes value %d from process %d\n", rank, count, rr[ProcessSenderOfMessage][0], ProcessSenderOfMessage);
                for (int k = 0; k < count; k++)
                    cout << rr[ProcessSenderOfMessage][k] << "|";
                cout << endl;
            }

            Counter++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (PrintBool == true)
    {
        cout << "Process rank = " << rank << " got message" << endl;
        for (const auto& r1 : rr)
        {
            cout << rank << " ";
            for (const auto& r2 : r1)
                cout << r2 << ",";

            cout << endl;
        }
    }
}

void MPIMessagesTest6_2_1_1(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8][9] = { { 8,8,8, 8,8,8, 8,8,8 }, { 1,1,1, 1,1,1, 1,1,1 }, { 2,2,2, 2,2,2, 2,2,2 }, { 3,3,3, 3,3,3, 3,3,3 }, { 4,4,4, 4,4,4, 4,4,4 }, { 5,5,5, 5,5,5, 5,5,5 }, { 6,6,6, 6,6,6, 6,6,6 }, { 7,7,7, 7,7,7, 7,7,7 } };

    int rr[8][9];

    for (auto& r1 : rr)
        for (auto& r2 : r1)
            r2 = 0;

    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            r[p][0] = rank;
            MPI_Isend(&r[p][0], (p + 1) * sizeof(int), MPI_CHAR, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, p + 1);
        }
        else
        {
            r[0][0] = rank;
            MPI_Isend(&r[0][0], (p + 1) * sizeof(int), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);

        if (flag == true)
        {
            int rrr[9];
            MPI_Recv(&rrr, count, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            int ProcessSenderOfMessage = rrr[0];
            for (int k = 0; k < count / 4; k++)
                rr[ProcessSenderOfMessage][k] = rrr[k];

            if (PrintBool == true)
            {
                printf("process %d got %d bytes value %d from process %d\n", rank, count, rr[ProcessSenderOfMessage][0], ProcessSenderOfMessage);
                for (int k = 0; k < count / 4; k++)
                    cout << rr[ProcessSenderOfMessage][k] << "|";
                cout << endl;
            }

            Counter++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (PrintBool == true)
    {
        cout << "Process rank = " << rank << " got message" << endl;
        for (const auto& r1 : rr)
        {
            cout << rank << " ";
            for (const auto& r2 : r1)
                cout << r2 << ",";

            cout << endl;
        }
    }
}

void MPIMessagesTest6_2_1_1_1(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    struct ALFA{ int a1; int a2; int a3; int a4; int a5; int a6; int a7; int a8; int a9; };

    ALFA r[8] = { { 8,8,8, 8,8,8, 8,8,8 }, { 1,1,1, 1,1,1, 1,1,1 }, { 2,2,2, 2,2,2, 2,2,2 }, { 3,3,3, 3,3,3, 3,3,3 }, { 4,4,4, 4,4,4, 4,4,4 }, { 5,5,5, 5,5,5, 5,5,5 }, { 6,6,6, 6,6,6, 6,6,6 }, { 7,7,7, 7,7,7, 7,7,7 } };

    ALFA rr[8];

    for (auto& r1 : rr)
        r1 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            r[p].a1 = rank;
            MPI_Isend(&r[p], (p + 1) * sizeof(int), MPI_CHAR, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, p + 1);
        }
        else
        {
            r[0].a1 = rank;
            MPI_Isend(&r[0], (p + 1) * sizeof(int), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);

        if (flag == true)
        {
            ALFA rrr;
            MPI_Recv(&rrr, count, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            int ProcessSenderOfMessage = rrr.a1;

            memcpy(&rr[ProcessSenderOfMessage], &rrr, count);

            if (PrintBool == true)
            {
                printf("process %d got %d bytes value %d from process %d\n", rank, count, rr[ProcessSenderOfMessage].a1, ProcessSenderOfMessage);

                cout << rr[ProcessSenderOfMessage].a1 << "|";
                cout << rr[ProcessSenderOfMessage].a2 << "|";
                cout << rr[ProcessSenderOfMessage].a3 << "|";
                cout << rr[ProcessSenderOfMessage].a4 << "|";
                cout << rr[ProcessSenderOfMessage].a5 << "|";
                cout << rr[ProcessSenderOfMessage].a6 << "|";
                cout << rr[ProcessSenderOfMessage].a7 << "|";
                cout << rr[ProcessSenderOfMessage].a8 << "|";
                cout << rr[ProcessSenderOfMessage].a9 << "|";
                cout << endl;
            }

            Counter++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (PrintBool == true)
    {
        cout << "Process rank = " << rank << " got message" << endl;
        for (const auto& r1 : rr)
        {
            cout << rank << " ";

            cout << r1.a1 << ",";
            cout << r1.a2 << ",";
            cout << r1.a3 << ",";
            cout << r1.a4 << ",";
            cout << r1.a5 << ",";
            cout << r1.a6 << ",";
            cout << r1.a7 << ",";
            cout << r1.a8 << ",";
            cout << r1.a9 << ",";

            cout << endl;
        }
    }
}

void MPIMessagesTest6_2_2(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8][9] = { { 8,8,8, 8,8,8, 8,8,8 }, { 1,1,1, 1,1,1, 1,1,1 }, { 2,2,2, 2,2,2, 2,2,2 }, { 3,3,3, 3,3,3, 3,3,3 }, { 4,4,4, 4,4,4, 4,4,4 }, { 5,5,5, 5,5,5, 5,5,5 }, { 6,6,6, 6,6,6, 6,6,6 }, { 7,7,7, 7,7,7, 7,7,7 } };

    int rr[8][9];

    for (auto& r1 : rr)
        for (auto& r2 : r1)
            r2 = 0;

    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            r[p][0] = rank;
            MPI_Isend(&r[p][0], p + 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, p + 1);
        }
        else
        {
            r[0][0] = rank;
            MPI_Isend(&r[0][0], p + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);

        int count;
        MPI_Get_count(&status, MPI_INT, &count);

        if (flag == true)
        {
            int rrr[9];
            MPI_Irecv(&rrr, count, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request);

            int ProcessSenderOfMessage = rrr[0];
            for (int k = 0; k < count; k++)
                rr[ProcessSenderOfMessage][k] = rrr[k];

            if (PrintBool == true)
            {
                printf("process %d got %d bytes value %d from process %d\n", rank, count, rr[ProcessSenderOfMessage][0], ProcessSenderOfMessage);
                for (int k = 0; k < count; k++)
                    cout << rr[ProcessSenderOfMessage][k] << "|";
                cout << endl;
            }

            Counter++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (PrintBool == true)
    {
        cout << "Process rank = " << rank << " got message" << endl;
        for (const auto& r1 : rr)
        {
            cout << rank << " ";
            for (const auto& r2 : r1)
                cout << r2 << ",";

            cout << endl;
        }
    }
}

void MPIMessagesTest6_3(const bool PrintBool)
{
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Request request = MPI_REQUEST_NULL;

    int r[8][9] = { { 8,8,8, 8,8,8, 8,8,8 }, { 1,1,1, 1,1,1, 1,1,1 }, { 2,2,2, 2,2,2, 2,2,2 }, { 3,3,3, 3,3,3, 3,3,3 }, { 4,4,4, 4,4,4, 4,4,4 }, { 5,5,5, 5,5,5, 5,5,5 }, { 6,6,6, 6,6,6, 6,6,6 }, { 7,7,7, 7,7,7, 7,7,7 } };

    int rr[8][9];
    for (auto& r1 : rr)
        for (auto& r2 : r1)
            r2 = 0;

    for (int p = 0; p < 8; p++)
    {
        if (p < 7)
        {
            r[p][0] = rank;
            MPI_Send(&r[p][0], p + 1, MPI_INT, p + 1, 0, MPI_COMM_WORLD);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, p + 1);
        }
        else
        {
            r[0][0] = rank;
            MPI_Send(&r[0][0], p + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            if (PrintBool == true)
                printf("process %d sent %d to dest %d\n", rank, p + 1, 0);
        }
    }

    int Counter = 0;
    while (Counter < 8)
    {
        MPI_Status status;

        MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,  &status);

        int count;
        MPI_Get_count(&status, MPI_INT, &count);

        int rrr[9];
        MPI_Recv(&rrr, count, MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        int ProcessSenderOfMessage = rrr[0];
        for (int k = 0; k < count; k++)
            rr[ProcessSenderOfMessage][k] = rrr[k];

        if (PrintBool == true)
        {
            printf("process %d got %d bytes value %d from process %d\n", rank, count, rr[ProcessSenderOfMessage][0], ProcessSenderOfMessage);
            for (int k = 0; k < count; k++)
                cout << rr[ProcessSenderOfMessage][k] << "|";
            cout << endl;
        }

        Counter++;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (PrintBool == true)
    {
        cout << "Process rank = " << rank << " got message" << endl;
        for (const auto& r1 : rr)
        {
            cout << rank << " ";
            for (const auto& r2 : r1)
                cout << r2 << ",";

            cout << endl;
        }
    }
}

void TestUnblockingMessagesTime()
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest6(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time1 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of generating mpi messages unblocking has taken time: ","Execution in threads") << endl;
}

void TestUnblockingMessagesTime_1()
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest6_1(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time1 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of generating mpi messages unblocking full has taken time: ","Execution in threads") << endl;
}

void TestBlockingMessagesTime()
{
    const auto start_time2 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest7(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time2 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time2, stop_time2, "Execution of generating mpi messages blocking has taken time: ","Execution in threads") << endl;
}

void TestUnblockingMessagesLongTime()
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest6_2_1(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time1 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of generating mpi messages unblocking full has taken time: ","Execution in threads") << endl;
}

void TestUnblockingMessagesLongTime_1()
{
    const auto start_time1 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest6_2_2(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time1 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time1, stop_time1, "Execution of generating mpi messages unblocking full full has taken time: ","Execution in threads") << endl;
}

void TestBlockingMessagesLongTime()
{
    const auto start_time2 = chrono::high_resolution_clock::now();

    for (int l = 0; l < 1000; l++)
    {
        MPIMessagesTest6_3(false);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    const auto stop_time2 = chrono::high_resolution_clock::now();

    cout << GetDurationTimeInOneLineStr(start_time2, stop_time2, "Execution of generating mpi messages blocking has taken time: ","Execution in threads") << endl;
}

void TestExecutionsShortMessagesTimes()
{
    TestUnblockingMessagesTime();
    TestUnblockingMessagesTime_1();
    TestBlockingMessagesTime();
}

void TestExecutionsLongMessagesTimes()
{
    TestUnblockingMessagesLongTime();
    TestUnblockingMessagesLongTime_1();
    TestBlockingMessagesLongTime();
}

void MPIAllTests()
{
    StartMPI();

    MPIMessagesTest6(true);

    MPIMessagesTest6_1(true);

    MPIMessagesTest7(true);

    MPIMessagesTest6_2_1(true);

    MPIMessagesTest6_2_1_1(true);

    MPIMessagesTest6_2_1_1_1(true);

    MPIMessagesTest6_2_2(true);

    MPIMessagesTest6_3(true);

    TestExecutionsLongMessagesTimes();

    MPI_Finalize();
}
