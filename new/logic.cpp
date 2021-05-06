#include <mpi.h>
#include <pthread.h>
#include <stdio.h>
#include "Graph.hpp"
#include "EBStop.hpp"

#define COMP_INSTANCES 4
#define STOR_INSTANCES 4

int main(int argc, char **argv)
{
    int num_vertex;
    int num_sampling;
    int my_rank;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE)
    {
        printf("MPI do not Support Multiple thread\n");
        exit(0);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
    {
    /*
    Main process,
    send (sampling) Work No. to each Computation process
    receive results from every work and determine whether stop
    */
        int work_no = 1;
        int workmap[COMP_INSTANCES + 1];
        int dst = 1;
        int estimation;
        EBStop::get_instance().init();
        for (; work_no <= COMP_INSTANCES; work_no++)
        {
            workmap[dst++] = work_no;
            MPI_Send(dst);
        }
        while (1)
        {
            MPI_Recv(result, source, MPI_ANY_SOURCE);
            dst = result[0];
            estimation = result[1];
            if (EBStop::get_instance().add(workmap[dst], estimation) == 0)
            {
                EBStop::get_instance().print_res();
                MPI_Bcast(work_no = 0);
                MPI_Finalize();
                exit(0);
            }
            workmap[dst] = work_no++;
            MPI_Send(dst);
        }
    }
    else if (my_rank <= COMP_INSTANCES)
    {
    /*
    computation process
    receive work no. from main process
    send sample size to each storage process
    receive sample, combine them and count the estimated result
    */
        int arr[STOR_INSTANCES] = {0};
        graph = Graph();
        while (1)
        {
            MPI_Recv(work_no, source = 0);
            if (work_no == 0)
            {
                MPI_Finalize();
                exit(0);
            }
            memset(arr, 0, STOR_INSTANCES * sizeof(int));
            srand(time(0));
            for (int i = 0; i < num_sampling; i++)
            {
                arr[rand() % STOR_INSTANCES]++;
            }
            for (int i = 0; i < STOR_INSTANCES; i++)
            {
                MPI_Send(arr[i], my_rank, COMP_INSTANCES + i + 1);
            }
            for (int i = 0; i < STOR_INSTANCES; i++)
            {
                MPI_Status status;
                int size;
                MPI_Probe(MPI_ANY_SOURCE, tag, status);
                MPI_Get_count(status, size);
                buf[i] = malloc(size);
                MPI_Recv(buf[i], size, status.source);
            }
            result = count(buf);
            MPI_Send(0, result, my_rank);
        }
    }
    else
    {
    /*
    storage process
    read from file of splited graph
    receive sample size and do sampling, then send back
    in different thread
    */
        graph = Graph();
        graph.init();
        while (1)
        {
            MPI_Recv(m, source, MPI_ANY_SOURCE);
            if (source == 0)
            {
                MPI_Finalize();
                exit(0);
            }
            subgraph = thread(m, source);
        }
    }
    return 0;
}

void *sampling(void *Param)
{
    MPI_Isend(subgraph, source);
    return;
}