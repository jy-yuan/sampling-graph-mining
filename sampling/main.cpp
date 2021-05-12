#include <mpi.h>
#include <pthread.h>
#include <stdio.h>

#include "EBStop.hpp"
#include "Graph.hpp"

#define COMP_INSTANCES 4
#define STOR_INSTANCES 4

#define TASK_TAG 0
#define SAMPLING_TAG 1
#define ESTIMATION_TAG 2

/*
pthread
sampling and send (non-blocking) zipped graph to computation process
*/
void *sampling(void *Param) {
    Graph *g = (Graph *)Param;
    int *subgraph = g->sample(g->m);
    int m = subgraph[0];
    int n = subgraph[1];
    int size = 2 * m + n + 3;
    MPI_Send(subgraph, sizeof(subgraph), MPI_INT, g->source, SAMPLING_TAG, MPI_COMM_WORLD);
    // return (void *)&m;
}

int main(int argc, char **argv) {
    int num_vertex; // need init
    int num_sampling; // need init
    int my_rank;
    int provided;
    bool single_thread = false;
    MPI_Request request;
    MPI_Status status;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("MPI do not Support Multiple thread\n");
        single_thread = true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0) {
        /*
        Main process,
        send (sampling) Work No. to each Computation process
        receive results from every work and determine whether stop
        */
        bool ended = false;
        int work_no = 1;
        int workmap[COMP_INSTANCES + 1];
        int dst = 1;
        int estimation;
        int stopbuf[2] = {0};
        EBStop::get_instance().init(0.05, 0.05);
        for (; work_no <= COMP_INSTANCES; work_no++) {
            workmap[dst++] = work_no;
            MPI_Send(&work_no, 1, MPI_INT, dst, TASK_TAG, MPI_COMM_WORLD);
        }
        while (1) {
            int result[2];
            MPI_Recv(&result, 2, MPI_INT, MPI_ANY_SOURCE, ESTIMATION_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dst = result[0];
            estimation = result[1];
            if (EBStop::get_instance().add(workmap[dst], estimation) == 0) {
                EBStop::get_instance().print_res();
                for (int i = COMP_INSTANCES + 1;
                     i <= COMP_INSTANCES + STOR_INSTANCES; i++) {
                    MPI_Isend(stopbuf, 2, MPI_INT, i, SAMPLING_TAG, MPI_COMM_WORLD, &request);
                }
                ended = true;
            }
            if (ended) {
                workmap[dst] = 0;
            } else {
                workmap[dst] = work_no++;
            }
            bool flag = true;  // all exit
            for (int i = 1; i <= COMP_INSTANCES; i++) {
                if (workmap[i] != 0) {
                    flag = false;
                }
            }
            MPI_Send(&workmap[dst], 1, MPI_INT, dst, TASK_TAG, MPI_COMM_WORLD);
            if (flag) {
                MPI_Finalize();
                exit(0);
            }
        }
    } else if (my_rank <= COMP_INSTANCES) {
        /*
        computation process
        receive work no. from main process
        send sample size to each storage process
        receive sample, combine them and count the estimated result
        */
        int arr[STOR_INSTANCES][2] = {0};  // random sizes
        int resultbuf[2] = {0};
        Graph graph = Graph();
        
        while (1) {
            int work_no;
            MPI_Recv(&work_no, 1, MPI_INT, 0, TASK_TAG, MPI_COMM_WORLD, &status);
            if (work_no == 0) {
                MPI_Finalize();
                exit(0);
            }
            memset(arr, 0, 2 * STOR_INSTANCES * sizeof(int));
            srand(time(0));
            for (int i = 0; i < num_sampling; i++) {
                arr[rand() % STOR_INSTANCES][0]++;
            }
            for (int i = 0; i < STOR_INSTANCES; i++) {
                arr[i][1] = my_rank;
                MPI_Isend(arr[i], 2, MPI_INT, COMP_INSTANCES + i + 1, SAMPLING_TAG, MPI_COMM_WORLD, &request);
            }
            for (int i = 0; i < STOR_INSTANCES; i++) {
                int size;
                MPI_Probe(MPI_ANY_SOURCE, SAMPLING_TAG, MPI_COMM_WORLD,
                          &status);
                MPI_Get_count(&status, MPI_INT, &size);
                int *buf = (int *)malloc(size * sizeof(int));
                MPI_Recv(buf, size, MPI_INT, MPI_ANY_SOURCE, SAMPLING_TAG,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // maybe problem
                graph.join(buf);
            }
            int result = graph.count();
            resultbuf[0] = my_rank;
            resultbuf[1] = result;
            MPI_Send(resultbuf, 2,MPI_INT,0, ESTIMATION_TAG, MPI_COMM_WORLD);
        }
    } else {
        /*
        storage process
        read from file of splited graph
        receive sample size and do sampling, then send back
        in different thread
        */
        Graph graph = Graph();
        pthread_t threads[COMP_INSTANCES];
        graph.init();
        int samplingbuf[2] = {0};
        while (1) {
            MPI_Recv(samplingbuf, 2, MPI_INT,MPI_ANY_SOURCE,SAMPLING_TAG, MPI_COMM_WORLD, &status);
            int m = samplingbuf[0];
            int source = samplingbuf[1];
            if (source == 0) {
                MPI_Finalize();
                exit(0);
            }
            graph.source = source;
            graph.m = m;
            if (single_thread) {
                sampling((void *)&graph);
            } else {
                pthread_create(&threads[source - 1], NULL, sampling, &graph);
            }
        }
    }
    return 0;
}