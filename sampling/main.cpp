#include <mpi.h>
#include <pthread.h>
#include <stdio.h>

#include "EBStop.hpp"
#include "Graph.hpp"
#include "IEStop.hpp"

#define COMP_INSTANCES 2
#define STOR_INSTANCES 2
#define ALPHA 0.05
#define DELTA 0.05
#define GRAPH_DIR "graph"

#define TASK_TAG 0
#define SAMPLING_TAG 1
#define ESTIMATION_TAG 2

#define DEBUG

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
#ifdef DEBUG
    printf("Sampling process send graph (of size %d) to compute process no %d:",
           size, g->source);
    for (int i = 0; i < size; i++) {
        printf(" %d", subgraph[i]);
    }
    printf("\n");
#endif
    MPI_Send(subgraph, size, MPI_INT, g->source, SAMPLING_TAG, MPI_COMM_WORLD);
    // return (void *)&m;
}

int main(int argc, char **argv) {
    int num_vertex = 30;   // need init
    int num_sampling = 5;  // need init
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
        IEStop::get_instance().init(ALPHA, DELTA);
        for (; work_no <= COMP_INSTANCES; work_no++) {
            workmap[dst] = work_no;
            MPI_Send(&work_no, 1, MPI_INT, dst, TASK_TAG, MPI_COMM_WORLD);
            dst++;
        }
        while (1) {
            int result[2];
            MPI_Recv(&result, 2, MPI_INT, MPI_ANY_SOURCE, ESTIMATION_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dst = result[0];
            estimation = result[1];
            if (IEStop::get_instance().add(workmap[dst], estimation) == 0) {
                IEStop::get_instance().print_res();
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
#ifdef DEBUG
            printf("Main process send task %d to %d process.\n", workmap[dst],
                   dst);
#endif
            MPI_Send(&workmap[dst], 1, MPI_INT, dst, TASK_TAG, MPI_COMM_WORLD);
            if (flag) {
                for (int i = COMP_INSTANCES + 1;
                     i <= COMP_INSTANCES + STOR_INSTANCES; i++) {
                    MPI_Send(stopbuf, 2, MPI_INT, i, SAMPLING_TAG,
                             MPI_COMM_WORLD);
                }
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
        while (1) {
            int work_no;
            MPI_Recv(&work_no, 1, MPI_INT, 0, TASK_TAG, MPI_COMM_WORLD,
                     &status);
#ifdef DEBUG
            printf("Compute process %d received task %d.\n", my_rank, work_no);
#endif
            if (work_no == 0) {
                MPI_Finalize();
                exit(0);
            }
            memset(arr, 0, 2 * STOR_INSTANCES * sizeof(int));
            srand(rand());
            for (int i = 0; i < num_sampling; i++) {
                arr[rand() % STOR_INSTANCES][0]++;
            }
#ifdef DEBUG
            printf("Compute process %d send sampling sizes:", my_rank);
            for (int i = 0; i < STOR_INSTANCES; i++) {
                printf(" %d", arr[i][0]);
            }
            printf("\n");
#endif
            for (int i = 0; i < STOR_INSTANCES; i++) {
                arr[i][1] = my_rank;
                MPI_Isend(arr[i], 2, MPI_INT, COMP_INSTANCES + i + 1,
                          SAMPLING_TAG, MPI_COMM_WORLD, &request);
            }
            Graph graph = Graph(); // new sampling graph
            graph.init(num_vertex);
            for (int i = 0; i < STOR_INSTANCES; i++) {
                int size;
                MPI_Probe(MPI_ANY_SOURCE, SAMPLING_TAG, MPI_COMM_WORLD,
                          &status);
                MPI_Get_count(&status, MPI_INT, &size);
#ifdef DEBUG
                printf("Compute process %d probed sampling of size %d.\n",
                       my_rank, size);
#endif
                int *buf = (int *)malloc(size * sizeof(int));
                MPI_Recv(buf, size, MPI_INT, status.MPI_SOURCE, SAMPLING_TAG,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef DEBUG
                printf("Compute process %d received Zipped graph: ", my_rank);
                for (int i = 0; i < size; i++) {
                    printf("%d ", buf[i]);
                }
                printf("\n");
#endif
                if (graph.join(buf) != 0) {
                    printf("Join graph failed.\n");
                }
            }
            int result = graph.count();
            resultbuf[0] = my_rank;
            resultbuf[1] = result;
#ifdef DEBUG
            printf("Compute process %d: estimation result %d.\n", my_rank,
                   result);
#endif
            MPI_Send(resultbuf, 2, MPI_INT, 0, ESTIMATION_TAG, MPI_COMM_WORLD);
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
        bool threadinit[COMP_INSTANCES] = {false};
        graph.init_from_file(GRAPH_DIR);
        int samplingbuf[2] = {0};
        while (1) {
            MPI_Recv(samplingbuf, 2, MPI_INT, MPI_ANY_SOURCE, SAMPLING_TAG,
                     MPI_COMM_WORLD, &status);
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
                if (threadinit[source - 1]) {
                    pthread_join(threads[source - 1], NULL);
                } else {
                    threadinit[source - 1] = true;
                }
                pthread_create(&threads[source - 1], NULL, sampling, &graph);
            }
        }
    }
    return 0;
}