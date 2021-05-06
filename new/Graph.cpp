#include "Graph.hpp"

/*
read a partial graph from doc like:
M N
u1 v1
u2 v2
...
u_n v_n

where M is # of vertexes and N is # of edges
*/
int Graph::init(const std::string dir)
{
    FILE *pFile = fopen(dir.c_str(), "r");
    fscanf(pFile, "%d %d", &M, &N);
    verExi.resize(M);
    csrInd.resize(M + 1);
    csrList.resize(N);
    int u, v;
    for (int i = 0; i < N; i++) {
        fscanf(pFile, "%d %d", &u, &v);
        verExi[u] = 1;
        csrList[csrInd[u]] = v;
        for (int j = u + 1; j <= M; j++) {
            csrInd[j]++;
        }
    }
    assert(csrInd[M] == N);
    return 0;
}

/*
join another partial graph into current graph
zipgraph is a data structure of graph (easy to trans by MPI)
contains,
int: num of vertex
int: num of edge
array: whether the vertex exits
array: csr index
array: csr list
*/
int Graph::join(int *zipgraph)
{
    int m = zipgraph[0];
    if (m != M) {
        printf("ERROR: Vertex of two graph are not the same!\n");
        return -1;
    }
    int n = zipgraph[1];
    int *newExi = zipgraph + 2;
    int *newInd = zipgraph + 2 + m;
    int *newList = zipgraph + 2 + m + m + 1;
    for (int i = 0; i < m; i++) {
        verExi[i] = (int)(newExi[i] && verExi[i]);
    }
    for (int i = 0; i < m; i++) {
        for (int j = newInd[i]; j < newInd[i + 1]; j++) {
            csrList.insert(csrList.begin() + csrInd[i + 1]++, newList[j]);
        }
    }
    assert(csrList.size() == N + n);
    assert(csrInd[m] == csrList.size());
    N += n;
    return 0;
}

/*
sample [num] vertexes from graph
in the structure of zipgraph
*/
int *Graph::sample(int num)
{
    int a[];
    return a;
}

//naive implementation of sprecise count
int Graph::count()
{
    int count = 0;
    for (int i = 0; i < M; i++) {
        if (verExi[i])
            count++;
    }
    return count;
}