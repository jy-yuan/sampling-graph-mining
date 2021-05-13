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
int Graph::init(const std::string dir) {
    FILE *pFile = fopen(dir.c_str(), "r");
    fscanf(pFile, "%d %d", &M, &N);
    verExi.resize(M);
    verDeg.resize(M);
    csrInd.resize(M + 1);
    csrList.reserve(N);
    int u, v;
    for (int i = 0; i < N; i++) {
        fscanf(pFile, "%d %d", &u, &v);
        verExi[u] = 1;
        verDeg[u]++;
        csrList.insert(csrList.begin() + csrInd[u], v);
        for (int j = u + 1; j <= M; j++) {
            csrInd[j]++;
        }
    }
    for (int i = 0; i < M; i++) {
        if (verExi[i]) {
            vertexes.push_back(i);
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
array: whether the vertex exits TODO:optimize
array: csr index
array: csr list
*/
int Graph::join(int *zipgraph) {
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
        csrInd[i + 1] = csrInd[i] + verDeg[i];
        for (int j = newInd[i]; j < newInd[i + 1]; j++) {
            csrList.insert(csrList.begin() + csrInd[i + 1]++, newList[j]);
        }
        verDeg[i] = verDeg[i] + newInd[i + 1] - newInd[i];
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
int *Graph::sample(int num) {
    std::vector<int> newExi;
    std::vector<int> newInd;
    std::vector<int> newList;
    newExi.resize(M);
    newInd.resize(M + 1);
    for (int i = 0; i < num; i++) {
        srand(time(0));
        int a = rand() % vertexes.size();
        while (newExi[vertexes[a]]) {
            a = rand() % vertexes.size();
        }
        newExi[vertexes[a]] = 1;
    }
    for (int i = 0; i < M; i++) {
        newInd[i + 1] = newInd[i];
        if (newExi[i]) {
            for (int j = csrInd[i]; j < csrInd[i + 1]; j++) {
                newList.insert(newList.begin() + newInd[i + 1]++, csrInd[j]);
            }
        }
    }
    int *sam = (int *)malloc((2 + M + M + 1 + newList.size()) * sizeof(int));
    sam[0] = M;
    sam[1] = newList.size();
    std::copy(newExi.begin(), newExi.end(), &sam[2]);
    std::copy(newInd.begin(), newInd.end(), &sam[2 + M]);
    memcpy(&sam[2 + M + M + 1], newList.data(), sizeof(int) * newList.size());
    return sam;
}

// naive implementation of precise count
int Graph::count() {
    int count = 0;
    for (int i = 0; i < M; i++) {
        if (verExi[i]) count++;
    }
    return count;
}