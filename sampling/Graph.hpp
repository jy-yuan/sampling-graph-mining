#ifndef GRAPH_HPP
#define GRAPH_HPP

#pragma once
#include <cstdio>
#include <string>
#include <vector>

class Graph {
    std::vector<int> vertexes;
    std::vector<int> verExi;
    std::vector<int> verDeg;
    std::vector<int> csrInd;
    std::vector<int> csrList;

   public:
    Graph(){};
    ~Graph(){};
    int init(int m);
    int init_from_file(const std::string dir = "graph");
    int join(int* zipgraph);
    int* sample(int num);
    int count();
    int m, source;  // for sampling
    int M, N;
};

#endif