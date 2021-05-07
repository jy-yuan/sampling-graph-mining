#ifndef GRAPH_HPP
#define GRAPH_HPP

#pragma once
#include <cstdio>
#include <string>
#include <vector>

class Graph {
    std::vector<int> verExi;
    std::vector<int> vertexes;
    std::vector<int> csrInd;
    std::vector<int> csrList;
    int M, N;

   public:
    Graph(){};
    ~Graph(){};
    int init(const std::string dir = "graph");
    int join(int* zipgraph);
    int* sample(int num);
    int count();
    int m, source;  // for sampling
};

#endif