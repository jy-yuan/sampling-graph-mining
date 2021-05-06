#pragma once
#include <vector>
#include <string>
#include <cstdio>

class Graph {
    std::vector<int> verExi;
    std::vector<int> csrInd;
    std::vector<int> csrList;
    int M, N;
  public:
    Graph() {};
    ~Graph() {};
    int init(const std::string dir = "graph");
    int join(int* zipgraph);
    int* sample(int num);
    int count();
};