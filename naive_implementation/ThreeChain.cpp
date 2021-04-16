#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>

using namespace std;

typedef unsigned int VertexId;
typedef unsigned long long EdgeId;

struct Edge {
	VertexId dst;
	Edge * next;
};

VertexId vertices;
EdgeId edges;
Edge ** edge_list;

inline void add_edge(VertexId src, VertexId dst) {
	Edge * edge_unit = new Edge();
	edge_unit->dst = dst;
	edge_unit->next = edge_list[src];
	edge_list[src] = edge_unit;
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << "./ThreeChain [graph]" << endl;
		exit(-1);
	}
	char * graph_name = argv[1];
	ifstream fin(graph_name);
	fin >> vertices >> edges; 
	edge_list = new Edge*[vertices];
	for (VertexId e_i = 0; e_i < edges; ++ e_i) {
		VertexId src, dst;
		fin >> src >> dst;
		add_edge(src, dst);
		add_edge(dst, src);
	}
	fin.close();
	EdgeId three_chain_count = 0;
	for (VertexId v_i = 0; v_i < vertices; ++ v_i) {
		EdgeId degree = 0;
		for (Edge * p_i = edge_list[v_i]; p_i; p_i = p_i->next) {
			++ degree;
		}
		three_chain_count += degree * (degree - 1) / 2;
	}
	cout << "count: " << three_chain_count << endl;
	return 0;
}

