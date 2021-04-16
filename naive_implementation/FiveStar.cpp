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
		cout << "./FiveStar [graph]" << endl;
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

	cout << "start counting..." << endl;
	EdgeId count = 0;
	for (VertexId v_0 = 0; v_0 < vertices; ++ v_0) { // the centric vertex
		if (v_0%10000==0) {
		cout << "counting... " << double(v_0) / double(vertices) << endl;
		cout << "\033[F";}
		EdgeId tmp = 0;
		for (Edge * p_0 = edge_list[v_0]; p_0; p_0 = p_0->next) {
			++ tmp;
		}
		if (tmp >= 4) {
			EdgeId delta = tmp * (tmp - 1) * (tmp - 2) * (tmp - 3);
			delta /= (4 * 3 * 2 * 1);
			count += delta;
		}
	}
	cout << endl;
	cout << "count: " << count << endl;
	return 0;
}

