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
		cout << "./FourChain [graph]" << endl;
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
	EdgeId four_chain_count = 0;
	for (VertexId v_0 = 0; v_0 < vertices; ++ v_0) { // the first vertex
		cout << "counting... " << double(v_0) / double(vertices) << endl;
		cout << "\033[F";
		for (Edge * p_0 = edge_list[v_0]; p_0; p_0 = p_0->next) {
			VertexId v_1 = p_0->dst; // the second vertex
			for (Edge * p_1 = edge_list[v_1]; p_1; p_1 = p_1->next) {
				VertexId v_2 = p_1->dst; // the third vertex
				if (v_2 != v_0) {
					for (Edge * p_2 = edge_list[v_2]; p_2; p_2 = p_2->next) {
						VertexId v_3 = p_2->dst; // the fourth vertex
						if (v_3 != v_0 && v_3 != v_1) {
							++ four_chain_count;
						}
					}
				}
			}
		}
	}
	cout << endl;
	cout << "count: " << four_chain_count / 2 << endl;
	return 0;
}

