#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <stdio.h>

using namespace std;

typedef unsigned int VertexId;
typedef unsigned long long EdgeId;

struct Edge {
	VertexId dst;
	VertexId value;
	Edge * next;
};

VertexId vertices;
EdgeId edges;
Edge ** edge_list;

inline void add_edge(VertexId src, VertexId dst, VertexId value) {
	Edge * edge_unit = new Edge();
	edge_unit->dst = dst;
	edge_unit->value = value;
	edge_unit->next = edge_list[src];
	edge_list[src] = edge_unit;
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << "./TriangleCounting [graph]" << endl;
		exit(-1);
	}
	char * graph_name = argv[1];
	ifstream fin(graph_name);
	fin >> vertices >> edges; 
	edge_list = new Edge*[vertices];
	for (VertexId e_i = 0; e_i < edges; ++ e_i) {
		VertexId src, dst, value;
		fin >> src >> dst >> value;
		add_edge(src, dst, value);
		add_edge(dst, src, value);
	}
	fin.close();
	EdgeId triangel_count = 0;
	for (VertexId v_i = 0; v_i < vertices; ++ v_i) {
		cout << "counting... " << double(v_i) / double(vertices) << endl;
		cout << "\033[F";
		for (Edge * e_i = edge_list[v_i]; e_i; e_i = e_i->next) {
			VertexId v_j = e_i->dst;
			if (v_j > v_i) {
				for (Edge * e_j = edge_list[v_j]; e_j; e_j = e_j->next) {
					VertexId v_k = e_j->dst;
					if (v_k > v_j) {
						bool flag = false;
						for (Edge * e_k = edge_list[v_k]; e_k; e_k = e_k->next) {
							if (e_i->value >= 50 || e_j->value >= 50 || e_k->value >= 50) {
								if (e_k->dst == v_i) {
									flag = true;
									break;
								}
							}
						}
						if (flag) {
							++ triangel_count;
						}
					}
				}
			}
		}
	}
	cout << endl;
	cout << "count: " << triangel_count << endl;
	return 0;
}

