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
		cout << "./FourClique [graph]" << endl;
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
	EdgeId count = 0;
	for (VertexId v_i = 0; v_i < vertices; ++ v_i) {
		if (v_i%1000 == 0){
		cout << "counting... " << double(v_i) / double(vertices) << endl;
		cout << "\033[F";
		}
		for (Edge * e_i = edge_list[v_i]; e_i; e_i = e_i->next) {
			VertexId v_j = e_i->dst;
			if (v_j > v_i) {
				for (Edge * e_j = edge_list[v_j]; e_j; e_j = e_j->next) {
					VertexId v_k = e_j->dst;
					if (v_k > v_j) {
						bool flag = false;
						for (Edge * e_k = edge_list[v_k]; e_k; e_k = e_k->next) {
							if (e_k->dst == v_i) {
								flag = true;
								break;
							}
						}
						if (flag) {
							for (VertexId v_z = v_k + 1; v_z < vertices; ++ v_z) {
								int ct = 0;
								for (Edge * e_z = edge_list[v_z]; e_z; e_z = e_z->next) {
									if (e_z->dst == v_i) ++ ct;
									if (e_z->dst == v_j) ++ ct;
									if (e_z->dst == v_k) ++ ct;
								}
								if (ct == 3) ++ count;
							}
						}
					}
				}
			}
		}
	}
	cout << endl;
	cout << "count: " << count << endl;
	return 0;
}

