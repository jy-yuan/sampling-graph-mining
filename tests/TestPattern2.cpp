#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <set>
#include <utility>

#include "SamplerGenerator.hpp"
#include "NetworkInterface.hpp"
#include "Debug.hpp"

using namespace std;

VertexId num_vertices;
EdgeId num_edges;
Edge<Empty> * edges;
int curr_pos;
VertexId * mapping;
bool * flag;
Pattern ** patterns;

void dfs(VertexId idx) {
	if (idx == num_vertices) {
		patterns[curr_pos ++]->setup_vertex_id_mapping(mapping);
		//for (VertexId i = 0; i < num_vertices; ++ i) {
		//	printf("%u ", mapping[i]);
		//}
		//printf("\n");
		return ;
	}
	for (VertexId i = 0; i < num_vertices; ++ i) {
		if (! flag[i]) {
			flag[i] = true;
			mapping[idx] = i;
			dfs(idx + 1);
			flag[i] = false;
		}
	}
}

void verify(bool output_patterns = true) {
	int num_patterns = 1;
	for (VertexId i = 1; i <= num_vertices; ++ i) {
		num_patterns *= i;
	}
	patterns = new Pattern* [num_patterns];
	for (int i = 0; i < num_patterns; ++ i) {
		patterns[i] = new Pattern(num_vertices, num_edges, edges);
	}
	mapping = new VertexId [num_vertices];
	flag = new bool [num_vertices];
	memset(flag, false, sizeof(bool) * num_vertices);

	curr_pos = 0;
	dfs(0);

	bool *is_duplicated = new bool[num_patterns];
	memset(is_duplicated, false, sizeof(bool) * num_patterns);

	int num_unique_patterns = 1;
	for (int i = 1; i < num_patterns; ++ i) {
		for (int j = 0; j < i; ++ j) {
			if (patterns[i]->is_automorphic(*patterns[j])) {
				is_duplicated[i] = true;
				break;
			}
		}
		if (! is_duplicated[i]) {
			++ num_unique_patterns;
		}
	}
	if (output_patterns) {
		for (int i = 0; i < num_patterns; ++ i) {
			if (! is_duplicated[i]) {
				cout << *patterns[i] << endl;
			}
		}
	}
	cout << num_unique_patterns << "/"  << num_patterns << endl;

	Graph<Empty> graph;
	SamplerGenerator<Empty> * sampler_generator = new SamplerGenerator<Empty>(&graph, patterns[0]);
	delete sampler_generator;

	for (int i = 0; i < num_patterns; ++ i) {
		delete patterns[i];
	}
	delete [] is_duplicated;
	delete [] patterns;
	delete [] mapping;
	delete [] flag;
}

int main(int argc, char ** argv) {
	NetworkInterface::init_network_interface(&argc, &argv);

	// Four-Chain Pattern
	num_vertices = 4;
	num_edges = 3;
	edges = new Edge<Empty>[num_edges];
	for (EdgeId i = 0; i < num_edges; ++ i) {
		edges[i].src = i;
		edges[i].dst = i + 1;
		edges[i].edge_id = i;
	}
	verify();
	delete [] edges;

	// Four-Clique Pattern
	
	num_vertices = 4;
	num_edges = num_vertices * (num_vertices - 1) / 2;
	edges = new Edge<Empty>[num_edges];
	EdgeId cnt = 0;
	for (VertexId i = 0; i < num_vertices; ++ i) {
		for (VertexId j = 0; j < i; ++ j) {
			edges[cnt].src = i;
			edges[cnt].dst = j;
			edges[cnt].edge_id = cnt;
			++ cnt;
		}
	}
	verify();
	delete [] edges;
	
	// Five-Star Pattern
	num_vertices = 5;
	num_edges = num_vertices - 1;
	edges = new Edge<Empty>[num_edges];
	for (EdgeId e_i = 0; e_i < num_edges; ++ e_i) {
		edges[e_i].src = 0;
		edges[e_i].dst = e_i + 1;
		edges[e_i].edge_id = e_i;
	}
	verify();
	delete [] edges;

	//// random
	//srand(time(NULL));
	//for (int i = 0; i < 10; ++ i) {
	//	num_vertices = 6;
	//	num_edges = 12;
	//	set<pair<VertexId, VertexId>> edge_set;
	//	edges = new Edge<Empty> [num_edges];
	//	for (EdgeId e_i = 0; e_i < num_edges; ++ e_i) {
	//		VertexId src, dst;
	//		while (true) {
	//			src = rand() % num_vertices;
	//			dst = rand() % num_vertices;
	//			if (src == dst) continue;
	//			if (src > dst) swap(src, dst);
	//			if (edge_set.find(make_pair(src, dst)) == edge_set.end()) break;
	//		}
	//		edge_set.insert(make_pair(src, dst));
	//		edges[e_i].src = src;
	//		edges[e_i].dst = dst;
	//		edges[e_i].edge_id = e_i;
	//	}
	//	verify(false);
	//	delete [] edges;
	//}
	
	NetworkInterface::finalize_instance();
	return 0;
}
