#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <iostream>

#include "SamplerGenerator.hpp"

using namespace std;

Pattern ** patterns;
VertexId num_vertices;
EdgeId num_edges;
Edge<Empty> *edges;

EdgeId *mapping;
bool *flag;
EdgeId curr_pos;

void dfs(EdgeId idx) {
	if (idx == num_edges) {
		//for (EdgeId e_i = 0; e_i < num_edges; ++ e_i) {
		//	printf("%lu ", mapping[e_i]);
		//}
		//printf("\n");
		patterns[curr_pos] = new Pattern(num_vertices, num_edges, edges);
		patterns[curr_pos]->setup_edge_id_mapping(mapping);
		++ curr_pos;
		return ;
	}
	for (EdgeId e_i = 0; e_i < num_edges; ++ e_i) {
		if (flag[e_i] == false) {
			flag[e_i] = true;
			mapping[idx] = e_i;
			dfs(idx + 1);
			flag[e_i] = false;
		}
	}
}

void verify() {
	mapping = new EdgeId [num_edges];
	EdgeId num_patterns = 1;
	for (EdgeId i = 1; i <= num_edges; ++ i) {
		num_patterns *= i;
	}
	printf("num_patterns = %lu\n", num_patterns);
	patterns = new Pattern* [num_patterns];
	flag = new bool [num_edges];
	memset(flag, false, sizeof(bool) * num_edges);
	curr_pos = 0;
	dfs(0);
	assert(curr_pos == num_patterns);
	//for (EdgeId p_i = 0; p_i < num_patterns; ++ p_i) {
	//	cout << *patterns[p_i] << endl;
	//}
	printf("cleaning duplicated patterns...\n");
	bool *is_duplicated = new bool[num_patterns];
	memset(is_duplicated, false, sizeof(bool) * num_patterns);
	for (EdgeId p_i = 1; p_i < num_patterns; ++ p_i) {
		for (EdgeId p_j = 0; p_j < p_i; ++ p_j) {
			if (is_duplicated[p_j] == false) {
				if (patterns[p_i]->has_same_edge_ordering(*patterns[p_j])) {
					is_duplicated[p_i] = true;
					break;
				}
			}
		}
	}
	EdgeId cnt = 0;
	for (EdgeId p_i = 0; p_i < num_patterns; ++ p_i) {
		if (! is_duplicated[p_i]) {
			cout << *patterns[p_i] << endl;
			++ cnt;
		}
	}
	printf("%lu/%lu\n", cnt, num_patterns);
	delete [] flag;
	delete [] mapping;
	delete [] is_duplicated;
	for (EdgeId p_i = 0; p_i < num_patterns; ++ p_i) {
		delete patterns[p_i];
	}
	delete [] patterns;
}

int main(int argc, char ** argv) {
	// five-chain
	printf("********** five-chain patterns: **********\n");
	num_vertices = 5;
	num_edges = 4;
	edges = new Edge<Empty> [num_edges];
	for (VertexId v_i = 1; v_i < num_vertices; ++ v_i) {
		edges[v_i - 1].src = v_i - 1;
		edges[v_i - 1].dst = v_i;
		edges[v_i - 1].edge_id = v_i - 1;
	}
	for (EdgeId e_i = 0; e_i < num_edges; ++ e_i) {
		printf("(%u, %u: %lu)\n", edges[e_i].src, edges[e_i].dst, edges[e_i].edge_id);
	}
	verify();
	delete [] edges;

	// 4-clique
	printf("********** four-clique patterns: **********\n");
	num_vertices = 4;
	num_edges = 4 * 3 / 2;
	edges = new Edge<Empty> [num_edges];
	{
		EdgeId tmp = 0;
		for (VertexId v_i = 0; v_i < num_vertices; ++ v_i) {
			for (VertexId v_j = 0; v_j < v_i; ++ v_j) {
				edges[tmp].src = v_i;
				edges[tmp].dst = v_j;
				edges[tmp].edge_id = tmp;
				++ tmp;
			}
		}
		assert(tmp == num_edges);
	}
	verify();
	delete [] edges;

	return 0;
}
