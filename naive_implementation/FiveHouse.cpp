#include <iostream>
#include <fstream>
#include <set>
#include <utility>

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <numa.h>
#include <assert.h>

using namespace std;

typedef unsigned int VertexId;
typedef unsigned long long EdgeId;

VertexId vertices;
EdgeId edges;
pair<VertexId, VertexId> *edge_list;
//set< pair<VertexId, VertexId> > * edge_sets;

VertexId * csr_list;
EdgeId * csr_idx;

EdgeId check(VertexId v_left, VertexId v_right, VertexId v_other) {
	EdgeId cnt = 0;
	for (EdgeId idx_left = csr_idx[v_left]; idx_left < csr_idx[v_left + 1]; ++ idx_left) {
		VertexId v_left_extend = csr_list[idx_left];
		if (v_left_extend != v_left && v_left_extend != v_right && v_left_extend != v_other) {
#pragma omp parallel for reduction(+:cnt)
			for (EdgeId idx_right = csr_idx[v_right]; idx_right < csr_idx[v_right + 1]; ++ idx_right) {
				EdgeId local_cnt = 0;
				VertexId v_right_extend = csr_list[idx_right];
				if (v_right_extend != v_left && v_right_extend != v_right && v_right_extend != v_other && v_right_extend != v_left_extend) {
					bool found = false;
					for (EdgeId idx = csr_idx[v_right_extend]; idx < csr_idx[v_right_extend + 1]; ++ idx) {
						if (csr_list[idx] == v_left_extend) {
							found = true;
							break;
						}
					}
					if (found) ++ local_cnt;
				}
				cnt += local_cnt;
			}
		}
	}
	return cnt;
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << "./FiveHouse [graph]" << endl;
		exit(-1);
	}
	int num_threads = numa_num_configured_cpus(); 
	omp_set_dynamic(0);
	omp_set_num_threads(num_threads);
	//edge_sets = new set< pair<VertexId, VertexId> > [num_threads];
	char * graph_name = argv[1];
	ifstream fin(graph_name);
	fin >> vertices >> edges; 
	edge_list = new pair<VertexId, VertexId>[edges];
	csr_idx = new EdgeId[vertices + 1];
	memset(csr_idx, 0, sizeof(EdgeId) * vertices);
	for (VertexId e_i = 0; e_i < edges; ++ e_i) {
		VertexId src, dst;
		fin >> src >> dst;
		edge_list[e_i].first = src;
		edge_list[e_i].second = dst;
		csr_idx[src + 1] ++;
		csr_idx[dst + 1] ++;
	}
	//cout << "A" << endl;
	for (VertexId v_i = 1; v_i <= vertices; ++ v_i) {
		csr_idx[v_i] += csr_idx[v_i - 1];
	}
	EdgeId * curr_pos = new EdgeId[vertices + 1];
	for (VertexId v_i = 0; v_i <= vertices; ++ v_i) {
		curr_pos[v_i] = csr_idx[v_i];
	}
	//cout << "B" << endl;
	csr_list = new VertexId[edges * 2];
	for (EdgeId e_i = 0; e_i < edges; ++ e_i) {
		VertexId src = edge_list[e_i].first;
		VertexId dst = edge_list[e_i].second;
		EdgeId pos = curr_pos[src] ++;
		csr_list[pos] = dst;
		pos = curr_pos[dst] ++;
		csr_list[pos] = src;
	}
	//cout << "C" << endl;
	for (VertexId v_i = 0; v_i < vertices; ++ v_i) {
		assert(curr_pos[v_i] == csr_idx[v_i + 1]);
	}
	fin.close();
	EdgeId reducer = 0;
	cout << "num_threads: " << num_threads << endl;
	for (VertexId v_i = 0; v_i < vertices; ++ v_i)
	{
			//if (v_i % 1 == 0 || v_i == vertices - 1) {
			      cout << "counting... " << double(v_i) / double(vertices) << endl;
			      cout << "\033[F";
			//}
			for (EdgeId idx_i = csr_idx[v_i]; idx_i < csr_idx[v_i + 1]; ++ idx_i) {
				VertexId v_j = csr_list[idx_i];
				if (v_j > v_i) {
					for (EdgeId idx_j = csr_idx[v_j]; idx_j < csr_idx[v_j + 1]; ++ idx_j) {
						EdgeId local_reducer = 0;
						VertexId v_k = csr_list[idx_j];
						if (v_k > v_j) {
							bool found = false;
							for (EdgeId idx_k = csr_idx[v_k]; idx_k < csr_idx[v_k + 1]; ++ idx_k) {
								if (csr_list[idx_k] == v_i) {
									found = true;
									break;
								}
							}
							if (found) {
								reducer += check(v_i, v_j, v_k);
								reducer += check(v_i, v_k, v_j);
								reducer += check(v_j, v_k, v_i);
							}
						}
					}
				}
			}
	}
	cout << endl;
	cout << "count: " << reducer << endl;
	return 0;
}

