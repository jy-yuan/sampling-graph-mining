#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <string.h>

using namespace std;

typedef unsigned int VertexId;
typedef unsigned long long EdgeId;

struct EdgeUnit {
	VertexId src;
	VertexId dst;
	EdgeUnit() {
	}
	EdgeUnit(VertexId src_, VertexId dst_): src(src_), dst(dst_) {
	}
};

EdgeUnit * read_edge_list;

struct Edge {
	VertexId dst;
	Edge * next;
	EdgeId edge_num;
};

const int num_partitions = 4;

vector<EdgeUnit> subgraph_edge_list;
VertexId subgraph_vertices;
EdgeId subgraph_edges;

VertexId vertices;
EdgeId edges;
Edge ** edge_list;

VertexId * mapping;
VertexId * partition_offset;

inline int get_partition_id(VertexId v_i) {
	for (int p_i = 0; p_i < num_partitions; ++ p_i) {
		if (v_i >= partition_offset[p_i] && v_i < partition_offset[p_i + 1]) {
			return p_i;
		}
	}
	assert(false);
}

inline void add_edge(VertexId src, VertexId dst, EdgeId edge_num) {
	Edge * edge_unit = new Edge();
	edge_unit->dst = dst;
	edge_unit->next = edge_list[src];
	edge_unit->edge_num = edge_num;
	edge_list[src] = edge_unit;
}

double sample_edge(EdgeId &sampled_edge_num) { 
	if (subgraph_edges == 0) {
		return -1;
	}
	sampled_edge_num = rand() % subgraph_edges;
	return 1. / double(subgraph_edges);
}

double conditional_sample_edge(vector<EdgeId> &sampled_edges_num, EdgeId &sampled_edge_num) { 
	EdgeId last_edge_num = 0;
	set<VertexId> sampled_vertices_set;
	vector<VertexId> sampled_vertices;
	for (EdgeId i = 0, j = sampled_edges_num.size(); i < j; ++ i) {
		EdgeId edge_num = sampled_edges_num[i];
		if (edge_num > last_edge_num) {
			last_edge_num = edge_num;
		}
		VertexId sampled_vertex = subgraph_edge_list[edge_num].src;
		if (sampled_vertices_set.find(sampled_vertex) == sampled_vertices_set.end()) {
			sampled_vertices_set.insert(sampled_vertex);
			sampled_vertices.push_back(sampled_vertex);
		}
		sampled_vertex = subgraph_edge_list[edge_num].dst;
		if (sampled_vertices_set.find(sampled_vertex) == sampled_vertices_set.end()) {
			sampled_vertices_set.insert(sampled_vertex);
			sampled_vertices.push_back(sampled_vertex);
		}
	}

	set<EdgeId> edges_to_be_sampled_set;
	vector<EdgeId> edges_to_be_sampled;
	for (VertexId i = 0, j = sampled_vertices.size(); i < j; ++ i) {
		VertexId sampled_vertex = sampled_vertices[i];
		for (Edge * p = edge_list[sampled_vertex]; p; p = p->next) {
			if (p->edge_num > last_edge_num) {
				if (edges_to_be_sampled_set.find(p->edge_num) == edges_to_be_sampled_set.end()) {
					edges_to_be_sampled_set.insert(p->edge_num);
					edges_to_be_sampled.push_back(p->edge_num);
				}
			}
		}
	}
	
	if (edges_to_be_sampled.empty()) {
		return -1;
	}
	sampled_edge_num = edges_to_be_sampled[rand() % edges_to_be_sampled.size()];
	return 1. / double(edges_to_be_sampled.size());
}

double estimate() {
	EdgeId sampled_edge;
	vector<EdgeId> sampled_edges;
	double p0 = sample_edge(sampled_edge);
	if (p0 < 0) return 0;
	sampled_edges.push_back(sampled_edge);
	double p1 = conditional_sample_edge(sampled_edges, sampled_edge);
	if (p1 < 0) return 0;
	return 1. / (p0 * p1);
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << "./ThreeChainApprox [graph]" << endl;
		exit(-1);
	}

	srand(time(NULL));
	char * graph_name = argv[1];
	ifstream fin(graph_name);
	fin >> vertices >> edges; 
	edge_list = new Edge*[vertices];
	read_edge_list = new EdgeUnit[edges];
	mapping = new VertexId[vertices];
	partition_offset = new VertexId[num_partitions + 1];

	// shuffle the graph
	for (VertexId v_i = 0; v_i < vertices; ++ v_i) {
		mapping[v_i] = v_i;
	}
	for (VertexId v_i = vertices - 1; v_i; -- v_i) {
		swap(mapping[v_i], mapping[rand() % (v_i + 1)]);
	}
	for (int p_i = 0; p_i < num_partitions; ++ p_i) {
		partition_offset[p_i] = vertices / num_partitions * p_i;
	}
	partition_offset[num_partitions] = vertices;

	for (EdgeId e_i = 0; e_i < edges; ++ e_i) {
		VertexId src, dst;
		fin >> src >> dst;
		//add_edge(src, dst, e_i);
		//add_edge(dst, src, e_i);
		read_edge_list[e_i].src = mapping[src];
		read_edge_list[e_i].dst = mapping[dst];
	}
	fin.close();

	EdgeId num_estimators = 2 << 18;
	double global_reducer = 0;
	for (int p_i = 0; p_i < num_partitions; ++ p_i) {
		memset(edge_list, 0, sizeof(Edge*) * vertices);
		subgraph_edge_list.clear();
		subgraph_edges = 0;
		subgraph_vertices = partition_offset[p_i + 1] - partition_offset[p_i];
		for (EdgeId e_i = 0; e_i < edges; ++ e_i) {
			VertexId src = read_edge_list[e_i].src;
			VertexId dst = read_edge_list[e_i].dst;
			if (get_partition_id(src) == p_i && get_partition_id(dst) == p_i) {
				add_edge(src, dst, subgraph_edges);
				add_edge(dst, src, subgraph_edges);
				++ subgraph_edges;
				subgraph_edge_list.push_back(EdgeUnit(src, dst));
			}
		}
		double reducer = 0;
		//EdgeId local_num_estimators = num_estimators / num_partitions;
		//if (p_i == num_partitions - 1) {
		//	local_num_estimators += num_estimators % num_partitions;
		//}
		cout << "partition " << p_i << endl;
		cout << "owned edges: " << subgraph_edge_list.size() << endl; 
		for (EdgeId i = 0; i < num_estimators; ++ i) {
			cout << "estimating " << double(i) / double(num_estimators) << endl;
			cout << "\033[F";
			reducer += estimate();
		}
		cout << endl;
		cout << "estimation: " << EdgeId(reducer / num_estimators) << endl;
		cout << endl;
		global_reducer += reducer;
	}
	cout << "estimated number of three-chain: " << EdgeId(global_reducer / num_estimators * num_partitions * num_partitions) << endl;

	//double cnt = 0;
	//for (EdgeId i = 0; i < num_estimators; ++ i) {
	//	cout << "estimating " << double(i) / double(num_estimators) << endl;
	//	cout << "\033[F";
	//	cnt += estimate();
	//}
	//cout << endl;
	//cnt /= double(num_estimators);
	//
	//cout << "estimated number of three-chain: " << EdgeId(cnt) << endl;

	delete [] edge_list;
	delete [] read_edge_list;
	delete [] mapping;
	delete [] partition_offset;
	return 0;
}

