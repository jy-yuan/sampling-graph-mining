#include <assert.h>
#include <numa.h>
#include <omp.h>

#include <algorithm>

#include "ZGraphInstance.hpp"
#include "SharedMemSys.hpp"
#include "Time.hpp"

#define BINARY_SEARCH_THRESHOLD 32

// Estimation

Estimation::Estimation() {
	reset();
}

void Estimation::reset() {
	num_finished_estimators = 0;
	current_estimation = 0;
}

void Estimation::update(unsigned long long new_num_estimators, double new_estimation) {
	current_estimation = (1. * num_finished_estimators / (num_finished_estimators + new_num_estimators)) * current_estimation + (1. * new_num_estimators / (num_finished_estimators + new_num_estimators)) * new_estimation;
	num_finished_estimators += new_num_estimators;
}

unsigned long long Estimation::get_num_finished_estimators() {
	return num_finished_estimators;
}

double Estimation::get_currnet_estimation() {
	return current_estimation;
}

// SubGraph

template<typename EdgeData>
SubGraph<EdgeData>::SubGraph() {
	min_valid_edge_id = 0;
	min_valid_vertex_id = 0;
}

template<typename EdgeData>
bool SubGraph<EdgeData>::is_in_vertices_list(VertexId vertex) {
	for (int i = 0, j = vertices.list_size; i < j; ++ i) {
		if (vertices.list[i] == vertex) {
			return true;
		}
	}
	return false;
}

template<typename EdgeData> 
void SubGraph<EdgeData>::add_edge(const Edge<EdgeData>& edge) {
	edges.push_back(edge); 
	if (edge.edge_id >= min_valid_edge_id) {
		min_valid_edge_id = edge.edge_id + 1;
	}
	add_vertex(edge.src);
	add_vertex(edge.dst);
}

template<typename EdgeData>
void SubGraph<EdgeData>::add_vertex(VertexId vertex) {
	if (! is_in_vertices_list(vertex)) {
		vertices.push_back(vertex); 
		if (vertex >= min_valid_vertex_id) {
			min_valid_vertex_id = vertex + 1;
		}
	}
}

template<typename EdgeData>
const MyList<Edge<EdgeData>>& SubGraph<EdgeData>::get_edges() const {
	return edges;
}

template<typename EdgeData>
const MyList<VertexId>& SubGraph<EdgeData>::get_vertices() const {
	return vertices;
}

template<typename EdgeData>
EdgeId SubGraph<EdgeData>::get_min_valid_edge_id() const {
	return min_valid_edge_id;
}

template<typename EdgeData>
VertexId SubGraph<EdgeData>::get_min_valid_vertex_id() const {
	return min_valid_vertex_id;
}

template<typename EdgeData> 
SubGraph<EdgeData> SubGraph<EdgeData>::operator-(const SubGraph<EdgeData> &sub_graph) {
	SubGraph<EdgeData> result;
	for (int i = 0; i < edges.list_size; ++ i) {
		const Edge<EdgeData> &edge_0 = edges.list[i];
		bool is_delete = false;
		for (int j = 0; j < sub_graph.edges.list_size; ++ j) {
			const Edge<EdgeData> &edge_1 = sub_graph.edges.list[j];
			if ((edge_0.src == edge_1.src && edge_0.dst == edge_1.dst) || (edge_0.src == edge_1.dst && edge_0.dst == edge_1.src)) {
				is_delete = true;
				break;
			}
		}
		if (! is_delete) {
			result.add_edge(edge_0);
		}
	}
	return result;
}

// ZGraphInstance

template<typename EdgeData>
double ZGraphInstance<EdgeData>::get_scaling_factor(VertexId num_vertices_in_pattern) {
	double factor = 1.;
	double num_workers = graph->p.num_sockets * graph->p.num_partitions;
	for (VertexId i = 1; i < num_vertices_in_pattern; ++ i) {
		factor *= num_workers;
	}
	return factor;
}

template<typename EdgeData>
std::pair<Edge<EdgeData>, double> ZGraphInstance<EdgeData>::core_conditional_sample_edge(const MyList<Edge<EdgeData>> &sampled_edges, const MyList<VertexId> &sampled_vertices) {
	int t_i = SharedMemSys::get_instance()->get_current_thread_id();
	int s_i = SharedMemSys::get_instance()->get_socket_id(t_i);

	EdgeId min_valid_edge_id = 0;
	for (int e_i = 0; e_i < sampled_edges.list_size; ++ e_i) {
		if (sampled_edges.list[e_i].edge_id + 1 > min_valid_edge_id) {
			min_valid_edge_id = sampled_edges.list[e_i].edge_id + 1;
		}
	}

	Edge<EdgeData> sampled_edge;
	double probability;

	EdgeId sum = 0;
	MyList<VertexId> current_pos;
	for (int i = 0; i < sampled_vertices.list_size; ++ i) {
		VertexId v_i = sampled_vertices.list[i];
		if (graph->outgoing_adj_bitmap[s_i]->get_bit(v_i)) {
			EdgeId begin_idx = graph->outgoing_adj_index[s_i][v_i];
			EdgeId end_idx = graph->outgoing_adj_index[s_i][v_i + 1];
			if (min_valid_edge_id <= graph->outgoing_adj_list[s_i][begin_idx].edge_id) { // all associated edges are valid
				sum += (end_idx - begin_idx);
				current_pos.push_back(begin_idx);
			} else if (min_valid_edge_id > graph->outgoing_adj_list[s_i][end_idx - 1].edge_id) { // all associated edges are invalid
				current_pos.push_back(end_idx);
			} else {
				if (end_idx - begin_idx <= BINARY_SEARCH_THRESHOLD) { // brute-force
					for (EdgeId idx = begin_idx; idx < end_idx; ++ idx) {
						if (min_valid_edge_id <= graph->outgoing_adj_list[s_i][idx].edge_id) {
							sum += (end_idx - idx);
							current_pos.push_back(idx);
							break;
						}
					}
				} else { // binary search
					EdgeId left_idx = begin_idx; // invariance: outgoing_adj_list[s_i][left_idx].edge_id < min_valid_edge_id
					EdgeId right_idx = end_idx - 1; // invariance: outgoing_adj_list[s_i][right_id].edge_id >= min_valid_edge_id
					EdgeId mid_idx;
					while (right_idx - left_idx > 1) {
						mid_idx = (left_idx + right_idx) >> 1;
						if (graph->outgoing_adj_list[s_i][mid_idx].edge_id < min_valid_edge_id) {
							left_idx = mid_idx;
						} else {
							right_idx = mid_idx;
						}
					}
					sum += (end_idx - right_idx);
					current_pos.push_back(right_idx);
				}
			}
		}
	}

	if (sum == 0) { // no valid edges
		return std::make_pair(sampled_edge, (double) -1);
	}

	EdgeId sampled_number = randoms[t_i].rand_int(0, sum - 1);
	int current_pos_idx = 0;
	bool is_successfully_sampled = false;
	for (int i = 0; i < sampled_vertices.list_size; ++ i) {
		VertexId v_i = sampled_vertices.list[i];
		if (graph->outgoing_adj_bitmap[s_i]->get_bit(v_i)) {
			EdgeId curr_idx = current_pos.list[current_pos_idx ++];
			EdgeId end_idx = graph->outgoing_adj_index[s_i][v_i + 1];
			if (curr_idx + sampled_number < end_idx) { // found sampled edge
				curr_idx += sampled_number;
				sampled_edge.edge_id = graph->outgoing_adj_list[s_i][curr_idx].edge_id;
				sampled_edge.src = graph->outgoing_adj_list[s_i][curr_idx].src;
				sampled_edge.dst = graph->outgoing_adj_list[s_i][curr_idx].dst;
				if (! std::is_same<EdgeData, Empty>::value) {
					sampled_edge.edge_data = graph->outgoing_adj_list[s_i][curr_idx].edge_data;
				}
				is_successfully_sampled = true;
				break;
			} else {
				sampled_number -= (end_idx - curr_idx);
			}
		}
	}
	double edge_appear_count = 0;
	for (int i = 0; i < sampled_vertices.list_size; ++ i) {
		if (sampled_vertices.list[i] == sampled_edge.src) {
			++ edge_appear_count;
		}
		if (sampled_vertices.list[i] == sampled_edge.dst) {
			++ edge_appear_count;
		}
	}
	assert(is_successfully_sampled);
	assert(edge_appear_count == 1 || edge_appear_count == 2);
	probability = 1. * edge_appear_count / sum;
	return std::make_pair(sampled_edge, probability);
}

template<typename EdgeData>
ZGraphInstance<EdgeData>::ZGraphInstance() {
	graph = nullptr;
	randoms = new Random[SharedMemSys::get_instance()->get_num_threads()];
	is_counting_sample_edge_calls = false;
	is_counting_conditional_edge_calls = false;
	cnt = 0;
}

template<typename EdgeData>
ZGraphInstance<EdgeData>::ZGraphInstance(Graph<EdgeData> * graph_): graph(graph_) {
	randoms = new Random[SharedMemSys::get_instance()->get_num_threads()]; // I am not sure that they will produce perfect random intergers, should be tested
	is_counting_sample_edge_calls = false;
	is_counting_conditional_edge_calls = false;
	cnt = 0;
}

template<typename EdgeData>
ZGraphInstance<EdgeData>::~ZGraphInstance() {
	delete [] randoms;
}

template<typename EdgeData>
Graph<EdgeData>* ZGraphInstance<EdgeData>::get_graph() {
	return graph;
}

template<typename EdgeData>
void ZGraphInstance<EdgeData>::set_graph(Graph<EdgeData> *graph_) {
	graph = graph_;
}

template<typename EdgeData>
double ZGraphInstance<EdgeData>::calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation) {
	//Debug::get_instance()->enter_function("ZGraphInstance::calulate");
	
	if (init_estimation.get_num_finished_estimators() >= num_estmators) {
		//Debug::get_instance()->leave_function("ZGraphInstance::calulate");
		return init_estimation.get_currnet_estimation();
	}

	num_estmators -= init_estimation.get_num_finished_estimators();
	int num_threads = SharedMemSys::get_instance()->get_num_threads();
	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	int num_threads_per_socket = num_threads / num_sockets;

	//Debug::get_instance()->log("preparing to estimate...");

	double reducer = 0;
#pragma omp parallel reduction(+:reducer)
	{
		double local_reducer = 0;
		int t_i = SharedMemSys::get_instance()->get_current_thread_id();
		unsigned local_num_estimators = num_estmators / num_threads_per_socket;
		if (SharedMemSys::get_instance()->get_socket_offset(t_i) == num_threads_per_socket - 1) {
			local_num_estimators += num_estmators % num_threads_per_socket;
		}
		//Debug::get_instance()->log("local_num_estimators(t_i = ", t_i, "): ", local_num_estimators);
		for (unsigned long long i = 0; i < local_num_estimators; ++ i) {
			local_reducer += estimate();
		}
		local_reducer /= double(num_estmators);
		reducer += local_reducer;
	}

	double global_reducer;
	MPI_Allreduce(&reducer, &global_reducer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double estimation = global_reducer * get_scaling_factor(num_vertices_in_pattern);
	init_estimation.update(num_estmators, estimation);

	//Debug::get_instance()->leave_function("ZGraphInstance::calulate");

	return init_estimation.get_currnet_estimation();
}

template<typename EdgeData>
std::pair<VertexId, double> ZGraphInstance<EdgeData>::sample_vertex() {
	int t_i = SharedMemSys::get_instance()->get_current_thread_id(); // could be optimized TODO
	int s_i = SharedMemSys::get_instance()->get_socket_id(t_i);
	if (graph->p.local_partition_offset[s_i] == graph->p.local_partition_offset[s_i + 1]) { // which means that there is no vertices in this sub-partition
		return std::make_pair<VertexId, double>((VertexId) 0, (double) -1);
	}
	VertexId min_v = graph->p.local_partition_offset[s_i];
	VertexId max_v = graph->p.local_partition_offset[s_i + 1] - 1;
	return std::make_pair<VertexId, double>(
			(VertexId) randoms[t_i].rand_int(min_v, max_v),
			(double) 1. / (max_v - min_v + 1)
			);
}

template<typename EdgeData>
std::pair<Edge<EdgeData>, double> ZGraphInstance<EdgeData>::sample_edge() {
	if (is_counting_sample_edge_calls) {
		++ cnt;
	}
	int t_i = SharedMemSys::get_instance()->get_current_thread_id();
	int s_i = SharedMemSys::get_instance()->get_socket_id(t_i);
	Edge<EdgeData> edge;
	if (graph->num_local_edges[s_i] == 0) { // which means that there is no edges in this sub-partition
		return std::make_pair(edge, (double) -1);
	}
	EdgeId sampled_edge_id = randoms[t_i].rand_int(0, graph->num_local_edges[s_i] - 1);
	edge.edge_id = sampled_edge_id;
	edge.src = graph->local_edges[s_i][sampled_edge_id].src; 
	edge.dst = graph->local_edges[s_i][sampled_edge_id].dst;
	if (! std::is_same<EdgeData, Empty>::value) {
		edge.edge_data = graph->local_edges[s_i][sampled_edge_id].edge_data;
	}
	double probability = 1. / graph->num_local_edges[s_i];
	return std::make_pair(edge, probability);
}

template<typename EdgeData>
std::pair<VertexId, double> ZGraphInstance<EdgeData>::conditional_sample_vertex(const SubGraph<EdgeData> &sub_graph) {
	int t_i = SharedMemSys::get_instance()->get_current_thread_id();
	int s_i = SharedMemSys::get_instance()->get_socket_id(t_i);
	VertexId min_valid_vertex_id = sub_graph.get_min_valid_vertex_id();
	if (min_valid_vertex_id >= graph->p.local_partition_offset[s_i + 1]) { // there is no valid vertex
		return std::make_pair((VertexId) 0, (double) -1);
	}
	VertexId sampled_v = randoms[t_i].rand_int(
			min_valid_vertex_id,
			graph->p.local_partition_offset[s_i + 1] - 1
			);
	double probability = 1. / (graph->p.local_partition_offset[s_i + 1] - min_valid_vertex_id);
	return std::make_pair(sampled_v, probability);
}

template<typename EdgeData>
std::pair<Edge<EdgeData>, double> ZGraphInstance<EdgeData>::conditional_sample_edge(const SubGraph<EdgeData> &sub_graph) {
	if (is_counting_conditional_edge_calls) {
		++ cnt;
	}
	const MyList<Edge<EdgeData>>& sampled_edges = sub_graph.get_edges();
	const MyList<VertexId>& sampled_vertices = sub_graph.get_vertices();

	return core_conditional_sample_edge(sampled_edges, sampled_vertices);
}

template<typename EdgeData>
std::pair<Edge<EdgeData>, double> ZGraphInstance<EdgeData>::conditional_sample_edge(const MyList<Edge<EdgeData>> &edges, const MyList<VertexId> &vertices) {
	if (is_counting_conditional_edge_calls) {
		++ cnt;
	}
	return core_conditional_sample_edge(edges, vertices);
}

template<typename EdgeData>
bool ZGraphInstance<EdgeData>::conditional_close(const SubGraph<EdgeData> &sampled_subgraph, const SubGraph<EdgeData> &subgraph_to_wait) {
	int s_i = SharedMemSys::get_instance()->get_current_socket_id();
	EdgeId min_valid_edge_id = sampled_subgraph.get_min_valid_edge_id();
	const MyList<VertexId>& sampled_vertices = sampled_subgraph.get_vertices();
	const MyList<Edge<EdgeData>>& edges_to_wait = subgraph_to_wait.get_edges();

	for (int i = 0; i < edges_to_wait.list_size; ++ i) {
		const Edge<EdgeData>& edge_to_wait = edges_to_wait.list[i];
		bool found = false;
		for (int j = 0; j < sampled_vertices.list_size; ++ j) {
			VertexId v_i = sampled_vertices.list[j];
			if (graph->outgoing_adj_bitmap[s_i]->get_bit(v_i)) {
				EdgeId begin_idx = graph->outgoing_adj_index[s_i][v_i];
				EdgeId end_idx = graph->outgoing_adj_index[s_i][v_i + 1];
				for (EdgeId idx = begin_idx; idx < end_idx; ++ idx) {
					AdjUnit<EdgeData>& edge = graph->outgoing_adj_list[s_i][idx];
					if (edge.edge_id < min_valid_edge_id) {
						continue;
					}
					if ((edge.src == edge_to_wait.src && edge.dst == edge_to_wait.dst) || (edge.dst == edge_to_wait.src && edge.src == edge_to_wait.dst)) {
						found = true;
						break;
					}
				}
			}
			if (found) {
				break;
			}
		}
		if (! found) {
			return false;
		}
	}

	return true;
}

template<typename EdgeData>
bool ZGraphInstance<EdgeData>::conditional_close(EdgeId min_valid_edge_id, const MyList<VertexId> &sampled_vertices, MyList<Edge<EdgeData>> &edges_to_wait) {
	int s_i = SharedMemSys::get_instance()->get_current_socket_id();
	for (int i = 0; i < edges_to_wait.list_size; ++ i) {
		Edge<EdgeData> & edge_to_wait = edges_to_wait.list[i];
		bool found = false;
		for (int j = 0; j < sampled_vertices.list_size; ++ j) {
			VertexId v_i = sampled_vertices.list[j];
			if (graph->outgoing_adj_bitmap[s_i]->get_bit(v_i)) {
				EdgeId begin_idx = graph->outgoing_adj_index[s_i][v_i];
				EdgeId end_idx = graph->outgoing_adj_index[s_i][v_i + 1];
				for (EdgeId idx = begin_idx; idx < end_idx; ++ idx) {
					AdjUnit<EdgeData>& edge = graph->outgoing_adj_list[s_i][idx];
					if (edge.edge_id < min_valid_edge_id) {
						continue;
					}
					if ((edge.src == edge_to_wait.src && edge.dst == edge_to_wait.dst) || (edge.dst == edge_to_wait.src && edge.src == edge_to_wait.dst)) {
						edge_to_wait.edge_id = edge.edge_id;
						found = true;
						break;
					}
				}
			}
			if (found) {
				break;
			}
		}
		if (! found) {
			return false;
		}
	}
	return true;
}

template<typename EdgeData>
unsigned long long ZGraphInstance<EdgeData>::sample_interger(unsigned long long min_int, unsigned long long max_int) {
	int t_i = SharedMemSys::get_instance()->get_current_thread_id();
	return randoms[t_i].rand_int(min_int, max_int);
}

template<typename EdgeData>
double ZGraphInstance<EdgeData>::get_num_sample_edge_calls() {
	Debug::get_instance()->enter_function("ZGraphInstance::get_num_sample_edge_calls");
	return 1; // TODO
	is_counting_sample_edge_calls = true;
	cnt = 0;
	int num_estimation = 1024;
	int num_success = 0;
	int sum = 0;
	for (int i = 0; i < num_estimation; ++ i) {
		cnt = 0;
		if (this->estimate() > 1e-10) {
			num_success ++;
			sum += cnt;
		}
	}
	is_counting_sample_edge_calls = false;
	Debug::get_instance()->leave_function("ZGraphInstance::get_num_sample_edge_calls");
	return sum * 1. / num_success;
}

template<typename EdgeData>
double ZGraphInstance<EdgeData>::get_num_conditional_sample_edge_calls() {
	Debug::get_instance()->enter_function("ZGraphInstance::get_num_conditional_sample_edge_calls");
	is_counting_conditional_edge_calls = true;
	cnt = 0;
	int num_estimation = 1024;
	int num_success = 0;
	int sum = 0;
	for (int i = 0; i < num_estimation; ++ i) {
		cnt = 0;
		if (this->estimate() > 1e-10) {
			num_success ++;
			sum += cnt;
		}
	}
	is_counting_conditional_edge_calls = false;
	Debug::get_instance()->leave_function("ZGraphInstance::get_num_conditional_sample_edge_calls");
	return sum * 1. / num_success;
}

// AllPredicateZGraphInstance

template<typename EdgeData>
AllPredicateZGraphInstance<EdgeData>::AllPredicateZGraphInstance(Graph<EdgeData> * graph_, std::function<bool(EdgeData)> predicate_): original_graph(graph_), predicate(predicate_){
	csr_constructor = new CSRConstructor<EdgeData>();
	destroyer = new NumaAwareGraphDestroyer<EdgeData>();
	graph_filter = new GraphFilter<EdgeData>(csr_constructor, destroyer, predicate);
	graph_filter->generate_graph(*original_graph, filtered_graph);
	ZGraphInstance<EdgeData>::set_graph(&filtered_graph);
}

template<typename EdgeData> 
AllPredicateZGraphInstance<EdgeData>::~AllPredicateZGraphInstance() {
	graph_filter->destroy_graph(filtered_graph);
	delete graph_filter;
	delete csr_constructor;
	delete destroyer;
}

template<typename EdgeData>
Graph<EdgeData>* AllPredicateZGraphInstance<EdgeData>::get_graph() {
	return original_graph;
}

template<typename EdgeData>
void AllPredicateZGraphInstance<EdgeData>::set_graph(Graph<EdgeData> * graph_) {
	graph_filter->destroy_graph(filtered_graph);
	original_graph = graph_;
	graph_filter->generate_graph(*original_graph, filtered_graph);
	ZGraphInstance<EdgeData>::set_graph(&filtered_graph);
}

template<typename EdgeData> 
double AllPredicateZGraphInstance<EdgeData>::calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation) {
	return ZGraphInstance<EdgeData>::calculate(num_estmators, num_vertices_in_pattern, init_estimation);
}

//// AtLeastOnePredicateZGraphInstance
//
//template<typename EdgeData>
//AtLeastOnePredicateZGraphInstance<EdgeData>::AtLeastOnePredicateZGraphInstance(Graph<EdgeData> * graph_, std::function<bool(EdgeData)> predicate_): ZGraphInstance<EdgeData>(graph_), predicate(predicate_) {
//	assert((! std::is_same<EdgeData, Empty>::value));
//	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
//	Graph<EdgeData> * graph = graph_;
//	matched_list_size = new EdgeId[num_sockets];
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//		matched_list_size[s_i] = 0;
//	}
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//#pragma omp parallel for 
//		for (EdgeId e_i = 0; e_i < graph->num_local_edges[s_i]; ++ e_i) {
//			if (predicate(graph->local_edges[s_i][e_i].edge_data)) {
//				__sync_fetch_and_add(&matched_list_size[s_i], 1);
//			}
//		}
//	}
//	matched_list = new Edge<EdgeData>*[num_sockets];
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//		matched_list[s_i] = (Edge<EdgeData>*) numa_alloc_onnode(sizeof(Edge<EdgeData>) * matched_list_size[s_i], s_i);
//	}
//	EdgeId * curr_pos = new EdgeId[num_sockets];
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//#pragma omp parallel for 
//		for (EdgeId e_i = 0; e_i < graph->num_local_edges[s_i]; ++ e_i) {
//			if (predicate(graph->local_edges[s_i][e_i].edge_data)) {
//				EdgeId pos = __sync_fetch_and_add(&curr_pos[s_i], 1);
//				matched_list[s_i][pos].edge_id = e_i;
//				matched_list[s_i][pos].src = graph->local_edges[s_i][e_i].src;
//				matched_list[s_i][pos].dst = graph->local_edges[s_i][e_i].dst;
//				matched_list[s_i][pos].edge_data = graph->local_edges[s_i][e_i].edge_data;
//			}
//		}
//	}
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//		assert(curr_pos[s_i] == matched_list_size[s_i]);
//	}
//	delete [] curr_pos;
//}
//
//template<typename EdgeData>
//AtLeastOnePredicateZGraphInstance<EdgeData>::~AtLeastOnePredicateZGraphInstance() {
//	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
//	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
//		numa_free(matched_list[s_i], sizeof(Edge<EdgeData>) * matched_list_size[s_i]);
//	}
//	delete [] matched_list;
//	delete [] matched_list_size;
//}
//
//template<typename EdgeData>
//std::pair<Edge<EdgeData>, double> AtLeastOnePredicateZGraphInstance<EdgeData>::sample_edge() {
//	if (this->is_counting_sample_edge_calls) {
//		++ this->cnt;
//	}
//	int t_i = SharedMemSys::get_instance()->get_current_thread_id();
//	int s_i = SharedMemSys::get_instance()->get_socket_id(t_i);
//	Edge<EdgeData> edge;
//	if (matched_list_size[s_i] == 0) {
//		return std::make_pair(edge, (double) -1);
//	}
//	EdgeId sample_edge_idx = this->randoms[t_i].rand_int(0, matched_list_size[s_i] - 1);
//	edge = matched_list[s_i][sample_edge_idx];
//	double probability = 1. / matched_list_size[s_i];
//	return std::make_pair(edge, probability);
//}

// CachedSubPatterns

template<typename EdgeData>
CachedSubPatterns<EdgeData>::CachedSubPatterns(unsigned long long num_estmators_, unsigned long long ** num_saved_sub_patterns_): num_estmators(num_estmators_), num_saved_sub_patterns(num_saved_sub_patterns_) {
	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	int num_threads_per_socket = SharedMemSys::get_instance()->get_num_threads_per_socket();

	// initialize saved_sub_patterns && saved_probability
	saved_sub_patterns = new SubGraph<EdgeData> ** [num_sockets];
	saved_probability = new double ** [num_sockets];
	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		saved_sub_patterns[s_i] = new SubGraph<EdgeData> * [num_threads_per_socket];
		saved_probability[s_i] = new double* [num_threads_per_socket];
		for (int t_i = 0; t_i < num_threads_per_socket; ++ t_i) {
			saved_sub_patterns[s_i][t_i] = (SubGraph<EdgeData>*) numa_alloc_onnode(sizeof(SubGraph<EdgeData>) * num_saved_sub_patterns[s_i][t_i], s_i);
			saved_probability[s_i][t_i] = (double *) numa_alloc_onnode(sizeof(double) * num_saved_sub_patterns[s_i][t_i], s_i);
		}
	}
}

template<typename EdgeData>
CachedSubPatterns<EdgeData>::~CachedSubPatterns() {
	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	int num_threads_per_socket = SharedMemSys::get_instance()->get_num_threads_per_socket();

	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		for (int t_i = 0; t_i < num_threads_per_socket; ++ t_i) {
			numa_free(saved_sub_patterns[s_i][t_i], sizeof(SubGraph<EdgeData>) * num_saved_sub_patterns[s_i][t_i]);
			numa_free(saved_probability[s_i][t_i], sizeof(double) * num_saved_sub_patterns[s_i][t_i]);
		}
		delete [] saved_sub_patterns[s_i];
		delete [] saved_probability[s_i];
	}
	delete [] saved_sub_patterns;
	delete [] saved_probability;

	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		delete [] num_saved_sub_patterns[s_i];
	}
	delete [] num_saved_sub_patterns;
}

template<typename EdgeData>
void CachedSubPatterns<EdgeData>::add(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData> &sub_graph, double probability) {
	saved_sub_patterns[socket_id][socket_offset][estimator_id] = sub_graph;
	saved_probability[socket_id][socket_offset][estimator_id] = probability;
}

//template<typename EdgeData>
//void CachedSubPatterns<EdgeData>::get(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData> &sub_graph, double &probability) {
//	sub_graph = saved_sub_patterns[socket_id][socket_offset][estimator_id];
//	probability = saved_probability[socket_id][socket_offset][estimator_id];
//}

template<typename EdgeData>
void CachedSubPatterns<EdgeData>::get(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData>* &subgraph, double* &probability) {
	subgraph = &saved_sub_patterns[socket_id][socket_offset][estimator_id];
	probability = &saved_probability[socket_id][socket_offset][estimator_id];
}

template<typename EdgeData>
std::pair<SubGraph<EdgeData>*, double*> CachedSubPatterns<EdgeData>::get_all(int socket_id, int socket_offset, unsigned long long &num_patterns) {
	num_patterns = num_saved_sub_patterns[socket_id][socket_offset];
	return std::make_pair(saved_sub_patterns[socket_id][socket_offset], saved_probability[socket_id][socket_offset]);
}

template<typename EdgeData>
unsigned long long CachedSubPatterns<EdgeData>::get_num_patterns(int socket_id, int socket_offset) {
	return num_saved_sub_patterns[socket_id][socket_offset];
}

// ZGraphInstanceWithCache

template<typename EdgeData>
void ZGraphInstanceWithCache<EdgeData>::get_subgraph(bool &is_cached_pattern, SubGraph<EdgeData>* &subgraph, double* &p) {
	if (this->is_counting_sample_edge_calls || this->is_counting_conditional_edge_calls) {
		is_cached_pattern = false;
		subgraph = &fake_subgraph;
		subgraph->edges.list_size = 0;
		subgraph->vertices.list_size = 0;
		subgraph->min_valid_edge_id = 0;
		subgraph->min_valid_vertex_id = 0;
		p = &fake_p;
		return ;
	}
	int socket_id = SharedMemSys::get_instance()->get_current_socket_id();
	int socket_offset = SharedMemSys::get_instance()->get_current_socket_offset();
	unsigned long long estimator_id = curr_estimator_id[socket_id][socket_offset];

	// since exported cache already includes imported cache, we always return exported cache
	if (imported_cached_sub_patterns == nullptr) { // there is no imported cached patterns
		is_cached_pattern = false;
	} else {
		if (estimator_id < imported_cached_sub_patterns->get_num_patterns(socket_id, socket_offset)) { // the imported cache is avaible
			is_cached_pattern = true;
		} else { // there is no avalible imported cache
			is_cached_pattern = false;
		}
	}
	exported_cached_sub_patterns->get(socket_id, socket_offset, estimator_id, subgraph,  p); 
}

template<typename EdgeData> 
ZGraphInstanceWithCache<EdgeData>::ZGraphInstanceWithCache(Graph<EdgeData> * graph_): ZGraphInstance<EdgeData>(graph_) {
	imported_cached_sub_patterns = nullptr;
	exported_cached_sub_patterns = nullptr;

	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	int num_threads_per_socket = SharedMemSys::get_instance()->get_num_threads_per_socket();
	curr_estimator_id = new unsigned long long* [num_sockets];
	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		curr_estimator_id[s_i] = new unsigned long long [num_threads_per_socket];
	}
}

template<typename EdgeData>
ZGraphInstanceWithCache<EdgeData>::~ZGraphInstanceWithCache() {
	if (exported_cached_sub_patterns != nullptr) {
		delete exported_cached_sub_patterns;
	}

	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		delete [] curr_estimator_id[s_i];
	}
	delete [] curr_estimator_id;
}

//template<typename EdgeData>
//void ZGraphInstanceWithCache<EdgeData>::import_cached_sub_pattern(SubGraph<EdgeData> &sub_graph, double &probability) {
//	int socket_id = SharedMemSys::get_instance()->get_current_socket_id();
//	int socket_offset = SharedMemSys::get_instance()->get_current_socket_offset();
//	imported_cached_sub_patterns->get(socket_id, socket_offset, curr_estimator_id[socket_id][socket_offset], sub_graph, probability);
//}
//
//template<typename EdgeData>
//void ZGraphInstanceWithCache<EdgeData>::export_cached_sub_pattern(SubGraph<EdgeData> &sub_graph, double probability) {
//	int socket_id = SharedMemSys::get_instance()->get_current_socket_id();
//	int socket_offset = SharedMemSys::get_instance()->get_current_socket_offset();
//	exported_cached_sub_patterns->add(socket_id, socket_offset, curr_estimator_id[socket_id][socket_offset], sub_graph, probability);
//}

template<typename EdgeData>
CachedSubPatterns<EdgeData> * ZGraphInstanceWithCache<EdgeData>::get_cached_sub_patterns() {
	return exported_cached_sub_patterns;
}

template<typename EdgeData>
void ZGraphInstanceWithCache<EdgeData>::set_cached_sub_patterns(CachedSubPatterns<EdgeData> * cached_sub_patterns_) {
	imported_cached_sub_patterns = cached_sub_patterns_;
}

template<typename EdgeData>
double ZGraphInstanceWithCache<EdgeData>::calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation) {
	Debug::get_instance()->enter_function("ZGraphInstanceWithCache::calulate");
	
	if (init_estimation.get_num_finished_estimators() >= num_estmators) {
		Debug::get_instance()->leave_function("ZGraphInstanceWithCache::calulate");
		return init_estimation.get_currnet_estimation();
	}

	num_estmators -= init_estimation.get_num_finished_estimators();
	int num_threads = SharedMemSys::get_instance()->get_num_threads();
	int num_sockets = SharedMemSys::get_instance()->get_num_sockets();
	int num_threads_per_socket = num_threads / num_sockets;

	//Timer::timer_start("cache_allocate");
	unsigned long long ** num_saved_sub_patterns = new unsigned long long* [num_sockets];
	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		num_saved_sub_patterns[s_i] = new unsigned long long [num_threads_per_socket];
		for (int t_i = 0; t_i < num_threads_per_socket; ++ t_i) {
			num_saved_sub_patterns[s_i][t_i] = num_estmators / num_threads_per_socket;
			if (t_i == num_threads_per_socket - 1) {
				num_saved_sub_patterns[s_i][t_i] += num_estmators % num_threads_per_socket;
			}
		}
	}
	if (exported_cached_sub_patterns != nullptr) {
		delete exported_cached_sub_patterns;
	}
	exported_cached_sub_patterns = new CachedSubPatterns<EdgeData>(num_estmators, num_saved_sub_patterns);
	// merge exported cache with imported cache (could be optimized TODO)
	if (imported_cached_sub_patterns != nullptr) {
#pragma omp parallel 
		{
			int s_i = SharedMemSys::get_instance()->get_current_socket_id();
			int t_i = SharedMemSys::get_instance()->get_current_socket_offset();
			unsigned long long num_imported_patterns;
			unsigned long long num_exported_patterns;
			std::pair<SubGraph<EdgeData>*, double*> imported_cache = imported_cached_sub_patterns->get_all(s_i, t_i, num_imported_patterns);
			std::pair<SubGraph<EdgeData>*, double*> exported_cache = exported_cached_sub_patterns->get_all(s_i, t_i, num_exported_patterns);
			unsigned long long num_patterns_to_copy = std::min(num_imported_patterns, num_exported_patterns);
			Debug::get_instance()->log("s_i = ", s_i, " t_i = ", t_i, " copying ", num_patterns_to_copy, " patterns");
			memcpy(exported_cache.first, imported_cache.first, sizeof(SubGraph<EdgeData>) * num_patterns_to_copy);  
			memcpy(exported_cache.second, imported_cache.second, sizeof(double) * num_patterns_to_copy);
		}
	}
	//Timer::timer_stop("cache_allocate");

	for (int s_i = 0; s_i < num_sockets; ++ s_i) {
		for (int t_i = 0; t_i < num_threads_per_socket; ++ t_i) {
			curr_estimator_id[s_i][t_i] = 0;
		}
	}

	Debug::get_instance()->log("preparing to estimate...");

	double reducer = 0;
#pragma omp parallel reduction(+:reducer)
	{
		double local_reducer = 0;
		int s_i = SharedMemSys::get_instance()->get_current_socket_id();
		int t_i = SharedMemSys::get_instance()->get_current_thread_id();
		int socket_offset = SharedMemSys::get_instance()->get_socket_offset(t_i);
		unsigned local_num_estimators = num_estmators / num_threads_per_socket;
		if (socket_offset == num_threads_per_socket - 1) {
			local_num_estimators += num_estmators % num_threads_per_socket;
		}
		Debug::get_instance()->log("local_num_estimators(t_i = ", t_i, "): ", local_num_estimators);
		for (unsigned long long i = 0; i < local_num_estimators; ++ i) {
			curr_estimator_id[s_i][socket_offset] = i;
			local_reducer += this->estimate();
		}
		local_reducer /= double(num_estmators);
		reducer += local_reducer;
	}

	double global_reducer;
	MPI_Allreduce(&reducer, &global_reducer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double estimation = global_reducer * this->get_scaling_factor(num_vertices_in_pattern);
	init_estimation.update(num_estmators, estimation);

	Debug::get_instance()->leave_function("ZGraphInstanceWithCache::calulate");

	imported_cached_sub_patterns = nullptr;

	return init_estimation.get_currnet_estimation();
}

template class SubGraph<Empty>;
template class SubGraph<VertexId>;

template class ZGraphInstance<Empty>;
template class ZGraphInstance<VertexId>;

//template class AllPredicateZGraphInstance<Empty>;
template class AllPredicateZGraphInstance<VertexId>;

////template class AtLeastOnetPredicateZGraphInstance<Empty>;
//template class AtLeastOnePredicateZGraphInstance<VertexId>;

template class CachedSubPatterns<Empty>;
template class CachedSubPatterns<VertexId>;

template class ZGraphInstanceWithCache<Empty>;
template class ZGraphInstanceWithCache<VertexId>;
