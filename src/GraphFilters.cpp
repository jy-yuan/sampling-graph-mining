#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <numa.h>

#include "GraphFilters.hpp"
#include "SharedMemSys.hpp"
#include "Debug.hpp"

// GraphFilter

template<typename EdgeData> 
GraphFilter<EdgeData>::GraphFilter(CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_, std::function<bool(EdgeData)> predicate_): csr_constructor(csr_constructor_), destroyer(destroyer_), predicate(predicate_){
}

template<typename EdgeData>
void GraphFilter<EdgeData>::generate_graph(const Graph<EdgeData> &input_graph, Graph<EdgeData> &output_graph) {
	Debug::get_instance()->enter_function("GraphFilter::generate_graph");

	// meta data
	output_graph.p = input_graph.p;
	output_graph.edge_data_size = input_graph.edge_data_size;
	output_graph.edge_unit_size = input_graph.edge_unit_size;
	output_graph.num_vertices = input_graph.num_vertices;
	
	// local edges && num_local_edges
	Bitmap ** is_reserved_bitmap = new Bitmap*[input_graph.p.num_sockets];
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		is_reserved_bitmap[s_i] = new Bitmap(input_graph.num_local_edges[s_i]);
		is_reserved_bitmap[s_i]->clear();
	}

	output_graph.num_local_edges = new EdgeId[input_graph.p.num_sockets];
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		output_graph.num_local_edges[s_i] = 0;
	}

#pragma omp parallel
	{
		int t_i = SharedMemSys::get_instance()->get_current_thread_id();
		int s_i = SharedMemSys::get_instance()->get_current_socket_id();
		int num_threads_per_socket = SharedMemSys::get_instance()->get_num_threads_per_socket();
		EdgeId begin_e_i = input_graph.num_local_edges[s_i] / num_threads_per_socket * SharedMemSys::get_instance()->get_socket_offset(t_i);
		EdgeId end_e_i = begin_e_i + input_graph.num_local_edges[s_i] / num_threads_per_socket;
		if (SharedMemSys::get_instance()->get_socket_offset(t_i) == num_threads_per_socket - 1) {
			end_e_i = input_graph.num_local_edges[s_i];
		}
		for (EdgeId e_i = begin_e_i; e_i < end_e_i; ++ e_i) {
			if (predicate(input_graph.local_edges[s_i][e_i].edge_data)) {
				is_reserved_bitmap[s_i]->set_bit(e_i);
				__sync_fetch_and_add(&output_graph.num_local_edges[s_i], 1);
			}
		}
	}

	output_graph.local_edges = new EdgeUnit<EdgeData>* [input_graph.p.num_sockets];
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		output_graph.local_edges[s_i] = (EdgeUnit<EdgeData>*) numa_alloc_onnode(sizeof(EdgeUnit<EdgeData>) * output_graph.num_local_edges[s_i], s_i);
	}

	EdgeId * curr_pos = new EdgeId[input_graph.p.num_sockets];
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		curr_pos[s_i] = 0;
	}
#pragma omp parallel
	{
		SharedMemSys * share_mem_sys = SharedMemSys::get_instance();
		int t_i = share_mem_sys->get_current_thread_id();
		int s_i = share_mem_sys->get_current_socket_id();
		int num_threads_per_socket = share_mem_sys->get_num_threads_per_socket();

		EdgeId begin_e_i = input_graph.num_local_edges[s_i] / num_threads_per_socket * share_mem_sys->get_socket_offset(t_i);
		EdgeId end_e_i = begin_e_i + input_graph.num_local_edges[s_i] / num_threads_per_socket;
		if (share_mem_sys->get_socket_offset(t_i) == num_threads_per_socket - 1) {
			end_e_i = input_graph.num_local_edges[s_i];
		}
		for (EdgeId e_i = begin_e_i; e_i < end_e_i; ++ e_i) {
			if (is_reserved_bitmap[s_i]->get_bit(e_i)) {
				EdgeId pos = __sync_fetch_and_add(&curr_pos[s_i], 1);
				output_graph.local_edges[s_i][pos] = input_graph.local_edges[s_i][e_i];
			}
		}
	}
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		assert(curr_pos[s_i] == output_graph.num_local_edges[s_i]);
	}
	delete [] curr_pos;

	// estimate output_graph.num_edges
	EdgeId input_num_edges_after_partition = 0;
	EdgeId output_num_edges_after_partition = 0;
	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		input_num_edges_after_partition += input_graph.num_local_edges[s_i];
		output_num_edges_after_partition += output_graph.num_local_edges[s_i];
	}
	MPI_Allreduce(MPI_IN_PLACE, &input_num_edges_after_partition, 1, get_mpi_data_type<EdgeId>(), MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &output_num_edges_after_partition, 1, get_mpi_data_type<EdgeId>(), MPI_SUM, MPI_COMM_WORLD);
	double sample_rate_after_partition = (double) input_num_edges_after_partition / input_graph.num_edges;
	Debug::get_instance()->log("1 / sample_rate_after_partition = ", 1. / sample_rate_after_partition);
	output_graph.num_edges = (EdgeId)(1. * output_num_edges_after_partition / sample_rate_after_partition);
	Debug::get_instance()->log("num_edges of input_graph: ", input_graph.num_edges, " num_edges of output_graph: ", output_graph.num_edges);

	// degree
	output_graph.degree = (EdgeId*) numa_alloc_interleaved(sizeof(EdgeId) * output_graph.num_vertices);
	memset(output_graph.degree, 0, sizeof(EdgeId) * output_graph.num_vertices);
#pragma omp parallel
	{
		SharedMemSys * share_mem_sys = SharedMemSys::get_instance();
		int t_i = share_mem_sys->get_current_thread_id();
		int s_i = share_mem_sys->get_current_socket_id();
		int num_threads_per_socket = share_mem_sys->get_num_threads_per_socket();
		int socket_offset = share_mem_sys->get_socket_offset(t_i);

		EdgeId begin_e_i = output_graph.num_local_edges[s_i] / num_threads_per_socket * socket_offset;
		EdgeId end_e_i = begin_e_i + output_graph.num_local_edges[s_i] / num_threads_per_socket;
		if (socket_offset == num_threads_per_socket - 1) {
			end_e_i = output_graph.num_local_edges[s_i];
		}
		for (EdgeId e_i = begin_e_i; e_i < end_e_i; ++ e_i) {
			VertexId src = output_graph.local_edges[s_i][e_i].src;
			VertexId dst = output_graph.local_edges[s_i][e_i].dst;
			__sync_fetch_and_add(&output_graph.degree[src], 1);
			__sync_fetch_and_add(&output_graph.degree[dst], 1);
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, output_graph.degree, output_graph.num_vertices, get_mpi_data_type<EdgeId>(), MPI_SUM, MPI_COMM_WORLD);

	EdgeId max_degree = 0;
#pragma omp parallel for reduction(max: max_degree)
	for (VertexId v_i = 0; v_i < output_graph.num_vertices; ++ v_i) {
		if (output_graph.degree[v_i] > max_degree) {
			max_degree = output_graph.degree[v_i];
		}
	}
	output_graph.max_degree = max_degree / sample_rate_after_partition;
	Debug::get_instance()->log("max-degree of input_graph: ", input_graph.max_degree, " max-degree of output_graph: ", output_graph.max_degree);

	// bitmap-assisted CSR
	csr_constructor->construct_bitmap_assisted_csr(output_graph.local_edges, output_graph.num_local_edges, output_graph);

	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		delete is_reserved_bitmap[s_i];
	}
	delete [] is_reserved_bitmap;

	Debug::get_instance()->log("partition of output-graph: ", output_graph.p);

	//{ // printing input_graph local_edges  
	//	std::string log_msg = "input-local-edges:\n";
	//	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
	//		log_msg += "s_i = " + std::to_string(s_i) + ", num_local_edges = " + std::to_string(input_graph.num_local_edges[s_i]) + "\n";
	//		for (EdgeId e_i = 0; e_i < input_graph.num_local_edges[s_i]; ++ e_i) {
	//			log_msg += "	" + std::to_string(e_i) + ": " + std::to_string(input_graph.local_edges[s_i][e_i].src) + " <--> " + std::to_string(input_graph.local_edges[s_i][e_i].dst) + " " + std::to_string(input_graph.local_edges[s_i][e_i].edge_data);
	//			if (e_i != input_graph.num_local_edges[s_i] - 1 || s_i != input_graph.p.num_sockets - 1) {
	//				log_msg += "\n";
	//			}
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}

	//{ // printing output_graph local_edges  
	//	std::string log_msg = "output-local-edges:\n";
	//	for (int s_i = 0; s_i < output_graph.p.num_sockets; ++ s_i) {
	//		log_msg += "s_i = " + std::to_string(s_i) + ", num_local_edges = " + std::to_string(output_graph.num_local_edges[s_i]) + "\n";
	//		for (EdgeId e_i = 0; e_i < output_graph.num_local_edges[s_i]; ++ e_i) {
	//			log_msg += "	" + std::to_string(e_i) + ": " + std::to_string(output_graph.local_edges[s_i][e_i].src) + " <--> " + std::to_string(output_graph.local_edges[s_i][e_i].dst) + " " + std::to_string(output_graph.local_edges[s_i][e_i].edge_data);
	//			if (e_i != output_graph.num_local_edges[s_i] - 1 || s_i != output_graph.p.num_sockets - 1) {
	//				log_msg += "\n";
	//			}
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}

	Debug::get_instance()->leave_function("GraphFilter::generate_graph");
}

template<typename EdgeData>
void GraphFilter<EdgeData>::destroy_graph(Graph<EdgeData> &graph) {
	destroyer->destroy_graph(graph);
}

template class GraphFilter<VertexId>;
