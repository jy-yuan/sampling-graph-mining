#include <stdlib.h>
#include <omp.h>
#include <numa.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>

#include <thread>

#include "Graph.hpp"
#include "NetworkInterface.hpp"
#include "Debug.hpp"
#include "Utilities.hpp"
#include "Random.hpp"
#include "Bitmap.hpp"
#include "SharedMemSys.hpp"
#include "Time.hpp"

#define CHUNK_SIZE (1 << 20)

// UniformGraphPartitioner

UniformGraphPartitioner::UniformGraphPartitioner(VertexId num_vertices_): num_vertices(num_vertices_) {
}

void UniformGraphPartitioner::partition_graph(GraphPartition &partition) {
	Debug::get_instance()->enter_function("UniformGraphPartitioner::partition_graph");

	NetworkInterface * net = NetworkInterface::get_instance();
	partition.partition_id = net->get_partition_id();
	partition.num_partitions = net->get_num_partitions();
	partition.num_sockets = SharedMemSys::get_instance()->get_num_sockets();

	VertexId * partition_offset = new VertexId [partition.num_partitions + 1];

	// uniformly partition the graph
	partition_offset[0] = 0;
	for (int p_i = 1; p_i < partition.num_partitions; ++ p_i) {
		partition_offset[p_i] = num_vertices / partition.num_partitions * p_i;
	}
	partition_offset[partition.num_partitions] = num_vertices;
	partition.owned_num_vertices = partition_offset[partition.partition_id + 1] - partition_offset[partition.partition_id];

	// check consistency
	MPI_Datatype vid_t = get_mpi_data_type<VertexId>();
	VertexId * global_partition_offset = new VertexId[partition.num_partitions + 1];
	MPI_Allreduce(partition_offset, global_partition_offset, partition.num_partitions + 1, vid_t, MPI_MAX, MPI_COMM_WORLD);
	for (int p_i = 0; p_i <= partition.num_partitions; ++ p_i) {
		assert(partition_offset[p_i] == global_partition_offset[p_i]);
	}
	MPI_Allreduce(partition_offset, global_partition_offset, partition.num_partitions + 1, vid_t, MPI_MIN, MPI_COMM_WORLD);
	for (int p_i = 0; p_i <= partition.num_partitions; ++ p_i) {
		assert(partition_offset[p_i] == global_partition_offset[p_i]);
	}
	delete [] global_partition_offset;

	partition.partition_offset = partition_offset;

	// local partition
	VertexId * local_partition_offset = new VertexId[partition.num_sockets + 1];
	local_partition_offset[0] = partition_offset[partition.partition_id];
	for (int s_i = 1; s_i < partition.num_sockets; ++ s_i) {
		local_partition_offset[s_i] = partition_offset[partition.partition_id] + (partition_offset[partition.partition_id + 1] - partition_offset[partition.partition_id]) / partition.num_sockets * s_i;
	}
	local_partition_offset[partition.num_sockets] = partition_offset[partition.partition_id + 1];

	partition.local_partition_offset = local_partition_offset;

	VertexId ** all_local_partition_offset = new VertexId*[partition.num_partitions];
	for (int p_i = 0; p_i < partition.num_partitions; ++ p_i) {
		all_local_partition_offset[p_i] = new VertexId[partition.num_sockets + 1];
		all_local_partition_offset[p_i][0] = partition_offset[p_i];
		VertexId delta = (partition_offset[p_i + 1] - partition_offset[p_i]) / partition.num_sockets;
		for (int s_i = 1; s_i < partition.num_sockets; ++ s_i) {
			all_local_partition_offset[p_i][s_i] = partition_offset[p_i] + delta * s_i;
		}
		all_local_partition_offset[p_i][partition.num_sockets] = partition_offset[p_i + 1];
	}
	partition.all_local_partition_offset = all_local_partition_offset;

	// printing partition result
	Debug::get_instance()->log("partition.partition_id: ", partition.partition_id);
	Debug::get_instance()->log("partition.num_partitions: ", partition.num_partitions);
	Debug::get_instance()->log("partition.num_sockets: ", partition.num_sockets);
	Debug::get_instance()->log("partition.owned_vertices: ", partition.owned_num_vertices);
	{
		std::string log_msg = "global_partition_result: ";
		for (int p_i = 0; p_i < partition.num_partitions; ++ p_i) {
			log_msg += "part " + std::to_string(p_i) + ": [" + std::to_string(partition_offset[p_i]) + ", " + std::to_string(partition_offset[p_i + 1]) + ") = " + std::to_string(partition_offset[p_i + 1] - partition_offset[p_i]) + " ";
		}
		Debug::get_instance()->log(log_msg);
	}
	{
		std::string log_msg = "local_partition_result: ";
		for (int s_i = 0; s_i < partition.num_sockets; ++ s_i) {
			log_msg += "part " + std::to_string(s_i) + ": [" + std::to_string(local_partition_offset[s_i]) + ", " + std::to_string(local_partition_offset[s_i + 1]) + ") = " + std::to_string(local_partition_offset[s_i + 1] - local_partition_offset[s_i]) + " ";
		}
		Debug::get_instance()->log(log_msg);
	}
	{
		std::string log_msg = "all_local_partition_result: ";
		for (int p_i = 0; p_i < partition.num_partitions; ++ p_i) {
			for (int s_i = 0; s_i < partition.num_sockets; ++ s_i) {
				log_msg += "part " + std::to_string(p_i) + "-" + std::to_string(s_i) + ": [" + std::to_string(all_local_partition_offset[p_i][s_i]) + ", " + std::to_string(all_local_partition_offset[p_i][s_i + 1]) + ") = " + std::to_string(all_local_partition_offset[p_i][s_i + 1] - all_local_partition_offset[p_i][s_i]) + " ";
			}
		}
		Debug::get_instance()->log(log_msg);
	}

	Debug::get_instance()->leave_function("UniformGraphPartitioner::partition_graph");
}

// IdenticalMapping

VertexId IdenticalMapping::get_mapping(VertexId vtx) {
	return vtx;
}

VertexId IdenticalMapping::get_inv_mapping(VertexId vtx) {
	return vtx;
}

// RandomMapping

RandomMapping::RandomMapping(unsigned int seed, VertexId num_vertices) { 
	//Debug::get_instance()->enter_function("RandomMapping::RandomMapping");
	srand(seed);
	mapping = new VertexId[num_vertices];
	inv_mapping = new VertexId[num_vertices];
#pragma omp parallel for 
	for (VertexId v_i = 0; v_i < num_vertices; ++ v_i) {
		mapping[v_i] = v_i;
	}
	VertexId pos, tmp;
	for (VertexId v_i = num_vertices - 1; v_i; -- v_i) {
		pos = rand() % (v_i + 1);
		tmp = mapping[v_i];
		mapping[v_i] = mapping[pos];
		mapping[pos] = tmp;
	}
#pragma omp parallel for 
	for (VertexId v_i = 0; v_i < num_vertices; ++ v_i) {
		inv_mapping[mapping[v_i]] = v_i;
	}
	//{
	//	std::string log_msg = "Mapping:\n";
	//	for (VertexId v_i = 0; v_i < num_vertices; ++ v_i) {
	//		log_msg += std::to_string(v_i) + " -> " + std::to_string(mapping[v_i]);
	//		if (v_i != num_vertices - 1) {
	//			log_msg += "\n";
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}
	//Debug::get_instance()->leave_function("RandomMapping::RandomMapping");
}

RandomMapping::~RandomMapping() {
	delete [] mapping;
	delete [] inv_mapping;
}

VertexId RandomMapping::get_mapping(VertexId vtx) {
	return mapping[vtx];
}

VertexId RandomMapping::get_inv_mapping(VertexId vtx) {
	return inv_mapping[vtx];
}

// NumaAwareGraphDestroyer

template<typename EdgeData>
void NumaAwareGraphDestroyer<EdgeData>::destroy_graph(Graph<EdgeData> &graph) {
	Debug::get_instance()->enter_function("NumaAwareGraphDestroyer::destroy_graph");

	Debug::get_instance()->log("to free graph.degree...");
	numa_free(graph.degree, sizeof(EdgeId) * graph.num_vertices);

	Debug::get_instance()->log("to free CSR bitmap...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		delete graph.outgoing_adj_bitmap[s_i]; 
	}
	delete [] graph.outgoing_adj_bitmap;

	Debug::get_instance()->log("to free CSR list...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) { 
		EdgeId end_idx = graph.outgoing_adj_index[s_i][graph.p.local_partition_offset[s_i + 1]];
		if (end_idx) {
			numa_free(graph.outgoing_adj_list[s_i], sizeof(AdjUnit<EdgeData>) * end_idx);
		}
	}
	delete [] graph.outgoing_adj_list;

	Debug::get_instance()->log("to free CSR index...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		numa_free(graph.outgoing_adj_index[s_i] + graph.p.local_partition_offset[s_i], sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1));
	}
	delete [] graph.outgoing_adj_index;

	Debug::get_instance()->log("to free num_local_edges && local_edges...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		if (graph.num_local_edges[s_i]) {
			numa_free(graph.local_edges[s_i], sizeof(EdgeUnit<EdgeData>) * graph.num_local_edges[s_i]);
		}
	}
	delete [] graph.local_edges;
	delete [] graph.num_local_edges;
	Debug::get_instance()->leave_function("NumaAwareGraphDestroyer::destroy_graph");
}

// CSRConstructor

template<typename EdgeData>
void CSRConstructor<EdgeData>::construct_bitmap_assisted_csr(EdgeUnit<EdgeData> ** local_edges, EdgeId * num_local_edges, Graph<EdgeData> &graph) {
	Debug::get_instance()->enter_function("CSRConstructor::construct_bitmap_assisted_csr");

	// shuffling local edges
	Debug::get_instance()->log("shuffling local edges...");
#pragma omp parallel for
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		Random random;
		EdgeId pos;
		EdgeUnit<EdgeData> tmp;
		for (EdgeId e_i = num_local_edges[s_i]; e_i > 1; -- e_i) {
			pos = random.rand_int() % (e_i);
			tmp = local_edges[s_i][e_i - 1];
			local_edges[s_i][e_i - 1] = local_edges[s_i][pos];
			local_edges[s_i][pos] = tmp;
		}
	}

	// alloc numa-aware memory space 
	Debug::get_instance()->log("allocating numa-aware memory space...");
	Bitmap ** outgoing_adj_bitmap = new Bitmap* [graph.p.num_sockets];
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		outgoing_adj_bitmap[s_i] = new Bitmap(graph.num_vertices);
		outgoing_adj_bitmap[s_i]->clear();
	}
	EdgeId ** outgoing_adj_index = new EdgeId* [graph.p.num_sockets];
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		EdgeId * tmp = (EdgeId*) numa_alloc_onnode(sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1), s_i);
		memset(tmp, 0, sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1));
		outgoing_adj_index[s_i] = tmp - graph.p.local_partition_offset[s_i];
	}
	AdjUnit<EdgeData> ** outgoing_adj_list = new AdjUnit<EdgeData>* [graph.p.num_sockets];
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		AdjUnit<EdgeData> * tmp = (AdjUnit<EdgeData>*) numa_alloc_onnode(sizeof(AdjUnit<EdgeData>) * num_local_edges[s_i] * 2, s_i);
		outgoing_adj_list[s_i] = tmp;
	}

	// construct bitmap-assisted CSR graph representation
	Debug::get_instance()->log("getting bitmap-assisted CSR graph representation...");
	Debug::get_instance()->log("getting csr index...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
#pragma omp parallel for
		for (EdgeId e_i = 0; e_i < num_local_edges[s_i]; ++ e_i) {
			VertexId src = local_edges[s_i][e_i].src;
			VertexId dst = local_edges[s_i][e_i].dst;
			assert(src >= graph.p.local_partition_offset[s_i] && src < graph.p.local_partition_offset[s_i + 1]); // there may be a bug here since I am not sure that whether assert() is thread-safe or not
			assert(dst >= graph.p.local_partition_offset[s_i] && dst < graph.p.local_partition_offset[s_i + 1]);
			__sync_fetch_and_add(&outgoing_adj_index[s_i][src + 1], 1);
			__sync_fetch_and_add(&outgoing_adj_index[s_i][dst + 1], 1);
		}
	}
#pragma omp parallel for 
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		for (VertexId v_i = graph.p.local_partition_offset[s_i] + 1; v_i <= graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
			outgoing_adj_index[s_i][v_i] += outgoing_adj_index[s_i][v_i - 1];
		}
		assert(outgoing_adj_index[s_i][graph.p.local_partition_offset[s_i + 1]] == num_local_edges[s_i] * 2);
	}

	EdgeId ** curr_pos = new EdgeId* [graph.p.num_sockets];
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		EdgeId * tmp = (EdgeId*) numa_alloc_onnode(sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1), s_i);
		memset(tmp, 0, sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1));
		curr_pos[s_i] = tmp - graph.p.local_partition_offset[s_i];
#pragma omp parallel for 
		for (VertexId v_i = graph.p.local_partition_offset[s_i]; v_i <= graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
			curr_pos[s_i][v_i] = outgoing_adj_index[s_i][v_i];
		}
	}

	Debug::get_instance()->log("getting csr list...");
#pragma omp parallel for
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		for (EdgeId e_i = 0; e_i < num_local_edges[s_i]; ++ e_i) {
			VertexId src = local_edges[s_i][e_i].src;
			VertexId dst = local_edges[s_i][e_i].dst;

			outgoing_adj_bitmap[s_i]->set_bit(src);
			outgoing_adj_bitmap[s_i]->set_bit(dst);

			EdgeId pos = curr_pos[s_i][src] ++;
			outgoing_adj_list[s_i][pos].edge_id = e_i;
			outgoing_adj_list[s_i][pos].src = src;
			outgoing_adj_list[s_i][pos].dst = dst;
			if (! std::is_same<EdgeData, Empty>::value) {
				outgoing_adj_list[s_i][pos].edge_data = local_edges[s_i][e_i].edge_data;
			}

			pos = curr_pos[s_i][dst] ++;
			outgoing_adj_list[s_i][pos].edge_id = e_i;
			outgoing_adj_list[s_i][pos].src = dst;
			outgoing_adj_list[s_i][pos].dst = src;
			if (! std::is_same<EdgeData, Empty>::value) {
				outgoing_adj_list[s_i][pos].edge_data = local_edges[s_i][e_i].edge_data;
			}
		}
	}

	// check the correctness 
	Debug::get_instance()->log("checking correctness...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
#pragma omp parallel for
		for (VertexId v_i = graph.p.local_partition_offset[s_i]; v_i < graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
			if (! (curr_pos[s_i][v_i] == outgoing_adj_index[s_i][v_i + 1])) {
				Debug::get_instance()->log("pos && outgoing_adj_index not matched! ", "s_i = ", s_i, " v_i = ", v_i);
			}
			assert(curr_pos[s_i][v_i] == outgoing_adj_index[s_i][v_i + 1]);
		}
	}
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
#pragma omp parallel for 
		for (VertexId v_i = graph.p.local_partition_offset[s_i]; v_i < graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
			for (EdgeId e_i = outgoing_adj_index[s_i][v_i]; e_i + 1 < outgoing_adj_index[s_i][v_i + 1]; ++ e_i) {
				if (! (outgoing_adj_list[s_i][e_i].edge_id < outgoing_adj_list[s_i][e_i + 1].edge_id)) {
					Debug::get_instance()->log("edge_id not increasing! ", "s_i = ", s_i, " e_i = ", e_i);
				}
				assert(outgoing_adj_list[s_i][e_i].edge_id < outgoing_adj_list[s_i][e_i + 1].edge_id);
			}
		}
	}

	Debug::get_instance()->log("deallocating curr_pos...");
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		EdgeId * tmp = curr_pos[s_i] + graph.p.local_partition_offset[s_i];
		numa_free(tmp, sizeof(EdgeId) * (graph.p.local_partition_offset[s_i + 1] - graph.p.local_partition_offset[s_i] + 1));
	}
	delete [] curr_pos;

	graph.outgoing_adj_bitmap = outgoing_adj_bitmap;
	graph.outgoing_adj_index = outgoing_adj_index;
	graph.outgoing_adj_list = outgoing_adj_list;

	//{
	//	std::string log_msg = "Bitmap-assisted CSR:\n";
	//	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
	//		log_msg += "socket = " + std::to_string(s_i) + "\n";
	//		log_msg += "bitmap: ";
	//		for (VertexId v_i = graph.p.local_partition_offset[s_i]; v_i < graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
	//			if (outgoing_adj_bitmap[s_i]->get_bit(v_i)) {
	//				log_msg += "1";
	//			} else {
	//				log_msg += "0";
	//			}
	//		}
	//		log_msg += "\n";
	//		log_msg += "index: ";
	//		for (VertexId v_i = graph.p.local_partition_offset[s_i]; v_i <= graph.p.local_partition_offset[s_i + 1]; ++ v_i) {
	//			log_msg += std::to_string(outgoing_adj_index[s_i][v_i]) + " ";
	//		}
	//		log_msg += "\n";
	//		log_msg += "list: ";
	//		for (EdgeId e_i = 0; e_i < outgoing_adj_index[s_i][graph.p.local_partition_offset[s_i + 1]]; ++ e_i) {
	//			log_msg += "(edge_id=" + std::to_string(outgoing_adj_list[s_i][e_i].edge_id) + ", src=" + std::to_string(outgoing_adj_list[s_i][e_i].src) + ", dst=" + std::to_string(outgoing_adj_list[s_i][e_i].dst) + ") ";
	//		}
	//		if (s_i != graph.p.num_sockets) {
	//			log_msg += "\n";
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}

	Debug::get_instance()->leave_function("CSRConstructor::construct_bitmap_assisted_csr");
}

// DistGraphLoader

template<typename EdgeData>
DistGraphLoader<EdgeData>::DistGraphLoader(GraphPartitioner * partitioner_, Mapping * mapping_, CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_): partitioner(partitioner_), mapping(mapping_), csr_constructor(csr_constructor_), destroyer(destroyer_) {
}

template<typename EdgeData>
void DistGraphLoader<EdgeData>::load_graph(std::string graph_path, VertexId num_vertices, Graph<EdgeData> &graph) {
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("load_graph");
	Debug::get_instance()->enter_function("DistGraphLoader::load_graph");

	partitioner->partition_graph(graph.p); 

	EdgeId * num_local_edges = new EdgeId [graph.p.num_sockets];
	EdgeUnit<EdgeData> ** local_edges = new EdgeUnit<EdgeData>* [graph.p.num_sockets];

	graph.edge_data_size = std::is_same<EdgeData, Empty>::value ? 0: sizeof(EdgeData);
	graph.edge_unit_size = graph.edge_data_size + sizeof(VertexId) * 2;
	assert(graph.edge_unit_size == sizeof(EdgeUnit<EdgeData>));
	graph.num_vertices = num_vertices;
	long total_bytes = file_size(graph_path.c_str());
	assert(total_bytes % graph.edge_unit_size == 0);
	graph.num_edges = total_bytes / graph.edge_unit_size;

	EdgeId edges_to_read = graph.num_edges / graph.p.num_partitions;
	if (graph.p.partition_id == graph.p.num_partitions - 1) {
		edges_to_read += graph.num_edges % graph.p.num_partitions;
	}

	long bytes_to_read = edges_to_read * graph.edge_unit_size;
	long read_offset = (graph.num_edges / graph.p.num_partitions) * graph.p.partition_id * graph.edge_unit_size;
	long read_bytes;

	int fin = open(graph_path.c_str(), O_RDONLY);
	EdgeUnit<EdgeData> * read_edge_buffer = new EdgeUnit<EdgeData> [CHUNK_SIZE];

	// calculating degree && num_local_edges
	Debug::get_instance()->log("calculating degree && num_local_edges");
	graph.degree = (EdgeId *) numa_alloc_interleaved(sizeof(EdgeId) * graph.num_vertices);
	assert(graph.degree != NULL);
#pragma omp parallel for 
	for (VertexId v_i = 0; v_i < num_vertices; ++ v_i) {
		graph.degree[v_i] = 0;
	}

	assert(lseek(fin, read_offset, SEEK_SET) == read_offset);
	read_bytes = 0;
	EdgeId ** tmp_num_edges = new EdgeId*[graph.p.num_partitions];
	for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
		tmp_num_edges[p_i] = new EdgeId[graph.p.num_sockets];
		for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
			tmp_num_edges[p_i][s_i] = 0;
		}
	}
	while (read_bytes < bytes_to_read) {
		long curr_read_bytes;
		if (bytes_to_read - read_bytes > long(graph.edge_unit_size) * CHUNK_SIZE) {
			curr_read_bytes = read(fin, read_edge_buffer, graph.edge_unit_size * CHUNK_SIZE);
		} else {
			curr_read_bytes = read(fin, read_edge_buffer, bytes_to_read - read_bytes);
		}
		assert(curr_read_bytes >= 0);
		read_bytes += curr_read_bytes;
		assert(curr_read_bytes % graph.edge_unit_size == 0);
		EdgeId curr_read_edges = curr_read_bytes / graph.edge_unit_size;
		//for (EdgeId e_i = 0; e_i < curr_read_edges; ++ e_i) { 
		//	VertexId src = mapping->get_mapping(read_edge_buffer[e_i].src);
		//	VertexId dst = mapping->get_mapping(read_edge_buffer[e_i].dst);
		//	Debug::get_instance()->log(src, "(", read_edge_buffer[e_i].src, ")", " <--> ", dst, "(", read_edge_buffer[e_i].dst, ")");
		//}
#pragma omp parallel for
		for (EdgeId e_i = 0; e_i < curr_read_edges; ++ e_i) {
			VertexId src = mapping->get_mapping(read_edge_buffer[e_i].src);
			VertexId dst = mapping->get_mapping(read_edge_buffer[e_i].dst);
			assert(src != dst); // self-loops are not allowed in the input graph file
			__sync_fetch_and_add(&graph.degree[src], 1);
			__sync_fetch_and_add(&graph.degree[dst], 1);
			int p_i = graph.p.get_partition_id(src);
			if (p_i == graph.p.get_partition_id(dst)) {
				int s_i = graph.p.get_local_partition_id(p_i, src); 
				if (s_i == graph.p.get_local_partition_id(p_i, dst)) {
					__sync_fetch_and_add(&tmp_num_edges[p_i][s_i], 1);
				}
			}
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, graph.degree, graph.num_vertices, get_mpi_data_type<EdgeId>(), MPI_SUM, MPI_COMM_WORLD);
	EdgeId max_degree = 0;
#pragma omp parallel for reduction(max: max_degree)
	for (VertexId v_i = 0; v_i < graph.num_vertices; ++ v_i) {
		if (graph.degree[v_i] > max_degree) {
			max_degree = graph.degree[v_i];
		}
	}
	graph.max_degree = max_degree;
	for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
		MPI_Allreduce(MPI_IN_PLACE, tmp_num_edges[p_i], graph.p.num_sockets, get_mpi_data_type<EdgeId>(), MPI_SUM, MPI_COMM_WORLD);
		for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
			Debug::get_instance()->log("tmp_num_edges[", p_i, "-", s_i, "]=", tmp_num_edges[p_i][s_i]);
		}
	}
	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
		local_edges[s_i] = (EdgeUnit<EdgeData>*) numa_alloc_onnode(sizeof(EdgeUnit<EdgeData>) * tmp_num_edges[graph.p.partition_id][s_i], s_i);
		num_local_edges[s_i] = 0;
	}
	for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
		delete [] tmp_num_edges[p_i];
	}
	delete [] tmp_num_edges;

	// get local edges
	Debug::get_instance()->log("getting local edges");
	int * buffered_edges = new int[graph.p.num_partitions];
	EdgeUnit<EdgeData> ** send_buffer = new EdgeUnit<EdgeData>* [graph.p.num_partitions];
	for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
		buffered_edges[p_i] = 0;
		send_buffer[p_i] = new EdgeUnit<EdgeData> [CHUNK_SIZE];
	}
	EdgeUnit<EdgeData> * recv_buffer = new EdgeUnit<EdgeData> [CHUNK_SIZE];

	{
		std::thread recv_thread(
				[&]() {
					int finished_count = 0;
					MPI_Status recv_status;
					while (finished_count < graph.p.num_partitions) {
						MPI_Probe(MPI_ANY_SOURCE, ShuffleGraph, MPI_COMM_WORLD, &recv_status);
						int i = recv_status.MPI_SOURCE;
						assert(recv_status.MPI_TAG == ShuffleGraph && i >= 0 && i < graph.p.num_partitions);
						int recv_bytes;
						MPI_Get_count(&recv_status, MPI_CHAR, &recv_bytes);
						if (recv_bytes == 1) {
							finished_count += 1;
							char c;
							MPI_Recv(&c, 1, MPI_CHAR, i, ShuffleGraph, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							continue;
						}
						assert(recv_bytes % graph.edge_unit_size == 0);
						EdgeId recv_edges = recv_bytes / graph.edge_unit_size;
						MPI_Recv(recv_buffer, recv_bytes, MPI_CHAR, i, ShuffleGraph, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#pragma omp parallel for 
						for (EdgeId e_i = 0; e_i < recv_edges; ++ e_i) {
							VertexId src = recv_buffer[e_i].src;
							VertexId dst = recv_buffer[e_i].dst;
							int s_i = graph.p.get_local_partition_id(src);
							if (s_i == graph.p.get_local_partition_id(dst)) { // those edges which crossing different partition should be ignored
								int pos = __sync_fetch_and_add(&num_local_edges[s_i], 1);
								local_edges[s_i][pos] = recv_buffer[e_i];
							}
						}
					}
				}
				);
	
		assert(lseek(fin, read_offset, SEEK_SET) == read_offset);
		read_bytes = 0;
		while (read_bytes < bytes_to_read) {
			long curr_read_bytes;
			if (bytes_to_read - read_bytes > long(graph.edge_unit_size) * CHUNK_SIZE) {
				curr_read_bytes = read(fin, read_edge_buffer, graph.edge_unit_size * CHUNK_SIZE);
			} else {
				curr_read_bytes = read(fin, read_edge_buffer, bytes_to_read - read_bytes);
			}
			assert(curr_read_bytes >= 0);
			read_bytes += curr_read_bytes;
			assert(curr_read_bytes % graph.edge_unit_size == 0);
			EdgeId curr_read_edges = curr_read_bytes / graph.edge_unit_size;
			for (EdgeId e_i = 0; e_i < curr_read_edges; ++ e_i) {
				VertexId src = mapping->get_mapping(read_edge_buffer[e_i].src);
				VertexId dst = mapping->get_mapping(read_edge_buffer[e_i].dst);
				int p_i = graph.p.get_partition_id(src);
				if (p_i == graph.p.get_partition_id(dst)) { // those edges which crossing different partition should be ignored
					int pos = buffered_edges[p_i] ++;
					send_buffer[p_i][pos].src = src;
					send_buffer[p_i][pos].dst = dst;
					if (! std::is_same<EdgeData, Empty>::value) {
						send_buffer[p_i][pos].edge_data = read_edge_buffer[e_i].edge_data;
					}
					if (buffered_edges[p_i] == CHUNK_SIZE) {
						MPI_Send(send_buffer[p_i], graph.edge_unit_size * CHUNK_SIZE, MPI_CHAR, p_i, ShuffleGraph, MPI_COMM_WORLD);
						buffered_edges[p_i] = 0;
					}
				}
			}
		}
		for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
			if (buffered_edges[p_i] == 0) continue;
			//Debug::get_instance()->log("edge_unit_size:", graph.edge_unit_size, " buffered_edges:", buffered_edges[p_i], " send_buffer[0].src:", send_buffer[p_i][0].src, " send_buffer[-1].src:", send_buffer[p_i][buffered_edges[p_i] - 1].src);
			//Debug::get_instance()->log("sending msg to ", p_i, ", size: ", graph.edge_unit_size * buffered_edges[p_i]);
			MPI_Send(send_buffer[p_i], graph.edge_unit_size * buffered_edges[p_i], MPI_CHAR, p_i, ShuffleGraph, MPI_COMM_WORLD);
			buffered_edges[p_i] = 0;
		}
		for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
			char c = 0;
			MPI_Send(&c, 1, MPI_CHAR, p_i, ShuffleGraph, MPI_COMM_WORLD);
		}
	
		recv_thread.join();
	}

	delete [] buffered_edges;
	for (int p_i = 0; p_i < graph.p.num_partitions; ++ p_i) {
		delete [] send_buffer[p_i];
	}
	delete [] send_buffer;
	delete [] recv_buffer;

	csr_constructor->construct_bitmap_assisted_csr(local_edges, num_local_edges, graph);

	//{ // printing local_edges  
	//	std::string log_msg = "local-edges:\n";
	//	for (int s_i = 0; s_i < graph.p.num_sockets; ++ s_i) {
	//		log_msg += "s_i = " + std::to_string(s_i) + ", num_local_edges = " + std::to_string(num_local_edges[s_i]) + "\n";
	//		for (EdgeId e_i = 0; e_i < num_local_edges[s_i]; ++ e_i) {
	//			//VertexId src = mapping->get_inv_mapping(local_edges[s_i][e_i].src);
	//			//VertexId dst = mapping->get_inv_mapping(local_edges[s_i][e_i].dst);
	//			log_msg += "	" + std::to_string(e_i) + ": " + std::to_string(local_edges[s_i][e_i].src) + " <--> " + std::to_string(local_edges[s_i][e_i].dst);
	//			if (e_i != num_local_edges[s_i] - 1 || s_i != graph.p.num_sockets - 1) {
	//				log_msg += "\n";
	//			}
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}

	graph.num_local_edges = num_local_edges;
	graph.local_edges = local_edges;

	close(fin);
	delete [] read_edge_buffer;

	Debug::get_instance()->leave_function("DistGraphLoader::load_graph");

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_stop("load_graph");
}

template<typename EdgeData>
void DistGraphLoader<EdgeData>::destroy_graph(Graph<EdgeData> &graph) {
	destroyer->destroy_graph(graph);
}

// GraphSampler

template<typename EdgeData>
GraphSampler<EdgeData>::GraphSampler(CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_, double sample_rate_): csr_constructor(csr_constructor_), destroyer(destroyer_), sample_rate(sample_rate_) {
}

template<typename EdgeData>
void GraphSampler<EdgeData>::sample_graph(const Graph<EdgeData> &input_graph, Graph<EdgeData> &output_graph) {
	Debug::get_instance()->enter_function("GraphSampler::sample_graph");

	// meta data
	output_graph.p = input_graph.p;
	output_graph.edge_data_size = input_graph.edge_data_size;
	output_graph.edge_unit_size = input_graph.edge_unit_size;
	output_graph.num_vertices = input_graph.num_vertices;

	// local_edges && num_local_edges
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
		Random random;
		const unsigned long long max_rand_int = (unsigned long long) 1 << 30;
		unsigned long long threshold = sample_rate * max_rand_int;

		int t_i = SharedMemSys::get_instance()->get_current_thread_id();
		int s_i = SharedMemSys::get_instance()->get_current_socket_id();
		int num_threads_per_socket = SharedMemSys::get_instance()->get_num_threads_per_socket();
		EdgeId begin_e_i = input_graph.num_local_edges[s_i] / num_threads_per_socket * SharedMemSys::get_instance()->get_socket_offset(t_i);
		EdgeId end_e_i = begin_e_i + input_graph.num_local_edges[s_i] / num_threads_per_socket;
		if (SharedMemSys::get_instance()->get_socket_offset(t_i) == num_threads_per_socket - 1) {
			end_e_i = input_graph.num_local_edges[s_i];
		}
		for (EdgeId e_i = begin_e_i; e_i < end_e_i; ++ e_i) {
			if (random.rand_int(1, max_rand_int) <= threshold) {
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
	//NetworkInterface::get_instance()->pause();

	//{ // printing local_edges  
	//	std::string log_msg = "local-edges:\n";
	//	for (int s_i = 0; s_i < output_graph.p.num_sockets; ++ s_i) {
	//		log_msg += "s_i = " + std::to_string(s_i) + ", num_local_edges = " + std::to_string(output_graph.num_local_edges[s_i]) + "\n";
	//		for (EdgeId e_i = 0; e_i < output_graph.num_local_edges[s_i]; ++ e_i) {
	//			//VertexId src = mapping->get_inv_mapping(local_edges[s_i][e_i].src);
	//			//VertexId dst = mapping->get_inv_mapping(local_edges[s_i][e_i].dst);
	//			log_msg += "	" + std::to_string(e_i) + ": " + std::to_string(output_graph.local_edges[s_i][e_i].src) + " <--> " + std::to_string(output_graph.local_edges[s_i][e_i].dst);
	//			if (e_i != output_graph.num_local_edges[s_i] - 1 || s_i != output_graph.p.num_sockets - 1) {
	//				log_msg += "\n";
	//			}
	//		}
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}

	// bitmap-assisted CSR
	csr_constructor->construct_bitmap_assisted_csr(output_graph.local_edges, output_graph.num_local_edges, output_graph);

	for (int s_i = 0; s_i < input_graph.p.num_sockets; ++ s_i) {
		delete is_reserved_bitmap[s_i];
	}
	delete [] is_reserved_bitmap;

	Debug::get_instance()->log("partition of output-graph: ", output_graph.p);

	Debug::get_instance()->leave_function("GraphSampler::sample_graph");
}

template<typename EdgeData>
void GraphSampler<EdgeData>::destroy_graph(Graph<EdgeData> &graph) { 
	destroyer->destroy_graph(graph);
}

template class DistGraphLoader<Empty>;
template class DistGraphLoader<VertexId>;

template class CSRConstructor<Empty>;
template class CSRConstructor<VertexId>;

template class GraphSampler<Empty>;
template class GraphSampler<VertexId>;

template class NumaAwareGraphDestroyer<Empty>;
template class NumaAwareGraphDestroyer<VertexId>;


