#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <assert.h>

#include <string>

#include "Type.hpp"
#include "Bitmap.hpp"
#include "Debug.hpp"

struct GraphPartition {
	int partition_id;
	int num_partitions;
	int num_sockets;

	VertexId owned_num_vertices;

	VertexId * partition_offset; // VertexId[num_partitions + 1]
	VertexId * local_partition_offset; // VertexId[num_sockets + 1]
	VertexId ** all_local_partition_offset; // VertexId[num_partitions][num_sockets + 1]

	GraphPartition() {
		partition_offset = nullptr;
		local_partition_offset = nullptr;
		all_local_partition_offset = nullptr;
	}
	~GraphPartition() {
		if (partition_offset != nullptr) {
			delete [] partition_offset;
		}
		if (local_partition_offset != nullptr) {
			delete [] local_partition_offset;
		}
		if (all_local_partition_offset != nullptr) {
			for (int p_i = 0; p_i < num_partitions; ++ p_i) {
				delete [] all_local_partition_offset[p_i];
			}
			delete [] all_local_partition_offset;
		}
	}
	GraphPartition(const GraphPartition &partition) {
		//Debug::get_instance()->log("****** GraphPartition:copy constructor called");
		partition_id = partition.partition_id;
		num_partitions = partition.num_partitions;
		num_sockets = partition.num_sockets;

		owned_num_vertices = partition.owned_num_vertices;

		partition_offset = new VertexId[num_partitions + 1];
		for (int p_i = 0; p_i <= num_partitions; ++ p_i) {
			partition_offset[p_i] = partition.partition_offset[p_i];
		}
		local_partition_offset = new VertexId[num_sockets + 1];
		for (int s_i = 0; s_i <= num_sockets; ++ s_i) {
			local_partition_offset[s_i] = partition.local_partition_offset[s_i];
		}
		all_local_partition_offset = new VertexId*[num_partitions];
		for (int p_i = 0; p_i < num_partitions; ++ p_i) {
			all_local_partition_offset[p_i] = new VertexId[num_sockets + 1];
			for (int s_i = 0; s_i <= num_sockets; ++ s_i) {
				all_local_partition_offset[p_i][s_i] = partition.all_local_partition_offset[p_i][s_i];
			}
		}
	}
	GraphPartition& operator=(const GraphPartition &partition) {
		if (this != &partition) {
			//Debug::get_instance()->log("****** GraphPartition:copy assignment operator called");
			partition_id = partition.partition_id;
			num_partitions = partition.num_partitions;
			num_sockets = partition.num_sockets;
			
			owned_num_vertices = partition.owned_num_vertices;
			
			partition_offset = new VertexId[num_partitions + 1];
			for (int p_i = 0; p_i <= num_partitions; ++ p_i) {
				partition_offset[p_i] = partition.partition_offset[p_i];
			}
			local_partition_offset = new VertexId[num_sockets + 1];
			for (int s_i = 0; s_i <= num_sockets; ++ s_i) {
				local_partition_offset[s_i] = partition.local_partition_offset[s_i];
			}
			all_local_partition_offset = new VertexId*[num_partitions];
			for (int p_i = 0; p_i < num_partitions; ++ p_i) {
				all_local_partition_offset[p_i] = new VertexId[num_sockets + 1];
				for (int s_i = 0; s_i <= num_sockets; ++ s_i) {
					all_local_partition_offset[p_i][s_i] = partition.all_local_partition_offset[p_i][s_i];
				}
			}
		}
		return *this;
	}

	inline int get_partition_id(VertexId v_i) {
		for (int p_i = 0; p_i < num_partitions; ++ p_i) {
			if (v_i >= partition_offset[p_i] && v_i < partition_offset[p_i + 1]) {
				return p_i;
			}
		}
		assert(false);
	}
	inline int get_local_partition_id(VertexId v_i) {
		for (int s_i = 0; s_i < num_sockets; ++ s_i) {
			if (v_i >= local_partition_offset[s_i] && v_i < local_partition_offset[s_i + 1]) {
				return s_i;
			}
		}
		assert(false);
	}
	inline int get_local_partition_id(int p_i, VertexId v_i) {
		for (int s_i = 0; s_i < num_sockets; ++ s_i) {
			if (v_i >= all_local_partition_offset[p_i][s_i] && v_i < all_local_partition_offset[p_i][s_i + 1]) {
				return s_i;
			}
		}
		assert(false);
	}

	friend std::ostream& operator<<(std::ostream& os, const GraphPartition &partition) {
		os << "Graph Partition: ";
		os << "partition_offset: ";
		for (int p_i = 0; p_i <= partition.num_partitions; ++ p_i) {
			os << std::to_string(partition.partition_offset[p_i]) << " ";
		}
		os << "local_partition_offset: ";
		for (int s_i = 0; s_i <= partition.num_sockets; ++ s_i) {
			os << std::to_string(partition.local_partition_offset[s_i]) << " ";
		}
		return os;
	}
};

class GraphPartitioner {
	public:
		virtual ~GraphPartitioner() {}
		virtual void partition_graph(GraphPartition &partition) = 0;
};

class UniformGraphPartitioner: public GraphPartitioner {
	private:
		VertexId num_vertices;
	public: 
		~UniformGraphPartitioner() {}
		UniformGraphPartitioner(VertexId num_vertices_); 
		void partition_graph(GraphPartition &partition); 
};

template<typename EdgeData>
struct Graph {
	GraphPartition p;

	size_t edge_data_size; 
	size_t edge_unit_size;

	VertexId num_vertices; // graph meta-data
	EdgeId num_edges;
	EdgeId * degree; // EdgeId[num_vertices]
	EdgeId max_degree; 

	EdgeId * num_local_edges; // EdgeId[num_sockets]
	EdgeUnit<EdgeData> ** local_edges; // EdgeUnit<EdgeData>[num_sockets][...]; numa-aware

	Bitmap ** outgoing_adj_bitmap; // Bitmap*[num_sockets]
	EdgeId ** outgoing_adj_index; // EdgeId[num_sockets][local_partition_offset[s_i] ... local_partition_offset[s_i+1]]; numa-aware
	AdjUnit<EdgeData> ** outgoing_adj_list; // AdjUnit<EdgeData>[num_sockets][...]; numa-aware
};

class Mapping { // every vertex u in original graph will be mapped to mapping[u] in the loaded graph, all vertex numbers in this program except those loader directed from the graph file will be the mapped numbers
	public:
		virtual ~Mapping() {}
		virtual VertexId get_mapping(VertexId vtx) = 0;
		virtual VertexId get_inv_mapping(VertexId vtx) = 0;
};

class IdenticalMapping: public Mapping {
	public:
		~IdenticalMapping() {}
		VertexId get_mapping(VertexId vtx);
		VertexId get_inv_mapping(VertexId vtx);
};

class RandomMapping: public Mapping {
	private:
		VertexId * mapping; // VertexId[num_vertices]
		VertexId * inv_mapping; // VertexId[num_vertices]
	public:
		RandomMapping(unsigned int seed, VertexId num_vertices); // not thread-safe
		~RandomMapping();
		VertexId get_mapping(VertexId vtx);
		VertexId get_inv_mapping(VertexId vtx);
};

template<typename EdgeData>
class GraphDestroyer {
	public:
		virtual ~GraphDestroyer() {}
		virtual void destroy_graph(Graph<EdgeData> &graph) = 0;
};

template<typename EdgeData>
class NumaAwareGraphDestroyer: public GraphDestroyer<EdgeData> {
	public:
		void destroy_graph(Graph<EdgeData> &graph); 
};

template<typename EdgeData>
class CSRConstructor {
	public:
		void construct_bitmap_assisted_csr(EdgeUnit<EdgeData> ** local_edges, EdgeId * num_local_edges, Graph<EdgeData> &graph);
};

template<typename EdgeData>
class GraphLoader {
	public:
		virtual ~GraphLoader() {}
		virtual void load_graph(std::string graph_path, VertexId num_vertices, Graph<EdgeData> &graph) = 0;
		virtual void destroy_graph(Graph<EdgeData> &graph) = 0;
};

template<typename EdgeData>
class DistGraphLoader: public GraphLoader<EdgeData> {
	private:
		GraphPartitioner * partitioner;
		Mapping * mapping;
		CSRConstructor<EdgeData> * csr_constructor;
		GraphDestroyer<EdgeData> * destroyer;
	public:
		DistGraphLoader(GraphPartitioner * partitioner_, Mapping * mapping_, CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_);
		~DistGraphLoader() {}
		void load_graph(std::string graph_path, VertexId num_vertices, Graph<EdgeData> &graph);
		void destroy_graph(Graph<EdgeData> &graph);
};

template<typename EdgeData>
class GraphSampler {
	private:
		CSRConstructor<EdgeData> * csr_constructor;
		GraphDestroyer<EdgeData> * destroyer;
		double sample_rate;
	public:
		GraphSampler(CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_, double sample_rate_ = 0.05); 
		void sample_graph(const Graph<EdgeData> &input_graph, Graph<EdgeData> &output_graph);
		void destroy_graph(Graph<EdgeData> &graph);
};

#endif
