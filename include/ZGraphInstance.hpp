#ifndef ZGRAPH_INSTANCE_HPP
#define ZGRAPH_INSTANCE_HPP

#include <assert.h>

#include <utility>
#include <vector>
#include <iostream>
#include <functional>

#include "Random.hpp"
#include "Type.hpp"
#include "Graph.hpp"
#include "GraphFilters.hpp"

#define MAX_LIST_SIZE 16

class Estimation {
	private:
		unsigned long long num_finished_estimators;
		double current_estimation;
	public:
		Estimation();
		void reset();
		void update(unsigned long long new_num_estimators, double new_estimation);
		unsigned long long get_num_finished_estimators();
		double get_currnet_estimation();
};

template<typename EdgeData>
struct Edge {
	EdgeId edge_id;
	VertexId src;
	VertexId dst;
	EdgeData edge_data;
};

template<>
struct Edge<Empty> {
	EdgeId edge_id;
	VertexId src;
	union {
		VertexId dst;
		Empty edge_data;
	};
};

template<typename Data>
struct MyList {
	int list_size;
	Data list[MAX_LIST_SIZE];

	MyList() {
		list_size = 0;
	}
	inline void push_back(const Data &data) {
		assert(list_size < MAX_LIST_SIZE);
		list[list_size ++] = data;
	}
};

template<typename Data>
struct VectorList {
	int list_size;
	std::vector<Data> list;

	VectorList() {
		list_size = 0;
	}
	inline void push_back(const Data &data) {
		list.push_back(data);
		list_size += 1;
	}
};

template<typename EdgeData> 
struct SubGraph {
	MyList<Edge<EdgeData>> edges;
	MyList<VertexId> vertices;
	EdgeId min_valid_edge_id;
	VertexId min_valid_vertex_id;

	SubGraph();
	
	bool is_in_vertices_list(VertexId vertex);
	
	void add_edge(const Edge<EdgeData>& edge);
	void add_vertex(VertexId vertex);
	
	const MyList<Edge<EdgeData>>& get_edges() const;
	const MyList<VertexId>& get_vertices() const;
	
	EdgeId get_min_valid_edge_id() const;
	VertexId get_min_valid_vertex_id() const;
	
	SubGraph<EdgeData> operator-(const SubGraph<EdgeData> &sub_graph); 

	friend std::ostream& operator<<(std::ostream& os, const SubGraph &sub_graph) {
		os << "SubGraph: ";
		for (int i = 0; i < sub_graph.edges.list_size; ++ i) {
			os << "(" << "[" << sub_graph.edges.list[i].edge_id << "]: " << sub_graph.edges.list[i].src << "<->" << sub_graph.edges.list[i].dst << ") ";
		}
		return os;
	}
};

template<typename EdgeData>
class ZGraphInstance {
	private:
		Graph<EdgeData> * graph;
	protected:
		Random * randoms; // Random[num_threads]
		bool is_counting_sample_edge_calls;
		bool is_counting_conditional_edge_calls;
		int cnt;

		double get_scaling_factor(VertexId num_vertices_in_pattern);
		std::pair<Edge<EdgeData>, double> core_conditional_sample_edge(const MyList<Edge<EdgeData>> &sampled_edges, const MyList<VertexId> &sampled_vertices);
	public:
		ZGraphInstance(); // set_graph must be called if this instance is initiliazed with no init graph object
		ZGraphInstance(Graph<EdgeData> * graph_);
		virtual ~ZGraphInstance();
		virtual Graph<EdgeData>* get_graph();
		virtual void set_graph(Graph<EdgeData> * graph_);
		virtual double calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation);

		virtual std::pair<VertexId, double> sample_vertex(); // APIs mentioned in ASAP paper // TODO could be protected
		virtual std::pair<Edge<EdgeData>, double> sample_edge();
		virtual std::pair<VertexId, double> conditional_sample_vertex(const SubGraph<EdgeData> &sub_graph);
		virtual std::pair<Edge<EdgeData>, double> conditional_sample_edge(const SubGraph<EdgeData> &sub_graph); // check whether the edges following those in the first subgraph are able to form the second subgraph
		virtual std::pair<Edge<EdgeData>, double> conditional_sample_edge(const MyList<Edge<EdgeData>> &edges, const MyList<VertexId> &vertices);
		virtual bool conditional_close(const SubGraph<EdgeData> &sampled_subgraph, const SubGraph<EdgeData> &subgraph_to_wait);
		virtual bool conditional_close(EdgeId min_valid_edge_id, const MyList<VertexId> &sampled_vertices, MyList<Edge<EdgeData>> &edges_to_wait);
		virtual unsigned long long sample_interger(unsigned long long min_int, unsigned long long max_int);

		virtual double get_num_sample_edge_calls(); // not thread-safe
		virtual double get_num_conditional_sample_edge_calls(); // not thread-safe

		virtual double estimate() = 0; // describe the behaviors of user-defined estimator, a negative number indicates that the sampling fails
};

template<typename EdgeData>
class AllPredicateZGraphInstance: public ZGraphInstance<EdgeData>{
	private:
		Graph<EdgeData> * original_graph;
		Graph<EdgeData> filtered_graph;
		GraphFilter<EdgeData> * graph_filter;
		CSRConstructor<EdgeData> * csr_constructor;
		GraphDestroyer<EdgeData> * destroyer;
		std::function<bool(EdgeData)> predicate;
	public:
		AllPredicateZGraphInstance(Graph<EdgeData> * graph_, std::function<bool(EdgeData)> predicate_); 
		~AllPredicateZGraphInstance();
		Graph<EdgeData> * get_graph();
		void set_graph(Graph<EdgeData> * graph_);
		double calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation);
};

//template<typename EdgeData>
//class AtLeastOnePredicateZGraphInstance: public ZGraphInstance<EdgeData> { // should be rewrite: sample_edge(), estimate(), constructor, destructor
//	private:
//		Edge<EdgeData> ** matched_list; // Edge<EdgeData>[num_sockets][num_matched_edges], numa-aware
//		EdgeId * matched_list_size; // EdgeId[num_sockets]
//	protected:
//		std::function<bool(EdgeData)> predicate;
//	public:
//		AtLeastOnePredicateZGraphInstance(Graph<EdgeData> * graph_, std::function<bool(EdgeData)> predicate_);
//		~AtLeastOnePredicateZGraphInstance();
//		std::pair<Edge<EdgeData>, double> sample_edge(); // the semantic of this API has been changed to sample an edge from the matched list 
//};

template<typename EdgeData>
class CachedSubPatterns { 
	private:
		unsigned long long num_estmators;
		unsigned long long ** num_saved_sub_patterns; // unsigned long long[num_sockets][num_threads_per_socket]
		SubGraph<EdgeData> *** saved_sub_patterns; // SubGraph<EdgeData>[num_sockets][num_threads_per_socket][...], numa-aware
		double *** saved_probability; // double[num_sockets][num_threads_per_socket][...], numa-aware
	public:
		CachedSubPatterns(unsigned long long num_estmators_, unsigned long long ** num_saved_sub_patterns_);
		~CachedSubPatterns();
		void add(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData> &sub_graph, double probability);
		//void get(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData> &sub_graph, double &probability); 
		void get(int socket_id, int socket_offset, unsigned long long estimator_id, SubGraph<EdgeData>* &subgraph, double* &probability);
		std::pair<SubGraph<EdgeData>*, double*> get_all(int socket_id, int socket_offset, unsigned long long &num_patterns);
		unsigned long long get_num_patterns(int socket_id, int socket_offset);
};

// issues about ZGraphInstanceWithCache:
// 1. number_of_estimators could be different
// 2. avoid calling copy contructor of Sugraph class
// 3. there is two options:
// 	(1) merge imported cache with exported cache at the beginning of the calculate function (trying this one)
// 	(2) don't care about the correctness of exported cache

template<typename EdgeData>
class ZGraphInstanceWithCache: public ZGraphInstance<EdgeData> { 
	private:
		CachedSubPatterns<EdgeData> * imported_cached_sub_patterns;
		CachedSubPatterns<EdgeData> * exported_cached_sub_patterns;
		unsigned long long ** curr_estimator_id; // unsigned long long[num_sockets][num_threads_per_socket] 
		SubGraph<EdgeData> fake_subgraph;
		double fake_p;
	protected:
		void get_subgraph(bool &is_cached_pattern,SubGraph<EdgeData>* &subgraph, double* &p); // if is_cached_pattern = false, probability is not defined
		//void import_cached_sub_pattern(SubGraph<EdgeData> &sub_graph, double &probability);
		//void export_cached_sub_pattern(SubGraph<EdgeData> &sub_graph, double probability);
	public:
		ZGraphInstanceWithCache(Graph<EdgeData> * graph_);
		~ZGraphInstanceWithCache();
		CachedSubPatterns<EdgeData> * get_cached_sub_patterns();
		void set_cached_sub_patterns(CachedSubPatterns<EdgeData> * cached_sub_patterns_);
		double calculate(unsigned long long num_estmators, VertexId num_vertices_in_pattern, Estimation &init_estimation);
};

#endif
