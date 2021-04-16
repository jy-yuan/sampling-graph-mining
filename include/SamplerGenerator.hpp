#ifndef SAMPLER_GENERATOR_HPP
#define SAMPLER_GENERATOR_HPP

#include <iostream>

#include "ZGraphInstance.hpp"
#include "Debug.hpp"

class Pattern { 
	private:
		EdgeId * csr_idx;
		Edge<Empty> * csr_list;
		Edge<Empty> * edges;

		VertexId num_vertices;
		EdgeId num_edges;
	private:
		void construct_csr();
		EdgeId get_idx_by_edge_id(EdgeId edge_id) const;
		bool is_vertex_mapping_correct(VertexId * vertex_mapping, const Pattern &pattern);
		bool core_has_same_edge_ordering(VertexId src, VertexId dst, const Pattern &pattern);
	public:
		Pattern(VertexId num_vertices_, EdgeId num_edges_, Edge<Empty> *edges);
		~Pattern();

		EdgeId* get_csr_idx();
		Edge<Empty>* get_csr_list();
		Edge<Empty>* get_edges();
		VertexId get_num_vertices();
		EdgeId get_num_edges();

		bool has_same_edge_ordering(const Pattern &pattern);
		void setup_edge_id_mapping(EdgeId *mapping);

		bool is_in_edges_list(Edge<Empty> edge);
		bool is_automorphic(VertexId num_vertices_, EdgeId num_edges_, const Edge<Empty> *edges_);
		bool is_automorphic(const Pattern &pattern);
		void setup_vertex_id_mapping(VertexId *mapping);

		friend std::ostream& operator << (std::ostream &os, const Pattern &pattern) {
			os << "Pattern: ";
			os << "num_vertices: " << pattern.num_vertices << ",";
			os << " num_edges: " << pattern.num_edges << ", ";
			for (EdgeId e_i = 0; e_i < pattern.num_edges; ++ e_i) {
				os << "(" << pattern.edges[e_i].edge_id << "):" << pattern.edges[e_i].src << "<->" << pattern.edges[e_i].dst << " ";
			}
			return os;
		}
};

template<typename EdgeData>
class SamplerGenerator: public ZGraphInstance<EdgeData> {
	private:
		VertexId num_vertices;
		EdgeId num_edges;
		VertexId num_patterns;
		VertexId curr_pos;
		VertexId num_unique_patterns;

		Pattern * pattern;
		Pattern ** patterns;
		VertexId * mapping;
		bool * is_used;
		bool * is_duplicated;
		Pattern ** unique_patterns;

		std::function<bool(SubGraph<EdgeData>)> predicate_fun;
		bool used_predicate_fun;

		void generate_mapping(VertexId idx);
		bool canonicality_check(const MyList<Edge<EdgeData>> &edges, const Edge<EdgeData> new_edge);
		void set_predicate_fun(std::function<bool(SubGraph<EdgeData>)> predicate_fun_);
	public:
		SamplerGenerator(Graph<EdgeData> * graph_, Pattern * pattern_);
		~SamplerGenerator();
		double estimate();

		double get_num_sample_edge_calls(); 
		double get_num_conditional_sample_edge_calls(); 
};

#endif
