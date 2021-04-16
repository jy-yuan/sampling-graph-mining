#ifndef GRAPH_FILTERS_HPP
#define GRAPH_FILTERS_HPP

#include <functional>

#include "Graph.hpp"

template<typename EdgeData>
class GraphFilter { // TODO could be merged with GraphSampler class
	private:
		CSRConstructor<EdgeData> * csr_constructor;
		GraphDestroyer<EdgeData> * destroyer;
		std::function<bool(EdgeData)> predicate;
	public:
		GraphFilter(CSRConstructor<EdgeData> * csr_constructor_, GraphDestroyer<EdgeData> * destroyer_, std::function<bool(EdgeData)> predicate_);
		void generate_graph(const Graph<EdgeData> &input_graph, Graph<EdgeData> &output_graph);
		void destroy_graph(Graph<EdgeData> &graph);
};

#endif
