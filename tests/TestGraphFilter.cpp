#include <stdio.h>

#include <string>

#include "SharedMemSys.hpp"
#include "Debug.hpp"
#include "Graph.hpp"
#include "Utilities.hpp"
#include "GraphFilters.hpp"

int main(int argc, char ** argv) {
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");

	SharedMemSys::init_shared_mem_sys();
	//std::string graph_path = "../datasets/simple-30/dataset/simple-30.weighted";
	std::string graph_path = "../datasets/simple-10/dataset/simple-10.weighted";
	//std::string graph_path = "../datasets/soc-livejournal/dataset/live-journal.weighted";
	VertexId num_vertices;
	load_graph_meta_data(graph_path + ".csv", num_vertices);
	Debug::get_instance()->log("num_vertices = ", num_vertices);

	UniformGraphPartitioner * partitioner = new UniformGraphPartitioner(num_vertices);
	unsigned int seed = 1234;
	RandomMapping * mapping = new RandomMapping(seed, num_vertices);
	CSRConstructor<VertexId> * csr_constructor = new CSRConstructor<VertexId>();
	NumaAwareGraphDestroyer<VertexId> * graph_destroyer = new NumaAwareGraphDestroyer<VertexId>();
	DistGraphLoader<VertexId> * graph_loader = new DistGraphLoader<VertexId>(partitioner, mapping, csr_constructor, graph_destroyer);
	Graph<VertexId> graph;

	graph_loader->load_graph(graph_path + ".biedgelist", num_vertices, graph);

	auto predicate = [&](VertexId edge_data) {
		return edge_data >= 50;
	};
	GraphFilter<VertexId> * graph_filter = new GraphFilter<VertexId>(csr_constructor, graph_destroyer, predicate);
	Graph<VertexId> filtered_graph;
	graph_filter->generate_graph(graph, filtered_graph);
	graph_filter->destroy_graph(filtered_graph);

	graph_loader->destroy_graph(graph);
	delete partitioner;
	delete mapping;
	delete graph_loader;
	delete graph_destroyer;
	delete csr_constructor;
	delete graph_filter;

	Debug::get_instance()->enter_function("main");
	NetworkInterface::finalize_instance();
	return 0;
}
