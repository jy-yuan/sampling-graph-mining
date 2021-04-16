#include <stdio.h>

#include <string>

#include "SharedMemSys.hpp"
#include "Debug.hpp"
#include "Graph.hpp"
#include "Utilities.hpp"

int main(int argc, char ** argv) {
	//NetworkInterface * net_instance = NetworkInterface::get_instance(&argc, &argv);
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");

	//SharedMemSys * shared_mem_sys = SharedMemSys::get_instance();
	SharedMemSys::init_shared_mem_sys();

	std::string graph_path = "../datasets/soc-livejournal/dataset/live-journal";
	//std::string graph_path = "../datasets/wiki-vote/dataset/wiki-vote";
	//std::string graph_path = "../datasets/simple-10000/dataset/simple-10000";
	//std::string graph_path = "../datasets/simple-10/dataset/simple-10";
	//std::string graph_path = "../datasets/simple-30/dataset/simple-30";
	VertexId num_vertices;
	load_graph_meta_data(graph_path + ".csv", num_vertices);
	Debug::get_instance()->log("num_vertices = ", num_vertices);

	UniformGraphPartitioner * partitioner = new UniformGraphPartitioner(num_vertices);
	unsigned int seed = 1234;
	RandomMapping * mapping = new RandomMapping(seed, num_vertices);
	CSRConstructor<Empty> * csr_constructor = new CSRConstructor<Empty>();
	NumaAwareGraphDestroyer<Empty> * graph_destroyer = new NumaAwareGraphDestroyer<Empty>();
	DistGraphLoader<Empty> * graph_loader = new DistGraphLoader<Empty>(partitioner, mapping, csr_constructor, graph_destroyer);
	Graph<Empty> graph;

	graph_loader->load_graph(graph_path + ".biedgelist", num_vertices, graph);

	GraphSampler<Empty> * sampler = new GraphSampler<Empty>(csr_constructor, graph_destroyer);
	Graph<Empty> sampled_graph;
	sampler->sample_graph(graph, sampled_graph);

	sampler->destroy_graph(sampled_graph); 
	graph_loader->destroy_graph(graph);
	delete sampler;
	delete partitioner;
	delete mapping;
	delete graph_loader;
	delete graph_destroyer;
	delete csr_constructor;

	Debug::get_instance()->leave_function("main");
	NetworkInterface::finalize_instance();
	return 0;
}
