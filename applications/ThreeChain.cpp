#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <utility>
#include <thread>
#include <chrono>

#include "ZGraphInstance.hpp"
#include "Debug.hpp"
#include "Utilities.hpp"
#include "SharedMemSys.hpp"
#include "Time.hpp"

template<typename EdgeData>
class ThreeChainCountingInstance: public ZGraphInstance<EdgeData> {
	public:
		ThreeChainCountingInstance(Graph<EdgeData> *graph_): ZGraphInstance<EdgeData>(graph_) {
		}

		double estimate() {
			SubGraph<EdgeData> sampled_subgraph;
			double p = 1;
			std::pair<Edge<EdgeData>, double> tmp = this->sample_edge();
			if (tmp.second < 0) {
				return 0;
			}
			p *= tmp.second;
			sampled_subgraph.add_edge(tmp.first);
			tmp = this->conditional_sample_edge(sampled_subgraph); 
			if (tmp.second < 0) {
				return 0;
			}
			p *= tmp.second;
			return 1. / p;
		}
};

int main(int argc, char ** argv) {
	//NetworkInterface * net_instance = NetworkInterface::get_instance(&argc, &argv);
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	//SharedMemSys * shared_mem_sys = SharedMemSys::get_instance();
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./ThreeChain [graph-path]");
		return -1;
	}

	std::string graph_path(argv[1]);
	VertexId num_vertices;
	load_graph_meta_data(graph_path + ".csv", num_vertices);
	Debug::get_instance()->log("num_vertices = ", num_vertices);

	UniformGraphPartitioner * partitioner = new UniformGraphPartitioner(num_vertices); // loading graph...
	unsigned int seed = get_globally_consistent_seed();
	Debug::get_instance()->log("seed = ", seed);
	RandomMapping * mapping = new RandomMapping(seed, num_vertices);
	CSRConstructor<Empty> * csr_constructor = new CSRConstructor<Empty>();
	NumaAwareGraphDestroyer<Empty> * graph_destroyer = new NumaAwareGraphDestroyer<Empty>();
	DistGraphLoader<Empty> * graph_loader = new DistGraphLoader<Empty>(partitioner, mapping, csr_constructor, graph_destroyer);
	Graph<Empty> graph;
	graph_loader->load_graph(graph_path + ".biedgelist", num_vertices, graph);

	ThreeChainCountingInstance<Empty> * three_chain_counting_instance = new ThreeChainCountingInstance<Empty>(&graph);

	unsigned long long num_estmators = (unsigned long long) 1 << 25;
	VertexId num_vertices_in_pattern = 3;
	Estimation init_estimation;
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	double cnt = three_chain_counting_instance->calculate(num_estmators, num_vertices_in_pattern, init_estimation);
	Timer::timer_stop("estimation");
	Debug::get_instance()->log("estimated number of 3-chain: ", (unsigned long long) cnt);

	delete three_chain_counting_instance;

	graph_loader->destroy_graph(graph);
	delete partitioner;
	delete mapping;
	delete csr_constructor;
	delete graph_loader;
	delete graph_destroyer;

	Debug::get_instance()->leave_function("main");
	MPI_Barrier(MPI_COMM_WORLD);
	std::this_thread::sleep_for(std::chrono::seconds(1));
	Timer::report_timers();
	NetworkInterface::finalize_instance();
	return 0;
}
