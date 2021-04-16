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
#include "ErrorLatencyProfile.hpp"

template<typename EdgeData>
SubGraph<EdgeData> FourClique(const MyList<VertexId> &vertices) {
	//Debug::get_instance()->log("num vertices in sampled_subgraph: ", vertices.list_size);
	assert(vertices.list_size == 4);
	SubGraph<EdgeData> four_clique;
	Edge<EdgeData> edge;
	for (int i = 0; i < vertices.list_size; ++ i) {
		VertexId v_i = vertices.list[i];
		for (int j = i + 1; j < vertices.list_size; ++ j) {
			VertexId v_j = vertices.list[j];
			edge.src = v_i;
			edge.dst = v_j;
			four_clique.add_edge(edge);
		}
	}
	return four_clique;
}

template<typename EdgeData> // TODO to check maybe some bugs here
class FourCliqueCountingInstance: public ZGraphInstance<EdgeData> {
	public:
		FourCliqueCountingInstance(Graph<EdgeData> *graph_): ZGraphInstance<EdgeData>(graph_) {
		}

		double estimate() {
			SubGraph<EdgeData> sampled_subgraph;
			double p = 1;

			std::pair<Edge<EdgeData>, double> tmp_0 = this->sample_edge();
			if (tmp_0.second < 0) {
				return 0;
			}
			p *= tmp_0.second;
			sampled_subgraph.add_edge(tmp_0.first);

			std::pair<Edge<EdgeData>, double> tmp_1 = this->conditional_sample_edge(sampled_subgraph); 
			if (tmp_1.second < 0) {
				return 0;
			}
			p *= tmp_1.second;
			sampled_subgraph.add_edge(tmp_1.first);

			std::pair<Edge<EdgeData>, double> tmp_2 = this->conditional_sample_edge(sampled_subgraph);
			if (tmp_2.second < 0) {
				return 0;
			}
			p *= tmp_2.second;
			sampled_subgraph.add_edge(tmp_2.first);

			if (sampled_subgraph.get_vertices().list_size == 3) { // a triangle is sampled
				return 0;
			}

			//Debug::get_instance()->log("sampled_subgraph: ", sampled_subgraph);
			SubGraph<EdgeData> subgraph_to_wait = FourClique<Empty>(sampled_subgraph.get_vertices()) - sampled_subgraph;
			if (this->conditional_close(sampled_subgraph, subgraph_to_wait)) {
				return 1. / p;
			}
			return 0;
		}
};

int main(int argc, char ** argv) {
	//NetworkInterface * net_instance = NetworkInterface::get_instance(&argc, &argv);
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	//SharedMemSys * shared_mem_sys = SharedMemSys::get_instance();
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./FourClique [graph-path]");
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

	FourCliqueCountingInstance<Empty> * four_clique_counting_instance = new FourCliqueCountingInstance<Empty>(&graph);
	Estimation init_estimation;
	//TimeProfile<Empty> time_profile;

	VertexId num_vertices_in_pattern = 4;
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	unsigned long long num_estmators = (unsigned long long) 1 << 23;
	double cnt = four_clique_counting_instance->calculate(num_estmators, num_vertices_in_pattern, init_estimation);
	//unsigned int time_limit = 5000;
	//double cnt = time_profile.run(*four_clique_counting_instance, num_vertices_in_pattern, time_limit);
	Timer::timer_stop("estimation");
	Debug::get_instance()->log("estimated number of four_clique: ", (unsigned long long) cnt);

	delete four_clique_counting_instance;

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
