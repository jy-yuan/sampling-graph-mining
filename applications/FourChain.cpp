#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <utility>
#include <thread>
#include <chrono>
#include <string>

#include "ZGraphInstance.hpp"
#include "Debug.hpp"
#include "Utilities.hpp"
#include "SharedMemSys.hpp"
#include "Time.hpp"
#include "ErrorLatencyProfile.hpp"

template<typename EdgeData> 
class FourChainCountingInstance: public ZGraphInstance<EdgeData> {
	public:
		FourChainCountingInstance(Graph<EdgeData> * graph_): ZGraphInstance<EdgeData>(graph_) {
		}

		double estimate() {
			const int num_sampling_methods = 3;
			int selected_sampling_methods = this->sample_interger(1, num_sampling_methods);
			if (selected_sampling_methods == 1) { // (0) <-> (1) <-> (2) type
				MyList<Edge<EdgeData>> edges;
				MyList<VertexId> vertices;
				double p = 1. / num_sampling_methods;

				std::pair<Edge<EdgeData>, double> tmp_0 = this->sample_edge();
				if (tmp_0.second < 0) return 0;
				p *= tmp_0.second;
				edges.push_back(tmp_0.first);
				vertices.push_back(tmp_0.first.src);
				vertices.push_back(tmp_0.first.dst);

				std::pair<Edge<EdgeData>, double> tmp_1 = this->conditional_sample_edge(edges, vertices);
				if (tmp_1.second < 0) return 0;
				p *= tmp_1.second;
				edges.push_back(tmp_1.first);
				vertices.list_size = 1;
				if (tmp_0.first.src == tmp_1.first.src) {
					vertices.list[0] = tmp_1.first.dst;
				} else if (tmp_0.first.src == tmp_1.first.dst) {
					vertices.list[0] = tmp_1.first.src;
				} else if (tmp_0.first.dst == tmp_1.first.src) {
					vertices.list[0] = tmp_1.first.dst;
				} else if (tmp_0.first.dst == tmp_1.first.dst) {
					vertices.list[0] = tmp_1.first.src;
				} else {
					assert(false);
				}

				std::pair<Edge<EdgeData>, double> tmp_2 = this->conditional_sample_edge(edges, vertices);
				if (tmp_2.second < 0) return 0;
				p *= tmp_2.second;
				if (tmp_2.first.src == tmp_0.first.src || tmp_2.first.src == tmp_0.first.dst) return 0;
				if (tmp_2.first.dst == tmp_0.first.src || tmp_2.first.dst == tmp_0.first.dst) return 0;

				return 1. / p;
			} else if (selected_sampling_methods == 2) { // (0) <-> (2) <-> (1) type
				MyList<Edge<EdgeData>> edges;
				MyList<VertexId> vertices;
				double p = 1. / num_sampling_methods;

				std::pair<Edge<EdgeData>, double> tmp_0 = this->sample_edge();
				if (tmp_0.second < 0) return 0;
				p *= tmp_0.second;
				edges.push_back(tmp_0.first);
				vertices.push_back(tmp_0.first.src);
				vertices.push_back(tmp_0.first.dst);

				std::pair<Edge<EdgeData>, double> tmp_1 = this->conditional_sample_edge(edges, vertices);
				if (tmp_1.second < 0) return 0;
				p *= tmp_1.second;
				vertices.list_size = 1;
				if (tmp_0.first.src == tmp_1.first.src) {
					vertices.list[0] = tmp_1.first.dst;
				} else if (tmp_0.first.src == tmp_1.first.dst) {
					vertices.list[0] = tmp_1.first.src;
				} else if (tmp_0.first.dst == tmp_1.first.src) {
					vertices.list[0] = tmp_1.first.dst;
				} else if (tmp_0.first.dst == tmp_1.first.dst) {
					vertices.list[0] = tmp_1.first.src;
				} else {
					assert(false);
				}

				std::pair<Edge<EdgeData>, double> tmp_2 = this->conditional_sample_edge(edges, vertices);
				if (tmp_2.second < 0) return 0;
				p *= tmp_2.second;
				if (tmp_2.first.src == tmp_0.first.src || tmp_2.first.src == tmp_0.first.dst) return 0;
				if (tmp_2.first.dst == tmp_0.first.src || tmp_2.first.dst == tmp_0.first.dst) return 0;
				if (tmp_2.first.edge_id >= tmp_1.first.edge_id) return 0;

				return 1. / p;
			} else { // (1) <-> (0) <-> (2) type
				MyList<Edge<EdgeData>> edges;
				MyList<VertexId> vertices;
				double p = 1. / num_sampling_methods;

				std::pair<Edge<EdgeData>, double> tmp_0 = this->sample_edge();
				if (tmp_0.second < 0) return 0;
				p *= tmp_0.second;
				edges.push_back(tmp_0.first);

				vertices.push_back(tmp_0.first.src);
				std::pair<Edge<EdgeData>, double> tmp_1 = this->conditional_sample_edge(edges, vertices);
				if (tmp_1.second < 0) return 0;
				p *= tmp_1.second;

				vertices.list[0] = tmp_0.first.dst;
				std::pair<Edge<EdgeData>, double> tmp_2 = this->conditional_sample_edge(edges, vertices);
				if (tmp_2.second < 0) return 0;
				p *= tmp_2.second;

				if (tmp_2.first.src == tmp_1.first.src || tmp_2.first.src == tmp_1.first.dst) return 0;
				if (tmp_2.first.dst == tmp_1.first.src || tmp_2.first.dst == tmp_1.first.dst) return 0;

				return 1. / p;
			} 
		}
};

int main(int argc, char ** argv) {

	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./FourChain [graph-path]");
		return -1;
	}

	Debug::get_instance()->print("loading graph...");
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
	Debug::get_instance()->print("graph loaded.");

	FourChainCountingInstance<Empty> * counting_instance = new FourChainCountingInstance<Empty>(&graph);
	double cnt;

	ErrorProfile<Empty> error_profile;
	double required_error_rate = 5. / 100.;
	double required_confidence = 0.95;
	Estimation init_estimation;

	//TimeProfile<Empty> time_profile;
	//unsigned int time_limit = 1000;
	
	//unsigned long long num_estimators = 1 << 16;
	//Estimation init_estimation;

	VertexId num_vertices_in_pattern = 4;
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	
	cnt = error_profile.run(*counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	//cnt = time_profile.run(*counting_instance, num_vertices_in_pattern, time_limit);
	//cnt = counting_instance->calculate(num_estimators, num_vertices_in_pattern, init_estimation);
	Debug::get_instance()->print("estimated number of 4-chains: ", (unsigned long long) cnt);

	Timer::timer_stop("estimation");

	delete counting_instance;

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
