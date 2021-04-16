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
class ChainMiningInstance: public ZGraphInstance<EdgeData> {
	private:
		int order;

		bool canonicality_check(const MyList<Edge<EdgeData>> &edges, const Edge<EdgeData> new_edge) { // increasement canonicality check
			if (edges.list_size <= 1) return true;
			bool found_neighbour = false;
			for (int i = 0; i < edges.list_size; ++ i) {
				if (found_neighbour == false) {
					if (new_edge.src == edges.list[i].src || new_edge.src == edges.list[i].dst 
							|| new_edge.dst == edges.list[i].src || new_edge.dst == edges.list[i].dst) {
						found_neighbour = true;
					}
				} else {
					if (edges.list[i].edge_id >= new_edge.edge_id) {
						//Debug::get_instance()->log("****** failed ******");
						return false;
					}
				}
			}
			return true;
		}
	public:
		ChainMiningInstance(Graph<EdgeData> * graph_, int order_): ZGraphInstance<EdgeData>(graph_), order(order_) {
		}

		double estimate() {
			MyList<Edge<EdgeData>> edges;
			MyList<Edge<EdgeData>> first_edge;
			MyList<VertexId> vertices;
			SubGraph<EdgeData> subgraph;
			//const MyList<VertexId> &vertices = subgraph.get_vertices();
			long double p = 1.;

			std::pair<Edge<EdgeData>, double> tmp_0 = this->sample_edge();
			if (tmp_0.second < 0) return 0;
			p *= tmp_0.second;
			edges.push_back(tmp_0.first);
			first_edge.push_back(tmp_0.first);
			subgraph.add_edge(tmp_0.first);
			vertices.push_back(tmp_0.first.src);
			vertices.push_back(tmp_0.first.dst);

			for (int i = 0; i < order - 2; ++ i) { // two edges remained
				std::pair<Edge<EdgeData>, double> tmp = this->conditional_sample_edge(first_edge, vertices);
				if (tmp.second < 0) return 0;
				p *= tmp.second;
				if (canonicality_check(edges, tmp.first) == false) return 0;
				edges.push_back(tmp.first);
				subgraph.add_edge(tmp.first);
				if (subgraph.get_vertices().list_size != i + 3) return 0;
				if (vertices.list[0] == tmp.first.src) {
					vertices.list[0] = tmp.first.dst;
				} else if (vertices.list[0] == tmp.first.dst) {
					vertices.list[0] = tmp.first.src;
				} else if (vertices.list[1] == tmp.first.src) {
					vertices.list[1] = tmp.first.dst;
				} else if (vertices.list[1] == tmp.first.dst) {
					vertices.list[1] = tmp.first.src;
				} else {
					assert(false);
				}
			}
			//Debug::get_instance()->log(subgraph);
			return 1. / p;
		}
};

int main(int argc, char ** argv) {
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2 && argc != 3) {
		Debug::get_instance()->print("usage: ./ChainMining [graph-path] [order, default = 4, should not smaller than 3]");
		return -1;
	}

	int order;
	if (argc == 2) {
		order = 4;
	} else if (argc == 3) {
		order = std::atoi(argv[2]);
		if (order < 3) {
			Debug::get_instance()->print("order should not be smaller than 3.");
			return -1;
		}
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

	ChainMiningInstance<Empty> * counting_instance = new ChainMiningInstance<Empty>(&graph, order);
	double cnt;

	ErrorProfile<Empty> error_profile;
	double required_error_rate = 5. / 100.;
	double required_confidence = 0.95;
	Estimation init_estimation;

	//TimeProfile<Empty> time_profile;
	//unsigned int time_limit = 3600000;
	
	//unsigned long long num_estimators = 1 << 16;
	//Estimation init_estimation;

	VertexId num_vertices_in_pattern = order;
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	
	cnt = error_profile.run(*counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	//cnt = time_profile.run(*counting_instance, num_vertices_in_pattern, time_limit);
	//cnt = counting_instance->calculate(num_estimators, num_vertices_in_pattern, init_estimation);
	Debug::get_instance()->print("estimated number of ", order, "-chains: ", cnt);

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
