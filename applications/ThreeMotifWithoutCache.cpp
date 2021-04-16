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
SubGraph<EdgeData> Triangle(const Edge<EdgeData> &edge_0, const Edge<EdgeData> &edge_1) {
	SubGraph<EdgeData> triangle;
	triangle.add_edge(edge_0);
	triangle.add_edge(edge_1);
	Edge<EdgeData> edge;
	if (edge_0.src == edge_1.src) {
		edge.src = edge_0.dst;
		edge.dst = edge_1.dst;
	} else if (edge_0.src == edge_1.dst) {
		edge.src = edge_0.dst;
		edge.dst = edge_1.src;
	} else if (edge_0.dst == edge_1.src) {
		edge.src = edge_0.src;
		edge.dst = edge_1.dst;
	} else if (edge_0.dst == edge_1.dst) {
		edge.src = edge_0.src;
		edge.dst = edge_1.src;
	} else {
		assert(false);
	}
	triangle.add_edge(edge);
	return triangle;
}

template<typename EdgeData>
class TriangleCoutingInstance: public ZGraphInstance<EdgeData> {
	public:
		TriangleCoutingInstance(Graph<EdgeData> *graph_): ZGraphInstance<EdgeData>(graph_) {
		}

		double estimate() {
			SubGraph<EdgeData> sampled_subgraph;
			double p = 1;

			std::pair<Edge<EdgeData>, double> tmp0 = this->sample_edge();
			if (tmp0.second < 0) { 
				return 0;
			}
			p *= tmp0.second;
			sampled_subgraph.add_edge(tmp0.first);

			std::pair<Edge<EdgeData>, double> tmp1 = this->conditional_sample_edge(sampled_subgraph);
			if (tmp1.second < 0) {
				return 0;
			}
			p *= tmp1.second;
			sampled_subgraph.add_edge(tmp1.first);

			SubGraph<EdgeData> subgraph_to_wait = Triangle(tmp0.first, tmp1.first) - sampled_subgraph;
			if (this->conditional_close(sampled_subgraph, subgraph_to_wait)) {
				return 1. / p;
			}
			return 0;
		}
};

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
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./ThreeMotifWithoutCache [graph-path]");
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
	TriangleCoutingInstance<Empty> * triangle_counting_instance = new TriangleCoutingInstance<Empty>(&graph);
	VertexId num_vertices_in_pattern = 3;

	ErrorProfile<Empty> error_profile;
	double required_error_rate = 5. / 100;
	double required_confidence = 0.95;
	//unsigned long long num_estmators_three_chain = 1 << 10;
	//unsigned long long num_estmators_triangle = 1 << 20;
	Timer::timer_start("estimation");
	Timer::timer_start("error_profile");
	unsigned long long num_estmators_three_chain = error_profile.get_num_estimators(*three_chain_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence);
	unsigned long long num_estmators_triangle = error_profile.get_num_estimators(*triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence);
	Timer::timer_stop("error_profile");

	//unsigned long long num_estmators = (unsigned long long) 1 << 20;
	Estimation init_estimation_three_chain;
	Estimation init_estimation_triangle;

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation(3-chain)");
	double cnt_three_chain = three_chain_counting_instance->calculate(num_estmators_three_chain, num_vertices_in_pattern, init_estimation_three_chain);
	Timer::timer_stop("estimation(3-chain)");

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation(triangle)");
	double cnt_triangle = triangle_counting_instance->calculate(num_estmators_triangle, num_vertices_in_pattern, init_estimation_triangle);
	Timer::timer_stop("estimation(triangle)");
	Timer::timer_stop("estimation");

	Debug::get_instance()->print("estimated number of 3-chains: ", (unsigned long long) cnt_three_chain);
	Debug::get_instance()->print("estimated number of triangles: ", (unsigned long long) cnt_triangle);

	delete three_chain_counting_instance; 
	delete triangle_counting_instance;

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
