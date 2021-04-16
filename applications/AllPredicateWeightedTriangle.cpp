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
class WeightedTriangleCoutingInstance: public AllPredicateZGraphInstance<EdgeData> {
	public:
		WeightedTriangleCoutingInstance(Graph<EdgeData> *graph_, std::function<bool(EdgeData)> predicate_): AllPredicateZGraphInstance<EdgeData>(graph_, predicate_) {
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

int main(int argc, char ** argv) {
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./AllPredicateWeightedTriangle [graph-path]");
		return -1;
	}

	// getting graph meta data
	std::string graph_path(argv[1]);
	VertexId num_vertices;
	load_graph_meta_data(graph_path + ".csv", num_vertices);
	Debug::get_instance()->log("num_vertices = ", num_vertices);

	// loading graph
	UniformGraphPartitioner * partitioner = new UniformGraphPartitioner(num_vertices); 
	unsigned int seed = get_globally_consistent_seed();
	Debug::get_instance()->log("seed = ", seed);
	RandomMapping * mapping = new RandomMapping(seed, num_vertices);
	CSRConstructor<VertexId> * csr_constructor = new CSRConstructor<VertexId>();
	NumaAwareGraphDestroyer<VertexId> * graph_destroyer = new NumaAwareGraphDestroyer<VertexId>();
	DistGraphLoader<VertexId> * graph_loader = new DistGraphLoader<VertexId>(partitioner, mapping, csr_constructor, graph_destroyer);
	Graph<VertexId> graph;
	graph_loader->load_graph(graph_path + ".biedgelist", num_vertices, graph);

	auto predicate = [&](VertexId edge_data) {
		return edge_data >= 50;
	};
	WeightedTriangleCoutingInstance<VertexId> * weighted_triangle_counting_instance = new WeightedTriangleCoutingInstance<VertexId>(&graph, predicate);
	VertexId num_vertices_in_pattern = 3;
	double cnt;

	//TimeProfile<VertexId> time_profile;
	//unsigned int time_limit = 1000;

	ErrorProfile<VertexId> error_profile;
	double required_error_rate = 5 / 100.;
	double required_confidence = 0.95;
	Estimation init_estimation;

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");

	//cnt = time_profile.run(*weighted_triangle_counting_instance, num_vertices_in_pattern, time_limit);
	cnt = error_profile.run(*weighted_triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);

	Timer::timer_stop("estimation");

	if (NetworkInterface::get_instance()->get_partition_id() == 0) {
		Debug::get_instance()->print("estimated number of triangles: ", (unsigned long long) cnt);
	}

	delete weighted_triangle_counting_instance;
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
