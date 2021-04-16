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
			//Debug::get_instance()->log(Triangle(tmp0.first, tmp1.first), sampled_subgraph, subgraph_to_wait); 
			if (this->conditional_close(sampled_subgraph, subgraph_to_wait)) {
				//Debug::get_instance()->log("Found!"); 
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
		Debug::get_instance()->log("usage: ./Triangle [graph-path]");
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

	TriangleCoutingInstance<Empty> * triangle_counting_instance = new TriangleCoutingInstance<Empty>(&graph);
	//TimeProfile<Empty> time_profile;
	ErrorProfile<Empty> error_profile;

	VertexId num_vertices_in_pattern = 3;
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	Estimation init_estimation;
	//unsigned long long num_estmators = (unsigned long long) 1 << 23;
	//double cnt = triangle_counting_instance->calculate(num_estmators, num_vertices_in_pattern, init_estimation);
	
	//unsigned int time_limit = 5000;
	//double cnt = time_profile.run(*triangle_counting_instance, num_vertices_in_pattern, time_limit);
	
	double required_error_rate = 5 / 100.;
	double required_confidence = 0.95;
	double cnt;


	cnt = error_profile.run(*triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	Debug::get_instance()->print("estimated number of triangles: ", (unsigned long long) cnt);
	
	//const unsigned long long num_data_points = 1 << 12;
	//unsigned long long num_estimators_per_data_points = 4;
	//unsigned long long num_estimators = 0;
	//double estimation[num_data_points];
	//double accurate_result = 7515023;
	//for (unsigned long long i = 0; i < num_data_points; ++ i) {
	//	num_estimators += num_estimators_per_data_points;
	//	estimation[i] = triangle_counting_instance->calculate(num_estimators, num_vertices_in_pattern, init_estimation);   
	//	if (NetworkInterface::get_instance()->get_partition_id() == 0) {
	//		double error = (estimation[i] - accurate_result) / estimation[i];
	//		if (error < 0) error = -error;
	//		printf("[%llu, %.4f]\n", num_estimators, error * 100);
	//	}
	//}


	Timer::timer_stop("estimation");

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
