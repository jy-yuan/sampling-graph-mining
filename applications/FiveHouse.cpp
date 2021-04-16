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
#include "SamplerGenerator.hpp"

Pattern* get_pattern() {
	Edge<Empty> edges[6];
	edges[0].src = 0, edges[0].dst = 1, edges[0].edge_id = 0;
	edges[1].src = 0, edges[1].dst = 2, edges[1].edge_id = 1;
	edges[2].src = 1, edges[2].dst = 2, edges[2].edge_id = 2;
	edges[3].src = 1, edges[3].dst = 3, edges[3].edge_id = 3;
	edges[4].src = 2, edges[4].dst = 4, edges[4].edge_id = 4;
	edges[5].src = 3, edges[5].dst = 4, edges[5].edge_id = 5;
	return new Pattern(5, 6, edges);
}

int main(int argc, char ** argv) { 
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./FiveHouse [graph-path]");
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

	Pattern * pattern = get_pattern();
	SamplerGenerator<Empty> * counting_instance = new SamplerGenerator<Empty>(&graph, pattern);
	ErrorProfile<Empty> error_profile;

	Estimation init_estimation;
	VertexId num_vertices_in_pattern = 5;
	double cnt;

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
//	unsigned long long num_estmators = (unsigned long long) 1 << 25; 
//	for (unsigned int i = 0; i < (1  << 20); ++ i) {
//		num_estmators += (unsigned long long)(1 << 25);
//		cnt = counting_instance->calculate(num_estmators, num_vertices_in_pattern, init_estimation);
//		double ground_truth = 2134981304;
//		double error = cnt / ground_truth - 1.;
//		if (error < 0) error = -error;
//		Debug::get_instance()->print("estimated number of 5-house: ", cnt, " error:", error);
//	}
	
	//TimeProfile<Empty> time_profile;
	//unsigned int time_limit = 28800000;
	//cnt = time_profile.run(*counting_instance, num_vertices_in_pattern, time_limit);
	
	double required_error_rate = 5 / 100.;
	double required_confidence = 0.95;

	cnt = error_profile.run(*counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	//Debug::get_instance()->log("estimated number of triangles: ", (unsigned long long) cnt);

	//cnt = error_profile.run(*triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	Debug::get_instance()->print("estimated number of 5-house: ", cnt);

	Timer::timer_stop("estimation");

	delete counting_instance;

	graph_loader->destroy_graph(graph);
	delete partitioner;
	delete mapping;
	delete csr_constructor;
	delete graph_loader;
	delete graph_destroyer;
	delete pattern;

	Debug::get_instance()->leave_function("main");
	MPI_Barrier(MPI_COMM_WORLD);
	std::this_thread::sleep_for(std::chrono::seconds(1));
	Timer::report_timers();
	NetworkInterface::finalize_instance();
	return 0;
}
