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

Pattern* get_triangle_pattern() {
	Edge<Empty> edges[3];
	edges[0].src = 0, edges[0].dst = 1;
	edges[1].src = 0, edges[1].dst = 2;
	edges[2].src = 1, edges[2].dst = 2;
	return new Pattern(3, 3, edges);
}

int main(int argc, char ** argv) { 
	NetworkInterface::init_network_interface(&argc, &argv);
	Debug::get_instance()->enter_function("main");
	SharedMemSys::init_shared_mem_sys();

	if (argc != 2) {
		Debug::get_instance()->log("usage: ./Triangle2 [graph-path]");
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

	Pattern * triangle_pattern = get_triangle_pattern();
	SamplerGenerator<Empty> * triangle_counting_instance = new SamplerGenerator<Empty>(&graph, triangle_pattern);
	ErrorProfile<Empty> error_profile;

	Estimation init_estimation;
	VertexId num_vertices_in_pattern = 3;
	double cnt;

	MPI_Barrier(MPI_COMM_WORLD);
	Timer::timer_start("estimation");
	//unsigned long long num_estmators = (unsigned long long) 3000000; 
	//cnt = triangle_counting_instance->calculate(num_estmators, num_vertices_in_pattern, init_estimation);
	
	//TimeProfile<Empty> time_profile;
	//unsigned int time_limit = 100000;
	//cnt = time_profile.run(*triangle_counting_instance, num_vertices_in_pattern, time_limit);
	
	double required_error_rate = 5 / 100.;
	double required_confidence = 0.95;

	cnt = error_profile.run(*triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	//Debug::get_instance()->log("estimated number of triangles: ", (unsigned long long) cnt);

	//cnt = error_profile.run(*triangle_counting_instance, num_vertices_in_pattern, required_error_rate, required_confidence, init_estimation);
	Debug::get_instance()->print("estimated number of triangles: ", (unsigned long long) cnt);

	Timer::timer_stop("estimation");

	delete triangle_counting_instance;

	graph_loader->destroy_graph(graph);
	delete partitioner;
	delete mapping;
	delete csr_constructor;
	delete graph_loader;
	delete graph_destroyer;
	delete triangle_pattern;

	Debug::get_instance()->leave_function("main");
	MPI_Barrier(MPI_COMM_WORLD);
	std::this_thread::sleep_for(std::chrono::seconds(1));
	Timer::report_timers();
	NetworkInterface::finalize_instance();
	return 0;
}
