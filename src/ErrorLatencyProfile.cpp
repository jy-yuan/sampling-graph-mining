#include <math.h>

#include "ErrorLatencyProfile.hpp"
#include "Time.hpp"
#include "SharedMemSys.hpp"
#include "Regression.hpp"

#define BASIC_NUM_ESTIMATORS_PER_THREAD 64

// TimeProfile

template<typename EdgeData>
void TimeProfile<EdgeData>::reset() {
	num_finished_estimators = 0;
	current_estimation = 0;
};

template<typename EdgeData>
void TimeProfile<EdgeData>::update(unsigned long long new_num_estimators, double new_estimation) {
	current_estimation = current_estimation * double(num_finished_estimators) / double(num_finished_estimators + new_num_estimators) + new_estimation * double(new_num_estimators) / double(num_finished_estimators + new_num_estimators);
	num_finished_estimators += new_num_estimators;
}

template<typename EdgeData>
double TimeProfile<EdgeData>::run(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, unsigned int time_limit) {
	reset();
	unsigned int used_time = 0;
	unsigned int curr_used_time = 0;
	unsigned long long next_num_estimators = BASIC_NUM_ESTIMATORS_PER_THREAD * SharedMemSys::get_instance()->get_num_threads_per_socket();
	while (used_time + curr_used_time * 2 <= time_limit) {
		MPI_Barrier(MPI_COMM_WORLD);
		double tmp = -get_time();
		Estimation init_estimation;
		double new_estimation = zgraph_instance.calculate(next_num_estimators, num_vertices_in_pattern, init_estimation);
		tmp += get_time();
		update(next_num_estimators, new_estimation);
		next_num_estimators <<= 1;
		curr_used_time = (unsigned int) (tmp * 1000);
		used_time += curr_used_time;
		Debug::get_instance()->log("current estimation: ", current_estimation);
	}
	Estimation init_estimation;
	unsigned long long remained_num_estimators = next_num_estimators * (time_limit - used_time) / (curr_used_time * 2);
	double new_estimation = zgraph_instance.calculate(remained_num_estimators, num_vertices_in_pattern, init_estimation);
	update(remained_num_estimators, new_estimation);
	return current_estimation;
}

// ErrorProfile

template<typename EdgeData>
double ErrorProfile<EdgeData>::calculate_bound(double epsilon, double delta, double num_sample_edge_calls, double num_conditional_sample_edge_calls, double estimation, EdgeId max_degree, EdgeId num_edges) { // this bound is not normalized
	// 1 / epsilon^2 * Delta^a * m^b / f(G) * ln(2 / delta)
	//double factor_0 = (2. + epsilon) / pow(epsilon, 2.);
	double factor_0 = 1. / pow(epsilon, 2.);
	double factor_1 = pow(double(max_degree), num_conditional_sample_edge_calls) * pow(double(num_edges), num_sample_edge_calls) / estimation;
	double factor_2 = log(2. / delta);
	return (factor_0 * factor_1 * factor_2);
}

template<typename EdgeData> 
double ErrorProfile<EdgeData>::run(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, double required_error_rate, double required_confidence, Estimation &init_estimation) { 
	Debug::get_instance()->enter_function("ErrorProfile::run");
	unsigned long long num_estimators_original_graph = get_num_estimators(zgraph_instance, num_vertices_in_pattern, required_error_rate, required_confidence);
	double accurate_estimation_original_graph = zgraph_instance.calculate(num_estimators_original_graph, num_vertices_in_pattern, init_estimation);

	Debug::get_instance()->leave_function("ErrorProfile::run");
	return accurate_estimation_original_graph;
}

template<typename EdgeData>
unsigned long long ErrorProfile<EdgeData>::get_num_estimators(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, double required_error_rate, double required_confidence) {
	Debug::get_instance()->enter_function("ErrorProfile::get_num_estimators");
	// [1] sampling a smaller graph G' from the original graph 
	// [2] calculate the f(G')
	// [3] approximate f(G) with f(G') * (1 / P), P = 0.05 by default
	// [4] get num_conditional_sample_edges call && num_sample_edge_calls
	// [5] get num_edges && max_degree
	// [6] calculate the bounding && get the required num_estimators for original graph
	// [7] get more accurate f(G) 
	
	// sampling a smaller graph
	long double sample_rate = 0.1; // TODO should be adjusted
	double num_sample_edge_calls = zgraph_instance.get_num_sample_edge_calls(); 
	double num_conditional_sample_edge_calls = zgraph_instance.get_num_conditional_sample_edge_calls();
	MPI_Allreduce(MPI_IN_PLACE, &num_sample_edge_calls, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	num_sample_edge_calls /= NetworkInterface::get_instance()->get_num_partitions();
	MPI_Allreduce(MPI_IN_PLACE, &num_conditional_sample_edge_calls, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	num_conditional_sample_edge_calls /= NetworkInterface::get_instance()->get_num_partitions();
	//num_conditional_sample_edge_calls = 3;
	CSRConstructor<EdgeData> * csr_constructor = new CSRConstructor<EdgeData>();
	NumaAwareGraphDestroyer<EdgeData> * graph_destroyer = new NumaAwareGraphDestroyer<EdgeData>();
	GraphSampler<EdgeData> * graph_sampler = new GraphSampler<EdgeData>(csr_constructor, graph_destroyer, sample_rate);
	Graph<EdgeData>& original_graph = *zgraph_instance.get_graph();
	Graph<EdgeData> sampled_graph;
	graph_sampler->sample_graph(original_graph, sampled_graph);
	zgraph_instance.set_graph(&sampled_graph);

	// true_bound = K * not_normalized_bound
	//const unsigned long long num_pre_estimators = 64;
	const unsigned long long num_post_estimators = 1 << 23;
	const unsigned long long num_data_points = 1 << 14;
	const unsigned long long num_estimators_per_data_points = 16;
	//const unsigned long long num_estimators_sampled_graph = (unsigned long long) num_data_points * num_estimators_per_data_points * 2;
	long double X[num_data_points];
	long double Y_tmp[num_data_points];
	long double Y[num_data_points];
	long double estimations[num_data_points];
	Estimation estimation;

	//{
	//	double new_estimation = zgraph_instance.calculate(num_pre_estimators, num_vertices_in_pattern);
	//	estimation.update(num_pre_estimators, new_estimation);
	//}
	for (unsigned long long i = 0; i < num_data_points; ++ i) {
		Estimation local_init_estimation;
		long double new_estimation = zgraph_instance.calculate(num_estimators_per_data_points, num_vertices_in_pattern, local_init_estimation);
		estimation.update(num_estimators_per_data_points, new_estimation);
		Y_tmp[i] = estimation.get_num_finished_estimators();
		estimations[i] = estimation.get_currnet_estimation();
	}
	{
		Estimation local_init_estimation;
		long double new_estimation = zgraph_instance.calculate(num_post_estimators, num_vertices_in_pattern, local_init_estimation);
		estimation.update(num_post_estimators, new_estimation);
	}
	unsigned long long used_num_data_points = 0;
	long double error_moving_average = 100; // a very large number
	int start = 0;
	double error_times = 1;
	double alpha = 0.04;
	for (unsigned long long i = 0; i < num_data_points; ++ i) {
		//double error = abs((estimations[i] - estimation.get_currnet_estimation()) / estimation.get_currnet_estimation());
		long double error = (estimations[i] - estimation.get_currnet_estimation()) / estimation.get_currnet_estimation() * 5;
		if (error < 0) error = -error;
		//error_moving_average = error;
		error_moving_average = error_moving_average * (1-alpha)+ error * alpha; // moving average is very important!!!
		double X_tmp = calculate_bound(error_moving_average * error_times, 1 - required_confidence, num_sample_edge_calls, num_conditional_sample_edge_calls, estimation.get_currnet_estimation(), sampled_graph.max_degree, sampled_graph.num_edges) * pow(error_moving_average*error_times, 2);
		//if (error_moving_average >= 0.3 && error_moving_average <= 0.95) {
			X[used_num_data_points] = sqrt(X_tmp/((i+1)*num_estimators_per_data_points));
			Y[used_num_data_points] = error_moving_average;
			if (NetworkInterface::get_instance()->get_partition_id() == 0) {
				Debug::get_instance()->log(error_moving_average, " [", X[used_num_data_points], ", ", Y[used_num_data_points], "]");
			}
			++ used_num_data_points;
		//}
	}
	assert(used_num_data_points > 0);
	Debug::get_instance()->log("confidence: ", required_confidence, " sample_edge_calls: ", num_sample_edge_calls, " num_conditional_sample_edge_calls: ", num_conditional_sample_edge_calls, " estimation: ", estimation.get_currnet_estimation(), " max_degree: ", sampled_graph.max_degree, " num_edges: ", sampled_graph.num_edges);
	Debug::get_instance()->log("accurate estimation on sampled graph: ", estimation.get_currnet_estimation());
	UnbiasedLinearRegression regression;
	Debug::get_instance()->log("used_num_data_points: ", used_num_data_points);
	regression.solve(X, Y, used_num_data_points);
	//{
	//	std::string log_msg = "";
	//	for (unsigned long long i = 0; i < num_data_points; ++ i) {
	//		log_msg += "(" + std::to_string(X[i]) + ", " + std::to_string(Y[i]) + ") ";
	//	}
	//	Debug::get_instance()->log(log_msg);
	//}
	double K = regression.get_k();
	MPI_Allreduce(MPI_IN_PLACE, &K, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	K /= NetworkInterface::get_instance()->get_num_partitions();

	//unsigned long long num_estimators_sampled_graph = (unsigned long long) 1 << 23; // this parameter could be ajudsted 
	//Debug::get_instance()->log("running estimation on sampled graph...");
	//double estimation_sampled_graph = zgraph_instance.calculate(num_estimators_sampled_graph, num_vertices_in_pattern);
	//Debug::get_instance()->log("estimation_sampled_graph: ", (int) estimation_sampled_graph);
	//double rough_estimation_original_graph = estimation_sampled_graph / pow(sample_rate, num_sample_edge_calls + num_conditional_sample_edge_calls);
	long double rough_estimation_original_graph = estimation.get_currnet_estimation() / pow(sample_rate, num_sample_edge_calls + num_conditional_sample_edge_calls);
	Debug::get_instance()->log("rough_estimation_original_graph: ", rough_estimation_original_graph);
	delete csr_constructor;
	graph_sampler->destroy_graph(sampled_graph);
	delete graph_sampler;
	delete graph_destroyer;

	zgraph_instance.set_graph(&original_graph);
	//int num_sample_edge_calls = 1;
	//int num_conditional_sample_edge_calls = 1;

	EdgeId num_edges = original_graph.num_edges;
	EdgeId max_degree = original_graph.max_degree;

	double theoretic_k = 2 + required_error_rate;
	for (int i = 0; i < num_conditional_sample_edge_calls; ++ i) {
		theoretic_k *= (i + 2);
	}
	Debug::get_instance()->log("theoretic_k = ", theoretic_k);
	Debug::get_instance()->log("K^2 = ", K*K);

	double num_estimators_original_graph = calculate_bound(required_error_rate, 1. - required_confidence, num_sample_edge_calls, num_conditional_sample_edge_calls, rough_estimation_original_graph, max_degree, num_edges);
	Debug::get_instance()->log("num_estimators required: ", num_estimators_original_graph);
	//num_estimators_original_graph /= 610;
	//if (K*K < theoretic_k) {
		num_estimators_original_graph *= (K*K); 
	//} else {
	//	num_estimators_original_graph *= theoretic_k;
	//}
	MPI_Allreduce(MPI_IN_PLACE, &num_estimators_original_graph, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	num_estimators_original_graph /= NetworkInterface::get_instance()->get_num_partitions();
	Debug::get_instance()->print("num_estimators required: ", num_estimators_original_graph);
	Debug::get_instance()->log("running estimation on original graph...");

	Debug::get_instance()->leave_function("ErrorProfile::get_num_estimators");

	return (unsigned long long)num_estimators_original_graph;
}

template class TimeProfile<Empty>;
template class TimeProfile<VertexId>;

template class ErrorProfile<Empty>;
template class ErrorProfile<VertexId>;
