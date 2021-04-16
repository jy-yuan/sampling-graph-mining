#ifndef ERROR_LATENCY_PROFILE
#define ERROR_LATENCY_PROFILE

#include <math.h>
#include <assert.h>

#include <string>

#include "ZGraphInstance.hpp"
#include "Type.hpp"
#include "Time.hpp"

template<typename EdgeData>
class TimeProfile {
	private:
		unsigned long long num_finished_estimators;
		double current_estimation;

		void reset();
		void update(unsigned long long new_num_estimators, double new_estimation);
	public:
		double run(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, unsigned int time_limit); // time_limit(ms)
};

template<typename EdgeData>
class ErrorProfile {
	private:
		double calculate_bound(double epsilon, double delta, double num_sample_edge_calls, double num_conditional_sample_edge_call, double estimation, EdgeId max_degree, EdgeId num_edges);
	public:
		double run(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, double required_error_rate, double required_confidence, Estimation &init_estimation);
		unsigned long long get_num_estimators(ZGraphInstance<EdgeData> &zgraph_instance, VertexId num_vertices_in_pattern, double required_error_rate, double required_confidence);
};

#endif
