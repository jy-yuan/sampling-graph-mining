#include <assert.h>
#include <numa.h>
#include <omp.h>
#include <string.h>

#include <string>

#include "SharedMemSys.hpp"
#include "Debug.hpp"

SharedMemSys* SharedMemSys::instance = nullptr;

SharedMemSys::SharedMemSys() {
	assert(numa_available() != -1);
	assert(sizeof(unsigned long) == 8);

	num_threads = numa_num_configured_cpus(); 
	//num_threads = 12;
	num_sockets = numa_num_configured_nodes();
	//assert(num_sockets==1);
	//num_sockets = 1; 
	num_threads_per_socket = num_threads / num_sockets;

	char nodestring[num_sockets * 2 + 1];
	memset(nodestring, 0, sizeof(nodestring));
	nodestring[0] = '0';
	for (int s_i = 1; s_i < num_sockets; ++ s_i) {
		nodestring[s_i * 2 - 1] = ',';
		nodestring[s_i * 2] = '0' + s_i;
	}
	struct bitmask * nodemask = numa_parse_nodestring(nodestring);
	numa_set_interleave_mask(nodemask);

	omp_set_dynamic(0);
	omp_set_num_threads(num_threads);
#pragma omp parallel for 
	for (int t_i = 0; t_i < num_threads; ++ t_i) {
		int s_i = get_socket_id(t_i);
		assert(numa_run_on_node(s_i) == 0);
	}
}

void SharedMemSys::init_shared_mem_sys() {
	assert(instance == nullptr);
	instance = new SharedMemSys();
}

SharedMemSys * SharedMemSys::get_instance() {
	if (instance == nullptr) {
		instance = new SharedMemSys();
	}
	return instance;
}

int SharedMemSys::get_socket_id(int thread_id) {
	return thread_id / num_threads_per_socket;
}

int SharedMemSys::get_socket_offset(int thread_id) {
	return thread_id % num_threads_per_socket;
}

int SharedMemSys::get_num_sockets() {
	return num_sockets;
}

int SharedMemSys::get_num_threads() {
	return num_threads;
}

int SharedMemSys::get_current_thread_id() {
	return omp_get_thread_num();
}

int SharedMemSys::get_current_socket_id() {
	return get_socket_id(omp_get_thread_num());
}

int SharedMemSys::get_num_threads_per_socket() {
	return num_threads_per_socket;
}

int SharedMemSys::get_current_socket_offset() {
	return omp_get_thread_num() % num_threads_per_socket;
}
