#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#include <fstream>

#include "Utilities.hpp"
#include "NetworkInterface.hpp"

bool file_exists(std::string filename) {
	struct stat st;
	return stat(filename.c_str(), &st) == 0;
}

long file_size(std::string filename) {
	struct stat st;
	assert(stat(filename.c_str(), &st) == 0);
	return st.st_size;
}

void load_graph_meta_data(std::string meta_path, VertexId &num_vertices) {
	std::ifstream fin(meta_path);
	std::string str;
	getline(fin, str);
	getline(fin, str);
	std::string tmp_str = "";
	int pos = 0;
	for (; str[pos] < '0' || str[pos] > '9'; ++ pos);
	for (; str[pos] >= '0' && str[pos] <= '9'; tmp_str += str[pos++]);
	num_vertices = std::stoi(tmp_str);
	fin.close();
}

unsigned int get_globally_consistent_seed() {
	unsigned int seed;
	if (NetworkInterface::get_instance()->get_partition_id() == 0) {
		seed = time(NULL);
		int num_partitions = NetworkInterface::get_instance()->get_num_partitions();
		for (int p_i = 1; p_i < num_partitions; ++ p_i) {
			MPI_Send(&seed, 1, MPI_UNSIGNED, p_i, TransferSeed, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&seed, 1, MPI_UNSIGNED, 0, TransferSeed, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	return seed;
}

