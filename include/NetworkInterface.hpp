#ifndef NETWORK_INTERFACE
#define NETWORK_INTERFACE

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <stdint.h>

enum MessageTag{
	TransferSeed,
	ShuffleGraph,
	ReduceEstimation
};

template<typename T>
MPI_Datatype get_mpi_data_type() {
	if (std::is_same<T, char>::value) {
		return MPI_CHAR;
	} else if (std::is_same<T, unsigned char>::value) {
		return MPI_UNSIGNED_CHAR;
	} else if (std::is_same<T, int>::value) {
		return MPI_INT;
	} else if (std::is_same<T, unsigned>::value) {
		return MPI_UNSIGNED;
	} else if (std::is_same<T, long>::value) {
		return MPI_LONG;
	} else if (std::is_same<T, unsigned long>::value) {
		return MPI_UNSIGNED_LONG;
	} else if (std::is_same<T, float>::value) {
		return MPI_FLOAT;
	} else if (std::is_same<T, double>::value) {
		return MPI_DOUBLE;
	} else {
		printf("Type not supported\n");
		exit(-1);
	}
}

class NetworkInterface {
	private:
		static NetworkInterface * instance;
		int partition_id;
		int partitions;
		NetworkInterface(int * argc, char *** argv);
		~NetworkInterface();
	public:
		static void init_network_interface(int * argc, char *** argv);
		static NetworkInterface * get_instance(int * argc = nullptr, char *** argv = nullptr);
		static void finalize_instance();

		void pause();
		int get_partition_id();
		int get_num_partitions();
};

#endif
