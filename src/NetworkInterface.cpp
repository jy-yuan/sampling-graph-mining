#include "NetworkInterface.hpp"

NetworkInterface * NetworkInterface::instance = nullptr;

NetworkInterface::NetworkInterface(int * argc, char *** argv) {
	int provided;
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &partition_id);
	MPI_Comm_size(MPI_COMM_WORLD, &partitions);
}

NetworkInterface::~NetworkInterface() {
	//printf("[%d] mpi finalized called\n", partition_id);
	MPI_Finalize();
	//printf("[%d] mpi finalized returned\n", partition_id);
}

void NetworkInterface::pause() {
	if (partition_id == 0) {
		getchar();
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

int NetworkInterface::get_partition_id() {
	return partition_id;
}

int NetworkInterface::get_num_partitions() {
	return partitions;
}

void NetworkInterface::init_network_interface(int * argc, char *** argv) {
	assert(instance == nullptr);
	assert(argc != nullptr && argv != nullptr);
	instance = new NetworkInterface(argc, argv);
}

NetworkInterface * NetworkInterface::get_instance(int * argc, char *** argv) {
	if (instance == nullptr) {
		assert(argc != nullptr && argv != nullptr);
		instance = new NetworkInterface(argc, argv);
	}
	return instance;
}

void NetworkInterface::finalize_instance() {
	delete instance;
}
