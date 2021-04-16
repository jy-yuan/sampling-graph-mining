#ifndef SHARED_MEM_SYS
#define SHARED_MEM_SYS

class SharedMemSys {
	private:
		static SharedMemSys * instance;
		int num_threads;
		int num_sockets;
		int num_threads_per_socket;

		SharedMemSys();
	public:
		static void init_shared_mem_sys();
		static SharedMemSys * get_instance();
		int get_socket_id(int thread_id);
		int get_socket_offset(int thread_id);
		int get_num_sockets();
		int get_num_threads();
		int get_num_threads_per_socket();

		int get_current_thread_id();
		int get_current_socket_id();
		int get_current_socket_offset();
};

#endif
