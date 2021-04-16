#include <stdio.h>

#include "Random.hpp"
#include "Debug.hpp"

int main(int argc, char ** argv) {
	//NetworkInterface * net_instance = NetworkInterface::get_instance(&argc, &argv);
	NetworkInterface::init_network_interface(&argc, &argv);

	Debug::get_instance()->enter_function("main");
	Debug::get_instance()->log("going to generate 10 random int between 0 and 9...");
	Random random;
	for (int i = 0; i < 10; ++ i) {
		Debug::get_instance()->log("rand_int[", i, "]: ", random.rand_int(0, 9));
	}
	Debug::get_instance()->leave_function("main");

	NetworkInterface::finalize_instance();
	return 0;
}


