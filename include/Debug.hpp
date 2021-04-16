#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <iostream>
#include <string>

#include "NetworkInterface.hpp"

//#define DEBUG_MODE

class Debug {
	private:
		static Debug * instance;

		Debug() {
		}

		template<typename T>
			void rec_log(T t) {
				std::cout << t << std::endl;
			}
		template<typename T, typename... Args>
			void rec_log(T t, Args... args) {
				std::cout << t;
				rec_log(args...);
			}
	public:
		static Debug * get_instance();

		template<typename... Args>
			void log(Args... args) {
#ifdef DEBUG_MODE
				std::cout << "[" << NetworkInterface::get_instance(nullptr, nullptr)->get_partition_id() << "] ";
				rec_log(args...);
#endif
			}
		template<typename ...Args>
			void print(Args... args) {
				std::cout << "[" << NetworkInterface::get_instance(nullptr, nullptr)->get_partition_id() << "] ";
				rec_log(args...);
			}
		void enter_function(std::string func_name);
		void leave_function(std::string func_name);
};

#endif
