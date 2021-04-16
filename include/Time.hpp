#ifndef TIME_HPP
#define TIME_HPP

#include <string>
#include <utility>
#include <string>
#include <vector>

double get_time();

class Timer {
	private:
		static std::vector<std::pair<std::string, double>> timers;
	public:
		static void timer_start(std::string timer_name);
		static void timer_stop(std::string timer_name);
		static void report_timers();
};

#endif
