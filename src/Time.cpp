#include <assert.h>

#include <sys/time.h>

#include "Time.hpp"
#include "NetworkInterface.hpp"
#include "Debug.hpp"

std::vector<std::pair<std::string, double>> Timer::timers;

double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1e6);
}

void Timer::timer_start(std::string timer_name) {
	for (int i = 0, j = timers.size(); i < j; ++ i) {
		if (timers[i].first == timer_name) {
			timers[i].second -= get_time();
			return ;
		}
	}
	timers.push_back(make_pair(timer_name, -get_time()));
}

void Timer::timer_stop(std::string timer_name) {
	for (int i = 0, j = timers.size(); i < j; ++ i) {
		if (timers[i].first == timer_name) {
			timers[i].second += get_time();
			return ;
		}
	}
	assert(false);
}

void Timer::report_timers() {
	if (NetworkInterface::get_instance()->get_partition_id() == 0) {
		Debug::get_instance()->print("timer_name,		time(ms)");
		for (std::pair<std::string, double> i: timers) {
			Debug::get_instance()->print(i.first, ",		", int(i.second * 1000));
		}
	}
}
