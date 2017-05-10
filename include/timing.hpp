#pragma once

#include <chrono>
#include <string>

class ScopeTimer
{
	using Clock = std::chrono::steady_clock;
	using TimePoint = std::chrono::time_point<Clock>;
public:
	ScopeTimer(const std::string& name) : start(Clock::now()), name(name) {};
	~ScopeTimer();
private:
	TimePoint start;
	std::string name;
};