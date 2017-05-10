#include "../include/timing.hpp"

#include <iostream>

ScopeTimer::~ScopeTimer()
{
	auto diff = Clock::now() - start;
	std::cout << name <<  " timer alive for " << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count() << " ms\n";
}
