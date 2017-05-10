#pragma once

#include "types.hpp"
#include <vector>
#include <limits>
#include <algorithm>

namespace brute_force
{
	template<typename Point, double(*Distance)(const Point&, const Point&)>
	void brutedistk(const std::vector<Point>& dataPts, int const point, int const k, Point& result, double& dist) {
		std::priority_queue<std::pair<double, int> > knn;
		knn.push(std::pair<double, int>(std::numeric_limits<double>::max(), -1));

		//int closestIdx = -1;
		const auto numPts = dataPts.size();
		for(size_t i = 0; i < numPts; ++i) {
			double dist = Distance(dataPts[point], dataPts[i]);
			if(dist < knn.top().first) {
				knn.push(std::pair<double, int>(dist, i));
				if((signed int)knn.size() > k)
					knn.pop();
			}
		}
		result = dataPts[knn.top().second];
		dist = knn.top().first;
	}

	template<typename Point, double(*Distance)(const Point&, const Point&)>
	void find_kth_neighbour(const std::vector<Point>& data, const int queryIdx, const int k,
		 Point& result, double& dist)
	{
		int N = (int)data.size();
		if(k == 1) {
			int best_idx = -1;
			double best_dist = std::numeric_limits<double>::max();
			for(int i = 0; i < N; ++i) {
				if(i == queryIdx) continue;
				double dist = Distance(data[queryIdx], data[i]);
				if(dist < best_dist) {
					best_dist = dist;
					best_idx = i;
				}
			}
			result = data[best_idx];
			dist = best_dist;
		}
		else {
			for(int i = 0; i < N; ++i)
				brutedistk<Point, Distance>(data, queryIdx, k + 1, result, dist);
		}
	}

	template<typename Point, double(*Distance)(const Point&, const Point&)>
	size_t fr_count(const std::vector<Point>& data, const int queryIdx, const double dist)
	{
		return std::count_if(data.begin(), data.end(),
							 [&data, queryIdx, &dist](auto& i) { return Distance(data[queryIdx], i) < dist; } );
	}

}
