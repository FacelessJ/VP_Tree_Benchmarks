#pragma once

#include <cmath>

namespace vp_detail
{
	constexpr int DIM = 2;

	struct VPNode
	{
		double threshold;
		int index;
		int left, right;
		//int padding; // Not needed. Can be replaced

		VPNode() : threshold(0.), index(0), left(-1), right(-1) {}
	};

	struct Point
	{
		double coords[DIM];
	};

	template <typename T>
	double EuclideanDist(const T &a, const T &b)
	{
		double total = 0.;
		for(size_t i = 0; i < DIM; ++i) {
			total += (b.coords[i] - a.coords[i]) * (b.coords[i] - a.coords[i]);
		}
		return std::sqrt(total);
	}

	template <typename T>
	double MaxNormDist(const T &a, const T &b)
	{
		double maxVal = 0.;
		for(size_t i = 0; i < DIM; ++i) {
			maxVal = std::max(maxVal, std::abs(b.coords[i] - a.coords[i]));
		}
		return maxVal;
	}
}
