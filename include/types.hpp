#pragma once


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
}