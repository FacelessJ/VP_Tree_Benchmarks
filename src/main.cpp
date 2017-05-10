#include <iostream>

#include "../include/basic_vp_tree.hpp"
#include "../include/brute_force.hpp"
#include "../include/types.hpp"

#include <vector>
#include <random>

int main()
{
	std::cout << "Hello, world!" << std::endl;

	constexpr size_t NUM_ELEMENTS = 1000;
	constexpr int k = 6;
	{
		using Point = vp_detail::Point;
		using VPTree = vp_basic::VpTree<Point, vp_detail::MaxNormDist<Point>>;

		std::uniform_real_distribution<double> distribution(0, 100);
		std::default_random_engine generator;

		// Create random point cloud
		std::vector<Point> data;
		data.reserve(NUM_ELEMENTS);
		std::generate_n(std::back_inserter(data), NUM_ELEMENTS, [distribution, &generator]() { 
			Point tmp;
			for(int i = 0; i < vp_detail::DIM; ++i)
				tmp.coords[i] = distribution(generator);
			return tmp; });

		
		// Create tree
		VPTree basicTree;
		basicTree.create(data);

		brute_force::fr_count<Point, vp_detail::MaxNormDist<Point>>(data, 0, 2.0);

		// find all knn for each point
		std::vector<Point> results;
		std::vector<double> dists;
		for(auto& pt : data) {
			basicTree.find_knn(pt, k, &results, &dists);
			// do something with results
		}

		// find kth neighbour for each point
		Point kth_neighbour;
		double dist_to_kth_neighbour;
		for(auto& pt : data) {
			basicTree.find_kth_neighbour(pt, k, kth_neighbour, dist_to_kth_neighbour);
			// do something with results
			std::cout << "Distance: " << dist_to_kth_neighbour << "\n";
		}
	}

	return 0;
}