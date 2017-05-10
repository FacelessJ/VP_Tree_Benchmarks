#include <iostream>

#include "../include/basic_vp_tree.hpp"
#include "../include/brute_force.hpp"
#include "../include/types.hpp"
#include "../include/timing.hpp"

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
		{
			ScopeTimer t("Create");
			basicTree.create(data);
		}

		brute_force::fr_count<Point, vp_detail::MaxNormDist<Point>>(data, 0, 2.0);

		// find kth neighbour for each point
		Point kth_neighbour;
		double dist_to_kth_neighbour[NUM_ELEMENTS];

		Point kth_neighbour_bf;
		double dist_to_kth_neighbour_bf;
		{
			ScopeTimer t("Dists");
			for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
				basicTree.find_kth_neighbour(data[i], k, kth_neighbour, dist_to_kth_neighbour[i]);
				//brute_force::find_kth_neighbour<Point, vp_detail::MaxNormDist<Point>>(data, i, k-1, kth_neighbour_bf, dist_to_kth_neighbour_bf);
				// do something with results
				//if(vp_detail::EuclideanDist<Point>(kth_neighbour, kth_neighbour_bf) > 0.000001)
				//{
				//	std::cout << "Diff for " << i << "\n";
				//}
			}
		}

		{
			ScopeTimer t("Count");
			for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
				auto count1 = basicTree.fr_count(data[i], dist_to_kth_neighbour[i] + 1);
				/*auto count2 = brute_force::fr_count<Point, vp_detail::MaxNormDist<Point>>(data, i, dist_to_kth_neighbour[i]+1);
				if(count1 != count2) {
					std::cout << "Count diff for " << i << "\n";
				}*/
			}
		}
		int a = 5;
	}

	return 0;
}