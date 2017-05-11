#include <iostream>

#include "../include/basic_vp_tree.hpp"
#include "../include/omp_vp_tree.hpp"
#include "../include/brute_force.hpp"
#include "../include/types.hpp"
#include "../include/timing.hpp"

#include <vector>
#include <random>

int main()
{
	std::cout << "Hello, world!" << std::endl;

	constexpr size_t NUM_ELEMENTS = 100000;
	constexpr int k = 6;
	{
		using Point = vp_detail::Point;
		using VPTree = vp_omp::VpTree<Point, vp_detail::MaxNormDist<Point>>;

		std::uniform_real_distribution<double> distribution(0, 100);
		std::default_random_engine generator;

		// Create random point cloud
		std::vector<Point> data;
		data.reserve(NUM_ELEMENTS);
		{
			//std::cout << "Generating" << std::endl;
			//ScopeTimer t("Generate");
			std::generate_n(std::back_inserter(data), NUM_ELEMENTS, [&distribution, &generator]() {
				Point tmp;
				for(int i = 0; i < vp_detail::DIM; ++i)
					tmp.coords[i] = distribution(generator);
				return tmp; });
		}
		
		// Create tree
		VPTree basicTree;
		{
			//std::cout << "Creating" << std::endl;
			ScopeTimer t("Create");
			basicTree.create(data);
		}

		brute_force::fr_count<Point, vp_detail::MaxNormDist<Point>>(data, 0, 2.0);

		// find kth neighbour for each point
		Point kth_neighbour;
		std::vector<double> dist_to_kth_neighbour(NUM_ELEMENTS);

		{
			ScopeTimer t("Dists Recursive");
			for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
				basicTree.find_kth_neighbour(data[i], k, kth_neighbour, dist_to_kth_neighbour[i]);
			}
		}

		std::vector<Point> kth_neighbours2(NUM_ELEMENTS);
		std::vector<double> dist_to_kth_neighbour2(NUM_ELEMENTS);

		{
			ScopeTimer t("Dists Batch");
			basicTree.batch_find_kth_neighbour(data, k, kth_neighbours2, dist_to_kth_neighbour2);
		}

		for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
			if(std::abs(dist_to_kth_neighbour[i] - dist_to_kth_neighbour2[i]) > 0.00001) {
				std::cout << "Diff\n";
			}
		}
		

		{
			ScopeTimer t("Count");
			for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
				auto count1 = basicTree.fr_count(data[i], dist_to_kth_neighbour[i]);
			}
		}
		{
			ScopeTimer t("Count Batch");
			auto count2 = basicTree.batch_fr_count(data, dist_to_kth_neighbour);
		}
	}
	system("pause");
	return 0;
}
