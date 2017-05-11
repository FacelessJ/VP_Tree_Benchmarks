#pragma once

// A VP-Tree implementation, by Steve Hanov. (steve.hanov@gmail.com)
// Released to the Public Domain
// Based on "Data Structures and Algorithms for Nearest Neighbor Search" by Peter N. Yianilos

#include "types.hpp"

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <cmath>
#include <stack>
#include <omp.h>

namespace vp_omp
{


	/**
	* Basic single-core implementation of Vantage Point trees.
	* Will perform knn_search and fixed_radius search.
	* Tree is created from points (point vector is shuffled to accomodate
	* tree). No mapping of new location to original location
	* is kept. Thus this data structure is best used for
	* types where order doesn't matter (i.e indistinguishable particles,
	* where we just need to know the average distance to kth nearest neighbour,
	* etc)
	*/
	template<typename T, double(*distance)(const T&, const T&)>
	class VpTree
	{
	public:

		VpTree() : items(nullptr), rootIdx(-1) { }

		~VpTree() {
		}

		/**
		* Creates a VP tree using data specified in new_items.
		* Tree is created in_place, that is, new_items is mutated,
		* and no record is kept of original indices.
		*/
		void create(std::vector<T>& new_items) {
			items = &new_items;
			nodes.reserve(items->size());
			rootIdx = buildFromPoints(0, items->size());
		}

		/**
		* Find the k nearest points to target and return them and their distances
		*/
		void find_knn(const T& target, const int k, std::vector<T>* results,
					  std::vector<double>* dists)
		{
			std::priority_queue<HeapItem> heap;

			double tau = std::numeric_limits<double>::max();
			find_knn_impl(rootIdx, target, k, heap, tau);

			results->clear(); dists->clear();

			while(!heap.empty()) {
				results->push_back((*items)[heap.top().index]);
				dists->push_back(heap.top().dist);
				heap.pop();
			}

			std::reverse(results->begin(), results->end());
			std::reverse(dists->begin(), dists->end());
		}

		/**
		* Find the kth nearest point to target and return it and its distances
		*/
		void find_kth_neighbour(const T& target, const int k, T& kth_neighbour,
								double& dist) const
		{
			std::priority_queue<HeapItem> heap;

			double tau = std::numeric_limits<double>::max();
			find_knn_impl(rootIdx, target, k, heap, tau);

			kth_neighbour = (*items)[heap.top().index];
			dist = heap.top().dist;
		}

		void batch_find_kth_neighbour(const std::vector<T>& queries, const int k, std::vector<T>& results,
									  std::vector<double>& dists)
		{
			results.resize(queries.size());
			dists.resize(queries.size());

			const auto num = static_cast<int>(queries.size());
#pragma omp parallel for
			for(auto i = 0; i < num; ++i) {
				find_kth_neighbour(queries[i], k, results[i], dists[i]);
			}
		}

		/**
		* Return count of how many points exist within dist of target point
		*/
		size_t fr_count(const T& target, const double dist) const
		{
			return fr_count_impl(rootIdx, target, dist);
		}

		std::vector<size_t> batch_fr_count(const std::vector<T>& queries,
										   const std::vector<double>& dists)
		{
			const auto num = static_cast<int>(queries.size());
			std::vector<size_t> ret(num);

#pragma omp parallel for
				for(auto i = 0; i < num; ++i) {
					ret[i] = fr_count(queries[i], dists[i]);
			}
			return ret;
		}

	private:
		std::vector<T>* items;

		struct Node
		{
			double threshold;
			int left;
			int right;
		};
		size_t rootIdx;

		std::vector<Node> nodes;

		struct HeapItem {
			HeapItem(int index, double dist) :
				index(index), dist(dist) {}
			int index;
			double dist;
			bool operator<(const HeapItem& o) const {
				return dist < o.dist;
			}
		};

		struct DistanceComparator
		{
			const T& item;
			DistanceComparator(const T& item) : item(item) {}
			bool operator()(const T& a, const T& b) {
				return distance(item, a) < distance(item, b);
			}
		};

		size_t buildFromPoints(int lower, int upper)
		{
			static_assert(std::is_pod<Node>::value, "Node isn't POD");
			if(upper == lower) {
				return -1;
			}

			auto nodeId = nodes.size();
			//auto nodeDefaultData = { 0., lower, -1, -1 };
			nodes.push_back(Node{ 0., -1, -1 });

			if(upper - lower > 1) {

				// choose an arbitrary point and move it to the start
				double m = (double)rand() / RAND_MAX;
				m = 0.5;
				int i = (int)(m * (upper - lower - 1)) + lower;
				std::swap((*items)[lower], (*items)[i]);

				int median = (upper + lower) / 2;

				// partitian around the median distance
				std::nth_element(
					items->begin() + lower + 1,
					items->begin() + median,
					items->begin() + upper,
					DistanceComparator((*items)[lower]));

				// what was the median?
				nodes[nodeId].threshold = distance((*items)[lower], (*items)[median]);

				nodes[nodeId].left = buildFromPoints(lower + 1, median);
				nodes[nodeId].right = buildFromPoints(median, upper);
			}

			return nodeId;
		}

		void find_knn_impl(const size_t node, const T& target, const int k,
						   std::priority_queue<HeapItem>& heap, double& tau) const
		{
			if(node == -1) return;

			double dist = distance((*items)[node], target);

			if(dist < tau) {
				if(heap.size() == k) heap.pop();
				heap.push(HeapItem(node, dist));
				if(heap.size() == k) tau = heap.top().dist;
			}

			if(nodes[node].left == -1 && nodes[node].right == -1) {
				return;
			}

			/**
			* Pick the more likely child first
			*/
			if(dist < nodes[node].threshold) {
				if(dist - tau <= nodes[node].threshold) {
					find_knn_impl(nodes[node].left, target, k, heap, tau);
				}

				if(dist + tau >= nodes[node].threshold) {
					find_knn_impl(nodes[node].right, target, k, heap, tau);
				}

			}
			else {
				if(dist + tau >= nodes[node].threshold) {
					find_knn_impl(nodes[node].right, target, k, heap, tau);
				}

				if(dist - tau <= nodes[node].threshold) {
					find_knn_impl(nodes[node].left, target, k, heap, tau);
				}
			}
		}

		size_t fr_count_impl(size_t nodeId, const T& target, const double max_dist) const
		{
			if(nodeId == -1) return 0;

			std::stack<int> workingSet;
			int currNode = -1;
			workingSet.push(nodeId);
			size_t count = 0;

			while(workingSet.empty() == false || currNode != -1) {
				if(currNode != -1) {
					double dist = distance((*items)[currNode], target);

					count += (dist < max_dist);
					if(nodes[currNode].left == -1 && nodes[currNode].right == -1) {
						currNode = -1;
						//continue;
					}
					else {
						// Pick the more likely child first
						if(dist < nodes[currNode].threshold) {
							if((dist - max_dist <= nodes[currNode].threshold))
								workingSet.push(nodes[currNode].left);
							if((dist + max_dist >= nodes[currNode].threshold))
								workingSet.push(nodes[currNode].right);
						}
						else {
							if((dist + max_dist >= nodes[currNode].threshold))
								workingSet.push(nodes[currNode].right);
							if((dist - max_dist <= nodes[currNode].threshold))
								workingSet.push(nodes[currNode].left);
						}
					}
				}
				if(workingSet.empty() == false) {
					currNode = workingSet.top();
					workingSet.pop();
				}
			}

			return count;
		}
	};
}