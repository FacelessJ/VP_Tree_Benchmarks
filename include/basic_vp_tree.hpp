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

namespace vp_basic
{
	template <typename T>
	double VDistance(const T &a, const T &b)
	{
		double total = 0.;
		for(size_t i = 0; i < 2; ++i) {
			total += (b.coords[i] - a.coords[i]) * (b.coords[i] - a.coords[i]);
		}
		return sqrt(total);
	}

	template <typename T>
	double VMaxNorm(const T &a, const T &b)
	{
		double maxVal = 0.;
		for(size_t i = 0; i < vp_detail::DIM; ++i) {
			maxVal = std::max(maxVal, std::abs(b.coords[i] - a.coords[i]));
		}
		return maxVal;
	}

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
		
		VpTree() : items(nullptr), tau(std::numeric_limits<double>::max()), rootIdx(-1) { }

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

		void find_knn(const T& target, int k, std::vector<T>* results,
					std::vector<double>* distances)
		{
			std::priority_queue<HeapItem> heap;

			tau = std::numeric_limits<double>::max();
			find_knn_impl(rootIdx, target, k, heap);

			results->clear(); distances->clear();

			while(!heap.empty()) {
				results->push_back((*items)[heap.top().index]);
				distances->push_back(heap.top().dist);
				heap.pop();
			}

			std::reverse(results->begin(), results->end());
			std::reverse(distances->begin(), distances->end());
		}

		void find_kth_neighbour(const T& target, int k, T& kth_neighbour,
								double& distance)
		{
			std::priority_queue<HeapItem> heap;

			tau = std::numeric_limits<double>::max();
			find_knn_impl(rootIdx, target, k, heap);

			kth_neighbour = (*items)[heap.top().index];
			distance = heap.top().dist;
		}

	private:
		std::vector<T>* items;
		double tau;

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

		void find_knn_impl(size_t node, const T& target, int k,
					std::priority_queue<HeapItem>& heap)
		{
			if(node == -1) return;

			double dist = distance((*items)[node], target);
			//printf("dist=%g tau=%gn", dist, tau );

			if(dist < tau) {
				if(heap.size() == k) heap.pop();
				heap.push(HeapItem(node, dist));
				if(heap.size() == k) tau = heap.top().dist;
			}

			if(nodes[node].left == -1 && nodes[node].right == -1) {
				return;
			}

			if(dist < nodes[node].threshold) {
				if(dist - tau <= nodes[node].threshold) {
					find_knn_impl(nodes[node].left, target, k, heap);
				}

				if(dist + tau >= nodes[node].threshold) {
					find_knn_impl(nodes[node].right, target, k, heap);
				}

			}
			else {
				if(dist + tau >= nodes[node].threshold) {
					find_knn_impl(nodes[node].right, target, k, heap);
				}

				if(dist - tau <= nodes[node].threshold) {
					find_knn_impl(nodes[node].left, target, k, heap);
				}
			}
		}
	};
}