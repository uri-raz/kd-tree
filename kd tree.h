#pragma once

#include <array>
#include <cmath>
#include <queue>
#include <vector>
#include <cassert>
#include <numeric>
#include <iostream>
#include <optional>
#include <algorithm>
#include <execution>

#include "kd box.h"

/************************************************************************************************************/
/* The number of points in a leaf is tunable, as a balance needs to be found - a large value will decrease  */
/* the tree height and mispredictions where regionCrossesRegion is called twice, but will increase the time */
/* spent processing leaves, and vice versa.                                                                 */
/* For my testing with k = 3, I've found a leaf size of 98 to be a good default value.                      */
/************************************************************************************************************/

template <unsigned int K, unsigned int leafSize = 98>
class kd_tree
{
public:
	using kd_point = std::array<float, K>;

private:
	static bool isValid(const kd_point &kdPoint) {
		return !std::any_of(kdPoint.cbegin(), kdPoint.cend(), [](float f) { return !std::isfinite(f); });
	}

protected:

	/************************************************************/
	/* There are two types of nodes - an internal one splits    */
	/* space across an axis, and a leaf node holds some points. */
	/************************************************************/

	class kd_node
	{
	public:
		virtual ~kd_node() = default;

		virtual void setExcluded(bool) = 0;

		[[nodiscard]] virtual bool isLeaf() const = 0;
		[[nodiscard]] virtual bool isInternal() const = 0;
		[[nodiscard]] virtual bool isExcluded() const = 0;
		[[nodiscard]] virtual uint64_t TreeHeight() const = 0;
		[[nodiscard]] virtual uint64_t getNumPoints() const = 0;
		[[nodiscard]] virtual const kd_box<K>& boundingBox() const = 0;

		[[nodiscard]] virtual uint64_t nodeCount(bool withInternalNodes) const = 0;

		virtual void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const = 0;

		virtual void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint, kd_box<K> &minRegion) const = 0;

		virtual void KNearestNeighbors(const kd_point &srcPoint, std::vector<std::pair<kd_point, float>>&nearPoints, kd_box<K> &minRegion) const = 0;
	};


	class kd_internal_node final : public kd_node
	{
	private:
		float m_splitVal;
		kd_box<K> m_boundingBox;

	public:
		kd_node *m_Left, *m_Right;
		struct {
			uint64_t m_numPoints : 56;
			size_t m_axis : 7;
			bool m_excluded : 1;
		};

		kd_internal_node(size_t axis, float splitVal, const kd_box<K>& boundingBox, kd_node* Left, kd_node* Right) :
			m_axis(axis), m_splitVal(splitVal), m_boundingBox(boundingBox),
			m_numPoints(Left->getNumPoints() + Right->getNumPoints()),
			m_Left(Left), m_Right(Right) {
			m_excluded = false;
		}

		void setExcluded(bool flag) override { m_excluded = flag; }

		void updateBox(const kd_point &Point) { m_boundingBox.updateBox(Point); }

		[[nodiscard]] size_t splitAxis() const { return m_axis; }
		[[nodiscard]] float splitVal() const { return m_splitVal; }
		[[nodiscard]] bool isLeaf() const override { return false; }
		[[nodiscard]] bool isInternal() const override { return true; }
		[[nodiscard]] bool isExcluded() const override { return m_excluded; }
		[[nodiscard]] const kd_box<K>& boundingBox() const override { return m_boundingBox; }
		[[nodiscard]] uint64_t getNumPoints() const override { return m_numPoints; }
		[[nodiscard]] uint64_t TreeHeight() const override { return 1 + std::max(m_Left->TreeHeight(), m_Right->TreeHeight()); }

		[[nodiscard]] uint64_t nodeCount(bool withInternalNodes) const override {
			return (withInternalNodes ? 1 : 0) + m_Left->nodeCount(withInternalNodes) + m_Right->nodeCount(withInternalNodes);
		}

		void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const override {
			if (regionCrossesRegion(searchBox, m_Left->boundingBox()))
				m_Left->pointsInBox(searchBox, Points);

			if (regionCrossesRegion(searchBox, m_Right->boundingBox()))
				m_Right->pointsInBox(searchBox, Points);
		}

		void KNearestNeighbors(const kd_point &srcPoint,
			                   std::vector<std::pair<kd_point, float>> &nearPoints, kd_box<K> &minRegion) const override
		{
			if (isExcluded())
				return;

			if (regionCrossesRegion(minRegion, m_Left->boundingBox()))
				m_Left->KNearestNeighbors(srcPoint, nearPoints, minRegion);

			if (regionCrossesRegion(minRegion, m_Right->boundingBox()))
				m_Right->KNearestNeighbors(srcPoint, nearPoints, minRegion);
		}

		void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint, kd_box<K> &minRegion) const override {
			if (regionCrossesRegion(minRegion, m_Left->boundingBox()))
				m_Left->nearestNeighbor(srcPoint, nearPoint, minRegion);

			if (regionCrossesRegion(minRegion, m_Right->boundingBox()))
				m_Right->nearestNeighbor(srcPoint, nearPoint, minRegion);
		}
	};


	class kd_leaf_node final : public kd_node
	{
	private:
		std::array<kd_point, leafSize> m_points;	// TODO: replace with inplace_vector
		kd_box<K> m_boundingBox;
		unsigned m_numPoints;
		bool m_excluded = false;

	public:
		kd_leaf_node(typename std::vector<kd_point>::iterator iterBegin, typename std::vector<kd_point>::iterator iterEnd) :
			m_numPoints(unsigned(iterEnd - iterBegin)), m_boundingBox(kd_box<K>(iterBegin, iterEnd))
		{
			std::copy_n(iterBegin, m_numPoints, m_points.begin());
		}

		void setExcluded(bool flag) override { m_excluded = flag; }

		[[nodiscard]] bool isLeaf() const override { return true; }
		[[nodiscard]] bool isInternal() const override { return false; }
		[[nodiscard]] bool isExcluded() const override { return false; }
		[[nodiscard]] uint64_t TreeHeight() const override { return 1; }
		[[nodiscard]] uint64_t nodeCount(bool) const override { return 1; }
		[[nodiscard]] uint64_t getNumPoints() const override { return m_numPoints; }
		[[nodiscard]] const kd_box<K>& boundingBox() const override { return m_boundingBox; }
		[[nodiscard]] const std::array<kd_point, leafSize>& getPoints() { return m_points; }

		[[nodiscard]] void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const override
		{
			for (unsigned i = 0; i < m_numPoints; i++)
				if (searchBox.pointIsInRegion(m_points[i]))
					Points.push_back(m_points[i]);
		}

		virtual void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint, kd_box<K> &minRegion) const override
		{
			float minDistanceSq = DistanceSq(srcPoint, nearPoint);

			for (size_t i = 0; i < m_numPoints; i++) {
				float currDistanceSq = DistanceSq(srcPoint, m_points[i]);

				if (currDistanceSq < minDistanceSq) {
					nearPoint = m_points[i];
					minDistanceSq = currDistanceSq;
					minRegion.set(srcPoint, std::sqrt(currDistanceSq));
				}
			}
		}

		virtual void KNearestNeighbors(const kd_point &srcPoint,
			                           std::vector<std::pair<kd_point, float>> &nearPoints, kd_box<K> &minRegion) const override
		{
			if (m_excluded) [[unlikely]]
				return;

			bool nearestChanged = false;
			size_t k = nearPoints.size();

			for (unsigned i = 0; i < m_numPoints; i++) {
				auto currDistance = DistanceSq(srcPoint, m_points[i]);
				if (currDistance >= nearPoints[k - 1].second * nearPoints[k - 1].second)
					continue;

				nearestChanged = true;

				nearPoints[k - 1].first = m_points[i];
				nearPoints[k - 1].second = std::sqrtf(currDistance);

				for (size_t i = k - 1; i > 0; i--)
					if (nearPoints[i - 1].second > nearPoints[i].second)
						std::swap(nearPoints[i - 1], nearPoints[i]);
					else
						break;
			}

			if (nearestChanged)
				minRegion.set(srcPoint, nearPoints[k - 1].second);
		}
	};

	kd_node* CreateTree(typename std::vector<kd_point>::iterator iterBegin,
		                typename std::vector<kd_point>::iterator iterEnd,
		                unsigned Depth = 0)
	{
		kd_node* retNode = nullptr;
		auto numPoints = iterEnd - iterBegin;

		if (numPoints <= leafSize)
			retNode = new kd_leaf_node(iterBegin, iterEnd);
		else {
			kd_box<K> boundingBox(iterBegin, iterEnd);

			unsigned splitAxis = diffIndex(boundingBox, Depth % K);
			assert(splitAxis != K);
			auto midIter = NthCoordMid(iterBegin, iterEnd, splitAxis);
			float splitValue = (*midIter)[splitAxis];

			retNode = new kd_internal_node(splitAxis, splitValue, boundingBox,
				                           CreateTree(iterBegin, midIter, Depth + 1),
				                           CreateTree(midIter, iterEnd, Depth + 1));
		}

		return retNode;
	}

	[[nodiscard]] std::tuple<kd_point,float> ApproxNearestNeighborPoint(const kd_point& srcPoint) const {
		kd_node* node = m_Root;

		while (node->isInternal()) {
			kd_internal_node* iNode = (kd_internal_node*)node;
			node = (srcPoint[iNode->splitAxis()] <= iNode->splitVal()) ? iNode->m_Left : iNode->m_Right;
		}

		kd_leaf_node* lNode = (kd_leaf_node*)node;

		const auto& points = lNode->getPoints();
		uint64_t numPoints = lNode->getNumPoints();

		kd_point candidatePoint{ points[0] };
		float candidateDistance = Distance(srcPoint, candidatePoint);

		for (unsigned i = 1; i < numPoints; i++) {
			float currDistance = Distance(srcPoint, points[i]);

			if (currDistance < candidateDistance) {
				candidateDistance = currDistance;
				candidatePoint = points[i];
			}
		}

		return { candidatePoint,candidateDistance };
	}

	// The routine has a desired side effect of partitioning Points
	std::vector<kd_point>::iterator NthCoordMid(typename std::vector<kd_point>::iterator iterBegin,
		                                        typename std::vector<kd_point>::iterator iterEnd,
		                                        const unsigned num)
	{
		// This function attempts to split the range iterBegin...iterEnd as close as possible to the middle.
		// As all points with the numth coordinate <= the chosen middle value go to the left son. Therefore,
		// the function has to take care of the corner case in which the points value on the numth axis are
		// identical from the middle on, e.g. Y₀, Y₁, Y₂, Y₂, Y₂, which would move all the points to the
		// left son, leaving the right son empty.

		size_t numPoints = iterEnd - iterBegin;

		std::nth_element(std::execution::par, iterBegin, iterBegin + numPoints / 2, iterEnd,
			             [num](const kd_point& A, const kd_point& B) { return A[num] < B[num]; });

		auto mIter = iterBegin + numPoints / 2;
		float val = (*mIter)[num];

		if ((*std::prev(iterEnd))[num] == val)
			// Decrease iterator until a different value is found.
			while ((*(--mIter))[num] == val);
		else {
			while ((*std::next(mIter))[num] == val)
				++mIter;

			// The implementation needs all points with the numth coordinate <= the chosen middle value to
			// go to the left son. If the value in the nth_element repeats multiple times, nth_element might
			// leave some of the appearances in random places to the right of the nth element.

			for (auto tIter = std::next(mIter); tIter != iterEnd; ++tIter)
				if ((*tIter)[num] == val)
					std::swap(*tIter, *(++mIter));
		}

		return mIter;
	}

	std::vector<std::pair<kd_point, float>> approximateKNearestNeigbors(const kd_point& srcPoint, unsigned k, kd_node*& exclNode) const
	{
		// Go down the tree as if to insert srcPoint until it finds a node with at least k points under it,
		// then pick the k nearest points to srcPoint of those as an initial estimate for the result.

		kd_node* currNode = m_Root;
		while (currNode->isInternal()) {
			kd_internal_node* iNode = (kd_internal_node*)currNode;
			kd_node* nextNode = (srcPoint[iNode->splitAxis()] <= iNode->splitVal()) ? iNode->m_Left : iNode->m_Right;

			if (currNode->getNumPoints() >= k && k > nextNode->getNumPoints()) [[unlikely]]
				break;

			currNode = nextNode;
		}

		exclNode = currNode;
		exclNode->setExcluded(true);

		std::vector<std::pair<kd_point, float>> nearPoints;
		nearPoints.reserve(currNode->getNumPoints());

		std::queue<kd_node*> nodes({ currNode });

		while (!nodes.empty()) {
			currNode = nodes.front();
			nodes.pop();

			if (currNode->isInternal()) {
				kd_internal_node* iNode = (kd_internal_node*)currNode;
				nodes.push(iNode->m_Left);
				nodes.push(iNode->m_Right);
			}
			else {
				kd_leaf_node* lNode = (kd_leaf_node*)currNode;
				const auto& points = lNode->getPoints();
				uint64_t numPoints = lNode->getNumPoints();

				for (unsigned i = 0; i < numPoints; i++)
					nearPoints.emplace_back(points[i], Distance(srcPoint, points[i]));
			}
		}

		if (nearPoints.size() >= k) {
			std::partial_sort(std::execution::par, nearPoints.begin(), nearPoints.begin() + k, nearPoints.end(),
				[](const std::pair<kd_point, float>& lhs, const std::pair<kd_point, float>& rhs) { return lhs.second < rhs.second; });
			nearPoints.resize(k);
		}

		return nearPoints;
	}

	void PrintTree(const kd_node* node, unsigned depth = 0) const {
		for (unsigned i = 0; i < depth; i++) std::cout << ' ';

		if (!node)
			std::cout << "null" << std::endl;
		else {
			if (node->isInternal()) {
				kd_internal_node* iNode = (kd_internal_node*)node;

				std::cout << "Split val is " << iNode->splitVal() << " for axis #" << iNode->splitAxis() << '\n';

				PrintTree(iNode->m_Left, depth + 1);
				PrintTree(iNode->m_Right, depth + 1);
			}
			else {
				kd_leaf_node* lNode = (kd_leaf_node*)node;
				auto numPoints = node->getNumPoints();
				const auto& points = lNode->getPoints();

				std::cout << "Points are\n";
				for (unsigned i = 0; i < numPoints; i++) {
					std::cout << " (";
					for (unsigned j = 0; j < K - 1; j++) {
						std::cout << points[i][j] << ' ';
					}
					std::cout << points[i][K - 1] << ")\n";
				}
			}
		}
	}

	kd_node *m_Root = nullptr;

public:
	kd_tree() = default;

	~kd_tree() {
		if (m_Root == nullptr) [[unlikely]]
			return;

		std::queue<kd_node*> nodes({ m_Root });
		m_Root = nullptr;

		while (!nodes.empty()) {
			kd_node* node = nodes.front();
			nodes.pop();

			if (node->isInternal()) {
				kd_internal_node* iNode = (kd_internal_node*)node;
				nodes.push(iNode->m_Left);
				nodes.push(iNode->m_Right);
			}

			delete node;
		}
	}

	kd_tree(const kd_tree& obj) = delete;
	bool operator=(const kd_tree<K>& rhs) = delete;
	bool operator==(const kd_tree<K> rhs) = delete;

	kd_tree(std::vector<kd_point> &Points) { insert(Points); }

	kd_tree& operator=(kd_tree&& other) noexcept {
		if (this == other) [[unlikely]]
			return this;

		std::swap(m_Root, other.m_Root);
		other->~kd_tree();

		return *this;
	}

	void insert(std::vector<kd_point> &Points) {
		this->~kd_tree();

		// Preamble: remove points with NaNs or infinity,  and duplicate points.
		std::erase_if(Points, [](const kd_point& point) {return !isValid(point); });

		if (Points.empty())
			return;

		std::sort(std::execution::par, Points.begin(), Points.end());
		auto newEnd = std::unique(std::execution::par, Points.begin(), Points.end());
		Points.erase(newEnd, Points.end());
		Points.reserve(Points.size());

		// Build tree with valid & unique points
		m_Root = CreateTree(Points.begin(), Points.end());
	}

	std::vector<kd_point> pointsInBox(const kd_box<K> &searchBox) const {
		std::vector<kd_point> Points;

		if (m_Root && regionCrossesRegion(searchBox, m_Root->boundingBox()))
			m_Root->pointsInBox(searchBox, Points);

		Points.shrink_to_fit();
		return Points;
	}

	[[nodiscard]] std::optional<kd_point> nearestNeighbor(const kd_point &srcPoint) const
	{
		if (!m_Root) [[unlikely]]
			return std::nullopt;

		auto [nearPoint, minDistance] = ApproxNearestNeighborPoint(srcPoint);
		kd_box<K> minBox(srcPoint, minDistance);
		m_Root->nearestNeighbor(srcPoint, nearPoint, minBox);

		return nearPoint;
	}

	std::optional<std::vector<kd_point>> KNearestNeighbors(const kd_point &srcPoint, const unsigned k) const
	{
		if (!m_Root) [[unlikely]]
			return std::nullopt;

		kd_node *exclNode = nullptr;
		auto nearPoints = approximateKNearestNeigbors(srcPoint, k, exclNode);

		if (m_Root->isInternal() && !m_Root->isExcluded() && nearPoints.size() == k) {
			unsigned tk = std::min(k, (unsigned)nearPoints.size());
			kd_box<K> minBox(srcPoint, Distance(srcPoint, nearPoints[tk - 1].first));

			m_Root->KNearestNeighbors(srcPoint, nearPoints, minBox);
		}

		exclNode->setExcluded(false);

		std::vector<kd_point> retVal;
		retVal.reserve(nearPoints.size());
		for (const auto& elem : nearPoints)
			retVal.push_back(elem.first);

		return retVal;
	}

	[[nodiscard]] auto nodeCount(bool withInternalNodes = false) const { return m_Root ? m_Root->nodeCount(withInternalNodes) : 0; }
	[[nodiscard]] auto getNumPoints() const { return m_Root ? m_Root->getNumPoints() : 0; }
	[[nodiscard]] auto TreeHeight() const { return this->m_Root ? m_Root->TreeHeight() : 0; }
	[[nodiscard]] std::optional<kd_box<K>> boundingBox() const { return m_Root ? std::optional(m_Root->boundingBox()) : std::nullopt; }

	void PrintTree() const { PrintTree(m_Root); }
};
