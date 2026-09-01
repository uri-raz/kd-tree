#pragma once

#include <array>
#include <cassert>

template <unsigned int K>
class kd_box
{
public:
	using kd_point = std::array<float, K>;
	kd_point first, second;

	kd_box() {
		first.fill(std::numeric_limits<float>::quiet_NaN());
		second.fill(std::numeric_limits<float>::quiet_NaN());
	}

	kd_box(const kd_point& Point1, const kd_point& Point2) { set(Point1, Point2); }

	void set(const kd_point& Point1, const kd_point& Point2) {
		for (unsigned i = 0; i < K; i++) {
			first[i] = std::min(Point1[i], Point2[i]);
			second[i] = std::max(Point1[i], Point2[i]);
		}
	}

	kd_box(const kd_box& lhs, const kd_box& rhs) {
		for (unsigned i = 0; i < K; i++) {
			first[i] = std::min(lhs.first[i], rhs.first[i]);
			second[i] = std::max(lhs.second[i], rhs.second[i]);
		}
	}

	kd_box(const kd_point& point, const float distance) { set(point, distance); }

	void set(const kd_point& point, const float distance) {
		for (unsigned i = 0; i < K; i++) {
			first[i] = point[i] - distance;
			second[i] = point[i] + distance;
		}
	}

	kd_box(const typename std::vector<kd_point>::iterator& iterBegin,
		   const typename std::vector<kd_point>::iterator& iterEnd)
	{
		assert(iterBegin != iterEnd);
		
		first = second = *iterBegin;
		for (auto currPoint = std::next(iterBegin); currPoint != iterEnd; currPoint++) {
			for (unsigned i = 0; i < K; i++) {
				first[i] = std::min(first[i], (*currPoint)[i]);
				second[i] = std::max(second[i], (*currPoint)[i]);
			}
		}
	}

	void updateBox(const kd_point& Point) {
		for (unsigned i = 0; i < K; i++) {
			first[i] = std::min(first[i], Point[i]);
			second[i] = std::max(second[i], Point[i]);
		}
	}

	bool pointIsInRegion(const kd_point& Point) const {
		for (unsigned i = 0; i < K; i++)
			if (first[i] > Point[i] || Point[i] > second[i])
				return false;

		return true;
	}
};

template <unsigned int K>
inline bool operator==(const kd_box<K>& lhs, const kd_box<K>& rhs) {
	return (lhs.first == rhs.first &&
		    lhs.second == rhs.second);
}

template <unsigned int K>
unsigned diffIndex(const std::array<float, K>& P, const std::array<float, K>& Q, unsigned s = 0) {
	for (unsigned i = 0; i < K; ++i, ++s %= K)
		if (P[s] != Q[s])
			return s;

	return K;
}

template <unsigned int K>
unsigned diffIndex(const kd_box<K>& box, unsigned s = 0) { return diffIndex(box.first, box.second, s); }

template <unsigned K>
bool regionCrossesRegion(const kd_box<K>& Region1, const kd_box<K>& Region2) {
	for (unsigned i = 0; i < K; i++)
		if (Region1.first[i] > Region2.second[i] || Region2.first[i] > Region1.second[i])
			return false;

	return true;
}

template <unsigned K>
float DistanceSq(const std::array<float,K>& P, const std::array<float, K>& Q) {
	float retVal = 0;
	for (unsigned i = 0; i < K; i++)
		retVal += (P[i] - Q[i]) * (P[i] - Q[i]);

	return retVal;
}

template <unsigned int K>
float Distance(const std::array<float, K>& P, const std::array<float, K>& Q) {
	return std::sqrtf(DistanceSq(P, Q));
}
