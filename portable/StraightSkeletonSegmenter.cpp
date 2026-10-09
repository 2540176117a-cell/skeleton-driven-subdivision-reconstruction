#include "StraightSkeletonSegmenter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <utility>


namespace
{
using Point = std::array<double, 3>;
using EdgeKey = std::pair<int, int>;

EdgeKey MakeEdgeKey(int a, int b)
{
	if (a > b) std::swap(a, b);
	return EdgeKey(a, b);
}

Point PointAt(const std::vector<double>& nodes, int index)
{
	return {
		nodes[3 * index],
		nodes[3 * index + 1],
		nodes[3 * index + 2]
	};
}

Point Add(const Point& a, const Point& b)
{
	return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

Point Subtract(const Point& a, const Point& b)
{
	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

Point Scale(const Point& v, double s)
{
	return { v[0] * s, v[1] * s, v[2] * s };
}

double Dot(const Point& a, const Point& b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Length(const Point& v)
{
	return std::sqrt(Dot(v, v));
}

double Distance(const Point& a, const Point& b)
{
	return Length(Subtract(a, b));
}

Point Lerp(const Point& a, const Point& b, double t)
{
	return Add(a, Scale(Subtract(b, a), t));
}

double Clamp01(double value)
{
	if (value < 0.0) return 0.0;
	if (value > 1.0) return 1.0;
	return value;
}

double ComputeBoundingBoxDiagonal(const std::vector<double>& nodes)
{
	if (nodes.size() < 3) {
		return 0.0;
	}

	double minP[3] = { nodes[0], nodes[1], nodes[2] };
	double maxP[3] = { nodes[0], nodes[1], nodes[2] };

	for (size_t i = 0; i < nodes.size() / 3; ++i) {
		for (int axis = 0; axis < 3; ++axis) {
			double value = nodes[3 * i + axis];
			minP[axis] = std::min(minP[axis], value);
			maxP[axis] = std::max(maxP[axis], value);
		}
	}

	double dx = maxP[0] - minP[0];
	double dy = maxP[1] - minP[1];
	double dz = maxP[2] - minP[2];
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double PolylineLength(const std::vector<Point>& polyline, int begin, int end)
{
	double length = 0.0;
	for (int i = begin + 1; i <= end; ++i) {
		length += Distance(polyline[i - 1], polyline[i]);
	}
	return length;
}

double DistanceToSurface(SurfaceDistance* locator, const Point& point)
{
    return locator && *locator ? (*locator)(point) : 0.0;
}

std::vector<Point> BuildPolyline(
	const std::vector<double>& nodes,
	const std::vector<int>& chain)
{
	std::vector<Point> polyline;
	polyline.reserve(chain.size());

	for (int node : chain) {
		polyline.push_back(PointAt(nodes, node));
	}

	return polyline;
}

Point PointAtDistance(const std::vector<Point>& polyline, double targetDistance)
{
	if (polyline.empty()) {
		return { 0.0, 0.0, 0.0 };
	}

	if (targetDistance <= 0.0) {
		return polyline.front();
	}

	double walked = 0.0;
	for (size_t i = 1; i < polyline.size(); ++i) {
		double edgeLength = Distance(polyline[i - 1], polyline[i]);
		if (edgeLength <= 1e-8) {
			continue;
		}

		if (walked + edgeLength >= targetDistance) {
			double t = (targetDistance - walked) / edgeLength;
			return Lerp(polyline[i - 1], polyline[i], Clamp01(t));
		}

		walked += edgeLength;
	}

	return polyline.back();
}

std::vector<Point> TrimPolyline(
	const std::vector<Point>& polyline,
	double trimStart,
	double trimEnd)
{
	std::vector<Point> trimmed;
	if (polyline.size() < 2) {
		return trimmed;
	}

	double totalLength = PolylineLength(polyline, 0, static_cast<int>(polyline.size()) - 1);
	double startDistance = std::min(trimStart, totalLength);
	double endDistance = std::max(startDistance, totalLength - trimEnd);

	if (endDistance - startDistance <= 1e-8) {
		return trimmed;
	}

	trimmed.push_back(PointAtDistance(polyline, startDistance));

	double walked = 0.0;
	for (size_t i = 1; i + 1 < polyline.size(); ++i) {
		walked += Distance(polyline[i - 1], polyline[i]);
		if (walked > startDistance && walked < endDistance) {
			trimmed.push_back(polyline[i]);
		}
	}

	trimmed.push_back(PointAtDistance(polyline, endDistance));
	return trimmed;
}

bool IsPolylineStraight(
	const std::vector<Point>& polyline,
	int begin,
	int end,
	double modelDiagonal)
{
	if (end - begin < 2) {
		return false;
	}

	const Point start = polyline[begin];
	const Point finish = polyline[end];
	const Point chord = Subtract(finish, start);
	const double chordLength = Length(chord);

	if (chordLength <= 1e-8) {
		return false;
	}

	const double arcLength = PolylineLength(polyline, begin, end);
	if (arcLength <= 1e-8 || arcLength / chordLength > 1.04) {
		return false;
	}

	const double invChordLength2 = 1.0 / (chordLength * chordLength);
	const double maxDeviation = std::max(chordLength * 0.035, modelDiagonal * 0.003);
	const double minDirectionCos = std::cos(12.0 * 3.14159265358979323846 / 180.0);

	for (int i = begin + 1; i < end; ++i) {
		const Point toPoint = Subtract(polyline[i], start);
		const double t = Dot(toPoint, chord) * invChordLength2;
		const Point projected = Add(start, Scale(chord, t));

		if (Distance(polyline[i], projected) > maxDeviation) {
			return false;
		}
	}

	for (int i = begin + 1; i <= end; ++i) {
		const Point localDirection = Subtract(polyline[i], polyline[i - 1]);
		const double localLength = Length(localDirection);

		if (localLength <= 1e-8) {
			continue;
		}

		const double cosValue = std::abs(Dot(localDirection, chord) / (localLength * chordLength));
		if (cosValue < minDirectionCos) {
			return false;
		}
	}

	return true;
}

void AppendStraightRuns(
	const std::vector<Point>& polyline,
	double modelDiagonal,
	int sourceChainId,
	std::vector<StraightSkeletonSegmenter::Segment>& output)
{
	if (polyline.size() < 2) {
		return;
	}

	const double minRunLength = std::max(modelDiagonal * 0.015, 1e-6);

	if (polyline.size() == 2) {
		if (Distance(polyline[0], polyline[1]) >= minRunLength) {
			StraightSkeletonSegmenter::Segment segment;
			segment.valid = true;
			segment.sourceChainId = sourceChainId;

			segment.points.push_back(polyline[0][0]);
			segment.points.push_back(polyline[0][1]);
			segment.points.push_back(polyline[0][2]);
			segment.points.push_back(polyline[1][0]);
			segment.points.push_back(polyline[1][1]);
			segment.points.push_back(polyline[1][2]);

			segment.originalPoints = segment.points;
			output.push_back(segment);
		}
		return;
	}

	const int pointCount = static_cast<int>(polyline.size());
	int start = 0;

	while (start + 2 < pointCount) {
		int bestEnd = -1;

		for (int end = start + 2; end < pointCount; ++end) {
			if (IsPolylineStraight(polyline, start, end, modelDiagonal)) {
				bestEnd = end;
			}
			else if (bestEnd >= 0) {
				break;
			}
		}

		if (bestEnd >= 0 &&
			PolylineLength(polyline, start, bestEnd) >= minRunLength) {
			StraightSkeletonSegmenter::Segment segment;
			segment.valid = true;
			segment.sourceChainId = sourceChainId;

			for (int i = start; i <= bestEnd; ++i) {
				segment.points.push_back(polyline[i][0]);
				segment.points.push_back(polyline[i][1]);
				segment.points.push_back(polyline[i][2]);
			}

			segment.originalPoints = segment.points;
			output.push_back(segment);

			start = bestEnd;
		}
		else {
			++start;
		}
	}
}
}

std::vector<StraightSkeletonSegmenter::Segment>
StraightSkeletonSegmenter::ExtractStraightSegments(
	const std::vector<double>& nodes,
	const std::vector<unsigned int>& edges,
	SurfaceDistance* surfaceLocator)
{
	std::vector<Segment> result;
	const int nodeCount = static_cast<int>(nodes.size() / 3);

	if (nodeCount <= 1 || edges.size() < 2) {
		return result;
	}

	std::vector<std::vector<int>> adjacency(nodeCount);

	for (size_t i = 0; i + 1 < edges.size(); i += 2) {
		int a = static_cast<int>(edges[i]);
		int b = static_cast<int>(edges[i + 1]);

		if (a < 0 || b < 0 || a >= nodeCount || b >= nodeCount || a == b) {
			continue;
		}

		adjacency[a].push_back(b);
		adjacency[b].push_back(a);
	}

	std::set<EdgeKey> visitedEdges;
	std::vector<std::vector<int>> branchChains;

	auto traceFrom = [&](int start, int next) {
		std::vector<int> chain;
		chain.push_back(start);
		chain.push_back(next);

		visitedEdges.insert(MakeEdgeKey(start, next));

		int prev = start;
		int curr = next;

		while (curr >= 0 &&
			curr < nodeCount &&
			adjacency[curr].size() == 2) {
			int candidate = adjacency[curr][0] == prev
				? adjacency[curr][1]
				: adjacency[curr][0];

			EdgeKey key = MakeEdgeKey(curr, candidate);
			if (visitedEdges.find(key) != visitedEdges.end()) {
				break;
			}

			chain.push_back(candidate);
			visitedEdges.insert(key);

			prev = curr;
			curr = candidate;
		}

		if (chain.size() >= 2) {
			branchChains.push_back(chain);
		}
	};

	for (int node = 0; node < nodeCount; ++node) {
		if (adjacency[node].size() == 2) {
			continue;
		}

		for (int next : adjacency[node]) {
			EdgeKey key = MakeEdgeKey(node, next);
			if (visitedEdges.find(key) != visitedEdges.end()) {
				continue;
			}

			traceFrom(node, next);
		}
	}

	const double modelDiagonal = ComputeBoundingBoxDiagonal(nodes);
	const double fallbackBlendDistance = modelDiagonal * 0.035;
	const double blendFactor = 2.5;

	for (size_t chainIndex = 0; chainIndex < branchChains.size(); ++chainIndex) {
		const std::vector<int>& chain = branchChains[chainIndex];

		if (chain.size() < 2) {
			continue;
		}

		std::vector<Point> polyline = BuildPolyline(nodes, chain);
		if (polyline.size() < 2) {
			continue;
		}

		const double totalLength = PolylineLength(
			polyline,
			0,
			static_cast<int>(polyline.size()) - 1);

		if (totalLength <= 1e-8) {
			continue;
		}

		int firstNode = chain.front();
		int lastNode = chain.back();

		const bool topologicallyClosedLoop =
			chain.size() > 2 &&
			firstNode == lastNode;

		const double endDistance = Distance(polyline.front(), polyline.back());
		const bool geometricallyClosedLoop =
			chain.size() > 4 &&
			endDistance < std::max(modelDiagonal * 0.006, totalLength * 0.08) &&
			totalLength > modelDiagonal * 0.02;

		if (topologicallyClosedLoop || geometricallyClosedLoop) {
			continue;
		}

		double trimStart = 0.0;
		double trimEnd = 0.0;

		if (adjacency[firstNode].size() > 2) {
			double radius = DistanceToSurface(surfaceLocator, PointAt(nodes, firstNode));
			if (radius <= 1e-8) {
				radius = fallbackBlendDistance;
			}
			trimStart = radius * blendFactor;
		}

		if (adjacency[lastNode].size() > 2) {
			double radius = DistanceToSurface(surfaceLocator, PointAt(nodes, lastNode));
			if (radius <= 1e-8) {
				radius = fallbackBlendDistance;
			}
			trimEnd = radius * blendFactor;
		}

		const bool bothEndsAreJunctions =
			adjacency[firstNode].size() > 2 &&
			adjacency[lastNode].size() > 2;

		const double maxTrimPerSide = totalLength * (bothEndsAreJunctions ? 0.25 : 0.35);

		trimStart = std::min(trimStart, maxTrimPerSide);
		trimEnd = std::min(trimEnd, maxTrimPerSide);

		std::vector<Point> trimmed = TrimPolyline(polyline, trimStart, trimEnd);

		AppendStraightRuns(
			trimmed,
			modelDiagonal,
			static_cast<int>(chainIndex),
			result);
	}

	return result;
}
