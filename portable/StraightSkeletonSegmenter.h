#pragma once
#include <vector>
#include <array>
#include <functional>
using SurfaceDistance = std::function<double(const std::array<double,3>&)>;
class StraightSkeletonSegmenter {
public:
    struct Segment {
        std::vector<int> nodes;
        std::vector<double> points, originalPoints;
        bool valid = true;
        int sourceChainId = -1;
    };
    static std::vector<Segment> ExtractStraightSegments(
        const std::vector<double>& nodes,
        const std::vector<unsigned int>& edges,
        SurfaceDistance* surfaceLocator = nullptr);
};
