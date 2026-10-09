#pragma once

#include <vector>

#include <vtkActor.h>
#include <vtkCellLocator.h>
#include <vtkSmartPointer.h>

class StraightSkeletonSegmenter
{
public:
	struct Segment
	{
		std::vector<int> nodes;

		// 当前有效的直线段点，用于显示和后续建笼
		std::vector<double> points;

		// 初始识别结果，用于后续需要恢复时使用
		std::vector<double> originalPoints;

		// 是否保留该直线段
		bool valid = true;

		// 该直线段来源于哪一条骨架 branch chain。
		// 后续自动控制笼共享端面时，只允许 sourceChainId 相同的 cage 互相共享。
		int sourceChainId = -1;
	};

	static std::vector<Segment> ExtractStraightSegments(
		const std::vector<double>& nodes,
		const std::vector<unsigned int>& edges,
		vtkCellLocator* surfaceLocator = nullptr);

	static vtkSmartPointer<vtkActor> CreateSegmentActor(
		const std::vector<double>& nodes,
		const std::vector<Segment>& segments);
};