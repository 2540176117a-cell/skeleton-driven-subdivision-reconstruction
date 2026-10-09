#pragma once

#include <vector>
#include <string>

class vtkPolyData;

class SkeletonExtractor 
{
public:
	// 定义骨架数据结构，必须放在类内部
	struct SkeletonData 
	{
		std::vector<double> nodes;           // 坐标 [x1, y1, z1, x2, y2, z2...]
		std::vector<unsigned int> edges;     // 边索引 [u1, v1, u2, v2...]
		std::vector<double> radii;           // 每个骨架节点对应的局部背景半径
		std::vector<unsigned int> supportOffsets; // 每个骨架节点支撑点在 supportPoints 中的起始偏移
		std::vector<double> supportPoints;   // CGAL 骨架节点对应的原始背景点 [x,y,z...]
	};

	SkeletonExtractor();
	~SkeletonExtractor();

	// 保留之前的各向同性重构函数
	bool remeshIsotropic(const std::string& input_stl_file,
		double target_edge_length,
		std::vector<double>& out_vertices,
		std::vector<unsigned int>& out_faces);

	// 新增的骨架提取函数
	bool extractSkeleton(const std::string& input_path,
		std::vector<double>& mesh_v,
		std::vector<unsigned int>& mesh_f,
		SkeletonData& skel_out);

	// 根据封闭 non-design 表面剔除内部骨架点，并同步重映射边、半径和支撑点。
	static SkeletonData filterSkeletonInsideClosedSurface(
		const SkeletonData& input,
		vtkPolyData* closedSurface,
		double tolerance = 1e-6);

};
