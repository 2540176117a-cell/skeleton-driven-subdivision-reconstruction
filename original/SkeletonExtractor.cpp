#ifndef NOMINMAX
#define NOMINMAX
#endif

// 保持这个头文件顺序，防止环境冲突
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
// --- 新增：骨架提取核心头文件 ---
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>

// --- 新增：Boost 图形处理头文件 (解决 vertices/edges 报错) ---
#include <boost/foreach.hpp>
#include <boost/graph/graph_traits.hpp>

#include "SkeletonExtractor.h"
#include <QDebug>
#include <QString>
#include <algorithm>
#include <map>

#include <set>
#include <cmath>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkMath.h>
#include <vtkSelectEnclosedPoints.h>
#include <QDebug>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                         TriangleMesh;
typedef boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

static double percentileValue(std::vector<double> values, double percentile)
{
	if (values.empty()) {
		return 0.0;
	}

	std::sort(values.begin(), values.end());

	if (percentile < 0.0) percentile = 0.0;
	if (percentile > 1.0) percentile = 1.0;

	const double rawIndex = percentile * static_cast<double>(values.size() - 1);
	const size_t lo = static_cast<size_t>(std::floor(rawIndex));
	const size_t hi = std::min(lo + 1, values.size() - 1);
	const double t = rawIndex - static_cast<double>(lo);

	return values[lo] * (1.0 - t) + values[hi] * t;
}

// 必须包含这两个空函数，否则链接器找不到它们
SkeletonExtractor::SkeletonExtractor() 
{
	// 这里可以为空，但必须有定义
}

SkeletonExtractor::~SkeletonExtractor()
{
	// 这里可以为空，但必须有定义
}
bool SkeletonExtractor::remeshIsotropic(const std::string& input_stl_file,
	double target_edge_length,
	std::vector<double>& out_vertices,
	std::vector<unsigned int>& out_faces)
{
	TriangleMesh mesh;

	// 1. 读取 STL
	if (!PMP::IO::read_polygon_mesh(input_stl_file, mesh) || mesh.is_empty()) {
		qDebug() << "Error: Cannot read STL file.";
		return false;
	}

	// 2. 自动计算边长 (如果传入 0)
	if (target_edge_length <= 0) {
		auto bbox = CGAL::bbox_3(mesh.points().begin(), mesh.points().end());
		double dx = bbox.xmax() - bbox.xmin();
		double dy = bbox.ymax() - bbox.ymin();
		target_edge_length = std::sqrt(dx*dx + dy * dy) * 0.015;
	}

	// 3. 特征保护：锁定锐利边和边界边
	// 注意：这里使用 bool 类型的 Property Map
	auto eif = mesh.add_property_map<TriangleMesh::Edge_index, bool>("e:is_feature", false).first;

	// 标记夹角大于 60 度的边
	PMP::detect_sharp_edges(mesh, 60.0, eif);

	// 手动标记所有边界边 (如果是开孔模型)
	for (auto e : mesh.edges()) {
		if (CGAL::is_border(e, mesh)) {
			put(eif, e, true);
		}
	}

	// 4. 执行重构
	// 如果 edge_is_feature_map 报错，请使用下方的 edge_is_constrained_map
	PMP::isotropic_remeshing(
		faces(mesh),
		target_edge_length,
		mesh,
		CGAL::parameters::number_of_iterations(5)
		.edge_is_constrained_map(eif)  // 改为 constrained_map 以兼容旧版本
		.protect_constraints(true)     // 开启约束保护，防止锯齿
	);

	// 5. 提取数据回传给外部 std::vector
	out_vertices.clear();
	out_faces.clear();
	std::map<vertex_descriptor, unsigned int> vid_map;
	unsigned int v_idx = 0;

	for (auto v : mesh.vertices()) {
		vid_map[v] = v_idx++;
		const Point_3& p = mesh.point(v);
		out_vertices.push_back(CGAL::to_double(p.x()));
		out_vertices.push_back(CGAL::to_double(p.y()));
		out_vertices.push_back(CGAL::to_double(p.z()));
	}

	for (auto f : mesh.faces()) {
		for (auto h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
			out_faces.push_back(vid_map[target(h, mesh)]);
		}
	}

	return true;
}

bool SkeletonExtractor::extractSkeleton(const std::string& input_path,
	std::vector<double>& mesh_v,
	std::vector<unsigned int>& mesh_f,
	SkeletonData& skel_out)
{
	// 1. 直接读取原始 STL 数据，完全不经过 Remesh 处理
	TriangleMesh mesh;
	if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(input_path, mesh) || mesh.is_empty()) {
		std::cout << "读取 STL 失败！" << std::endl;
		return false;
	}

	// 导出原始网格数据，用于主程序显示透明背景
	mesh_v.clear();
	mesh_f.clear();
	std::map<TriangleMesh::Vertex_index, unsigned int> vid_map;
	unsigned int v_idx = 0;
	for (auto v : mesh.vertices()) {
		vid_map[v] = v_idx++;
		const auto& p = mesh.point(v);
		mesh_v.push_back(CGAL::to_double(p.x()));
		mesh_v.push_back(CGAL::to_double(p.y()));
		mesh_v.push_back(CGAL::to_double(p.z()));
	}
	for (auto f : mesh.faces()) {
		for (auto h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
			mesh_f.push_back(vid_map[target(h, mesh)]);
	}

	// 2. 骨架提取逻辑
	typedef CGAL::Mean_curvature_flow_skeletonization<TriangleMesh> Skeletonization;
	typedef Skeletonization::Skeleton Skeleton;

	Skeletonization mcs(mesh);

	// 执行网格收缩
	mcs.contract_until_convergence();

	// --- 核心修复：CGAL 官方 API 名是 convert_to_skeleton ---
	Skeleton skeleton;
	mcs.convert_to_skeleton(skeleton);

	// 3. 将骨架数据转存到输出结构中
	skel_out.nodes.clear();
	skel_out.edges.clear();
	skel_out.radii.clear();
	skel_out.supportOffsets.clear();
	skel_out.supportPoints.clear();
	skel_out.supportOffsets.push_back(0);

	// 提取节点
	BOOST_FOREACH(Skeleton::vertex_descriptor v, vertices(skeleton)) {
		const auto& p = skeleton[v].point;
		skel_out.nodes.push_back(CGAL::to_double(p.x()));
		skel_out.nodes.push_back(CGAL::to_double(p.y()));
		skel_out.nodes.push_back(CGAL::to_double(p.z()));

		std::vector<double> localDistances;
		for (auto meshVertex : skeleton[v].vertices) {
			const auto& q = mesh.point(meshVertex);

			const double dx = CGAL::to_double(q.x()) - CGAL::to_double(p.x());
			const double dy = CGAL::to_double(q.y()) - CGAL::to_double(p.y());
			const double dz = CGAL::to_double(q.z()) - CGAL::to_double(p.z());

			localDistances.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
			skel_out.supportPoints.push_back(CGAL::to_double(q.x()));
			skel_out.supportPoints.push_back(CGAL::to_double(q.y()));
			skel_out.supportPoints.push_back(CGAL::to_double(q.z()));
		}

		skel_out.radii.push_back(percentileValue(localDistances, 0.75));
		skel_out.supportOffsets.push_back(
			static_cast<unsigned int>(skel_out.supportPoints.size() / 3));
	}

	// 提取边连接关系
	BOOST_FOREACH(Skeleton::edge_descriptor e, edges(skeleton)) {
		skel_out.edges.push_back((unsigned int)source(e, skeleton));
		skel_out.edges.push_back((unsigned int)target(e, skeleton));
	}

	qDebug() << "[Skeleton] 提取完成，节点数: " << skel_out.nodes.size() / 3 ;

	return true;
}


SkeletonExtractor::SkeletonData SkeletonExtractor::filterSkeletonInsideClosedSurface(
	const SkeletonData& input,
	vtkPolyData* closedSurface,
	double tolerance)
{
	if (!closedSurface ||
		closedSurface->GetNumberOfPoints() == 0 ||
		closedSurface->GetNumberOfCells() == 0 ||
		input.nodes.empty())
	{
		return input;
	}

	const size_t oldNodeCount = input.nodes.size() / 3;
	if (oldNodeCount == 0)
	{
		return input;
	}

	vtkSmartPointer<vtkPoints> queryPoints = vtkSmartPointer<vtkPoints>::New();
	queryPoints->SetNumberOfPoints(static_cast<vtkIdType>(oldNodeCount));

	for (size_t i = 0; i < oldNodeCount; ++i)
	{
		queryPoints->SetPoint(
			static_cast<vtkIdType>(i),
			input.nodes[3 * i],
			input.nodes[3 * i + 1],
			input.nodes[3 * i + 2]);
	}

	vtkSmartPointer<vtkPolyData> queryPolyData = vtkSmartPointer<vtkPolyData>::New();
	queryPolyData->SetPoints(queryPoints);

	vtkSmartPointer<vtkSelectEnclosedPoints> selector =
		vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selector->SetInputData(queryPolyData);
	selector->SetSurfaceData(closedSurface);
	selector->SetTolerance(tolerance);
	selector->Update();

	std::vector<char> keepNode(oldNodeCount, 1);
	std::vector<int> oldToNew(oldNodeCount, -1);

	for (size_t i = 0; i < oldNodeCount; ++i)
	{
		if (selector->IsInside(static_cast<vtkIdType>(i)))
		{
			keepNode[i] = 0;
		}
	}

	SkeletonData output;
	output.supportOffsets.push_back(0);

	size_t newIndex = 0;
	for (size_t oldIndex = 0; oldIndex < oldNodeCount; ++oldIndex)
	{
		if (!keepNode[oldIndex])
		{
			continue;
		}

		oldToNew[oldIndex] = static_cast<int>(newIndex++);

		output.nodes.push_back(input.nodes[3 * oldIndex]);
		output.nodes.push_back(input.nodes[3 * oldIndex + 1]);
		output.nodes.push_back(input.nodes[3 * oldIndex + 2]);

		output.radii.push_back(oldIndex < input.radii.size() ? input.radii[oldIndex] : 0.0);

		if (input.supportOffsets.size() == oldNodeCount + 1)
		{
			const unsigned int begin = input.supportOffsets[oldIndex];
			const unsigned int end = input.supportOffsets[oldIndex + 1];

			for (unsigned int k = begin; k < end; ++k)
			{
				const size_t p = static_cast<size_t>(3) * static_cast<size_t>(k);
				if (p + 2 < input.supportPoints.size())
				{
					output.supportPoints.push_back(input.supportPoints[p]);
					output.supportPoints.push_back(input.supportPoints[p + 1]);
					output.supportPoints.push_back(input.supportPoints[p + 2]);
				}
			}
		}

		output.supportOffsets.push_back(static_cast<unsigned int>(output.supportPoints.size() / 3));
	}

	std::set<std::pair<unsigned int, unsigned int>> uniqueEdges;
	for (size_t i = 0; i + 1 < input.edges.size(); i += 2)
	{
		const unsigned int oldA = input.edges[i];
		const unsigned int oldB = input.edges[i + 1];

		if (oldA >= oldNodeCount || oldB >= oldNodeCount)
		{
			continue;
		}

		const int newA = oldToNew[oldA];
		const int newB = oldToNew[oldB];

		if (newA < 0 || newB < 0 || newA == newB)
		{
			continue;
		}

		const unsigned int a = static_cast<unsigned int>(std::min(newA, newB));
		const unsigned int b = static_cast<unsigned int>(std::max(newA, newB));
		if (uniqueEdges.insert({ a, b }).second)
		{
			output.edges.push_back(static_cast<unsigned int>(newA));
			output.edges.push_back(static_cast<unsigned int>(newB));
		}
	}

	qDebug() << "[SkeletonFilter] non-design 内部骨架剔除完成："
		<< "原节点 =" << static_cast<int>(oldNodeCount)
		<< "剩余节点 =" << static_cast<int>(output.nodes.size() / 3)
		<< "原边 =" << static_cast<int>(input.edges.size() / 2)
		<< "剩余边 =" << static_cast<int>(output.edges.size() / 2);

	return output;
}
