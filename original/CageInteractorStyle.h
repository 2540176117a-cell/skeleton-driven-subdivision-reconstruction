class vtkPolyData;
#ifndef CAGE_INTERACTOR_STYLE_H
#define CAGE_INTERACTOR_STYLE_H

#include "StraightSkeletonSegmenter.h"
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkImplicitBoolean.h>
#include <vtkCylinder.h>
#include <vtkSphereSource.h>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <array>
#include <functional>
#include <QString>

// 引入 OpenSubdiv 核心引用
#include <opensubdiv/far/topologyRefiner.h>

struct CageSegment {
	double pStart[4][3];
	double pEnd[4][3];
};

// 专门为 OpenSubdiv 准备的轻量化计算点结构体
struct OsdVertex {
	float x, y, z;
	void Clear() { x = y = z = 0.0f; }
	void AddWithWeight(OsdVertex const & src, float weight) {
		x += weight * src.x;
		y += weight * src.y;
		z += weight * src.z;
	}
};

class CageInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
	enum class StraightCageTemplateMode {
		ThinRod,
		LargePlate
	};

	static CageInteractorStyle* New();
	vtkTypeMacro(CageInteractorStyle, vtkInteractorStyleTrackballCamera);

	void SetReferenceData(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkCellLocator> locator);
	void StandardizeCagePoints(double points[4][3], double axis[3], double up[3], double side[3]);

	virtual void OnMouseMove() override;
	virtual void OnLeftButtonDown() override;
	virtual void OnRightButtonDown() override;
	virtual void OnMiddleButtonDown() override;
	virtual void OnKeyPress() override;

	vtkSmartPointer<vtkActor> GetPreviewBoxActor() { return m_PreviewBoxActor; }
	vtkSmartPointer<vtkActor> GetFixedBoxActor() { return m_FixedBoxActor; }
	vtkSmartPointer<vtkActor> GetLinesActor() { return m_LinesActor; }

	void outfile();

	// 导出时用于融合原始 STP 中的 non-design 保留实体。
	// 这里只保存匹配所需的源数据，outfile() 内部会在用户确认后尝试 BooleanFuse。
	void SetNonDesignExportSource(
		const QString& stepPath,
		const std::vector<std::array<double, 3>>& elementCenters,
		const std::vector<char>& elementIsNonDesign,
		bool allStepSolidsAreNonDesign = false);
	void SetNonDesignReferenceSurface(vtkPolyData* nonDesignSurface);

	void SetSkeletonData(
		const std::vector<double>& nodes,
		const std::vector<unsigned int>& edges,
		const std::vector<double>& radii,
		const std::vector<unsigned int>& supportOffsets,
		const std::vector<double>& supportPoints,
		vtkSmartPointer<vtkPolyData> skelPolyData);
	const std::vector<double>& GetSkeletonNodes() const { return m_skelNodes; }
	const std::vector<unsigned int>& GetSkeletonEdges() const { return m_skelEdges; }
	const std::vector<double>& GetSkeletonRadii() const { return m_skelRadii; }

	// 骨架编辑模式：节点合并 / 删除骨架边 / 延伸或从已有节点插入骨架分支。
	void EnableSkeletonNodeEditMode(bool enable);
	void EnableSkeletonDeleteMode(bool enable);
	void EnableSkeletonExtendMode(bool enable);
	bool IsSkeletonNodeEditMode() const { return m_isSkeletonEditMode; }
	bool IsSkeletonDeleteMode() const { return m_isSkeletonDeleteMode; }
	bool IsSkeletonExtendMode() const { return m_isSkeletonExtendMode; }

	int BuildStraightCagesFromSegments(
		const std::vector<StraightSkeletonSegmenter::Segment>& segments,
		bool replaceExisting = true,
		StraightCageTemplateMode templateMode = StraightCageTemplateMode::ThinRod);

	void SetStraightSegmentsForEditing(
		std::vector<StraightSkeletonSegmenter::Segment>* segments,
		vtkSmartPointer<vtkActor> segmentActor,
		const std::vector<double>& skeletonNodes);

	void EnableStraightSegmentEditMode(bool enable);
	bool IsStraightSegmentEditMode() const { return m_isStraightSegmentEditMode; }
	void EnableStraightSegmentSelectMode(bool enable);
	bool IsStraightSegmentSelectMode() const { return m_isStraightSegmentSelectMode; }
	void EnableStraightCageEditMode(bool enable);
	bool IsStraightCageEditMode() const { return m_isStraightCageEditMode; }
	void EnableStraightCageRotateMode(bool enable);
	bool IsStraightCageRotateMode() const { return m_isStraightCageRotateMode; }
	void SetEditModeChangedCallback(std::function<void()> callback) { m_editModeChangedCallback = callback; }

	// 🚨 终极多通缝合接口，取代原有的 GenerateBranch 和 GenerateHBranch
	void GenerateNBranch();
	void GenerateNBranch3DSkeleton();
	// 自动拓扑按钮接口
	void RecognizeAutoTopology(vtkRenderer* ren);
	void BuildAutoTopology(vtkRenderer* ren);

protected:
	CageInteractorStyle();
	~CageInteractorStyle();

	void UpdatePreviewBox(double center[3], double axis[3], double up[3], double side[3], double width, double depth);
	void GetTransformCorners(vtkTransform* t, double corners[4][3]);

	// ==== 🚨 RMF 最小旋转标架核心算法 ====
	void ComputeRMF(const double x0[3], const double t0[3], const double u0[3],
		const double x1[3], const double t1[3], double u1[3], double v1[3]);

private:
	vtkSmartPointer<vtkPolyData> m_RefData;
	vtkSmartPointer<vtkCellLocator> m_RefLocator;
	vtkSmartPointer<vtkActor> m_PreviewBoxActor;

	int m_ClickState;
	vtkSmartPointer<vtkActor> m_FixedBoxActor;
	vtkSmartPointer<vtkTransform> m_StartTransform;

	vtkSmartPointer<vtkActor> m_LinesActor;
	vtkSmartPointer<vtkPolyData> m_LinesData;

	vtkSmartPointer<vtkPolyData> m_subD_Mesh;
	vtkSmartPointer<vtkActor> m_GlobalSmoothActor;

	// ==== 全局控制笼数据容器 ====
	std::vector<CageSegment> m_allCages;                                 // 直管身控制笼
	std::vector<vtkSmartPointer<vtkPolyData>> m_hubCages;                // 🚨 泛型多通接头 (Quad Hub)

	struct SuppressedCageSideFace {
		int cageIndex;
		int sideFaceIndex;
	};
	std::vector<SuppressedCageSideFace> m_suppressedCageSideFaces;       // 侧向接管时需要打开的 cage 侧面


	// ==== NBranch 二维核心复用接口 ====
	// GenerateNBranch() 原有外部入口不变；三维 1+2+2 只准备 ring，再调用这个核心。
	struct NBranchInputPort {
		int sourceCageIndex = -1;
		bool sourceIsStart = false;
		bool isVirtual = false;
		double ring[4][3];
		double preferredDir[3] = { 0.0, 0.0, 0.0 };
		bool hasPreferredDir = false;
	};

	struct NBranchBuildOptions {
		bool addSdf = true;
		bool rebuildNow = false;
		bool debugLog = true;

		bool forceHubCenter = false;
		double hubCenter[3] = { 0.0, 0.0, 0.0 };

		bool forceNormalUp = false;
		double normalUp[3] = { 0.0, 0.0, 1.0 };

		double radiusScale = 0.8;
		double minRadiusScale = 0.045;
	};

	bool GenerateNBranchFromRings(
		const std::vector<NBranchInputPort>& inputPorts,
		const NBranchBuildOptions& options);

	// N=3 二维分支专用：局部三角 hub。
	// A 端只作为底边方向参考，不再复制 A 输入端面；
	// B0/B1 决定尖端方向，骨架交点仅作为局部锚点。
	bool GenerateNBranchTriHub2D();

	void RebuildGlobalMesh();
	vtkSmartPointer<vtkActor> m_HistoryCagesActor;
	vtkSmartPointer<vtkActor> m_CageSolidActor;
	int m_ManualRotationOffset = 0;

	// 导出 / 拓扑问题线调试显示
	vtkSmartPointer<vtkActor> m_DebugExportWireActor;

	// 原始 STP non-design 保留区导出融合数据。
	// 注意：这只用于 outfile() 的最终导出阶段，不影响原有重建/细分/显示流程。
	QString m_nonDesignExportStepPath;
	std::vector<std::array<double, 3>> m_nonDesignExportElementCenters;
	std::vector<char> m_nonDesignExportIsNonDesign;
	bool m_nonDesignExportAllStepSolids = false;
	vtkSmartPointer<vtkPolyData> m_nonDesignReferenceSurface;


	// ==== OpenSubdiv 数据 ====
	std::vector<OsdVertex> m_controlVerts;
	std::vector<int> m_faceVerts;
	std::vector<int> m_vertsPerFace;
	std::vector<OsdVertex> m_refinedVerts;
	OpenSubdiv::v3_7_0::Far::TopologyRefiner* m_refiner;

	bool m_IsSnapped = false;
	bool m_SnappedIsStart = false;
	int m_SnappedCageIndex = -1;
	int m_SnappedFaceOffset = -1; // 0/1 are caps, 2~5 are side faces.
	double m_SnappedPoints[4][3];

	struct SelectedFace {
		int cageIndex;
		bool isStart;
	};

	std::vector<SelectedFace> m_selectedFacesForBranch;
	bool m_isBranchMode = false;
	vtkSmartPointer<vtkActor> m_BranchHighlightActor;
	vtkSmartPointer<vtkPolyData> m_BranchHighlightData;
	void UpdateBranchHighlight(vtkRenderer* ren);

	// ==== 骨架寻路与编辑数据 ====
	std::vector<double> m_skelNodes;
	std::vector<unsigned int> m_skelEdges;
	std::vector<double> m_skelRadii;
	std::vector<unsigned int> m_skelSupportOffsets;
	std::vector<double> m_skelSupportPoints;
	vtkSmartPointer<vtkActor> m_PathHighlightActor;
	vtkSmartPointer<vtkActor> m_AutoComplexTraceActor;
	vtkSmartPointer<vtkActor> m_AutoComplex3DTraceActor;
	vtkSmartPointer<vtkActor> m_AutoSingleLoopTraceActor;
	vtkSmartPointer<vtkActor> m_AutoMultiLoopTraceActor;
	bool m_useForcedLoopTopology = false;
	std::vector<int> m_forcedLoopCycleNodes;
	std::vector<std::vector<int>> m_forcedLoopConnectorPaths;
	struct LoopAnchorInfo {
		int node;
		double center[3];
		double normal[3];
	};
	std::vector<std::vector<int>> BuildSkeletonAdjacency() const;
	std::vector<int> FindPathBFS(int startIdx, int endIdx, const std::vector<std::vector<int>>& adjList);
	std::vector<int> FindPathBFSIgnoringEdge(int startIdx, int endIdx,
		const std::vector<std::vector<int>>& adjList,
		int blockedU,
		int blockedV);
	void ShowSkeletonPath(vtkRenderer* ren);
	int FindNearestSkeletonNodeToSelectedFace(const SelectedFace& sel,
		double outCenter[3],
		double outNormal[3]) const;


	std::vector<int> FindBestSkeletonCycle(const std::vector<LoopAnchorInfo>& anchors,
		const std::vector<std::vector<int>>& adjList,
		std::vector<std::vector<int>>& connectorPaths,
		bool& isMultiBranchLoop);
	bool DetectAndShowLoopTopology(vtkRenderer* ren);
	void ShowSkeletonCycle(vtkRenderer* ren,
		const std::vector<int>& cycleNodes,
		const std::vector<std::vector<int>>& connectorPaths,
		bool isMultiBranchLoop);
	bool GenerateSingleBranchLoopTopology(vtkRenderer* ren);
	void AutoDetectAndShowTopology(vtkRenderer* ren, bool generateTopology = false);
	void ShowAutoTopologyTrace(vtkRenderer* ren,
		vtkSmartPointer<vtkActor>& actor,
		const std::vector<std::vector<int>>& paths,
		double r,
		double g,
		double b,
		double lineWidth);

	bool m_isSkeletonEditMode = false;
	bool m_isSkeletonDeleteMode = false;
	bool m_isSkeletonExtendMode = false;
	int m_selectedSkelNode1 = -1;
	int m_selectedExtendSkelNode = -1;
	vtkSmartPointer<vtkPolyData> m_skelPolyData;
	vtkSmartPointer<vtkActor> m_skelNodeHighlightActor;
	vtkSmartPointer<vtkActor> m_skeletonExtendPreviewActor;
	int PickSkeletonNode(int x, int y, vtkRenderer* ren);
	int PickSkeletonEdge(int x, int y, vtkRenderer* ren) const;
	void UpdateSkeletonHighlight(vtkRenderer* ren);
	void RefreshSkeletonPolyData();
	void MergeSkeletonNodes(int n1, int n2);
	bool DeleteSkeletonEdge(int edgeIndex);
	bool ExtendSkeletonFromSelectedTail(int x, int y, vtkRenderer* ren);
	bool UpdateSkeletonExtendPreview(int x, int y, vtkRenderer* ren);
	bool GetSkeletonExtensionTargetPoint(int x, int y, vtkRenderer* ren, int tail, double outPoint[3]) const;
	void SmoothSkeletonAroundNode(int centerNode, int iterations = 3);
	void ClearSkeletonExtendSelection(vtkRenderer* ren = nullptr);
	void HideSkeletonExtendPreview();
	int GetSkeletonNodeDegree(int nodeIndex) const;
	int GetSingleSkeletonNeighbor(int nodeIndex) const;
	double EstimateSkeletonMeanEdgeLength() const;
	void SmoothSkeleton(int iterations = 10);
	void NotifyEditModeChanged();
	void HideManualPreviewActors();
	void ExitAllEditModes();


	// ==== 直线段识别结果编辑 ====
	bool m_isStraightSegmentEditMode = false;
	bool m_isStraightSegmentSelectMode = false;
	std::vector<StraightSkeletonSegmenter::Segment>* m_editStraightSegments = nullptr;
	std::vector<double> m_editStraightSkeletonNodes;
	vtkSmartPointer<vtkActor> m_editStraightSegmentActor;

	bool HandleStraightSegmentTrimClick(int x, int y, vtkRenderer* ren);
	bool HandleStraightSegmentSelectClick(int x, int y, vtkRenderer* ren);
	bool BuildStraightSegmentFromSkeletonEdge(int edgeIndex, StraightSkeletonSegmenter::Segment& segment) const;
	void RefreshStraightSegmentEditActor();

	enum class StraightCageEditTargetType {
		None,
		Vertex,
		Face,
		Edge,
		HubVertex,
		HubFace,
		HubEdge
	};

	struct StraightCageEditTarget {
		StraightCageEditTargetType type = StraightCageEditTargetType::None;
		int cageIndex = -1;
		int faceIndex = -1;
		int vertexIndex = -1;
		int edgeVertexA = -1;
		int edgeVertexB = -1;
		int hubIndex = -1;
		int hubFaceIndex = -1;
		vtkIdType hubPointId = -1;
		vtkIdType hubPointA = -1;
		vtkIdType hubPointB = -1;
		std::vector<vtkIdType> hubFacePointIds;
		double center[3] = { 0.0, 0.0, 0.0 };
		double normal[3] = { 0.0, 0.0, 1.0 };
		double alternateNormal[3] = { 0.0, 0.0, 1.0 };
		double edgeDirection[3] = { 0.0, 0.0, 0.0 };
		bool hasAlternateNormal = false;
		double originalA[3] = { 0.0, 0.0, 0.0 };
		double originalB[3] = { 0.0, 0.0, 0.0 };
		double minOffset = -1e100;
		double maxOffset = 1e100;
	};

	bool m_isStraightCageEditMode = false;
	bool m_straightCageEditActive = false;
	StraightCageEditTarget m_straightCageEditTarget;
	int m_straightCageEditStartMouse[2] = { 0, 0 };
	double m_straightCageEditCurrentOffset = 0.0;
	double m_straightCageEditCurrentDelta[3] = { 0.0, 0.0, 0.0 };
	double m_straightCageEditReferenceStep = 1.0;
	vtkSmartPointer<vtkActor> m_straightCageEditHighlightActor;
	vtkSmartPointer<vtkActor> m_straightCageEditPreviewActor;
	bool m_isStraightCageRotateMode = false;
	bool m_straightCageRotateActive = false;
	StraightCageEditTarget m_straightCageRotateTarget;
	int m_straightCageRotateStartMouse[2] = { 0, 0 };
	double m_straightCageRotateStartDisplayAngle = 0.0;
	double m_straightCageRotateCurrentAngle = 0.0;
	std::function<void()> m_editModeChangedCallback;

	bool PickStraightCageEditTarget(int x, int y, vtkRenderer* ren, StraightCageEditTarget& target);
	bool HandleStraightCageEditLeftClick(int x, int y, vtkRenderer* ren);
	bool UpdateStraightCageEditPreview(int x, int y, vtkRenderer* ren);
	void ApplyStraightCageEdit();
	bool RotateStraightCageEditTarget(double angleRadians);
	bool HandleStraightCageRotateLeftClick(int x, int y, vtkRenderer* ren);
	bool UpdateStraightCageRotatePreview(int x, int y, vtkRenderer* ren);
	void ApplyStraightCageRotation();
	void CancelStraightCageRotation(vtkRenderer* ren = nullptr);
	void UpdateStraightCageRotateActors(vtkRenderer* ren, bool showPreview);
	void CancelStraightCageEdit(vtkRenderer* ren = nullptr);
	void UpdateStraightCageEditActors(vtkRenderer* ren, bool showPreview);

	// ==== 🚨 论文创新点：混合隐函数投影数据结构 ====
	struct SDFTube
	{
		double center[3];
		double axis[3];   // 必须是归一化的方向向量
		double radius;
	};

	struct HubSDFInfo
	{
		double hubCenter[3];
		double maxRadius;
		double coreRadius = 0.0;
		std::vector<SDFTube> tubes;
	};

	std::vector<HubSDFInfo> m_sdfHubs; // 存储场景中所有多通接头的数学隐函数场参数
	bool m_use3DHubSdfProjection = false;

	// 隐函数牛顿投影算子
	void ProjectToSmoothSDF(vtkSmartPointer<vtkPolyData> mesh);

};

#endif
