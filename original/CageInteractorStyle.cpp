#include "CageInteractorStyle.h"
#include<QDebug>
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkRenderWindow.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkStripper.h>
#include <vtkAppendPolyData.h>
#include <vtkIdList.h>
#include <QFileDialog>
#include <QMessageBox>
#include <vtkModifiedBSPTree.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <array>
#include <algorithm>
#include <limits>
#include <numeric>
#include <tuple>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkTubeFilter.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <windows.h>
#include <vtkOBBTree.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkGenericCell.h>


// ============================================================
// NBranch 局部 hubRing 生成期让步控制
// 只在 GenerateNBranch3DSkeleton() 的 1+2+2 main-only 模式中短暂启用。
// 目的：不要先生成 hub/bridge 再后处理点位，避免 hubPoly 和 bridge.pStart 分裂成两根线。
// ============================================================
namespace
{
	struct NBranchHubRingRelaxContext
	{
		bool enabled = false;
		int targetInputIndex = -1;
		double axisAtoB[3] = { 1.0, 0.0, 0.0 };
		double aCenter[3] = { 0.0, 0.0, 0.0 };
		double requestMove = 0.0;
		double safeGap = 0.0;
		double maxPortDim = 1.0;
		std::vector<std::array<double, 3>> boundaryPoints;
		std::array<std::vector<std::array<double, 3>>, 4> boundaryPointGroups;
	};

	static NBranchHubRingRelaxContext g_NBranchHubRingRelax;

	static double Dot3Local(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	static double Distance2PointArrayLocal(const double p[3], const std::array<double, 3>& q)
	{
		double dx = p[0] - q[0];
		double dy = p[1] - q[1];
		double dz = p[2] - q[2];
		return dx * dx + dy * dy + dz * dz;
	}

	static double AxisTFromALocal(const double p[3], const NBranchHubRingRelaxContext& ctx)
	{
		double v[3] = {
			p[0] - ctx.aCenter[0],
			p[1] - ctx.aCenter[1],
			p[2] - ctx.aCenter[2]
		};
		return Dot3Local(v, ctx.axisAtoB);
	}

	static void PerpFromAxisLocal(
		const double p[3],
		const NBranchHubRingRelaxContext& ctx,
		double out[3])
	{
		double t = AxisTFromALocal(p, ctx);
		out[0] = p[0] - (ctx.aCenter[0] + ctx.axisAtoB[0] * t);
		out[1] = p[1] - (ctx.aCenter[1] + ctx.axisAtoB[1] * t);
		out[2] = p[2] - (ctx.aCenter[2] + ctx.axisAtoB[2] * t);
	}

	static double BoundaryPointTLocal(
		const std::array<double, 3>& p,
		const NBranchHubRingRelaxContext& ctx)
	{
		double v[3] = {
			p[0] - ctx.aCenter[0],
			p[1] - ctx.aCenter[1],
			p[2] - ctx.aCenter[2]
		};
		return Dot3Local(v, ctx.axisAtoB);
	}


	static void ClampNBranchHubPolyAgainstBEndNodesInPlace(vtkPoints* pts)
	{
		NBranchHubRingRelaxContext& ctx = g_NBranchHubRingRelax;
		if (!pts || !ctx.enabled || ctx.boundaryPoints.empty()) {
			return;
		}

		// v8：真正按“靠近 2 分支的 hub 节点”和“2 分支端面的 8 个节点”做硬比较。
		// 不再主动把 hub 往 B 方向推，也不再用插值边界/统一截面。
		// 对每一个 hubPoly 顶点：
		//   1) 找到 B0/B1 两个输入端面 8 个角点中，在非主方向坐标上最近的那个角点；
		//   2) 只比较主方向坐标 S；
		//   3) 如果 hub 顶点越过这个 B 端角点，就沿真实 A->B 轴反向拉回。
		double axis[3] = { ctx.axisAtoB[0], ctx.axisAtoB[1], ctx.axisAtoB[2] };
		double axisNorm = std::sqrt(Dot3Local(axis, axis));
		if (axisNorm < 1e-8) {
			axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
			axisNorm = 1.0;
		}
		axis[0] /= axisNorm;
		axis[1] /= axisNorm;
		axis[2] /= axisNorm;

		int majorDim = 0;
		if (std::abs(axis[1]) > std::abs(axis[majorDim])) majorDim = 1;
		if (std::abs(axis[2]) > std::abs(axis[majorDim])) majorDim = 2;
		const double sign = (axis[majorDim] >= 0.0) ? 1.0 : -1.0;
		const double axisMajorAbs = std::max(1e-8, std::abs(axis[majorDim]));

		auto SCoord = [&](const double p[3]) -> double {
			return sign * p[majorDim];
		};
		auto SCoordArray = [&](const std::array<double, 3>& p) -> double {
			return sign * p[majorDim];
		};
		auto PerpDist2ByMajorCoord = [&](const double p[3], const std::array<double, 3>& q) -> double {
			double d2 = 0.0;
			for (int d = 0; d < 3; ++d) {
				if (d == majorDim) continue;
				double dv = p[d] - q[d];
				d2 += dv * dv;
			}
			return d2;
		};

		const double safeGap = std::max(ctx.safeGap, ctx.maxPortDim * 0.060);
		int clampedCount = 0;
		double minGap = std::numeric_limits<double>::max();

		vtkIdType nPts = pts->GetNumberOfPoints();
		for (vtkIdType i = 0; i < nPts; ++i) {
			double p[3];
			pts->GetPoint(i, p);

			int bestNode = -1;
			double bestD2 = std::numeric_limits<double>::max();
			for (int j = 0; j < static_cast<int>(ctx.boundaryPoints.size()); ++j) {
				double d2 = PerpDist2ByMajorCoord(p, ctx.boundaryPoints[j]);
				if (d2 < bestD2) {
					bestD2 = d2;
					bestNode = j;
				}
			}
			if (bestNode < 0) continue;

			double pS = SCoord(p);
			double entryS = SCoordArray(ctx.boundaryPoints[bestNode]);
			double limitS = entryS - safeGap;

			if (pS > limitS) {
				double deltaS = pS - limitS;
				double pullBack = deltaS / axisMajorAbs;
				p[0] -= axis[0] * pullBack;
				p[1] -= axis[1] * pullBack;
				p[2] -= axis[2] * pullBack;
				pts->SetPoint(i, p);
				pS = limitS;
				++clampedCount;
			}

			minGap = std::min(minGap, entryS - pS);
			if (false) qDebug() << "[NBranchHubRingRelax][v8Exact8NodeClamp]"
				<< "hubPt =" << i
				<< "bestBNode =" << bestNode
				<< "majorDim =" << majorDim
				<< "sign =" << sign
				<< "hubS =" << pS
				<< "entryS =" << entryS
				<< "limitS =" << limitS
				<< "gap =" << (entryS - pS)
				<< "clamped =" << (pS >= limitS - 1e-8);
		}

		if (clampedCount > 0) {
			pts->Modified();
		}
		if (false) qDebug() << "[NBranchHubRingRelax] applied v8 exact 8-node clamp"
			<< "hubPointCount =" << nPts
			<< "bNodeCount =" << static_cast<int>(ctx.boundaryPoints.size())
			<< "majorDim =" << majorDim
			<< "axis =" << axis[0] << axis[1] << axis[2]
			<< "sign =" << sign
			<< "safeGap =" << safeGap
			<< "clampedPoints =" << clampedCount
			<< "minGapToB =" << minGap;
	}

	static void ApplyNBranchHubRingRelaxInPlace(
		double hubRing[4][3],
		const int sourceInputIndex)
	{
		NBranchHubRingRelaxContext& ctx = g_NBranchHubRingRelax;
		if (!ctx.enabled || sourceInputIndex != ctx.targetInputIndex) {
			return;
		}
		if (ctx.boundaryPoints.empty()) {
			return;
		}

		// v11：最直接的 B 端面硬限制。
		// 不再猜 hubRing[k] 和 B.ring[k] 的编号对应关系；不做采样；不做半空间；不横向推点。
		// 对每一个 hub 顶点：
		//   1) 先按原本的 A->B 主骨架方向做避让移动；
		//   2) 把该 hub 点和 B0/B1 两个输入端面的全部原始顶点，在主方向坐标 S 上逐一比较；
		//   3) B 端面所有顶点中 S 最小的那个，就是最靠近 hub 的入口边界；
		//   4) 如果当前 hub 点超过这个入口边界，就只把当前这个 hub 点沿 -axisAtoB 回缩。
		double axis[3] = { ctx.axisAtoB[0], ctx.axisAtoB[1], ctx.axisAtoB[2] };
		double axisNorm = std::sqrt(Dot3Local(axis, axis));
		if (axisNorm < 1e-8) {
			axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
			axisNorm = 1.0;
		}
		axis[0] /= axisNorm;
		axis[1] /= axisNorm;
		axis[2] /= axisNorm;

		int majorDim = 0;
		if (std::abs(axis[1]) > std::abs(axis[majorDim])) majorDim = 1;
		if (std::abs(axis[2]) > std::abs(axis[majorDim])) majorDim = 2;
		const double sign = (axis[majorDim] >= 0.0) ? 1.0 : -1.0;
		const double axisMajorAbs = std::max(1e-8, std::abs(axis[majorDim]));

		auto SCoord = [&](const double p[3]) -> double
		{
			return sign * p[majorDim];
		};
		auto SCoordArray = [&](const std::array<double, 3>& p) -> double
		{
			return sign * p[majorDim];
		};

		double bMinS = std::numeric_limits<double>::max();
		int bMinIndex = -1;
		for (int i = 0; i < static_cast<int>(ctx.boundaryPoints.size()); ++i) {
			const double s = SCoordArray(ctx.boundaryPoints[i]);
			if (s < bMinS) {
				bMinS = s;
				bMinIndex = i;
			}
		}
		if (bMinIndex < 0 || bMinS >= std::numeric_limits<double>::max() * 0.25) {
			return;
		}

		const double requestMove = std::max(ctx.requestMove, ctx.maxPortDim * 0.18);
		const double safeGap = std::max(ctx.safeGap, ctx.maxPortDim * 0.055);
		const double limitS = bMinS - safeGap;

		int clampedCount = 0;
		double minGapToB = std::numeric_limits<double>::max();
		double minMoveS = std::numeric_limits<double>::max();
		double maxMoveS = -std::numeric_limits<double>::max();

		if (false) qDebug() << "[NBranchHubRingRelax][v11DirectAllBVerticesClamp] use direct selected B face vertices only"
			<< "bVertexCount =" << static_cast<int>(ctx.boundaryPoints.size())
			<< "bMinIndex =" << bMinIndex
			<< "bMinS =" << bMinS
			<< "limitS =" << limitS
			<< "majorDim =" << majorDim
			<< "sign =" << sign;

		for (int k = 0; k < 4; ++k) {
			double original[3] = { hubRing[k][0], hubRing[k][1], hubRing[k][2] };
			double moved[3] = {
				original[0] + axis[0] * requestMove,
				original[1] + axis[1] * requestMove,
				original[2] + axis[2] * requestMove
			};

			const double originalS = SCoord(original);
			double movedS = SCoord(moved);
			bool clamped = false;

			// 每个 hub 点都直接和 B 端面全部顶点形成的最靠近入口限制比较。
			// 只要避让后超过 B 端面最近入口，就沿 -axisAtoB 回退当前点。
			if (movedS > limitS) {
				const double deltaS = movedS - limitS;
				const double pullBack = deltaS / axisMajorAbs;
				moved[0] -= axis[0] * pullBack;
				moved[1] -= axis[1] * pullBack;
				moved[2] -= axis[2] * pullBack;
				movedS = limitS;
				clamped = true;
				++clampedCount;
			}

			hubRing[k][0] = moved[0];
			hubRing[k][1] = moved[1];
			hubRing[k][2] = moved[2];

			minGapToB = std::min(minGapToB, bMinS - movedS);
			minMoveS = std::min(minMoveS, movedS - originalS);
			maxMoveS = std::max(maxMoveS, movedS - originalS);

			if (false) qDebug() << "[NBranchHubRingRelax][v11DirectAllBVerticesClamp]"
				<< "k =" << k
				<< "originalS =" << originalS
				<< "movedS =" << movedS
				<< "bMinS =" << bMinS
				<< "limitS =" << limitS
				<< "gap =" << (bMinS - movedS)
				<< "clamped =" << clamped;
		}

		if (false) qDebug() << "[NBranchHubRingRelax] applied v11 direct all-B-vertices clamp"
			<< "sourceInputIndex =" << sourceInputIndex
			<< "axis =" << axis[0] << axis[1] << axis[2]
			<< "majorDim =" << majorDim
			<< "sign =" << sign
			<< "requestMove =" << requestMove
			<< "safeGap =" << safeGap
			<< "bMinS =" << bMinS
			<< "limitS =" << limitS
			<< "clampedPoints =" << clampedCount
			<< "minMoveS =" << minMoveS
			<< "maxMoveS =" << maxMoveS
			<< "minGapToB =" << minGapToB;
	}

	static void ClampNBranchHubBsideVerticesAgainstBInPlace(
		vtkPoints* pts,
		const int targetSortedK)
	{
		NBranchHubRingRelaxContext& ctx = g_NBranchHubRingRelax;
		if (!pts || !ctx.enabled || ctx.boundaryPoints.empty()) {
			return;
		}
		const vtkIdType nPts = pts->GetNumberOfPoints();
		if (nPts < 6 || targetSortedK < 0) {
			return;
		}

		// N=3 主 fan 中，靠近 B0/B1 两个端面的 hub 节点不是 A 侧接口四边形的 4 个点，
		// 而是剩下的那个 sector 的上下 2 个点：sector = targetSortedK + 1。
		// 这里不再给这两个点施加“避让前推”，只做 B 端面硬限制：
		// 如果它们沿主骨架方向超过 B0/B1 输入端面 8 个原始顶点形成的最近入口线，
		// 就只沿 -axisAtoB 回退当前这个点。
		const int N = static_cast<int>(nPts / 2);
		if (N != 3) {
			if (false) qDebug() << "[NBranchHubRingRelax][v12BsideTwoNodeClamp] skipped because hub point count is not 6"
				<< "nPts =" << nPts
				<< "N =" << N;
			return;
		}

		const int bSideSector = (targetSortedK + 1) % N;
		const vtkIdType bSidePt0 = static_cast<vtkIdType>(2 * bSideSector);
		const vtkIdType bSidePt1 = static_cast<vtkIdType>(2 * bSideSector + 1);

		double axis[3] = { ctx.axisAtoB[0], ctx.axisAtoB[1], ctx.axisAtoB[2] };
		double axisNorm = std::sqrt(Dot3Local(axis, axis));
		if (axisNorm < 1e-8) {
			axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
			axisNorm = 1.0;
		}
		axis[0] /= axisNorm;
		axis[1] /= axisNorm;
		axis[2] /= axisNorm;

		int majorDim = 0;
		if (std::abs(axis[1]) > std::abs(axis[majorDim])) majorDim = 1;
		if (std::abs(axis[2]) > std::abs(axis[majorDim])) majorDim = 2;
		const double sign = (axis[majorDim] >= 0.0) ? 1.0 : -1.0;
		const double axisMajorAbs = std::max(1e-8, std::abs(axis[majorDim]));

		auto SCoord = [&](const double p[3]) -> double
		{
			return sign * p[majorDim];
		};
		auto SCoordArray = [&](const std::array<double, 3>& p) -> double
		{
			return sign * p[majorDim];
		};

		double bMinS = std::numeric_limits<double>::max();
		int bMinIndex = -1;
		for (int i = 0; i < static_cast<int>(ctx.boundaryPoints.size()); ++i) {
			const double s = SCoordArray(ctx.boundaryPoints[i]);
			if (s < bMinS) {
				bMinS = s;
				bMinIndex = i;
			}
		}
		if (bMinIndex < 0 || bMinS >= std::numeric_limits<double>::max() * 0.25) {
			return;
		}

		const double safeGap = std::max(ctx.safeGap, ctx.maxPortDim * 0.055);
		const double sectionGap = std::max(ctx.safeGap * 1.35, ctx.maxPortDim * 0.075);
		const double limitS = bMinS - safeGap;

		// 重要补丁：B 侧两个 hub 节点被 B 端硬限制回退时，不能压到 A 侧已经避让后的四点截面上。
		// 先统计除 B 侧两个点以外的四个 hub 点在 A->B 主方向上的最大 S。
		// 正常情况下应满足：B侧节点S >= A侧截面maxS + sectionGap，同时 B侧节点S <= B端入口limitS。
		// 如果空间不足，则优先保持 B 端安全距离，把 A 侧四点整体沿 -axisAtoB 再回退，让两截面重新拉开。
		double otherMaxS = -std::numeric_limits<double>::max();
		std::vector<vtkIdType> otherIds;
		otherIds.reserve(static_cast<size_t>(std::max<vtkIdType>(0, nPts - 2)));
		for (vtkIdType oid = 0; oid < nPts; ++oid) {
			if (oid == bSidePt0 || oid == bSidePt1) {
				continue;
			}
			double op[3];
			pts->GetPoint(oid, op);
			otherMaxS = std::max(otherMaxS, SCoord(op));
			otherIds.push_back(oid);
		}

		if (otherMaxS <= -std::numeric_limits<double>::max() * 0.25) {
			otherMaxS = limitS - sectionGap;
		}

		double minAllowedBsideS = otherMaxS + sectionGap;
		int sectionPullbackCount = 0;
		double sectionPullbackS = 0.0;
		if (minAllowedBsideS > limitS) {
			sectionPullbackS = minAllowedBsideS - limitS;
			const double pullBack = sectionPullbackS / axisMajorAbs;
			for (vtkIdType oid : otherIds) {
				double op[3];
				pts->GetPoint(oid, op);
				op[0] -= axis[0] * pullBack;
				op[1] -= axis[1] * pullBack;
				op[2] -= axis[2] * pullBack;
				pts->SetPoint(oid, op);
				++sectionPullbackCount;
			}
			otherMaxS -= sectionPullbackS;
			minAllowedBsideS = otherMaxS + sectionGap;
		}

		int clampedCount = 0;
		double minGapToB = std::numeric_limits<double>::max();
		double minGapToASection = std::numeric_limits<double>::max();

		if (false) qDebug() << "[NBranchHubRingRelax][v13BsideSectionGapClamp] check the remaining B-side two hub vertices"
			<< "targetSortedK =" << targetSortedK
			<< "bSideSector =" << bSideSector
			<< "bSidePt0 =" << bSidePt0
			<< "bSidePt1 =" << bSidePt1
			<< "bVertexCount =" << static_cast<int>(ctx.boundaryPoints.size())
			<< "bMinIndex =" << bMinIndex
			<< "bMinS =" << bMinS
			<< "limitS =" << limitS
			<< "majorDim =" << majorDim
			<< "sign =" << sign;

		const vtkIdType ids[2] = { bSidePt0, bSidePt1 };
		for (int k = 0; k < 2; ++k) {
			const vtkIdType pid = ids[k];
			if (pid < 0 || pid >= nPts) {
				continue;
			}

			double p[3];
			pts->GetPoint(pid, p);
			double pS = SCoord(p);
			const double beforeS = pS;
			bool clamped = false;

			if (pS > limitS) {
				const double deltaS = pS - limitS;
				const double pullBack = deltaS / axisMajorAbs;
				p[0] -= axis[0] * pullBack;
				p[1] -= axis[1] * pullBack;
				p[2] -= axis[2] * pullBack;
				pS = limitS;
				clamped = true;
			}

			if (pS < minAllowedBsideS) {
				const double deltaS = minAllowedBsideS - pS;
				const double pushForward = deltaS / axisMajorAbs;
				p[0] += axis[0] * pushForward;
				p[1] += axis[1] * pushForward;
				p[2] += axis[2] * pushForward;
				pS = minAllowedBsideS;
				clamped = true;
			}

			pts->SetPoint(pid, p);
			if (clamped) {
				++clampedCount;
			}

			minGapToB = std::min(minGapToB, bMinS - pS);
			minGapToASection = std::min(minGapToASection, pS - otherMaxS);
			if (false) qDebug() << "[NBranchHubRingRelax][v13BsideSectionGapClamp]"
				<< "pid =" << pid
				<< "beforeS =" << beforeS
				<< "afterS =" << pS
				<< "bMinS =" << bMinS
				<< "limitS =" << limitS
				<< "otherMaxS =" << otherMaxS
				<< "minAllowedBsideS =" << minAllowedBsideS
				<< "gapToB =" << (bMinS - pS)
				<< "gapToASection =" << (pS - otherMaxS)
				<< "clamped =" << clamped;
		}

		if (clampedCount > 0) {
			pts->Modified();
		}

		if (false) qDebug() << "[NBranchHubRingRelax] applied v13 B-side two-node clamp with A-section safety gap"
			<< "targetSortedK =" << targetSortedK
			<< "bSideSector =" << bSideSector
			<< "axis =" << axis[0] << axis[1] << axis[2]
			<< "safeGap =" << safeGap
			<< "sectionGap =" << sectionGap
			<< "bMinS =" << bMinS
			<< "limitS =" << limitS
			<< "otherMaxS =" << otherMaxS
			<< "minAllowedBsideS =" << minAllowedBsideS
			<< "sectionPullbackCount =" << sectionPullbackCount
			<< "sectionPullbackS =" << sectionPullbackS
			<< "clampedPoints =" << clampedCount
			<< "minGapToB =" << minGapToB
			<< "minGapToASection =" << minGapToASection;
	}

}


// --- OCCT ---
#include <gp_Pnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRep_Tool.hxx>
#include <ShapeFix_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <TopExp.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <STEPControl_Writer.hxx>
#include <STEPControl_Reader.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <Geom_BezierSurface.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp_Explorer.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
// --- OpenSubdiv ---
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <set>

static void AppendTopoDSWireToDebugPolyData(
	const TopoDS_Wire& wire,
	vtkSmartPointer<vtkPoints> pts,
	vtkSmartPointer<vtkCellArray> lines,
	int samplesPerEdge = 12)
{
	TopExp_Explorer edgeExp(wire, TopAbs_EDGE);

	while (edgeExp.More()) {
		TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

		Standard_Real first = 0.0;
		Standard_Real last = 0.0;

		Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);

		if (!curve.IsNull()) {
			vtkIdType firstId = -1;
			vtkIdType prevId = -1;

			for (int i = 0; i <= samplesPerEdge; ++i) {
				double t = static_cast<double>(i) / static_cast<double>(samplesPerEdge);
				Standard_Real u = first + (last - first) * t;

				gp_Pnt p = curve->Value(u);

				vtkIdType id = pts->InsertNextPoint(p.X(), p.Y(), p.Z());

				if (i == 0) {
					firstId = id;
				}

				if (prevId >= 0) {
					vtkNew<vtkLine> line;
					line->GetPointIds()->SetId(0, prevId);
					line->GetPointIds()->SetId(1, id);
					lines->InsertNextCell(line);
				}

				prevId = id;
			}
		}

		edgeExp.Next();
	}
}
static vtkSmartPointer<vtkActor> MakeDebugWireActor(vtkPolyData* wirePoly)
{
	vtkNew<vtkTubeFilter> tube;
	tube->SetInputData(wirePoly);
	tube->SetRadius(0.35);
	tube->SetNumberOfSides(8);
	tube->CappingOn();
	tube->Update();

	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(tube->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	actor->GetProperty()->SetOpacity(1.0);
	actor->GetProperty()->SetLighting(false);
	actor->SetPickable(0);

	return actor;
}
static void DebugPolyDataBoundary(vtkPolyData* poly, const char* tag)
{
	if (!poly) {
		if (false) qDebug() << "[BoundaryDebug]" << tag << "poly is null";
		return;
	}

	vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(poly);
	cleaner->SetTolerance(1e-7);
	cleaner->Update();

	vtkPolyData* cleanPoly = cleaner->GetOutput();

	vtkSmartPointer<vtkFeatureEdges> boundary =
		vtkSmartPointer<vtkFeatureEdges>::New();
	boundary->SetInputData(cleanPoly);
	boundary->BoundaryEdgesOn();
	boundary->FeatureEdgesOff();
	boundary->ManifoldEdgesOff();
	boundary->NonManifoldEdgesOn();
	boundary->Update();

	vtkPolyData* edgePoly = boundary->GetOutput();

	if (false) qDebug() << "[BoundaryDebug]" << tag
		<< "points =" << cleanPoly->GetNumberOfPoints()
		<< "cells =" << cleanPoly->GetNumberOfCells()
		<< "boundary/nonmanifold edge cells =" << edgePoly->GetNumberOfCells()
		<< "edge points =" << edgePoly->GetNumberOfPoints();

	if (edgePoly->GetNumberOfCells() > 0) {
		double b[6];
		edgePoly->GetBounds(b);

		if (false) qDebug() << "[BoundaryDebug]" << tag
			<< "edge bounds x=" << b[0] << b[1]
			<< "y=" << b[2] << b[3]
			<< "z=" << b[4] << b[5];
	}
}
static void DebugPolyDataFaceTypes(vtkPolyData* poly, const char* tag)
{
	if (!poly) {
		if (false) qDebug() << "[FaceTypeDebug]" << tag << "poly is null";
		return;
	}

	vtkIdType triCount = 0;
	vtkIdType quadCount = 0;
	vtkIdType ngonCount = 0;
	vtkIdType lineOrVertexCount = 0;

	vtkCellArray* polys = poly->GetPolys();
	if (!polys) {
		if (false) qDebug() << "[FaceTypeDebug]" << tag << "no polys";
		return;
	}

	vtkSmartPointer<vtkIdList> idList =
		vtkSmartPointer<vtkIdList>::New();

	polys->InitTraversal();

	while (polys->GetNextCell(idList)) {
		vtkIdType npts = idList->GetNumberOfIds();

		if (npts == 3) {
			triCount++;
		}
		else if (npts == 4) {
			quadCount++;
		}
		else if (npts > 4) {
			ngonCount++;
		}
		else {
			lineOrVertexCount++;
		}
	}

	if (false) qDebug() << "[FaceTypeDebug]" << tag
		<< "tri =" << triCount
		<< "quad =" << quadCount
		<< "ngon =" << ngonCount
		<< "line/vertex =" << lineOrVertexCount;
}

using namespace OpenSubdiv;

vtkStandardNewMacro(CageInteractorStyle);
namespace
{
	double StraightEditDistance3(const double a[3], const double b[3])
	{
		const double dx = a[0] - b[0];
		const double dy = a[1] - b[1];
		const double dz = a[2] - b[2];
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	double StraightEditDot3(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	double StraightEditPolylineLength(const std::vector<double>& pts)
	{
		if (pts.size() < 6) {
			return 0.0;
		}

		double length = 0.0;

		for (size_t i = 3; i + 2 < pts.size(); i += 3) {
			double a[3] = { pts[i - 3], pts[i - 2], pts[i - 1] };
			double b[3] = { pts[i], pts[i + 1], pts[i + 2] };
			length += StraightEditDistance3(a, b);
		}

		return length;
	}

	double StraightEditSkeletonDiagonal(const std::vector<double>& nodes)
	{
		if (nodes.size() < 3) {
			return 1.0;
		}

		double minP[3] = { nodes[0], nodes[1], nodes[2] };
		double maxP[3] = { nodes[0], nodes[1], nodes[2] };

		for (size_t i = 0; i + 2 < nodes.size(); i += 3) {
			for (int k = 0; k < 3; ++k) {
				minP[k] = std::min(minP[k], nodes[i + k]);
				maxP[k] = std::max(maxP[k], nodes[i + k]);
			}
		}

		const double dx = maxP[0] - minP[0];
		const double dy = maxP[1] - minP[1];
		const double dz = maxP[2] - minP[2];

		const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
		return diag > 1e-8 ? diag : 1.0;
	}

	bool StraightEditTrimSegmentAtPoint(
		StraightSkeletonSegmenter::Segment& segment,
		const double pickedPoint[3],
		double modelDiagonal)
	{
		if (!segment.valid || segment.points.size() < 6) {
			return false;
		}

		std::vector<double>& pts = segment.points;

		double bestDistance = std::numeric_limits<double>::max();
		double bestAlong = 0.0;
		double bestProjected[3] = { 0.0, 0.0, 0.0 };

		double walked = 0.0;

		for (size_t i = 3; i + 2 < pts.size(); i += 3) {
			double a[3] = { pts[i - 3], pts[i - 2], pts[i - 1] };
			double b[3] = { pts[i], pts[i + 1], pts[i + 2] };

			double ab[3] = {
				b[0] - a[0],
				b[1] - a[1],
				b[2] - a[2]
			};

			const double abLen2 = StraightEditDot3(ab, ab);
			if (abLen2 <= 1e-12) {
				continue;
			}

			double ap[3] = {
				pickedPoint[0] - a[0],
				pickedPoint[1] - a[1],
				pickedPoint[2] - a[2]
			};

			double t = StraightEditDot3(ap, ab) / abLen2;
			t = std::max(0.0, std::min(1.0, t));

			double projected[3] = {
				a[0] + ab[0] * t,
				a[1] + ab[1] * t,
				a[2] + ab[2] * t
			};

			const double distance = StraightEditDistance3(projected, pickedPoint);
			const double edgeLength = std::sqrt(abLen2);

			if (distance < bestDistance) {
				bestDistance = distance;
				bestAlong = walked + edgeLength * t;
				bestProjected[0] = projected[0];
				bestProjected[1] = projected[1];
				bestProjected[2] = projected[2];
			}

			walked += edgeLength;
		}

		const double totalLength = StraightEditPolylineLength(pts);
		if (totalLength <= 1e-8) {
			segment.valid = false;
			return true;
		}

		const double pickTolerance = std::max(modelDiagonal * 0.03, totalLength * 0.20);
		if (bestDistance > pickTolerance) {
			return false;
		}

		const bool trimStart = bestAlong <= totalLength * 0.5;

		std::vector<double> newPts;

		if (trimStart) {
			newPts.push_back(bestProjected[0]);
			newPts.push_back(bestProjected[1]);
			newPts.push_back(bestProjected[2]);

			double walked2 = 0.0;

			for (size_t i = 3; i + 2 < pts.size(); i += 3) {
				double a[3] = { pts[i - 3], pts[i - 2], pts[i - 1] };
				double b[3] = { pts[i], pts[i + 1], pts[i + 2] };
				const double edgeLength = StraightEditDistance3(a, b);

				if (walked2 + edgeLength > bestAlong + 1e-8) {
					newPts.push_back(b[0]);
					newPts.push_back(b[1]);
					newPts.push_back(b[2]);
				}

				walked2 += edgeLength;
			}
		}
		else {
			double walked2 = 0.0;

			newPts.push_back(pts[0]);
			newPts.push_back(pts[1]);
			newPts.push_back(pts[2]);

			for (size_t i = 3; i + 2 < pts.size(); i += 3) {
				double a[3] = { pts[i - 3], pts[i - 2], pts[i - 1] };
				double b[3] = { pts[i], pts[i + 1], pts[i + 2] };
				const double edgeLength = StraightEditDistance3(a, b);

				if (walked2 + edgeLength < bestAlong - 1e-8) {
					newPts.push_back(b[0]);
					newPts.push_back(b[1]);
					newPts.push_back(b[2]);
				}
				else {
					break;
				}

				walked2 += edgeLength;
			}

			newPts.push_back(bestProjected[0]);
			newPts.push_back(bestProjected[1]);
			newPts.push_back(bestProjected[2]);
		}

		if (newPts.size() < 6 || StraightEditPolylineLength(newPts) < std::max(modelDiagonal * 0.008, 1e-6)) {
			segment.valid = false;
		}
		else {
			segment.points = newPts;
			segment.valid = true;
		}

		return true;
	}
}
namespace
{
	struct AutoStraightCageSegment
	{
		std::vector<double> points;
		bool valid = true;

		// 来源骨架 chain id。
		// 只有 sourceChainId 相同的自动 cage 才允许共享中间端面。
		int sourceChainId = -1;
	};

	static double AutoStraightDot3(const double a[3], const double b[3])
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	static double AutoStraightDistance3(const double a[3], const double b[3])
	{
		const double dx = a[0] - b[0];
		const double dy = a[1] - b[1];
		const double dz = a[2] - b[2];
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	static bool AutoStraightNormalize(double v[3])
	{
		const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

		if (n <= 1e-12) {
			return false;
		}

		v[0] /= n;
		v[1] /= n;
		v[2] /= n;

		return true;
	}

	static double AutoStraightMedian(std::vector<double> values)
	{
		if (values.empty()) {
			return 0.0;
		}

		std::sort(values.begin(), values.end());

		const size_t n = values.size();

		if (n % 2 == 1) {
			return values[n / 2];
		}

		return 0.5 * (values[n / 2 - 1] + values[n / 2]);
	}

	static double AutoStraightPercentile(std::vector<double> values, double percentile)
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

	static double AutoStraightReferenceDiagonal(vtkPolyData* refData)
	{
		if (!refData) {
			return 1.0;
		}

		double bounds[6];
		refData->GetBounds(bounds);

		const double dx = bounds[1] - bounds[0];
		const double dy = bounds[3] - bounds[2];
		const double dz = bounds[5] - bounds[4];

		const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);

		return diag > 1e-8 ? diag : 1.0;
	}

	static bool AutoStraightGetPoint(
		const std::vector<double>& pts,
		int index,
		double p[3])
	{
		const int pointCount = static_cast<int>(pts.size() / 3);

		if (index < 0 || index >= pointCount) {
			return false;
		}

		p[0] = pts[index * 3 + 0];
		p[1] = pts[index * 3 + 1];
		p[2] = pts[index * 3 + 2];

		return true;
	}

	static int AutoStraightPointCount(const std::vector<double>& pts)
	{
		return static_cast<int>(pts.size() / 3);
	}

	static double AutoStraightPolylineLength(const std::vector<double>& pts)
	{
		const int pointCount = AutoStraightPointCount(pts);

		if (pointCount < 2) {
			return 0.0;
		}

		double length = 0.0;

		for (int i = 1; i < pointCount; ++i) {
			double a[3];
			double b[3];

			AutoStraightGetPoint(pts, i - 1, a);
			AutoStraightGetPoint(pts, i, b);

			length += AutoStraightDistance3(a, b);
		}

		return length;
	}

	static bool AutoStraightPointAtDistance(
		const std::vector<double>& pts,
		double targetDistance,
		double outPoint[3])
	{
		const int pointCount = AutoStraightPointCount(pts);

		if (pointCount <= 0) {
			return false;
		}

		if (targetDistance <= 0.0) {
			return AutoStraightGetPoint(pts, 0, outPoint);
		}

		double walked = 0.0;

		for (int i = 1; i < pointCount; ++i) {
			double a[3];
			double b[3];

			AutoStraightGetPoint(pts, i - 1, a);
			AutoStraightGetPoint(pts, i, b);

			const double edgeLength = AutoStraightDistance3(a, b);

			if (edgeLength <= 1e-8) {
				continue;
			}

			if (walked + edgeLength >= targetDistance) {
				const double t = (targetDistance - walked) / edgeLength;

				for (int d = 0; d < 3; ++d) {
					outPoint[d] = a[d] * (1.0 - t) + b[d] * t;
				}

				return true;
			}

			walked += edgeLength;
		}

		return AutoStraightGetPoint(pts, pointCount - 1, outPoint);
	}

	static bool AutoStraightBuildAxisFromPoints(
		const std::vector<double>& pts,
		double outStart[3],
		double outEnd[3],
		double outAxis[3])
	{
		const int pointCount = AutoStraightPointCount(pts);

		if (pointCount < 2) {
			return false;
		}

		double first[3];
		double last[3];

		if (!AutoStraightGetPoint(pts, 0, first)) {
			return false;
		}

		if (!AutoStraightGetPoint(pts, pointCount - 1, last)) {
			return false;
		}

		outAxis[0] = last[0] - first[0];
		outAxis[1] = last[1] - first[1];
		outAxis[2] = last[2] - first[2];

		if (!AutoStraightNormalize(outAxis)) {
			return false;
		}

		double center[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < pointCount; ++i) {
			double p[3];
			AutoStraightGetPoint(pts, i, p);

			center[0] += p[0] / static_cast<double>(pointCount);
			center[1] += p[1] / static_cast<double>(pointCount);
			center[2] += p[2] / static_cast<double>(pointCount);
		}

		double minT = std::numeric_limits<double>::max();
		double maxT = -std::numeric_limits<double>::max();

		for (int i = 0; i < pointCount; ++i) {
			double p[3];
			AutoStraightGetPoint(pts, i, p);

			double v[3] = {
				p[0] - center[0],
				p[1] - center[1],
				p[2] - center[2]
			};

			const double t = AutoStraightDot3(v, outAxis);

			minT = std::min(minT, t);
			maxT = std::max(maxT, t);
		}

		if (maxT - minT <= 1e-8) {
			return false;
		}

		for (int d = 0; d < 3; ++d) {
			outStart[d] = center[d] + outAxis[d] * minT;
			outEnd[d] = center[d] + outAxis[d] * maxT;
		}

		return true;
	}

	static void AutoStraightReversePoints(std::vector<double>& pts)
	{
		const int pointCount = AutoStraightPointCount(pts);

		for (int i = 0; i < pointCount / 2; ++i) {
			const int j = pointCount - 1 - i;

			for (int d = 0; d < 3; ++d) {
				std::swap(pts[i * 3 + d], pts[j * 3 + d]);
			}
		}
	}

	static double AutoStraightEndpointDistance(
		const std::vector<double>& a,
		bool aStart,
		const std::vector<double>& b,
		bool bStart)
	{
		double pa[3];
		double pb[3];

		const int ac = AutoStraightPointCount(a);
		const int bc = AutoStraightPointCount(b);

		if (ac <= 0 || bc <= 0) {
			return std::numeric_limits<double>::max();
		}

		AutoStraightGetPoint(a, aStart ? 0 : ac - 1, pa);
		AutoStraightGetPoint(b, bStart ? 0 : bc - 1, pb);

		return AutoStraightDistance3(pa, pb);
	}

	static std::vector<AutoStraightCageSegment> AutoStraightSplitSegmentAtTurns(
		const std::vector<double>& pts,
		int sourceChainId,
		double modelDiagonal)
	{
		std::vector<AutoStraightCageSegment> result;

		const int pointCount = AutoStraightPointCount(pts);

		if (pointCount < 2) {
			return result;
		}

		const double minPieceLength = std::max(modelDiagonal * 0.012, 1e-6);
		const double turnCos = std::cos(28.0 * 3.14159265358979323846 / 180.0);

		std::vector<double> current;
		current.push_back(pts[0]);
		current.push_back(pts[1]);
		current.push_back(pts[2]);

		for (int i = 1; i < pointCount; ++i) {
			current.push_back(pts[i * 3 + 0]);
			current.push_back(pts[i * 3 + 1]);
			current.push_back(pts[i * 3 + 2]);

			if (i >= pointCount - 1) {
				continue;
			}

			double prev[3];
			double curr[3];
			double next[3];

			AutoStraightGetPoint(pts, i - 1, prev);
			AutoStraightGetPoint(pts, i, curr);
			AutoStraightGetPoint(pts, i + 1, next);

			double d0[3] = {
				curr[0] - prev[0],
				curr[1] - prev[1],
				curr[2] - prev[2]
			};

			double d1[3] = {
				next[0] - curr[0],
				next[1] - curr[1],
				next[2] - curr[2]
			};

			if (!AutoStraightNormalize(d0) || !AutoStraightNormalize(d1)) {
				continue;
			}

			const bool isTurn = AutoStraightDot3(d0, d1) < turnCos;

			if (!isTurn) {
				continue;
			}

			if (AutoStraightPolylineLength(current) >= minPieceLength) {
				AutoStraightCageSegment piece;
				piece.points = current;
				piece.valid = true;
				piece.sourceChainId = sourceChainId;

				result.push_back(piece);

				current.clear();
				current.push_back(curr[0]);
				current.push_back(curr[1]);
				current.push_back(curr[2]);
			}
		}

		if (AutoStraightPolylineLength(current) >= minPieceLength) {
			AutoStraightCageSegment piece;
			piece.points = current;
			piece.valid = true;
			piece.sourceChainId = sourceChainId;

			result.push_back(piece);
		}

		return result;
	}

	static std::vector<AutoStraightCageSegment> AutoStraightOrderAndOrientSegmentsForCage(
		const std::vector<AutoStraightCageSegment>& input,
		double modelDiagonal)
	{
		std::map<int, std::vector<AutoStraightCageSegment>> groups;

		for (const AutoStraightCageSegment& item : input) {
			groups[item.sourceChainId].push_back(item);
		}

		std::vector<AutoStraightCageSegment> ordered;
		const double connectGap = std::max(modelDiagonal * 0.035, 1e-6);

		for (auto& groupPair : groups) {
			std::vector<AutoStraightCageSegment>& group = groupPair.second;
			std::vector<bool> used(group.size(), false);

			while (true) {
				size_t seed = static_cast<size_t>(-1);

				for (size_t i = 0; i < group.size(); ++i) {
					if (!used[i]) {
						seed = i;
						break;
					}
				}

				if (seed == static_cast<size_t>(-1)) {
					break;
				}

				AutoStraightCageSegment current = group[seed];
				used[seed] = true;

				ordered.push_back(current);

				while (true) {
					size_t best = static_cast<size_t>(-1);
					bool bestReverse = false;
					double bestGap = std::numeric_limits<double>::max();

					for (size_t i = 0; i < group.size(); ++i) {
						if (used[i]) {
							continue;
						}

						const double gapStart =
							AutoStraightEndpointDistance(
								ordered.back().points,
								false,
								group[i].points,
								true);

						const double gapEnd =
							AutoStraightEndpointDistance(
								ordered.back().points,
								false,
								group[i].points,
								false);

						if (gapStart < bestGap) {
							bestGap = gapStart;
							best = i;
							bestReverse = false;
						}

						if (gapEnd < bestGap) {
							bestGap = gapEnd;
							best = i;
							bestReverse = true;
						}
					}

					if (best == static_cast<size_t>(-1) || bestGap > connectGap) {
						break;
					}

					AutoStraightCageSegment next = group[best];

					if (bestReverse) {
						AutoStraightReversePoints(next.points);
					}

					used[best] = true;
					ordered.push_back(next);
				}
			}
		}

		return ordered;
	}

	static void AutoStraightAppendWithoutDuplicate(
		std::vector<double>& dst,
		const std::vector<double>& src)
	{
		if (src.empty()) {
			return;
		}

		if (dst.empty()) {
			dst = src;
			return;
		}

		double lastDst[3] = {
			dst[dst.size() - 3],
			dst[dst.size() - 2],
			dst[dst.size() - 1]
		};

		double firstSrc[3] = {
			src[0],
			src[1],
			src[2]
		};

		int startIndex = 0;

		if (AutoStraightDistance3(lastDst, firstSrc) < 1e-6) {
			startIndex = 1;
		}

		const int pointCount = AutoStraightPointCount(src);

		for (int i = startIndex; i < pointCount; ++i) {
			dst.push_back(src[i * 3 + 0]);
			dst.push_back(src[i * 3 + 1]);
			dst.push_back(src[i * 3 + 2]);
		}
	}

	static bool AutoStraightTryMergeTwoSegments(
		AutoStraightCageSegment& a,
		AutoStraightCageSegment& b,
		double modelDiagonal)
	{
		if (!a.valid || !b.valid) {
			return false;
		}

		// 不同原始骨架 branch chain 的直线段不能合并。
		if (a.sourceChainId != b.sourceChainId) {
			return false;
		}

		if (a.points.size() < 6 || b.points.size() < 6) {
			return false;
		}

		double aStart[3], aEnd[3], aAxis[3];
		double bStart[3], bEnd[3], bAxis[3];

		if (!AutoStraightBuildAxisFromPoints(a.points, aStart, aEnd, aAxis)) {
			return false;
		}

		if (!AutoStraightBuildAxisFromPoints(b.points, bStart, bEnd, bAxis)) {
			return false;
		}

		const double cosAngle = std::abs(AutoStraightDot3(aAxis, bAxis));

		if (cosAngle < std::cos(15.0 * 3.14159265358979323846 / 180.0)) {
			return false;
		}

		const double lenA = AutoStraightPolylineLength(a.points);
		const double lenB = AutoStraightPolylineLength(b.points);

		if (lenA <= 1e-8 || lenB <= 1e-8) {
			return false;
		}

		const double mergeGap =
			std::max(modelDiagonal * 0.025, std::min(lenA, lenB) * 0.35);

		struct Candidate
		{
			int mode;
			double gap;
		};

		Candidate candidates[4] = {
			{ 0, AutoStraightDistance3(aEnd,   bStart) },
			{ 1, AutoStraightDistance3(aStart, bEnd) },
			{ 2, AutoStraightDistance3(aStart, bStart) },
			{ 3, AutoStraightDistance3(aEnd,   bEnd) }
		};

		int bestMode = -1;
		double bestGap = std::numeric_limits<double>::max();

		for (const Candidate& c : candidates) {
			if (c.gap < bestGap) {
				bestGap = c.gap;
				bestMode = c.mode;
			}
		}

		if (bestGap > mergeGap) {
			return false;
		}

		std::vector<double> merged;

		if (bestMode == 0) {
			merged = a.points;
			AutoStraightAppendWithoutDuplicate(merged, b.points);
		}
		else if (bestMode == 1) {
			merged = b.points;
			AutoStraightAppendWithoutDuplicate(merged, a.points);
		}
		else if (bestMode == 2) {
			std::vector<double> ra = a.points;
			AutoStraightReversePoints(ra);

			merged = ra;
			AutoStraightAppendWithoutDuplicate(merged, b.points);
		}
		else {
			std::vector<double> rb = b.points;
			AutoStraightReversePoints(rb);

			merged = a.points;
			AutoStraightAppendWithoutDuplicate(merged, rb);
		}

		if (merged.size() < 6) {
			return false;
		}

		a.points = merged;
		b.valid = false;

		return true;
	}

	static std::vector<AutoStraightCageSegment> AutoStraightSubdivideSegmentsForCage(
		const std::vector<AutoStraightCageSegment>& input,
		double modelDiagonal,
		vtkCellLocator* locator)
	{
		(void)modelDiagonal;
		(void)locator;

		// 当前先关闭“根据截面尺寸变化自动细分”。
		// 同一条直线段统一用平均截面建立控制笼，稳定优先。
		return input;
	}

	static std::vector<AutoStraightCageSegment> AutoStraightMergeSegmentsForCage(
		const std::vector<StraightSkeletonSegmenter::Segment>& segments,
		double modelDiagonal,
		vtkCellLocator* locator,
		bool splitAtTurns)
	{
		std::vector<AutoStraightCageSegment> work;

		for (const StraightSkeletonSegmenter::Segment& s : segments) {
			if (!s.valid) {
				continue;
			}

			if (s.points.size() < 6) {
				continue;
			}

			const double len = AutoStraightPolylineLength(s.points);

			if (len < std::max(modelDiagonal * 0.002, 1e-6)) {
				continue;
			}

			if (splitAtTurns) {
				std::vector<AutoStraightCageSegment> pieces =
					AutoStraightSplitSegmentAtTurns(
						s.points,
						s.sourceChainId,
						modelDiagonal);

				work.insert(work.end(), pieces.begin(), pieces.end());
			}
			else {
				AutoStraightCageSegment item;
				item.points = s.points;
				item.valid = true;
				item.sourceChainId = s.sourceChainId;
				work.push_back(item);
			}
		}

		bool changed = true;

		while (changed) {
			changed = false;

			for (size_t i = 0; i < work.size(); ++i) {
				if (!work[i].valid) {
					continue;
				}

				for (size_t j = i + 1; j < work.size(); ++j) {
					if (!work[j].valid) {
						continue;
					}

					if (AutoStraightTryMergeTwoSegments(
						work[i],
						work[j],
						modelDiagonal))
					{
						changed = true;
					}
				}
			}
		}

		std::vector<AutoStraightCageSegment> merged;

		for (AutoStraightCageSegment& item : work) {
			if (!item.valid) {
				continue;
			}

			const double len = AutoStraightPolylineLength(item.points);

			if (len < std::max(modelDiagonal * 0.002, 1e-6)) {
				continue;
			}

			merged.push_back(item);
		}

		std::vector<AutoStraightCageSegment> ordered =
			AutoStraightOrderAndOrientSegmentsForCage(
				merged,
				modelDiagonal);

		return AutoStraightSubdivideSegmentsForCage(
			ordered,
			modelDiagonal,
			locator);
	}

	static bool AutoStraightFindClosestVectorToSurface(
		vtkCellLocator* locator,
		const double p[3],
		double outVector[3],
		double& outDistance)
	{
		outVector[0] = 0.0;
		outVector[1] = 0.0;
		outVector[2] = 0.0;
		outDistance = 0.0;

		if (!locator) {
			return false;
		}

		double closest[3];
		vtkIdType cellId = -1;
		int subId = -1;
		double dist2 = 0.0;

		locator->FindClosestPoint(
			const_cast<double*>(p),
			closest,
			cellId,
			subId,
			dist2);

		if (cellId < 0 || dist2 <= 1e-12) {
			return false;
		}

		outVector[0] = closest[0] - p[0];
		outVector[1] = closest[1] - p[1];
		outVector[2] = closest[2] - p[2];

		outDistance = std::sqrt(dist2);

		return true;
	}

	static void AutoStraightBuildStableFrame(
		vtkCellLocator* locator,
		const double mid[3],
		const double axis[3],
		double up[3],
		double side[3])
	{
		double closestVector[3];
		double closestDistance = 0.0;

		bool hasSurfaceVector =
			AutoStraightFindClosestVectorToSurface(
				locator,
				mid,
				closestVector,
				closestDistance);

		if (hasSurfaceVector) {
			const double projection = AutoStraightDot3(closestVector, axis);

			up[0] = closestVector[0] - projection * axis[0];
			up[1] = closestVector[1] - projection * axis[1];
			up[2] = closestVector[2] - projection * axis[2];

			if (!AutoStraightNormalize(up)) {
				hasSurfaceVector = false;
			}
		}

		if (!hasSurfaceVector) {
			double worldUp[3] = { 0.0, 0.0, 1.0 };

			if (std::abs(AutoStraightDot3(worldUp, axis)) > 0.90) {
				worldUp[0] = 0.0;
				worldUp[1] = 1.0;
				worldUp[2] = 0.0;
			}

			const double projection = AutoStraightDot3(worldUp, axis);

			up[0] = worldUp[0] - projection * axis[0];
			up[1] = worldUp[1] - projection * axis[1];
			up[2] = worldUp[2] - projection * axis[2];

			AutoStraightNormalize(up);
		}

		vtkMath::Cross(axis, up, side);
		AutoStraightNormalize(side);

		vtkMath::Cross(side, axis, up);
		AutoStraightNormalize(up);
	}

	static bool AutoStraightRayDistance(
		vtkCellLocator* locator,
		const double p[3],
		const double dir[3],
		double farDistance,
		double& outDistance)
	{
		outDistance = 0.0;

		if (!locator || farDistance <= 1e-8) {
			return false;
		}

		double p1[3] = {
			p[0] + dir[0] * farDistance,
			p[1] + dir[1] * farDistance,
			p[2] + dir[2] * farDistance
		};

		double t = 0.0;
		double x[3] = { 0.0, 0.0, 0.0 };
		double pcoords[3] = { 0.0, 0.0, 0.0 };
		int subId = -1;
		vtkIdType cellId = -1;

		const int hit = locator->IntersectWithLine(
			const_cast<double*>(p),
			p1,
			1e-6,
			t,
			x,
			pcoords,
			subId,
			cellId);

		if (!hit || cellId < 0) {
			return false;
		}

		outDistance = AutoStraightDistance3(p, x);

		return outDistance > 1e-8;
	}

	static double AutoStraightEstimateSkeletonRadiusForSegment(
		const std::vector<double>& skeletonNodes,
		const std::vector<double>& skeletonRadii,
		const double start[3],
		const double end[3],
		double modelDiagonal)
	{
		const int nodeCount =
			std::min(
				static_cast<int>(skeletonNodes.size() / 3),
				static_cast<int>(skeletonRadii.size()));

		if (nodeCount <= 0) {
			return 0.0;
		}

		double axis[3] = {
			end[0] - start[0],
			end[1] - start[1],
			end[2] - start[2]
		};

		const double length = AutoStraightDistance3(start, end);
		if (length <= 1e-8 || !AutoStraightNormalize(axis)) {
			return 0.0;
		}

		const double axialMargin =
			std::max(modelDiagonal * 0.008, length * 0.10);
		const double perpendicularLimit =
			std::max(modelDiagonal * 0.018, length * 0.18);

		std::vector<double> values;

		for (int i = 0; i < nodeCount; ++i) {
			const double radius = skeletonRadii[i];
			if (radius <= 1e-8) {
				continue;
			}

			double p[3] = {
				skeletonNodes[3 * i],
				skeletonNodes[3 * i + 1],
				skeletonNodes[3 * i + 2]
			};

			double v[3] = {
				p[0] - start[0],
				p[1] - start[1],
				p[2] - start[2]
			};

			const double t = AutoStraightDot3(v, axis);
			if (t < -axialMargin || t > length + axialMargin) {
				continue;
			}

			double projected[3] = {
				start[0] + axis[0] * t,
				start[1] + axis[1] * t,
				start[2] + axis[2] * t
			};

			if (AutoStraightDistance3(p, projected) > perpendicularLimit) {
				continue;
			}

			values.push_back(radius);
		}

		if (!values.empty()) {
			return AutoStraightMedian(values);
		}

		double mid[3] = {
			0.5 * (start[0] + end[0]),
			0.5 * (start[1] + end[1]),
			0.5 * (start[2] + end[2])
		};

		double bestDistance2 = std::numeric_limits<double>::max();
		double bestRadius = 0.0;

		for (int i = 0; i < nodeCount; ++i) {
			const double radius = skeletonRadii[i];
			if (radius <= 1e-8) {
				continue;
			}

			double p[3] = {
				skeletonNodes[3 * i],
				skeletonNodes[3 * i + 1],
				skeletonNodes[3 * i + 2]
			};

			const double d2 = vtkMath::Distance2BetweenPoints(mid, p);
			if (d2 < bestDistance2) {
				bestDistance2 = d2;
				bestRadius = radius;
			}
		}

		return bestRadius;
	}

	static bool AutoStraightEstimateSupportHalfExtentsAtPoint(
		const std::vector<double>& skeletonNodes,
		const std::vector<unsigned int>& skeletonEdges,
		const std::vector<unsigned int>& supportOffsets,
		const std::vector<double>& supportPoints,
		const double center[3],
		const double segmentStart[3],
		const double segmentEnd[3],
		const double axis[3],
		double modelDiagonal,
		const double up[3],
		const double side[3],
		double& halfWidth,
		double& halfDepth)
	{
		const int nodeCount = static_cast<int>(skeletonNodes.size() / 3);
		if (nodeCount <= 0 ||
			supportOffsets.size() < static_cast<size_t>(nodeCount + 1) ||
			supportPoints.size() < 3)
		{
			return false;
		}

		std::vector<int> degree(nodeCount, 0);
		for (size_t ei = 0; ei + 1 < skeletonEdges.size(); ei += 2) {
			const int a = static_cast<int>(skeletonEdges[ei]);
			const int b = static_cast<int>(skeletonEdges[ei + 1]);
			if (a >= 0 && a < nodeCount && b >= 0 && b < nodeCount && a != b) {
				++degree[a];
				++degree[b];
			}
		}

		const double segmentLength =
			AutoStraightDistance3(segmentStart, segmentEnd);

		const bool hasAxis =
			segmentLength > 1e-8 &&
			std::abs(AutoStraightDot3(axis, axis) - 1.0) < 1e-4;

		int bestNode = -1;
		double bestDistance2 = std::numeric_limits<double>::max();
		int bestNodeAny = -1;
		double bestDistance2Any = std::numeric_limits<double>::max();

		for (int i = 0; i < nodeCount; ++i) {
			const unsigned int begin = supportOffsets[i];
			const unsigned int end = supportOffsets[i + 1];

			if (end <= begin || static_cast<size_t>(end) * 3 > supportPoints.size()) {
				continue;
			}

			double p[3] = {
				skeletonNodes[3 * i],
				skeletonNodes[3 * i + 1],
				skeletonNodes[3 * i + 2]
			};

			double score = vtkMath::Distance2BetweenPoints(center, p);

			if (hasAxis) {
				double vStart[3] = {
					p[0] - segmentStart[0],
					p[1] - segmentStart[1],
					p[2] - segmentStart[2]
				};

				const double t = AutoStraightDot3(vStart, axis);

				double projected[3] = {
					segmentStart[0] + axis[0] * t,
					segmentStart[1] + axis[1] * t,
					segmentStart[2] + axis[2] * t
				};

				const double offAxis2 = vtkMath::Distance2BetweenPoints(p, projected);
				const double outside =
					t < 0.0 ? -t :
					t > segmentLength ? t - segmentLength :
					0.0;

				score =
					offAxis2 +
					outside * outside * 9.0 +
					vtkMath::Distance2BetweenPoints(center, projected) * 0.04;
			}

			if (score < bestDistance2Any) {
				bestDistance2Any = score;
				bestNodeAny = i;
			}

			if (degree[i] > 2) {
				continue;
			}

			if (score < bestDistance2) {
				bestDistance2 = score;
				bestNode = i;
			}
		}

		if (bestNode < 0) {
			bestNode = bestNodeAny;
		}

		if (bestNode < 0) {
			return false;
		}

		const unsigned int begin = supportOffsets[bestNode];
		const unsigned int end = supportOffsets[bestNode + 1];

		std::vector<double> sideValues;
		std::vector<double> upValues;
		sideValues.reserve(end - begin);
		upValues.reserve(end - begin);

		for (unsigned int pi = begin; pi < end; ++pi) {
			if (static_cast<size_t>(pi) * 3 + 2 >= supportPoints.size()) {
				continue;
			}

			double q[3] = {
				supportPoints[3 * pi],
				supportPoints[3 * pi + 1],
				supportPoints[3 * pi + 2]
			};

			double v[3] = {
				q[0] - center[0],
				q[1] - center[1],
				q[2] - center[2]
			};

			if (hasAxis) {
				const double axial = std::abs(AutoStraightDot3(v, axis));
				const double sideProjection = AutoStraightDot3(v, side);
				const double upProjection = AutoStraightDot3(v, up);
				const double radiusHint =
					std::sqrt(
						sideProjection * sideProjection +
						upProjection * upProjection);

				if (axial > std::max(modelDiagonal * 0.004, radiusHint * 1.25)) {
					continue;
				}
			}

			sideValues.push_back(std::abs(AutoStraightDot3(v, side)));
			upValues.push_back(std::abs(AutoStraightDot3(v, up)));
		}

		if (sideValues.size() < 4 || upValues.size() < 4) {
			return false;
		}

		halfWidth = AutoStraightPercentile(sideValues, 0.92) * 1.02;
		halfDepth = AutoStraightPercentile(upValues, 0.92) * 1.02;

		return halfWidth > 1e-8 && halfDepth > 1e-8;
	}

	static bool AutoStraightEstimateBroadSupportHalfExtentsAtPoint(
		const std::vector<double>& skeletonNodes,
		const std::vector<unsigned int>& skeletonEdges,
		const std::vector<unsigned int>& supportOffsets,
		const std::vector<double>& supportPoints,
		const double center[3],
		const double segmentStart[3],
		const double segmentEnd[3],
		const double axis[3],
		double modelDiagonal,
		const double up[3],
		const double side[3],
		double& halfWidth,
		double& halfDepth)
	{
		const int nodeCount = static_cast<int>(skeletonNodes.size() / 3);
		if (nodeCount <= 0 ||
			supportOffsets.size() < static_cast<size_t>(nodeCount + 1) ||
			supportPoints.size() < 3)
		{
			return false;
		}

		std::vector<int> degree(nodeCount, 0);
		for (size_t ei = 0; ei + 1 < skeletonEdges.size(); ei += 2) {
			const int a = static_cast<int>(skeletonEdges[ei]);
			const int b = static_cast<int>(skeletonEdges[ei + 1]);
			if (a >= 0 && a < nodeCount && b >= 0 && b < nodeCount && a != b) {
				++degree[a];
				++degree[b];
			}
		}

		const double segmentLength =
			AutoStraightDistance3(segmentStart, segmentEnd);

		const bool hasAxis =
			segmentLength > 1e-8 &&
			std::abs(AutoStraightDot3(axis, axis) - 1.0) < 1e-4;

		struct CandidateNode
		{
			int index = -1;
			double score = 0.0;
			bool branch = false;
		};

		std::vector<CandidateNode> candidates;
		candidates.reserve(static_cast<size_t>(nodeCount));

		for (int i = 0; i < nodeCount; ++i) {
			const unsigned int begin = supportOffsets[i];
			const unsigned int end = supportOffsets[i + 1];

			if (end <= begin || static_cast<size_t>(end) * 3 > supportPoints.size()) {
				continue;
			}

			double p[3] = {
				skeletonNodes[3 * i],
				skeletonNodes[3 * i + 1],
				skeletonNodes[3 * i + 2]
			};

			double score = vtkMath::Distance2BetweenPoints(center, p);

			if (hasAxis) {
				double vStart[3] = {
					p[0] - segmentStart[0],
					p[1] - segmentStart[1],
					p[2] - segmentStart[2]
				};

				const double t = AutoStraightDot3(vStart, axis);

				double projected[3] = {
					segmentStart[0] + axis[0] * t,
					segmentStart[1] + axis[1] * t,
					segmentStart[2] + axis[2] * t
				};

				const double offAxis2 = vtkMath::Distance2BetweenPoints(p, projected);
				const double outside =
					t < 0.0 ? -t :
					t > segmentLength ? t - segmentLength :
					0.0;

				score =
					offAxis2 +
					outside * outside * 9.0 +
					vtkMath::Distance2BetweenPoints(center, projected) * 0.06;
			}

			CandidateNode c;
			c.index = i;
			c.score = score;
			c.branch = degree[i] > 2;
			candidates.push_back(c);
		}

		if (candidates.empty()) {
			return false;
		}

		std::sort(
			candidates.begin(),
			candidates.end(),
			[](const CandidateNode& a, const CandidateNode& b) {
			return a.score < b.score;
		});

		const double bestScore = std::max(candidates.front().score, 1e-12);
		const double collectRadius =
			std::max(modelDiagonal * 0.026, segmentLength * 0.42);
		const double collectScoreLimit =
			std::max(bestScore * 7.0, collectRadius * collectRadius);

		std::vector<double> sideValues;
		std::vector<double> upValues;
		sideValues.reserve(256);
		upValues.reserve(256);

		int usedNodes = 0;
		int usedNonBranchNodes = 0;

		for (const CandidateNode& candidate : candidates) {
			if (candidate.index < 0) {
				continue;
			}

			if (candidate.branch && usedNonBranchNodes >= 3) {
				continue;
			}

			if (usedNodes >= 12 && candidate.score > collectScoreLimit) {
				break;
			}

			if (candidate.score > collectScoreLimit && usedNodes >= 4) {
				continue;
			}

			const int nodeIndex = candidate.index;
			const unsigned int begin = supportOffsets[nodeIndex];
			const unsigned int end = supportOffsets[nodeIndex + 1];
			bool thisNodeContributed = false;

			for (unsigned int pi = begin; pi < end; ++pi) {
				if (static_cast<size_t>(pi) * 3 + 2 >= supportPoints.size()) {
					continue;
				}

				double q[3] = {
					supportPoints[3 * pi],
					supportPoints[3 * pi + 1],
					supportPoints[3 * pi + 2]
				};

				double v[3] = {
					q[0] - center[0],
					q[1] - center[1],
					q[2] - center[2]
				};

				if (hasAxis) {
					const double axial = std::abs(AutoStraightDot3(v, axis));
					const double sideProjection = AutoStraightDot3(v, side);
					const double upProjection = AutoStraightDot3(v, up);
					const double radiusHint =
						std::sqrt(
							sideProjection * sideProjection +
							upProjection * upProjection);

					const double axialLimit =
						std::max(
							std::max(modelDiagonal * 0.014, segmentLength * 0.24),
							radiusHint * 1.70);

					if (axial > axialLimit) {
						continue;
					}
				}

				const double sideValue = std::abs(AutoStraightDot3(v, side));
				const double upValue = std::abs(AutoStraightDot3(v, up));

				if (sideValue <= 1e-9 && upValue <= 1e-9) {
					continue;
				}

				sideValues.push_back(sideValue);
				upValues.push_back(upValue);
				thisNodeContributed = true;
			}

			if (thisNodeContributed) {
				++usedNodes;
				if (!candidate.branch) {
					++usedNonBranchNodes;
				}
			}
		}

		if (sideValues.size() < 4 || upValues.size() < 4) {
			return false;
		}

		halfWidth = AutoStraightPercentile(sideValues, 0.94) * 1.06;
		halfDepth = AutoStraightPercentile(upValues, 0.94) * 1.06;

		return halfWidth > 1e-8 && halfDepth > 1e-8;
	}

	static double AutoStraightEstimateHalfExtentInDirection(
		vtkCellLocator* locator,
		const double start[3],
		const double end[3],
		const double dir[3],
		double modelDiagonal,
		double fallbackHalfSize,
		double minHalfSize,
		double maxHalfSize)
	{
		std::vector<double> values;

		double negDir[3] = {
			-dir[0],
			-dir[1],
			-dir[2]
		};

		for (int i = 0; i < 9; ++i) {
			const double t = 0.18 + 0.64 * static_cast<double>(i) / 8.0;

			double p[3] = {
				start[0] * (1.0 - t) + end[0] * t,
				start[1] * (1.0 - t) + end[1] * t,
				start[2] * (1.0 - t) + end[2] * t
			};

			double dPlus = 0.0;
			double dMinus = 0.0;

			const bool hitPlus =
				AutoStraightRayDistance(
					locator,
					p,
					dir,
					modelDiagonal,
					dPlus);

			const bool hitMinus =
				AutoStraightRayDistance(
					locator,
					p,
					negDir,
					modelDiagonal,
					dMinus);

			if (hitPlus && hitMinus) {
				const double nearSide = std::min(dPlus, dMinus);
				const double farSide = std::max(dPlus, dMinus);

				if (nearSide > 1e-8 && farSide / nearSide > 2.2) {
					values.push_back(std::min(farSide, nearSide * 1.90));
				}
				else {
					values.push_back(0.5 * (dPlus + dMinus));
				}
			}
			else if (hitPlus) {
				values.push_back(dPlus);
			}
			else if (hitMinus) {
				values.push_back(dMinus);
			}
		}

		double halfExtent = fallbackHalfSize;

		if (!values.empty()) {
			halfExtent = AutoStraightMedian(values);
		}

		halfExtent = std::max(minHalfSize, halfExtent);
		halfExtent = std::min(maxHalfSize, halfExtent);

		return halfExtent;
	}

	static double AutoStraightClampDoubleLocal(double value, double lo, double hi)
	{
		if (value < lo) return lo;
		if (value > hi) return hi;
		return value;
	}

	static bool AutoStraightFitAxisToBoundary(
		vtkCellLocator* locator,
		double center[3],
		const double dir[3],
		double modelDiagonal,
		double& halfExtent,
		double& maxHalfExtent)
	{
		if (!locator || modelDiagonal <= 1e-8) {
			return false;
		}

		double negDir[3] = { -dir[0], -dir[1], -dir[2] };
		double dPlus = 0.0;
		double dMinus = 0.0;
		const double farDistance = std::max(modelDiagonal * 1.5, 1e-6);

		const bool hitPlus =
			AutoStraightRayDistance(locator, center, dir, farDistance, dPlus);
		const bool hitMinus =
			AutoStraightRayDistance(locator, center, negDir, farDistance, dMinus);

		if (hitPlus && hitMinus) {
			halfExtent = std::max(
				halfExtent,
				std::max(dPlus, dMinus) * 1.005);
			maxHalfExtent = std::max(maxHalfExtent, halfExtent * 1.08);
			return halfExtent > 1e-8;
		}

		if (hitPlus) {
			halfExtent = std::max(halfExtent, dPlus * 1.005);
			maxHalfExtent = std::max(maxHalfExtent, halfExtent * 1.08);
			return halfExtent > 1e-8;
		}

		if (hitMinus) {
			halfExtent = std::max(halfExtent, dMinus * 1.005);
			maxHalfExtent = std::max(maxHalfExtent, halfExtent * 1.08);
			return halfExtent > 1e-8;
		}

		return false;
	}

	static bool AutoStraightEstimateSectionHalfExtentsAtPoint(
		vtkCellLocator* locator,
		const std::vector<double>& skeletonNodes,
		const std::vector<unsigned int>& skeletonEdges,
		const std::vector<unsigned int>& supportOffsets,
		const std::vector<double>& supportPoints,
		const double center[3],
		const double segmentStart[3],
		const double segmentEnd[3],
		const double axis[3],
		double modelDiagonal,
		const double up[3],
		const double side[3],
		double fallbackHalfSize,
		double minHalfSize,
		double maxHalfSize,
		double& halfWidth,
		double& halfDepth,
		bool& usedSupport)
	{
		usedSupport = false;

		double rayWidth =
			AutoStraightEstimateHalfExtentInDirection(
				locator,
				center,
				center,
				side,
				modelDiagonal,
				fallbackHalfSize,
				minHalfSize,
				maxHalfSize);

		double rayDepth =
			AutoStraightEstimateHalfExtentInDirection(
				locator,
				center,
				center,
				up,
				modelDiagonal,
				fallbackHalfSize,
				minHalfSize,
				maxHalfSize);

		double supportWidth = 0.0;
		double supportDepth = 0.0;

		const bool hasSupport =
			AutoStraightEstimateBroadSupportHalfExtentsAtPoint(
				skeletonNodes,
				skeletonEdges,
				supportOffsets,
				supportPoints,
				center,
				segmentStart,
				segmentEnd,
				axis,
				modelDiagonal,
				up,
				side,
				supportWidth,
				supportDepth);

		if (hasSupport) {
			usedSupport = true;
			halfWidth = std::max(rayWidth, supportWidth);
			halfDepth = std::max(rayDepth, supportDepth);
		}
		else {
			halfWidth = rayWidth;
			halfDepth = rayDepth;
		}

		halfWidth = AutoStraightClampDoubleLocal(halfWidth, minHalfSize, maxHalfSize);
		halfDepth = AutoStraightClampDoubleLocal(halfDepth, minHalfSize, maxHalfSize);

		return halfWidth > 1e-8 && halfDepth > 1e-8;
	}

	static void AutoStraightMakeRing(
		const double center[3],
		const double up[3],
		const double side[3],
		double halfWidth,
		double halfDepth,
		double ring[4][3])
	{
		for (int d = 0; d < 3; ++d) {
			ring[0][d] = center[d] - side[d] * halfWidth - up[d] * halfDepth;
			ring[1][d] = center[d] + side[d] * halfWidth - up[d] * halfDepth;
			ring[2][d] = center[d] + side[d] * halfWidth + up[d] * halfDepth;
			ring[3][d] = center[d] - side[d] * halfWidth + up[d] * halfDepth;
		}
	}


	static bool AutoStraightPointInsideClosedSurface(
		vtkPolyData* surface,
		const double p[3],
		double tolerance = 1e-6)
	{
		if (!surface || surface->GetNumberOfPoints() == 0 || surface->GetNumberOfCells() == 0) {
			return false;
		}

		double b[6];
		surface->GetBounds(b);
		const double diag = std::sqrt(
			(b[1] - b[0]) * (b[1] - b[0]) +
			(b[3] - b[2]) * (b[3] - b[2]) +
			(b[5] - b[4]) * (b[5] - b[4]));
		const double eps = std::max(diag * 1e-6, tolerance);
		if (p[0] < b[0] - eps || p[0] > b[1] + eps ||
			p[1] < b[2] - eps || p[1] > b[3] + eps ||
			p[2] < b[4] - eps || p[2] > b[5] + eps)
		{
			return false;
		}

		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		pts->InsertNextPoint(p);
		vtkSmartPointer<vtkPolyData> query = vtkSmartPointer<vtkPolyData>::New();
		query->SetPoints(pts);

		vtkSmartPointer<vtkSelectEnclosedPoints> selector =
			vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selector->SetInputData(query);
		selector->SetSurfaceData(surface);
		selector->SetTolerance(tolerance);
		selector->Update();
		return selector->IsInside(0) != 0;
	}

	static bool AutoStraightRayDistanceToSurface(
		vtkOBBTree* tree,
		const double origin[3],
		const double dirIn[3],
		double maxLen,
		double& outDist)
	{
		outDist = 0.0;
		if (!tree || maxLen <= 1e-9) return false;

		double dir[3] = { dirIn[0], dirIn[1], dirIn[2] };
		if (vtkMath::Normalize(dir) < 1e-10) return false;

		double p1[3] = { origin[0], origin[1], origin[2] };
		double p2[3] = {
			origin[0] + dir[0] * maxLen,
			origin[1] + dir[1] * maxLen,
			origin[2] + dir[2] * maxLen
		};

		vtkSmartPointer<vtkPoints> hits = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
		if (tree->IntersectWithLine(p1, p2, hits, ids) == 0 || hits->GetNumberOfPoints() == 0) {
			return false;
		}

		double best = std::numeric_limits<double>::max();
		for (vtkIdType i = 0; i < hits->GetNumberOfPoints(); ++i) {
			double h[3];
			hits->GetPoint(i, h);
			double d = std::sqrt(vtkMath::Distance2BetweenPoints(origin, h));
			if (d > maxLen * 1e-6) best = std::min(best, d);
		}
		if (!std::isfinite(best) || best == std::numeric_limits<double>::max()) {
			return false;
		}
		outDist = best;
		return true;
	}

	static bool AutoStraightRayFarthestDistanceToSurface(
		vtkOBBTree* tree,
		const double origin[3],
		const double dirIn[3],
		double maxLen,
		double& outDist)
	{
		outDist = 0.0;
		if (!tree || maxLen <= 1e-9) return false;

		double dir[3] = { dirIn[0], dirIn[1], dirIn[2] };
		if (vtkMath::Normalize(dir) < 1e-10) return false;

		double end[3] = {
			origin[0] + dir[0] * maxLen,
			origin[1] + dir[1] * maxLen,
			origin[2] + dir[2] * maxLen
		};
		vtkSmartPointer<vtkPoints> hits = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
		if (tree->IntersectWithLine(origin, end, hits, ids) == 0 || hits->GetNumberOfPoints() == 0) {
			return false;
		}

		double farthest = 0.0;
		for (vtkIdType i = 0; i < hits->GetNumberOfPoints(); ++i) {
			double hit[3];
			hits->GetPoint(i, hit);
			double along[3] = { hit[0] - origin[0], hit[1] - origin[1], hit[2] - origin[2] };
			double distance = vtkMath::Dot(along, dir);
			if (distance > maxLen * 1e-6 && distance <= maxLen * (1.0 + 1e-8)) {
				farthest = std::max(farthest, distance);
			}
		}
		if (farthest <= 0.0) return false;
		outDist = farthest;
		return true;
	}

	static bool AutoStraightClampRingToNonDesignSurface(
		vtkPolyData* surface,
		const double center[3],
		const double up[3],
		const double side[3],
		double modelDiagonal,
		double& halfWidth,
		double& halfDepth,
		double ring[4][3])
	{
		if (!surface || surface->GetNumberOfPoints() == 0 || surface->GetNumberOfCells() == 0) {
			return false;
		}

		bool touchesNonDesign = AutoStraightPointInsideClosedSurface(surface, center, 1e-6);
		for (int i = 0; i < 4 && !touchesNonDesign; ++i) {
			if (AutoStraightPointInsideClosedSurface(surface, ring[i], 1e-6)) {
				touchesNonDesign = true;
			}
		}
		if (!touchesNonDesign) {
			return false;
		}

		vtkSmartPointer<vtkOBBTree> tree = vtkSmartPointer<vtkOBBTree>::New();
		tree->SetDataSet(surface);
		tree->BuildLocator();

		const double maxLen = std::max(modelDiagonal * 2.0, std::max(halfWidth, halfDepth) * 20.0);
		double sideNeg[3] = { -side[0], -side[1], -side[2] };
		double upNeg[3] = { -up[0], -up[1], -up[2] };

		double dsP = 0.0, dsN = 0.0, duP = 0.0, duN = 0.0;
		bool hsP = AutoStraightRayDistanceToSurface(tree, center, side, maxLen, dsP);
		bool hsN = AutoStraightRayDistanceToSurface(tree, center, sideNeg, maxLen, dsN);
		bool huP = AutoStraightRayDistanceToSurface(tree, center, up, maxLen, duP);
		bool huN = AutoStraightRayDistanceToSurface(tree, center, upNeg, maxLen, duN);

		const double oldW = halfWidth;
		const double oldD = halfDepth;
		const double safety = 0.88;

		if (hsP && hsN) {
			halfWidth = std::min(halfWidth, std::max(1e-6, safety * std::min(dsP, dsN)));
		}
		if (huP && huN) {
			halfDepth = std::min(halfDepth, std::max(1e-6, safety * std::min(duP, duN)));
		}

		const bool changed =
			(halfWidth < oldW * 0.995) ||
			(halfDepth < oldD * 0.995);

		if (changed) {
			AutoStraightMakeRing(center, up, side, halfWidth, halfDepth, ring);
		}
		return changed;
	}

	struct AutoStraightCageBuildItem
	{
		CageSegment cage;

		double startCenter[3];
		double endCenter[3];

		double axis[3];
		double halfWidth = 0.0;
		double halfDepth = 0.0;

		double length = 0.0;

		// 来源骨架 chain id。
		// 共享端面必须要求这个 id 相同。
		int sourceChainId = -1;
	};

	static void AutoStraightRingCenter(
		const double ring[4][3],
		double center[3])
	{
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;

		for (int i = 0; i < 4; ++i) {
			center[0] += ring[i][0] * 0.25;
			center[1] += ring[i][1] * 0.25;
			center[2] += ring[i][2] * 0.25;
		}
	}

	static void AutoStraightCopyRingLocal(
		double dst[4][3],
		const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	}

	static void AutoStraightGetEndpointRing(
		const AutoStraightCageBuildItem& item,
		bool isStart,
		double ring[4][3])
	{
		if (isStart) {
			AutoStraightCopyRingLocal(ring, item.cage.pStart);
		}
		else {
			AutoStraightCopyRingLocal(ring, item.cage.pEnd);
		}
	}

	static void AutoStraightSetEndpointRing(
		AutoStraightCageBuildItem& item,
		bool isStart,
		const double ring[4][3])
	{
		if (isStart) {
			AutoStraightCopyRingLocal(item.cage.pStart, ring);
		}
		else {
			AutoStraightCopyRingLocal(item.cage.pEnd, ring);
		}
	}

	static void AutoStraightGetEndpointCenter(
		const AutoStraightCageBuildItem& item,
		bool isStart,
		double center[3])
	{
		if (isStart) {
			center[0] = item.startCenter[0];
			center[1] = item.startCenter[1];
			center[2] = item.startCenter[2];
		}
		else {
			center[0] = item.endCenter[0];
			center[1] = item.endCenter[1];
			center[2] = item.endCenter[2];
		}
	}

	static void AutoStraightSetEndpointCenter(
		AutoStraightCageBuildItem& item,
		bool isStart,
		const double center[3])
	{
		if (isStart) {
			item.startCenter[0] = center[0];
			item.startCenter[1] = center[1];
			item.startCenter[2] = center[2];
		}
		else {
			item.endCenter[0] = center[0];
			item.endCenter[1] = center[1];
			item.endCenter[2] = center[2];
		}
	}

	static double AutoStraightRingAverageRadius(const double ring[4][3])
	{
		double center[3];
		AutoStraightRingCenter(ring, center);

		double radius = 0.0;

		for (int i = 0; i < 4; ++i) {
			radius += AutoStraightDistance3(center, ring[i]) * 0.25;
		}

		return radius;
	}

	static void AutoStraightScaleEndpointRing(
		AutoStraightCageBuildItem& item,
		bool isStart,
		double scale)
	{
		double ring[4][3];
		AutoStraightGetEndpointRing(item, isStart, ring);

		double center[3];
		AutoStraightRingCenter(ring, center);

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				ring[i][d] = center[d] + (ring[i][d] - center[d]) * scale;
			}
		}

		AutoStraightSetEndpointRing(item, isStart, ring);
	}

	static void AutoStraightResolveCrossChainEndpointConflicts(
		std::vector<AutoStraightCageBuildItem>& items,
		double modelDiagonal)
	{
		if (items.size() < 2) {
			return;
		}

		struct EndpointRef
		{
			size_t itemIndex;
			bool isStart;
		};

		std::vector<EndpointRef> endpoints;

		for (size_t i = 0; i < items.size(); ++i) {
			endpoints.push_back({ i, true });
			endpoints.push_back({ i, false });
		}

		auto ExtractRingFrame =
			[](const double ring[4][3], double side[3], double up[3]) -> bool
		{
			double leftMid[3] = {
				0.5 * (ring[0][0] + ring[3][0]),
				0.5 * (ring[0][1] + ring[3][1]),
				0.5 * (ring[0][2] + ring[3][2])
			};

			double rightMid[3] = {
				0.5 * (ring[1][0] + ring[2][0]),
				0.5 * (ring[1][1] + ring[2][1]),
				0.5 * (ring[1][2] + ring[2][2])
			};

			double bottomMid[3] = {
				0.5 * (ring[0][0] + ring[1][0]),
				0.5 * (ring[0][1] + ring[1][1]),
				0.5 * (ring[0][2] + ring[1][2])
			};

			double topMid[3] = {
				0.5 * (ring[2][0] + ring[3][0]),
				0.5 * (ring[2][1] + ring[3][1]),
				0.5 * (ring[2][2] + ring[3][2])
			};

			side[0] = rightMid[0] - leftMid[0];
			side[1] = rightMid[1] - leftMid[1];
			side[2] = rightMid[2] - leftMid[2];

			up[0] = topMid[0] - bottomMid[0];
			up[1] = topMid[1] - bottomMid[1];
			up[2] = topMid[2] - bottomMid[2];

			if (!AutoStraightNormalize(side)) {
				return false;
			}

			const double upProjection = AutoStraightDot3(up, side);

			up[0] -= upProjection * side[0];
			up[1] -= upProjection * side[1];
			up[2] -= upProjection * side[2];

			if (!AutoStraightNormalize(up)) {
				return false;
			}

			return true;
		};

		auto RingHalfExtentAlong =
			[](const double ring[4][3], const double center[3], const double dir[3]) -> double
		{
			double value = 0.0;

			for (int i = 0; i < 4; ++i) {
				double v[3] = {
					ring[i][0] - center[0],
					ring[i][1] - center[1],
					ring[i][2] - center[2]
				};

				value = std::max(value, std::abs(AutoStraightDot3(v, dir)));
			}

			return value;
		};

		auto RingSupportAlong =
			[](const double ring[4][3], const double center[3], const double dir[3]) -> double
		{
			double support = 0.0;

			for (int i = 0; i < 4; ++i) {
				double v[3] = {
					ring[i][0] - center[0],
					ring[i][1] - center[1],
					ring[i][2] - center[2]
				};

				support = std::max(support, std::abs(AutoStraightDot3(v, dir)));
			}

			return support;
		};

		auto EndpointNormal =
			[](const AutoStraightCageBuildItem& item, bool isStart, double n[3])
		{
			if (isStart) {
				n[0] = -item.axis[0];
				n[1] = -item.axis[1];
				n[2] = -item.axis[2];
			}
			else {
				n[0] = item.axis[0];
				n[1] = item.axis[1];
				n[2] = item.axis[2];
			}

			AutoStraightNormalize(n);
		};

		auto ShouldResolveRealOverlap =
			[&](
				const AutoStraightCageBuildItem& itemA,
				bool aStart,
				const AutoStraightCageBuildItem& itemB,
				bool bStart,
				double detectPadding) -> bool
		{
			double ringA[4][3];
			double ringB[4][3];

			AutoStraightGetEndpointRing(itemA, aStart, ringA);
			AutoStraightGetEndpointRing(itemB, bStart, ringB);

			double centerA[3];
			double centerB[3];

			AutoStraightRingCenter(ringA, centerA);
			AutoStraightRingCenter(ringB, centerB);

			double deltaAB[3] = {
				centerB[0] - centerA[0],
				centerB[1] - centerA[1],
				centerB[2] - centerA[2]
			};

			double sideA[3];
			double upA[3];
			double sideB[3];
			double upB[3];

			if (!ExtractRingFrame(ringA, sideA, upA) ||
				!ExtractRingFrame(ringB, sideB, upB))
			{
				return false;
			}

			double normalA[3];
			double normalB[3];

			EndpointNormal(itemA, aStart, normalA);
			EndpointNormal(itemB, bStart, normalB);

			// ============================================================
			// 关键：
			// 不再只看中心距离 + 平均半径。
			// 必须在 A 的 side/up 平面方向上发生重叠，
			// 同时也必须在 B 的 side/up 平面方向上发生重叠。
			// 否则说明只是附近端点，不是真正碰撞。
			// ============================================================

			const double aSide =
				RingHalfExtentAlong(ringA, centerA, sideA);
			const double aUp =
				RingHalfExtentAlong(ringA, centerA, upA);

			const double bSupportOnASide =
				RingSupportAlong(ringB, centerB, sideA);
			const double bSupportOnAUp =
				RingSupportAlong(ringB, centerB, upA);

			const double gapOnASide =
				std::abs(AutoStraightDot3(deltaAB, sideA));

			const double gapOnAUp =
				std::abs(AutoStraightDot3(deltaAB, upA));

			const bool overlapInAFrame =
				gapOnASide <= aSide + bSupportOnASide + detectPadding &&
				gapOnAUp <= aUp + bSupportOnAUp + detectPadding;

			if (!overlapInAFrame) {
				return false;
			}

			double deltaBA[3] = {
				-centerB[0] + centerA[0],
				-centerB[1] + centerA[1],
				-centerB[2] + centerA[2]
			};

			const double bSide =
				RingHalfExtentAlong(ringB, centerB, sideB);
			const double bUp =
				RingHalfExtentAlong(ringB, centerB, upB);

			const double aSupportOnBSide =
				RingSupportAlong(ringA, centerA, sideB);
			const double aSupportOnBUp =
				RingSupportAlong(ringA, centerA, upB);

			const double gapOnBSide =
				std::abs(AutoStraightDot3(deltaBA, sideB));

			const double gapOnBUp =
				std::abs(AutoStraightDot3(deltaBA, upB));

			const bool overlapInBFrame =
				gapOnBSide <= bSide + aSupportOnBSide + detectPadding &&
				gapOnBUp <= bUp + aSupportOnBUp + detectPadding;

			if (!overlapInBFrame) {
				return false;
			}

			// 法向距离过大时也不认为是端面碰撞。
			// 这可以避免一个端面只是处在附近节点区域，却被错误缩小。
			const double radiusA =
				AutoStraightRingAverageRadius(ringA);
			const double radiusB =
				AutoStraightRingAverageRadius(ringB);

			const double normalGapA =
				std::abs(AutoStraightDot3(deltaAB, normalA));

			const double normalGapB =
				std::abs(AutoStraightDot3(deltaAB, normalB));

			const double normalLimit =
				std::max(
					modelDiagonal * 0.018,
					0.45 * std::max(radiusA, radiusB));

			if (normalGapA > normalLimit && normalGapB > normalLimit) {
				return false;
			}

			return true;
		};

		const double detectPadding =
			std::max(modelDiagonal * 0.006, 1e-6);

		for (int iter = 0; iter < 3; ++iter) {
			for (size_t a = 0; a < endpoints.size(); ++a) {
				for (size_t b = a + 1; b < endpoints.size(); ++b) {
					AutoStraightCageBuildItem& itemA =
						items[endpoints[a].itemIndex];

					AutoStraightCageBuildItem& itemB =
						items[endpoints[b].itemIndex];

					if (endpoints[a].itemIndex == endpoints[b].itemIndex) {
						continue;
					}

					if (itemA.sourceChainId == itemB.sourceChainId) {
						continue;
					}

					double centerA[3];
					double centerB[3];
					double ringA[4][3];
					double ringB[4][3];

					AutoStraightGetEndpointCenter(
						itemA,
						endpoints[a].isStart,
						centerA);

					AutoStraightGetEndpointCenter(
						itemB,
						endpoints[b].isStart,
						centerB);

					AutoStraightGetEndpointRing(
						itemA,
						endpoints[a].isStart,
						ringA);

					AutoStraightGetEndpointRing(
						itemB,
						endpoints[b].isStart,
						ringB);

					const double centerGap =
						AutoStraightDistance3(centerA, centerB);

					if (centerGap <= 1e-8) {
						continue;
					}

					const double radiusA =
						AutoStraightRingAverageRadius(ringA);

					const double radiusB =
						AutoStraightRingAverageRadius(ringB);

					const double currentSum =
						radiusA + radiusB;

					const double allowedSum =
						centerGap * 1.08 - detectPadding;

					if (currentSum <= 1e-8 || currentSum <= allowedSum) {
						continue;
					}

					// 新增真实重叠判断：
					// 不满足真实重叠，则不缩放。
					if (!ShouldResolveRealOverlap(
						itemA,
						endpoints[a].isStart,
						itemB,
						endpoints[b].isStart,
						detectPadding))
					{
						continue;
					}

					const double scale =
						std::max(
							0.82,
							std::min(
								0.98,
								allowedSum / currentSum));

					AutoStraightScaleEndpointRing(
						itemA,
						endpoints[a].isStart,
						scale);

					AutoStraightScaleEndpointRing(
						itemB,
						endpoints[b].isStart,
						scale);
				}
			}
		}
	}

	static bool AutoStraightEndpointCanShareFace(
		const AutoStraightCageBuildItem& a,
		bool aStart,
		const AutoStraightCageBuildItem& b,
		bool bStart,
		double modelDiagonal)
	{
		// 关键限制：
		// 只有同一条原始骨架 branch chain 内拆出来的 cage，
		// 才允许共享中间端面。
		if (a.sourceChainId < 0 ||
			b.sourceChainId < 0 ||
			a.sourceChainId != b.sourceChainId)
		{
			return false;
		}

		double ca[3];
		double cb[3];

		AutoStraightGetEndpointCenter(a, aStart, ca);
		AutoStraightGetEndpointCenter(b, bStart, cb);

		const double gap = AutoStraightDistance3(ca, cb);

		const double lengthScale =
			std::max(1e-8, std::min(a.length, b.length));

		const double maxGap =
			std::max(modelDiagonal * 0.018, lengthScale * 0.35);

		if (gap > maxGap) {
			return false;
		}

		const double axisDot = std::abs(AutoStraightDot3(a.axis, b.axis));
		const double minAxisDot =
			std::cos(50.0 * 3.14159265358979323846 / 180.0);

		if (axisDot < minAxisDot) {
			return false;
		}

		return true;
	}

	static void AutoStraightBuildSharedRing(
		const double ringA[4][3],
		const double ringB[4][3],
		double sharedRing[4][3])
	{
		double bestScore = std::numeric_limits<double>::max();
		int bestShift = 0;

		// 只允许 cyclic shift，不做镜像反转。
		for (int shift = 0; shift < 4; ++shift) {
			double score = 0.0;

			for (int i = 0; i < 4; ++i) {
				const int j = (i + shift) % 4;

				const double dx = ringA[i][0] - ringB[j][0];
				const double dy = ringA[i][1] - ringB[j][1];
				const double dz = ringA[i][2] - ringB[j][2];

				score += dx * dx + dy * dy + dz * dz;
			}

			if (score < bestScore) {
				bestScore = score;
				bestShift = shift;
			}
		}

		double alignedB[4][3];

		for (int i = 0; i < 4; ++i) {
			const int j = (i + bestShift) % 4;

			for (int d = 0; d < 3; ++d) {
				alignedB[i][d] = ringB[j][d];
			}
		}

		double centerA[3];
		double centerB[3];
		double sharedCenter[3];

		AutoStraightRingCenter(ringA, centerA);
		AutoStraightRingCenter(alignedB, centerB);

		for (int d = 0; d < 3; ++d) {
			sharedCenter[d] = 0.5 * (centerA[d] + centerB[d]);
		}

		double radiusA = AutoStraightRingAverageRadius(ringA);
		double radiusB = AutoStraightRingAverageRadius(alignedB);

		if (radiusA > radiusB * 1.18) {
			for (int i = 0; i < 4; ++i) {
				for (int d = 0; d < 3; ++d) {
					sharedRing[i][d] =
						sharedCenter[d] + (ringA[i][d] - centerA[d]);
				}
			}
			return;
		}

		if (radiusB > radiusA * 1.18) {
			for (int i = 0; i < 4; ++i) {
				for (int d = 0; d < 3; ++d) {
					sharedRing[i][d] =
						sharedCenter[d] + (alignedB[i][d] - centerB[d]);
				}
			}
			return;
		}

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				sharedRing[i][d] = 0.5 * (ringA[i][d] + alignedB[i][d]);
			}
		}
	}

	static void AutoStraightSnapAdjacentCageFaces(
		std::vector<AutoStraightCageBuildItem>& items,
		double modelDiagonal)
	{
		if (items.size() < 2) {
			return;
		}

		struct EndpointRef
		{
			size_t itemIndex;
			bool isStart;
		};

		std::vector<EndpointRef> endpoints;

		for (size_t i = 0; i < items.size(); ++i) {
			endpoints.push_back({ i, true });
			endpoints.push_back({ i, false });
		}

		std::vector<bool> used(endpoints.size(), false);

		for (size_t a = 0; a < endpoints.size(); ++a) {
			if (used[a]) {
				continue;
			}

			const EndpointRef& ea = endpoints[a];

			size_t bestIndex = static_cast<size_t>(-1);
			double bestGap = std::numeric_limits<double>::max();

			for (size_t b = a + 1; b < endpoints.size(); ++b) {
				if (used[b]) {
					continue;
				}

				const EndpointRef& eb = endpoints[b];

				if (ea.itemIndex == eb.itemIndex) {
					continue;
				}

				if (!AutoStraightEndpointCanShareFace(
					items[ea.itemIndex],
					ea.isStart,
					items[eb.itemIndex],
					eb.isStart,
					modelDiagonal))
				{
					continue;
				}

				double ca[3];
				double cb[3];

				AutoStraightGetEndpointCenter(
					items[ea.itemIndex],
					ea.isStart,
					ca);

				AutoStraightGetEndpointCenter(
					items[eb.itemIndex],
					eb.isStart,
					cb);

				const double gap = AutoStraightDistance3(ca, cb);

				if (gap < bestGap) {
					bestGap = gap;
					bestIndex = b;
				}
			}

			if (bestIndex == static_cast<size_t>(-1)) {
				continue;
			}

			const EndpointRef& eb = endpoints[bestIndex];

			double ringA[4][3];
			double ringB[4][3];
			double sharedRing[4][3];

			AutoStraightGetEndpointRing(
				items[ea.itemIndex],
				ea.isStart,
				ringA);

			AutoStraightGetEndpointRing(
				items[eb.itemIndex],
				eb.isStart,
				ringB);

			AutoStraightBuildSharedRing(
				ringA,
				ringB,
				sharedRing);

			double sharedCenter[3];
			AutoStraightRingCenter(sharedRing, sharedCenter);

			AutoStraightSetEndpointRing(
				items[ea.itemIndex],
				ea.isStart,
				sharedRing);

			AutoStraightSetEndpointRing(
				items[eb.itemIndex],
				eb.isStart,
				sharedRing);

			AutoStraightSetEndpointCenter(
				items[ea.itemIndex],
				ea.isStart,
				sharedCenter);

			AutoStraightSetEndpointCenter(
				items[eb.itemIndex],
				eb.isStart,
				sharedCenter);

			used[a] = true;
			used[bestIndex] = true;

			if (false) qDebug() << "[自动直线控制笼] 共享端面:"
				<< "cage" << static_cast<int>(ea.itemIndex)
				<< (ea.isStart ? "start" : "end")
				<< "<->"
				<< "cage" << static_cast<int>(eb.itemIndex)
				<< (eb.isStart ? "start" : "end")
				<< "gap =" << bestGap;
		}
	}
}

CageInteractorStyle::CageInteractorStyle()
{
	m_ClickState = 0;
	m_refiner = nullptr; // 

	vtkNew<vtkPlaneSource> plane;
	plane->SetOrigin(-0.5, -0.5, 0);
	plane->SetPoint1(0.5, -0.5, 0);
	plane->SetPoint2(-0.5, 0.5, 0);
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(plane->GetOutputPort());
	m_PreviewBoxActor = vtkSmartPointer<vtkActor>::New();
	m_PreviewBoxActor->SetMapper(mapper);
	m_PreviewBoxActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
	m_PreviewBoxActor->GetProperty()->SetRepresentationToWireframe();
	m_PreviewBoxActor->GetProperty()->SetLineWidth(2);
	m_PreviewBoxActor->VisibilityOff();

	vtkNew<vtkPolyDataMapper> fixedMapper;
	fixedMapper->SetInputConnection(plane->GetOutputPort());
	m_FixedBoxActor = vtkSmartPointer<vtkActor>::New();
	m_FixedBoxActor->SetMapper(fixedMapper);
	m_FixedBoxActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
	m_FixedBoxActor->GetProperty()->SetRepresentationToWireframe();
	m_FixedBoxActor->GetProperty()->SetLineWidth(2);
	m_FixedBoxActor->VisibilityOff();

	m_LinesData = vtkSmartPointer<vtkPolyData>::New();
	vtkNew<vtkPoints> pts;
	pts->SetNumberOfPoints(8);
	m_LinesData->SetPoints(pts);

	vtkNew<vtkCellArray> lines;
	for (int i = 0; i < 4; i++) {
		vtkNew<vtkLine> line;
		line->GetPointIds()->SetId(0, i);
		line->GetPointIds()->SetId(1, i + 4);
		lines->InsertNextCell(line);
	}
	m_LinesData->SetLines(lines);

	vtkNew<vtkPolyDataMapper> lineMapper;
	lineMapper->SetInputData(m_LinesData);
	m_LinesActor = vtkSmartPointer<vtkActor>::New();
	m_LinesActor->SetMapper(lineMapper);
	m_LinesActor->GetProperty()->SetColor(0.78, 0.78, 0.74);
	m_LinesActor->GetProperty()->SetLineWidth(2);
	m_LinesActor->VisibilityOff();

	m_StartTransform = vtkSmartPointer<vtkTransform>::New();
}

CageInteractorStyle::~CageInteractorStyle()
{
	if (m_refiner)
	{
		delete m_refiner; // 
	}
}
void CageInteractorStyle::SetReferenceData(vtkSmartPointer<vtkPolyData> data, vtkSmartPointer<vtkCellLocator> locator)
{
	m_RefData = data;
	m_RefLocator = locator;
}

void CageInteractorStyle::SetNonDesignExportSource(
	const QString& stepPath,
	const std::vector<std::array<double, 3>>& elementCenters,
	const std::vector<char>& elementIsNonDesign,
	bool allStepSolidsAreNonDesign)
{
	m_nonDesignExportStepPath = stepPath;
	m_nonDesignExportElementCenters = elementCenters;
	m_nonDesignExportIsNonDesign = elementIsNonDesign;
	m_nonDesignExportAllStepSolids = allStepSolidsAreNonDesign;

	qDebug() << "[导出融合] non-design 源数据已传入 CageInteractorStyle:"
		<< "step =" << m_nonDesignExportStepPath
		<< "elements =" << static_cast<int>(m_nonDesignExportElementCenters.size())
		<< "flags =" << static_cast<int>(m_nonDesignExportIsNonDesign.size())
		<< "all STEP solids =" << m_nonDesignExportAllStepSolids;
}

void CageInteractorStyle::SetNonDesignReferenceSurface(vtkPolyData* nonDesignSurface)
{
	m_nonDesignReferenceSurface = nonDesignSurface;
	qDebug() << "[non-design reference] surface passed to CageInteractorStyle:"
		<< (m_nonDesignReferenceSurface ? m_nonDesignReferenceSurface->GetNumberOfPoints() : 0)
		<< "points"
		<< (m_nonDesignReferenceSurface ? m_nonDesignReferenceSurface->GetNumberOfCells() : 0)
		<< "cells";
}
void CageInteractorStyle::StandardizeCagePoints(double points[4][3], double axis[3], double up[3], double side[3])
{
	double center[3] = { 0, 0, 0 };
	for (int i = 0; i < 4; i++) {
		center[0] += points[i][0] / 4.0;
		center[1] += points[i][1] / 4.0;
		center[2] += points[i][2] / 4.0;
	}

	struct PolarPoint { int id; double angle; };
	std::vector<PolarPoint> sorted;
	for (int i = 0; i < 4; i++) {
		double v[3] = { points[i][0] - center[0], points[i][1] - center[1], points[i][2] - center[2] };
		double x = vtkMath::Dot(v, side);
		double y = vtkMath::Dot(v, up);
		sorted.push_back({ i, std::atan2(y, x) });
	}
	std::sort(sorted.begin(), sorted.end(), [](const PolarPoint& a, const PolarPoint& b) {
		return a.angle < b.angle;
	});

	double temp[4][3];
	for (int i = 0; i < 4; i++) memcpy(temp[i], points[sorted[i].id], sizeof(double) * 3);
	memcpy(points, temp, sizeof(double) * 12);
}
void CageInteractorStyle::OnMouseMove()
{
	if (!m_RefLocator || !m_RefData) return;
	vtkRenderWindowInteractor* interactor = this->GetInteractor();
	if (!interactor) return;

	int x = interactor->GetEventPosition()[0];
	int y = interactor->GetEventPosition()[1];

	vtkRenderer* ren = interactor->FindPokedRenderer(x, y);
	if (!ren) {
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isStraightCageEditMode) {
		HideManualPreviewActors();
		if (m_straightCageEditActive) {
			UpdateStraightCageEditPreview(x, y, ren);
			return;
		}
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isStraightCageRotateMode) {
		HideManualPreviewActors();
		if (m_straightCageRotateActive) {
			UpdateStraightCageRotatePreview(x, y, ren);
			return;
		}
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isStraightSegmentEditMode) {
		HideManualPreviewActors();
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isStraightSegmentSelectMode) {
		HideManualPreviewActors();
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isSkeletonExtendMode) {
		if (m_PreviewBoxActor) m_PreviewBoxActor->VisibilityOff();
		if (m_selectedExtendSkelNode >= 0) {
			UpdateSkeletonExtendPreview(x, y, ren);
		}
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_isSkeletonEditMode || m_isSkeletonDeleteMode) {
		if (m_PreviewBoxActor) m_PreviewBoxActor->VisibilityOff();
		vtkInteractorStyleTrackballCamera::OnMouseMove();
		return;
	}

	if (m_PreviewBoxActor) m_PreviewBoxActor->SetPickable(0);
	if (m_FixedBoxActor) m_FixedBoxActor->SetPickable(0);
	if (m_LinesActor) m_LinesActor->SetPickable(0);
	if (m_HistoryCagesActor) m_HistoryCagesActor->SetPickable(0);
	if (m_GlobalSmoothActor) m_GlobalSmoothActor->SetPickable(0);

	vtkActor* pickedActor = nullptr;
	vtkIdType pickedCellId = -1;

	if (m_CageSolidActor) {
		vtkNew<vtkCellPicker> cagePicker;
		cagePicker->SetTolerance(0.005);
		m_CageSolidActor->SetPickable(1);
		m_CageSolidActor->SetVisibility(true);
		cagePicker->AddPickList(m_CageSolidActor);
		cagePicker->PickFromListOn();

		cagePicker->Pick(x, y, 0, ren);

		if (cagePicker->GetActor() == m_CageSolidActor) {
			pickedActor = cagePicker->GetActor();
			pickedCellId = cagePicker->GetCellId();
		}

		m_CageSolidActor->SetVisibility(false);
		m_CageSolidActor->SetPickable(0);
	}

	if (pickedActor == nullptr) {
		vtkNew<vtkCellPicker> bgPicker;
		bgPicker->SetTolerance(0.001);
		bgPicker->Pick(x, y, 0, ren);
		pickedActor = bgPicker->GetActor();
		pickedCellId = bgPicker->GetCellId();
	}

	// ========================================================
	// 1. 鎷惧彇鍒板凡缁忕敓鎴愮殑瀹炰綋鎺у埗绗?
	// ========================================================
	if (pickedActor != nullptr && pickedActor == m_CageSolidActor) {
		vtkPolyData* poly = vtkPolyDataMapper::SafeDownCast(m_CageSolidActor->GetMapper())->GetInput();
		vtkCell* cell = poly->GetCell(pickedCellId);
		if (cell && cell->GetNumberOfPoints() == 4) {
			double p0[3], p1[3], p2[3], p3[3];
			cell->GetPoints()->GetPoint(0, p0);
			cell->GetPoints()->GetPoint(1, p1);
			cell->GetPoints()->GetPoint(2, p2);
			cell->GetPoints()->GetPoint(3, p3);

			double trueCenter[3];
			for (int i = 0; i < 3; i++) trueCenter[i] = (p0[i] + p1[i] + p2[i] + p3[i]) / 4.0;

			double v10[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
			double v20[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
			double faceNormal[3];
			vtkMath::Cross(v10, v20, faceNormal);
			vtkMath::Normalize(faceNormal);

			int cageIndex = pickedCellId / 6;
			if (cageIndex >= 0 && cageIndex < m_allCages.size()) {
				double volCenter[3] = { 0, 0, 0 };
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 3; j++) {
						volCenter[j] += m_allCages[cageIndex].pStart[i][j] / 8.0;
						volCenter[j] += m_allCages[cageIndex].pEnd[i][j] / 8.0;
					}
				}
				double outDir[3] = { trueCenter[0] - volCenter[0], trueCenter[1] - volCenter[1], trueCenter[2] - volCenter[2] };

				if (vtkMath::Dot(faceNormal, outDir) < 0) {
					for (int i = 0; i < 3; i++) faceNormal[i] = -faceNormal[i];
				}
			}

			double camPos[3], camFocal[3], camDir[3];
			ren->GetActiveCamera()->GetPosition(camPos);
			ren->GetActiveCamera()->GetFocalPoint(camFocal);
			for (int i = 0; i < 3; i++) camDir[i] = camFocal[i] - camPos[i];

			if (vtkMath::Dot(faceNormal, camDir) > 0) {
				m_PreviewBoxActor->VisibilityOff();
				this->GetInteractor()->Render();
				return;
			}

			// 馃毃 鏍稿績淇锛氬潥鍐虫憭寮冮暱鐭竟鍒ゆ柇锛佷弗鏍兼牴鎹潰椤剁偣鐨勬嫇鎵戦『搴忛攣瀹氭湰鍦板潗鏍囩郴
			double axis[3], side[3], normalUp[3];
			for (int i = 0; i < 3; i++) axis[i] = faceNormal[i];

			for (int i = 0; i < 3; i++) side[i] = p1[i] - p0[i];
			vtkMath::Normalize(side);

			vtkMath::Cross(axis, side, normalUp);
			vtkMath::Normalize(normalUp);

			double minX = 1e9, maxX = -1e9, minY = 1e9, maxY = -1e9;
			double pts[4][3] = { {p0[0],p0[1],p0[2]}, {p1[0],p1[1],p1[2]}, {p2[0],p2[1],p2[2]}, {p3[0],p3[1],p3[2]} };

			for (int i = 0; i < 4; i++) {
				double dx = pts[i][0] - trueCenter[0];
				double dy = pts[i][1] - trueCenter[1];
				double dz = pts[i][2] - trueCenter[2];
				double projX = dx * side[0] + dy * side[1] + dz * side[2];
				double projY = dx * normalUp[0] + dy * normalUp[1] + dz * normalUp[2];

				if (projX < minX) minX = projX;
				if (projX > maxX) maxX = projX;
				if (projY < minY) minY = projY;
				if (projY > maxY) maxY = projY;
			}
			UpdatePreviewBox(trueCenter, axis, normalUp, side, maxX - minX, maxY - minY);
			m_PreviewBoxActor->VisibilityOn();
			this->GetInteractor()->Render();
			return;
		}
	}

	// ========================================================
	// 2. 鎷惧彇鍒拌儗鏅弬鑰冩ā鍨?
	// ========================================================
	if (pickedCellId != -1 && pickedActor != m_CageSolidActor)
	{
		vtkCell* cell = m_RefData->GetCell(pickedCellId);
		if (cell->GetNumberOfPoints() < 3) return;

		double p0[3], p1[3], p2[3], p3[3];
		cell->GetPoints()->GetPoint(0, p0);
		cell->GetPoints()->GetPoint(1, p1);
		cell->GetPoints()->GetPoint(2, p2);

		bool isQuad = (cell->GetNumberOfPoints() >= 4);
		if (isQuad) cell->GetPoints()->GetPoint(3, p3);

		double cellCenter[3];
		if (isQuad) {
			for (int i = 0; i < 3; i++) cellCenter[i] = (p0[i] + p1[i] + p2[i] + p3[i]) / 4.0;
		}
		else {
			for (int i = 0; i < 3; i++) cellCenter[i] = (p0[i] + p1[i] + p2[i]) / 3.0;
		}

		double v10[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
		double normalUp[3];
		if (isQuad) {
			double v30[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
			vtkMath::Cross(v10, v30, normalUp);
		}
		else {
			double v20[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
			vtkMath::Cross(v10, v20, normalUp);
		}
		vtkMath::Normalize(normalUp);

		double axis[3] = { 1, 0, 0 };
		double dirU[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
		vtkMath::Normalize(dirU);

		static double lastAxis[3] = { 0, 0, 0 };
		if (isQuad) {
			double dirV[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
			vtkMath::Normalize(dirV);
			if (vtkMath::Norm(lastAxis) < 0.1) {
				if (vtkMath::Distance2BetweenPoints(p0, p1) > vtkMath::Distance2BetweenPoints(p0, p3)) {
					for (int i = 0; i < 3; i++) lastAxis[i] = dirU[i];
				}
				else {
					for (int i = 0; i < 3; i++) lastAxis[i] = dirV[i];
				}
			}
			double dotU = std::abs(vtkMath::Dot(dirU, lastAxis));
			double dotV = std::abs(vtkMath::Dot(dirV, lastAxis));
			if (dotU > dotV) {
				double sign = (vtkMath::Dot(dirU, lastAxis) > 0) ? 1.0 : -1.0;
				for (int i = 0; i < 3; i++) axis[i] = dirU[i] * sign;
			}
			else {
				double sign = (vtkMath::Dot(dirV, lastAxis) > 0) ? 1.0 : -1.0;
				for (int i = 0; i < 3; i++) axis[i] = dirV[i] * sign;
			}
		}
		else {
			if (vtkMath::Norm(lastAxis) < 0.1) {
				for (int i = 0; i < 3; i++) lastAxis[i] = dirU[i];
			}
			double dotL = vtkMath::Dot(lastAxis, normalUp);
			for (int i = 0; i < 3; i++) axis[i] = lastAxis[i] - dotL * normalUp[i];
			if (vtkMath::Norm(axis) < 1e-5) {
				for (int i = 0; i < 3; i++) axis[i] = dirU[i];
			}
			vtkMath::Normalize(axis);
		}
		for (int i = 0; i < 3; i++) lastAxis[i] = axis[i];

		double side[3];
		if (m_ClickState == 1 && m_PreviewBoxActor->GetVisibility()) {
			vtkTransform* prevTransform = vtkTransform::SafeDownCast(m_PreviewBoxActor->GetUserTransform());
			if (prevTransform) {
				double prevX[3] = { prevTransform->GetMatrix()->GetElement(0,0),
									prevTransform->GetMatrix()->GetElement(1,0),
									prevTransform->GetMatrix()->GetElement(2,0) };
				double prevZ[3] = { prevTransform->GetMatrix()->GetElement(0,2),
									prevTransform->GetMatrix()->GetElement(1,2),
									prevTransform->GetMatrix()->GetElement(2,2) };

				if (vtkMath::Dot(axis, prevX) < 0) {
					for (int i = 0; i < 3; i++) axis[i] = -axis[i];
				}
				if (vtkMath::Dot(normalUp, prevZ) < 0) {
					for (int i = 0; i < 3; i++) normalUp[i] = -normalUp[i];
				}
			}
		}

		double t, pcoords[3];
		int subId;
		vtkIdType oppositeCellId;

		double baseSize = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		if (!isQuad) {
			double dist12 = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
			double dist20 = sqrt(vtkMath::Distance2BetweenPoints(p2, p0));
			if (dist12 > baseSize) baseSize = dist12;
			if (dist20 > baseSize) baseSize = dist20;
		}

		double rayStartDown[3], rayDownEnd[3];
		for (int i = 0; i < 3; i++) {
			rayStartDown[i] = cellCenter[i] - normalUp[i] * 0.1;
			rayDownEnd[i] = cellCenter[i] - normalUp[i] * 500.0;
		}

		double pDown[3];
		int hitDown = m_RefLocator->IntersectWithLine(rayStartDown, rayDownEnd, 0.001, t, pDown, pcoords, subId, oppositeCellId);

		double trueCenter[3];
		double depth = baseSize;

		if (hitDown) {
			depth = sqrt(vtkMath::Distance2BetweenPoints(cellCenter, pDown));
			for (int i = 0; i < 3; i++) trueCenter[i] = (cellCenter[i] + pDown[i]) / 2.0;
		}
		else {
			for (int i = 0; i < 3; i++) trueCenter[i] = cellCenter[i] - normalUp[i] * (depth / 2.0);
		}

		// 杩欓噷蹇呴』瑕佸厛鐢ㄦ櫘閫氳绠楄繃涓€閬嶏紝鑾峰彇鍒板ぇ姒傜殑鍒濆 side锛屾柟渚挎眰瀹藉害
		vtkMath::Cross(normalUp, axis, side);
		vtkMath::Normalize(side);
		vtkMath::Cross(side, normalUp, axis);
		vtkMath::Normalize(axis);

		double rayRightEnd[3], rayLeftEnd[3];
		for (int i = 0; i < 3; i++) {
			rayRightEnd[i] = trueCenter[i] + side[i] * 500.0;
			rayLeftEnd[i] = trueCenter[i] - side[i] * 500.0;
		}

		double pRight[3], pLeft[3];
		int hitRight = m_RefLocator->IntersectWithLine(trueCenter, rayRightEnd, 0.001, t, pRight, pcoords, subId, oppositeCellId);
		int hitLeft = m_RefLocator->IntersectWithLine(trueCenter, rayLeftEnd, 0.001, t, pLeft, pcoords, subId, oppositeCellId);

		double width = baseSize;
		if (hitRight && hitLeft) width = sqrt(vtkMath::Distance2BetweenPoints(pRight, pLeft));
		else if (hitRight) width = 2.0 * sqrt(vtkMath::Distance2BetweenPoints(trueCenter, pRight));
		else if (hitLeft) width = 2.0 * sqrt(vtkMath::Distance2BetweenPoints(trueCenter, pLeft));

		if (width < 0.1 * baseSize) width = baseSize;
		if (depth < 0.1 * baseSize) depth = baseSize;


		// =======================================================
		// 馃毃 鏍稿績淇锛氬紩鍏?RMF 寮鸿娑堥櫎寤朵几鏃剁殑鎷撴墤鎵洸
		// =======================================================
		if (m_ClickState == 1 && m_StartTransform)
		{
			double x0[3], t0[3], u0[3];
			vtkMatrix4x4* mat = m_StartTransform->GetMatrix();
			for (int i = 0; i < 3; i++) {
				u0[i] = mat->GetElement(i, 1); // Previous Up (Y)
				t0[i] = mat->GetElement(i, 2); // Previous Axis (Z)
				x0[i] = mat->GetElement(i, 3); // Previous Center
			}

			// 鏂扮殑鍒囩嚎鏂瑰悜鎸囧悜褰撳墠榧犳爣鍦ㄨ儗鏅笂姹傚緱鐨勪綅缃?
			double t1[3] = { trueCenter[0] - x0[0], trueCenter[1] - x0[1], trueCenter[2] - x0[2] };
			if (vtkMath::Norm(t1) > 1e-6) vtkMath::Normalize(t1);
			else for (int i = 0; i < 3; i++) t1[i] = t0[i];

			// RMF 璁＄畻鍑哄畬缇庢棤鎵洸鐨?normalUp (U) 鍜?side (V)
			ComputeRMF(x0, t0, u0, trueCenter, t1, normalUp, side);
			for (int i = 0; i < 3; i++) axis[i] = t1[i];
		}


		bool isPoison = false;
		for (int i = 0; i < 3; i++) {
			if (std::isnan(trueCenter[i]) || std::isnan(axis[i]) || std::isnan(normalUp[i]) || std::isnan(side[i])) isPoison = true;
		}
		if (std::isnan(width) || std::isnan(depth)) isPoison = true;

		if (!isPoison) {
			UpdatePreviewBox(trueCenter, axis, normalUp, side, width * 1.1, depth * 1.1);
			m_PreviewBoxActor->VisibilityOn();
		}
		else {
			m_PreviewBoxActor->VisibilityOff();
		}
		this->GetInteractor()->Render();
	}
	else {
		m_PreviewBoxActor->VisibilityOff();
		this->GetInteractor()->Render();
	}

	vtkInteractorStyleTrackballCamera::OnMouseMove();
}
void CageInteractorStyle::UpdatePreviewBox(double center[3], double axis[3], double up[3], double side[3], double width, double depth)
{
	vtkNew<vtkMatrix4x4> matrix;
	matrix->Identity();
	for (int i = 0; i < 3; i++) {
		matrix->SetElement(i, 0, side[i]);
		matrix->SetElement(i, 1, up[i]);
		matrix->SetElement(i, 2, axis[i]);
		matrix->SetElement(i, 3, center[i]);
	}
	vtkNew<vtkTransform> transform;
	transform->PostMultiply();
	transform->Scale(width, depth, 1.0);
	transform->Concatenate(matrix);
	m_PreviewBoxActor->SetUserTransform(transform);

	if (m_ClickState == 1) {
		double endCorners[4][3], startCorners[4][3];
		GetTransformCorners(transform, endCorners);
		for (int i = 0; i < 4; i++) m_LinesData->GetPoints()->GetPoint(i, startCorners[i]);

		// 馃毃 缁堟瀬鎶楁壄鏇诧細2D 鏃堕挓鎶曞奖瀵归綈
		auto MatchTwistFree = [](double A[4][3], double B[4][3], int mapA2B[4]) {
			double cA[3] = { 0 }, cB[3] = { 0 };
			for (int i = 0; i < 4; i++) for (int d = 0; d < 3; d++) { cA[d] += A[i][d] / 4.0; cB[d] += B[i][d] / 4.0; }

			double V[3] = { cB[0] - cA[0], cB[1] - cA[1], cB[2] - cA[2] };
			if (vtkMath::Norm(V) < 1e-6) {
				double v1[3] = { A[1][0] - A[0][0], A[1][1] - A[0][1], A[1][2] - A[0][2] };
				double v2[3] = { A[2][0] - A[0][0], A[2][1] - A[0][1], A[2][2] - A[0][2] };
				vtkMath::Cross(v1, v2, V);
			}
			vtkMath::Normalize(V);

			double Y[3] = { 0,0,1 };
			if (std::abs(V[2]) > 0.99) { Y[0] = 1; Y[1] = 0; Y[2] = 0; }
			double X[3];
			vtkMath::Cross(Y, V, X); vtkMath::Normalize(X);
			vtkMath::Cross(V, X, Y); vtkMath::Normalize(Y);

			struct PolarNode { int id; double angle; };
			std::vector<PolarNode> sortA(4), sortB(4);
			for (int i = 0; i < 4; i++) {
				double vA[3] = { A[i][0] - cA[0], A[i][1] - cA[1], A[i][2] - cA[2] };
				sortA[i] = { i, std::atan2(vtkMath::Dot(vA, Y), vtkMath::Dot(vA, X)) };
				double vB[3] = { B[i][0] - cB[0], B[i][1] - cB[1], B[i][2] - cB[2] };
				sortB[i] = { i, std::atan2(vtkMath::Dot(vB, Y), vtkMath::Dot(vB, X)) };
			}
			std::sort(sortA.begin(), sortA.end(), [](const PolarNode& a, const PolarNode& b) { return a.angle < b.angle; });
			std::sort(sortB.begin(), sortB.end(), [](const PolarNode& a, const PolarNode& b) { return a.angle < b.angle; });

			for (int i = 0; i < 4; i++) mapA2B[sortA[i].id] = sortB[i].id;
		};

		int mapA2B[4];
		MatchTwistFree(startCorners, endCorners, mapA2B);

		for (int i = 0; i < 4; i++) {
			int shifted_idx = (i + m_ManualRotationOffset) % 4; // 鐢ㄦ埛婊氳疆渚濈劧鍙互鎺у埗寰皟
			m_LinesData->GetPoints()->SetPoint(i + 4, endCorners[mapA2B[shifted_idx]]);
		}
		m_LinesData->GetPoints()->Modified();
	}
}

namespace
{
	static void CopyRing4(double dst[4][3], const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	}


	static void CopyRing4MaybeReverse(double dst[4][3], const double src[4][3], bool reverse)
	{
		static const int normalOrder[4] = { 0, 1, 2, 3 };
		static const int reverseOrder[4] = { 0, 3, 2, 1 };

		const int* order = reverse ? reverseOrder : normalOrder;

		for (int i = 0; i < 4; ++i) {
			int id = order[i];
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[id][d];
			}
		}
	}

	// m_CageSolidActor 閲屾瘡涓帶鍒剁鏈?6 涓崟鍏冿細
	// 0 = 璧峰绔洊
	// 1 = 缁撴潫绔洊
	// 2~5 = side faces
	static bool DecodePickedCageCap(
		vtkIdType pickedCellId,
		int cageCount,
		int& cageIndex,
		bool& isStartCap)
	{
		if (pickedCellId < 0) return false;

		cageIndex = static_cast<int>(pickedCellId / 6);
		int faceOffset = static_cast<int>(pickedCellId % 6);

		if (cageIndex < 0 || cageIndex >= cageCount) return false;

		if (faceOffset == 0) {
			isStartCap = true;
			return true;
		}

		if (faceOffset == 1) {
			isStartCap = false;
			return true;
		}

		return false;
	}

	static bool DecodePickedCageFace(
		vtkIdType pickedCellId,
		int cageCount,
		int& cageIndex,
		int& faceOffset)
	{
		if (pickedCellId < 0) return false;

		cageIndex = static_cast<int>(pickedCellId / 6);
		faceOffset = static_cast<int>(pickedCellId % 6);

		if (cageIndex < 0 || cageIndex >= cageCount) return false;
		if (faceOffset < 0 || faceOffset > 5) return false;

		return true;
	}

	static bool GetCageFaceLocalIds(int faceOffset, int ids[4])
	{
		static const int faceTable[6][4] = {
			{ 0, 1, 2, 3 },
			{ 4, 5, 6, 7 },
			{ 0, 1, 5, 4 },
			{ 1, 2, 6, 5 },
			{ 2, 3, 7, 6 },
			{ 3, 0, 4, 7 }
		};

		if (faceOffset < 0 || faceOffset > 5) return false;
		for (int i = 0; i < 4; ++i) {
			ids[i] = faceTable[faceOffset][i];
		}

		return true;
	}

	static bool GetCanonicalCageFaceRing(
		const std::vector<CageSegment>& cages,
		int cageIndex,
		int faceOffset,
		double outRing[4][3])
	{
		if (cageIndex < 0 || cageIndex >= static_cast<int>(cages.size())) {
			return false;
		}

		int ids[4];
		if (!GetCageFaceLocalIds(faceOffset, ids)) {
			return false;
		}

		const CageSegment& cage = cages[cageIndex];
		for (int i = 0; i < 4; ++i) {
			const int id = ids[i];
			const double* src = id < 4 ? cage.pStart[id] : cage.pEnd[id - 4];
			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = src[d];
			}
		}

		return true;
	}

	static vtkSmartPointer<vtkPolyData> BuildRingQuadData(const double ring[4][3])
	{
		vtkNew<vtkPoints> points;
		for (int i = 0; i < 4; ++i) {
			points->InsertNextPoint(ring[i]);
		}

		vtkNew<vtkCellArray> polys;
		polys->InsertNextCell(4);
		for (vtkIdType i = 0; i < 4; ++i) {
			polys->InsertCellPoint(i);
		}

		vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
		data->SetPoints(points);
		data->SetPolys(polys);
		return data;
	}

	// 涓€瀹氳浠?m_allCages 鍙栫闈㈢殑鍘熷鎷撴墤椤哄簭銆?
	// 涓嶈浠?vtkCell 鍙栵紝鍥犱负 solid actor 鐨勮捣濮嬬鐩?鏄弽鐫€鍐欑殑銆?
	static bool GetCanonicalCageCapRing(
		const std::vector<CageSegment>& cages,
		int cageIndex,
		bool isStartCap,
		double outRing[4][3])
	{
		if (cageIndex < 0 || cageIndex >= static_cast<int>(cages.size())) {
			return false;
		}

		const CageSegment& cage = cages[cageIndex];

		for (int i = 0; i < 4; ++i) {
			const double* src = isStartCap ? cage.pStart[i] : cage.pEnd[i];

			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = src[d];
			}
		}

		return true;
	}

	// 宸叉湁 cage 鐨勭闈紝鍦ㄥ畠鑷繁鐨勪晶闈㈤噷鍛堢幇鐨勮竟鐣屾柟鍚戯細
	//
	// existing pStart = +1
	// existing pEnd   = -1
	static int ExistingCapBoundarySign(bool existingCapIsStart)
	{
		return existingCapIsStart ? +1 : -1;
	}

	// 鏂?cage 鐨勭鐐癸紝鍦ㄥ畠鑷繁鐨勪晶闈㈤噷澶╃劧鍛堢幇鐨勮竟鐣屾柟鍚戯細
	//
	// new pStart = +1
	// new pEnd   = -1
	static int NewEndpointBoundarySign(bool newEndpointIsStart)
	{
		return newEndpointIsStart ? +1 : -1;
	}

	// 鍒ゆ柇锛氫负浜嗚鏂?cage 鍜屽凡鏂?cage 鍦ㄥ叕鍏辫竟涓婃柟鍚戠浉鍙嶏紝
	// 鏂?cage 鐨勮繖涓闈㈡槸鍚﹂渶瑕佸弽杞洓鐐圭幆銆?
	static bool NeedReverseForConnectedEndpoint(
		bool existingCapIsStart,
		bool newEndpointIsStart)
	{
		int existingSign = ExistingCapBoundarySign(existingCapIsStart);
		int desiredNewSign = -existingSign;
		int naturalNewSign = NewEndpointBoundarySign(newEndpointIsStart);

		return naturalNewSign != desiredNewSign;
	}

	// 鍙厑璁稿惊鐜钩绉伙紝涓嶅厑璁搁暅鍍忋€?
	//
	// 杩欐牱鍙互鍑忓皯鎵浆锛屼絾涓嶄細鏀瑰彉鐜悜銆?
	static void AlignRingByCyclicShift(
		const double fixedRing[4][3],
		const double movingRing[4][3],
		double outMovingRing[4][3])
	{
		double bestScore = 1e100;
		int bestShift = 0;

		for (int shift = 0; shift < 4; ++shift) {
			double score = 0.0;

			for (int i = 0; i < 4; ++i) {
				int j = (i + shift) % 4;
				score += vtkMath::Distance2BetweenPoints(fixedRing[i], movingRing[j]);
			}

			if (score < bestScore) {
				bestScore = score;
				bestShift = shift;
			}
		}

		for (int i = 0; i < 4; ++i) {
			int j = (i + bestShift) % 4;
			for (int d = 0; d < 3; ++d) {
				outMovingRing[i][d] = movingRing[j][d];
			}
		}
	}


	// 手动连接/自动连接统一使用的 8 候选 ring 对齐：
	// 4 种循环平移 + 4 种反向循环平移，选择连接长度最短且四边面最不扭的对应。
	// 这比旧的 AlignRingByCyclicShift 更安全，因为旧逻辑禁止反向，遇到端面环方向不一致时只会打印日志继续扭着连。
	struct RingBridgeAlignScore
	{
		double totalLen2 = 0.0;
		double maxLen2 = 0.0;
		double diagDiff = 0.0;
		double edgeParallelPenalty = 0.0;
	};

	static double RingBridgeDotNormalizedLocal(
		const double a0[3], const double a1[3],
		const double b0[3], const double b1[3])
	{
		double ea[3] = { a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2] };
		double eb[3] = { b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2] };
		double la = vtkMath::Norm(ea);
		double lb = vtkMath::Norm(eb);
		if (la < 1e-10 || lb < 1e-10) return 0.0;
		return vtkMath::Dot(ea, eb) / (la * lb);
	}

	static RingBridgeAlignScore EvaluateRingBridgeAlign(
		const double a[4][3],
		const double b[4][3])
	{
		RingBridgeAlignScore score;
		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) & 3;
			double l2 = vtkMath::Distance2BetweenPoints(a[i], b[i]);
			score.totalLen2 += l2;
			score.maxLen2 = std::max(score.maxLen2, l2);

			double diagA = std::sqrt(vtkMath::Distance2BetweenPoints(a[i], b[j]));
			double diagB = std::sqrt(vtkMath::Distance2BetweenPoints(a[j], b[i]));
			score.diagDiff = std::max(score.diagDiff, std::abs(diagA - diagB));

			double dot = RingBridgeDotNormalizedLocal(a[i], a[j], b[i], b[j]);
			score.edgeParallelPenalty += (1.0 - dot);
		}
		return score;
	}

	static bool IsBetterRingBridgeAlign(
		const RingBridgeAlignScore& a,
		const RingBridgeAlignScore& b)
	{
		// 先保证对应边最短，避免手动重连时远角对远角造成麻花；
		// 再用对角差和边平行度消除局部交叉候选。
		if (a.totalLen2 < b.totalLen2 - 1e-8) return true;
		if (a.totalLen2 > b.totalLen2 + 1e-8) return false;
		if (a.maxLen2 < b.maxLen2 - 1e-8) return true;
		if (a.maxLen2 > b.maxLen2 + 1e-8) return false;
		if (a.diagDiff < b.diagDiff - 1e-7) return true;
		if (a.diagDiff > b.diagDiff + 1e-7) return false;
		return a.edgeParallelPenalty < b.edgeParallelPenalty;
	}

	static void BuildRingBridgeCandidate(
		const double src[4][3],
		bool reverse,
		int shift,
		double out[4][3])
	{
		static const int forwardOrder[4] = { 0, 1, 2, 3 };
		static const int reverseOrder[4] = { 0, 3, 2, 1 };
		const int* order = reverse ? reverseOrder : forwardOrder;
		for (int i = 0; i < 4; ++i) {
			int id = order[(i + shift) & 3];
			for (int d = 0; d < 3; ++d) {
				out[i][d] = src[id][d];
			}
		}
	}

	static void AlignRingByEightCandidates(
		const double fixedRing[4][3],
		const double movingRing[4][3],
		double outMovingRing[4][3],
		const char* debugTag)
	{
		double raw[4][3];
		CopyRing4(raw, movingRing);
		RingBridgeAlignScore rawScore = EvaluateRingBridgeAlign(fixedRing, raw);

		double best[4][3];
		CopyRing4(best, raw);
		RingBridgeAlignScore bestScore = rawScore;
		bool bestReverse = false;
		int bestShift = 0;

		for (int family = 0; family < 2; ++family) {
			bool reverse = (family == 1);
			for (int shift = 0; shift < 4; ++shift) {
				double cand[4][3];
				BuildRingBridgeCandidate(movingRing, reverse, shift, cand);
				RingBridgeAlignScore sc = EvaluateRingBridgeAlign(fixedRing, cand);
				if (IsBetterRingBridgeAlign(sc, bestScore)) {
					bestScore = sc;
					CopyRing4(best, cand);
					bestReverse = reverse;
					bestShift = shift;
				}
			}
		}

		CopyRing4(outMovingRing, best);
		const bool changed = bestReverse || bestShift != 0 ||
			(bestScore.totalLen2 < rawScore.totalLen2 - 1e-8) ||
			(bestScore.diagDiff < rawScore.diagDiff - 1e-7);

		if (changed || bestScore.diagDiff > rawScore.diagDiff + 1e-6) {
			if (false) qDebug() << (debugTag ? debugTag : "[RingBridgeAlign8]")
				<< "changed =" << changed
				<< "reverse =" << bestReverse
				<< "shift =" << bestShift
				<< "rawTotalLen2 =" << rawScore.totalLen2
				<< "bestTotalLen2 =" << bestScore.totalLen2
				<< "rawMaxLen2 =" << rawScore.maxLen2
				<< "bestMaxLen2 =" << bestScore.maxLen2
				<< "rawDiagDiff =" << rawScore.diagDiff
				<< "bestDiagDiff =" << bestScore.diagDiff
				<< "rawParallelPenalty =" << rawScore.edgeParallelPenalty
				<< "bestParallelPenalty =" << bestScore.edgeParallelPenalty;
		}
	}
}

void CageInteractorStyle::OnLeftButtonDown()
{
	vtkRenderWindowInteractor* interactor = this->GetInteractor();
	if (!interactor) return;

	int x = interactor->GetEventPosition()[0];
	int y = interactor->GetEventPosition()[1];

	vtkRenderer* ren = interactor->FindPokedRenderer(x, y);
	if (!ren) {
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	if (m_isStraightCageEditMode) {
		HideManualPreviewActors();
		if (HandleStraightCageEditLeftClick(x, y, ren)) {
			return;
		}
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	if (m_isStraightCageRotateMode) {
		HideManualPreviewActors();
		if (HandleStraightCageRotateLeftClick(x, y, ren)) {
			return;
		}
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	if (m_isStraightSegmentEditMode) {
		HideManualPreviewActors();
		if (HandleStraightSegmentTrimClick(x, y, ren)) {
			return;
		}

		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	if (m_isStraightSegmentSelectMode) {
		HideManualPreviewActors();
		if (HandleStraightSegmentSelectClick(x, y, ren)) {
			return;
		}

		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	// ============================================================
	// 1A. 骨架边删除模式：点击最近的红色骨架边，删除该边。
	// ============================================================
	if (m_isSkeletonDeleteMode) {
		const int edgeIndex = PickSkeletonEdge(x, y, ren);
		if (edgeIndex >= 0) {
			DeleteSkeletonEdge(edgeIndex);
			UpdateSkeletonHighlight(ren);
			qDebug() << "[骨架编辑] 已删除骨架边" << edgeIndex;
			return;
		}
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	// ============================================================
	// 1B. 骨架延伸/插入模式：
	//     第一次左键选择任意骨架节点；移动鼠标显示预览线；第二次左键确认新增分支。
	//     如果第二次点击命中已有骨架节点/边，则直接接入已有骨架。
	// ============================================================
	if (m_isSkeletonExtendMode) {
		if (m_selectedExtendSkelNode < 0) {
			const int pickedNode = PickSkeletonNode(x, y, ren);
			if (pickedNode >= 0) {
				m_selectedExtendSkelNode = pickedNode;
				m_selectedSkelNode1 = pickedNode;
				UpdateSkeletonHighlight(ren);
				UpdateSkeletonExtendPreview(x, y, ren);
				qDebug() << "[骨架延伸] 已选择骨架节点" << pickedNode
					<< "，当前节点度数 =" << GetSkeletonNodeDegree(pickedNode)
					<< "。移动鼠标预览，再次左键确认；右键取消。";
				return;
			}
			vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
			return;
		}
		else {
			if (ExtendSkeletonFromSelectedTail(x, y, ren)) {
				qDebug() << "[骨架延伸] 延伸/插入完成，已解除当前节点锁定，可重新选择其它节点。";
			}
			ClearSkeletonExtendSelection(ren);
		}
		return;
	}

	// ============================================================
	// 1. 骨架节点编辑模式：点击两个骨架节点，合并两点之间的骨架路径。
	// ============================================================
	if (m_isSkeletonEditMode) {
		int pickedNode = PickSkeletonNode(x, y, ren);

		if (pickedNode != -1) {
			if (m_selectedSkelNode1 == -1) {
				m_selectedSkelNode1 = pickedNode;
				UpdateSkeletonHighlight(ren);
				qDebug() << "[骨架编辑] 已选择第一个骨架点";
			}
			else if (pickedNode != m_selectedSkelNode1) {
				MergeSkeletonNodes(m_selectedSkelNode1, pickedNode);
				m_selectedSkelNode1 = -1;
				UpdateSkeletonHighlight(ren);
				qDebug() << "[骨架编辑] 骨架点合并完成";
			}
			return;
		}

		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	// ============================================================
	// 2. 浼樺厛鎷惧彇宸叉湁鎺у埗绗肩闈?
	// ============================================================
	vtkActor* pickedActor = nullptr;
	vtkIdType pickedCellId = -1;

	if (m_CageSolidActor) {
		vtkNew<vtkCellPicker> cagePicker;
		cagePicker->SetTolerance(0.005);

		m_CageSolidActor->SetPickable(1);
		m_CageSolidActor->SetVisibility(true);

		cagePicker->AddPickList(m_CageSolidActor);
		cagePicker->PickFromListOn();

		cagePicker->Pick(x, y, 0, ren);

		if (cagePicker->GetActor() == m_CageSolidActor) {
			pickedActor = cagePicker->GetActor();
			pickedCellId = cagePicker->GetCellId();
		}

		m_CageSolidActor->SetVisibility(false);
		m_CageSolidActor->SetPickable(0);
	}

	// 濡傛灉娌℃湁鎷惧彇鍒板凡鏈夋帶鍒剁锛屽啀鎷惧彇鑳屾櫙妯″瀷
	if (pickedActor == nullptr) {
		vtkNew<vtkCellPicker> bgPicker;
		bgPicker->SetTolerance(0.001);
		bgPicker->Pick(x, y, 0, ren);

		pickedActor = bgPicker->GetActor();
		pickedCellId = bgPicker->GetCellId();
	}

	if (pickedCellId == -1) {
		if (m_PreviewBoxActor) {
			m_PreviewBoxActor->VisibilityOff();
		}

		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
		return;
	}

	// ============================================================
	// 3. 绗竴娆＄偣鍑伙細纭畾鏂版帶鍒剁鐨勮捣鐐?
	// ============================================================
	if (m_ClickState == 0)
	{
		m_IsSnapped = false;
		m_SnappedIsStart = false;
		m_SnappedCageIndex = -1;
		m_SnappedFaceOffset = -1;

		int snappedCageIndex = -1;
		bool snappedIsStart = false;
		int snappedFaceOffset = -1;
		double snappedRing[4][3];

		// ------------------------------------------------------------
		// 濡傛灉绗竴娆＄偣涓凡鏈?cage 鐨?璧峰/缁撴潫绔洊銆?
		// 蹇呴』浠?m_allCages 鍙栧師濮嬫嫇鎵戦『搴忋€?
		// ------------------------------------------------------------
		if (pickedActor == m_CageSolidActor) {
			if (DecodePickedCageFace(
				pickedCellId,
				static_cast<int>(m_allCages.size()),
				snappedCageIndex,
				snappedFaceOffset))
			{
				const bool snappedIsCap = (snappedFaceOffset == 0 || snappedFaceOffset == 1);
				snappedIsStart = (snappedFaceOffset == 0);

				if ((!m_isBranchMode || snappedIsCap) &&
					GetCanonicalCageFaceRing(
					m_allCages,
					snappedCageIndex,
					snappedFaceOffset,
					snappedRing))
				{
					m_IsSnapped = true;
					m_SnappedIsStart = snappedIsStart;
					m_SnappedCageIndex = snappedCageIndex;
					m_SnappedFaceOffset = snappedFaceOffset;

					CopyRing4(m_SnappedPoints, snappedRing);
				}
			}
		}

		// ------------------------------------------------------------
		// 鍒嗘敮閫夋嫨妯″紡淇濇寔鍘熼€昏緫锛屽彧鏄?snapped 鏁版嵁鏉ユ簮鍙樻纭€?
		// ------------------------------------------------------------
		if (m_isBranchMode && m_IsSnapped) {
			bool alreadySelected = false;

			for (auto& face : m_selectedFacesForBranch) {
				if (face.cageIndex == snappedCageIndex &&
					face.isStart == snappedIsStart)
				{
					alreadySelected = true;
					break;
				}
			}

			if (alreadySelected) return;

			m_selectedFacesForBranch.push_back({ snappedCageIndex, snappedIsStart });

			int currentN = static_cast<int>(m_selectedFacesForBranch.size());
			if (false) qDebug() << "[分支模式] 当前已选择端面数量 =" << currentN;

			UpdateBranchHighlight(ren);

			if (currentN >= 2) {
				ShowSkeletonPath(ren);
			}

			if (currentN < 3) {
				if (false) qDebug() << "[分支模式] 端面数量还不足，当前数量 =" << currentN;
			}
			else {
				if (false) qDebug() << "[分支模式] 已满足复杂分支生成条件，端面数量 =" << currentN;
			}

			this->GetInteractor()->GetRenderWindow()->Render();
			return;
		}

		// ------------------------------------------------------------
		// 璁板綍璧风偣 transform锛岀敤浜庢樉绀哄浐瀹氭銆?
		// ------------------------------------------------------------
		m_StartTransform->DeepCopy(m_PreviewBoxActor->GetUserTransform());

		double startCorners[4][3];

		if (m_IsSnapped) {
			CopyRing4(startCorners, m_SnappedPoints);
		}
		else {
			GetTransformCorners(m_StartTransform, startCorners);
		}

		if (m_FixedBoxActor) {
			m_FixedBoxActor->SetUserTransform(nullptr);
			vtkPolyDataMapper::SafeDownCast(m_FixedBoxActor->GetMapper())->SetInputData(
				BuildRingQuadData(startCorners));
			m_FixedBoxActor->VisibilityOn();
		}

		for (int i = 0; i < 4; ++i) {
			m_LinesData->GetPoints()->SetPoint(i, startCorners[i]);
		}

		m_LinesActor->VisibilityOn();
		m_ClickState = 1;
	}

	// ============================================================
	// 4. 绗簩娆＄偣鍑伙細纭畾鏂版帶鍒剁鐨勭粓鐐癸紝骞剁湡姝ｇ敓鎴?CageSegment
	// ============================================================
	else if (m_ClickState == 1)
	{
		// ------------------------------------------------------------
		// 绗簩娆＄偣鍑讳篃瑕佹娴嬫槸鍚︾偣涓簡宸叉湁 cage 鐨勭闈€?
		// ------------------------------------------------------------
		bool endIsSnapped = false;
		bool endSnappedIsStart = false;
		int endSnappedCageIndex = -1;
		int endSnappedFaceOffset = -1;
		double endSnappedRing[4][3];

		if (pickedActor == m_CageSolidActor) {
			if (DecodePickedCageFace(
				pickedCellId,
				static_cast<int>(m_allCages.size()),
				endSnappedCageIndex,
				endSnappedFaceOffset))
			{
				endSnappedIsStart = (endSnappedFaceOffset == 0);
				if (GetCanonicalCageFaceRing(
					m_allCages,
					endSnappedCageIndex,
					endSnappedFaceOffset,
					endSnappedRing))
				{
					endIsSnapped = true;
				}
			}
		}

		// ------------------------------------------------------------
		// 璇诲彇褰撳墠棰勮 涓殑璧风偣鍜岀粓鐐广€?
		// ------------------------------------------------------------
		double rawStart[4][3];
		double rawEnd[4][3];

		for (int i = 0; i < 4; ++i) {
			m_LinesData->GetPoints()->GetPoint(i, rawStart[i]);
			m_LinesData->GetPoints()->GetPoint(i + 4, rawEnd[i]);
		}

		if (endIsSnapped) {
			CopyRing4(rawEnd, endSnappedRing);
		}

		// ------------------------------------------------------------
		// 鍏抽敭淇锛?
		//
		// 涓婁竴鐗堝彧鍙嶈浆浜嗚繛鎺ョ锛屾墍浠ュ嚭鐜伴夯鑺便€?
		//
		// 杩欎竴鐗堟敼鎴愶細
		// 濡傛灉鍙湁涓€绔?snap锛屽苟涓旇繖绔负浜嗘帴鍙ｅ钩婊戦渶瑕佸弽杞紝
		// 閭ｄ箞鍙︿竴绔篃涓€璧峰弽杞€?
		//
		// 杩欐牱锛?
		// 1. 鎺ュ彛澶勫叡浜竟鏂瑰悜姝ｇ‘銆?
		// 2. 鎺у埗绗兼湰浣撲笉浼氭壄鎴?X銆?
		// ------------------------------------------------------------
		bool startNeedsReverse = false;
		bool endNeedsReverse = false;

		if (m_IsSnapped && (m_SnappedFaceOffset == 0 || m_SnappedFaceOffset == 1)) {
			startNeedsReverse = NeedReverseForConnectedEndpoint(
				m_SnappedIsStart,
				true   // 鏂?cage 璧风偣鏄?pStart
			);
		}

		if (endIsSnapped && (endSnappedFaceOffset == 0 || endSnappedFaceOffset == 1)) {
			endNeedsReverse = NeedReverseForConnectedEndpoint(
				endSnappedIsStart,
				false  // 鏂?cage 缁堢偣銆俻End
			);
		}

		bool reverseStartRing = false;
		bool reverseEndRing = false;

		if (m_IsSnapped && !endIsSnapped) {
			// 鍙湁璧风偣杩炴帴宸叉湁绔潰銆?
			// 濡傛灉璧风偣闇€瑕佸弽杞紝缁堢偣涔熷繀椤讳竴璧峰弽杞紝閬垮厤绠¤韩楹昏姳銆?
			reverseStartRing = startNeedsReverse;
			reverseEndRing = startNeedsReverse;
		}
		else if (!m_IsSnapped && endIsSnapped) {
			// 鍙湁缁堢偣杩炴帴宸叉湁绔潰銆?
			// 濡傛灉缁堢偣闇€瑕佸弽杞紝璧风偣涔熷繀椤讳竴璧峰弽杞紝閬垮厤绠¤韩楹昏姳銆?
			reverseStartRing = endNeedsReverse;
			reverseEndRing = endNeedsReverse;
		}
		else if (m_IsSnapped && endIsSnapped) {
			// 涓ょ閮借繛鎺ュ凡鏈夌闈€?
			// 杩欑鎯呭喌姣旇緝灏戣銆?
			// 鍏堝垎鍒弧瓒充袱绔帴鍙ｃ€?
			// 濡傛灉涓ょ闇€姹備笉涓€鑷达紝璇存槑杩欎袱涓闈箣闂存湰韬瓨鍦ㄥ崐鎵浆鍏崇郴銆?
			reverseStartRing = startNeedsReverse;
			reverseEndRing = endNeedsReverse;

			if (reverseStartRing != reverseEndRing) {
				if (false) qDebug() << "[控制笼连接] 端面环方向不一致";
			}
		}
		else {
			// 涓ょ閮戒笉鏄凡鏈夌闈紝淇濇寔棰勮缁撴灉銆?
			reverseStartRing = false;
			reverseEndRing = false;
		}

		double finalStartBeforeShift[4][3];
		double finalEndBeforeShift[4][3];

		CopyRing4MaybeReverse(finalStartBeforeShift, rawStart, reverseStartRing);
		CopyRing4MaybeReverse(finalEndBeforeShift, rawEnd, reverseEndRing);

		// ------------------------------------------------------------
		// 最终安全对齐：允许 4 种循环平移 + 4 种反向循环平移。
		// 旧逻辑在两端都 snap 且 reverseStartRing != reverseEndRing 时只打印
		// “[控制笼连接] 端面环方向不一致”，但仍然只做 cyclic shift，
		// 这会导致整体重建后局部麻花、sewing 留下 free wire。
		// 这里必须把手动重连入口也接入同一套 8 候选 ring 对齐。
		// ------------------------------------------------------------
		double finalEnd[4][3];
		AlignRingByEightCandidates(
			finalStartBeforeShift,
			finalEndBeforeShift,
			finalEnd,
			"[ManualRingBridgeAlign8]");

		// ------------------------------------------------------------
		// 鍐欏叆鏂?cage銆?
		// ------------------------------------------------------------
		CageSegment newCage;

		CopyRing4(newCage.pStart, finalStartBeforeShift);
		CopyRing4(newCage.pEnd, finalEnd);

		m_allCages.push_back(newCage);

		auto suppressSnappedSideFace = [&](int cageIndex, int faceOffset)
		{
			if (cageIndex < 0 || faceOffset < 2 || faceOffset > 5) {
				return;
			}

			const int sideFaceIndex = faceOffset - 2;
			for (const auto& item : m_suppressedCageSideFaces) {
				if (item.cageIndex == cageIndex && item.sideFaceIndex == sideFaceIndex) {
					return;
				}
			}

			m_suppressedCageSideFaces.push_back({ cageIndex, sideFaceIndex });
		};

		if (m_IsSnapped) {
			suppressSnappedSideFace(m_SnappedCageIndex, m_SnappedFaceOffset);
		}
		if (endIsSnapped) {
			suppressSnappedSideFace(endSnappedCageIndex, endSnappedFaceOffset);
		}

		RebuildGlobalMesh();

		if (m_FixedBoxActor) {
			m_FixedBoxActor->VisibilityOff();
		}

		if (m_LinesActor) {
			m_LinesActor->VisibilityOff();
		}

		m_ClickState = 0;
		m_ManualRotationOffset = 0;
		m_IsSnapped = false;
		m_SnappedIsStart = false;
		m_SnappedCageIndex = -1;
		m_SnappedFaceOffset = -1;
	}

	this->GetInteractor()->GetRenderWindow()->Render();
}
void CageInteractorStyle::OnRightButtonDown()
{
	if (m_isStraightCageEditMode && m_straightCageEditActive) {
		vtkRenderer* ren = nullptr;
		if (this->GetInteractor()) {
			int x = this->GetInteractor()->GetEventPosition()[0];
			int y = this->GetInteractor()->GetEventPosition()[1];
			ren = this->GetInteractor()->FindPokedRenderer(x, y);
		}
		CancelStraightCageEdit(ren);
		return;
	}

	if (m_isStraightCageRotateMode && m_straightCageRotateActive) {
		vtkRenderer* ren = nullptr;
		if (this->GetInteractor()) {
			int x = this->GetInteractor()->GetEventPosition()[0];
			int y = this->GetInteractor()->GetEventPosition()[1];
			ren = this->GetInteractor()->FindPokedRenderer(x, y);
		}
		CancelStraightCageRotation(ren);
		return;
	}

	if (m_isSkeletonExtendMode && m_selectedExtendSkelNode >= 0) {
		vtkRenderer* ren = nullptr;
		if (this->GetInteractor()) {
			int x = this->GetInteractor()->GetEventPosition()[0];
			int y = this->GetInteractor()->GetEventPosition()[1];
			ren = this->GetInteractor()->FindPokedRenderer(x, y);
		}
		ClearSkeletonExtendSelection(ren);
		qDebug() << "[骨架延伸] 已取消当前节点选择。";
		return;
	}

	if (m_ClickState == 1) {
		m_ClickState = 0;
		if (m_FixedBoxActor) m_FixedBoxActor->VisibilityOff();
		if (m_LinesActor) m_LinesActor->VisibilityOff();
		m_ManualRotationOffset = 0;
		this->GetInteractor()->GetRenderWindow()->Render();
		return;
	}

	vtkRenderWindowInteractor* interactor = this->GetInteractor();
	if (!interactor) return;
	int x = interactor->GetEventPosition()[0];
	int y = interactor->GetEventPosition()[1];

	vtkRenderer* ren = interactor->FindPokedRenderer(x, y);
	if (!ren) {
		vtkInteractorStyleTrackballCamera::OnRightButtonDown();
		return;
	}

	// 
	if (m_CageSolidActor) {
		vtkNew<vtkCellPicker> cagePicker;
		cagePicker->SetTolerance(0.005);
		m_CageSolidActor->SetPickable(1);
		m_CageSolidActor->SetVisibility(true);
		cagePicker->AddPickList(m_CageSolidActor);
		cagePicker->PickFromListOn();

		cagePicker->Pick(x, y, 0, ren);

		if (cagePicker->GetActor() == m_CageSolidActor) {
			vtkIdType cellId = cagePicker->GetCellId();
			if (cellId != -1) {
				int cageIndex = cellId / 6;
				if (cageIndex >= 0 && cageIndex < m_allCages.size()) {
					m_allCages.erase(m_allCages.begin() + cageIndex);
					for (auto it = m_suppressedCageSideFaces.begin(); it != m_suppressedCageSideFaces.end();) {
						if (it->cageIndex == cageIndex) {
							it = m_suppressedCageSideFaces.erase(it);
						}
						else {
							if (it->cageIndex > cageIndex) {
								--it->cageIndex;
							}
							++it;
						}
					}
					RebuildGlobalMesh();
					if (m_PreviewBoxActor) m_PreviewBoxActor->VisibilityOff();

					m_CageSolidActor->SetVisibility(false);
					m_CageSolidActor->SetPickable(0);

					this->GetInteractor()->GetRenderWindow()->Render();
					return;
				}
			}
		}
		m_CageSolidActor->SetVisibility(false);
		m_CageSolidActor->SetPickable(0);
	}

	vtkInteractorStyleTrackballCamera::OnRightButtonDown();
}
void CageInteractorStyle::OnMiddleButtonDown()
{
	if (m_ClickState == 1 && m_PreviewBoxActor && m_PreviewBoxActor->GetVisibility()) {
		m_ManualRotationOffset = (m_ManualRotationOffset + 1) % 4;
		this->OnMouseMove();
		return;
	}
	vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
}

void CageInteractorStyle::GetTransformCorners(vtkTransform* t, double corners[4][3])
{
	double local[4][4] = {
		{-0.5, -0.5, 0, 1}, { 0.5, -0.5, 0, 1},
		{ 0.5,  0.5, 0, 1}, {-0.5,  0.5, 0, 1}
	};
	for (int i = 0; i < 4; i++) {
		double world[4];
		t->MultiplyPoint(local[i], world);
		corners[i][0] = world[0] / world[3];
		corners[i][1] = world[1] / world[3];
		corners[i][2] = world[2] / world[3];
	}
}
void CageInteractorStyle::RebuildGlobalMesh()
{
	if (m_allCages.empty() && m_hubCages.empty()) return;

	m_controlVerts.clear();
	m_faceVerts.clear();
	m_vertsPerFace.clear();

	// ============================================================
	// 1. 全局点复用 / 焊接
	//    关键目的：
	//    hubPoly 的插入孔边界点 和 侧枝 CageSegment 的 pStart 点
	//    只要空间位置足够接近，就必须复用同一个 control vertex index。
	// ============================================================

	const double weldTol = 1e-3;
	const double weldTol2 = weldTol * weldTol;

	auto addOrFindVertex = [&](const double p[3]) -> int
	{
		int bestId = -1;
		double bestD2 = weldTol2;

		for (size_t i = 0; i < m_controlVerts.size(); ++i)
		{
			double dx = static_cast<double>(m_controlVerts[i].x) - p[0];
			double dy = static_cast<double>(m_controlVerts[i].y) - p[1];
			double dz = static_cast<double>(m_controlVerts[i].z) - p[2];

			double d2 = dx * dx + dy * dy + dz * dz;

			if (d2 <= bestD2) {
				bestD2 = d2;
				bestId = static_cast<int>(i);
			}
		}

		if (bestId >= 0) {
			return bestId;
		}

		m_controlVerts.push_back({
			static_cast<float>(p[0]),
			static_cast<float>(p[1]),
			static_cast<float>(p[2])
			});

		return static_cast<int>(m_controlVerts.size()) - 1;
	};

	auto getVertPoint = [&](int id, double p[3])
	{
		p[0] = m_controlVerts[id].x;
		p[1] = m_controlVerts[id].y;
		p[2] = m_controlVerts[id].z;
	};

	// ============================================================
	// 2. 面方向、重复面、退化面处理
	// ============================================================

	auto orientFaceOutward = [&](std::vector<int>& face, const double volumeCenter[3])
	{
		if (face.size() < 3) return;

		double p0[3], p1[3], p2[3];
		getVertPoint(face[0], p0);
		getVertPoint(face[1], p1);

		bool foundNonCollinear = false;

		for (size_t k = 2; k < face.size(); ++k) {
			getVertPoint(face[k], p2);

			double e1[3] = {
				p1[0] - p0[0],
				p1[1] - p0[1],
				p1[2] - p0[2]
			};

			double e2[3] = {
				p2[0] - p0[0],
				p2[1] - p0[1],
				p2[2] - p0[2]
			};

			double n[3];
			vtkMath::Cross(e1, e2, n);

			if (vtkMath::Norm(n) > 1e-10) {
				foundNonCollinear = true;
				break;
			}
		}

		if (!foundNonCollinear) return;

		double e1[3] = {
			p1[0] - p0[0],
			p1[1] - p0[1],
			p1[2] - p0[2]
		};

		double e2[3] = {
			p2[0] - p0[0],
			p2[1] - p0[1],
			p2[2] - p0[2]
		};

		double normal[3];
		vtkMath::Cross(e1, e2, normal);

		if (vtkMath::Norm(normal) < 1e-10) return;
		vtkMath::Normalize(normal);

		double faceCenter[3] = { 0.0, 0.0, 0.0 };

		for (int id : face) {
			double p[3];
			getVertPoint(id, p);

			faceCenter[0] += p[0] / static_cast<double>(face.size());
			faceCenter[1] += p[1] / static_cast<double>(face.size());
			faceCenter[2] += p[2] / static_cast<double>(face.size());
		}

		double outDir[3] = {
			faceCenter[0] - volumeCenter[0],
			faceCenter[1] - volumeCenter[1],
			faceCenter[2] - volumeCenter[2]
		};

		if (vtkMath::Dot(normal, outDir) < 0.0) {
			std::reverse(face.begin(), face.end());
		}
	};

	auto makeFaceKey = [](const std::vector<int>& face) -> std::set<int>
	{
		std::set<int> key;

		for (int id : face) {
			key.insert(id);
		}

		return key;
	};

	auto cleanFaceIds = [](const std::vector<int>& inFace) -> std::vector<int>
	{
		std::vector<int> out;

		for (int id : inFace) {
			if (id < 0) continue;

			if (!out.empty() && out.back() == id) {
				continue;
			}

			out.push_back(id);
		}

		if (out.size() >= 2 && out.front() == out.back()) {
			out.pop_back();
		}

		std::set<int> uniqueIds(out.begin(), out.end());

		if (uniqueIds.size() < 3) {
			out.clear();
		}

		return out;
	};

	struct PendingFace
	{
		std::vector<int> ids;
		std::set<int> key;
	};

	std::vector<PendingFace> pendingFaces;
	std::map<std::set<int>, int> faceUseCount;

	auto shouldSuppressCageSideFace = [&](int cageIndex, int sideFaceIndex) -> bool
	{
		for (const auto& item : m_suppressedCageSideFaces) {
			if (item.cageIndex == cageIndex && item.sideFaceIndex == sideFaceIndex) {
				return true;
			}
		}

		return false;
	};

	auto enqueueFace = [&](const std::vector<int>& rawFace)
	{
		std::vector<int> face = cleanFaceIds(rawFace);

		if (face.size() < 3) return;

		PendingFace pf;
		pf.ids = face;
		pf.key = makeFaceKey(face);

		if (pf.key.size() < 3) return;

		pendingFaces.push_back(pf);
		faceUseCount[pf.key]++;
	};

	// ============================================================
	// 3. 普通控制笼侧面加入待细分面
	// ============================================================

	for (size_t cageIndex = 0; cageIndex < m_allCages.size(); ++cageIndex)
	{
		const auto& cage = m_allCages[cageIndex];

		int v[8];

		for (int i = 0; i < 4; ++i) {
			v[i] = addOrFindVertex(cage.pStart[i]);
		}

		for (int i = 0; i < 4; ++i) {
			v[i + 4] = addOrFindVertex(cage.pEnd[i]);
		}

		double cageCenter[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				cageCenter[d] += cage.pStart[i][d] / 8.0;
				cageCenter[d] += cage.pEnd[i][d] / 8.0;
			}
		}

		std::vector<std::vector<int>> sideFaces;

		sideFaces.push_back({ v[0], v[1], v[5], v[4] });
		sideFaces.push_back({ v[1], v[2], v[6], v[5] });
		sideFaces.push_back({ v[2], v[3], v[7], v[6] });
		sideFaces.push_back({ v[3], v[0], v[4], v[7] });

		for (int sideFaceIndex = 0; sideFaceIndex < static_cast<int>(sideFaces.size()); ++sideFaceIndex) {
			if (shouldSuppressCageSideFace(static_cast<int>(cageIndex), sideFaceIndex)) {
				continue;
			}

			auto f = sideFaces[sideFaceIndex];
			orientFaceOutward(f, cageCenter);
			enqueueFace(f);
		}
	}

	// ============================================================
	// 4. 接头 patch / hubPoly 加入待细分面
	//    这里同样必须使用 addOrFindVertex，否则 hub 与 branch 不会共享拓扑点。
	// ============================================================

	for (const auto& hub : m_hubCages)
	{
		if (!hub || !hub->GetPoints() || !hub->GetPolys()) continue;

		vtkPoints* pts = hub->GetPoints();

		double hubCenter[3] = { 0.0, 0.0, 0.0 };
		int nPts = static_cast<int>(pts->GetNumberOfPoints());

		if (nPts <= 0) continue;

		std::vector<int> hubPtMap(nPts, -1);

		for (int i = 0; i < nPts; ++i) {
			double p[3];
			pts->GetPoint(i, p);

			hubCenter[0] += p[0] / static_cast<double>(nPts);
			hubCenter[1] += p[1] / static_cast<double>(nPts);
			hubCenter[2] += p[2] / static_cast<double>(nPts);

			hubPtMap[i] = addOrFindVertex(p);
		}

		vtkCellArray* polys = hub->GetPolys();
		polys->InitTraversal();

		vtkNew<vtkIdList> idList;

		while (polys->GetNextCell(idList))
		{
			int nV = idList->GetNumberOfIds();
			if (nV < 3) continue;

			std::vector<int> face;
			face.reserve(nV);

			bool valid = true;

			for (int i = 0; i < nV; ++i) {
				vtkIdType localId = idList->GetId(i);

				if (localId < 0 || localId >= static_cast<vtkIdType>(hubPtMap.size())) {
					valid = false;
					break;
				}

				int gid = hubPtMap[static_cast<int>(localId)];

				if (gid < 0) {
					valid = false;
					break;
				}

				face.push_back(gid);
			}

			if (!valid) continue;

			orientFaceOutward(face, hubCenter);
			enqueueFace(face);
		}
	}

	// ============================================================
	// 5. 删除完全重复的内部面，只把外露面送入 OpenSubdiv
	// ============================================================

	for (const auto& pf : pendingFaces)
	{
		auto it = faceUseCount.find(pf.key);
		if (it == faceUseCount.end()) continue;

		if (it->second != 1) {
			continue;
		}

		std::vector<int> cleanFace = cleanFaceIds(pf.ids);

		if (cleanFace.size() < 3) {
			continue;
		}

		m_vertsPerFace.push_back(static_cast<int>(cleanFace.size()));

		for (int id : cleanFace) {
			m_faceVerts.push_back(id);
		}
	}

	if (m_faceVerts.empty()) return;

	// ============================================================
// Debug: 显示进入 OpenSubdiv 前的非四边 coarse face
// 黄色线 = outfile 会跳过的非四边面边界
// 如果黄色线出现在插入区域，那么红线来源就是 skipped non-quad face。
// ============================================================
	{
		vtkNew<vtkPoints> nonQuadPts;
		vtkNew<vtkCellArray> nonQuadLines;

		size_t cursor = 0;

		int triFaceCount = 0;
		int quadFaceCount = 0;
		int ngonFaceCount = 0;
		int otherFaceCount = 0;

		for (size_t f = 0; f < m_vertsPerFace.size(); ++f) {
			int n = m_vertsPerFace[f];

			if (cursor + n > m_faceVerts.size()) {
				break;
			}

			if (n == 3) {
				triFaceCount++;
			}
			else if (n == 4) {
				quadFaceCount++;
			}
			else if (n > 4) {
				ngonFaceCount++;
			}
			else {
				otherFaceCount++;
			}

			if (n != 4 && n >= 3) {
				std::vector<vtkIdType> localIds;
				localIds.reserve(n);

				for (int i = 0; i < n; ++i) {
					int vid = m_faceVerts[cursor + i];

					if (vid < 0 || vid >= static_cast<int>(m_controlVerts.size())) {
						continue;
					}

					vtkIdType pid = nonQuadPts->InsertNextPoint(
						m_controlVerts[vid].x,
						m_controlVerts[vid].y,
						m_controlVerts[vid].z
					);

					localIds.push_back(pid);
				}

				if (static_cast<int>(localIds.size()) == n) {
					for (int i = 0; i < n; ++i) {
						vtkIdType line[2] = {
							localIds[i],
							localIds[(i + 1) % n]
						};

						nonQuadLines->InsertNextCell(2, line);
					}
				}
			}

			cursor += n;
		}

		if (false) qDebug() << "[TopologyDebug] coarse face type before OpenSubdiv:"
			<< "tri =" << triFaceCount
			<< "quad =" << quadFaceCount
			<< "ngon =" << ngonFaceCount
			<< "other =" << otherFaceCount;

		if (nonQuadPts->GetNumberOfPoints() > 0) {
			vtkNew<vtkPolyData> nonQuadPoly;
			nonQuadPoly->SetPoints(nonQuadPts);
			nonQuadPoly->SetLines(nonQuadLines);

			vtkNew<vtkTubeFilter> tube;
			tube->SetInputData(nonQuadPoly);
			tube->SetRadius(0.08);
			tube->SetNumberOfSides(8);
			tube->CappingOn();
			tube->Update();

			vtkNew<vtkPolyDataMapper> mapper;
			mapper->SetInputConnection(tube->GetOutputPort());

			static vtkSmartPointer<vtkActor> nonQuadDebugActor = nullptr;

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(1.0, 1.0, 0.0);
			actor->GetProperty()->SetOpacity(1.0);
			actor->GetProperty()->SetLighting(false);
			actor->SetPickable(0);

			if (this->GetCurrentRenderer()) {
				if (nonQuadDebugActor) {
					this->GetCurrentRenderer()->RemoveActor(nonQuadDebugActor);
				}

				nonQuadDebugActor = actor;
				this->GetCurrentRenderer()->AddActor(nonQuadDebugActor);

				if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
					this->GetInteractor()->GetRenderWindow()->Render();
				}
			}

			if (false) qDebug() << "[TopologyDebug] non-quad coarse faces displayed in yellow."
				<< "points =" << nonQuadPts->GetNumberOfPoints()
				<< "lines =" << nonQuadLines->GetNumberOfCells();
		}
	}

	// ============================================================
	// 5.5 进入 OpenSubdiv 前的控制网格拓扑检查
	//     boundary edges 包括真实管口边界，但如果插入区没有焊接，也会出现在这里。
	// ============================================================

	{
		std::map<std::pair<int, int>, int> edgeUseCount;

		size_t cursor = 0;

		for (size_t f = 0; f < m_vertsPerFace.size(); ++f) {
			int n = m_vertsPerFace[f];

			if (n < 2 || cursor + n > m_faceVerts.size()) {
				cursor += n;
				continue;
			}

			for (int i = 0; i < n; ++i) {
				int a = m_faceVerts[cursor + i];
				int b = m_faceVerts[cursor + ((i + 1) % n)];

				if (a == b) continue;

				if (a > b) {
					std::swap(a, b);
				}

				edgeUseCount[{a, b}]++;
			}

			cursor += n;
		}

		int boundaryEdgeCount = 0;
		int nonManifoldEdgeCount = 0;

		for (const auto& kv : edgeUseCount) {
			if (kv.second == 1) {
				boundaryEdgeCount++;
			}
			else if (kv.second > 2) {
				nonManifoldEdgeCount++;
			}
		}

		if (false) qDebug() << "[TopologyDebug] RebuildGlobalMesh final control topology:"
			<< "control points =" << m_controlVerts.size()
			<< "faces =" << m_vertsPerFace.size()
			<< "boundary edges =" << boundaryEdgeCount
			<< "nonmanifold edges =" << nonManifoldEdgeCount
			<< "weldTol =" << weldTol;
	}


	// ============================================================
	// 6. OpenSubdiv 细分
	// ============================================================

	Far::TopologyDescriptor desc;
	desc.numVertices = static_cast<int>(m_controlVerts.size());
	desc.numFaces = static_cast<int>(m_vertsPerFace.size());
	desc.numVertsPerFace = m_vertsPerFace.data();
	desc.vertIndicesPerFace = m_faceVerts.data();

	if (m_refiner) {
		delete m_refiner;
		m_refiner = nullptr;
	}

	Sdc::Options osOptions;
	osOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);

	m_refiner = Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(
		desc,
		Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options(
			Sdc::SCHEME_CATMARK,
			osOptions
		)
	);

	// OpenSubdiv subdivision level.
	int maxLevel = 4;

	if (maxLevel > 0) {
		m_refiner->RefineUniform(
			Far::TopologyRefiner::UniformOptions(maxLevel)
		);
	}

	int numBaseVerts = m_refiner->GetLevel(0).GetNumVertices();

	m_refinedVerts.resize(m_refiner->GetNumVerticesTotal());

	for (int i = 0; i < numBaseVerts; ++i) {
		m_refinedVerts[i] = m_controlVerts[i];
	}

	if (maxLevel > 0) {
		Far::StencilTable const* st =
			Far::StencilTableFactory::Create(
				*m_refiner,
				Far::StencilTableFactory::Options()
			);

		if (st) {
			st->UpdateValues(
				m_controlVerts.data(),
				m_refinedVerts.data() + numBaseVerts
			);

			delete st;
		}
	}

	// ============================================================
	// 7. 生成细分后的 VTK 网格
	// ============================================================

	vtkNew<vtkPoints> outPts;

	int lastLevelOffset = 0;

	for (int i = 0; i < maxLevel; ++i) {
		lastLevelOffset += m_refiner->GetLevel(i).GetNumVertices();
	}

	Far::TopologyLevel const& lastLevel = m_refiner->GetLevel(maxLevel);

	for (int i = 0; i < lastLevel.GetNumVertices(); ++i) {
		OsdVertex vtx = m_refinedVerts[lastLevelOffset + i];
		outPts->InsertNextPoint(vtx.x, vtx.y, vtx.z);
	}

	vtkNew<vtkCellArray> outPolys;

	for (int f = 0; f < lastLevel.GetNumFaces(); ++f) {
		Far::ConstIndexArray fv = lastLevel.GetFaceVertices(f);

		outPolys->InsertNextCell(fv.size());

		for (int i = 0; i < fv.size(); ++i) {
			outPolys->InsertCellPoint(fv[i]);
		}
	}

	vtkNew<vtkPolyData> rawMesh;
	rawMesh->SetPoints(outPts);
	rawMesh->SetPolys(outPolys);

	// ============================================================
	// 8. 根据边界环生成显示用端盖
	//    这些端盖只用于 m_subD_Mesh 显示，不进入 OpenSubdiv coarse topology。
	// ============================================================

	vtkNew<vtkFeatureEdges> edges;
	edges->SetInputData(rawMesh);
	edges->BoundaryEdgesOn();
	edges->FeatureEdgesOff();
	edges->NonManifoldEdgesOff();
	edges->ManifoldEdgesOff();
	edges->Update();

	vtkNew<vtkStripper> stripper;
	stripper->SetInputConnection(edges->GetOutputPort());
	stripper->Update();

	vtkPolyData* boundaryLines = stripper->GetOutput();

	vtkNew<vtkPolyData> capPolyData;
	capPolyData->SetPoints(boundaryLines->GetPoints());

	vtkNew<vtkCellArray> capPolys;
	vtkCellArray* lines = boundaryLines->GetLines();

	if (lines) {
		lines->InitTraversal();

		vtkNew<vtkIdList> idListCap;

		while (lines->GetNextCell(idListCap)) {
			if (idListCap->GetNumberOfIds() >= 3) {
				capPolys->InsertNextCell(idListCap);
			}
		}
	}

	capPolyData->SetPolys(capPolys);

	vtkNew<vtkAppendPolyData> append;
	append->AddInputData(rawMesh);
	append->AddInputData(capPolyData);
	append->Update();

	vtkNew<vtkCleanPolyData> clean;
	clean->SetInputConnection(append->GetOutputPort());
	clean->Update();

	if (!m_subD_Mesh) {
		m_subD_Mesh = vtkSmartPointer<vtkPolyData>::New();
	}

	m_subD_Mesh->DeepCopy(clean->GetOutput());

	vtkNew<vtkTriangleFilter> tri;
	tri->SetInputData(m_subD_Mesh);

	vtkNew<vtkPolyDataNormals> norms;
	norms->SetInputConnection(tri->GetOutputPort());
	norms->ComputePointNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();

	if (!m_GlobalSmoothActor) {
		m_GlobalSmoothActor = vtkSmartPointer<vtkActor>::New();
		m_GlobalSmoothActor->SetMapper(vtkNew<vtkPolyDataMapper>());
		m_GlobalSmoothActor->GetProperty()->SetColor(0.82, 0.82, 0.78);
		m_GlobalSmoothActor->GetProperty()->SetSpecular(0.35);
		m_GlobalSmoothActor->GetProperty()->SetSpecularPower(28.0);

		if (this->GetCurrentRenderer()) {
			this->GetCurrentRenderer()->AddActor(m_GlobalSmoothActor);
		}
	}

	vtkPolyDataMapper::SafeDownCast(
		m_GlobalSmoothActor->GetMapper()
	)->SetInputData(norms->GetOutput());

	// ============================================================
	// 9. 绘制历史控制笼线
	// ============================================================

	vtkNew<vtkPoints> cagePts;
	vtkNew<vtkCellArray> cageLines;

	vtkIdType ptIdx = 0;

	for (const auto& cage : m_allCages)
	{
		for (int i = 0; i < 4; ++i) {
			cagePts->InsertNextPoint(cage.pStart[i]);
		}

		for (int i = 0; i < 4; ++i) {
			cagePts->InsertNextPoint(cage.pEnd[i]);
		}

		for (int i = 0; i < 4; ++i) {
			vtkIdType l1[2] = {
				ptIdx + i,
				ptIdx + (i + 1) % 4
			};

			vtkIdType l2[2] = {
				ptIdx + 4 + i,
				ptIdx + 4 + (i + 1) % 4
			};

			vtkIdType l3[2] = {
				ptIdx + i,
				ptIdx + 4 + i
			};

			cageLines->InsertNextCell(2, l1);
			cageLines->InsertNextCell(2, l2);
			cageLines->InsertNextCell(2, l3);
		}

		ptIdx += 8;
	}

	for (const auto& hub : m_hubCages)
	{
		if (!hub || !hub->GetPoints() || !hub->GetPolys()) continue;

		vtkPoints* hpts = hub->GetPoints();
		vtkCellArray* hpolys = hub->GetPolys();

		vtkIdType hubStartIdx = cagePts->GetNumberOfPoints();

		for (int i = 0; i < hpts->GetNumberOfPoints(); ++i) {
			double p[3];
			hpts->GetPoint(i, p);
			cagePts->InsertNextPoint(p);
		}

		hpolys->InitTraversal();

		vtkNew<vtkIdList> cidList;

		while (hpolys->GetNextCell(cidList)) {
			int nV = cidList->GetNumberOfIds();

			for (int i = 0; i < nV; ++i) {
				vtkIdType line[2] = {
					hubStartIdx + cidList->GetId(i),
					hubStartIdx + cidList->GetId((i + 1) % nV)
				};

				cageLines->InsertNextCell(2, line);
			}
		}
	}

	vtkNew<vtkPolyData> cageData;
	cageData->SetPoints(cagePts);
	cageData->SetLines(cageLines);

	if (!m_HistoryCagesActor) {
		m_HistoryCagesActor = vtkSmartPointer<vtkActor>::New();
		m_HistoryCagesActor->SetMapper(vtkNew<vtkPolyDataMapper>());
		m_HistoryCagesActor->GetProperty()->SetColor(0.78, 0.78, 0.74);
		m_HistoryCagesActor->GetProperty()->SetLineWidth(2.0);

		if (this->GetCurrentRenderer()) {
			this->GetCurrentRenderer()->AddActor(m_HistoryCagesActor);
		}
	}

	vtkPolyDataMapper::SafeDownCast(
		m_HistoryCagesActor->GetMapper()
	)->SetInputData(cageData);

	// ============================================================
	// 10. 控制笼实体预览，仅用于调试显示，不参与 OpenSubdiv。
	// ============================================================

	vtkNew<vtkPoints> solidPts;
	vtkNew<vtkCellArray> solidPolys;

	auto orientLocalSolidFace = [](
		std::vector<int>& face,
		const double localPts[8][3],
		const double center[3])
	{
		if (face.size() < 3) return;

		const double* p0 = localPts[face[0]];
		const double* p1 = localPts[face[1]];
		const double* p2 = localPts[face[2]];

		double e1[3] = {
			p1[0] - p0[0],
			p1[1] - p0[1],
			p1[2] - p0[2]
		};

		double e2[3] = {
			p2[0] - p0[0],
			p2[1] - p0[1],
			p2[2] - p0[2]
		};

		double normal[3];
		vtkMath::Cross(e1, e2, normal);

		if (vtkMath::Norm(normal) < 1e-10) return;
		vtkMath::Normalize(normal);

		double fc[3] = { 0.0, 0.0, 0.0 };

		for (int id : face) {
			fc[0] += localPts[id][0] / static_cast<double>(face.size());
			fc[1] += localPts[id][1] / static_cast<double>(face.size());
			fc[2] += localPts[id][2] / static_cast<double>(face.size());
		}

		double outDir[3] = {
			fc[0] - center[0],
			fc[1] - center[1],
			fc[2] - center[2]
		};

		if (vtkMath::Dot(normal, outDir) < 0.0) {
			std::reverse(face.begin(), face.end());
		}
	};

	vtkIdType off = 0;

	for (const auto& cage : m_allCages)
	{
		double localPts[8][3];

		for (int i = 0; i < 4; ++i) {
			memcpy(localPts[i], cage.pStart[i], sizeof(double) * 3);
			solidPts->InsertNextPoint(localPts[i]);
		}

		for (int i = 0; i < 4; ++i) {
			memcpy(localPts[i + 4], cage.pEnd[i], sizeof(double) * 3);
			solidPts->InsertNextPoint(localPts[i + 4]);
		}

		double c[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 8; ++i) {
			c[0] += localPts[i][0] / 8.0;
			c[1] += localPts[i][1] / 8.0;
			c[2] += localPts[i][2] / 8.0;
		}

		std::vector<std::vector<int>> solidFaces;

		solidFaces.push_back({ 0, 1, 2, 3 });
		solidFaces.push_back({ 4, 5, 6, 7 });
		solidFaces.push_back({ 0, 1, 5, 4 });
		solidFaces.push_back({ 1, 2, 6, 5 });
		solidFaces.push_back({ 2, 3, 7, 6 });
		solidFaces.push_back({ 3, 0, 4, 7 });

		for (auto& f : solidFaces) {
			orientLocalSolidFace(f, localPts, c);

			solidPolys->InsertNextCell(static_cast<vtkIdType>(f.size()));

			for (int id : f) {
				solidPolys->InsertCellPoint(off + id);
			}
		}

		off += 8;
	}

	vtkNew<vtkPolyData> sd;
	sd->SetPoints(solidPts);
	sd->SetPolys(solidPolys);

	DebugPolyDataFaceTypes(sd, "RebuildGlobalMesh cage solid preview");
	DebugPolyDataBoundary(sd, "RebuildGlobalMesh cage solid preview");

	if (!m_CageSolidActor) {
		m_CageSolidActor = vtkSmartPointer<vtkActor>::New();
		m_CageSolidActor->SetMapper(vtkNew<vtkPolyDataMapper>());
		m_CageSolidActor->GetProperty()->SetColor(0.82, 0.82, 0.78);
		m_CageSolidActor->GetProperty()->SetOpacity(0.15);
		m_CageSolidActor->GetProperty()->SetRepresentationToSurface();

		if (this->GetCurrentRenderer()) {
			this->GetCurrentRenderer()->AddActor(m_CageSolidActor);
		}
	}

	vtkPolyDataMapper::SafeDownCast(
		m_CageSolidActor->GetMapper()
	)->SetInputData(sd);

	m_CageSolidActor->VisibilityOff();
	m_CageSolidActor->SetPickable(0);

}

static int HbOccAccessViolationFilter(DWORD code, DWORD* sehCode)
{
	if (code == EXCEPTION_ACCESS_VIOLATION || code == EXCEPTION_IN_PAGE_ERROR) {
		if (sehCode) {
			*sehCode = code;
		}
		return EXCEPTION_EXECUTE_HANDLER;
	}

	return EXCEPTION_CONTINUE_SEARCH;
}

static bool SafeBuildFillingWithSEH(BRepOffsetAPI_MakeFilling& filler, DWORD& sehCode)
{
	sehCode = 0;

	__try
	{
		filler.Build();
		return true;
	}
	__except (HbOccAccessViolationFilter(GetExceptionCode(), &sehCode))
	{
		return false;
	}
}

static QString HbSehCodeText(DWORD sehCode)
{
	return QString("0x%1")
		.arg(static_cast<unsigned int>(sehCode), 8, 16, QChar('0'));
}

static bool SafeAddFillingEdgeWithSEH(
	BRepOffsetAPI_MakeFilling* filler,
	const TopoDS_Edge* edge,
	DWORD& sehCode)
{
	sehCode = 0;

	__try
	{
		filler->Add(*edge, GeomAbs_C0);
		return true;
	}
	__except (HbOccAccessViolationFilter(GetExceptionCode(), &sehCode))
	{
		return false;
	}
}

static Handle(Geom_BSplineSurface) HbBuildClampedBSplineSurface(
	const TColgp_Array2OfPnt& poles,
	int uPoleCount,
	int vPoleCount)
{
	const int uDegree = std::min(3, std::max(1, uPoleCount - 1));
	const int vDegree = std::min(3, std::max(1, vPoleCount - 1));

	if (uPoleCount < uDegree + 1 || vPoleCount < vDegree + 1) {
		return Handle(Geom_BSplineSurface)();
	}

	const int uKnotCount = uPoleCount - uDegree + 1;
	const int vKnotCount = vPoleCount - vDegree + 1;

	TColStd_Array1OfReal uKnots(1, uKnotCount);
	TColStd_Array1OfReal vKnots(1, vKnotCount);
	TColStd_Array1OfInteger uMults(1, uKnotCount);
	TColStd_Array1OfInteger vMults(1, vKnotCount);

	for (int i = 1; i <= uKnotCount; ++i) {
		const double t = (uKnotCount == 1)
			? 0.0
			: static_cast<double>(i - 1) / static_cast<double>(uKnotCount - 1);
		uKnots.SetValue(i, t);
		uMults.SetValue(i, (i == 1 || i == uKnotCount) ? uDegree + 1 : 1);
	}

	for (int i = 1; i <= vKnotCount; ++i) {
		const double t = (vKnotCount == 1)
			? 0.0
			: static_cast<double>(i - 1) / static_cast<double>(vKnotCount - 1);
		vKnots.SetValue(i, t);
		vMults.SetValue(i, (i == 1 || i == vKnotCount) ? vDegree + 1 : 1);
	}

	try {
		return new Geom_BSplineSurface(
			poles,
			uKnots,
			vKnots,
			uMults,
			vMults,
			uDegree,
			vDegree,
			Standard_False,
			Standard_False);
	}
	catch (...) {
		return Handle(Geom_BSplineSurface)();
	}
}

namespace
{
	static double HbExportBoundsDiag(const double b[6])
	{
		const double dx = b[1] - b[0];
		const double dy = b[3] - b[2];
		const double dz = b[5] - b[4];
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	static bool HbExportPointInBounds(const std::array<double, 3>& p, const double b[6])
	{
		return p[0] >= b[0] && p[0] <= b[1] &&
			p[1] >= b[2] && p[1] <= b[3] &&
			p[2] >= b[4] && p[2] <= b[5];
	}

	static bool HbExportGetExpandedShapeBounds(const TopoDS_Shape& shape, double outB[6])
	{
		if (shape.IsNull()) {
			return false;
		}

		Bnd_Box box;
		BRepBndLib::Add(shape, box);

		if (box.IsVoid()) {
			return false;
		}

		double xmin = 0.0, ymin = 0.0, zmin = 0.0;
		double xmax = 0.0, ymax = 0.0, zmax = 0.0;
		box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

		outB[0] = xmin;
		outB[1] = xmax;
		outB[2] = ymin;
		outB[3] = ymax;
		outB[4] = zmin;
		outB[5] = zmax;

		const double eps = std::max(1e-6, HbExportBoundsDiag(outB) * 1e-4);
		outB[0] -= eps;
		outB[1] += eps;
		outB[2] -= eps;
		outB[3] += eps;
		outB[4] -= eps;
		outB[5] += eps;

		return true;
	}

	static std::vector<TopoDS_Shape> HbExportLoadMatchedNonDesignSolids(
		const QString& stepPath,
		const std::vector<std::array<double, 3>>& elementCenters,
		const std::vector<char>& elementIsNonDesign,
		bool allStepSolidsAreNonDesign)
	{
		std::vector<TopoDS_Shape> matchedSolids;

		if (stepPath.isEmpty() || (!allStepSolidsAreNonDesign &&
			(elementCenters.empty() || elementCenters.size() != elementIsNonDesign.size()))) {
			qWarning() << "[导出融合] non-design 匹配输入为空或尺寸不一致。";
			return matchedSolids;
		}

		STEPControl_Reader reader;
		IFSelect_ReturnStatus status = reader.ReadFile(stepPath.toStdString().c_str());
		if (status != IFSelect_RetDone) {
			qWarning() << "[导出融合] STEP 读取失败:" << stepPath;
			return matchedSolids;
		}

		reader.TransferRoots();
		TopoDS_Shape root = reader.OneShape();
		if (root.IsNull()) {
			qWarning() << "[导出融合] STEP root shape 为空。";
			return matchedSolids;
		}

		std::vector<TopoDS_Shape> solids;
		for (TopExp_Explorer exp(root, TopAbs_SOLID); exp.More(); exp.Next()) {
			solids.push_back(exp.Current());
		}

		qDebug() << "[导出融合] STEP solid 数量 =" << static_cast<int>(solids.size());
		if (allStepSolidsAreNonDesign) {
			qDebug() << "[导出融合] 当前 STEP 仅包含 non-design，保留全部 solid。";
			return solids;
		}

		for (size_t si = 0; si < solids.size(); ++si) {
			double b[6];
			if (!HbExportGetExpandedShapeBounds(solids[si], b)) {
				continue;
			}

			int ndInside = 0;
			int designInside = 0;
			for (size_t i = 0; i < elementCenters.size(); ++i) {
				if (!HbExportPointInBounds(elementCenters[i], b)) {
					continue;
				}

				if (elementIsNonDesign[i]) {
					++ndInside;
				}
				else {
					++designInside;
				}
			}

			const int totalInside = ndInside + designInside;
			const double ndRatio = totalInside > 0
				? static_cast<double>(ndInside) / static_cast<double>(totalInside)
				: 0.0;

			const bool selected = (ndInside > 0 && ndRatio >= 0.50);

			qDebug() << "[导出融合] solid" << static_cast<int>(si)
				<< "ndInside =" << ndInside
				<< "designInside =" << designInside
				<< "ndRatio =" << ndRatio
				<< "selected =" << selected;

			if (selected) {
				matchedSolids.push_back(solids[si]);
			}
		}

		qDebug() << "[导出融合] 匹配到 non-design solid 数量 =" << static_cast<int>(matchedSolids.size());
		return matchedSolids;
	}

	static TopoDS_Shape HbExportMakeCompound(
		const TopoDS_Shape& rebuiltShape,
		const std::vector<TopoDS_Shape>& nonDesignSolids)
	{
		BRep_Builder builder;
		TopoDS_Compound compound;
		builder.MakeCompound(compound);

		if (!rebuiltShape.IsNull()) {
			builder.Add(compound, rebuiltShape);
		}

		for (const TopoDS_Shape& nd : nonDesignSolids) {
			if (!nd.IsNull()) {
				builder.Add(compound, nd);
			}
		}

		return compound;
	}

	static double HbExportShapeVolume(const TopoDS_Shape& shape)
	{
		if (shape.IsNull()) {
			return -1.0;
		}

		try {
			GProp_GProps props;
			BRepGProp::VolumeProperties(shape, props);
			const double v = props.Mass();
			return std::isfinite(v) ? std::abs(v) : -1.0;
		}
		catch (...) {
			return -1.0;
		}
	}

	static bool HbExportShapeBoundsContainsShape(
		const TopoDS_Shape& containerShape,
		const TopoDS_Shape& childShape,
		double relativeTol = 1e-4)
	{
		if (containerShape.IsNull() || childShape.IsNull()) {
			return false;
		}

		double cb[6];
		double nb[6];
		if (!HbExportGetExpandedShapeBounds(containerShape, cb) ||
			!HbExportGetExpandedShapeBounds(childShape, nb)) {
			return false;
		}

		const double tol = std::max({ 1e-6, HbExportBoundsDiag(cb) * relativeTol, HbExportBoundsDiag(nb) * relativeTol });
		return cb[0] <= nb[0] + tol && cb[1] >= nb[1] - tol &&
			cb[2] <= nb[2] + tol && cb[3] >= nb[3] - tol &&
			cb[4] <= nb[4] + tol && cb[5] >= nb[5] - tol;
	}

	static TopoDS_Shape HbExportHealShapeForBoolean(const TopoDS_Shape& shape)
	{
		if (shape.IsNull()) {
			return shape;
		}

		try {
			ShapeFix_Shape fixer(shape);
			fixer.Perform();
			TopoDS_Shape fixed = fixer.Shape();
			return fixed.IsNull() ? shape : fixed;
		}
		catch (...) {
			return shape;
		}
	}

	static bool HbExportTryFuseOnce(
		const TopoDS_Shape& currentShape,
		const TopoDS_Shape& ndShape,
		double fuzzyValue,
		bool healInputs,
		TopoDS_Shape& outShape,
		QString& failReason)
	{
		outShape.Nullify();
		failReason.clear();

		if (currentShape.IsNull() || ndShape.IsNull()) {
			failReason = QString::fromUtf8(u8"输入 shape 为空");
			return false;
		}

		TopoDS_Shape a = healInputs ? HbExportHealShapeForBoolean(currentShape) : currentShape;
		TopoDS_Shape b = healInputs ? HbExportHealShapeForBoolean(ndShape) : ndShape;

		try {
			BRepAlgoAPI_Fuse fuse(a, b);
			if (fuzzyValue > 0.0) {
				fuse.SetFuzzyValue(fuzzyValue);
			}
			fuse.Build();

			if (!fuse.IsDone()) {
				failReason = QString::fromUtf8(u8"BRepAlgoAPI_Fuse::IsDone=false");
				return false;
			}

			TopoDS_Shape result = fuse.Shape();
			if (result.IsNull()) {
				failReason = QString::fromUtf8(u8"Fuse 返回空 shape");
				return false;
			}

			// 对布尔结果再做一次轻量修复。注意：这里不重新提取某个 solid，避免把 non-design 子实体丢掉。
			result = HbExportHealShapeForBoolean(result);
			if (result.IsNull()) {
				failReason = QString::fromUtf8(u8"Fuse 后 ShapeFix 返回空 shape");
				return false;
			}

			outShape = result;
			return true;
		}
		catch (const Standard_Failure& e) {
			failReason = QString::fromUtf8(u8"OCC 异常：%1").arg(QString::fromLocal8Bit(e.GetMessageString()));
			return false;
		}
		catch (...) {
			failReason = QString::fromUtf8(u8"未知异常");
			return false;
		}
	}

	static bool HbExportAcceptTrueUnionResult(
		const TopoDS_Shape& candidate,
		const TopoDS_Shape& currentBefore,
		const TopoDS_Shape& ndShape,
		double currentVolumeBefore,
		double ndVolume,
		QString& rejectReason)
	{
		rejectReason.clear();

		if (candidate.IsNull()) {
			rejectReason = QString::fromUtf8(u8"candidate 为空");
			return false;
		}

		BRepCheck_Analyzer analyzer(candidate);
		const bool validShape = analyzer.IsValid();

		// 并集结果必须覆盖被融合的 non-design 的空间包围盒。
		// 这条检查直接针对你遇到的“插入部分被吞掉，但 non-design 外壳没保留”的问题。
		const bool boundsKeepNd = HbExportShapeBoundsContainsShape(candidate, ndShape, 2e-3);

		const double fusedVolume = HbExportShapeVolume(candidate);

		bool volumeOk = true;
		if (currentVolumeBefore > 0.0 && fusedVolume > 0.0) {
			// 真正的并集不应该比当前实体明显变小。
			// 允许 1% 的 OCC 体积计算波动，但不接受明显丢体。
			volumeOk = fusedVolume >= currentVolumeBefore * 0.99;
		}

		bool hasReasonableVolume = true;
		if (ndVolume > 0.0 && currentVolumeBefore > 0.0 && fusedVolume > 0.0) {
			// 如果 non-design 完全没有进入结果，通常 boundsKeepNd 会失败；这里再加一道体积下限。
			// 由于设计域和 non-design 可能有较大重叠，不能要求 +ndVolume，只要求不能明显低于当前体积。
			hasReasonableVolume = fusedVolume >= std::max(currentVolumeBefore, ndVolume) * 0.99;
		}

		if (!validShape || !boundsKeepNd || !volumeOk || !hasReasonableVolume) {
			rejectReason = QString::fromUtf8(
				u8"valid=%1, boundsKeepNd=%2, volumeOk=%3, volumeReasonable=%4, volBefore=%5, ndVol=%6, fusedVol=%7")
				.arg(validShape)
				.arg(boundsKeepNd)
				.arg(volumeOk)
				.arg(hasReasonableVolume)
				.arg(currentVolumeBefore)
				.arg(ndVolume)
				.arg(fusedVolume);
			return false;
		}

		return true;
	}

	static TopoDS_Shape HbExportTryBooleanFuseWithNonDesign(
		const TopoDS_Shape& rebuiltShape,
		const QString& stepPath,
		const std::vector<std::array<double, 3>>& elementCenters,
		const std::vector<char>& elementIsNonDesign,
		bool allStepSolidsAreNonDesign,
		bool& usedFuse,
		bool& usedCompoundFallback,
		QString& report)
	{
		usedFuse = false;
		usedCompoundFallback = false;
		report.clear();

		if (rebuiltShape.IsNull()) {
			report = QString::fromUtf8(u8"重建实体为空，不能融合。");
			usedCompoundFallback = true;
			return rebuiltShape;
		}

		std::vector<TopoDS_Shape> nonDesignSolids = HbExportLoadMatchedNonDesignSolids(
			stepPath,
			elementCenters,
			elementIsNonDesign,
			allStepSolidsAreNonDesign);

		if (nonDesignSolids.empty()) {
			report = QString::fromUtf8(u8"没有匹配到 STP non-design 实体，保持原重建实体导出。");
			return rebuiltShape;
		}

		TopoDS_Shape current = HbExportHealShapeForBoolean(rebuiltShape);
		if (current.IsNull()) {
			current = rebuiltShape;
		}

		int successCount = 0;
		QStringList failureReports;

		for (size_t i = 0; i < nonDesignSolids.size(); ++i) {
			TopoDS_Shape ndShape = HbExportHealShapeForBoolean(nonDesignSolids[i]);
			if (ndShape.IsNull()) {
				failureReports << QString::fromUtf8(u8"index=%1：non-design shape 为空").arg(static_cast<int>(i));
				usedCompoundFallback = true;
				break;
			}

			const double currentVolumeBefore = HbExportShapeVolume(current);
			const double ndVolume = HbExportShapeVolume(ndShape);

			double cb[6];
			double nb[6];
			const bool hasCurrentBounds = HbExportGetExpandedShapeBounds(current, cb);
			const bool hasNdBounds = HbExportGetExpandedShapeBounds(ndShape, nb);
			const double diag = (hasCurrentBounds && hasNdBounds)
				? std::max(HbExportBoundsDiag(cb), HbExportBoundsDiag(nb))
				: 1.0;

			std::vector<std::pair<double, bool>> attempts;
			attempts.push_back({ 0.0, false });
			attempts.push_back({ std::max(1e-7, diag * 1e-7), false });
			attempts.push_back({ std::max(1e-6, diag * 1e-6), false });
			attempts.push_back({ std::max(5e-6, diag * 5e-6), true });
			attempts.push_back({ std::max(1e-5, diag * 1e-5), true });

			bool accepted = false;
			TopoDS_Shape acceptedShape;
			QStringList attemptRejects;

			for (size_t ai = 0; ai < attempts.size(); ++ai) {
				TopoDS_Shape candidate;
				QString failReason;
				const double fuzzy = attempts[ai].first;
				const bool healInputs = attempts[ai].second;

				const bool fuseOk = HbExportTryFuseOnce(
					current,
					ndShape,
					fuzzy,
					healInputs,
					candidate,
					failReason);

				if (!fuseOk) {
					attemptRejects << QString::fromUtf8(u8"attempt %1 fuzzy=%2 heal=%3：%4")
						.arg(static_cast<int>(ai))
						.arg(fuzzy)
						.arg(healInputs)
						.arg(failReason);
					continue;
				}

				QString rejectReason;
				const bool accept = HbExportAcceptTrueUnionResult(
					candidate,
					current,
					ndShape,
					currentVolumeBefore,
					ndVolume,
					rejectReason);

				qDebug() << "[导出融合] Fuse 尝试 index =" << static_cast<int>(i)
					<< "attempt =" << static_cast<int>(ai)
					<< "fuzzy =" << fuzzy
					<< "heal =" << healInputs
					<< "accepted =" << accept
					<< rejectReason;

				if (accept) {
					accepted = true;
					acceptedShape = candidate;
					break;
				}

				attemptRejects << QString::fromUtf8(u8"attempt %1 fuzzy=%2 heal=%3：%4")
					.arg(static_cast<int>(ai))
					.arg(fuzzy)
					.arg(healInputs)
					.arg(rejectReason);
			}

			if (!accepted || acceptedShape.IsNull()) {
				usedCompoundFallback = true;
				failureReports << QString::fromUtf8(u8"index=%1：所有 BooleanFuse 尝试均未得到可靠真并集。\n%2")
					.arg(static_cast<int>(i))
					.arg(attemptRejects.join(QString::fromUtf8(u8"\n")));

				qWarning() << "[导出融合] 严格融合失败，停止导出，non-design index =" << static_cast<int>(i);
				break;
			}

			current = acceptedShape;
			++successCount;
			usedFuse = true;
			qDebug() << "[导出融合] Fuse 接受，non-design index =" << static_cast<int>(i);
		}

		if (usedCompoundFallback) {
			report = QString::fromUtf8(
				u8"BooleanFuse 严格融合失败，已停止本次融合导出。\n\n"
				u8"这次不会再自动退回 Compound 装配体，因为你需要的是单一融合实体。\n"
				u8"已成功融合数量：%1 / %2。\n\n"
				u8"失败原因：\n%3\n\n"
				u8"含义：OCC 返回的 Fuse 结果没有同时满足“结果有效、包围盒保留 non-design、总体积不异常变小”。\n"
				u8"这说明当前接口处的布尔结果仍然会吞掉 non-design 外壳，不能作为最终实体导出。")
				.arg(successCount)
				.arg(static_cast<int>(nonDesignSolids.size()))
				.arg(failureReports.join(QString::fromUtf8(u8"\n\n")));

			return TopoDS_Shape();
		}

		report = QString::fromUtf8(u8"BooleanFuse 成功：已严格融合 %1 个 non-design 实体，未使用 Compound 回退。")
			.arg(successCount);

		return current.IsNull() ? rebuiltShape : current;
	}

}

void CageInteractorStyle::outfile()
{
	if (!m_refiner || m_faceVerts.empty()) {
		qDebug() << "[导出] m_refiner 为空或 m_faceVerts 为空，无法导出。";
		return;
	}

	QString savePath = QFileDialog::getSaveFileName(
		nullptr,
		"导出实体",
		"C:/",
		"STEP Files (*.stp *.step)"
	);

	if (savePath.isEmpty()) {
		return;
	}

	qDebug() << ">>> 开始执行全局连续放样导出...";
	qDebug() << ">>> [正式导出+SEH防崩定位版] 保留原有 MakeFilling 自动补洞；仅在 OCC 访问冲突时显示当前 wire 并停止。";

	try {
		double safeTol = 1e-3;
		const double faceTol = 1e-6;
		BRepBuilderAPI_Sewing sewer(safeTol);

		// ============================================================
		// 调试显示：只在 MakeFilling 真的触发 SEH 崩溃时使用
		// ============================================================

		vtkSmartPointer<vtkPoints> debugWirePts = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkCellArray> debugWireLines = vtkSmartPointer<vtkCellArray>::New();

		auto AppendSingleWireToDebugPolyData = [&](
			const TopoDS_Wire& wire,
			int samplesPerEdge = 8)
		{
			debugWirePts->Reset();
			debugWireLines->Reset();

			TopExp_Explorer edgeExp(wire, TopAbs_EDGE);

			while (edgeExp.More()) {
				TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

				Standard_Real first = 0.0;
				Standard_Real last = 0.0;

				Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);

				if (!curve.IsNull()) {
					vtkIdType prevId = -1;

					for (int k = 0; k <= samplesPerEdge; ++k) {
						double t = static_cast<double>(k) / static_cast<double>(samplesPerEdge);
						Standard_Real u = first + (last - first) * t;

						gp_Pnt p = curve->Value(u);

						vtkIdType id = debugWirePts->InsertNextPoint(
							p.X(),
							p.Y(),
							p.Z()
						);

						if (prevId >= 0) {
							vtkNew<vtkLine> line;
							line->GetPointIds()->SetId(0, prevId);
							line->GetPointIds()->SetId(1, id);
							debugWireLines->InsertNextCell(line);
						}

						prevId = id;
					}
				}

				edgeExp.Next();
			}
		};

		auto ShowDebugWiresNow = [&]()
		{
			if (debugWirePts->GetNumberOfPoints() <= 0) {
				if (false) qDebug() << "[ExportDebug] no debug wire points to display.";
				return;
			}

			double bounds[6];
			debugWirePts->GetBounds(bounds);

			double center[3] = {
				0.5 * (bounds[0] + bounds[1]),
				0.5 * (bounds[2] + bounds[3]),
				0.5 * (bounds[4] + bounds[5])
			};

			double dx = bounds[1] - bounds[0];
			double dy = bounds[3] - bounds[2];
			double dz = bounds[5] - bounds[4];
			double diag = std::sqrt(dx * dx + dy * dy + dz * dz);

			if (false) qDebug() << "[ExportDebug] crash wire display bounds:"
				<< "x =" << bounds[0] << bounds[1]
				<< "y =" << bounds[2] << bounds[3]
				<< "z =" << bounds[4] << bounds[5]
				<< "center =" << center[0] << center[1] << center[2]
				<< "diag =" << diag;

			vtkSmartPointer<vtkPolyData> debugWirePoly =
				vtkSmartPointer<vtkPolyData>::New();

			debugWirePoly->SetPoints(debugWirePts);
			debugWirePoly->SetLines(debugWireLines);

			vtkNew<vtkTubeFilter> tube;
			tube->SetInputData(debugWirePoly);
			tube->SetRadius(0.25);
			tube->SetNumberOfSides(12);
			tube->CappingOn();
			tube->Update();

			vtkNew<vtkSphereSource> sphere;
			sphere->SetCenter(center);
			sphere->SetRadius(std::max(1.0, diag * 0.08));
			sphere->SetThetaResolution(24);
			sphere->SetPhiResolution(12);
			sphere->Update();

			vtkNew<vtkAppendPolyData> appendDebug;
			appendDebug->AddInputConnection(tube->GetOutputPort());
			appendDebug->AddInputConnection(sphere->GetOutputPort());
			appendDebug->Update();

			vtkNew<vtkPolyDataMapper> mapper;
			mapper->SetInputConnection(appendDebug->GetOutputPort());

			vtkSmartPointer<vtkActor> actor =
				vtkSmartPointer<vtkActor>::New();

			actor->SetMapper(mapper);
			actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
			actor->GetProperty()->SetOpacity(1.0);
			actor->GetProperty()->SetLighting(false);
			actor->SetPickable(0);
			actor->VisibilityOn();

			vtkRenderer* ren = this->GetCurrentRenderer();

			if (!ren && this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
				vtkRendererCollection* rens =
					this->GetInteractor()->GetRenderWindow()->GetRenderers();

				if (rens) {
					rens->InitTraversal();
					ren = rens->GetNextItem();
				}
			}

			if (!ren) {
				if (false) qDebug() << "[ExportDebug] cannot display crash wire: renderer is null.";
				return;
			}

			static vtkSmartPointer<vtkActor> lastCrashWireActor = nullptr;

			if (lastCrashWireActor) {
				ren->RemoveActor(lastCrashWireActor);
			}

			lastCrashWireActor = actor;
			ren->AddActor(lastCrashWireActor);
			ren->ResetCameraClippingRange();

			if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
				this->GetInteractor()->GetRenderWindow()->Render();
			}

			if (false) qDebug() << "[ExportDebug] OCC crash wire displayed in red with center sphere."
				<< "points =" << debugWirePts->GetNumberOfPoints()
				<< "lines =" << debugWireLines->GetNumberOfCells();
		};

		auto PrintWireInfo = [&](const TopoDS_Wire& wire, int wireIndex, int wireEdgeCount)
		{
			Bnd_Box wireBox;
			BRepBndLib::Add(wire, wireBox);

			if (!wireBox.IsVoid()) {
				double xmin, ymin, zmin, xmax, ymax, zmax;
				wireBox.Get(xmin, ymin, zmin, xmax, ymax, zmax);

				qDebug() << "[ExportDebug] wire" << wireIndex
					<< "edge count =" << wireEdgeCount
					<< "bounds x=" << xmin << xmax
					<< "y=" << ymin << ymax
					<< "z=" << zmin << zmax;
			}
			else {
				qDebug() << "[ExportDebug] wire" << wireIndex
					<< "edge count =" << wireEdgeCount
					<< "bounds = void";
			}

			bool wireClosed = false;

			try {
				wireClosed = BRep_Tool::IsClosed(wire);
			}
			catch (...) {
				wireClosed = false;
			}

			qDebug() << "[ExportDebug] wire" << wireIndex
				<< "closed =" << wireClosed;

			std::vector<gp_Pnt> starts;
			std::vector<gp_Pnt> ends;
			std::vector<TopoDS_Edge> orderedEdges;

			try {
				for (BRepTools_WireExplorer wx(wire); wx.More(); wx.Next()) {
					TopoDS_Edge edge = TopoDS::Edge(wx.Current());
					TopoDS_Vertex v0;
					TopoDS_Vertex v1;
					TopExp::Vertices(edge, v0, v1);

					if (v0.IsNull() || v1.IsNull()) {
						continue;
					}

					if (edge.Orientation() == TopAbs_REVERSED) {
						std::swap(v0, v1);
					}

					gp_Pnt p0 = BRep_Tool::Pnt(v0);
					gp_Pnt p1 = BRep_Tool::Pnt(v1);

					starts.push_back(p0);
					ends.push_back(p1);
					orderedEdges.push_back(edge);
				}
			}
			catch (...) {
				qWarning() << "[ExportDebug] failed to inspect ordered wire edges"
					<< "wire =" << wireIndex;
			}

			double maxGap = 0.0;
			for (size_t ei = 0; ei < starts.size(); ++ei) {
				const gp_Pnt& a = starts[ei];
				const gp_Pnt& b = ends[ei];
				const double edgeLen = a.Distance(b);

				double gapToPrev = 0.0;
				if (!ends.empty()) {
					const gp_Pnt& prev = ends[(ei + ends.size() - 1) % ends.size()];
					gapToPrev = prev.Distance(a);
					maxGap = std::max(maxGap, gapToPrev);
				}

				qDebug() << "[ExportDebug] wire edge"
					<< "wire =" << wireIndex
					<< "edge =" << static_cast<int>(ei)
					<< "len =" << edgeLen
					<< "gapFromPrev =" << gapToPrev
					<< "start =" << a.X() << a.Y() << a.Z()
					<< "end =" << b.X() << b.Y() << b.Z();
			}

			if (!starts.empty()) {
				qWarning() << "[ExportDebug] wire closure detail"
					<< "wire =" << wireIndex
					<< "ordered edges =" << static_cast<int>(starts.size())
					<< "max gap =" << maxGap;
			}

			for (size_t a = 0; a < orderedEdges.size(); ++a) {
				for (size_t b = a + 1; b < orderedEdges.size(); ++b) {
					const bool oppositeEndpoints =
						starts[a].Distance(ends[b]) <= 1e-6 &&
						ends[a].Distance(starts[b]) <= 1e-6;

					if (!oppositeEndpoints) {
						continue;
					}

					Standard_Real aFirst = 0.0;
					Standard_Real aLast = 0.0;
					Standard_Real bFirst = 0.0;
					Standard_Real bLast = 0.0;
					Handle(Geom_Curve) ca = BRep_Tool::Curve(orderedEdges[a], aFirst, aLast);
					Handle(Geom_Curve) cb = BRep_Tool::Curve(orderedEdges[b], bFirst, bLast);

					if (ca.IsNull() || cb.IsNull()) {
						continue;
					}

					gp_Pnt aMid = ca->Value(0.5 * (aFirst + aLast));
					gp_Pnt bMid = cb->Value(0.5 * (bFirst + bLast));
					const double midGap = aMid.Distance(bMid);

					qWarning() << "[ExportDebug] opposite duplicate free edges"
						<< "wire =" << wireIndex
						<< "edgeA =" << static_cast<int>(a)
						<< "edgeB =" << static_cast<int>(b)
						<< "midGap =" << midGap;
				}
			}
		};

		auto ShowCrashWireAndStop = [&](const TopoDS_Wire& wire,
			int wireIndex,
			int wireEdgeCount,
			const QString& reason)
		{
			AppendSingleWireToDebugPolyData(wire, 8);
			ShowDebugWiresNow();

			qWarning() << "[ExportDebug] STOP after intercepted OCC crash."
				<< "reason =" << reason
				<< "wire index =" << wireIndex
				<< "edge count =" << wireEdgeCount
				<< "This wire is displayed in red.";
		};

		// ============================================================
		// 1. 四边 coarse face 放样
		// ============================================================

		auto AppendControlFaceToDebug = [&](const std::vector<int>& face)
		{
			if (face.size() < 2) return;
			std::vector<vtkIdType> ids;
			ids.reserve(face.size());
			for (int vid : face) {
				if (vid < 0 || vid >= static_cast<int>(m_controlVerts.size())) continue;
				const OsdVertex& v = m_controlVerts[vid];
				ids.push_back(debugWirePts->InsertNextPoint(v.x, v.y, v.z));
			}
			if (ids.size() < 2) return;
			for (size_t i = 0; i < ids.size(); ++i) {
				vtkIdType line[2] = { ids[i], ids[(i + 1) % ids.size()] };
				debugWireLines->InsertNextCell(2, line);
			}
		};

		auto AppendControlEdgeToDebug = [&](int a, int b)
		{
			if (a < 0 || b < 0 ||
				a >= static_cast<int>(m_controlVerts.size()) ||
				b >= static_cast<int>(m_controlVerts.size()))
			{
				return;
			}
			const OsdVertex& va = m_controlVerts[a];
			const OsdVertex& vb = m_controlVerts[b];
			vtkIdType p0 = debugWirePts->InsertNextPoint(va.x, va.y, va.z);
			vtkIdType p1 = debugWirePts->InsertNextPoint(vb.x, vb.y, vb.z);
			vtkIdType line[2] = { p0, p1 };
			debugWireLines->InsertNextCell(2, line);
		};

		auto ValidateControlMeshBeforeExport = [&]() -> bool
		{
			debugWirePts->Reset();
			debugWireLines->Reset();

			if (m_controlVerts.empty() || m_faceVerts.empty() || m_vertsPerFace.empty()) {
				qWarning() << "[导出保护] 控制网格为空，停止导出。";
				return false;
			}

			double bounds[6] = {
				m_controlVerts[0].x, m_controlVerts[0].x,
				m_controlVerts[0].y, m_controlVerts[0].y,
				m_controlVerts[0].z, m_controlVerts[0].z
			};
			for (const OsdVertex& v : m_controlVerts) {
				bounds[0] = std::min(bounds[0], static_cast<double>(v.x));
				bounds[1] = std::max(bounds[1], static_cast<double>(v.x));
				bounds[2] = std::min(bounds[2], static_cast<double>(v.y));
				bounds[3] = std::max(bounds[3], static_cast<double>(v.y));
				bounds[4] = std::min(bounds[4], static_cast<double>(v.z));
				bounds[5] = std::max(bounds[5], static_cast<double>(v.z));
			}
			const double dx = bounds[1] - bounds[0];
			const double dy = bounds[3] - bounds[2];
			const double dz = bounds[5] - bounds[4];
			const double diag = std::max(1e-6, std::sqrt(dx * dx + dy * dy + dz * dz));
			const double pointTol = std::max(1e-7, diag * 1e-6);
			const double edgeTol2 = pointTol * pointTol;
			const double areaTol = std::max(1e-12, diag * diag * 1e-10);

			auto getPoint = [&](int vid, double p[3]) -> bool
			{
				if (vid < 0 || vid >= static_cast<int>(m_controlVerts.size())) return false;
				const OsdVertex& v = m_controlVerts[vid];
				if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) return false;
				p[0] = v.x;
				p[1] = v.y;
				p[2] = v.z;
				return true;
			};

			auto dist2PointSegment = [&](const double p[3], const double a[3], const double b[3], double& t) -> double
			{
				double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
				double ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
				const double len2 = vtkMath::Dot(ab, ab);
				if (len2 <= 1e-20) {
					t = 0.0;
					return vtkMath::Distance2BetweenPoints(p, a);
				}
				t = vtkMath::Dot(ap, ab) / len2;
				t = std::max(0.0, std::min(1.0, t));
				double q[3] = { a[0] + ab[0] * t, a[1] + ab[1] * t, a[2] + ab[2] * t };
				return vtkMath::Distance2BetweenPoints(p, q);
			};

			int nonQuadCount = 0;
			int invalidFaceCount = 0;
			int degenerateQuadCount = 0;
			int nonManifoldEdgeCount = 0;
			int tJunctionCount = 0;
			size_t cursor = 0;
			std::map<std::pair<int, int>, int> edgeUseCount;

			for (size_t f = 0; f < m_vertsPerFace.size(); ++f) {
				const int n = m_vertsPerFace[f];
				if (n <= 0 || cursor + static_cast<size_t>(n) > m_faceVerts.size()) {
					++invalidFaceCount;
					break;
				}

				std::vector<int> face;
				face.reserve(n);
				for (int i = 0; i < n; ++i) {
					face.push_back(m_faceVerts[cursor + i]);
				}
				cursor += n;

				std::set<int> uniqueIds(face.begin(), face.end());
				if (uniqueIds.size() < 3) {
					++invalidFaceCount;
					AppendControlFaceToDebug(face);
					continue;
				}

				std::vector<std::array<double, 3>> pts(n);
				bool validPoint = true;
				for (int i = 0; i < n; ++i) {
					validPoint = validPoint && getPoint(face[i], pts[i].data());
				}
				if (!validPoint) {
					++invalidFaceCount;
					AppendControlFaceToDebug(face);
					continue;
				}

				bool tinyEdge = false;
				for (int i = 0; i < n; ++i) {
					if (vtkMath::Distance2BetweenPoints(pts[i].data(), pts[(i + 1) % n].data()) <= edgeTol2) {
						tinyEdge = true;
						break;
					}
				}

				double area2 = 0.0;
				for (int i = 1; i + 1 < n; ++i) {
					double e0[3] = {
						pts[i][0] - pts[0][0],
						pts[i][1] - pts[0][1],
						pts[i][2] - pts[0][2]
					};
					double e1[3] = {
						pts[i + 1][0] - pts[0][0],
						pts[i + 1][1] - pts[0][1],
						pts[i + 1][2] - pts[0][2]
					};
					double cross[3];
					vtkMath::Cross(e0, e1, cross);
					area2 += vtkMath::Norm(cross);
				}

				if (tinyEdge || area2 <= areaTol) {
					++degenerateQuadCount;
					AppendControlFaceToDebug(face);
					continue;
				}

				if (n != 4) {
					++nonQuadCount;
					continue;
				}

				for (int i = 0; i < 4; ++i) {
					int a = face[i];
					int b = face[(i + 1) % 4];
					if (a > b) std::swap(a, b);
					edgeUseCount[{a, b}]++;
				}
			}

			for (const auto& kv : edgeUseCount) {
				if (kv.second > 2) {
					++nonManifoldEdgeCount;
					AppendControlEdgeToDebug(kv.first.first, kv.first.second);
				}
			}

			for (const auto& kv : edgeUseCount) {
				const int aId = kv.first.first;
				const int bId = kv.first.second;
				double a[3], b[3];
				if (!getPoint(aId, a) || !getPoint(bId, b)) continue;
				for (int vid = 0; vid < static_cast<int>(m_controlVerts.size()); ++vid) {
					if (vid == aId || vid == bId) continue;
					double p[3];
					if (!getPoint(vid, p)) continue;
					double t = 0.0;
					if (dist2PointSegment(p, a, b, t) <= edgeTol2 &&
						t > 1e-4 &&
						t < 1.0 - 1e-4)
					{
						++tJunctionCount;
						AppendControlEdgeToDebug(aId, bId);
						break;
					}
				}
			}

			const bool ok =
				invalidFaceCount == 0 &&
				degenerateQuadCount == 0 &&
				nonManifoldEdgeCount == 0 &&
				tJunctionCount == 0;

			if (!ok) {
				ShowDebugWiresNow();
				{
					QString msg = QString::fromUtf8(
						u8"导出前拓扑保护已拦截。\n\n"
						u8"检测到真正危险的控制网格拓扑：\n"
						u8"非法面=%1，退化面=%2，非流形边=%3，T-junction=%4。\n"
						u8"允许跳过的非四边面=%5。\n\n"
						u8"红色高亮处是需要修正的位置。正常三角 hub 面不会被拦截；但开孔边界不能被额外切碎或形成 T-junction。")
						.arg(invalidFaceCount)
						.arg(degenerateQuadCount)
						.arg(nonManifoldEdgeCount)
						.arg(tJunctionCount)
						.arg(nonQuadCount);
					QMessageBox::critical(nullptr, QString::fromUtf8(u8"导出被保护拦截"), msg);
					qWarning() << "[ExportGuard]" << msg;
				}
				return ok;
				QString msg = QString::fromUtf8(
					u8"导出前拓扑保护已拦截。\n\n"
					u8"当前控制网格存在会导致 OCC 导出/缝合崩溃的开孔或接头拓扑：\n"
					u8"非四边面=%1，非法面=%2，退化四边面=%3，非流形边=%4，T-junction=%5。\n\n"
					u8"红色高亮处就是需要修正的位置。开孔边界必须保持四点四边，并且不能被额外切碎。")
					.arg(nonQuadCount)
					.arg(invalidFaceCount)
					.arg(degenerateQuadCount)
					.arg(nonManifoldEdgeCount)
					.arg(tJunctionCount);
				QMessageBox::critical(nullptr, QString::fromUtf8(u8"导出被保护拦截"), msg);
				qWarning() << "[导出保护]" << msg;
			}

			return ok;
		};

		if (!ValidateControlMeshBeforeExport()) {
			return;
		}

		int maxLevel = m_refiner->GetMaxLevel();
		int numCoarseFaces = m_refiner->GetLevel(0).GetNumFaces();
		int gridSteps = 1 << maxLevel;
		int ptsGridSize = gridSteps + 1;

		if (false) qDebug() << "[ExportDebug] maxLevel =" << maxLevel
			<< "numCoarseFaces =" << numCoarseFaces
			<< "gridSteps =" << gridSteps;

		int coarseTriCount = 0;
		int coarseQuadCount = 0;
		int coarseNgonCount = 0;
		int coarseOtherCount = 0;

		for (int cf = 0; cf < numCoarseFaces; ++cf) {
			Far::ConstIndexArray cv =
				m_refiner->GetLevel(0).GetFaceVertices(cf);

			if (cv.size() == 3) {
				coarseTriCount++;
			}
			else if (cv.size() == 4) {
				coarseQuadCount++;
			}
			else if (cv.size() > 4) {
				coarseNgonCount++;
			}
			else {
				coarseOtherCount++;
			}
		}

		if (false) qDebug() << "[ExportDebug] coarse faces:"
			<< "total =" << numCoarseFaces
			<< "tri =" << coarseTriCount
			<< "quad =" << coarseQuadCount
			<< "ngon =" << coarseNgonCount
			<< "other =" << coarseOtherCount;

		int lastLevelOffset = 0;

		for (int i = 0; i < maxLevel; ++i) {
			lastLevelOffset += m_refiner->GetLevel(i).GetNumVertices();
		}

		int exportedFaceCount = 0;
		int skippedNonQuadCount = 0;
		int skippedInvalidGridCount = 0;
		int failedFitCount = 0;
		int failedMakeFaceCount = 0;

		for (int coarseFace = 0; coarseFace < numCoarseFaces; ++coarseFace) {
			Far::ConstIndexArray coarseVerts =
				m_refiner->GetLevel(0).GetFaceVertices(coarseFace);

			if (coarseVerts.size() != 4) {
				skippedNonQuadCount++;
				continue;
			}

			std::vector<std::vector<int>> faceGrid(
				gridSteps,
				std::vector<int>(gridSteps, -1)
			);

			std::function<void(int, int, int, int, int)> fillFaceGrid;

			fillFaceGrid = [&](int level, int f, int r, int c, int size)
			{
				if (level == maxLevel) {
					faceGrid[r][c] = f;
					return;
				}

				int half = size / 2;

				Far::ConstIndexArray children =
					m_refiner->GetLevel(level).GetFaceChildFaces(f);

				if (children.size() == 4) {
					fillFaceGrid(level + 1, children[0], r, c, half);
					fillFaceGrid(level + 1, children[1], r, c + half, half);
					fillFaceGrid(level + 1, children[2], r + half, c + half, half);
					fillFaceGrid(level + 1, children[3], r + half, c, half);
				}
			};

			fillFaceGrid(0, coarseFace, 0, 0, gridSteps);

			TColgp_Array2OfPnt poles(1, ptsGridSize, 1, ptsGridSize);

			bool isValidGrid = true;

			for (int r = 0; r < gridSteps; ++r) {
				for (int c = 0; c < gridSteps; ++c) {
					int microFace = faceGrid[r][c];

					if (microFace == -1) {
						isValidGrid = false;
						break;
					}

					Far::ConstIndexArray fverts =
						m_refiner->GetLevel(maxLevel).GetFaceVertices(microFace);

					if (fverts.size() != 4) {
						isValidGrid = false;
						break;
					}

					auto getPnt = [&](int vIdx) -> gp_Pnt
					{
						OsdVertex ov =
							m_refinedVerts[lastLevelOffset + vIdx];

						return gp_Pnt(ov.x, ov.y, ov.z);
					};

					poles.SetValue(r + 1, c + 1, getPnt(fverts[0]));
					poles.SetValue(r + 1, c + 2, getPnt(fverts[1]));
					poles.SetValue(r + 2, c + 2, getPnt(fverts[2]));
					poles.SetValue(r + 2, c + 1, getPnt(fverts[3]));
				}

				if (!isValidGrid) {
					break;
				}
			}

			if (!isValidGrid) {
				skippedInvalidGridCount++;

				if (false) qDebug() << "[ExportDebug] skip invalid subdiv grid coarseFace ="
					<< coarseFace;

				continue;
			}

			try {
				Handle(Geom_BSplineSurface) surface =
					HbBuildClampedBSplineSurface(poles, ptsGridSize, ptsGridSize);

				if (surface.IsNull()) {
					failedFitCount++;

					qWarning() << "[ExportDebug] BSpline construction failed coarseFace ="
						<< coarseFace;

					continue;
				}

				BRepBuilderAPI_MakeFace mkFace(surface, faceTol);

				if (mkFace.IsDone() && !mkFace.Face().IsNull()) {
					sewer.Add(mkFace.Face());
					exportedFaceCount++;
				}
				else {
					failedMakeFaceCount++;

					if (false) qDebug() << "[ExportDebug] MakeFace from BSpline failed coarseFace ="
						<< coarseFace;
				}
			}
			catch (const Standard_Failure& e) {
				failedFitCount++;

				if (false) qDebug() << "[ExportDebug] OCC exception while exporting coarseFace ="
					<< coarseFace
					<< e.GetMessageString();
			}
			catch (...) {
				failedFitCount++;

				if (false) qDebug() << "[ExportDebug] unknown exception while exporting coarseFace ="
					<< coarseFace;
			}
		}

		if (false) qDebug() << "[ExportDebug] surface export summary:"
			<< "exported =" << exportedFaceCount
			<< "skipped non-quad =" << skippedNonQuadCount
			<< "skipped invalid grid =" << skippedInvalidGridCount
			<< "failed fit =" << failedFitCount
			<< "failed makeface =" << failedMakeFaceCount;

		if (exportedFaceCount <= 0) {
			qDebug() << "[导出] 没有成功生成任何曲面，导出中止。";
			return;
		}

		// ============================================================
		// 2. 第一次 sewing + ShapeFix
		// ============================================================

		qWarning() << "[ExportDebug] entering first sewer.Perform."
			<< "exported faces =" << exportedFaceCount
			<< "skipped non-quad =" << skippedNonQuadCount
			<< "failed fit =" << failedFitCount;

		try {
			sewer.Perform();
		}
		catch (const Standard_Failure& e) {
			qWarning() << "[ExportDebug] OCC exception in first sewer.Perform."
				<< "exported faces =" << exportedFaceCount
				<< "skipped non-quad =" << skippedNonQuadCount
				<< "failed fit =" << failedFitCount
				<< e.GetMessageString();
			return;
		}
		catch (...) {
			qWarning() << "[ExportDebug] unknown C++ exception in first sewer.Perform."
				<< "exported faces =" << exportedFaceCount
				<< "skipped non-quad =" << skippedNonQuadCount
				<< "failed fit =" << failedFitCount;
			return;
		}

		TopoDS_Shape sewed = sewer.SewedShape();

		if (sewed.IsNull()) {
			qDebug() << "[导出] 第一次 sewing 结果为空。";
			return;
		}

		ShapeFix_Shape healer(sewed);
		healer.Perform();

		TopoDS_Shape healedShape = healer.Shape();

		if (healedShape.IsNull()) {
			qDebug() << "[导出] ShapeFix 后结果为空。";
			return;
		}

		// ============================================================
		// 3. 自由边补洞
		// ============================================================

		BRepBuilderAPI_Sewing finalSewer(safeTol);
		finalSewer.Add(healedShape);

		TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;

		TopExp::MapShapesAndAncestors(
			healedShape,
			TopAbs_EDGE,
			TopAbs_FACE,
			edgeFaceMap
		);

		Handle(TopTools_HSequenceOfShape) freeEdges =
			new TopTools_HSequenceOfShape();

		for (int i = 1; i <= edgeFaceMap.Extent(); i++) {
			if (edgeFaceMap(i).Extent() == 1) {
				freeEdges->Append(TopoDS::Edge(edgeFaceMap.FindKey(i)));
			}
		}

		if (false) qDebug() << "[ExportDebug] free edge count after first sewing ="
			<< freeEdges->Length();

		int connectedWireCount = 0;
		int makeFaceCapCount = 0;
		int fillingCapCount = 0;
		int failedFillingCount = 0;
		int sehCrashCount = 0;

		if (freeEdges->Length() > 0)
		{
			Handle(TopTools_HSequenceOfShape) closedWires =
				new TopTools_HSequenceOfShape();

			ShapeAnalysis_FreeBounds::ConnectEdgesToWires(
				freeEdges,
				safeTol,
				Standard_False,
				closedWires
			);

			connectedWireCount = closedWires->Length();

			if (false) qDebug() << "[ExportDebug] connected free wire count ="
				<< connectedWireCount;

			for (int i = 1; i <= closedWires->Length(); i++)
			{
				TopoDS_Wire wire =
					TopoDS::Wire(closedWires->Value(i));

				int wireEdgeCount = 0;

				TopExp_Explorer countExp(wire, TopAbs_EDGE);

				while (countExp.More()) {
					wireEdgeCount++;
					countExp.Next();
				}

				PrintWireInfo(wire, i, wireEdgeCount);

				bool capOk = false;
				TopoDS_Shape capShape;

				try {
					BRepBuilderAPI_MakeFace capMaker(wire);

					if (capMaker.IsDone() && !capMaker.Face().IsNull())
					{
						capOk = true;
						capShape = capMaker.Face();
						makeFaceCapCount++;

						if (false) qDebug() << "[ExportDebug] MakeFace cap success on wire" << i;
					}
					else {
						if (false) qDebug() << "[ExportDebug] MakeFace cap not done on wire" << i;
					}
				}
				catch (const Standard_Failure& e) {
					if (false) qDebug() << "[ExportDebug] MakeFace cap exception on wire" << i
						<< e.GetMessageString();
				}
				catch (...) {
					if (false) qDebug() << "[ExportDebug] MakeFace cap unknown exception on wire" << i;
				}

				if (capOk && !capShape.IsNull())
				{
					finalSewer.Add(capShape);
					continue;
				}

				if (wireEdgeCount > 4) {
					failedFillingCount++;
					qWarning() << "[ExportDebug] skip MakeFilling for complex free wire."
						<< "wire index =" << i
						<< "edge count =" << wireEdgeCount
						<< "This is a generation/sewing gap, not a normal cap hole.";

					ShowCrashWireAndStop(
						wire,
						i,
						wireEdgeCount,
						QString("Complex free wire with %1 edges; skipped unsafe MakeFilling")
						.arg(wireEdgeCount)
					);

					return;
				}

				// ====================================================
				// 关键：不再按 wireEdgeCount 提前停。
				// 所有 wire 都恢复原来的 MakeFilling 自动补洞。
				// 只有 filler.Build() 真的 SEH 崩溃时，才显示当前 wire 并停止。
				// ====================================================

				if (false) qDebug() << ">>> MakeFace 失败，启用 MakeFilling 自动补洞。"
					<< "wire index =" << i
					<< "edge count =" << wireEdgeCount;

				try {
					BRepOffsetAPI_MakeFilling filler;

					TopExp_Explorer edgeExp(wire, TopAbs_EDGE);

					int addedEdgeCount = 0;

					while (edgeExp.More())
					{
						DWORD addSehCode = 0;
						const TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
						if (!SafeAddFillingEdgeWithSEH(&filler, &edge, addSehCode)) {
							sehCrashCount++;

							const QString codeText = HbSehCodeText(addSehCode);
							qWarning() << "[ExportDebug] OCC SEH in MakeFilling.Add edge."
								<< "wire index =" << i
								<< "wire edge count =" << wireEdgeCount
								<< "added edge index =" << addedEdgeCount
								<< "SEH code =" << codeText;

							ShowCrashWireAndStop(
								wire,
								i,
								wireEdgeCount,
								QString("MakeFilling.Add edge SEH crash, code = %1").arg(codeText)
							);

							return;
						}
						addedEdgeCount++;
						edgeExp.Next();
					}

					qDebug() << "[ExportDebug] filling added edge count ="
						<< addedEdgeCount
						<< "wire index =" << i;

					DWORD sehCode = 0;
					bool buildReturned = SafeBuildFillingWithSEH(filler, sehCode);

					if (!buildReturned) {
						sehCrashCount++;

						QString codeText = HbSehCodeText(sehCode);

						qWarning() << "[ExportDebug] MakeFilling SEH crash intercepted."
							<< "wire index =" << i
							<< "edge count =" << wireEdgeCount
							<< "SEH code =" << codeText;

						ShowCrashWireAndStop(
							wire,
							i,
							wireEdgeCount,
							QString("MakeFilling SEH crash intercepted, code = %1").arg(codeText)
						);

						return;
					}

					qDebug() << "[ExportDebug] reading MakeFilling result"
						<< "wire index =" << i
						<< "edge count =" << wireEdgeCount;

					bool fillingDone = filler.IsDone();
					TopoDS_Shape fillingShape;
					if (fillingDone) {
						fillingShape = filler.Shape();
					}

					qDebug() << "[ExportDebug] filling result:"
						<< "wire index =" << i
						<< "IsDone =" << fillingDone
						<< "ShapeNull =" << fillingShape.IsNull();

					if (fillingDone && !fillingShape.IsNull())
					{
						finalSewer.Add(fillingShape);
						fillingCapCount++;

						qDebug() << "[ExportDebug] MakeFilling cap success on wire" << i;
					}
					else {
						failedFillingCount++;

						qWarning() << "[ExportDebug] MakeFilling cap failed on wire" << i;
					}
				}
				catch (const Standard_Failure& e) {
					failedFillingCount++;

					qWarning() << "[ExportDebug] MakeFilling Standard_Failure on wire" << i
						<< e.GetMessageString();

					ShowCrashWireAndStop(
						wire,
						i,
						wireEdgeCount,
						QString("MakeFilling Standard_Failure: %1")
						.arg(QString::fromUtf8(e.GetMessageString()))
					);

					return;
				}
				catch (...) {
					failedFillingCount++;

					qWarning() << "[ExportDebug] MakeFilling C++ unknown exception on wire" << i;

					ShowCrashWireAndStop(
						wire,
						i,
						wireEdgeCount,
						QString("MakeFilling unknown C++ exception")
					);

					return;
				}
			}
		}

		qDebug() << "[ExportDebug] cap summary:"
			<< "connected wires =" << connectedWireCount
			<< "MakeFace caps =" << makeFaceCapCount
			<< "MakeFilling caps =" << fillingCapCount
			<< "failed fillings =" << failedFillingCount
			<< "SEH crashes =" << sehCrashCount;

		// ============================================================
		// 4. final sewing
		// ============================================================

		qWarning() << "[ExportDebug] entering finalSewer.Perform."
			<< "connected wires =" << connectedWireCount
			<< "MakeFace caps =" << makeFaceCapCount
			<< "MakeFilling caps =" << fillingCapCount
			<< "failed fillings =" << failedFillingCount;

		try {
			finalSewer.Perform();
		}
		catch (const Standard_Failure& e) {
			qWarning() << "[ExportDebug] OCC exception in finalSewer.Perform."
				<< "connected wires =" << connectedWireCount
				<< "MakeFace caps =" << makeFaceCapCount
				<< "MakeFilling caps =" << fillingCapCount
				<< "failed fillings =" << failedFillingCount
				<< e.GetMessageString();
			ShowDebugWiresNow();
			return;
		}
		catch (...) {
			qWarning() << "[ExportDebug] unknown C++ exception in finalSewer.Perform."
				<< "connected wires =" << connectedWireCount
				<< "MakeFace caps =" << makeFaceCapCount
				<< "MakeFilling caps =" << fillingCapCount
				<< "failed fillings =" << failedFillingCount;
			ShowDebugWiresNow();
			return;
		}

		TopoDS_Shape finalClosedShape = finalSewer.SewedShape();

		if (finalClosedShape.IsNull()) {
			qDebug() << "[导出] final sewing 结果为空。";
			return;
		}

		TopTools_IndexedDataMapOfShapeListOfShape finalEdgeFaceMap;

		TopExp::MapShapesAndAncestors(
			finalClosedShape,
			TopAbs_EDGE,
			TopAbs_FACE,
			finalEdgeFaceMap
		);

		int finalFreeEdgeCount = 0;
		int finalNonManifoldEdgeCount = 0;

		for (int i = 1; i <= finalEdgeFaceMap.Extent(); ++i) {
			int faceCount = finalEdgeFaceMap(i).Extent();

			if (faceCount == 1) {
				finalFreeEdgeCount++;
			}
			else if (faceCount > 2) {
				finalNonManifoldEdgeCount++;
			}
		}

		if (false) qDebug() << "[ExportDebug] after final sewing:"
			<< "free edges =" << finalFreeEdgeCount
			<< "nonmanifold edges =" << finalNonManifoldEdgeCount;

		// ============================================================
		// 5. MakeSolid
		// ============================================================

		BRepBuilderAPI_MakeSolid solidMaker;

		TopExp_Explorer shellExp(finalClosedShape, TopAbs_SHELL);

		int shellCount = 0;

		while (shellExp.More()) {
			solidMaker.Add(TopoDS::Shell(shellExp.Current()));
			shellExp.Next();
			shellCount++;
		}

		if (false) qDebug() << "[ExportDebug] shell count =" << shellCount;

		if (shellCount > 0) {
			try {
				solidMaker.Build();

				if (solidMaker.IsDone() && !solidMaker.Solid().IsNull()) {
					finalClosedShape = solidMaker.Solid();

					qDebug() << "[导出] 已生成封闭实体";
				}
				else {
					if (false) qDebug() << "[ExportDebug] MakeSolid not done, keep sewed shape.";
				}
			}
			catch (const Standard_Failure& e) {
				if (false) qDebug() << "[ExportDebug] MakeSolid exception:"
					<< e.GetMessageString()
					<< "keep sewed shape.";
			}
			catch (...) {
				if (false) qDebug() << "[ExportDebug] MakeSolid unknown exception, keep sewed shape.";
			}
		}

		// ============================================================
		// 6. 可选：与原始 STP 中的 non-design 保留实体做融合导出
		//    不破坏原 outfile：用户选择“否”时仍然只导出原重建实体。
		//    用户选择“是”时，先尝试 BooleanFuse；失败则回退 Compound 导出。
		// ============================================================

		TopoDS_Shape shapeToWrite = finalClosedShape;

		if (!m_nonDesignExportStepPath.isEmpty() &&
			(m_nonDesignExportAllStepSolids ||
			(!m_nonDesignExportElementCenters.empty() &&
				m_nonDesignExportElementCenters.size() == m_nonDesignExportIsNonDesign.size())))
		{
			QMessageBox::StandardButton fuseReply = QMessageBox::question(
				nullptr,
				QString::fromUtf8(u8"导出选项"),
				QString::fromUtf8(
					u8"检测到原始 STP non-design 保留区域数据。\n"
					u8"是否尝试将当前重建实体与 non-design 区域进行布尔融合后再导出？\n\n"
					u8"选择 Yes：尝试 BooleanFuse 布尔融合。若某个 non-design 融合结果不可靠，程序会自动以 Compound 方式保底保留，避免非设计区域丢失。\n"
					u8"选择 No：保持原 outfile 逻辑，只导出当前重建实体。\n"
					u8"选择 Cancel：取消本次导出。"),
				QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel,
				QMessageBox::Yes);

			if (fuseReply == QMessageBox::Cancel) {
				qDebug() << "[导出融合] 用户取消导出。";
				return;
			}

			if (fuseReply == QMessageBox::Yes) {
				bool usedFuse = false;
				bool usedCompoundFallback = false;
				QString fuseReport;

				shapeToWrite = HbExportTryBooleanFuseWithNonDesign(
					finalClosedShape,
					m_nonDesignExportStepPath,
					m_nonDesignExportElementCenters,
					m_nonDesignExportIsNonDesign,
					m_nonDesignExportAllStepSolids,
					usedFuse,
					usedCompoundFallback,
					fuseReport);

				qDebug() << "[导出融合]" << fuseReport;

				if (shapeToWrite.IsNull() || usedCompoundFallback) {
					QMessageBox::critical(
						nullptr,
						QString::fromUtf8(u8"non-design 融合失败"),
						fuseReport);
					qWarning() << "[导出融合] 严格融合失败，未导出 Compound 装配体。";
					return;
				}

				QMessageBox::information(
					nullptr,
					QString::fromUtf8(u8"non-design 融合结果"),
					fuseReport);
			}
		}
		else {
			qDebug() << "[导出融合] 没有可用 non-design 源数据，本次保持原 outfile 导出。";
		}

		// ============================================================
		// 7. STEP write
		// ============================================================

		STEPControl_Writer stepWriter;
		Interface_Static::SetIVal("write.step.schema", 4);

		IFSelect_ReturnStatus transferStatus =
			stepWriter.Transfer(shapeToWrite, STEPControl_AsIs);

		if (false) qDebug() << "[ExportDebug] STEP transfer status ="
			<< static_cast<int>(transferStatus);

		if (transferStatus != IFSelect_RetDone) {
			qDebug() << "[导出] STEP Transfer 失败。";
			return;
		}

		IFSelect_ReturnStatus writeStatus =
			stepWriter.Write(savePath.toStdString().c_str());

		if (false) qDebug() << "[ExportDebug] STEP write status ="
			<< static_cast<int>(writeStatus);

		if (writeStatus != IFSelect_RetDone) {
			qDebug() << "[导出] STEP Write 失败。";
			return;
		}

		qDebug() << ">>> 导出完毕，实体编辑已解锁。";

	}
	catch (const Standard_Failure& e) {
		if (false) qDebug() << ">>> [拦截] OpenCascade 异常:"
			<< e.GetMessageString();
	}
	catch (const std::exception& e) {
		if (false) qDebug() << ">>> [拦截] std::exception:"
			<< e.what();
	}
	catch (...) {
		if (false) qDebug() << ">>> [拦截] 未知错误已拦截";
	}
}

void CageInteractorStyle::OnKeyPress()
{
	vtkRenderWindowInteractor* rwi = this->GetInteractor();
	std::string key = rwi->GetKeySym();

	if (key == "Escape" || key == "Esc") {
		ExitAllEditModes();
		if (rwi) {
			rwi->Render();
		}
		return;
	}

	// 1. 
	if ((key == "z" || key == "Z") && rwi->GetControlKey()) {
		if (!m_allCages.empty()) {
			m_allCages.pop_back();
			for (auto it = m_suppressedCageSideFaces.begin(); it != m_suppressedCageSideFaces.end();) {
				if (it->cageIndex >= static_cast<int>(m_allCages.size())) {
					it = m_suppressedCageSideFaces.erase(it);
				}
				else {
					++it;
				}
			}
			RebuildGlobalMesh();
			if (m_allCages.empty() && m_PreviewBoxActor) {
				m_PreviewBoxActor->VisibilityOff();
			}
			rwi->Render();
		}
		return;
	}

	// 2. B 
	if (key == "b" || key == "B") {
		vtkRenderer* ren = rwi->FindPokedRenderer(rwi->GetEventPosition()[0], rwi->GetEventPosition()[1]);
		if (!ren) ren = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

		m_isBranchMode = !m_isBranchMode;
		if (m_isBranchMode) {
			qDebug() << "[分支模式] 进入：1 个端面生成单分支圆环，2 个端面生成双分支圆环，3 个及以上端面生成复杂分支";
			m_selectedFacesForBranch.clear();
			if (m_PathHighlightActor) m_PathHighlightActor->VisibilityOff();
		}
		else {
			qDebug() << "[分支模式] 退出";
			m_selectedFacesForBranch.clear();
			if (m_PathHighlightActor) m_PathHighlightActor->VisibilityOff();
		}

		if (ren) UpdateBranchHighlight(ren);
		rwi->Render();
		return;
	}

	if (key == "s" || key == "S") {
		vtkRenderer* ren = rwi->FindPokedRenderer(rwi->GetEventPosition()[0], rwi->GetEventPosition()[1]);
		if (!ren) ren = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

		AutoDetectAndShowTopology(ren);
		rwi->Render();
		return;
	}

	if (key == "a" || key == "A") {
		vtkRenderer* ren = rwi->FindPokedRenderer(rwi->GetEventPosition()[0], rwi->GetEventPosition()[1]);
		if (!ren) ren = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

		AutoDetectAndShowTopology(ren, true);
		rwi->Render();
		return;
	}

	if (key == "Return" || key == "KP_Enter")
	{
		if (m_isBranchMode) {
			int N = static_cast<int>(m_selectedFacesForBranch.size());
			vtkRenderer* ren = rwi->FindPokedRenderer(rwi->GetEventPosition()[0], rwi->GetEventPosition()[1]);
			if (!ren) ren = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();

			auto computeSelectedFaceFrame = [&](const SelectedFace& sel, double center[3], double normal[3]) -> bool
			{
				if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
					return false;
				}

				const CageSegment& cage = m_allCages[sel.cageIndex];
				const double(*ring)[3] = sel.isStart ? cage.pStart : cage.pEnd;
				const double(*oppositeRing)[3] = sel.isStart ? cage.pEnd : cage.pStart;

				center[0] = center[1] = center[2] = 0.0;
				double oppositeCenter[3] = { 0.0, 0.0, 0.0 };

				for (int i = 0; i < 4; ++i) {
					for (int d = 0; d < 3; ++d) {
						center[d] += ring[i][d] * 0.25;
						oppositeCenter[d] += oppositeRing[i][d] * 0.25;
					}
				}

				for (int d = 0; d < 3; ++d) {
					normal[d] = center[d] - oppositeCenter[d];
				}

				double nLen = vtkMath::Norm(normal);
				if (nLen < 1e-8) {
					double e01[3] = {
						ring[1][0] - ring[0][0],
						ring[1][1] - ring[0][1],
						ring[1][2] - ring[0][2]
					};
					double e03[3] = {
						ring[3][0] - ring[0][0],
						ring[3][1] - ring[0][1],
						ring[3][2] - ring[0][2]
					};
					vtkMath::Cross(e01, e03, normal);
					nLen = vtkMath::Norm(normal);
				}

				if (nLen < 1e-8) {
					normal[0] = 1.0;
					normal[1] = 0.0;
					normal[2] = 0.0;
				}
				else {
					vtkMath::Normalize(normal);
				}

				return true;
			};

			auto shouldUse2DCenterHubForManualBranch = [&]() -> bool
			{
				if (N == 3) {
					if (false) qDebug() << "[ManualAutoTopology] selected faces = 3, force 2D center hub topology.";
					return true;
				}

				if (N < 4) {
					if (false) qDebug() << "[ManualAutoTopology] selected faces =" << N << ", loop topology path is handled separately.";
					return false;
				}

				std::vector<std::array<double, 3>> centers;
				std::vector<std::array<double, 3>> normals;
				centers.reserve(N);
				normals.reserve(N);

				for (const SelectedFace& sel : m_selectedFacesForBranch) {
					double c[3];
					double n[3];
					if (!computeSelectedFaceFrame(sel, c, n)) {
						if (false) qDebug() << "[ManualAutoTopology] invalid selected face, fallback to 3D branch. cageIndex =" << sel.cageIndex;
						return false;
					}

					centers.push_back({ c[0], c[1], c[2] });
					normals.push_back({ n[0], n[1], n[2] });
				}

				// N >= 4 的手动自动拓扑只先判断“端面中心 + 端面法向”是否近似共面。
				// 1+2+2 不再在这里提前抢走流程；只有判定不共面后，才进入
				// GenerateNBranch3DSkeleton()，再由三维内部自己的分组逻辑处理。

				double centroid[3] = { 0.0, 0.0, 0.0 };
				for (const auto& c : centers) {
					centroid[0] += c[0] / static_cast<double>(N);
					centroid[1] += c[1] / static_cast<double>(N);
					centroid[2] += c[2] / static_cast<double>(N);
				}

				double scale = 0.0;
				for (const auto& c : centers) {
					double p[3] = { c[0], c[1], c[2] };
					scale = std::max(scale, std::sqrt(vtkMath::Distance2BetweenPoints(p, centroid)));
				}
				scale = std::max(scale, 1.0);

				double planeNormal[3] = { 0.0, 0.0, 0.0 };
				double bestCrossLen = -1.0;

				for (int i = 0; i < N; ++i) {
					double a[3] = {
						centers[i][0] - centroid[0],
						centers[i][1] - centroid[1],
						centers[i][2] - centroid[2]
					};

					for (int j = i + 1; j < N; ++j) {
						double b[3] = {
							centers[j][0] - centroid[0],
							centers[j][1] - centroid[1],
							centers[j][2] - centroid[2]
						};

						double cr[3];
						vtkMath::Cross(a, b, cr);
						double crLen = vtkMath::Norm(cr);

						if (crLen > bestCrossLen) {
							bestCrossLen = crLen;
							planeNormal[0] = cr[0];
							planeNormal[1] = cr[1];
							planeNormal[2] = cr[2];
						}
					}
				}

				if (vtkMath::Norm(planeNormal) < 1e-8) {
					for (int i = 0; i < N; ++i) {
						double ni[3] = { normals[i][0], normals[i][1], normals[i][2] };
						for (int j = i + 1; j < N; ++j) {
							double nj[3] = { normals[j][0], normals[j][1], normals[j][2] };
							double cr[3];
							vtkMath::Cross(ni, nj, cr);
							if (vtkMath::Norm(cr) > vtkMath::Norm(planeNormal)) {
								planeNormal[0] = cr[0];
								planeNormal[1] = cr[1];
								planeNormal[2] = cr[2];
							}
						}
					}
				}

				if (vtkMath::Norm(planeNormal) < 1e-8) {
					if (false) qDebug() << "[ManualAutoTopology] cannot estimate common plane, fallback to 3D branch.";
					return false;
				}

				vtkMath::Normalize(planeNormal);

				double maxCenterPlaneDist = 0.0;
				for (const auto& c : centers) {
					double v[3] = {
						c[0] - centroid[0],
						c[1] - centroid[1],
						c[2] - centroid[2]
					};
					maxCenterPlaneDist = std::max(maxCenterPlaneDist, std::abs(vtkMath::Dot(v, planeNormal)));
				}

				double maxNormalOutOfPlane = 0.0;
				double avgNormalOutOfPlane = 0.0;
				for (const auto& n : normals) {
					double nv[3] = { n[0], n[1], n[2] };
					double v = std::abs(vtkMath::Dot(nv, planeNormal));
					maxNormalOutOfPlane = std::max(maxNormalOutOfPlane, v);
					avgNormalOutOfPlane += v / static_cast<double>(N);
				}

				const double centerPlaneTol = std::max(1e-4, scale * 0.22);
				const double normalPlaneTol = 0.38;

				const bool centersNearPlane = maxCenterPlaneDist <= centerPlaneTol;
				const bool normalsNearPlane = maxNormalOutOfPlane <= normalPlaneTol;
				const bool use2D = centersNearPlane && normalsNearPlane;

				if (false) qDebug() << "[ManualAutoTopology] selected faces =" << N
					<< "planeNormal =" << planeNormal[0] << planeNormal[1] << planeNormal[2]
					<< "centerMaxDist =" << maxCenterPlaneDist
					<< "centerTol =" << centerPlaneTol
					<< "normalOutMax =" << maxNormalOutOfPlane
					<< "normalOutAvg =" << avgNormalOutOfPlane
					<< "normalTol =" << normalPlaneTol
					<< "centersNearPlane =" << centersNearPlane
					<< "normalsNearPlane =" << normalsNearPlane
					<< "result =" << (use2D ? "2D center hub" : "3D side-branch");

				return use2D;
			};

			if (N >= 3)
			{
				const bool use2D = shouldUse2DCenterHubForManualBranch();

				if (use2D) {
					if (false) qDebug() << "[ManualAutoTopology] ENTER -> GenerateNBranch(), use 2D center hub topology.";
					GenerateNBranch();
				}
				else {
					if (false) qDebug() << "[ManualAutoTopology] ENTER -> GenerateNBranch3DSkeleton(), use 3D side-branch topology.";
					GenerateNBranch3DSkeleton();
				}

				m_isBranchMode = false;
				m_selectedFacesForBranch.clear();
				UpdateBranchHighlight(ren);
				rwi->Render();
			}
			else if (N == 1 || N == 2)
			{
				bool ok = GenerateSingleBranchLoopTopology(ren);
				if (ok) {
					m_isBranchMode = false;
					m_selectedFacesForBranch.clear();
					UpdateBranchHighlight(ren);
				}
				rwi->Render();
			}
			else
			{
				if (false) qDebug() << "[ManualAutoTopology] ENTER ignored, no selected branch face.";
				return;
			}
		}
		return;
	}


	if (key == "e" || key == "E")
	{
		EnableSkeletonNodeEditMode(!m_isSkeletonEditMode);
		return;
	}


	vtkInteractorStyleTrackballCamera::OnKeyPress();
}
void CageInteractorStyle::UpdateBranchHighlight(vtkRenderer* ren)
{
	if (!m_BranchHighlightData) {
		m_BranchHighlightData = vtkSmartPointer<vtkPolyData>::New();
	}

	vtkNew<vtkPoints> pts;
	vtkNew<vtkCellArray> polys;
	vtkIdType ptIdx = 0;


	for (const auto& sel : m_selectedFacesForBranch) {
		if (sel.cageIndex >= 0 && sel.cageIndex < m_allCages.size()) {
			for (int i = 0; i < 4; i++) {
				if (sel.isStart) pts->InsertNextPoint(m_allCages[sel.cageIndex].pStart[i]);
				else pts->InsertNextPoint(m_allCages[sel.cageIndex].pEnd[i]);
			}
			polys->InsertNextCell(4);
			polys->InsertCellPoint(ptIdx);
			polys->InsertCellPoint(ptIdx + 1);
			polys->InsertCellPoint(ptIdx + 2);
			polys->InsertCellPoint(ptIdx + 3);
			ptIdx += 4;
		}
	}

	m_BranchHighlightData->SetPoints(pts);
	m_BranchHighlightData->SetPolys(polys);


	if (!m_BranchHighlightActor) {
		vtkNew<vtkPolyDataMapper> mapper;
		m_BranchHighlightActor = vtkSmartPointer<vtkActor>::New();
		m_BranchHighlightActor->SetMapper(mapper);


		m_BranchHighlightActor->GetProperty()->SetColor(0.0, 1.0, 0.0);
		m_BranchHighlightActor->GetProperty()->SetRepresentationToWireframe();
		m_BranchHighlightActor->GetProperty()->SetLineWidth(5.0);
		m_BranchHighlightActor->GetProperty()->SetLighting(false);
		m_BranchHighlightActor->SetPickable(0);

		if (ren) ren->AddActor(m_BranchHighlightActor);
	}

	vtkPolyDataMapper::SafeDownCast(m_BranchHighlightActor->GetMapper())->SetInputData(m_BranchHighlightData);


	if (m_selectedFacesForBranch.empty()) {
		m_BranchHighlightActor->VisibilityOff();
	}
	else {
		m_BranchHighlightActor->VisibilityOn();
	}

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::SetSkeletonData(
	const std::vector<double>& nodes,
	const std::vector<unsigned int>& edges,
	const std::vector<double>& radii,
	const std::vector<unsigned int>& supportOffsets,
	const std::vector<double>& supportPoints,
	vtkSmartPointer<vtkPolyData> skelPolyData)
{
	m_skelNodes = nodes;
	m_skelEdges = edges;
	m_skelRadii = radii;
	m_skelSupportOffsets = supportOffsets;
	m_skelSupportPoints = supportPoints;
	m_skelPolyData = skelPolyData;
}

int CageInteractorStyle::BuildStraightCagesFromSegments(
	const std::vector<StraightSkeletonSegmenter::Segment>& segments,
	bool replaceExisting,
	StraightCageTemplateMode templateMode)
{
	if (segments.empty()) {
		if (false) qDebug() << "[自动直线控制笼] 没有可用直线段。";
		return 0;
	}

	if (!m_RefLocator || !m_RefData) {
		if (false) qDebug() << "[自动直线控制笼] 缺少背景模型 locator，无法估计截面。";
		return 0;
	}

	const double modelDiagonal = AutoStraightReferenceDiagonal(m_RefData);

	std::vector<AutoStraightCageSegment> mergedSegments =
		AutoStraightMergeSegmentsForCage(
			segments,
			modelDiagonal,
			m_RefLocator,
			templateMode == StraightCageTemplateMode::ThinRod);

	if (mergedSegments.empty()) {
		if (false) qDebug() << "[自动直线控制笼] 合并后没有有效直线段。";
		return 0;
	}

	std::vector<AutoStraightCageBuildItem> buildItems;

	struct AutoStraightFrameState
	{
		bool valid = false;

		double endCenter[3] = { 0.0, 0.0, 0.0 };
		double axis[3] = { 0.0, 0.0, 0.0 };
		double up[3] = { 0.0, 0.0, 0.0 };
		double side[3] = { 0.0, 0.0, 0.0 };

		double halfWidth = 0.0;
		double halfDepth = 0.0;

		int lastBuildItemIndex = -1;
	};

	struct AutoStraightSharedPair
	{
		int previousIndex = -1;
		int currentIndex = -1;
	};

	std::map<int, AutoStraightFrameState> chainFrameStates;
	std::vector<AutoStraightSharedPair> sharedPairs;

	auto RingHalfExtentAlong =
		[](const double ring[4][3], const double center[3], const double dir[3]) -> double
	{
		double value = 0.0;

		for (int i = 0; i < 4; ++i) {
			double v[3] = {
				ring[i][0] - center[0],
				ring[i][1] - center[1],
				ring[i][2] - center[2]
			};

			value = std::max(value, std::abs(AutoStraightDot3(v, dir)));
		}

		return value;
	};

	auto ExtractRingUpSide =
		[](const double ring[4][3], double side[3], double up[3]) -> bool
	{
		double leftMid[3] = {
			0.5 * (ring[0][0] + ring[3][0]),
			0.5 * (ring[0][1] + ring[3][1]),
			0.5 * (ring[0][2] + ring[3][2])
		};

		double rightMid[3] = {
			0.5 * (ring[1][0] + ring[2][0]),
			0.5 * (ring[1][1] + ring[2][1]),
			0.5 * (ring[1][2] + ring[2][2])
		};

		double bottomMid[3] = {
			0.5 * (ring[0][0] + ring[1][0]),
			0.5 * (ring[0][1] + ring[1][1]),
			0.5 * (ring[0][2] + ring[1][2])
		};

		double topMid[3] = {
			0.5 * (ring[2][0] + ring[3][0]),
			0.5 * (ring[2][1] + ring[3][1]),
			0.5 * (ring[2][2] + ring[3][2])
		};

		side[0] = rightMid[0] - leftMid[0];
		side[1] = rightMid[1] - leftMid[1];
		side[2] = rightMid[2] - leftMid[2];

		up[0] = topMid[0] - bottomMid[0];
		up[1] = topMid[1] - bottomMid[1];
		up[2] = topMid[2] - bottomMid[2];

		return AutoStraightNormalize(side) && AutoStraightNormalize(up);
	};

	auto ProjectVectorToPlane =
		[](const double v[3], const double normal[3], double out[3]) -> bool
	{
		const double projection = AutoStraightDot3(v, normal);

		out[0] = v[0] - projection * normal[0];
		out[1] = v[1] - projection * normal[1];
		out[2] = v[2] - projection * normal[2];

		return AutoStraightNormalize(out);
	};

	auto BuildCornerSharedRing =
		[&](
			AutoStraightCageBuildItem& previousItem,
			AutoStraightCageBuildItem& currentItem)
	{
		double previousCenter[3];
		double currentCenter[3];

		AutoStraightRingCenter(previousItem.cage.pEnd, previousCenter);
		AutoStraightRingCenter(currentItem.cage.pStart, currentCenter);

		double sharedCenter[3] = {
			0.5 * (previousCenter[0] + currentCenter[0]),
			0.5 * (previousCenter[1] + currentCenter[1]),
			0.5 * (previousCenter[2] + currentCenter[2])
		};

		// 拐角共享截面的法向使用两段轴向的角平分方向。
		// 这样上一段 end 和下一段 start 是同一个面，
		// 但这个面不会被强行压成上一段或下一段的方向。
		double sharedAxis[3] = {
			previousItem.axis[0] + currentItem.axis[0],
			previousItem.axis[1] + currentItem.axis[1],
			previousItem.axis[2] + currentItem.axis[2]
		};

		if (!AutoStraightNormalize(sharedAxis)) {
			sharedAxis[0] = previousItem.axis[0];
			sharedAxis[1] = previousItem.axis[1];
			sharedAxis[2] = previousItem.axis[2];

			if (!AutoStraightNormalize(sharedAxis)) {
				return;
			}
		}

		double previousSide[3];
		double previousUp[3];

		double currentSide[3];
		double currentUp[3];

		const bool hasPreviousFrame =
			ExtractRingUpSide(
				previousItem.cage.pEnd,
				previousSide,
				previousUp);

		const bool hasCurrentFrame =
			ExtractRingUpSide(
				currentItem.cage.pStart,
				currentSide,
				currentUp);

		double sharedUp[3];
		double sharedSide[3];

		bool hasSharedUp = false;

		if (hasPreviousFrame) {
			hasSharedUp =
				ProjectVectorToPlane(
					previousUp,
					sharedAxis,
					sharedUp);
		}

		if (!hasSharedUp && hasCurrentFrame) {
			hasSharedUp =
				ProjectVectorToPlane(
					currentUp,
					sharedAxis,
					sharedUp);
		}

		if (!hasSharedUp) {
			AutoStraightBuildStableFrame(
				m_RefLocator,
				sharedCenter,
				sharedAxis,
				sharedUp,
				sharedSide);
		}
		else {
			vtkMath::Cross(sharedAxis, sharedUp, sharedSide);

			if (!AutoStraightNormalize(sharedSide)) {
				AutoStraightBuildStableFrame(
					m_RefLocator,
					sharedCenter,
					sharedAxis,
					sharedUp,
					sharedSide);
			}
			else {
				vtkMath::Cross(sharedSide, sharedAxis, sharedUp);
				AutoStraightNormalize(sharedUp);
			}
		}

		// 共享截面尺寸取两侧端面的较大值，避免拐角内圈/外圈被挤瘪。
		const double previousWidth =
			RingHalfExtentAlong(
				previousItem.cage.pEnd,
				previousCenter,
				sharedSide);

		const double previousDepth =
			RingHalfExtentAlong(
				previousItem.cage.pEnd,
				previousCenter,
				sharedUp);

		const double currentWidth =
			RingHalfExtentAlong(
				currentItem.cage.pStart,
				currentCenter,
				sharedSide);

		const double currentDepth =
			RingHalfExtentAlong(
				currentItem.cage.pStart,
				currentCenter,
				sharedUp);

		double sharedHalfWidth =
			std::max(previousWidth, currentWidth);

		double sharedHalfDepth =
			std::max(previousDepth, currentDepth);

		if (sharedHalfWidth <= 1e-8) {
			sharedHalfWidth =
				std::max(previousItem.halfWidth, currentItem.halfWidth);
		}

		if (sharedHalfDepth <= 1e-8) {
			sharedHalfDepth =
				std::max(previousItem.halfDepth, currentItem.halfDepth);
		}

		double rawSharedRing[4][3];
		double alignedSharedRing[4][3];

		AutoStraightMakeRing(
			sharedCenter,
			sharedUp,
			sharedSide,
			sharedHalfWidth,
			sharedHalfDepth,
			rawSharedRing);

		StandardizeCagePoints(
			rawSharedRing,
			sharedAxis,
			sharedUp,
			sharedSide);

		// 保持上一段 end 的点序，避免内外圈点顺序跳变。
		AlignRingByCyclicShift(
			previousItem.cage.pEnd,
			rawSharedRing,
			alignedSharedRing);

		AutoStraightCopyRingLocal(
			previousItem.cage.pEnd,
			alignedSharedRing);

		AutoStraightCopyRingLocal(
			currentItem.cage.pStart,
			alignedSharedRing);

		for (int d = 0; d < 3; ++d) {
			previousItem.endCenter[d] = sharedCenter[d];
			currentItem.startCenter[d] = sharedCenter[d];
		}
	};

	for (size_t i = 0; i < mergedSegments.size(); ++i) {
		const AutoStraightCageSegment& segment = mergedSegments[i];

		if (!segment.valid || segment.points.size() < 6) {
			continue;
		}

		double startCenter[3];
		double endCenter[3];
		double axis[3];

		if (!AutoStraightBuildAxisFromPoints(
			segment.points,
			startCenter,
			endCenter,
			axis))
		{
			continue;
		}

		AutoStraightFrameState& frameState =
			chainFrameStates[segment.sourceChainId];

		if (frameState.valid) {
			const double distToStart =
				AutoStraightDistance3(frameState.endCenter, startCenter);

			const double distToEnd =
				AutoStraightDistance3(frameState.endCenter, endCenter);

			if (distToEnd < distToStart) {
				for (int d = 0; d < 3; ++d) {
					std::swap(startCenter[d], endCenter[d]);
					axis[d] = -axis[d];
				}
			}
		}

		const double length =
			AutoStraightDistance3(startCenter, endCenter);

		if (length < std::max(modelDiagonal * 0.002, 1e-6)) {
			continue;
		}

		double mid[3] = {
			0.5 * (startCenter[0] + endCenter[0]),
			0.5 * (startCenter[1] + endCenter[1]),
			0.5 * (startCenter[2] + endCenter[2])
		};

		double up[3];
		double side[3];

		AutoStraightBuildStableFrame(
			m_RefLocator,
			mid,
			axis,
			up,
			side);

		bool connectedToPreviousFrame = false;

		if (frameState.valid) {
			const double connectDistance =
				AutoStraightDistance3(frameState.endCenter, startCenter);

			const double frameConnectLimit =
				std::max(modelDiagonal * 0.075, length * 0.80);

			if (connectDistance <= frameConnectLimit) {
				connectedToPreviousFrame = true;

				double propagatedUp[3] = {
					frameState.up[0],
					frameState.up[1],
					frameState.up[2]
				};

				const double projection =
					AutoStraightDot3(propagatedUp, axis);

				for (int d = 0; d < 3; ++d) {
					propagatedUp[d] -= projection * axis[d];
				}

				if (AutoStraightNormalize(propagatedUp)) {
					up[0] = propagatedUp[0];
					up[1] = propagatedUp[1];
					up[2] = propagatedUp[2];

					vtkMath::Cross(axis, up, side);

					if (AutoStraightNormalize(side)) {
						vtkMath::Cross(side, axis, up);
						AutoStraightNormalize(up);
					}
				}
				else {
					connectedToPreviousFrame = false;
				}
			}
		}

		double startHalfWidth = 0.0;
		double startHalfDepth = 0.0;
		double endHalfWidth = 0.0;
		double endHalfDepth = 0.0;

		if (templateMode == StraightCageTemplateMode::ThinRod) {
			double fallbackHalfSize =
				std::max(modelDiagonal * 0.012, length * 0.10);

			double minHalfSize =
				std::max(modelDiagonal * 0.004, length * 0.020);

			double maxHalfSize =
				std::max(
					minHalfSize * 2.0,
					std::min(
						modelDiagonal * 0.055,
						std::max(fallbackHalfSize * 2.2, length * 0.25)));

			double halfWidth =
				AutoStraightEstimateHalfExtentInDirection(
					m_RefLocator,
					startCenter,
					endCenter,
					side,
					modelDiagonal,
					fallbackHalfSize,
					minHalfSize,
					maxHalfSize);

			double halfDepth =
				AutoStraightEstimateHalfExtentInDirection(
					m_RefLocator,
					startCenter,
					endCenter,
					up,
					modelDiagonal,
					fallbackHalfSize,
					minHalfSize,
					maxHalfSize);

			if (connectedToPreviousFrame &&
				frameState.halfWidth > 1e-8 &&
				frameState.halfDepth > 1e-8)
			{
				halfWidth = std::max(halfWidth, frameState.halfWidth * 0.88);
				halfWidth = std::min(halfWidth, frameState.halfWidth * 1.18);

				halfDepth = std::max(halfDepth, frameState.halfDepth * 0.88);
				halfDepth = std::min(halfDepth, frameState.halfDepth * 1.18);
			}

			const double sectionScale = 1.08;

			halfWidth *= sectionScale;
			halfDepth *= sectionScale;

			startHalfWidth = halfWidth;
			startHalfDepth = halfDepth;
			endHalfWidth = halfWidth;
			endHalfDepth = halfDepth;
		}
		else {
			double fallbackHalfSize =
				std::max(modelDiagonal * 0.018, length * 0.16);

			double minHalfSize =
				std::max(modelDiagonal * 0.0035, length * 0.012);

			double maxHalfSize =
				std::max(
					minHalfSize * 2.0,
					std::max(fallbackHalfSize * 3.0, modelDiagonal * 0.35));

			double midHalfWidth = fallbackHalfSize;
			double midHalfDepth = fallbackHalfSize;
			startHalfWidth = fallbackHalfSize;
			startHalfDepth = fallbackHalfSize;
			endHalfWidth = fallbackHalfSize;
			endHalfDepth = fallbackHalfSize;

			bool startUsedSupport = false;
			bool midUsedSupport = false;
			bool endUsedSupport = false;

			AutoStraightEstimateSectionHalfExtentsAtPoint(
				m_RefLocator,
				m_skelNodes,
				m_skelEdges,
				m_skelSupportOffsets,
				m_skelSupportPoints,
				startCenter,
				startCenter,
				endCenter,
				axis,
				modelDiagonal,
				up,
				side,
				fallbackHalfSize,
				minHalfSize,
				maxHalfSize,
				startHalfWidth,
				startHalfDepth,
				startUsedSupport);

			AutoStraightEstimateSectionHalfExtentsAtPoint(
				m_RefLocator,
				m_skelNodes,
				m_skelEdges,
				m_skelSupportOffsets,
				m_skelSupportPoints,
				mid,
				startCenter,
				endCenter,
				axis,
				modelDiagonal,
				up,
				side,
				fallbackHalfSize,
				minHalfSize,
				maxHalfSize,
				midHalfWidth,
				midHalfDepth,
				midUsedSupport);

			AutoStraightEstimateSectionHalfExtentsAtPoint(
				m_RefLocator,
				m_skelNodes,
				m_skelEdges,
				m_skelSupportOffsets,
				m_skelSupportPoints,
				endCenter,
				startCenter,
				endCenter,
				axis,
				modelDiagonal,
				up,
				side,
				fallbackHalfSize,
				minHalfSize,
				maxHalfSize,
				endHalfWidth,
				endHalfDepth,
				endUsedSupport);

			startHalfWidth = std::max(startHalfWidth, midHalfWidth * 0.72);
			startHalfDepth = std::max(startHalfDepth, midHalfDepth * 0.72);
			endHalfWidth = std::max(endHalfWidth, midHalfWidth * 0.72);
			endHalfDepth = std::max(endHalfDepth, midHalfDepth * 0.72);

			if (connectedToPreviousFrame &&
				frameState.halfWidth > 1e-8 &&
				frameState.halfDepth > 1e-8)
			{
				const bool hasSupportDrivenSize =
					startUsedSupport || midUsedSupport || endUsedSupport;

				if (!hasSupportDrivenSize) {
					startHalfWidth = std::max(startHalfWidth, frameState.halfWidth * 0.88);
					startHalfWidth = std::min(startHalfWidth, frameState.halfWidth * 1.18);
					endHalfWidth = std::max(endHalfWidth, frameState.halfWidth * 0.88);
					endHalfWidth = std::min(endHalfWidth, frameState.halfWidth * 1.18);

					startHalfDepth = std::max(startHalfDepth, frameState.halfDepth * 0.88);
					startHalfDepth = std::min(startHalfDepth, frameState.halfDepth * 1.18);
					endHalfDepth = std::max(endHalfDepth, frameState.halfDepth * 0.88);
					endHalfDepth = std::min(endHalfDepth, frameState.halfDepth * 1.18);
				}
			}

			AutoStraightFitAxisToBoundary(
				m_RefLocator,
				startCenter,
				side,
				modelDiagonal,
				startHalfWidth,
				maxHalfSize);

			AutoStraightFitAxisToBoundary(
				m_RefLocator,
				startCenter,
				up,
				modelDiagonal,
				startHalfDepth,
				maxHalfSize);

			AutoStraightFitAxisToBoundary(
				m_RefLocator,
				endCenter,
				side,
				modelDiagonal,
				endHalfWidth,
				maxHalfSize);

			AutoStraightFitAxisToBoundary(
				m_RefLocator,
				endCenter,
				up,
				modelDiagonal,
				endHalfDepth,
				maxHalfSize);

			const double sectionScale = 1.0;

			startHalfWidth = AutoStraightClampDoubleLocal(startHalfWidth * sectionScale, minHalfSize, maxHalfSize);
			startHalfDepth = AutoStraightClampDoubleLocal(startHalfDepth * sectionScale, minHalfSize, maxHalfSize);
			endHalfWidth = AutoStraightClampDoubleLocal(endHalfWidth * sectionScale, minHalfSize, maxHalfSize);
			endHalfDepth = AutoStraightClampDoubleLocal(endHalfDepth * sectionScale, minHalfSize, maxHalfSize);
		}

		double startRing[4][3];
		double endRingRaw[4][3];
		double endRing[4][3];

		// 注意：
		// 这里不再直接复制上一段 endRing。
		// 每段先生成规则端面，后面统一用 BuildCornerSharedRing
		// 把连续处改成角平分方向共享端面。

		AutoStraightMakeRing(
			startCenter,
			up,
			side,
			startHalfWidth,
			startHalfDepth,
			startRing);

		StandardizeCagePoints(
			startRing,
			axis,
			up,
			side);

		if (AutoStraightClampRingToNonDesignSurface(
			m_nonDesignReferenceSurface,
			startCenter,
			up,
			side,
			modelDiagonal,
			startHalfWidth,
			startHalfDepth,
			startRing))
		{
			StandardizeCagePoints(startRing, axis, up, side);
		}

		AutoStraightMakeRing(
			endCenter,
			up,
			side,
			endHalfWidth,
			endHalfDepth,
			endRingRaw);

		StandardizeCagePoints(
			endRingRaw,
			axis,
			up,
			side);

		AlignRingByCyclicShift(
			startRing,
			endRingRaw,
			endRing);

		if (AutoStraightClampRingToNonDesignSurface(
			m_nonDesignReferenceSurface,
			endCenter,
			up,
			side,
			modelDiagonal,
			endHalfWidth,
			endHalfDepth,
			endRing))
		{
			AlignRingByCyclicShift(startRing, endRing, endRingRaw);
			AutoStraightCopyRingLocal(endRing, endRingRaw);
		}

		AutoStraightCageBuildItem item;

		AutoStraightCopyRingLocal(item.cage.pStart, startRing);
		AutoStraightCopyRingLocal(item.cage.pEnd, endRing);

		for (int d = 0; d < 3; ++d) {
			item.startCenter[d] = startCenter[d];
			item.endCenter[d] = endCenter[d];
			item.axis[d] = axis[d];
		}

		item.length = length;
		item.sourceChainId = segment.sourceChainId;
		item.halfWidth = std::max(startHalfWidth, endHalfWidth);
		item.halfDepth = std::max(startHalfDepth, endHalfDepth);

		const int previousIndex = frameState.lastBuildItemIndex;
		const int currentIndex = static_cast<int>(buildItems.size());

		buildItems.push_back(item);

		if (connectedToPreviousFrame &&
			previousIndex >= 0 &&
			previousIndex < currentIndex)
		{
			AutoStraightSharedPair pair;
			pair.previousIndex = previousIndex;
			pair.currentIndex = currentIndex;
			sharedPairs.push_back(pair);
		}

		frameState.valid = true;

		for (int d = 0; d < 3; ++d) {
			frameState.endCenter[d] = endCenter[d];
			frameState.axis[d] = axis[d];
			frameState.up[d] = up[d];
			frameState.side[d] = side[d];
		}

		frameState.halfWidth = item.halfWidth;
		frameState.halfDepth = item.halfDepth;
		frameState.lastBuildItemIndex = currentIndex;
	}

	if (buildItems.empty()) {
		if (false) qDebug() << "[自动直线控制笼] 没有生成任何有效控制笼。";
		return 0;
	}

	// 先处理不同 chain 的端点冲突。
	AutoStraightResolveCrossChainEndpointConflicts(
		buildItems,
		modelDiagonal);

	// 再强制恢复同一 sourceChainId 内的连续性。
	// 每个连续处使用角平分方向共享 ring。
	for (const AutoStraightSharedPair& pair : sharedPairs) {
		if (pair.previousIndex < 0 ||
			pair.currentIndex < 0 ||
			pair.previousIndex >= static_cast<int>(buildItems.size()) ||
			pair.currentIndex >= static_cast<int>(buildItems.size()))
		{
			continue;
		}

		AutoStraightCageBuildItem& previousItem =
			buildItems[static_cast<size_t>(pair.previousIndex)];

		AutoStraightCageBuildItem& currentItem =
			buildItems[static_cast<size_t>(pair.currentIndex)];

		if (previousItem.sourceChainId != currentItem.sourceChainId) {
			continue;
		}

		BuildCornerSharedRing(previousItem, currentItem);
	}

	std::vector<CageSegment> newCages;
	newCages.reserve(buildItems.size());

	for (const AutoStraightCageBuildItem& item : buildItems) {
		newCages.push_back(item.cage);
	}

	if (newCages.empty()) {
		if (false) qDebug() << "[自动直线控制笼] 处理后没有有效控制笼。";
		return 0;
	}

	if (replaceExisting) {
		m_allCages.clear();
		m_hubCages.clear();
		m_suppressedCageSideFaces.clear();
		m_sdfHubs.clear();
		m_use3DHubSdfProjection = false;
		m_selectedFacesForBranch.clear();

		m_useForcedLoopTopology = false;
		m_forcedLoopCycleNodes.clear();
		m_forcedLoopConnectorPaths.clear();
	}

	for (const CageSegment& cage : newCages) {
		m_allCages.push_back(cage);
	}

	m_ClickState = 0;
	m_ManualRotationOffset = 0;
	m_IsSnapped = false;
	m_SnappedIsStart = false;
	m_SnappedCageIndex = -1;
	m_SnappedFaceOffset = -1;

	if (m_FixedBoxActor) {
		m_FixedBoxActor->VisibilityOff();
	}

	if (m_LinesActor) {
		m_LinesActor->VisibilityOff();
	}

	if (m_PreviewBoxActor) {
		m_PreviewBoxActor->VisibilityOff();
	}

	RebuildGlobalMesh();

	if (this->GetInteractor() &&
		this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}

	if (false) qDebug() << "[自动直线控制笼] 原始直线段数量 =" << segments.size()
		<< "合并后数量 =" << mergedSegments.size()
		<< "生成 cage 数量 =" << newCages.size()
		<< "共享连续端面数量 =" << sharedPairs.size()
		<< "当前 m_allCages 数量 =" << m_allCages.size();

	return static_cast<int>(newCages.size());
}

void CageInteractorStyle::SetStraightSegmentsForEditing(
	std::vector<StraightSkeletonSegmenter::Segment>* segments,
	vtkSmartPointer<vtkActor> segmentActor,
	const std::vector<double>& skeletonNodes)
{
	m_editStraightSegments = segments;
	m_editStraightSegmentActor = segmentActor;
	m_editStraightSkeletonNodes = skeletonNodes;
}
void CageInteractorStyle::EnableStraightSegmentEditMode(bool enable)
{
	if (m_isStraightSegmentEditMode == enable) {
		return;
	}

	m_isStraightSegmentEditMode = enable;
	if (enable) {
		m_isStraightSegmentSelectMode = false;
		m_isStraightCageEditMode = false;
		m_isStraightCageRotateMode = false;
		CancelStraightCageEdit(nullptr);
		CancelStraightCageRotation(nullptr);
		HideManualPreviewActors();
	}

	if (m_editStraightSegmentActor) {
		m_editStraightSegmentActor->SetPickable(enable ? 1 : 0);
	}

	if (!enable) {
		m_ClickState = 0;

		if (m_FixedBoxActor) {
			m_FixedBoxActor->VisibilityOff();
		}

		if (m_LinesActor) {
			m_LinesActor->VisibilityOff();
		}
	}

	if (this->GetInteractor() &&
		this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}

	NotifyEditModeChanged();
	if (false) qDebug() << "[直线段编辑]" << (enable ? "进入编辑模式" : "退出编辑模式");
}
void CageInteractorStyle::EnableStraightSegmentSelectMode(bool enable)
{
	if (m_isStraightSegmentSelectMode == enable) {
		return;
	}

	m_isStraightSegmentSelectMode = enable;
	if (enable) {
		m_isStraightSegmentEditMode = false;
		m_isStraightCageEditMode = false;
		m_isStraightCageRotateMode = false;
		CancelStraightCageEdit(nullptr);
		CancelStraightCageRotation(nullptr);
		HideManualPreviewActors();
	}

	if (m_editStraightSegmentActor) {
		m_editStraightSegmentActor->SetPickable(0);
	}

	if (this->GetInteractor() &&
		this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}

	NotifyEditModeChanged();
}

void CageInteractorStyle::RefreshStraightSegmentEditActor()
{
	if (!m_editStraightSegmentActor || !m_editStraightSegments) {
		return;
	}

	vtkSmartPointer<vtkActor> newActor =
		StraightSkeletonSegmenter::CreateSegmentActor(
			m_editStraightSkeletonNodes,
			*m_editStraightSegments);

	if (!newActor || !newActor->GetMapper()) {
		return;
	}

	m_editStraightSegmentActor->SetMapper(newActor->GetMapper());
	m_editStraightSegmentActor->GetProperty()->SetColor(0.0, 1.0, 0.25);
	m_editStraightSegmentActor->GetProperty()->SetLineWidth(8.0);
	m_editStraightSegmentActor->GetProperty()->SetRenderLinesAsTubes(true);
	m_editStraightSegmentActor->GetProperty()->SetLighting(false);
	m_editStraightSegmentActor->SetPickable(m_isStraightSegmentEditMode ? 1 : 0);
}
bool CageInteractorStyle::HandleStraightSegmentTrimClick(int x, int y, vtkRenderer* ren)
{
	if (!ren || !m_editStraightSegments || !m_editStraightSegmentActor) {
		return false;
	}

	if (m_editStraightSegments->empty()) {
		return false;
	}

	vtkNew<vtkCellPicker> picker;
	picker->SetTolerance(0.01);

	m_editStraightSegmentActor->SetPickable(1);

	picker->AddPickList(m_editStraightSegmentActor);
	picker->PickFromListOn();

	const int picked = picker->Pick(x, y, 0, ren);

	if (!picked || picker->GetActor() != m_editStraightSegmentActor) {
		return false;
	}

	double pickedPoint[3];
	picker->GetPickPosition(pickedPoint);

	const double modelDiagonal = StraightEditSkeletonDiagonal(m_editStraightSkeletonNodes);

	int bestSegmentIndex = -1;
	double bestDistance = std::numeric_limits<double>::max();

	for (int si = 0; si < static_cast<int>(m_editStraightSegments->size()); ++si) {
		const StraightSkeletonSegmenter::Segment& segment = (*m_editStraightSegments)[si];

		if (!segment.valid || segment.points.size() < 6) {
			continue;
		}

		const std::vector<double>& pts = segment.points;

		for (size_t i = 3; i + 2 < pts.size(); i += 3) {
			double a[3] = { pts[i - 3], pts[i - 2], pts[i - 1] };
			double b[3] = { pts[i], pts[i + 1], pts[i + 2] };

			double ab[3] = {
				b[0] - a[0],
				b[1] - a[1],
				b[2] - a[2]
			};

			const double abLen2 = StraightEditDot3(ab, ab);
			if (abLen2 <= 1e-12) {
				continue;
			}

			double ap[3] = {
				pickedPoint[0] - a[0],
				pickedPoint[1] - a[1],
				pickedPoint[2] - a[2]
			};

			double t = StraightEditDot3(ap, ab) / abLen2;
			t = std::max(0.0, std::min(1.0, t));

			double projected[3] = {
				a[0] + ab[0] * t,
				a[1] + ab[1] * t,
				a[2] + ab[2] * t
			};

			const double distance = StraightEditDistance3(projected, pickedPoint);

			if (distance < bestDistance) {
				bestDistance = distance;
				bestSegmentIndex = si;
			}
		}
	}

	if (bestSegmentIndex < 0) {
		return false;
	}

	bool changed = StraightEditTrimSegmentAtPoint(
		(*m_editStraightSegments)[bestSegmentIndex],
		pickedPoint,
		modelDiagonal);

	if (!changed) {
		return false;
	}

	RefreshStraightSegmentEditActor();

	if (this->GetInteractor() &&
		this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}


	return true;
}

bool CageInteractorStyle::BuildStraightSegmentFromSkeletonEdge(
	int edgeIndex,
	StraightSkeletonSegmenter::Segment& segment) const
{
	segment = StraightSkeletonSegmenter::Segment();
	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
	const int edgeCount = static_cast<int>(m_skelEdges.size() / 2);
	if (edgeIndex < 0 || edgeIndex >= edgeCount || nodeCount < 2) {
		return false;
	}

	const int seedA = static_cast<int>(m_skelEdges[2 * edgeIndex]);
	const int seedB = static_cast<int>(m_skelEdges[2 * edgeIndex + 1]);
	if (seedA < 0 || seedB < 0 || seedA >= nodeCount || seedB >= nodeCount || seedA == seedB) {
		return false;
	}

	std::vector<std::vector<int>> adjacency(nodeCount);
	for (int ei = 0; ei < edgeCount; ++ei) {
		const int a = static_cast<int>(m_skelEdges[2 * ei]);
		const int b = static_cast<int>(m_skelEdges[2 * ei + 1]);
		if (a < 0 || b < 0 || a >= nodeCount || b >= nodeCount || a == b) {
			continue;
		}
		adjacency[a].push_back(b);
		adjacency[b].push_back(a);
	}

	auto getPoint = [&](int node, double p[3])
	{
		p[0] = m_skelNodes[3 * node];
		p[1] = m_skelNodes[3 * node + 1];
		p[2] = m_skelNodes[3 * node + 2];
	};

	double pA[3], pB[3];
	getPoint(seedA, pA);
	getPoint(seedB, pB);
	double seedDir[3] = {
		pB[0] - pA[0],
		pB[1] - pA[1],
		pB[2] - pA[2]
	};
	if (vtkMath::Normalize(seedDir) == 0.0) {
		return false;
	}

	const double minCos = std::cos(15.0 * 3.14159265358979323846 / 180.0);
	auto alignedWithSeed = [&](int from, int to) -> bool
	{
		double pf[3], pt[3];
		getPoint(from, pf);
		getPoint(to, pt);
		double dir[3] = {
			pt[0] - pf[0],
			pt[1] - pf[1],
			pt[2] - pf[2]
		};
		if (vtkMath::Normalize(dir) == 0.0) {
			return false;
		}
		return std::abs(vtkMath::Dot(dir, seedDir)) >= minCos;
	};

	auto walkAway = [&](int curr, int prev)
	{
		std::vector<int> walked;
		walked.push_back(curr);
		std::set<int> visited;
		visited.insert(prev);
		visited.insert(curr);

		while (curr >= 0 &&
			curr < nodeCount &&
			adjacency[curr].size() == 2)
		{
			int next = -1;
			for (int candidate : adjacency[curr]) {
				if (candidate != prev) {
					next = candidate;
					break;
				}
			}

			if (next < 0 ||
				visited.find(next) != visited.end() ||
				!alignedWithSeed(curr, next))
			{
				break;
			}

			walked.push_back(next);
			visited.insert(next);
			prev = curr;
			curr = next;
		}

		return walked;
	};

	std::vector<int> left = walkAway(seedA, seedB);
	std::vector<int> right = walkAway(seedB, seedA);
	std::reverse(left.begin(), left.end());

	segment.nodes = left;
	segment.nodes.insert(segment.nodes.end(), right.begin(), right.end());
	if (segment.nodes.size() < 2) {
		return false;
	}

	for (int node : segment.nodes) {
		double p[3];
		getPoint(node, p);
		segment.points.push_back(p[0]);
		segment.points.push_back(p[1]);
		segment.points.push_back(p[2]);
	}

	if (StraightEditPolylineLength(segment.points) < std::max(StraightEditSkeletonDiagonal(m_skelNodes) * 0.004, 1e-6)) {
		return false;
	}

	segment.originalPoints = segment.points;
	segment.valid = true;
	segment.sourceChainId = 100000 + edgeIndex;
	return true;
}

bool CageInteractorStyle::HandleStraightSegmentSelectClick(int x, int y, vtkRenderer* ren)
{
	if (!ren || !m_editStraightSegments) {
		return false;
	}

	const int edgeIndex = PickSkeletonEdge(x, y, ren);
	if (edgeIndex < 0) {
		return false;
	}

	StraightSkeletonSegmenter::Segment newSegment;
	if (!BuildStraightSegmentFromSkeletonEdge(edgeIndex, newSegment)) {
		return false;
	}

	const double modelDiagonal = StraightEditSkeletonDiagonal(m_skelNodes);
	const double duplicateTol2 = std::pow(std::max(modelDiagonal * 0.006, 1e-6), 2.0);
	auto samePoint = [&](const double* a, const double* b)
	{
		return vtkMath::Distance2BetweenPoints(a, b) <= duplicateTol2;
	};

	const double newStart[3] = { newSegment.points[0], newSegment.points[1], newSegment.points[2] };
	const size_t n = newSegment.points.size();
	const double newEnd[3] = { newSegment.points[n - 3], newSegment.points[n - 2], newSegment.points[n - 1] };

	for (const StraightSkeletonSegmenter::Segment& existing : *m_editStraightSegments) {
		if (!existing.valid || existing.points.size() < 6) {
			continue;
		}
		const double oldStart[3] = { existing.points[0], existing.points[1], existing.points[2] };
		const size_t m = existing.points.size();
		const double oldEnd[3] = { existing.points[m - 3], existing.points[m - 2], existing.points[m - 1] };

		const bool sameDirection = samePoint(newStart, oldStart) && samePoint(newEnd, oldEnd);
		const bool reverseDirection = samePoint(newStart, oldEnd) && samePoint(newEnd, oldStart);
		if (sameDirection || reverseDirection) {
			return true;
		}
	}

	m_editStraightSegments->push_back(newSegment);
	RefreshStraightSegmentEditActor();

	if (this->GetInteractor() &&
		this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}

	return true;
}

void CageInteractorStyle::NotifyEditModeChanged()
{
	if (m_editModeChangedCallback) {
		m_editModeChangedCallback();
	}
}

void CageInteractorStyle::HideManualPreviewActors()
{
	if (m_PreviewBoxActor) m_PreviewBoxActor->VisibilityOff();
	if (m_FixedBoxActor) m_FixedBoxActor->VisibilityOff();
	if (m_LinesActor) m_LinesActor->VisibilityOff();
}

void CageInteractorStyle::ExitAllEditModes()
{
	m_isSkeletonEditMode = false;
	m_isSkeletonDeleteMode = false;
	m_isSkeletonExtendMode = false;
	m_isStraightSegmentEditMode = false;
	m_isStraightSegmentSelectMode = false;
	m_isStraightCageEditMode = false;
	m_isStraightCageRotateMode = false;
	m_selectedSkelNode1 = -1;
	m_selectedExtendSkelNode = -1;
	m_ClickState = 0;
	HideSkeletonExtendPreview();
	HideManualPreviewActors();
	CancelStraightCageEdit(nullptr);
	CancelStraightCageRotation(nullptr);
	if (m_editStraightSegmentActor) {
		m_editStraightSegmentActor->SetPickable(0);
	}
	NotifyEditModeChanged();
}

void CageInteractorStyle::EnableStraightCageEditMode(bool enable)
{
	if (m_isStraightCageEditMode == enable) {
		return;
	}

	m_isStraightCageEditMode = enable;
	if (enable) {
		m_isSkeletonEditMode = false;
		m_isSkeletonDeleteMode = false;
		m_isSkeletonExtendMode = false;
		m_isStraightSegmentEditMode = false;
		m_isStraightSegmentSelectMode = false;
		m_isStraightCageRotateMode = false;
		CancelStraightCageRotation(nullptr);
		m_selectedSkelNode1 = -1;
		m_selectedExtendSkelNode = -1;
		HideSkeletonExtendPreview();
		if (m_editStraightSegmentActor) {
			m_editStraightSegmentActor->SetPickable(0);
		}
	}
	else {
		CancelStraightCageEdit(nullptr);
	}
	HideManualPreviewActors();
	NotifyEditModeChanged();
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::EnableStraightCageRotateMode(bool enable)
{
	if (m_isStraightCageRotateMode == enable) {
		return;
	}

	m_isStraightCageRotateMode = enable;
	if (enable) {
		m_isSkeletonEditMode = false;
		m_isSkeletonDeleteMode = false;
		m_isSkeletonExtendMode = false;
		m_isStraightSegmentEditMode = false;
		m_isStraightSegmentSelectMode = false;
		m_isStraightCageEditMode = false;
		CancelStraightCageEdit(nullptr);
		m_selectedSkelNode1 = -1;
		m_selectedExtendSkelNode = -1;
		HideSkeletonExtendPreview();
		if (m_editStraightSegmentActor) {
			m_editStraightSegmentActor->SetPickable(0);
		}
	}
	else {
		CancelStraightCageRotation(nullptr);
	}
	HideManualPreviewActors();
	NotifyEditModeChanged();
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

namespace
{
	static void StraightCageEditFaceIds(int faceIndex, int ids[4])
	{
		static const int faceTable[6][4] = {
			{ 0, 1, 2, 3 },
			{ 4, 5, 6, 7 },
			{ 0, 1, 5, 4 },
			{ 1, 2, 6, 5 },
			{ 2, 3, 7, 6 },
			{ 3, 0, 4, 7 }
		};

		const int safeFace = std::max(0, std::min(5, faceIndex));
		for (int i = 0; i < 4; ++i) {
			ids[i] = faceTable[safeFace][i];
		}
	}

	static void StraightCageEditGetPoint(const CageSegment& cage, int localId, double p[3])
	{
		const int id = std::max(0, std::min(7, localId));
		const double* src = id < 4 ? cage.pStart[id] : cage.pEnd[id - 4];
		p[0] = src[0];
		p[1] = src[1];
		p[2] = src[2];
	}

	static void StraightCageEditSetPoint(CageSegment& cage, int localId, const double p[3])
	{
		const int id = std::max(0, std::min(7, localId));
		double* dst = id < 4 ? cage.pStart[id] : cage.pEnd[id - 4];
		dst[0] = p[0];
		dst[1] = p[1];
		dst[2] = p[2];
	}

	static void StraightCageEditCageCenter(const CageSegment& cage, double center[3])
	{
		center[0] = center[1] = center[2] = 0.0;
		for (int i = 0; i < 8; ++i) {
			double p[3];
			StraightCageEditGetPoint(cage, i, p);
			center[0] += p[0] / 8.0;
			center[1] += p[1] / 8.0;
			center[2] += p[2] / 8.0;
		}
	}

	static double StraightCageEditPointSegmentDistance2D(
		double px,
		double py,
		double ax,
		double ay,
		double bx,
		double by)
	{
		const double vx = bx - ax;
		const double vy = by - ay;
		const double len2 = vx * vx + vy * vy;
		if (len2 <= 1e-12) {
			const double dx = px - ax;
			const double dy = py - ay;
			return std::sqrt(dx * dx + dy * dy);
		}

		double t = ((px - ax) * vx + (py - ay) * vy) / len2;
		t = std::max(0.0, std::min(1.0, t));
		const double qx = ax + vx * t;
		const double qy = ay + vy * t;
		const double dx = px - qx;
		const double dy = py - qy;
		return std::sqrt(dx * dx + dy * dy);
	}

	static void StraightCageEditCentersAxis(
		const CageSegment& cage,
		double startCenter[3],
		double endCenter[3],
		double axis[3])
	{
		startCenter[0] = startCenter[1] = startCenter[2] = 0.0;
		endCenter[0] = endCenter[1] = endCenter[2] = 0.0;
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				startCenter[d] += cage.pStart[i][d] * 0.25;
				endCenter[d] += cage.pEnd[i][d] * 0.25;
			}
		}
		for (int d = 0; d < 3; ++d) {
			axis[d] = endCenter[d] - startCenter[d];
		}
		vtkMath::Normalize(axis);
	}

	static void StraightCageRotatePoint(
		double p[3],
		const double center[3],
		const double axis[3],
		double angleRadians)
	{
		const double c = std::cos(angleRadians);
		const double s = std::sin(angleRadians);
		double v[3] = {
			p[0] - center[0],
			p[1] - center[1],
			p[2] - center[2]
		};
		const double parallel = vtkMath::Dot(v, axis);
		double vParallel[3] = {
			axis[0] * parallel,
			axis[1] * parallel,
			axis[2] * parallel
		};
		double vPerp[3] = {
			v[0] - vParallel[0],
			v[1] - vParallel[1],
			v[2] - vParallel[2]
		};
		double cross[3];
		vtkMath::Cross(axis, vPerp, cross);
		for (int d = 0; d < 3; ++d) {
			p[d] = center[d] + vParallel[d] + vPerp[d] * c + cross[d] * s;
		}
	}

	static CageSegment StraightCageRotatedCopy(const CageSegment& cage, double angleRadians)
	{
		CageSegment rotated = cage;
		double startCenter[3], endCenter[3], axis[3];
		StraightCageEditCentersAxis(cage, startCenter, endCenter, axis);
		if (vtkMath::Norm(axis) <= 1e-8) {
			return rotated;
		}
		for (int i = 0; i < 4; ++i) {
			StraightCageRotatePoint(rotated.pStart[i], startCenter, axis, angleRadians);
			StraightCageRotatePoint(rotated.pEnd[i], endCenter, axis, angleRadians);
		}
		return rotated;
	}

	static vtkSmartPointer<vtkPolyData> StraightCageBuildWholeData(const CageSegment& cage)
	{
		static const int faceTable[6][4] = {
			{ 0, 1, 2, 3 },
			{ 4, 7, 6, 5 },
			{ 0, 4, 5, 1 },
			{ 1, 5, 6, 2 },
			{ 2, 6, 7, 3 },
			{ 3, 7, 4, 0 }
		};

		vtkNew<vtkPoints> points;
		for (int i = 0; i < 4; ++i) {
			points->InsertNextPoint(cage.pStart[i]);
		}
		for (int i = 0; i < 4; ++i) {
			points->InsertNextPoint(cage.pEnd[i]);
		}

		vtkNew<vtkCellArray> polys;
		for (const auto& face : faceTable) {
			polys->InsertNextCell(4);
			for (int id : face) {
				polys->InsertCellPoint(id);
			}
		}

		vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
		data->SetPoints(points);
		data->SetPolys(polys);
		return data;
	}
}

bool CageInteractorStyle::PickStraightCageEditTarget(
	int x,
	int y,
	vtkRenderer* ren,
	StraightCageEditTarget& target)
{
	target = StraightCageEditTarget();
	if (!ren || (m_allCages.empty() && m_hubCages.empty())) {
		return false;
	}

	static const int edgeTable[12][2] = {
		{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
		{ 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
		{ 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
	};

	double bestEdgeDistance = std::numeric_limits<double>::max();
	int bestCage = -1;
	int bestA = -1;
	int bestB = -1;
	double bestVertexDistance = std::numeric_limits<double>::max();
	int bestVertexCage = -1;
	int bestVertexId = -1;

	for (int ci = 0; ci < static_cast<int>(m_allCages.size()); ++ci) {
		const CageSegment& cage = m_allCages[ci];
		for (int vi = 0; vi < 8; ++vi) {
			double p[3];
			StraightCageEditGetPoint(cage, vi, p);
			ren->SetWorldPoint(p[0], p[1], p[2], 1.0);
			ren->WorldToDisplay();
			double* display = ren->GetDisplayPoint();
			const double dx = display[0] - static_cast<double>(x);
			const double dy = display[1] - static_cast<double>(y);
			const double distance = std::sqrt(dx * dx + dy * dy);
			if (distance < bestVertexDistance) {
				bestVertexDistance = distance;
				bestVertexCage = ci;
				bestVertexId = vi;
			}
		}
		for (const auto& edge : edgeTable) {
			double a[3], b[3];
			StraightCageEditGetPoint(cage, edge[0], a);
			StraightCageEditGetPoint(cage, edge[1], b);

			ren->SetWorldPoint(a[0], a[1], a[2], 1.0);
			ren->WorldToDisplay();
			double* da = ren->GetDisplayPoint();
			double ax = da[0];
			double ay = da[1];

			ren->SetWorldPoint(b[0], b[1], b[2], 1.0);
			ren->WorldToDisplay();
			double* db = ren->GetDisplayPoint();
			double bx = db[0];
			double by = db[1];

			const double distance =
				StraightCageEditPointSegmentDistance2D(
					static_cast<double>(x),
					static_cast<double>(y),
					ax,
					ay,
					bx,
					by);

			if (distance < bestEdgeDistance) {
				bestEdgeDistance = distance;
				bestCage = ci;
				bestA = edge[0];
				bestB = edge[1];
			}
		}
	}

	double bestHubDistance = std::numeric_limits<double>::max();
	int bestHub = -1;
	int bestHubFace = -1;
	vtkIdType bestHubA = -1;
	vtkIdType bestHubB = -1;
	double bestHubVertexDistance = std::numeric_limits<double>::max();
	int bestHubVertex = -1;
	vtkIdType bestHubPointId = -1;
	double bestHubFaceDepth = std::numeric_limits<double>::max();
	int bestHubFaceHit = -1;
	int bestHubFaceHitCell = -1;

	for (int hi = 0; hi < static_cast<int>(m_hubCages.size()); ++hi) {
		vtkPolyData* hub = m_hubCages[hi];
		if (!hub || !hub->GetPoints()) {
			continue;
		}

		for (vtkIdType pi = 0; pi < hub->GetNumberOfPoints(); ++pi) {
			double p[3];
			hub->GetPoint(pi, p);
			ren->SetWorldPoint(p[0], p[1], p[2], 1.0);
			ren->WorldToDisplay();
			double* display = ren->GetDisplayPoint();
			const double dx = display[0] - static_cast<double>(x);
			const double dy = display[1] - static_cast<double>(y);
			const double distance = std::sqrt(dx * dx + dy * dy);
			if (distance < bestHubVertexDistance) {
				bestHubVertexDistance = distance;
				bestHubVertex = hi;
				bestHubPointId = pi;
			}
		}

		for (vtkIdType fi = 0; fi < hub->GetNumberOfCells(); ++fi) {
			vtkCell* cell = hub->GetCell(fi);
			if (!cell || cell->GetNumberOfPoints() < 3) {
				continue;
			}

			const vtkIdType n = cell->GetNumberOfPoints();
			std::vector<std::array<double, 3>> displayPolygon;
			displayPolygon.reserve(static_cast<size_t>(n));
			double averageDepth = 0.0;
			for (vtkIdType k = 0; k < n; ++k) {
				double p[3];
				hub->GetPoint(cell->GetPointId(k), p);
				ren->SetWorldPoint(p[0], p[1], p[2], 1.0);
				ren->WorldToDisplay();
				double* display = ren->GetDisplayPoint();
				displayPolygon.push_back({ display[0], display[1], display[2] });
				averageDepth += display[2] / static_cast<double>(n);
			}

			bool inside = false;
			for (vtkIdType k = 0, j = n - 1; k < n; j = k++) {
				const double xi = displayPolygon[static_cast<size_t>(k)][0];
				const double yi = displayPolygon[static_cast<size_t>(k)][1];
				const double xj = displayPolygon[static_cast<size_t>(j)][0];
				const double yj = displayPolygon[static_cast<size_t>(j)][1];
				const bool crosses =
					((yi > static_cast<double>(y)) != (yj > static_cast<double>(y))) &&
					(static_cast<double>(x) <
					(xj - xi) * (static_cast<double>(y) - yi) /
						((yj - yi) == 0.0 ? 1e-12 : (yj - yi)) + xi);
				if (crosses) {
					inside = !inside;
				}
			}
			if (inside && averageDepth < bestHubFaceDepth) {
				bestHubFaceDepth = averageDepth;
				bestHubFaceHit = hi;
				bestHubFaceHitCell = static_cast<int>(fi);
			}

			for (vtkIdType k = 0; k < n; ++k) {
				const vtkIdType aId = cell->GetPointId(k);
				const vtkIdType bId = cell->GetPointId((k + 1) % n);
				double a[3], b[3];
				hub->GetPoint(aId, a);
				hub->GetPoint(bId, b);

				ren->SetWorldPoint(a[0], a[1], a[2], 1.0);
				ren->WorldToDisplay();
				double* da = ren->GetDisplayPoint();
				double ax = da[0];
				double ay = da[1];

				ren->SetWorldPoint(b[0], b[1], b[2], 1.0);
				ren->WorldToDisplay();
				double* db = ren->GetDisplayPoint();
				double bx = db[0];
				double by = db[1];

				const double distance =
					StraightCageEditPointSegmentDistance2D(
						static_cast<double>(x),
						static_cast<double>(y),
						ax,
						ay,
						bx,
						by);

				if (distance < bestHubDistance) {
					bestHubDistance = distance;
					bestHub = hi;
					bestHubFace = static_cast<int>(fi);
					bestHubA = aId;
					bestHubB = bId;
				}
			}
		}

		// N 边柱 hub 的上下截面分别存为多边形 cell，连接上下截面的边
		// (0,1)、(2,3)... 不属于任何 cell 周边，必须显式加入拾取候选。
		const vtkIdType pointCount = hub->GetNumberOfPoints();
		if (pointCount >= 6 && pointCount % 2 == 0) {
			for (vtkIdType aId = 0; aId + 1 < pointCount; aId += 2) {
				const vtkIdType bId = aId + 1;
				double a[3], b[3];
				hub->GetPoint(aId, a);
				hub->GetPoint(bId, b);

				ren->SetWorldPoint(a[0], a[1], a[2], 1.0);
				ren->WorldToDisplay();
				double* da = ren->GetDisplayPoint();
				double ax = da[0], ay = da[1];
				ren->SetWorldPoint(b[0], b[1], b[2], 1.0);
				ren->WorldToDisplay();
				double* db = ren->GetDisplayPoint();

				double distance = StraightCageEditPointSegmentDistance2D(
					static_cast<double>(x), static_cast<double>(y),
					ax, ay, db[0], db[1]);
				if (distance < bestHubDistance) {
					bestHubDistance = distance;
					bestHub = hi;
					bestHubFace = -1;
					bestHubA = aId;
					bestHubB = bId;
				}
			}
		}
	}

	if (bestHubVertex >= 0 && bestHubVertexDistance <= 9.0) {
		vtkPolyData* hub = m_hubCages[bestHubVertex];
		if (!hub || !hub->GetPoints() ||
			bestHubPointId < 0 ||
			bestHubPointId >= hub->GetNumberOfPoints())
		{
			return false;
		}

		double p[3];
		hub->GetPoint(bestHubPointId, p);
		target.type = StraightCageEditTargetType::HubVertex;
		target.hubIndex = bestHubVertex;
		target.hubPointId = bestHubPointId;
		for (int d = 0; d < 3; ++d) {
			target.center[d] = p[d];
			target.originalA[d] = p[d];
		}
		return true;
	}

	if (bestHub >= 0 && bestHubDistance <= 8.0 && bestHubDistance <= bestEdgeDistance + 1.5) {
		vtkPolyData* hub = m_hubCages[bestHub];
		vtkCell* cell = (hub && bestHubFace >= 0) ? hub->GetCell(bestHubFace) : nullptr;
		if (!hub) {
			return false;
		}

		double a[3], b[3];
		hub->GetPoint(bestHubA, a);
		hub->GetPoint(bestHubB, b);

		target.type = StraightCageEditTargetType::HubEdge;
		target.hubIndex = bestHub;
		target.hubFaceIndex = bestHubFace;
		target.hubPointA = bestHubA;
		target.hubPointB = bestHubB;
		for (int d = 0; d < 3; ++d) {
			target.center[d] = 0.5 * (a[d] + b[d]);
			target.originalA[d] = a[d];
			target.originalB[d] = b[d];
		}

		double edgeDir[3] = {
			b[0] - a[0],
			b[1] - a[1],
			b[2] - a[2]
		};
		if (vtkMath::Normalize(edgeDir) == 0.0) {
			return false;
		}
		for (int d = 0; d < 3; ++d) {
			target.edgeDirection[d] = edgeDir[d];
		}

		double moveDir[3] = { 0.0, 0.0, 0.0 };
		double positiveLimit = 1e100;
		if (cell && cell->GetNumberOfPoints() >= 3) {
			double faceCenter[3] = { 0.0, 0.0, 0.0 };
			for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
				double p[3];
				hub->GetPoint(cell->GetPointId(k), p);
				for (int d = 0; d < 3; ++d) faceCenter[d] += p[d] / static_cast<double>(cell->GetNumberOfPoints());
			}

			double p0[3], p1[3], p2[3];
			hub->GetPoint(cell->GetPointId(0), p0);
			hub->GetPoint(cell->GetPointId(1), p1);
			hub->GetPoint(cell->GetPointId(2), p2);
			double e01[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
			double e02[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
			double faceNormal[3];
			vtkMath::Cross(e01, e02, faceNormal);
			if (vtkMath::Normalize(faceNormal) == 0.0) return false;
			vtkMath::Cross(faceNormal, edgeDir, moveDir);
			if (vtkMath::Normalize(moveDir) == 0.0) return false;

			double toFaceCenter[3] = {
				faceCenter[0] - target.center[0],
				faceCenter[1] - target.center[1],
				faceCenter[2] - target.center[2]
			};
			if (vtkMath::Dot(moveDir, toFaceCenter) < 0.0) {
				for (int d = 0; d < 3; ++d) moveDir[d] = -moveDir[d];
			}

			for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
				const vtkIdType id = cell->GetPointId(k);
				if (id == bestHubA || id == bestHubB) continue;
				double p[3];
				hub->GetPoint(id, p);
				double ap[3] = { p[0] - target.center[0], p[1] - target.center[1], p[2] - target.center[2] };
				const double dist = vtkMath::Dot(ap, moveDir);
				if (dist > 1e-8) positiveLimit = std::min(positiveLimit, dist);
			}
		}
		else {
			// 上下截面连接边沿“hub 中心 -> 边中心”的径向移动，
			// 并去掉沿边自身的分量，保持边长和原有上下对应关系。
			double hubCenter[3] = { 0.0, 0.0, 0.0 };
			for (vtkIdType i = 0; i < hub->GetNumberOfPoints(); ++i) {
				double p[3];
				hub->GetPoint(i, p);
				for (int d = 0; d < 3; ++d) hubCenter[d] += p[d] / static_cast<double>(hub->GetNumberOfPoints());
			}
			for (int d = 0; d < 3; ++d) moveDir[d] = target.center[d] - hubCenter[d];
			double alongEdge = vtkMath::Dot(moveDir, edgeDir);
			for (int d = 0; d < 3; ++d) moveDir[d] -= edgeDir[d] * alongEdge;
			if (vtkMath::Normalize(moveDir) == 0.0) return false;
			positiveLimit = std::max(
				std::sqrt(vtkMath::Distance2BetweenPoints(target.center, hubCenter)) * 1.5,
				std::sqrt(vtkMath::Distance2BetweenPoints(a, b)));
		}
		if (positiveLimit >= 1e90) {
			positiveLimit = std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
		}

		for (int d = 0; d < 3; ++d) {
			target.normal[d] = moveDir[d];
		}
		target.minOffset = -positiveLimit * 0.85;
		target.maxOffset = positiveLimit * 0.85;

		return true;
	}

	if (bestHubFaceHit >= 0 && bestHubFaceHitCell >= 0) {
		vtkPolyData* hub = m_hubCages[bestHubFaceHit];
		vtkCell* cell = hub ? hub->GetCell(bestHubFaceHitCell) : nullptr;
		if (!hub || !cell || cell->GetNumberOfPoints() < 3) {
			return false;
		}

		target.type = StraightCageEditTargetType::HubFace;
		target.hubIndex = bestHubFaceHit;
		target.hubFaceIndex = bestHubFaceHitCell;
		target.center[0] = target.center[1] = target.center[2] = 0.0;
		target.hubFacePointIds.reserve(
			static_cast<size_t>(cell->GetNumberOfPoints()));

		for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
			const vtkIdType pointId = cell->GetPointId(k);
			double p[3];
			hub->GetPoint(pointId, p);
			target.hubFacePointIds.push_back(pointId);
			for (int d = 0; d < 3; ++d) {
				target.center[d] +=
					p[d] / static_cast<double>(cell->GetNumberOfPoints());
			}
		}

		target.normal[0] = target.normal[1] = target.normal[2] = 0.0;
		for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
			double p[3], q[3];
			hub->GetPoint(cell->GetPointId(k), p);
			hub->GetPoint(
				cell->GetPointId((k + 1) % cell->GetNumberOfPoints()), q);
			target.normal[0] += (p[1] - q[1]) * (p[2] + q[2]);
			target.normal[1] += (p[2] - q[2]) * (p[0] + q[0]);
			target.normal[2] += (p[0] - q[0]) * (p[1] + q[1]);
		}
		if (vtkMath::Normalize(target.normal) == 0.0) {
			return false;
		}

		double hubCenter[3] = { 0.0, 0.0, 0.0 };
		for (vtkIdType i = 0; i < hub->GetNumberOfPoints(); ++i) {
			double p[3];
			hub->GetPoint(i, p);
			for (int d = 0; d < 3; ++d) {
				hubCenter[d] +=
					p[d] / static_cast<double>(hub->GetNumberOfPoints());
			}
		}
		double outward[3] = {
			target.center[0] - hubCenter[0],
			target.center[1] - hubCenter[1],
			target.center[2] - hubCenter[2]
		};
		if (vtkMath::Dot(target.normal, outward) < 0.0) {
			for (int d = 0; d < 3; ++d) {
				target.normal[d] = -target.normal[d];
			}
		}

		return true;
	}

	if (bestVertexCage >= 0 && bestVertexDistance <= 9.0) {
		const CageSegment& cage = m_allCages[bestVertexCage];
		double p[3];
		StraightCageEditGetPoint(cage, bestVertexId, p);
		target.type = StraightCageEditTargetType::Vertex;
		target.cageIndex = bestVertexCage;
		target.vertexIndex = bestVertexId;
		for (int d = 0; d < 3; ++d) {
			target.center[d] = p[d];
			target.originalA[d] = p[d];
		}
		return true;
	}

	if (bestCage >= 0 && bestEdgeDistance <= 10.0) {
		const CageSegment& cage = m_allCages[bestCage];
		double a[3], b[3], cageCenter[3];
		StraightCageEditGetPoint(cage, bestA, a);
		StraightCageEditGetPoint(cage, bestB, b);
		StraightCageEditCageCenter(cage, cageCenter);

		target.type = StraightCageEditTargetType::Edge;
		target.cageIndex = bestCage;
		target.edgeVertexA = bestA;
		target.edgeVertexB = bestB;
		for (int d = 0; d < 3; ++d) {
			target.center[d] = 0.5 * (a[d] + b[d]);
		}

		double edgeDir[3] = {
			b[0] - a[0],
			b[1] - a[1],
			b[2] - a[2]
		};
		if (vtkMath::Normalize(edgeDir) == 0.0) {
			return false;
		}
		for (int d = 0; d < 3; ++d) {
			target.edgeDirection[d] = edgeDir[d];
		}

		double moveDirs[2][3] = {
			{ 0.0, 0.0, 1.0 },
			{ 0.0, 0.0, 1.0 }
		};
		int moveDirCount = 0;

		for (int faceIndex = 0; faceIndex < 6; ++faceIndex) {
			int faceIds[4];
			StraightCageEditFaceIds(faceIndex, faceIds);

			bool hasA = false;
			bool hasB = false;
			for (int k = 0; k < 4; ++k) {
				hasA = hasA || (faceIds[k] == bestA);
				hasB = hasB || (faceIds[k] == bestB);
			}
			if (!hasA || !hasB) {
				continue;
			}

			double fp[4][3];
			double faceCenter[3] = { 0.0, 0.0, 0.0 };
			for (int k = 0; k < 4; ++k) {
				StraightCageEditGetPoint(cage, faceIds[k], fp[k]);
				for (int d = 0; d < 3; ++d) {
					faceCenter[d] += fp[k][d] * 0.25;
				}
			}

			double e01[3] = {
				fp[1][0] - fp[0][0],
				fp[1][1] - fp[0][1],
				fp[1][2] - fp[0][2]
			};
			double e02[3] = {
				fp[2][0] - fp[0][0],
				fp[2][1] - fp[0][1],
				fp[2][2] - fp[0][2]
			};
			double faceNormal[3];
			vtkMath::Cross(e01, e02, faceNormal);
			if (vtkMath::Normalize(faceNormal) == 0.0) {
				continue;
			}

			double moveDir[3];
			vtkMath::Cross(faceNormal, edgeDir, moveDir);
			if (vtkMath::Normalize(moveDir) == 0.0) {
				continue;
			}

			double outOnFace[3] = {
				target.center[0] - faceCenter[0],
				target.center[1] - faceCenter[1],
				target.center[2] - faceCenter[2]
			};
			if (vtkMath::Dot(moveDir, outOnFace) < 0.0) {
				for (int d = 0; d < 3; ++d) {
					moveDir[d] = -moveDir[d];
				}
			}

			bool duplicateDirection = false;
			for (int i = 0; i < moveDirCount; ++i) {
				if (std::abs(vtkMath::Dot(moveDirs[i], moveDir)) > 0.985) {
					duplicateDirection = true;
					break;
				}
			}
			if (!duplicateDirection && moveDirCount < 2) {
				for (int d = 0; d < 3; ++d) {
					moveDirs[moveDirCount][d] = moveDir[d];
				}
				++moveDirCount;
			}
		}

		if (moveDirCount <= 0) {
			double out[3] = {
				target.center[0] - cageCenter[0],
				target.center[1] - cageCenter[1],
				target.center[2] - cageCenter[2]
			};
			const double along = vtkMath::Dot(out, edgeDir);
			for (int d = 0; d < 3; ++d) {
				out[d] -= along * edgeDir[d];
			}
			if (vtkMath::Normalize(out) != 0.0) {
				for (int d = 0; d < 3; ++d) {
					moveDirs[0][d] = out[d];
				}
				moveDirCount = 1;
			}
		}

		if (moveDirCount <= 0) {
			return false;
		}

		for (int d = 0; d < 3; ++d) {
			target.normal[d] = moveDirs[0][d];
			target.alternateNormal[d] = moveDirs[moveDirCount > 1 ? 1 : 0][d];
		}
		target.hasAlternateNormal = moveDirCount > 1;

		return true;
	}

	if (!m_CageSolidActor) {
		return false;
	}

	vtkNew<vtkCellPicker> picker;
	picker->SetTolerance(0.005);
	m_CageSolidActor->SetPickable(1);
	m_CageSolidActor->SetVisibility(true);
	picker->AddPickList(m_CageSolidActor);
	picker->PickFromListOn();
	picker->Pick(x, y, 0, ren);
	const bool picked = picker->GetActor() == m_CageSolidActor && picker->GetCellId() >= 0;
	vtkIdType cellId = picker->GetCellId();
	m_CageSolidActor->SetVisibility(false);
	m_CageSolidActor->SetPickable(0);

	if (!picked) {
		return false;
	}

	const int cageIndex = static_cast<int>(cellId / 6);
	const int faceIndex = static_cast<int>(cellId % 6);
	if (cageIndex < 0 || cageIndex >= static_cast<int>(m_allCages.size())) {
		return false;
	}

	const CageSegment& cage = m_allCages[cageIndex];
	int ids[4];
	StraightCageEditFaceIds(faceIndex, ids);

	double p[4][3];
	for (int i = 0; i < 4; ++i) {
		StraightCageEditGetPoint(cage, ids[i], p[i]);
	}

	target.type = StraightCageEditTargetType::Face;
	target.cageIndex = cageIndex;
	target.faceIndex = faceIndex;
	target.center[0] = target.center[1] = target.center[2] = 0.0;
	for (int i = 0; i < 4; ++i) {
		target.center[0] += p[i][0] * 0.25;
		target.center[1] += p[i][1] * 0.25;
		target.center[2] += p[i][2] * 0.25;
	}

	double e01[3] = {
		p[1][0] - p[0][0],
		p[1][1] - p[0][1],
		p[1][2] - p[0][2]
	};
	double e02[3] = {
		p[2][0] - p[0][0],
		p[2][1] - p[0][1],
		p[2][2] - p[0][2]
	};
	vtkMath::Cross(e01, e02, target.normal);

	if (vtkMath::Normalize(target.normal) == 0.0) {
		return false;
	}

	double cageCenter[3];
	StraightCageEditCageCenter(cage, cageCenter);
	double outDir[3] = {
		target.center[0] - cageCenter[0],
		target.center[1] - cageCenter[1],
		target.center[2] - cageCenter[2]
	};
	if (vtkMath::Dot(target.normal, outDir) < 0.0) {
		for (int d = 0; d < 3; ++d) {
			target.normal[d] = -target.normal[d];
		}
	}

	return true;
}

bool CageInteractorStyle::HandleStraightCageEditLeftClick(int x, int y, vtkRenderer* ren)
{
	if (!m_isStraightCageEditMode || !ren) {
		return false;
	}

	if (m_straightCageEditActive) {
		ApplyStraightCageEdit();
		CancelStraightCageEdit(ren);
		if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
			this->GetInteractor()->GetRenderWindow()->Render();
		}
		return true;
	}

	StraightCageEditTarget target;
	if (!PickStraightCageEditTarget(x, y, ren, target)) {
		return false;
	}

	m_straightCageEditTarget = target;
	m_straightCageEditActive = true;
	m_straightCageEditStartMouse[0] = x;
	m_straightCageEditStartMouse[1] = y;
	m_straightCageEditCurrentOffset = 0.0;
	m_straightCageEditCurrentDelta[0] = 0.0;
	m_straightCageEditCurrentDelta[1] = 0.0;
	m_straightCageEditCurrentDelta[2] = 0.0;

	m_straightCageEditReferenceStep = 1e-3;
	if (target.type == StraightCageEditTargetType::HubVertex) {
		m_straightCageEditReferenceStep = std::max(1e-3, 1.0);
	}
	else if (target.type == StraightCageEditTargetType::HubEdge) {
		m_straightCageEditReferenceStep =
			std::max(
				1e-3,
				std::sqrt(vtkMath::Distance2BetweenPoints(target.originalA, target.originalB)) * 0.50);
	}
	else if (target.type == StraightCageEditTargetType::HubFace) {
		vtkPolyData* hub = m_hubCages[target.hubIndex];
		for (vtkIdType pointId : target.hubFacePointIds) {
			double p[3];
			hub->GetPoint(pointId, p);
			m_straightCageEditReferenceStep =
				std::max(
					m_straightCageEditReferenceStep,
					std::sqrt(vtkMath::Distance2BetweenPoints(p, target.center)) * 0.35);
		}
	}
	else {
		const CageSegment& cage = m_allCages[target.cageIndex];
		double cageCenter[3];
		StraightCageEditCageCenter(cage, cageCenter);
		for (int i = 0; i < 8; ++i) {
			double p[3];
			StraightCageEditGetPoint(cage, i, p);
			m_straightCageEditReferenceStep =
				std::max(
					m_straightCageEditReferenceStep,
					std::sqrt(vtkMath::Distance2BetweenPoints(p, cageCenter)) * 0.35);
		}
	}

	UpdateStraightCageEditActors(ren, true);
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

bool CageInteractorStyle::UpdateStraightCageEditPreview(int x, int y, vtkRenderer* ren)
{
	if (!m_straightCageEditActive || !ren) {
		return false;
	}

	const double* c = m_straightCageEditTarget.center;

	ren->SetWorldPoint(c[0], c[1], c[2], 1.0);
	ren->WorldToDisplay();
	double* dc = ren->GetDisplayPoint();
	const double cx = dc[0];
	const double cy = dc[1];

	const double dx = static_cast<double>(x - m_straightCageEditStartMouse[0]);
	const double dy = static_cast<double>(y - m_straightCageEditStartMouse[1]);

	if (m_straightCageEditTarget.type == StraightCageEditTargetType::Vertex ||
		m_straightCageEditTarget.type == StraightCageEditTargetType::Edge ||
		m_straightCageEditTarget.type == StraightCageEditTargetType::HubVertex ||
		m_straightCageEditTarget.type == StraightCageEditTargetType::HubEdge)
	{
		auto worldPointAtMouseOnCenterPlane =
			[&](int px, int py, double out[3]) -> bool
		{
			ren->SetDisplayPoint(static_cast<double>(px), static_cast<double>(py), dc[2]);
			ren->DisplayToWorld();
			double* world = ren->GetWorldPoint();
			if (!world || std::abs(world[3]) < 1e-12) {
				return false;
			}
			out[0] = world[0] / world[3];
			out[1] = world[1] / world[3];
			out[2] = world[2] / world[3];
			return true;
		};

		double startWorld[3];
		double currentWorld[3];
		if (!worldPointAtMouseOnCenterPlane(
			m_straightCageEditStartMouse[0],
			m_straightCageEditStartMouse[1],
			startWorld) ||
			!worldPointAtMouseOnCenterPlane(x, y, currentWorld))
		{
			return false;
		}

		double delta[3] = {
			currentWorld[0] - startWorld[0],
			currentWorld[1] - startWorld[1],
			currentWorld[2] - startWorld[2]
		};

		double deltaLen = vtkMath::Norm(delta);
		if (m_straightCageEditTarget.type == StraightCageEditTargetType::HubEdge) {
			const double maxLen = std::max(
				std::abs(m_straightCageEditTarget.minOffset),
				std::abs(m_straightCageEditTarget.maxOffset));
			if (std::isfinite(maxLen) && maxLen > 1e-8 && deltaLen > maxLen) {
				const double scale = maxLen / deltaLen;
				for (int d = 0; d < 3; ++d) {
					delta[d] *= scale;
				}
				deltaLen = maxLen;
			}
		}

		for (int d = 0; d < 3; ++d) {
			m_straightCageEditCurrentDelta[d] = delta[d];
		}
		m_straightCageEditCurrentOffset = deltaLen;

		UpdateStraightCageEditActors(ren, true);
		if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
			this->GetInteractor()->GetRenderWindow()->Render();
		}
		return true;
	}

	auto projectOffsetForDirection =
		[&](const double dir[3], double& offset, double& score)
	{
		double p1[3] = {
			c[0] + dir[0] * m_straightCageEditReferenceStep,
			c[1] + dir[1] * m_straightCageEditReferenceStep,
			c[2] + dir[2] * m_straightCageEditReferenceStep
		};
		ren->SetWorldPoint(p1[0], p1[1], p1[2], 1.0);
		ren->WorldToDisplay();
		double* dn = ren->GetDisplayPoint();

		double sx = dn[0] - cx;
		double sy = dn[1] - cy;
		double len2 = sx * sx + sy * sy;
		if (len2 <= 1e-8) {
			sx = 0.0;
			sy = 1.0;
			len2 = 1.0;
		}

		const double projectedPixels = (dx * sx + dy * sy) / len2;
		offset = projectedPixels * m_straightCageEditReferenceStep;
		score = std::abs(dx * sx + dy * sy) / std::sqrt(len2);
	};

	double selectedOffset = 0.0;
	double selectedScore = 0.0;
	projectOffsetForDirection(
		m_straightCageEditTarget.normal,
		selectedOffset,
		selectedScore);

	if (m_straightCageEditTarget.type == StraightCageEditTargetType::Edge &&
		m_straightCageEditTarget.hasAlternateNormal)
	{
		double alternateOffset = 0.0;
		double alternateScore = 0.0;
		projectOffsetForDirection(
			m_straightCageEditTarget.alternateNormal,
			alternateOffset,
			alternateScore);

		if (alternateScore > selectedScore + 1e-6) {
			for (int d = 0; d < 3; ++d) {
				std::swap(
					m_straightCageEditTarget.normal[d],
					m_straightCageEditTarget.alternateNormal[d]);
			}
			selectedOffset = alternateOffset;
		}
	}

	m_straightCageEditCurrentOffset =
		selectedOffset;
	m_straightCageEditCurrentOffset =
		std::max(
			m_straightCageEditTarget.minOffset,
			std::min(
				m_straightCageEditTarget.maxOffset,
				m_straightCageEditCurrentOffset));

	for (int d = 0; d < 3; ++d) {
		m_straightCageEditCurrentDelta[d] =
			m_straightCageEditTarget.normal[d] *
			m_straightCageEditCurrentOffset;
	}

	UpdateStraightCageEditActors(ren, true);
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

void CageInteractorStyle::ApplyStraightCageEdit()
{
	const StraightCageEditTarget& target = m_straightCageEditTarget;
	if (target.type == StraightCageEditTargetType::HubVertex) {
		if (target.hubIndex < 0 ||
			target.hubIndex >= static_cast<int>(m_hubCages.size()) ||
			!m_hubCages[target.hubIndex] ||
			!m_hubCages[target.hubIndex]->GetPoints() ||
			target.hubPointId < 0 ||
			target.hubPointId >= m_hubCages[target.hubIndex]->GetNumberOfPoints() ||
			vtkMath::Norm(m_straightCageEditCurrentDelta) <= 1e-9)
		{
			return;
		}

		double moved[3] = {
			target.originalA[0] + m_straightCageEditCurrentDelta[0],
			target.originalA[1] + m_straightCageEditCurrentDelta[1],
			target.originalA[2] + m_straightCageEditCurrentDelta[2]
		};

		const double tol2 = 1e-6;
		for (CageSegment& cage : m_allCages) {
			for (int i = 0; i < 8; ++i) {
				double p[3];
				StraightCageEditGetPoint(cage, i, p);
				if (vtkMath::Distance2BetweenPoints(p, target.originalA) <= tol2) {
					StraightCageEditSetPoint(cage, i, moved);
				}
			}
		}

		for (vtkSmartPointer<vtkPolyData>& hub : m_hubCages) {
			if (!hub || !hub->GetPoints()) {
				continue;
			}
			vtkPoints* points = hub->GetPoints();
			for (vtkIdType pi = 0; pi < points->GetNumberOfPoints(); ++pi) {
				double p[3];
				points->GetPoint(pi, p);
				if (vtkMath::Distance2BetweenPoints(p, target.originalA) <= tol2) {
					points->SetPoint(pi, moved);
				}
			}
			points->Modified();
			hub->Modified();
		}

		RebuildGlobalMesh();
		return;
	}

	if (target.type == StraightCageEditTargetType::HubFace) {
		if (target.hubIndex < 0 ||
			target.hubIndex >= static_cast<int>(m_hubCages.size()) ||
			!m_hubCages[target.hubIndex] ||
			!m_hubCages[target.hubIndex]->GetPoints() ||
			target.hubFacePointIds.empty() ||
			std::abs(m_straightCageEditCurrentOffset) <= 1e-9)
		{
			return;
		}

		struct HubPointMove {
			vtkIdType pointId = -1;
			double from[3] = { 0.0, 0.0, 0.0 };
			double to[3] = { 0.0, 0.0, 0.0 };
		};

		vtkPolyData* sourceHub = m_hubCages[target.hubIndex];
		std::vector<HubPointMove> moves;
		moves.reserve(target.hubFacePointIds.size());
		for (vtkIdType pointId : target.hubFacePointIds) {
			if (pointId < 0 || pointId >= sourceHub->GetNumberOfPoints()) {
				continue;
			}
			HubPointMove move;
			move.pointId = pointId;
			sourceHub->GetPoint(pointId, move.from);
			for (int d = 0; d < 3; ++d) {
				move.to[d] =
					move.from[d] +
					m_straightCageEditCurrentDelta[d];
			}
			moves.push_back(move);
		}
		if (moves.empty()) {
			return;
		}

		for (const HubPointMove& move : moves) {
			sourceHub->GetPoints()->SetPoint(move.pointId, move.to);
		}
		sourceHub->GetPoints()->Modified();
		sourceHub->Modified();

		const double tol2 = 1e-6;
		for (CageSegment& cage : m_allCages) {
			for (int i = 0; i < 8; ++i) {
				double p[3];
				StraightCageEditGetPoint(cage, i, p);
				for (const HubPointMove& move : moves) {
					if (vtkMath::Distance2BetweenPoints(p, move.from) <= tol2) {
						StraightCageEditSetPoint(cage, i, move.to);
						break;
					}
				}
			}
		}

		for (int hi = 0; hi < static_cast<int>(m_hubCages.size()); ++hi) {
			if (hi == target.hubIndex ||
				!m_hubCages[hi] ||
				!m_hubCages[hi]->GetPoints())
			{
				continue;
			}
			vtkPoints* points = m_hubCages[hi]->GetPoints();
			for (vtkIdType pi = 0; pi < points->GetNumberOfPoints(); ++pi) {
				double p[3];
				points->GetPoint(pi, p);
				for (const HubPointMove& move : moves) {
					if (vtkMath::Distance2BetweenPoints(p, move.from) <= tol2) {
						points->SetPoint(pi, move.to);
						break;
					}
				}
			}
			points->Modified();
			m_hubCages[hi]->Modified();
		}

		RebuildGlobalMesh();
		return;
	}

	if (target.type == StraightCageEditTargetType::HubEdge) {
		if (target.hubIndex < 0 ||
			target.hubIndex >= static_cast<int>(m_hubCages.size()) ||
			!m_hubCages[target.hubIndex] ||
			!m_hubCages[target.hubIndex]->GetPoints())
		{
			return;
		}

		vtkPolyData* hub = m_hubCages[target.hubIndex];
		double movedA[3];
		double movedB[3];
		for (int d = 0; d < 3; ++d) {
			movedA[d] = target.originalA[d] + m_straightCageEditCurrentDelta[d];
			movedB[d] = target.originalB[d] + m_straightCageEditCurrentDelta[d];
		}

		auto MoveMatchingCagePoint = [&](const double originalPoint[3], const double movedPoint[3])
		{
			const double tol2 = 1e-6;
			for (CageSegment& cage : m_allCages) {
				for (int i = 0; i < 8; ++i) {
					double p[3];
					StraightCageEditGetPoint(cage, i, p);
					if (vtkMath::Distance2BetweenPoints(p, originalPoint) <= tol2) {
						StraightCageEditSetPoint(cage, i, movedPoint);
					}
				}
			}
		};

		hub->GetPoints()->SetPoint(target.hubPointA, movedA);
		hub->GetPoints()->SetPoint(target.hubPointB, movedB);
		hub->GetPoints()->Modified();
		hub->Modified();

		MoveMatchingCagePoint(target.originalA, movedA);
		MoveMatchingCagePoint(target.originalB, movedB);

		RebuildGlobalMesh();
		return;
	}

	if (target.cageIndex < 0 || target.cageIndex >= static_cast<int>(m_allCages.size())) {
		return;
	}

	int ids[4] = { -1, -1, -1, -1 };
	int count = 0;
	if (target.type == StraightCageEditTargetType::Vertex) {
		ids[0] = target.vertexIndex;
		count = 1;
	}
	else if (target.type == StraightCageEditTargetType::Face) {
		StraightCageEditFaceIds(target.faceIndex, ids);
		count = 4;
	}
	else if (target.type == StraightCageEditTargetType::Edge) {
		ids[0] = target.edgeVertexA;
		ids[1] = target.edgeVertexB;
		count = 2;
	}

	if (count <= 0) {
		return;
	}

	if (vtkMath::Norm(m_straightCageEditCurrentDelta) <= 1e-9) {
		return;
	}

	struct PointMove {
		double from[3];
		double to[3];
	};

	std::vector<PointMove> moves;
	moves.reserve(count);

	const CageSegment& sourceCage = m_allCages[target.cageIndex];

	for (int i = 0; i < count; ++i) {
		PointMove move;
		StraightCageEditGetPoint(sourceCage, ids[i], move.from);
		for (int d = 0; d < 3; ++d) {
			move.to[d] = move.from[d] + m_straightCageEditCurrentDelta[d];
		}

		moves.push_back(move);
	}

	auto MoveLocalCagePointIfMatches = [&](CageSegment& cage, const PointMove& move, double tol2)
	{
		for (int i = 0; i < 8; ++i) {
			double p[3];
			StraightCageEditGetPoint(cage, i, p);
			if (vtkMath::Distance2BetweenPoints(p, move.from) <= tol2) {
				StraightCageEditSetPoint(cage, i, move.to);
			}
		}
	};

	const double tol2 = 1e-6;

	if (target.type == StraightCageEditTargetType::Face ||
		target.type == StraightCageEditTargetType::Vertex) {
		// 平面编辑只平移所选面的四个逻辑节点。相邻 cage 中坐标相同的点
		// 是同一逻辑节点的存储副本，需要同步移动以保持共享边连续。
		for (CageSegment& cage : m_allCages) {
			for (const PointMove& move : moves) {
				MoveLocalCagePointIfMatches(cage, move, tol2);
			}
		}

		// hub 中若共享了同一节点，也只同步这些匹配点；不移动整个 hub，
		// 更不会在原面和新面之间创建额外控制笼段。
		for (vtkSmartPointer<vtkPolyData>& hub : m_hubCages) {
			if (!hub || !hub->GetPoints()) {
				continue;
			}

			vtkPoints* points = hub->GetPoints();
			const vtkIdType pointCount = points->GetNumberOfPoints();

			for (vtkIdType pi = 0; pi < pointCount; ++pi) {
				double p[3];
				points->GetPoint(pi, p);
				for (const PointMove& move : moves) {
					if (vtkMath::Distance2BetweenPoints(p, move.from) <= tol2) {
						points->SetPoint(pi, move.to);
						break;
					}
				}
			}

			points->Modified();
			hub->Modified();
		}
	}
	else {
		for (const PointMove& move : moves) {
			for (CageSegment& cage : m_allCages) {
				MoveLocalCagePointIfMatches(cage, move, tol2);
			}
		}
	}

	RebuildGlobalMesh();
}

bool CageInteractorStyle::RotateStraightCageEditTarget(double angleRadians)
{
	const StraightCageEditTarget& target = m_straightCageEditTarget;
	if (target.type == StraightCageEditTargetType::HubEdge ||
		target.type == StraightCageEditTargetType::HubFace)
	{
		return false;
	}
	if (target.cageIndex < 0 || target.cageIndex >= static_cast<int>(m_allCages.size())) {
		return false;
	}

	CageSegment& cage = m_allCages[target.cageIndex];
	double startCenter[3], endCenter[3], axis[3];
	StraightCageEditCentersAxis(cage, startCenter, endCenter, axis);
	if (vtkMath::Norm(axis) <= 1e-8) {
		return false;
	}

	for (int i = 0; i < 4; ++i) {
		StraightCageRotatePoint(cage.pStart[i], startCenter, axis, angleRadians);
		StraightCageRotatePoint(cage.pEnd[i], endCenter, axis, angleRadians);
	}

	CancelStraightCageEdit(nullptr);
	RebuildGlobalMesh();
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

bool CageInteractorStyle::HandleStraightCageRotateLeftClick(int x, int y, vtkRenderer* ren)
{
	if (!m_isStraightCageRotateMode || !ren) {
		return false;
	}

	if (m_straightCageRotateActive) {
		ApplyStraightCageRotation();
		CancelStraightCageRotation(ren);
		return true;
	}

	StraightCageEditTarget target;
	if (!PickStraightCageEditTarget(x, y, ren, target)) {
		return false;
	}
	if (target.type == StraightCageEditTargetType::HubEdge ||
		target.type == StraightCageEditTargetType::HubFace ||
		target.cageIndex < 0 ||
		target.cageIndex >= static_cast<int>(m_allCages.size()))
	{
		return false;
	}

	const CageSegment& cage = m_allCages[target.cageIndex];
	double center[3];
	StraightCageEditCageCenter(cage, center);
	ren->SetWorldPoint(center[0], center[1], center[2], 1.0);
	ren->WorldToDisplay();
	double* display = ren->GetDisplayPoint();
	if (!display) {
		return false;
	}

	m_straightCageRotateTarget = target;
	m_straightCageRotateActive = true;
	m_straightCageRotateStartMouse[0] = x;
	m_straightCageRotateStartMouse[1] = y;
	m_straightCageRotateStartDisplayAngle =
		std::atan2(static_cast<double>(y) - display[1], static_cast<double>(x) - display[0]);
	m_straightCageRotateCurrentAngle = 0.0;

	UpdateStraightCageRotateActors(ren, true);
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

bool CageInteractorStyle::UpdateStraightCageRotatePreview(int x, int y, vtkRenderer* ren)
{
	if (!m_straightCageRotateActive || !ren) {
		return false;
	}
	if (m_straightCageRotateTarget.cageIndex < 0 ||
		m_straightCageRotateTarget.cageIndex >= static_cast<int>(m_allCages.size()))
	{
		return false;
	}

	const CageSegment& cage = m_allCages[m_straightCageRotateTarget.cageIndex];
	double center[3];
	StraightCageEditCageCenter(cage, center);
	ren->SetWorldPoint(center[0], center[1], center[2], 1.0);
	ren->WorldToDisplay();
	double* display = ren->GetDisplayPoint();
	if (!display) {
		return false;
	}

	const double dx = static_cast<double>(x) - display[0];
	const double dy = static_cast<double>(y) - display[1];
	if (dx * dx + dy * dy > 4.0) {
		double angle =
			std::atan2(dy, dx) - m_straightCageRotateStartDisplayAngle;
		const double twoPi = 2.0 * 3.14159265358979323846;
		while (angle > 3.14159265358979323846) angle -= twoPi;
		while (angle < -3.14159265358979323846) angle += twoPi;
		m_straightCageRotateCurrentAngle = angle;
	}

	UpdateStraightCageRotateActors(ren, true);
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

void CageInteractorStyle::ApplyStraightCageRotation()
{
	if (m_straightCageRotateTarget.cageIndex < 0 ||
		m_straightCageRotateTarget.cageIndex >= static_cast<int>(m_allCages.size()))
	{
		return;
	}

	StraightCageEditTarget previous = m_straightCageEditTarget;
	m_straightCageEditTarget = m_straightCageRotateTarget;
	RotateStraightCageEditTarget(m_straightCageRotateCurrentAngle);
	m_straightCageEditTarget = previous;
}

void CageInteractorStyle::CancelStraightCageRotation(vtkRenderer* ren)
{
	m_straightCageRotateActive = false;
	m_straightCageRotateTarget = StraightCageEditTarget();
	m_straightCageRotateCurrentAngle = 0.0;
	if (m_straightCageEditHighlightActor) {
		m_straightCageEditHighlightActor->VisibilityOff();
	}
	if (m_straightCageEditPreviewActor) {
		m_straightCageEditPreviewActor->VisibilityOff();
	}
	if (ren && this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::UpdateStraightCageRotateActors(vtkRenderer* ren, bool showPreview)
{
	if (!ren ||
		m_straightCageRotateTarget.cageIndex < 0 ||
		m_straightCageRotateTarget.cageIndex >= static_cast<int>(m_allCages.size()))
	{
		return;
	}

	auto ensureActor =
		[&](vtkSmartPointer<vtkActor>& actor, double r, double g, double b, double lineWidth)
	{
		if (!actor) {
			vtkNew<vtkPolyDataMapper> mapper;
			actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->SetPickable(0);
			actor->GetProperty()->SetLighting(false);
			ren->AddActor(actor);
		}
		actor->GetProperty()->SetColor(r, g, b);
		actor->GetProperty()->SetRepresentationToWireframe();
		actor->GetProperty()->SetLineWidth(lineWidth);
		actor->GetProperty()->SetPointSize(lineWidth * 2.0);
	};

	ensureActor(m_straightCageEditHighlightActor, 0.0, 1.0, 0.0, 5.0);
	ensureActor(m_straightCageEditPreviewActor, 1.0, 1.0, 0.0, 4.0);

	const CageSegment& cage = m_allCages[m_straightCageRotateTarget.cageIndex];
	vtkPolyDataMapper::SafeDownCast(
		m_straightCageEditHighlightActor->GetMapper())->SetInputData(
			StraightCageBuildWholeData(cage));
	m_straightCageEditHighlightActor->VisibilityOn();

	if (showPreview) {
		const CageSegment rotated =
			StraightCageRotatedCopy(cage, m_straightCageRotateCurrentAngle);
		vtkPolyDataMapper::SafeDownCast(
			m_straightCageEditPreviewActor->GetMapper())->SetInputData(
				StraightCageBuildWholeData(rotated));
		m_straightCageEditPreviewActor->VisibilityOn();
	}
	else if (m_straightCageEditPreviewActor) {
		m_straightCageEditPreviewActor->VisibilityOff();
	}
}

void CageInteractorStyle::CancelStraightCageEdit(vtkRenderer* ren)
{
	m_straightCageEditActive = false;
	m_straightCageEditTarget = StraightCageEditTarget();
	m_straightCageEditCurrentOffset = 0.0;
	m_straightCageEditCurrentDelta[0] = 0.0;
	m_straightCageEditCurrentDelta[1] = 0.0;
	m_straightCageEditCurrentDelta[2] = 0.0;
	if (m_straightCageEditHighlightActor) {
		m_straightCageEditHighlightActor->VisibilityOff();
	}
	if (m_straightCageEditPreviewActor) {
		m_straightCageEditPreviewActor->VisibilityOff();
	}
	if (ren && this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::UpdateStraightCageEditActors(vtkRenderer* ren, bool showPreview)
{
	if (!ren) {
		return;
	}

	if (m_straightCageEditTarget.type == StraightCageEditTargetType::HubVertex ||
		m_straightCageEditTarget.type == StraightCageEditTargetType::HubEdge ||
		m_straightCageEditTarget.type == StraightCageEditTargetType::HubFace)
	{
		if (m_straightCageEditTarget.hubIndex < 0 ||
			m_straightCageEditTarget.hubIndex >= static_cast<int>(m_hubCages.size()))
		{
			return;
		}
	}
	else if (m_straightCageEditTarget.cageIndex < 0 ||
		m_straightCageEditTarget.cageIndex >= static_cast<int>(m_allCages.size()))
	{
		return;
	}

	auto ensureActor =
		[&](vtkSmartPointer<vtkActor>& actor, double r, double g, double b, double lineWidth)
	{
		if (!actor) {
			vtkNew<vtkPolyDataMapper> mapper;
			actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->SetPickable(0);
			actor->GetProperty()->SetLighting(false);
			ren->AddActor(actor);
		}
		actor->GetProperty()->SetColor(r, g, b);
		actor->GetProperty()->SetRepresentationToWireframe();
		actor->GetProperty()->SetLineWidth(lineWidth);
		actor->GetProperty()->SetPointSize(lineWidth * 2.0);
	};

	ensureActor(m_straightCageEditHighlightActor, 0.0, 1.0, 0.0, 5.0);
	ensureActor(m_straightCageEditPreviewActor, 1.0, 1.0, 0.0, 4.0);

	auto buildData =
		[&](const double delta[3]) -> vtkSmartPointer<vtkPolyData>
	{
		vtkNew<vtkPoints> points;
		vtkNew<vtkCellArray> verts;
		vtkNew<vtkCellArray> polys;
		vtkNew<vtkCellArray> lines;

		if (m_straightCageEditTarget.type == StraightCageEditTargetType::HubVertex) {
			double p[3];
			for (int d = 0; d < 3; ++d) {
				p[d] = m_straightCageEditTarget.originalA[d] + delta[d];
			}
			const vtkIdType id = points->InsertNextPoint(p);
			verts->InsertNextCell(1, &id);
		}
		else if (m_straightCageEditTarget.type == StraightCageEditTargetType::HubEdge) {
			double a[3];
			double b[3];
			for (int d = 0; d < 3; ++d) {
				a[d] = m_straightCageEditTarget.originalA[d] +
					delta[d];
				b[d] = m_straightCageEditTarget.originalB[d] +
					delta[d];
			}
			points->InsertNextPoint(a);
			points->InsertNextPoint(b);
			vtkIdType ids[2] = { 0, 1 };
			lines->InsertNextCell(2, ids);
		}
		else if (m_straightCageEditTarget.type == StraightCageEditTargetType::HubFace) {
			vtkPolyData* hub = m_hubCages[m_straightCageEditTarget.hubIndex];
			const std::vector<vtkIdType>& pointIds =
				m_straightCageEditTarget.hubFacePointIds;
			for (vtkIdType pointId : pointIds) {
				double p[3];
				hub->GetPoint(pointId, p);
				for (int d = 0; d < 3; ++d) {
					p[d] += delta[d];
				}
				points->InsertNextPoint(p);
			}
			polys->InsertNextCell(static_cast<vtkIdType>(pointIds.size()));
			for (vtkIdType i = 0; i < static_cast<vtkIdType>(pointIds.size()); ++i) {
				polys->InsertCellPoint(i);
			}
		}
		else if (m_straightCageEditTarget.type == StraightCageEditTargetType::Face) {
			const CageSegment& cage = m_allCages[m_straightCageEditTarget.cageIndex];
			int ids[4];
			StraightCageEditFaceIds(m_straightCageEditTarget.faceIndex, ids);
			for (int i = 0; i < 4; ++i) {
				double p[3];
				StraightCageEditGetPoint(cage, ids[i], p);
				for (int d = 0; d < 3; ++d) {
					p[d] += delta[d];
				}
				points->InsertNextPoint(p);
			}
			polys->InsertNextCell(4);
			for (vtkIdType i = 0; i < 4; ++i) {
				polys->InsertCellPoint(i);
			}
		}
		else if (m_straightCageEditTarget.type == StraightCageEditTargetType::Vertex) {
			const CageSegment& cage = m_allCages[m_straightCageEditTarget.cageIndex];
			double p[3];
			StraightCageEditGetPoint(cage, m_straightCageEditTarget.vertexIndex, p);
			for (int d = 0; d < 3; ++d) {
				p[d] += delta[d];
			}
			const vtkIdType id = points->InsertNextPoint(p);
			verts->InsertNextCell(1, &id);
		}
		else if (m_straightCageEditTarget.type == StraightCageEditTargetType::Edge) {
			const CageSegment& cage = m_allCages[m_straightCageEditTarget.cageIndex];
			for (int i = 0; i < 2; ++i) {
				const int id =
					i == 0
					? m_straightCageEditTarget.edgeVertexA
					: m_straightCageEditTarget.edgeVertexB;
				double p[3];
				StraightCageEditGetPoint(cage, id, p);
				for (int d = 0; d < 3; ++d) {
					p[d] += delta[d];
				}
				points->InsertNextPoint(p);
			}
			vtkIdType ids[2] = { 0, 1 };
			lines->InsertNextCell(2, ids);
		}

		vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
		data->SetPoints(points);
		data->SetVerts(verts);
		data->SetPolys(polys);
		data->SetLines(lines);
		return data;
	};

	double zeroDelta[3] = { 0.0, 0.0, 0.0 };
	vtkPolyDataMapper::SafeDownCast(
		m_straightCageEditHighlightActor->GetMapper())->SetInputData(buildData(zeroDelta));
	m_straightCageEditHighlightActor->VisibilityOn();

	if (showPreview) {
		vtkPolyDataMapper::SafeDownCast(
			m_straightCageEditPreviewActor->GetMapper())->SetInputData(
				buildData(m_straightCageEditCurrentDelta));
		m_straightCageEditPreviewActor->VisibilityOn();
	}
	else {
		m_straightCageEditPreviewActor->VisibilityOff();
	}
}

std::vector<int> CageInteractorStyle::FindPathBFS(int startIdx, int endIdx, const std::vector<std::vector<int>>& adjList) {
	std::vector<int> parent(adjList.size(), -1);
	std::vector<bool> visited(adjList.size(), false);
	std::queue<int> q;

	q.push(startIdx);
	visited[startIdx] = true;

	while (!q.empty()) {
		int curr = q.front();
		q.pop();

		if (curr == endIdx) break;

		for (int neighbor : adjList[curr]) {
			if (!visited[neighbor]) {
				visited[neighbor] = true;
				parent[neighbor] = curr;
				q.push(neighbor);
			}
		}
	}

	std::vector<int> path;
	if (!visited[endIdx]) return path;

	int curr = endIdx;
	while (curr != -1) {
		path.push_back(curr);
		curr = parent[curr];
	}
	std::reverse(path.begin(), path.end());
	return path;
}
std::vector<std::vector<int>> CageInteractorStyle::BuildSkeletonAdjacency() const
{
	int numNodes = static_cast<int>(m_skelNodes.size() / 3);
	std::vector<std::vector<int>> adjList(numNodes);

	for (size_t i = 0; i < m_skelEdges.size() / 2; ++i) {
		int u = static_cast<int>(m_skelEdges[2 * i]);
		int v = static_cast<int>(m_skelEdges[2 * i + 1]);

		if (u < 0 || v < 0 || u >= numNodes || v >= numNodes || u == v) {
			continue;
		}

		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}

	return adjList;
}
std::vector<int> CageInteractorStyle::FindPathBFSIgnoringEdge(
	int startIdx,
	int endIdx,
	const std::vector<std::vector<int>>& adjList,
	int blockedU,
	int blockedV)
{
	if (startIdx < 0 || endIdx < 0 ||
		startIdx >= static_cast<int>(adjList.size()) ||
		endIdx >= static_cast<int>(adjList.size()))
	{
		return std::vector<int>();
	}

	std::vector<int> parent(adjList.size(), -1);
	std::vector<bool> visited(adjList.size(), false);
	std::queue<int> q;

	q.push(startIdx);
	visited[startIdx] = true;

	while (!q.empty()) {
		int curr = q.front();
		q.pop();

		if (curr == endIdx) break;

		for (int neighbor : adjList[curr]) {
			bool isBlocked =
				(curr == blockedU && neighbor == blockedV) ||
				(curr == blockedV && neighbor == blockedU);

			if (isBlocked || visited[neighbor]) {
				continue;
			}

			visited[neighbor] = true;
			parent[neighbor] = curr;
			q.push(neighbor);
		}
	}

	if (!visited[endIdx]) {
		return std::vector<int>();
	}

	std::vector<int> path;
	int curr = endIdx;

	while (curr != -1) {
		path.push_back(curr);
		curr = parent[curr];
	}

	std::reverse(path.begin(), path.end());
	return path;
}

int CageInteractorStyle::FindNearestSkeletonNodeToSelectedFace(const SelectedFace& sel,
	double outCenter[3],
	double outNormal[3]) const
{
	if (m_skelNodes.empty()) return -1;
	if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) return -1;

	double startCenter[3] = { 0.0, 0.0, 0.0 };
	double endCenter[3] = { 0.0, 0.0, 0.0 };
	const CageSegment& cage = m_allCages[sel.cageIndex];

	for (int i = 0; i < 4; ++i) {
		for (int d = 0; d < 3; ++d) {
			startCenter[d] += cage.pStart[i][d] * 0.25;
			endCenter[d] += cage.pEnd[i][d] * 0.25;
		}
	}

	for (int d = 0; d < 3; ++d) {
		outCenter[d] = sel.isStart ? startCenter[d] : endCenter[d];
		outNormal[d] = sel.isStart ? (startCenter[d] - endCenter[d]) : (endCenter[d] - startCenter[d]);
	}

	if (vtkMath::Norm(outNormal) < 1e-8) {
		outNormal[0] = 1.0;
		outNormal[1] = 0.0;
		outNormal[2] = 0.0;
	}
	else {
		vtkMath::Normalize(outNormal);
	}

	int bestNode = -1;
	double bestDist2 = 1e100;

	for (int i = 0; i < static_cast<int>(m_skelNodes.size() / 3); ++i) {
		double p[3] = {
			m_skelNodes[3 * i],
			m_skelNodes[3 * i + 1],
			m_skelNodes[3 * i + 2]
		};

		double d2 = vtkMath::Distance2BetweenPoints(outCenter, p);

		if (d2 < bestDist2) {
			bestDist2 = d2;
			bestNode = i;
		}
	}

	return bestNode;
}



std::vector<int> CageInteractorStyle::FindBestSkeletonCycle(
	const std::vector<LoopAnchorInfo>& anchors,
	const std::vector<std::vector<int>>& adjList,
	std::vector<std::vector<int>>& connectorPaths,
	bool& isMultiBranchLoop)
{
	isMultiBranchLoop = false;
	connectorPaths.clear();

	if (anchors.empty() || adjList.empty()) {
		return std::vector<int>();
	}

	auto computeGraphDistances = [&](int startNode) {
		std::vector<int> dist(adjList.size(), -1);
		std::queue<int> q;

		if (startNode < 0 || startNode >= static_cast<int>(adjList.size())) {
			return dist;
		}

		dist[startNode] = 0;
		q.push(startNode);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			for (int next : adjList[curr]) {
				if (dist[next] != -1) continue;

				dist[next] = dist[curr] + 1;
				q.push(next);
			}
		}

		return dist;
	};

	auto buildPathFromParent = [](int endNode, const std::vector<int>& parent) {
		std::vector<int> path;
		int curr = endNode;
		while (curr != -1) {
			path.push_back(curr);
			curr = parent[curr];
		}
		std::reverse(path.begin(), path.end());
		return path;
	};

	auto findNearestPathToCycle = [&](int startNode, const std::set<int>& cycleKey) {
		std::vector<int> parent(adjList.size(), -1);
		std::vector<bool> visited(adjList.size(), false);
		std::queue<int> q;

		q.push(startNode);
		visited[startNode] = true;

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (cycleKey.find(curr) != cycleKey.end()) {
				return buildPathFromParent(curr, parent);
			}

			for (int next : adjList[curr]) {
				if (visited[next]) continue;

				visited[next] = true;
				parent[next] = curr;
				q.push(next);
			}
		}

		return std::vector<int>();
	};

	std::vector<std::vector<int>> anchorDistances;
	for (const auto& anchor : anchors) {
		anchorDistances.push_back(computeGraphDistances(anchor.node));
	}

	std::set<std::set<int>> visitedCycleKeys;
	std::vector<int> bestCycle;
	std::vector<std::vector<int>> bestConnectorPaths;
	double bestScore = 1e100;
	int bestBranchCount = 0;

	for (int u = 0; u < static_cast<int>(adjList.size()); ++u) {
		for (int v : adjList[u]) {
			if (u >= v) continue;

			std::vector<int> cycle = FindPathBFSIgnoringEdge(u, v, adjList, u, v);
			if (cycle.size() < 3) continue;

			std::set<int> cycleKey(cycle.begin(), cycle.end());
			if (visitedCycleKeys.find(cycleKey) != visitedCycleKeys.end()) {
				continue;
			}
			visitedCycleKeys.insert(cycleKey);

			double score = 0.0;
			std::vector<std::vector<int>> candidateConnectorPaths;
			bool allAnchorsInFront = true;

			for (size_t a = 0; a < anchorDistances.size(); ++a) {
				int minHop = 1000000;
				int closestCycleNode = -1;

				for (int node : cycle) {
					int d = anchorDistances[a][node];
					if (d >= 0 && d < minHop) {
						minHop = d;
						closestCycleNode = node;
					}
				}

				if (minHop == 1000000 || closestCycleNode < 0) {
					score += 1e9;
					allAnchorsInFront = false;
				}
				else {
					double closestPoint[3] = {
						m_skelNodes[3 * closestCycleNode],
						m_skelNodes[3 * closestCycleNode + 1],
						m_skelNodes[3 * closestCycleNode + 2]
					};

					double toCycle[3] = {
						closestPoint[0] - anchors[a].center[0],
						closestPoint[1] - anchors[a].center[1],
						closestPoint[2] - anchors[a].center[2]
					};

					double len = vtkMath::Norm(toCycle);
					double frontDot = -1.0;

					if (len > 1e-8) {
						toCycle[0] /= len;
						toCycle[1] /= len;
						toCycle[2] /= len;
						frontDot = vtkMath::Dot(toCycle, anchors[a].normal);
					}

					if (frontDot <= 0.05) {
						allAnchorsInFront = false;
						score += 1e8 + (0.05 - frontDot) * 1e6;
					}
					else {
						score += (1.0 - frontDot) * 5000.0;
					}

					std::vector<int> connectorPath = findNearestPathToCycle(anchors[a].node, cycleKey);
					if (connectorPath.empty()) {
						allAnchorsInFront = false;
						score += 1e9;
					}
					else {
						score += static_cast<double>(connectorPath.size()) * 10000.0;
						candidateConnectorPaths.push_back(connectorPath);
					}
				}
			}

			if (!allAnchorsInFront) {
				continue;
			}

			int branchCount = 0;
			for (int node : cycle) {
				if (adjList[node].size() >= 3) {
					++branchCount;
				}
			}

			score += static_cast<double>(cycle.size());
			score -= static_cast<double>(branchCount) * 0.25;

			if (score < bestScore) {
				bestScore = score;
				bestCycle = cycle;
				bestConnectorPaths = candidateConnectorPaths;
				bestBranchCount = branchCount;
			}
		}
	}

	if (!bestCycle.empty()) {
		connectorPaths = bestConnectorPaths;
		isMultiBranchLoop = (anchors.size() >= 2 || bestBranchCount >= 2);
	}

	return bestCycle;
}

bool CageInteractorStyle::DetectAndShowLoopTopology(vtkRenderer* ren)
{
	int selectedCount = static_cast<int>(m_selectedFacesForBranch.size());

	if (selectedCount != 1 && selectedCount != 2) {
		if (false) qDebug() << "[圆环检测] 请选择 1 个或 2 个控制笼端面。";
		return false;
	}

	if (m_skelNodes.empty() || m_skelEdges.empty()) {
		if (false) qDebug() << "[圆环检测] 骨架数据为空，请先使用带骨架的入口生成骨架。";
		return false;
	}

	std::vector<std::vector<int>> adjList = BuildSkeletonAdjacency();
	std::vector<LoopAnchorInfo> anchors;

	for (const auto& sel : m_selectedFacesForBranch) {
		LoopAnchorInfo anchor;
		anchor.node = FindNearestSkeletonNodeToSelectedFace(sel, anchor.center, anchor.normal);

		int node = anchor.node;
		if (node < 0) {
			if (false) qDebug() << "[圆环检测] 选中端面绑定到骨架失败。";
			return false;
		}

		anchors.push_back(anchor);
	}

	bool isMultiBranchLoop = false;
	std::vector<std::vector<int>> connectorPaths;
	std::vector<int> cycleNodes = FindBestSkeletonCycle(anchors, adjList, connectorPaths, isMultiBranchLoop);

	if (cycleNodes.empty()) {
		if (false) qDebug() << "[圆环检测] 未在选中端面前方找到闭合骨架圆环。";
		return false;
	}

	ShowSkeletonCycle(ren, cycleNodes, connectorPaths, isMultiBranchLoop);

	if (false) qDebug() << "[圆环检测] 选中端面数量 =" << selectedCount
		<< "锚点0 =" << anchors[0].node
		<< "锚点1 =" << (anchors.size() > 1 ? anchors[1].node : -1)
		<< "圆环骨架点数量 =" << static_cast<int>(cycleNodes.size())
		<< "连接路径数量 =" << static_cast<int>(connectorPaths.size())
		<< "类型 =" << (isMultiBranchLoop ? "多分支圆环" : "单分支圆环");

	return true;
}

void CageInteractorStyle::ShowSkeletonCycle(
	vtkRenderer* ren,
	const std::vector<int>& cycleNodes,
	const std::vector<std::vector<int>>& connectorPaths,
	bool isMultiBranchLoop)
{
	if (cycleNodes.size() < 3) return;

	vtkNew<vtkPoints> pts;
	for (size_t i = 0; i < m_skelNodes.size() / 3; ++i) {
		pts->InsertNextPoint(
			m_skelNodes[3 * i],
			m_skelNodes[3 * i + 1],
			m_skelNodes[3 * i + 2]
		);
	}

	vtkNew<vtkCellArray> lines;

	for (size_t i = 1; i < cycleNodes.size(); ++i) {
		vtkIdType ids[2] = {
			static_cast<vtkIdType>(cycleNodes[i - 1]),
			static_cast<vtkIdType>(cycleNodes[i])
		};
		lines->InsertNextCell(2, ids);
	}

	vtkIdType closeIds[2] = {
		static_cast<vtkIdType>(cycleNodes.back()),
		static_cast<vtkIdType>(cycleNodes.front())
	};
	lines->InsertNextCell(2, closeIds);

	for (const auto& path : connectorPaths) {
		for (size_t i = 1; i < path.size(); ++i) {
			vtkIdType ids[2] = {
				static_cast<vtkIdType>(path[i - 1]),
				static_cast<vtkIdType>(path[i])
			};
			lines->InsertNextCell(2, ids);
		}
	}

	vtkNew<vtkPolyData> pd;
	pd->SetPoints(pts);
	pd->SetLines(lines);

	if (!m_PathHighlightActor) {
		vtkNew<vtkPolyDataMapper> mapper;
		m_PathHighlightActor = vtkSmartPointer<vtkActor>::New();
		m_PathHighlightActor->SetMapper(mapper);
		m_PathHighlightActor->GetProperty()->SetRenderLinesAsTubes(true);
		m_PathHighlightActor->GetProperty()->SetLighting(false);
		m_PathHighlightActor->SetPickable(0);

		if (ren) {
			ren->AddActor(m_PathHighlightActor);
		}
	}

	vtkPolyDataMapper::SafeDownCast(m_PathHighlightActor->GetMapper())->SetInputData(pd);

	if (isMultiBranchLoop) {
		m_PathHighlightActor->GetProperty()->SetColor(1.0, 0.45, 0.0);
	}
	else {
		m_PathHighlightActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
	}

	m_PathHighlightActor->GetProperty()->SetLineWidth(10.0);
	m_PathHighlightActor->VisibilityOn();

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::ShowAutoTopologyTrace(
	vtkRenderer* ren,
	vtkSmartPointer<vtkActor>& actor,
	const std::vector<std::vector<int>>& paths,
	double r,
	double g,
	double b,
	double lineWidth)
{
	if (!ren || m_skelNodes.empty()) return;

	vtkNew<vtkPoints> pts;
	for (size_t i = 0; i < m_skelNodes.size() / 3; ++i) {
		pts->InsertNextPoint(
			m_skelNodes[3 * i],
			m_skelNodes[3 * i + 1],
			m_skelNodes[3 * i + 2]
		);
	}

	vtkNew<vtkCellArray> lines;

	for (const auto& path : paths) {
		if (path.size() < 2) continue;

		for (size_t i = 1; i < path.size(); ++i) {
			vtkIdType ids[2] = {
				static_cast<vtkIdType>(path[i - 1]),
				static_cast<vtkIdType>(path[i])
			};
			lines->InsertNextCell(2, ids);
		}
	}

	vtkNew<vtkPolyData> pd;
	pd->SetPoints(pts);
	pd->SetLines(lines);

	if (!actor) {
		vtkNew<vtkPolyDataMapper> mapper;
		actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetRenderLinesAsTubes(true);
		actor->GetProperty()->SetLighting(false);
		actor->SetPickable(0);
		ren->AddActor(actor);
	}

	vtkPolyDataMapper::SafeDownCast(actor->GetMapper())->SetInputData(pd);
	actor->GetProperty()->SetColor(r, g, b);
	actor->GetProperty()->SetLineWidth(lineWidth);

	if (lines->GetNumberOfCells() > 0) {
		actor->VisibilityOn();
	}
	else {
		actor->VisibilityOff();
	}
}

void CageInteractorStyle::AutoDetectAndShowTopology(vtkRenderer* ren, bool generateTopology)
{
	if (!ren || m_allCages.empty() || m_skelNodes.empty() || m_skelEdges.empty()) {
		if (false) qDebug() << "[自动识别] 缺少控制笼或骨架数据。";
		return;
	}

	std::vector<std::vector<int>> adjList = BuildSkeletonAdjacency();

	auto makeEdgeKey = [](int a, int b) {
		if (a > b) std::swap(a, b);
		return std::make_pair(a, b);
	};

	auto isValidNode = [&](int node) {
		return node >= 0 && node < static_cast<int>(adjList.size());
	};

	// 自动识别时，已经生成直管的骨架边相当于被端面占住，后续寻路不能穿过去。
	std::vector<std::pair<int, int>> cageEndNodes(m_allCages.size(), std::make_pair(-1, -1));
	std::set<int> allCageCapNodes;

	for (int cageIndex = 0; cageIndex < static_cast<int>(m_allCages.size()); ++cageIndex) {
		for (int capSide = 0; capSide < 2; ++capSide) {
			SelectedFace sel;
			sel.cageIndex = cageIndex;
			sel.isStart = (capSide == 0);

			double center[3];
			double normal[3];
			int node = FindNearestSkeletonNodeToSelectedFace(sel, center, normal);

			if (!isValidNode(node)) continue;

			if (sel.isStart) {
				cageEndNodes[cageIndex].first = node;
			}
			else {
				cageEndNodes[cageIndex].second = node;
			}

			allCageCapNodes.insert(node);
		}
	}

	std::set<std::pair<int, int>> occupiedEdges;
	std::set<int> occupiedInnerNodes;

	for (const auto& ends : cageEndNodes) {
		if (!isValidNode(ends.first) || !isValidNode(ends.second)) continue;

		std::vector<int> path = FindPathBFS(ends.first, ends.second, adjList);
		if (path.size() < 2) continue;

		for (size_t i = 1; i < path.size(); ++i) {
			occupiedEdges.insert(makeEdgeKey(path[i - 1], path[i]));
		}

		for (size_t i = 1; i + 1 < path.size(); ++i) {
			if (allCageCapNodes.find(path[i]) == allCageCapNodes.end()) {
				occupiedInnerNodes.insert(path[i]);
			}
		}
	}

	struct FreeCapInfo {
		SelectedFace sel;
		LoopAnchorInfo anchor;
		double ring[4][3];
	};

	auto capRingsCoincident = [](const double a[4][3], const double b[4][3]) -> bool {
		double maxMinD2 = 0.0;

		for (int i = 0; i < 4; ++i) {
			double bestD2 = 1e100;

			for (int j = 0; j < 4; ++j) {
				bestD2 = std::min(bestD2, vtkMath::Distance2BetweenPoints(a[i], b[j]));
			}

			maxMinD2 = std::max(maxMinD2, bestD2);
		}

		return maxMinD2 < 1e-5;
	};

	std::vector<FreeCapInfo> caps;

	for (int cageIndex = 0; cageIndex < static_cast<int>(m_allCages.size()); ++cageIndex) {
		for (int capSide = 0; capSide < 2; ++capSide) {
			SelectedFace sel;
			sel.cageIndex = cageIndex;
			sel.isStart = (capSide == 0);

			double ring[4][3];
			if (!GetCanonicalCageCapRing(m_allCages, sel.cageIndex, sel.isStart, ring)) {
				continue;
			}

			bool sharedWithOtherCap = false;

			for (int otherIndex = 0; otherIndex < static_cast<int>(m_allCages.size()) && !sharedWithOtherCap; ++otherIndex) {
				for (int otherSide = 0; otherSide < 2; ++otherSide) {
					if (otherIndex == cageIndex && otherSide == capSide) {
						continue;
					}

					double otherRing[4][3];
					if (!GetCanonicalCageCapRing(m_allCages, otherIndex, otherSide == 0, otherRing)) {
						continue;
					}

					if (capRingsCoincident(ring, otherRing)) {
						sharedWithOtherCap = true;
						break;
					}
				}
			}

			if (sharedWithOtherCap) {
				continue;
			}

			LoopAnchorInfo anchor;
			anchor.node = FindNearestSkeletonNodeToSelectedFace(sel, anchor.center, anchor.normal);
			if (anchor.node < 0) {
				continue;
			}

			FreeCapInfo info;
			info.sel = sel;
			info.anchor = anchor;
			CopyRing4(info.ring, ring);
			caps.push_back(info);
		}
	}

	std::set<int> freeCapNodes;
	for (const auto& cap : caps) {
		freeCapNodes.insert(cap.anchor.node);
	}

	auto isSearchStepBlocked = [&](int curr, int next, int startNode, int endNode, int blockedNode) {
		if (!isValidNode(next)) return true;
		if (next == blockedNode) return true;
		if (occupiedEdges.find(makeEdgeKey(curr, next)) != occupiedEdges.end()) return true;
		if (next != startNode && next != endNode &&
			occupiedInnerNodes.find(next) != occupiedInnerNodes.end()) {
			return true;
		}
		if (next != startNode && next != endNode &&
			freeCapNodes.find(next) != freeCapNodes.end()) {
			return true;
		}
		return false;
	};

	auto findConstrainedPath = [&](int startNode, int endNode, int blockedNode) {
		std::vector<int> emptyPath;

		if (startNode < 0 || endNode < 0 ||
			startNode >= static_cast<int>(adjList.size()) ||
			endNode >= static_cast<int>(adjList.size()) ||
			startNode == blockedNode || endNode == blockedNode) {
			return emptyPath;
		}

		std::vector<int> parent(adjList.size(), -1);
		std::vector<bool> visited(adjList.size(), false);
		std::queue<int> q;

		visited[startNode] = true;
		q.push(startNode);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (curr == endNode) break;

			for (int next : adjList[curr]) {
				if (visited[next]) continue;
				if (isSearchStepBlocked(curr, next, startNode, endNode, blockedNode)) continue;

				visited[next] = true;
				parent[next] = curr;
				q.push(next);
			}
		}

		if (!visited[endNode]) return emptyPath;

		std::vector<int> path;
		for (int curr = endNode; curr != -1; curr = parent[curr]) {
			path.push_back(curr);
		}

		std::reverse(path.begin(), path.end());
		return path;
	};

	auto makeSearchAdjacency = [&](const std::set<int>& allowedCapNodes) {
		std::vector<std::vector<int>> filtered(adjList.size());

		for (int u = 0; u < static_cast<int>(adjList.size()); ++u) {
			if (occupiedInnerNodes.find(u) != occupiedInnerNodes.end()) continue;
			if (freeCapNodes.find(u) != freeCapNodes.end() &&
				allowedCapNodes.find(u) == allowedCapNodes.end()) {
				continue;
			}

			for (int v : adjList[u]) {
				if (occupiedEdges.find(makeEdgeKey(u, v)) != occupiedEdges.end()) continue;
				if (occupiedInnerNodes.find(v) != occupiedInnerNodes.end()) continue;
				if (freeCapNodes.find(v) != freeCapNodes.end() &&
					allowedCapNodes.find(v) == allowedCapNodes.end()) {
					continue;
				}

				filtered[u].push_back(v);
			}
		}

		return filtered;
	};

	auto findNearestForwardJunction = [&](int capIndex, std::vector<int>& outPath) -> int {
		outPath.clear();

		const FreeCapInfo& cap = caps[capIndex];
		std::vector<int> parent(adjList.size(), -1);
		std::vector<bool> visited(adjList.size(), false);
		std::queue<int> q;

		visited[cap.anchor.node] = true;
		q.push(cap.anchor.node);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (curr != cap.anchor.node && adjList[curr].size() >= 3) {
				double p[3] = {
					m_skelNodes[3 * curr],
					m_skelNodes[3 * curr + 1],
					m_skelNodes[3 * curr + 2]
				};
				double v[3] = {
					p[0] - cap.anchor.center[0],
					p[1] - cap.anchor.center[1],
					p[2] - cap.anchor.center[2]
				};

				double len = vtkMath::Norm(v);
				double frontDot = -1.0;

				if (len > 1e-8) {
					v[0] /= len;
					v[1] /= len;
					v[2] /= len;
					frontDot = vtkMath::Dot(v, cap.anchor.normal);
				}

				if (frontDot > 0.0) {
					for (int node = curr; node != -1; node = parent[node]) {
						outPath.push_back(node);
					}

					std::reverse(outPath.begin(), outPath.end());
					return curr;
				}
			}

			for (int next : adjList[curr]) {
				if (visited[next]) continue;
				if (isSearchStepBlocked(curr, next, cap.anchor.node, -1, -1)) continue;

				visited[next] = true;
				parent[next] = curr;
				q.push(next);
			}
		}

		return -1;
	};

	auto traceTwoInputLoop = [&](int capA, int capB,
		std::vector<int>& outCycle,
		std::vector<std::vector<int>>& outConnectors,
		int& outJunction) -> bool {
		outCycle.clear();
		outConnectors.clear();
		outJunction = -1;

		std::vector<int> inputPath = findConstrainedPath(caps[capA].anchor.node, caps[capB].anchor.node, -1);
		if (inputPath.size() < 3) return false;

		int junctionIndex = -1;
		int bestDegree = -1;

		for (int i = 1; i + 1 < static_cast<int>(inputPath.size()); ++i) {
			int node = inputPath[i];
			int degree = static_cast<int>(adjList[node].size());

			if (degree >= 4 && degree > bestDegree) {
				bestDegree = degree;
				junctionIndex = i;
			}
		}

		if (junctionIndex < 0) {
			for (int i = 1; i + 1 < static_cast<int>(inputPath.size()); ++i) {
				int node = inputPath[i];
				int degree = static_cast<int>(adjList[node].size());

				if (degree >= 3 && degree > bestDegree) {
					bestDegree = degree;
					junctionIndex = i;
				}
			}
		}

		if (junctionIndex < 0) return false;

		int junction = inputPath[junctionIndex];
		std::set<int> inputNeighbors;
		inputNeighbors.insert(inputPath[junctionIndex - 1]);
		inputNeighbors.insert(inputPath[junctionIndex + 1]);

		std::vector<int> loopNeighbors;
		for (int n : adjList[junction]) {
			if (inputNeighbors.find(n) == inputNeighbors.end()) {
				loopNeighbors.push_back(n);
			}
		}

		std::vector<int> bestLoopPath;

		for (size_t i = 0; i < loopNeighbors.size(); ++i) {
			for (size_t j = i + 1; j < loopNeighbors.size(); ++j) {
				std::vector<int> candidatePath =
					findConstrainedPath(loopNeighbors[i], loopNeighbors[j], junction);

				if (candidatePath.empty()) continue;

				if (bestLoopPath.empty() || candidatePath.size() < bestLoopPath.size()) {
					bestLoopPath = candidatePath;
				}
			}
		}

		if (bestLoopPath.empty()) return false;

		outJunction = junction;
		outCycle.push_back(junction);
		outCycle.insert(outCycle.end(), bestLoopPath.begin(), bestLoopPath.end());
		outConnectors.push_back(findConstrainedPath(caps[capA].anchor.node, junction, -1));
		outConnectors.push_back(findConstrainedPath(caps[capB].anchor.node, junction, -1));

		if (outConnectors[0].empty() || outConnectors[1].empty()) {
			outCycle.clear();
			outConnectors.clear();
			outJunction = -1;
			return false;
		}

		return true;
	};

	auto makeCyclePath = [](const std::vector<int>& cycleNodes) {
		std::vector<int> path = cycleNodes;
		if (!path.empty()) {
			path.push_back(path.front());
		}
		return path;
	};

	auto computeAutoSelectedFaceFrame = [&](const SelectedFace& sel, double center[3], double normal[3]) -> bool
	{
		if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
			return false;
		}

		const CageSegment& cage = m_allCages[sel.cageIndex];
		const double(*ring)[3] = sel.isStart ? cage.pStart : cage.pEnd;
		const double(*oppositeRing)[3] = sel.isStart ? cage.pEnd : cage.pStart;

		center[0] = center[1] = center[2] = 0.0;
		double oppositeCenter[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				center[d] += ring[i][d] * 0.25;
				oppositeCenter[d] += oppositeRing[i][d] * 0.25;
			}
		}

		for (int d = 0; d < 3; ++d) {
			normal[d] = center[d] - oppositeCenter[d];
		}

		double nLen = vtkMath::Norm(normal);
		if (nLen < 1e-8) {
			double e01[3] = {
				ring[1][0] - ring[0][0],
				ring[1][1] - ring[0][1],
				ring[1][2] - ring[0][2]
			};
			double e03[3] = {
				ring[3][0] - ring[0][0],
				ring[3][1] - ring[0][1],
				ring[3][2] - ring[0][2]
			};
			vtkMath::Cross(e01, e03, normal);
			nLen = vtkMath::Norm(normal);
		}

		if (nLen < 1e-8) {
			normal[0] = 1.0;
			normal[1] = 0.0;
			normal[2] = 0.0;
		}
		else {
			vtkMath::Normalize(normal);
		}

		return true;
	};

	auto shouldUse2DCenterHubForAutoBranch = [&](const std::vector<SelectedFace>& faces) -> bool
	{
		const int count = static_cast<int>(faces.size());

		if (count == 3) {
			if (false) qDebug() << "[AutoTopology] selected faces = 3, use 2D center hub topology.";
			return true;
		}

		if (count < 4) {
			return false;
		}

		std::vector<std::array<double, 3>> centers;
		std::vector<std::array<double, 3>> normals;
		centers.reserve(count);
		normals.reserve(count);

		for (const SelectedFace& sel : faces) {
			double c[3];
			double n[3];
			if (!computeAutoSelectedFaceFrame(sel, c, n)) {
				if (false) qDebug() << "[AutoTopology] invalid complex face, fallback to 3D branch. cageIndex =" << sel.cageIndex;
				return false;
			}

			centers.push_back({ c[0], c[1], c[2] });
			normals.push_back({ n[0], n[1], n[2] });
		}

		double centroid[3] = { 0.0, 0.0, 0.0 };
		for (const auto& c : centers) {
			centroid[0] += c[0] / static_cast<double>(count);
			centroid[1] += c[1] / static_cast<double>(count);
			centroid[2] += c[2] / static_cast<double>(count);
		}

		double scale = 0.0;
		for (const auto& c : centers) {
			double p[3] = { c[0], c[1], c[2] };
			scale = std::max(scale, std::sqrt(vtkMath::Distance2BetweenPoints(p, centroid)));
		}
		scale = std::max(scale, 1.0);

		double planeNormal[3] = { 0.0, 0.0, 0.0 };
		double bestCrossLen = -1.0;

		for (int i = 0; i < count; ++i) {
			double a[3] = {
				centers[i][0] - centroid[0],
				centers[i][1] - centroid[1],
				centers[i][2] - centroid[2]
			};

			for (int j = i + 1; j < count; ++j) {
				double b[3] = {
					centers[j][0] - centroid[0],
					centers[j][1] - centroid[1],
					centers[j][2] - centroid[2]
				};

				double cr[3];
				vtkMath::Cross(a, b, cr);
				double crLen = vtkMath::Norm(cr);

				if (crLen > bestCrossLen) {
					bestCrossLen = crLen;
					planeNormal[0] = cr[0];
					planeNormal[1] = cr[1];
					planeNormal[2] = cr[2];
				}
			}
		}

		if (vtkMath::Norm(planeNormal) < 1e-8) {
			for (int i = 0; i < count; ++i) {
				double ni[3] = { normals[i][0], normals[i][1], normals[i][2] };
				for (int j = i + 1; j < count; ++j) {
					double nj[3] = { normals[j][0], normals[j][1], normals[j][2] };
					double cr[3];
					vtkMath::Cross(ni, nj, cr);
					if (vtkMath::Norm(cr) > vtkMath::Norm(planeNormal)) {
						planeNormal[0] = cr[0];
						planeNormal[1] = cr[1];
						planeNormal[2] = cr[2];
					}
				}
			}
		}

		if (vtkMath::Norm(planeNormal) < 1e-8) {
			if (false) qDebug() << "[AutoTopology] cannot estimate common plane, fallback to 3D branch.";
			return false;
		}

		vtkMath::Normalize(planeNormal);

		double maxCenterPlaneDist = 0.0;
		for (const auto& c : centers) {
			double v[3] = {
				c[0] - centroid[0],
				c[1] - centroid[1],
				c[2] - centroid[2]
			};
			maxCenterPlaneDist = std::max(maxCenterPlaneDist, std::abs(vtkMath::Dot(v, planeNormal)));
		}

		double maxNormalOutOfPlane = 0.0;
		double avgNormalOutOfPlane = 0.0;
		for (const auto& n : normals) {
			double nv[3] = { n[0], n[1], n[2] };
			double v = std::abs(vtkMath::Dot(nv, planeNormal));
			maxNormalOutOfPlane = std::max(maxNormalOutOfPlane, v);
			avgNormalOutOfPlane += v / static_cast<double>(count);
		}

		const double centerPlaneTol = std::max(1e-4, scale * 0.22);
		const double normalPlaneTol = 0.38;

		const bool centersNearPlane = maxCenterPlaneDist <= centerPlaneTol;
		const bool normalsNearPlane = maxNormalOutOfPlane <= normalPlaneTol;
		const bool use2D = centersNearPlane && normalsNearPlane;

		if (false) qDebug() << "[AutoTopology] complex branch faces =" << count
			<< "planeNormal =" << planeNormal[0] << planeNormal[1] << planeNormal[2]
			<< "centerMaxDist =" << maxCenterPlaneDist
			<< "centerTol =" << centerPlaneTol
			<< "normalOutMax =" << maxNormalOutOfPlane
			<< "normalOutAvg =" << avgNormalOutOfPlane
			<< "normalTol =" << normalPlaneTol
			<< "centersNearPlane =" << centersNearPlane
			<< "normalsNearPlane =" << normalsNearPlane
			<< "result =" << (use2D ? "2D center hub" : "3D side-branch");

		return use2D;
	};

	std::vector<bool> used(caps.size(), false);
	std::vector<std::vector<int>> complex2DPaths;
	std::vector<std::vector<int>> complex3DPaths;
	std::vector<std::vector<int>> singleLoopPaths;
	std::vector<std::vector<int>> multiLoopPaths;

	struct AutoLoopGroup {
		std::vector<SelectedFace> faces;
		std::vector<int> cycle;
		std::vector<std::vector<int>> connectors;
	};

	struct AutoComplexGroup {
		std::vector<SelectedFace> faces;
		std::vector<std::vector<int>> paths;
		bool use2D;
	};

	std::vector<AutoComplexGroup> complexFaceGroups;
	std::vector<AutoLoopGroup> singleLoopGroups;
	std::vector<AutoLoopGroup> multiLoopGroups;

	std::map<int, std::vector<std::pair<int, std::vector<int>>>> complexGroups;

	for (int i = 0; i < static_cast<int>(caps.size()); ++i) {
		std::vector<int> path;
		int junction = findNearestForwardJunction(i, path);

		if (junction >= 0 && !path.empty()) {
			complexGroups[junction].push_back(std::make_pair(i, path));
		}
	}

	for (const auto& item : complexGroups) {
		if (item.second.size() < 3) continue;

		AutoComplexGroup group;
		group.use2D = true;

		for (const auto& capPath : item.second) {
			used[capPath.first] = true;
			group.paths.push_back(capPath.second);
			group.faces.push_back(caps[capPath.first].sel);
		}

		group.use2D = shouldUse2DCenterHubForAutoBranch(group.faces);
		if (group.use2D) {
			for (const auto& path : group.paths) {
				complex2DPaths.push_back(path);
			}
		}
		else {
			for (const auto& path : group.paths) {
				complex3DPaths.push_back(path);
			}
		}

		complexFaceGroups.push_back(group);
	}

	struct PairLoopCandidate {
		int a;
		int b;
		int junction;
		std::vector<int> cycle;
		std::vector<std::vector<int>> connectors;
		double score;
	};

	std::vector<PairLoopCandidate> pairCandidates;

	for (int i = 0; i < static_cast<int>(caps.size()); ++i) {
		if (used[i]) continue;

		for (int j = i + 1; j < static_cast<int>(caps.size()); ++j) {
			if (used[j]) continue;

			PairLoopCandidate candidate;
			candidate.a = i;
			candidate.b = j;

			if (!traceTwoInputLoop(i, j, candidate.cycle, candidate.connectors, candidate.junction)) {
				continue;
			}

			candidate.score = static_cast<double>(candidate.cycle.size());
			for (const auto& p : candidate.connectors) {
				candidate.score += static_cast<double>(p.size()) * 0.5;
			}

			pairCandidates.push_back(candidate);
		}
	}

	std::sort(
		pairCandidates.begin(),
		pairCandidates.end(),
		[](const PairLoopCandidate& a, const PairLoopCandidate& b) {
		return a.score < b.score;
	}
	);

	for (const auto& candidate : pairCandidates) {
		if (used[candidate.a] || used[candidate.b]) {
			continue;
		}

		used[candidate.a] = true;
		used[candidate.b] = true;

		multiLoopPaths.push_back(makeCyclePath(candidate.cycle));
		for (const auto& p : candidate.connectors) {
			multiLoopPaths.push_back(p);
		}

		AutoLoopGroup group;
		group.faces.push_back(caps[candidate.a].sel);
		group.faces.push_back(caps[candidate.b].sel);
		group.cycle = candidate.cycle;
		group.connectors = candidate.connectors;
		multiLoopGroups.push_back(group);
	}

	for (int i = 0; i < static_cast<int>(caps.size()); ++i) {
		if (used[i]) continue;

		bool isMultiBranchLoop = false;
		std::vector<std::vector<int>> connectorPaths;
		std::vector<LoopAnchorInfo> anchors;
		anchors.push_back(caps[i].anchor);

		std::set<int> allowedCapNodes;
		allowedCapNodes.insert(caps[i].anchor.node);
		std::vector<std::vector<int>> filteredAdjList = makeSearchAdjacency(allowedCapNodes);

		std::vector<int> cycleNodes = FindBestSkeletonCycle(
			anchors,
			filteredAdjList,
			connectorPaths,
			isMultiBranchLoop
		);

		if (!cycleNodes.empty() && !connectorPaths.empty() && !isMultiBranchLoop) {
			used[i] = true;
			singleLoopPaths.push_back(makeCyclePath(cycleNodes));

			for (const auto& p : connectorPaths) {
				singleLoopPaths.push_back(p);
			}

			AutoLoopGroup group;
			group.faces.push_back(caps[i].sel);
			group.cycle = cycleNodes;
			group.connectors = connectorPaths;
			singleLoopGroups.push_back(group);
		}
	}

	if (!generateTopology) {
		// 自动识别路径颜色：二维复杂分支 = 黄色；三维复杂分支 = 紫色；
		// 单输入 loop = 青色；双输入 loop = 橙色。
		ShowAutoTopologyTrace(ren, m_AutoComplexTraceActor, complex2DPaths, 1.0, 1.0, 0.0, 9.0);
		ShowAutoTopologyTrace(ren, m_AutoComplex3DTraceActor, complex3DPaths, 0.75, 0.25, 1.0, 10.0);
		ShowAutoTopologyTrace(ren, m_AutoSingleLoopTraceActor, singleLoopPaths, 0.0, 0.9, 1.0, 8.0);
		ShowAutoTopologyTrace(ren, m_AutoMultiLoopTraceActor, multiLoopPaths, 1.0, 0.45, 0.0, 8.0);
	}

	if (generateTopology) {
		int generatedComplex = 0;
		int generatedMultiLoop = 0;
		int generatedSingleLoop = 0;

		bool oldBranchMode = m_isBranchMode;

		if (m_BranchHighlightActor) m_BranchHighlightActor->VisibilityOff();
		if (m_PathHighlightActor) m_PathHighlightActor->VisibilityOff();
		if (m_AutoComplexTraceActor) m_AutoComplexTraceActor->VisibilityOff();
		if (m_AutoComplex3DTraceActor) m_AutoComplex3DTraceActor->VisibilityOff();
		if (m_AutoSingleLoopTraceActor) m_AutoSingleLoopTraceActor->VisibilityOff();
		if (m_AutoMultiLoopTraceActor) m_AutoMultiLoopTraceActor->VisibilityOff();

		for (const auto& group : complexFaceGroups) {
			if (group.faces.size() < 3) continue;

			m_selectedFacesForBranch = group.faces;
			if (group.use2D) {
				if (false) qDebug() << "[AutoTopology] generate complex branch -> GenerateNBranch(), 2D center hub.";
				GenerateNBranch();
			}
			else {
				if (false) qDebug() << "[AutoTopology] generate complex branch -> GenerateNBranch3DSkeleton(), 3D side-branch.";
				GenerateNBranch3DSkeleton();
			}
			++generatedComplex;
		}

		auto generateLoopGroup = [&](const AutoLoopGroup& group) -> bool {
			if (group.faces.empty() || group.cycle.empty() || group.connectors.empty()) {
				return false;
			}

			m_selectedFacesForBranch = group.faces;
			m_useForcedLoopTopology = true;
			m_forcedLoopCycleNodes = group.cycle;
			m_forcedLoopConnectorPaths = group.connectors;

			bool ok = GenerateSingleBranchLoopTopology(ren);

			m_useForcedLoopTopology = false;
			m_forcedLoopCycleNodes.clear();
			m_forcedLoopConnectorPaths.clear();

			return ok;
		};

		for (const auto& group : multiLoopGroups) {
			if (generateLoopGroup(group)) {
				++generatedMultiLoop;
			}
		}

		for (const auto& group : singleLoopGroups) {
			if (generateLoopGroup(group)) {
				++generatedSingleLoop;
			}
		}

		m_selectedFacesForBranch.clear();
		m_isBranchMode = oldBranchMode;

		if (m_BranchHighlightActor) {
			m_BranchHighlightActor->VisibilityOff();
		}

		if (m_AutoComplexTraceActor) m_AutoComplexTraceActor->VisibilityOff();
		if (m_AutoComplex3DTraceActor) m_AutoComplex3DTraceActor->VisibilityOff();
		if (m_AutoSingleLoopTraceActor) m_AutoSingleLoopTraceActor->VisibilityOff();
		if (m_AutoMultiLoopTraceActor) m_AutoMultiLoopTraceActor->VisibilityOff();

	}

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::ShowSkeletonPath(vtkRenderer* ren) {
	if (m_skelNodes.empty()) return;

	// 1. 
	std::vector<int> targetNodes;
	for (const auto& sel : m_selectedFacesForBranch) {
		double center[3] = { 0,0,0 };
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 3; j++) {
				center[j] += (sel.isStart ? m_allCages[sel.cageIndex].pStart[i][j] : m_allCages[sel.cageIndex].pEnd[i][j]) / 4.0;
			}
		}
		int bestNode = 0;
		double minDist = 1e18;
		for (size_t i = 0; i < m_skelNodes.size() / 3; i++) {
			double d = vtkMath::Distance2BetweenPoints(center, new double[3]{ m_skelNodes[3 * i], m_skelNodes[3 * i + 1], m_skelNodes[3 * i + 2] });
			if (d < minDist) { minDist = d; bestNode = i; }
		}
		targetNodes.push_back(bestNode);
	}

	// 2. 
	int numNodes = m_skelNodes.size() / 3;
	std::vector<std::vector<int>> adjList(numNodes);
	for (size_t i = 0; i < m_skelEdges.size() / 2; i++) {
		int u = m_skelEdges[2 * i];
		int v = m_skelEdges[2 * i + 1];
		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}

	// 3.
	std::set<std::pair<int, int>> pathEdges;
	auto addPathEdges = [&](const std::vector<int>& p) {
		for (size_t i = 1; i < p.size(); i++) {
			int u = p[i - 1], v = p[i];
			pathEdges.insert({ std::min(u,v), std::max(u,v) });
		}
	};

	for (size_t i = 1; i < targetNodes.size(); i++) {
		std::vector<int> path = FindPathBFS(targetNodes[i], targetNodes[0], adjList);
		addPathEdges(path);
	}

	// 4. 
	vtkNew<vtkPoints> pts;
	for (size_t i = 0; i < m_skelNodes.size() / 3; i++) {
		pts->InsertNextPoint(m_skelNodes[3 * i], m_skelNodes[3 * i + 1], m_skelNodes[3 * i + 2]);
	}
	vtkNew<vtkCellArray> lines;
	for (auto edge : pathEdges) {
		vtkNew<vtkLine> line;
		line->GetPointIds()->SetId(0, edge.first);
		line->GetPointIds()->SetId(1, edge.second);
		lines->InsertNextCell(line);
	}

	vtkNew<vtkPolyData> pd; pd->SetPoints(pts); pd->SetLines(lines);

	if (!m_PathHighlightActor) {
		vtkNew<vtkPolyDataMapper> mapper;
		m_PathHighlightActor = vtkSmartPointer<vtkActor>::New();
		m_PathHighlightActor->SetMapper(mapper);
		m_PathHighlightActor->GetProperty()->SetColor(1.0, 1.0, 0.0); // 
		m_PathHighlightActor->GetProperty()->SetLineWidth(8.0);       // 
		m_PathHighlightActor->GetProperty()->SetRenderLinesAsTubes(true); // 
		m_PathHighlightActor->GetProperty()->SetLighting(false);
		if (ren) ren->AddActor(m_PathHighlightActor);
	}
	vtkPolyDataMapper::SafeDownCast(m_PathHighlightActor->GetMapper())->SetInputData(pd);
	m_PathHighlightActor->VisibilityOn();
}




void CageInteractorStyle::EnableSkeletonNodeEditMode(bool enable)
{
	if (m_isSkeletonEditMode == enable) {
		return;
	}

	m_isSkeletonEditMode = enable;
	if (enable)
	{
		m_isSkeletonDeleteMode = false;
		m_isSkeletonExtendMode = false;
		m_isStraightCageEditMode = false;
		m_isStraightCageRotateMode = false;
		m_isStraightSegmentEditMode = false;
		m_isStraightSegmentSelectMode = false;
		CancelStraightCageEdit(nullptr);
		CancelStraightCageRotation(nullptr);
		HideManualPreviewActors();
		m_selectedExtendSkelNode = -1;
		m_isBranchMode = false;
		m_selectedFacesForBranch.clear();
		if (m_BranchHighlightActor) m_BranchHighlightActor->VisibilityOff();
		qDebug() << "[骨架节点编辑] 进入：依次点击两个骨架节点进行合并。";
	}
	else
	{
		m_selectedSkelNode1 = -1;
		if (m_skelNodeHighlightActor) m_skelNodeHighlightActor->VisibilityOff();
		qDebug() << "[骨架节点编辑] 退出";
	}

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow())
	{
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	NotifyEditModeChanged();
}

void CageInteractorStyle::EnableSkeletonDeleteMode(bool enable)
{
	if (m_isSkeletonDeleteMode == enable) {
		return;
	}

	m_isSkeletonDeleteMode = enable;
	if (enable)
	{
		m_isSkeletonEditMode = false;
		m_isSkeletonExtendMode = false;
		m_isStraightCageEditMode = false;
		m_isStraightCageRotateMode = false;
		m_isStraightSegmentEditMode = false;
		m_isStraightSegmentSelectMode = false;
		CancelStraightCageEdit(nullptr);
		CancelStraightCageRotation(nullptr);
		HideManualPreviewActors();
		m_selectedSkelNode1 = -1;
		m_selectedExtendSkelNode = -1;
		if (m_skelNodeHighlightActor) m_skelNodeHighlightActor->VisibilityOff();
		qDebug() << "[骨架编辑] 进入删除模式：点击一条骨架边删除。";
	}
	else
	{
		qDebug() << "[骨架编辑] 退出删除模式";
	}

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow())
	{
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	NotifyEditModeChanged();
}

void CageInteractorStyle::EnableSkeletonExtendMode(bool enable)
{
	if (m_isSkeletonExtendMode == enable) {
		return;
	}

	m_isSkeletonExtendMode = enable;
	if (enable)
	{
		m_isSkeletonEditMode = false;
		m_isSkeletonDeleteMode = false;
		m_isStraightCageEditMode = false;
		m_isStraightCageRotateMode = false;
		m_isStraightSegmentEditMode = false;
		m_isStraightSegmentSelectMode = false;
		CancelStraightCageEdit(nullptr);
		CancelStraightCageRotation(nullptr);
		HideManualPreviewActors();
		m_selectedExtendSkelNode = -1;
		m_selectedSkelNode1 = -1;
		HideSkeletonExtendPreview();
		if (m_skelNodeHighlightActor) m_skelNodeHighlightActor->VisibilityOff();
		qDebug() << "[骨架延伸] 进入延伸/插入模式：左键选骨架节点，移动预览，再左键确认；右键取消。";
	}
	else
	{
		ClearSkeletonExtendSelection(nullptr);
		qDebug() << "[骨架延伸] 退出延伸模式";
	}

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow())
	{
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	NotifyEditModeChanged();
}

bool CageInteractorStyle::DeleteSkeletonEdge(int edgeIndex)
{
	const int edgeCount = static_cast<int>(m_skelEdges.size() / 2);
	if (edgeIndex < 0 || edgeIndex >= edgeCount)
	{
		return false;
	}

	std::vector<unsigned int> newEdges;
	newEdges.reserve(m_skelEdges.size() - 2);
	for (int i = 0; i < edgeCount; ++i)
	{
		if (i == edgeIndex)
		{
			continue;
		}
		newEdges.push_back(m_skelEdges[2 * i]);
		newEdges.push_back(m_skelEdges[2 * i + 1]);
	}

	m_skelEdges.swap(newEdges);
	RefreshSkeletonPolyData();
	return true;
}


void CageInteractorStyle::HideSkeletonExtendPreview()
{
	if (m_skeletonExtendPreviewActor) {
		m_skeletonExtendPreviewActor->VisibilityOff();
	}
	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

void CageInteractorStyle::ClearSkeletonExtendSelection(vtkRenderer* ren)
{
	m_selectedExtendSkelNode = -1;
	m_selectedSkelNode1 = -1;
	HideSkeletonExtendPreview();
	if (ren) {
		UpdateSkeletonHighlight(ren);
	}
}

bool CageInteractorStyle::GetSkeletonExtensionTargetPoint(
	int x,
	int y,
	vtkRenderer* ren,
	int tail,
	double outPoint[3]) const
{
	const int numNodes = static_cast<int>(m_skelNodes.size() / 3);
	if (!ren || tail < 0 || tail >= numNodes || !outPoint) {
		return false;
	}

	double pTail[3] = {
		m_skelNodes[3 * tail],
		m_skelNodes[3 * tail + 1],
		m_skelNodes[3 * tail + 2]
	};

	// 关键修改：延伸骨架时不能再直接拾取模型表面。
	// 否则鼠标如果落到灰色拓扑体或绿色 non-design 的上表面，
	// 延伸终点会被吸到该表面上，导致骨架线向上漂移。
	//
	// 当前采用“过尾端点、平行屏幕的工作平面”：
	// 1. 先把尾端点 P0 投影到屏幕，取得它的 display-z；
	// 2. 再把当前鼠标屏幕坐标 (x,y) 以同一个 display-z 反投影回世界坐标；
	// 3. 这个世界点作为延伸预览/确认点。
	//
	// 这样鼠标只控制延伸方向和长度，不会被模型表面拾取点带偏。
	ren->SetWorldPoint(pTail[0], pTail[1], pTail[2], 1.0);
	ren->WorldToDisplay();
	double* tailDisplay = ren->GetDisplayPoint();
	if (!tailDisplay) {
		return false;
	}

	const double fixedDisplayZ = tailDisplay[2];

	ren->SetDisplayPoint(static_cast<double>(x), static_cast<double>(y), fixedDisplayZ);
	ren->DisplayToWorld();
	double* world = ren->GetWorldPoint();
	if (!world || std::abs(world[3]) < 1e-12) {
		return false;
	}

	outPoint[0] = world[0] / world[3];
	outPoint[1] = world[1] / world[3];
	outPoint[2] = world[2] / world[3];

	return vtkMath::Distance2BetweenPoints(outPoint, pTail) > 1e-12;
}

bool CageInteractorStyle::UpdateSkeletonExtendPreview(int x, int y, vtkRenderer* ren)
{
	const int tail = m_selectedExtendSkelNode;
	const int neighbor = GetSingleSkeletonNeighbor(tail);
	const int numNodes = static_cast<int>(m_skelNodes.size() / 3);

	if (!ren || tail < 0 || tail >= numNodes) {
		HideSkeletonExtendPreview();
		return false;
	}

	double pTail[3] = {
		m_skelNodes[3 * tail],
		m_skelNodes[3 * tail + 1],
		m_skelNodes[3 * tail + 2]
	};
	double picked[3];
	if (!GetSkeletonExtensionTargetPoint(x, y, ren, tail, picked)) {
		HideSkeletonExtendPreview();
		return false;
	}

	double v[3] = { picked[0] - pTail[0], picked[1] - pTail[1], picked[2] - pTail[2] };
	double len = vtkMath::Normalize(v);
	const double meanLen = std::max(EstimateSkeletonMeanEdgeLength(), 1e-6);
	if (!std::isfinite(len) || len < meanLen * 0.10) {
		HideSkeletonExtendPreview();
		return false;
	}
	len = std::max(meanLen * 0.20, std::min(meanLen * 8.0, len));

	double oldDir[3] = { v[0], v[1], v[2] };
	bool hasTailDirection = false;
	if (neighbor >= 0 && neighbor < numNodes) {
		double pNei[3] = {
			m_skelNodes[3 * neighbor],
			m_skelNodes[3 * neighbor + 1],
			m_skelNodes[3 * neighbor + 2]
		};
		oldDir[0] = pTail[0] - pNei[0];
		oldDir[1] = pTail[1] - pNei[1];
		oldDir[2] = pTail[2] - pNei[2];
		hasTailDirection = vtkMath::Normalize(oldDir) > 1e-10;
	}

	const double dotOldNew = std::max(-1.0, std::min(1.0, vtkMath::Dot(oldDir, v)));
	const bool needBlendNode = hasTailDirection && dotOldNew < 0.94;

	// 预览使用一段短的过渡折线，而不是硬折角的一根线。
	double blendDir[3] = {
		oldDir[0] + v[0],
		oldDir[1] + v[1],
		oldDir[2] + v[2]
	};
	if (vtkMath::Normalize(blendDir) < 1e-10) {
		blendDir[0] = v[0]; blendDir[1] = v[1]; blendDir[2] = v[2];
	}

	double pMid[3] = {
		pTail[0] + blendDir[0] * len * 0.35,
		pTail[1] + blendDir[1] * len * 0.35,
		pTail[2] + blendDir[2] * len * 0.35
	};
	double pNew[3] = {
		pTail[0] + v[0] * len,
		pTail[1] + v[1] * len,
		pTail[2] + v[2] * len
	};

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(pTail);
	if (needBlendNode) {
		pts->InsertNextPoint(pMid);
	}
	pts->InsertNextPoint(pNew);

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	if (needBlendNode) {
		vtkIdType ids0[2] = { 0, 1 };
		vtkIdType ids1[2] = { 1, 2 };
		lines->InsertNextCell(2, ids0);
		lines->InsertNextCell(2, ids1);
	}
	else {
		vtkIdType ids[2] = { 0, 1 };
		lines->InsertNextCell(2, ids);
	}

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->SetPoints(pts);
	pd->SetLines(lines);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(pd);
	mapper->ScalarVisibilityOff();

	if (!m_skeletonExtendPreviewActor) {
		m_skeletonExtendPreviewActor = vtkSmartPointer<vtkActor>::New();
		m_skeletonExtendPreviewActor->SetPickable(0);
		m_skeletonExtendPreviewActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
		m_skeletonExtendPreviewActor->GetProperty()->SetLineWidth(4.0);
		m_skeletonExtendPreviewActor->GetProperty()->SetLighting(false);
		if (ren) ren->AddActor(m_skeletonExtendPreviewActor);
	}
	m_skeletonExtendPreviewActor->SetMapper(mapper);
	m_skeletonExtendPreviewActor->VisibilityOn();

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow()) {
		this->GetInteractor()->GetRenderWindow()->Render();
	}
	return true;
}

bool CageInteractorStyle::ExtendSkeletonFromSelectedTail(int x, int y, vtkRenderer* ren)
{
	const int tail = m_selectedExtendSkelNode;
	const int neighbor = GetSingleSkeletonNeighbor(tail);
	const int numNodes = static_cast<int>(m_skelNodes.size() / 3);

	if (!ren || tail < 0 || tail >= numNodes)
	{
		qDebug() << "[骨架延伸] 当前选择不是有效骨架节点，无法延伸/插入。";
		return false;
	}

	double pTail[3] = {
		m_skelNodes[3 * tail],
		m_skelNodes[3 * tail + 1],
		m_skelNodes[3 * tail + 2]
	};
	double picked[3];
	if (!GetSkeletonExtensionTargetPoint(x, y, ren, tail, picked)) {
		return false;
	}

	double targetDir[3] = { picked[0] - pTail[0], picked[1] - pTail[1], picked[2] - pTail[2] };
	double extendLen = vtkMath::Normalize(targetDir);
	const double meanLen = std::max(EstimateSkeletonMeanEdgeLength(), 1e-6);
	if (!std::isfinite(extendLen) || extendLen < meanLen * 0.10) {
		return false;
	}
	extendLen = std::max(meanLen * 0.20, std::min(meanLen * 8.0, extendLen));

	auto edgeExists = [&](int a, int b) -> bool
	{
		if (a < 0 || b < 0 || a == b) {
			return true;
		}
		for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2) {
			const int u = static_cast<int>(m_skelEdges[i]);
			const int v = static_cast<int>(m_skelEdges[i + 1]);
			if ((u == a && v == b) || (u == b && v == a)) {
				return true;
			}
		}
		return false;
	};

	auto addSkeletonEdge = [&](int a, int b) -> bool
	{
		if (a < 0 || b < 0 || a >= numNodes || b >= numNodes || edgeExists(a, b)) {
			return false;
		}
		m_skelEdges.push_back(static_cast<unsigned int>(a));
		m_skelEdges.push_back(static_cast<unsigned int>(b));
		return true;
	};

	const int clickedNode = PickSkeletonNode(x, y, ren);
	if (clickedNode >= 0 && clickedNode != tail) {
		if (!addSkeletonEdge(tail, clickedNode)) {
			return false;
		}
		RefreshSkeletonPolyData();
		if (m_editStraightSegments) {
			m_editStraightSkeletonNodes = m_skelNodes;
			RefreshStraightSegmentEditActor();
		}
		ClearSkeletonExtendSelection(ren);
		return true;
	}

	const int hitEdge = PickSkeletonEdge(x, y, ren);
	if (hitEdge >= 0 && 2 * hitEdge + 1 < static_cast<int>(m_skelEdges.size())) {
		const int edgeA = static_cast<int>(m_skelEdges[2 * hitEdge]);
		const int edgeB = static_cast<int>(m_skelEdges[2 * hitEdge + 1]);
		if (edgeA >= 0 && edgeB >= 0 &&
			edgeA < numNodes && edgeB < numNodes &&
			edgeA != edgeB && edgeA != tail && edgeB != tail)
		{
			double pA[3] = {
				m_skelNodes[3 * edgeA],
				m_skelNodes[3 * edgeA + 1],
				m_skelNodes[3 * edgeA + 2]
			};
			double pB[3] = {
				m_skelNodes[3 * edgeB],
				m_skelNodes[3 * edgeB + 1],
				m_skelNodes[3 * edgeB + 2]
			};
			double ab[3] = { pB[0] - pA[0], pB[1] - pA[1], pB[2] - pA[2] };
			const double abLen2 = vtkMath::Dot(ab, ab);
			if (abLen2 > 1e-12) {
				double ap[3] = { picked[0] - pA[0], picked[1] - pA[1], picked[2] - pA[2] };
				double t = vtkMath::Dot(ap, ab) / abLen2;
				t = std::max(0.0, std::min(1.0, t));

				const double snapT = 0.05;
				if (t <= snapT) {
					if (!addSkeletonEdge(tail, edgeA)) return false;
				}
				else if (t >= 1.0 - snapT) {
					if (!addSkeletonEdge(tail, edgeB)) return false;
				}
				else {
					double pSplit[3] = {
						pA[0] + ab[0] * t,
						pA[1] + ab[1] * t,
						pA[2] + ab[2] * t
					};
					const unsigned int splitNode = static_cast<unsigned int>(m_skelNodes.size() / 3);
					m_skelNodes.push_back(pSplit[0]);
					m_skelNodes.push_back(pSplit[1]);
					m_skelNodes.push_back(pSplit[2]);

					double splitRadius = tail < static_cast<int>(m_skelRadii.size()) ? m_skelRadii[tail] : 0.0;
					if (edgeA < static_cast<int>(m_skelRadii.size()) &&
						edgeB < static_cast<int>(m_skelRadii.size()))
					{
						splitRadius = (1.0 - t) * m_skelRadii[edgeA] + t * m_skelRadii[edgeB];
					}
					m_skelRadii.push_back(splitRadius);
					if (!m_skelSupportOffsets.empty()) m_skelSupportOffsets.push_back(m_skelSupportOffsets.back());

					m_skelEdges[2 * hitEdge] = static_cast<unsigned int>(edgeA);
					m_skelEdges[2 * hitEdge + 1] = splitNode;
					m_skelEdges.push_back(splitNode);
					m_skelEdges.push_back(static_cast<unsigned int>(edgeB));
					m_skelEdges.push_back(static_cast<unsigned int>(tail));
					m_skelEdges.push_back(splitNode);
				}

				RefreshSkeletonPolyData();
				if (m_editStraightSegments) {
					m_editStraightSkeletonNodes = m_skelNodes;
					RefreshStraightSegmentEditActor();
				}
				ClearSkeletonExtendSelection(ren);
				return true;
			}
		}
	}

	double oldDir[3] = { targetDir[0], targetDir[1], targetDir[2] };
	bool hasTailDirection = false;
	if (neighbor >= 0 && neighbor < numNodes) {
		double pNei[3] = {
			m_skelNodes[3 * neighbor],
			m_skelNodes[3 * neighbor + 1],
			m_skelNodes[3 * neighbor + 2]
		};
		oldDir[0] = pTail[0] - pNei[0];
		oldDir[1] = pTail[1] - pNei[1];
		oldDir[2] = pTail[2] - pNei[2];
		hasTailDirection = vtkMath::Normalize(oldDir) > 1e-10;
	}

	double blendDir[3] = {
		oldDir[0] + targetDir[0],
		oldDir[1] + targetDir[1],
		oldDir[2] + targetDir[2]
	};
	if (vtkMath::Normalize(blendDir) < 1e-10) {
		blendDir[0] = targetDir[0];
		blendDir[1] = targetDir[1];
		blendDir[2] = targetDir[2];
	}

	const double dotOldNew = std::max(-1.0, std::min(1.0, vtkMath::Dot(oldDir, targetDir)));
	const bool needBlendNode = hasTailDirection && (dotOldNew < 0.94);

	unsigned int lastNode = static_cast<unsigned int>(tail);
	int midNodeIndex = -1;

	if (needBlendNode) {
		double pMid[3] = {
			pTail[0] + blendDir[0] * extendLen * 0.35,
			pTail[1] + blendDir[1] * extendLen * 0.35,
			pTail[2] + blendDir[2] * extendLen * 0.35
		};
		const unsigned int midNode = static_cast<unsigned int>(m_skelNodes.size() / 3);
		midNodeIndex = static_cast<int>(midNode);
		m_skelNodes.push_back(pMid[0]);
		m_skelNodes.push_back(pMid[1]);
		m_skelNodes.push_back(pMid[2]);
		m_skelEdges.push_back(static_cast<unsigned int>(tail));
		m_skelEdges.push_back(midNode);
		m_skelRadii.push_back(tail < static_cast<int>(m_skelRadii.size()) ? m_skelRadii[tail] : 0.0);
		if (!m_skelSupportOffsets.empty()) m_skelSupportOffsets.push_back(m_skelSupportOffsets.back());
		lastNode = midNode;
	}

	double pNew[3] = {
		pTail[0] + targetDir[0] * extendLen,
		pTail[1] + targetDir[1] * extendLen,
		pTail[2] + targetDir[2] * extendLen
	};
	const unsigned int newNode = static_cast<unsigned int>(m_skelNodes.size() / 3);
	const int newNodeIndex = static_cast<int>(newNode);
	m_skelNodes.push_back(pNew[0]);
	m_skelNodes.push_back(pNew[1]);
	m_skelNodes.push_back(pNew[2]);
	m_skelEdges.push_back(lastNode);
	m_skelEdges.push_back(newNode);
	m_skelRadii.push_back(tail < static_cast<int>(m_skelRadii.size()) ? m_skelRadii[tail] : 0.0);
	if (!m_skelSupportOffsets.empty()) m_skelSupportOffsets.push_back(m_skelSupportOffsets.back());

	RefreshSkeletonPolyData();
	if (hasTailDirection) {
		SmoothSkeletonAroundNode(static_cast<int>(lastNode), 3);
	}

	// 方案 A：如果被延伸的尾端本来就是某条已识别直线段的端点，
	// 则把新增的过渡点/终点直接并入该直线段。
	// 这样手动延伸后，不需要重新依赖 StraightSkeletonSegmenter 的长度/角度阈值，
	// 新增尾段会继承原来的绿色直线段显示与后续建笼数据。
	if (hasTailDirection && m_editStraightSegments && !m_editStraightSegments->empty())
	{
		auto pointFromNode = [&](int nodeId, double out[3]) -> bool
		{
			const int n = static_cast<int>(m_skelNodes.size() / 3);
			if (nodeId < 0 || nodeId >= n) return false;
			out[0] = m_skelNodes[3 * nodeId];
			out[1] = m_skelNodes[3 * nodeId + 1];
			out[2] = m_skelNodes[3 * nodeId + 2];
			return true;
		};

		auto dist2 = [](const double a[3], const double b[3]) -> double
		{
			const double dx = a[0] - b[0];
			const double dy = a[1] - b[1];
			const double dz = a[2] - b[2];
			return dx * dx + dy * dy + dz * dz;
		};

		double currentTail[3];
		if (pointFromNode(tail, currentTail))
		{
			int bestSegment = -1;
			bool attachToBack = true;
			double bestD2 = std::numeric_limits<double>::max();

			for (int si = 0; si < static_cast<int>(m_editStraightSegments->size()); ++si)
			{
				StraightSkeletonSegmenter::Segment& seg = (*m_editStraightSegments)[si];
				if (!seg.valid || seg.points.size() < 6) continue;

				double front[3] = { seg.points[0], seg.points[1], seg.points[2] };
				double back[3] = {
					seg.points[seg.points.size() - 3],
					seg.points[seg.points.size() - 2],
					seg.points[seg.points.size() - 1]
				};

				const double dFront = dist2(front, currentTail);
				if (dFront < bestD2)
				{
					bestD2 = dFront;
					bestSegment = si;
					attachToBack = false;
				}

				const double dBack = dist2(back, currentTail);
				if (dBack < bestD2)
				{
					bestD2 = dBack;
					bestSegment = si;
					attachToBack = true;
				}
			}

			const double inheritTol = std::max(meanLen * 0.85, 1e-6);
			if (bestSegment >= 0 && bestD2 <= inheritTol * inheritTol)
			{
				StraightSkeletonSegmenter::Segment& seg = (*m_editStraightSegments)[bestSegment];

				double pTailNow[3];
				double pMidNow[3];
				double pNewNow[3];
				pointFromNode(tail, pTailNow);
				if (midNodeIndex >= 0) pointFromNode(midNodeIndex, pMidNow);
				pointFromNode(newNodeIndex, pNewNow);

				auto pushPointBack = [&](const double p[3])
				{
					seg.points.push_back(p[0]);
					seg.points.push_back(p[1]);
					seg.points.push_back(p[2]);
				};

				auto pushPointFront = [&](const double p[3])
				{
					seg.points.insert(seg.points.begin(), { p[0], p[1], p[2] });
				};

				if (attachToBack)
				{
					// 把直线段原尾点校正到当前骨架尾端，避免显示上出现小缝隙。
					seg.points[seg.points.size() - 3] = pTailNow[0];
					seg.points[seg.points.size() - 2] = pTailNow[1];
					seg.points[seg.points.size() - 1] = pTailNow[2];
					if (midNodeIndex >= 0) pushPointBack(pMidNow);
					pushPointBack(pNewNow);
				}
				else
				{
					// 直线段起点是当前尾端，则从前端继承延伸。
					seg.points[0] = pTailNow[0];
					seg.points[1] = pTailNow[1];
					seg.points[2] = pTailNow[2];
					if (midNodeIndex >= 0) pushPointFront(pMidNow);
					pushPointFront(pNewNow);
				}

				seg.originalPoints = seg.points;
				m_editStraightSkeletonNodes = m_skelNodes;
				RefreshStraightSegmentEditActor();

				qDebug() << "[骨架延伸] 新增尾段已继承到直线段" << bestSegment
					<< "新增节点" << newNodeIndex
					<< "过渡节点" << midNodeIndex;
			}
		}
	}

	ClearSkeletonExtendSelection(ren);
	return true;
}


void CageInteractorStyle::RefreshSkeletonPolyData()
{
	if (!m_skelPolyData)
	{
		m_skelPolyData = vtkSmartPointer<vtkPolyData>::New();
	}

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	const vtkIdType nodeCount = static_cast<vtkIdType>(m_skelNodes.size() / 3);
	pts->SetNumberOfPoints(nodeCount);
	for (vtkIdType i = 0; i < nodeCount; ++i)
	{
		pts->SetPoint(i,
			m_skelNodes[3 * i],
			m_skelNodes[3 * i + 1],
			m_skelNodes[3 * i + 2]);
	}

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2)
	{
		const unsigned int a = m_skelEdges[i];
		const unsigned int b = m_skelEdges[i + 1];
		if (a >= static_cast<unsigned int>(nodeCount) || b >= static_cast<unsigned int>(nodeCount) || a == b)
		{
			continue;
		}

		vtkIdType ids[2] = {
			static_cast<vtkIdType>(a),
			static_cast<vtkIdType>(b)
		};
		lines->InsertNextCell(2, ids);
	}

	m_skelPolyData->SetPoints(pts);
	m_skelPolyData->SetLines(lines);
	m_skelPolyData->Modified();

	if (this->GetInteractor() && this->GetInteractor()->GetRenderWindow())
	{
		this->GetInteractor()->GetRenderWindow()->Render();
	}
}

int CageInteractorStyle::PickSkeletonEdge(int x, int y, vtkRenderer* ren) const
{
	if (!ren || m_skelNodes.empty() || m_skelEdges.empty())
	{
		return -1;
	}

	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
	const int edgeCount = static_cast<int>(m_skelEdges.size() / 2);

	auto worldToDisplay = [&](int nodeId, double out[2]) -> bool
	{
		if (nodeId < 0 || nodeId >= nodeCount)
		{
			return false;
		}

		ren->SetWorldPoint(
			m_skelNodes[3 * nodeId],
			m_skelNodes[3 * nodeId + 1],
			m_skelNodes[3 * nodeId + 2],
			1.0);
		ren->WorldToDisplay();
		double* d = ren->GetDisplayPoint();
		out[0] = d[0];
		out[1] = d[1];
		return true;
	};

	int bestEdge = -1;
	double bestDist = 12.0; // 像素阈值，避免误删远处骨架边。

	for (int ei = 0; ei < edgeCount; ++ei)
	{
		const int a = static_cast<int>(m_skelEdges[2 * ei]);
		const int b = static_cast<int>(m_skelEdges[2 * ei + 1]);
		double pa[2], pb[2];
		if (!worldToDisplay(a, pa) || !worldToDisplay(b, pb))
		{
			continue;
		}

		const double vx = pb[0] - pa[0];
		const double vy = pb[1] - pa[1];
		const double wx = static_cast<double>(x) - pa[0];
		const double wy = static_cast<double>(y) - pa[1];
		const double len2 = vx * vx + vy * vy;
		if (len2 < 1e-12)
		{
			continue;
		}

		double t = (wx * vx + wy * vy) / len2;
		t = std::max(0.0, std::min(1.0, t));
		const double cx = pa[0] + t * vx;
		const double cy = pa[1] + t * vy;
		const double dx = static_cast<double>(x) - cx;
		const double dy = static_cast<double>(y) - cy;
		const double dist = std::sqrt(dx * dx + dy * dy);

		if (dist < bestDist)
		{
			bestDist = dist;
			bestEdge = ei;
		}
	}

	return bestEdge;
}

int CageInteractorStyle::GetSkeletonNodeDegree(int nodeIndex) const
{
	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
	if (nodeIndex < 0 || nodeIndex >= nodeCount)
	{
		return 0;
	}

	int degree = 0;
	for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2)
	{
		const int a = static_cast<int>(m_skelEdges[i]);
		const int b = static_cast<int>(m_skelEdges[i + 1]);
		if (a == nodeIndex || b == nodeIndex)
		{
			++degree;
		}
	}
	return degree;
}

int CageInteractorStyle::GetSingleSkeletonNeighbor(int nodeIndex) const
{
	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
	if (nodeIndex < 0 || nodeIndex >= nodeCount)
	{
		return -1;
	}

	int neighbor = -1;
	int degree = 0;
	for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2)
	{
		const int a = static_cast<int>(m_skelEdges[i]);
		const int b = static_cast<int>(m_skelEdges[i + 1]);
		if (a == nodeIndex)
		{
			neighbor = b;
			++degree;
		}
		else if (b == nodeIndex)
		{
			neighbor = a;
			++degree;
		}
	}

	return degree == 1 ? neighbor : -1;
}

double CageInteractorStyle::EstimateSkeletonMeanEdgeLength() const
{
	if (m_skelNodes.empty() || m_skelEdges.empty())
	{
		return 1.0;
	}

	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
	double sum = 0.0;
	int count = 0;

	for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2)
	{
		const int a = static_cast<int>(m_skelEdges[i]);
		const int b = static_cast<int>(m_skelEdges[i + 1]);
		if (a < 0 || b < 0 || a >= nodeCount || b >= nodeCount || a == b)
		{
			continue;
		}

		double pa[3] = {
			m_skelNodes[3 * a],
			m_skelNodes[3 * a + 1],
			m_skelNodes[3 * a + 2]
		};
		double pb[3] = {
			m_skelNodes[3 * b],
			m_skelNodes[3 * b + 1],
			m_skelNodes[3 * b + 2]
		};

		sum += std::sqrt(vtkMath::Distance2BetweenPoints(pa, pb));
		++count;
	}

	return count > 0 ? sum / static_cast<double>(count) : 1.0;
}

int CageInteractorStyle::PickSkeletonNode(int x, int y, vtkRenderer* ren)
{
	if (m_skelNodes.empty()) return -1;

	int bestNode = -1;
	double minDist = 15.0; // 

	for (size_t i = 0; i < m_skelNodes.size() / 3; i++) {
		ren->SetWorldPoint(m_skelNodes[3 * i], m_skelNodes[3 * i + 1], m_skelNodes[3 * i + 2], 1.0);
		ren->WorldToDisplay();
		double* display = ren->GetDisplayPoint();

		// 
		double dist = std::sqrt(std::pow(display[0] - x, 2) + std::pow(display[1] - y, 2));
		if (dist < minDist) {
			minDist = dist;
			bestNode = (int)i;
		}
	}
	return bestNode;
}
void CageInteractorStyle::UpdateSkeletonHighlight(vtkRenderer* ren) {
	if (!m_skelNodeHighlightActor) {
		vtkNew<vtkPolyDataMapper> mapper;
		m_skelNodeHighlightActor = vtkSmartPointer<vtkActor>::New();
		m_skelNodeHighlightActor->SetMapper(mapper);
		m_skelNodeHighlightActor->GetProperty()->SetColor(1.0, 1.0, 0.0); // 
		m_skelNodeHighlightActor->GetProperty()->SetLighting(false);      // 
		m_skelNodeHighlightActor->SetPickable(0);
		ren->AddActor(m_skelNodeHighlightActor);
	}

	if (m_selectedSkelNode1 == -1) {
		m_skelNodeHighlightActor->VisibilityOff();
	}
	else {

		vtkNew<vtkSphereSource> sphere;
		sphere->SetCenter(m_skelNodes[3 * m_selectedSkelNode1],
			m_skelNodes[3 * m_selectedSkelNode1 + 1],
			m_skelNodes[3 * m_selectedSkelNode1 + 2]);


		sphere->SetRadius(1.0);

		sphere->SetThetaResolution(16);
		sphere->SetPhiResolution(16);
		sphere->Update();

		vtkPolyDataMapper::SafeDownCast(m_skelNodeHighlightActor->GetMapper())->SetInputData(sphere->GetOutput());
		m_skelNodeHighlightActor->VisibilityOn();
	}
	this->GetInteractor()->GetRenderWindow()->Render();
}
void CageInteractorStyle::MergeSkeletonNodes(int n1, int n2) {
	if (!m_skelPolyData) return;


	int numNodes = m_skelNodes.size() / 3;
	std::vector<std::vector<int>> adjList(numNodes);
	for (size_t i = 0; i < m_skelEdges.size() / 2; i++) {
		int u = m_skelEdges[2 * i];
		int v = m_skelEdges[2 * i + 1];
		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}


	std::vector<int> path = FindPathBFS(n1, n2, adjList);
	std::set<int> pathNodes;
	for (int node : path) {
		pathNodes.insert(node);
	}


	if (pathNodes.empty()) {
		pathNodes.insert(n1);
		pathNodes.insert(n2);
	}

	// 2. 
	double p_new[3];
	for (int d = 0; d < 3; d++) {
		p_new[d] = (m_skelNodes[3 * n1 + d] + m_skelNodes[3 * n2 + d]) / 2.0;
		m_skelNodes[3 * n1 + d] = p_new[d];
	}

	// 3.
	std::vector<unsigned int> newEdges;
	std::set<std::pair<unsigned int, unsigned int>> uniqueEdges;

	for (size_t i = 0; i < m_skelEdges.size() / 2; i++) {
		unsigned int u = m_skelEdges[2 * i];
		unsigned int v = m_skelEdges[2 * i + 1];

		// 
		if (pathNodes.find(u) != pathNodes.end()) u = n1;
		if (pathNodes.find(v) != pathNodes.end()) v = n1;

		// 
		if (u == v) continue;

		// 
		unsigned int minNode = std::min(u, v);
		unsigned int maxNode = std::max(u, v);
		if (uniqueEdges.find({ minNode, maxNode }) == uniqueEdges.end()) {
			uniqueEdges.insert({ minNode, maxNode });
			newEdges.push_back(u);
			newEdges.push_back(v);
		}
	}

	// 
	m_skelEdges = newEdges;

	// 4. 
	m_skelPolyData->GetPoints()->SetPoint(n1, p_new);
	m_skelPolyData->GetPoints()->Modified(); // 

	// 
	vtkNew<vtkCellArray> newLines;
	for (size_t i = 0; i < m_skelEdges.size() / 2; i++) {
		vtkIdType lineIds[2] = { (vtkIdType)m_skelEdges[2 * i], (vtkIdType)m_skelEdges[2 * i + 1] };
		newLines->InsertNextCell(2, lineIds);
	}
	m_skelPolyData->SetLines(newLines);
	SmoothSkeleton(10);
}
void CageInteractorStyle::SmoothSkeleton(int iterations) {
	if (m_skelNodes.empty() || !m_skelPolyData) return;

	int numNodes = m_skelNodes.size() / 3;

	// 1.
	std::vector<std::vector<int>> adjList(numNodes);
	for (size_t i = 0; i < m_skelEdges.size() / 2; i++) {
		int u = m_skelEdges[2 * i];
		int v = m_skelEdges[2 * i + 1];
		if (u != v) {
			adjList[u].push_back(v);
			adjList[v].push_back(u);
		}
	}

	// 2. 
	std::vector<double> newNodes = m_skelNodes;
	double lambda = 0.5; // 

	for (int iter = 0; iter < iterations; ++iter) {
		std::vector<double> tempNodes = newNodes;

		for (int i = 0; i < numNodes; ++i) {
			// 
			if (adjList[i].size() == 2) {
				int n1 = adjList[i][0];
				int n2 = adjList[i][1];

				for (int d = 0; d < 3; ++d) {
					double avg = (tempNodes[3 * n1 + d] + tempNodes[3 * n2 + d]) / 2.0;
					newNodes[3 * i + d] = tempNodes[3 * i + d] + lambda * (avg - tempNodes[3 * i + d]);
				}
			}
		}
	}

	// 3. 
	m_skelNodes = newNodes;
	for (int i = 0; i < numNodes; ++i) {
		m_skelPolyData->GetPoints()->SetPoint(i, &m_skelNodes[3 * i]);
	}

	m_skelPolyData->GetPoints()->Modified(); // 
	m_skelPolyData->Modified();              // 
	this->GetInteractor()->GetRenderWindow()->Render();
}

void CageInteractorStyle::SmoothSkeletonAroundNode(int centerNode, int iterations)
{
	const int numNodes = static_cast<int>(m_skelNodes.size() / 3);
	if (centerNode < 0 || centerNode >= numNodes || iterations <= 0) {
		return;
	}

	std::vector<std::vector<int>> adjList(numNodes);
	for (size_t i = 0; i + 1 < m_skelEdges.size(); i += 2) {
		int u = static_cast<int>(m_skelEdges[i]);
		int v = static_cast<int>(m_skelEdges[i + 1]);
		if (u < 0 || v < 0 || u >= numNodes || v >= numNodes || u == v) continue;
		adjList[u].push_back(v);
		adjList[v].push_back(u);
	}

	std::set<int> localNodes;
	localNodes.insert(centerNode);
	for (int n : adjList[centerNode]) {
		localNodes.insert(n);
		for (int n2 : adjList[n]) {
			localNodes.insert(n2);
		}
	}

	std::vector<double> newNodes = m_skelNodes;
	const double lambda = 0.35;
	for (int iter = 0; iter < iterations; ++iter) {
		std::vector<double> temp = newNodes;
		for (int id : localNodes) {
			if (id <= 0 || id >= numNodes) continue;
			if (adjList[id].size() != 2) continue;
			const int a = adjList[id][0];
			const int b = adjList[id][1];
			for (int d = 0; d < 3; ++d) {
				double avg = 0.5 * (temp[3 * a + d] + temp[3 * b + d]);
				newNodes[3 * id + d] = temp[3 * id + d] + lambda * (avg - temp[3 * id + d]);
			}
		}
	}

	m_skelNodes.swap(newNodes);
	RefreshSkeletonPolyData();
}


static void DebugBridgeCageTopology(
	int bridgeIndex,
	int originalTubeIndex,
	bool tubeIsStart,
	const double tubeCanonical[4][3],
	const double hubCanonical[4][3],
	const CageSegment& bridge,
	const double normalUp[3])
{
	auto PrintPoint = [](const char* name, int i, const double p[3])
	{
		if (false) qDebug().nospace()
			<< name << "[" << i << "] = ("
			<< p[0] << ", "
			<< p[1] << ", "
			<< p[2] << ")";
	};

	auto Vec = [](const double a[3], const double b[3], double out[3])
	{
		out[0] = b[0] - a[0];
		out[1] = b[1] - a[1];
		out[2] = b[2] - a[2];
	};

	auto DotNormalized = [&](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double
	{
		double va[3], vb[3];

		Vec(a0, a1, va);
		Vec(b0, b1, vb);

		double na = vtkMath::Norm(va);
		double nb = vtkMath::Norm(vb);

		if (na < 1e-8 || nb < 1e-8) return 0.0;

		vtkMath::Normalize(va);
		vtkMath::Normalize(vb);

		return vtkMath::Dot(va, vb);
	};

	auto Distance = [](const double a[3], const double b[3]) -> double
	{
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto SignedArea = [](const double ring[4][3], const double n[3]) -> double
	{
		double area[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double c[3];
			vtkMath::Cross(ring[i], ring[j], c);

			area[0] += c[0];
			area[1] += c[1];
			area[2] += c[2];
		}

		return vtkMath::Dot(area, n);
	};

	auto FindNearestTubeIndex = [&](const double p[3]) -> int
	{
		int bestId = -1;
		double bestD2 = 1e100;

		for (int i = 0; i < 4; ++i) {
			double d2 = vtkMath::Distance2BetweenPoints(p, tubeCanonical[i]);

			if (d2 < bestD2) {
				bestD2 = d2;
				bestId = i;
			}
		}

		return bestId;
	};

	auto FindNearestHubIndex = [&](const double p[3]) -> int
	{
		int bestId = -1;
		double bestD2 = 1e100;

		for (int i = 0; i < 4; ++i) {
			double d2 = vtkMath::Distance2BetweenPoints(p, hubCanonical[i]);

			if (d2 < bestD2) {
				bestD2 = d2;
				bestId = i;
			}
		}

		return bestId;
	};

	if (false) qDebug() << "\n\n==================== [桥接调试开始] ====================";
	if (false) qDebug() << "桥接编号 =" << bridgeIndex
		<< "原始管口序号 =" << originalTubeIndex
		<< "是否起始端 =" << tubeIsStart;

	if (false) qDebug() << "上方向 normalUp =" << normalUp[0] << normalUp[1] << normalUp[2];

	if (false) qDebug() << "\n--- tubeCanonical  ---";
	for (int i = 0; i < 4; ++i) {
		PrintPoint("tubeCanonical", i, tubeCanonical[i]);
	}

	if (false) qDebug() << "\n--- hubCanonical 接头侧口四点，理论语义 H0,H1,H2,H3 ---";
	if (false) qDebug() << "H0 = 上侧(前一个)";
	if (false) qDebug() << "H1 = 下侧(前一个)";
	if (false) qDebug() << "H2 = 下侧(当前)";
	if (false) qDebug() << "H3 = 上侧(当前)";
	for (int i = 0; i < 4; ++i) {
		PrintPoint("hubCanonical", i, hubCanonical[i]);
	}

	if (false) qDebug() << "\n--- bridge.pStart  ---";
	for (int i = 0; i < 4; ++i) {
		PrintPoint("bridge.pStart", i, bridge.pStart[i]);
		if (false) qDebug() << "    最近的接头侧标准点序号 =" << FindNearestHubIndex(bridge.pStart[i]);
	}

	if (false) qDebug() << "\n--- bridge.pEnd ---";
	for (int i = 0; i < 4; ++i) {
		PrintPoint("bridge.pEnd", i, bridge.pEnd[i]);
		if (false) qDebug() << "    最近的管口标准点序号 =" << FindNearestTubeIndex(bridge.pEnd[i]);
	}

	if (false) qDebug() << "\n--- bridge pStart[i] -> pEnd[i] ---";
	double totalConnLen = 0.0;

	for (int i = 0; i < 4; ++i) {
		double len = Distance(bridge.pStart[i], bridge.pEnd[i]);
		totalConnLen += len;

		if (false) qDebug().nospace()
			<< "  i=" << i
			<< " : hubNearest=" << FindNearestHubIndex(bridge.pStart[i])
			<< " -> tubeNearest=" << FindNearestTubeIndex(bridge.pEnd[i])
			<< " , length=" << len;
	}

	if (false) qDebug() << "  总连接长度 =" << totalConnLen;

	if (false) qDebug() << "\n--- ring signed area 閻╃顕?normalUp ---";
	double tubeArea = SignedArea(tubeCanonical, normalUp);

	double bridgeStartArea = SignedArea(bridge.pStart, normalUp);
	double bridgeEndArea = SignedArea(bridge.pEnd, normalUp);
	double hubArea = SignedArea(hubCanonical, normalUp);

	if (false) qDebug() << "  管口标准环面积 =" << tubeArea;
	if (false) qDebug() << "  接头标准环面积 =" << hubArea;
	if (false) qDebug() << "  桥接起点环面积 =" << bridgeStartArea;
	if (false) qDebug() << "  桥接终点环面积 =" << bridgeEndArea;

	if (false) qDebug() << "  桥接起点环与接头环方向是否一致？"
		<< ((bridgeStartArea * hubArea >= 0.0) ? "YES" : "NO");

	if (false) qDebug() << "  桥接终点环与管口环方向是否一致？"
		<< ((bridgeEndArea * tubeArea >= 0.0) ? "YES" : "NO");

	if (false) qDebug() << "  F0 = pStart0, pStart1, pEnd1, pEnd0";
	if (false) qDebug() << "  F1 = pStart1, pStart2, pEnd2, pEnd1";
	if (false) qDebug() << "  F2 = pStart2, pStart3, pEnd3, pEnd2";
	if (false) qDebug() << "  F3 = pStart3, pStart0, pEnd0, pEnd3";

	for (int i = 0; i < 4; ++i)
	{
		int j = (i + 1) % 4;

		// 
		double edgeDot = DotNormalized(
			bridge.pStart[i],
			bridge.pStart[j],
			bridge.pEnd[i],
			bridge.pEnd[j]
		);

		// 
		double len0 = Distance(bridge.pStart[i], bridge.pEnd[i]);
		double len1 = Distance(bridge.pStart[j], bridge.pEnd[j]);

		// 
		double diagA = Distance(bridge.pStart[i], bridge.pEnd[j]);
		double diagB = Distance(bridge.pStart[j], bridge.pEnd[i]);

		if (false) qDebug().nospace()
			<< "  F" << i
			<< " edgeDot(startEdge,endEdge)=" << edgeDot
			<< "  len_i=" << len0
			<< "  len_j=" << len1
			<< "  diagA=" << diagA
			<< "  diagB=" << diagB;

		if (edgeDot < -0.3) {
			if (false) qDebug() << "    !!! WARNING: 璇ヤ晶闈㈢殑涓ょ妯悜杈硅繎浼煎弽鍚戯紝bridge 鍙兘鍦ㄨ繖閲屾嫥楹昏姳";
		}

		if (std::abs(diagA - diagB) > std::max(len0, len1) * 1.5) {
			if (false) qDebug() << "[桥接调试] 侧面对角线长度差异较大";
		}
	}

	if (false) qDebug() << "\n--- 接头上下共享边方向检查 ---";
	if (false) qDebug() << "接头上表面中，此口相关边通常是 H0 -> H3";
	if (false) qDebug() << "桥接起点环中，对应上边是 pStart3 -> pStart0";
	if (false) qDebug() << "接头上边应与桥接上边方向相反，dot 越接近 -1 越好";

	double hubTopVsBridgeTop = DotNormalized(
		hubCanonical[0],
		hubCanonical[3],
		bridge.pStart[3],
		bridge.pStart[0]
	);

	if (false) qDebug() << "  dot( hub H0->H3 , bridge S3->S0 ) =" << hubTopVsBridgeTop;

	if (false) qDebug() << "接头下表面中，此口相关边通常是 H2 -> H1";
	if (false) qDebug() << "桥接起点环中，对应下边是 pStart1 -> pStart2";
	if (false) qDebug() << "接头下边应与桥接下边方向相反，dot 越接近 -1 越好";

	double hubBottomVsBridgeBottom = DotNormalized(
		hubCanonical[2],
		hubCanonical[1],
		bridge.pStart[1],
		bridge.pStart[2]
	);

	if (false) qDebug() << "  dot( hub H2->H1 , bridge S1->S2 ) =" << hubBottomVsBridgeBottom;

	if (hubTopVsBridgeTop > 0.3 || hubBottomVsBridgeBottom > 0.3) {
		if (false) qDebug() << "[桥接调试] 接头和桥接的共享边方向可能相同";
	}

	if (false) qDebug() << "===================== [桥接调试结束] =====================\n\n";
}

bool CageInteractorStyle::GenerateNBranchFromRings(
	const std::vector<NBranchInputPort>& inputPorts,
	const NBranchBuildOptions& options)
{
	const int N = static_cast<int>(inputPorts.size());
	if (N < 3) {
		if (false) qDebug() << "[GenerateNBranchFromRings] input port count < 3.";
		return false;
	}

	auto CopyRing = [](double dst[4][3], const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	};

	auto RingCenter = [](const double ring[4][3], double c[3])
	{
		c[0] = c[1] = c[2] = 0.0;
		for (int i = 0; i < 4; ++i) {
			c[0] += ring[i][0] * 0.25;
			c[1] += ring[i][1] * 0.25;
			c[2] += ring[i][2] * 0.25;
		}
	};

	auto Distance = [](const double a[3], const double b[3]) -> double
	{
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto NormalizeSafe = [](double v[3], const double fallback[3]) -> bool
	{
		if (vtkMath::Norm(v) < 1e-8) {
			v[0] = fallback[0];
			v[1] = fallback[1];
			v[2] = fallback[2];
			if (vtkMath::Norm(v) < 1e-8) {
				v[0] = 1.0;
				v[1] = 0.0;
				v[2] = 0.0;
			}
			vtkMath::Normalize(v);
			return false;
		}
		vtkMath::Normalize(v);
		return true;
	};

	auto DotEdgeNormalized = [](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double
	{
		double va[3] = { a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2] };
		double vb[3] = { b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2] };
		double na = vtkMath::Norm(va);
		double nb = vtkMath::Norm(vb);
		if (na < 1e-10 || nb < 1e-10) return 0.0;
		va[0] /= na; va[1] /= na; va[2] /= na;
		vb[0] /= nb; vb[1] /= nb; vb[2] /= nb;
		return vtkMath::Dot(va, vb);
	};

	struct BridgeScore
	{
		double minEdgeDot;
		double avgEdgeDot;
		double totalLength;
		double maxDiagDiff;
	};

	auto EvaluateBridge = [&](const double startRing[4][3], const double endRing[4][3]) -> BridgeScore
	{
		BridgeScore score;
		score.minEdgeDot = 1.0;
		score.avgEdgeDot = 0.0;
		score.totalLength = 0.0;
		score.maxDiagDiff = 0.0;

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;
			double dot = DotEdgeNormalized(startRing[i], startRing[j], endRing[i], endRing[j]);
			score.minEdgeDot = std::min(score.minEdgeDot, dot);
			score.avgEdgeDot += dot * 0.25;
			score.totalLength += Distance(startRing[i], endRing[i]);
			double diagA = Distance(startRing[i], endRing[j]);
			double diagB = Distance(startRing[j], endRing[i]);
			score.maxDiagDiff = std::max(score.maxDiagDiff, std::abs(diagA - diagB));
		}

		return score;
	};

	auto IsBetterScore = [](const BridgeScore& a, const BridgeScore& b) -> bool
	{
		if (a.minEdgeDot > b.minEdgeDot + 1e-6) return true;
		if (a.minEdgeDot < b.minEdgeDot - 1e-6) return false;
		if (a.avgEdgeDot > b.avgEdgeDot + 1e-6) return true;
		if (a.avgEdgeDot < b.avgEdgeDot - 1e-6) return false;
		if (a.maxDiagDiff < b.maxDiagDiff - 1e-6) return true;
		if (a.maxDiagDiff > b.maxDiagDiff + 1e-6) return false;
		return a.totalLength < b.totalLength;
	};

	auto BuildTubeCandidateRing = [](const double tubeCanonical[4][3], bool reverse, int shift, double outRing[4][3])
	{
		static const int forwardBase[4] = { 0, 1, 2, 3 };
		static const int reverseBase[4] = { 0, 3, 2, 1 };
		const int* base = reverse ? reverseBase : forwardBase;
		for (int i = 0; i < 4; ++i) {
			int id = base[(i + shift) % 4];
			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = tubeCanonical[id][d];
			}
		}
	};

	struct TubeData
	{
		int index = -1;
		bool isVirtual = false;
		double center[3] = { 0.0, 0.0, 0.0 };
		double dirOut[3] = { 1.0, 0.0, 0.0 };
		double angle = 0.0;
		double w = 0.0;
		double h = 0.0;
		double distToHub = 0.0;
		double capRing[4][3];
	};

	std::vector<TubeData> tubes(N);
	double geomCenter[3] = { 0.0, 0.0, 0.0 };
	for (int k = 0; k < N; ++k) {
		tubes[k].index = k;
		tubes[k].isVirtual = inputPorts[k].isVirtual;
		CopyRing(tubes[k].capRing, inputPorts[k].ring);
		RingCenter(tubes[k].capRing, tubes[k].center);
		tubes[k].w = Distance(tubes[k].capRing[0], tubes[k].capRing[1]);
		tubes[k].h = Distance(tubes[k].capRing[1], tubes[k].capRing[2]);
		for (int d = 0; d < 3; ++d) {
			geomCenter[d] += tubes[k].center[d] / static_cast<double>(N);
		}
	}

	double hubCenter[3] = { geomCenter[0], geomCenter[1], geomCenter[2] };
	if (options.forceHubCenter) {
		hubCenter[0] = options.hubCenter[0];
		hubCenter[1] = options.hubCenter[1];
		hubCenter[2] = options.hubCenter[2];
	}

	double maxPortDim = 0.0;
	double minCenterDist = 1e100;
	for (int k = 0; k < N; ++k) {
		double dir[3] = {
			tubes[k].center[0] - hubCenter[0],
			tubes[k].center[1] - hubCenter[1],
			tubes[k].center[2] - hubCenter[2]
		};
		tubes[k].distToHub = vtkMath::Norm(dir);
		if (tubes[k].distToHub < 1e-8 && inputPorts[k].hasPreferredDir) {
			dir[0] = inputPorts[k].preferredDir[0];
			dir[1] = inputPorts[k].preferredDir[1];
			dir[2] = inputPorts[k].preferredDir[2];
		}
		NormalizeSafe(dir, inputPorts[k].hasPreferredDir ? inputPorts[k].preferredDir : tubes[k].capRing[1]);
		for (int d = 0; d < 3; ++d) tubes[k].dirOut[d] = dir[d];
		minCenterDist = std::min(minCenterDist, std::max(tubes[k].distToHub, 1e-6));
		maxPortDim = std::max(maxPortDim, std::max(tubes[k].w, tubes[k].h));
	}

	double normalUp[3] = { 0.0, 0.0, 1.0 };
	if (options.forceNormalUp) {
		normalUp[0] = options.normalUp[0];
		normalUp[1] = options.normalUp[1];
		normalUp[2] = options.normalUp[2];
		NormalizeSafe(normalUp, tubes[0].capRing[1]);
	}
	else {
		double bestCrossNorm = 0.0;
		for (int i = 0; i < N; ++i) {
			for (int j = i + 1; j < N; ++j) {
				double c[3];
				vtkMath::Cross(tubes[i].dirOut, tubes[j].dirOut, c);
				double n = vtkMath::Norm(c);
				if (n > bestCrossNorm) {
					bestCrossNorm = n;
					normalUp[0] = c[0];
					normalUp[1] = c[1];
					normalUp[2] = c[2];
				}
			}
		}
		NormalizeSafe(normalUp, tubes[0].capRing[1]);
	}

	double refAxis[3] = { tubes[0].dirOut[0], tubes[0].dirOut[1], tubes[0].dirOut[2] };
	double refSide[3];
	vtkMath::Cross(normalUp, refAxis, refSide);
	NormalizeSafe(refSide, tubes[0].capRing[1]);
	vtkMath::Cross(refSide, normalUp, refAxis);
	NormalizeSafe(refAxis, tubes[0].dirOut);

	for (TubeData& t : tubes) {
		t.angle = std::atan2(vtkMath::Dot(t.dirOut, refSide), vtkMath::Dot(t.dirOut, refAxis));
	}

	std::sort(tubes.begin(), tubes.end(), [](const TubeData& a, const TubeData& b) {
		return a.angle < b.angle;
	});

	if (minCenterDist < 1e-6 || minCenterDist > 1e90) {
		minCenterDist = std::max(maxPortDim, 1.0);
	}

	vtkSmartPointer<vtkPolyData> hubPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

	for (int k = 0; k < N; ++k) {
		int nextK = (k + 1) % N;
		double a1 = tubes[k].angle;
		double a2 = tubes[nextK].angle;
		if (a2 <= a1) a2 += 2.0 * vtkMath::Pi();
		double midAngle = 0.5 * (a1 + a2);

		double localW = 0.5 * (tubes[k].w + tubes[nextK].w);
		double localH = 0.5 * (tubes[k].h + tubes[nextK].h);
		double localR = std::min(localW * options.radiusScale, minCenterDist * 0.90);
		localR = std::max(localR, std::max(localW * options.minRadiusScale, maxPortDim * 0.020));

		double px = localR * std::cos(midAngle);
		double py = localR * std::sin(midAngle);
		double base[3] = {
			hubCenter[0] + px * refAxis[0] + py * refSide[0],
			hubCenter[1] + px * refAxis[1] + py * refSide[1],
			hubCenter[2] + px * refAxis[2] + py * refSide[2]
		};

		pts->InsertNextPoint(
			base[0] + normalUp[0] * localH * 0.5,
			base[1] + normalUp[1] * localH * 0.5,
			base[2] + normalUp[2] * localH * 0.5);
		pts->InsertNextPoint(
			base[0] - normalUp[0] * localH * 0.5,
			base[1] - normalUp[1] * localH * 0.5,
			base[2] - normalUp[2] * localH * 0.5);
	}

	polys->InsertNextCell(N);
	for (int k = 0; k < N; ++k) polys->InsertCellPoint(k * 2);
	polys->InsertNextCell(N);
	for (int k = N - 1; k >= 0; --k) polys->InsertCellPoint(k * 2 + 1);

	hubPoly->SetPoints(pts);
	hubPoly->SetPolys(polys);
	m_hubCages.push_back(hubPoly);

	// 重要：保持原先避让机制：只处理靠近 1 端面(targetInputIndex)的 hub 接口面，
	// 沿主骨架方向朝 2 分支让位，用来避让侧枝骨架和对应端面的交点。
	// 不再整体改整个 hubPoly；否则会把靠近 2 分支的一侧 hub 节点也拉乱。
	// 在 ApplyNBranchHubRingRelaxInPlace 内部再用 B0/B1 两个端面的真实 8 个节点做硬越界限制。
	if (g_NBranchHubRingRelax.enabled) {
		int targetSortedK = -1;
		for (int kk = 0; kk < N; ++kk) {
			if (tubes[kk].index == g_NBranchHubRingRelax.targetInputIndex) {
				targetSortedK = kk;
				break;
			}
		}

		if (targetSortedK >= 0) {
			int prevK = (targetSortedK - 1 + N) % N;
			double relaxedHubRing[4][3];
			pts->GetPoint(2 * prevK, relaxedHubRing[0]);
			pts->GetPoint(2 * prevK + 1, relaxedHubRing[1]);
			pts->GetPoint(2 * targetSortedK + 1, relaxedHubRing[2]);
			pts->GetPoint(2 * targetSortedK, relaxedHubRing[3]);

			ApplyNBranchHubRingRelaxInPlace(relaxedHubRing, tubes[targetSortedK].index);

			pts->SetPoint(2 * prevK, relaxedHubRing[0]);
			pts->SetPoint(2 * prevK + 1, relaxedHubRing[1]);
			pts->SetPoint(2 * targetSortedK + 1, relaxedHubRing[2]);
			pts->SetPoint(2 * targetSortedK, relaxedHubRing[3]);

			// 关键补充：A 侧接口四边形避让之后，还必须检查 hub 中靠近 B0/B1 的剩余两个节点。
			// 之前只检查 relaxedHubRing 的 4 个点，所以真正贴近 B 端面的两个 hub 节点没有被限制。
			ClampNBranchHubBsideVerticesAgainstBInPlace(pts, targetSortedK);

			pts->Modified();
			hubPoly->Modified();

			if (false) qDebug() << "[NBranchHubRingRelax] pre-applied targeted A-side concession and B-side two-node clamp"
				<< "targetSortedK =" << targetSortedK
				<< "targetInputIndex =" << g_NBranchHubRingRelax.targetInputIndex;
		}
	}

	for (int k = 0; k < N; ++k) {
		int prevK = (k - 1 + N) % N;
		double hubRing[4][3];
		pts->GetPoint(2 * prevK, hubRing[0]);
		pts->GetPoint(2 * prevK + 1, hubRing[1]);
		pts->GetPoint(2 * k + 1, hubRing[2]);
		pts->GetPoint(2 * k, hubRing[3]);


		double bestTubeRing[4][3];
		BridgeScore bestScore;
		bestScore.minEdgeDot = -1e100;
		bestScore.avgEdgeDot = -1e100;
		bestScore.totalLength = 1e100;
		bestScore.maxDiagDiff = 1e100;

		for (int family = 0; family < 2; ++family) {
			bool reverse = family == 1;
			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				BuildTubeCandidateRing(tubes[k].capRing, reverse, shift, candidate);
				BridgeScore score = EvaluateBridge(hubRing, candidate);
				if (IsBetterScore(score, bestScore)) {
					bestScore = score;
					CopyRing(bestTubeRing, candidate);
				}
			}
		}

		CageSegment bridge;
		CopyRing(bridge.pStart, hubRing);
		CopyRing(bridge.pEnd, bestTubeRing);
		m_allCages.push_back(bridge);

		if (options.debugLog && bestScore.minEdgeDot < 0.05) {
			if (false) qDebug() << "[GenerateNBranchFromRings] bridge twist risk"
				<< "sortedIndex =" << k
				<< "sourceIndex =" << tubes[k].index
				<< "isVirtual =" << tubes[k].isVirtual
				<< "minEdgeDot =" << bestScore.minEdgeDot;
		}
	}

	if (options.addSdf) {
		HubSDFInfo sdfInfo;
		for (int d = 0; d < 3; ++d) sdfInfo.hubCenter[d] = hubCenter[d];
		sdfInfo.maxRadius = 0.0;
		sdfInfo.coreRadius = maxPortDim * 0.35;

		for (int k = 0; k < N; ++k) {
			SDFTube tube;
			for (int d = 0; d < 3; ++d) {
				tube.center[d] = hubCenter[d];
				tube.axis[d] = tubes[k].dirOut[d];
			}
			NormalizeSafe(tube.axis, tubes[k].dirOut);
			tube.radius = std::min(tubes[k].w, tubes[k].h) * 0.5;
			sdfInfo.tubes.push_back(tube);
			sdfInfo.maxRadius = std::max(sdfInfo.maxRadius,
				tubes[k].distToHub + std::max(tubes[k].w, tubes[k].h));
		}

		m_sdfHubs.push_back(sdfInfo);
	}

	if (options.rebuildNow) {
		RebuildGlobalMesh();
	}

	if (false) qDebug() << "[GenerateNBranchFromRings] generated reusable 2D fan topology, N =" << N;
	return true;
}

void CageInteractorStyle::GenerateNBranch3DSkeleton()
{
	m_use3DHubSdfProjection = false;

	const int N = static_cast<int>(m_selectedFacesForBranch.size());
	if (N < 4) {
		if (false) qDebug() << "[3D branch] N < 4, fallback to GenerateNBranch().";
		GenerateNBranch();
		return;
	}

	if (m_skelNodes.empty() || m_skelEdges.empty()) {
		if (false) qDebug() << "[3D branch] skeleton is empty, fallback to GenerateNBranch().";
		GenerateNBranch();
		return;
	}

	auto CopyRing = [](double dst[4][3], const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	};

	auto Distance = [](const double a[3], const double b[3]) -> double
	{
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto NormalizeSafe = [](double v[3], const double fallback[3]) -> bool
	{
		if (vtkMath::Norm(v) < 1e-8) {
			v[0] = fallback[0];
			v[1] = fallback[1];
			v[2] = fallback[2];

			if (vtkMath::Norm(v) < 1e-8) {
				v[0] = 1.0;
				v[1] = 0.0;
				v[2] = 0.0;
			}

			vtkMath::Normalize(v);
			return false;
		}

		vtkMath::Normalize(v);
		return true;
	};

	auto ComputeRingCenter = [](const double ring[4][3], double c[3])
	{
		c[0] = c[1] = c[2] = 0.0;

		for (int i = 0; i < 4; ++i) {
			c[0] += ring[i][0] * 0.25;
			c[1] += ring[i][1] * 0.25;
			c[2] += ring[i][2] * 0.25;
		}
	};

	auto TranslateRing = [](const double src[4][3], const double dir[3], double dist, double dst[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d] + dir[d] * dist;
			}
		}
	};

	auto DotEdgeNormalized = [](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double
	{
		double va[3] = {
			a1[0] - a0[0],
			a1[1] - a0[1],
			a1[2] - a0[2]
		};

		double vb[3] = {
			b1[0] - b0[0],
			b1[1] - b0[1],
			b1[2] - b0[2]
		};

		double na = vtkMath::Norm(va);
		double nb = vtkMath::Norm(vb);

		if (na < 1e-10 || nb < 1e-10) {
			return 0.0;
		}

		va[0] /= na;
		va[1] /= na;
		va[2] /= na;

		vb[0] /= nb;
		vb[1] /= nb;
		vb[2] /= nb;

		return vtkMath::Dot(va, vb);
	};

	struct BridgeScore
	{
		double minEdgeDot;
		double avgEdgeDot;
		double totalLength;
		double maxDiagDiff;
	};

	auto EvaluateBridge = [&DotEdgeNormalized, &Distance](
		const double startRing[4][3],
		const double endRing[4][3]) -> BridgeScore
	{
		BridgeScore score;
		score.minEdgeDot = 1.0;
		score.avgEdgeDot = 0.0;
		score.totalLength = 0.0;
		score.maxDiagDiff = 0.0;

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double dot = DotEdgeNormalized(
				startRing[i],
				startRing[j],
				endRing[i],
				endRing[j]
			);

			score.minEdgeDot = std::min(score.minEdgeDot, dot);
			score.avgEdgeDot += dot * 0.25;
			score.totalLength += Distance(startRing[i], endRing[i]);

			double diagA = Distance(startRing[i], endRing[j]);
			double diagB = Distance(startRing[j], endRing[i]);

			score.maxDiagDiff = std::max(score.maxDiagDiff, std::abs(diagA - diagB));
		}

		return score;
	};

	auto IsBetterScore = [](const BridgeScore& a, const BridgeScore& b) -> bool
	{
		if (a.minEdgeDot > b.minEdgeDot + 1e-6) return true;
		if (a.minEdgeDot < b.minEdgeDot - 1e-6) return false;

		if (a.avgEdgeDot > b.avgEdgeDot + 1e-6) return true;
		if (a.avgEdgeDot < b.avgEdgeDot - 1e-6) return false;

		if (a.maxDiagDiff < b.maxDiagDiff - 1e-6) return true;
		if (a.maxDiagDiff > b.maxDiagDiff + 1e-6) return false;

		return a.totalLength < b.totalLength;
	};

	auto BuildTubeCandidateRing = [](
		const double tubeCanonical[4][3],
		bool reverse,
		int shift,
		double outRing[4][3])
	{
		static const int forwardBase[4] = { 0, 1, 2, 3 };
		static const int reverseBase[4] = { 0, 3, 2, 1 };

		const int* base = reverse ? reverseBase : forwardBase;

		for (int i = 0; i < 4; ++i) {
			int id = base[(i + shift) % 4];

			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = tubeCanonical[id][d];
			}
		}
	};

	auto AddSegmentBetweenRings = [&](const double a[4][3], const double b[4][3]) -> int
	{
		// 统一的安全 ring-to-ring 连接。
		// 不能直接相信 b[0..3] 与 a[0..3] 天然对应，否则在长段、回环段、
		// 自动/手动重连段中会偶发产生麻花面；麻花面视觉上可能仍像贴合，
		// 但导出时会导致 sewing 失败并留下 free wire。
		double bestB[4][3];
		CopyRing(bestB, b);

		BridgeScore rawScore = EvaluateBridge(a, b);
		BridgeScore bestScore = rawScore;

		for (int family = 0; family < 2; ++family) {
			const bool reverse = (family == 1);
			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				BuildTubeCandidateRing(b, reverse, shift, candidate);
				BridgeScore score = EvaluateBridge(a, candidate);

				if (IsBetterScore(score, bestScore)) {
					bestScore = score;
					CopyRing(bestB, candidate);
				}
			}
		}

		const bool changed =
			(bestScore.minEdgeDot > rawScore.minEdgeDot + 1e-6) ||
			(bestScore.avgEdgeDot > rawScore.avgEdgeDot + 1e-6) ||
			(bestScore.maxDiagDiff < rawScore.maxDiagDiff - 1e-6) ||
			(bestScore.totalLength < rawScore.totalLength - 1e-6);

		if (changed || bestScore.minEdgeDot < 0.05) {
			if (false) qDebug() << "[RingBridgeSafe]"
				<< "index =" << static_cast<int>(m_allCages.size())
				<< "changed =" << changed
				<< "rawMinDot =" << rawScore.minEdgeDot
				<< "bestMinDot =" << bestScore.minEdgeDot
				<< "rawAvgDot =" << rawScore.avgEdgeDot
				<< "bestAvgDot =" << bestScore.avgEdgeDot
				<< "rawDiagDiff =" << rawScore.maxDiagDiff
				<< "bestDiagDiff =" << bestScore.maxDiagDiff
				<< "rawLen =" << rawScore.totalLength
				<< "bestLen =" << bestScore.totalLength;
		}

		CageSegment segment;
		CopyRing(segment.pStart, a);
		CopyRing(segment.pEnd, bestB);

		int index = static_cast<int>(m_allCages.size());
		m_allCages.push_back(segment);
		return index;
	};

	auto MeasureRingExtent = [](const double ring[4][3], const double axis[3]) -> double
	{
		double minT = std::numeric_limits<double>::max();
		double maxT = -std::numeric_limits<double>::max();

		for (int i = 0; i < 4; ++i) {
			double t = vtkMath::Dot(ring[i], axis);
			minT = std::min(minT, t);
			maxT = std::max(maxT, t);
		}

		return std::max(0.0, maxT - minT);
	};

	struct BranchPort3D
	{
		SelectedFace sel;
		int anchorNode = -1;

		double ring[4][3];

		double center[3];
		double normal[3];

		double u[3];
		double v[3];

		double dirFromHub[3];

		double width = 0.0;
		double height = 0.0;
		double distToHub = 0.0;
	};

	std::vector<BranchPort3D> ports;
	ports.reserve(N);

	for (const SelectedFace& sel : m_selectedFacesForBranch) {
		if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
			if (false) qDebug() << "[3D branch] invalid selected cage index =" << sel.cageIndex;
			return;
		}

		BranchPort3D port;
		port.sel = sel;
		port.anchorNode = FindNearestSkeletonNodeToSelectedFace(sel, port.center, port.normal);

		if (port.anchorNode < 0) {
			if (false) qDebug() << "[3D branch] selected face cannot bind to skeleton.";
			return;
		}

		const CageSegment& cage = m_allCages[sel.cageIndex];

		for (int i = 0; i < 4; ++i) {
			const double* p = sel.isStart ? cage.pStart[i] : cage.pEnd[i];

			for (int d = 0; d < 3; ++d) {
				port.ring[i][d] = p[d];
			}
		}

		double rawU[3] = {
			port.ring[1][0] - port.ring[0][0],
			port.ring[1][1] - port.ring[0][1],
			port.ring[1][2] - port.ring[0][2]
		};

		double rawV[3] = {
			port.ring[3][0] - port.ring[0][0],
			port.ring[3][1] - port.ring[0][1],
			port.ring[3][2] - port.ring[0][2]
		};

		port.width = vtkMath::Norm(rawU);
		port.height = vtkMath::Norm(rawV);

		if (port.width < 1e-8 || port.height < 1e-8) {
			if (false) qDebug() << "[3D branch] selected face ring is degenerated.";
			return;
		}

		vtkMath::Normalize(rawU);
		vtkMath::Normalize(rawV);

		for (int d = 0; d < 3; ++d) {
			port.u[d] = rawU[d];
			port.v[d] = rawV[d];
		}

		ports.push_back(port);
	}

	std::vector<std::vector<int>> adjList = BuildSkeletonAdjacency();
	const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);

	if (adjList.empty() || nodeCount <= 0) {
		if (false) qDebug() << "[3D branch] skeleton graph is empty.";
		GenerateNBranch();
		return;
	}

	auto GetNodePoint = [&](int node, double p[3])
	{
		p[0] = m_skelNodes[3 * node];
		p[1] = m_skelNodes[3 * node + 1];
		p[2] = m_skelNodes[3 * node + 2];
	};

	auto ComputeGraphDistances = [&](int startNode)
	{
		std::vector<int> dist(nodeCount, -1);

		if (startNode < 0 || startNode >= nodeCount) {
			return dist;
		}

		std::queue<int> q;
		dist[startNode] = 0;
		q.push(startNode);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (curr < 0 || curr >= static_cast<int>(adjList.size())) continue;

			for (int next : adjList[curr]) {
				if (next < 0 || next >= nodeCount) continue;
				if (dist[next] >= 0) continue;

				dist[next] = dist[curr] + 1;
				q.push(next);
			}
		}

		return dist;
	};

	std::vector<std::vector<int>> allDistances;
	allDistances.reserve(ports.size());

	for (const BranchPort3D& port : ports) {
		allDistances.push_back(ComputeGraphDistances(port.anchorNode));
	}

	auto ScoreSkeletonNode = [&](int node, bool requireJunction, bool requireFront) -> double
	{
		if (node < 0 || node >= nodeCount) return 1e100;

		if (requireJunction) {
			if (node >= static_cast<int>(adjList.size()) || adjList[node].size() < 3) {
				return 1e100;
			}
		}

		double p[3];
		GetNodePoint(node, p);

		double score = 0.0;

		for (size_t i = 0; i < ports.size(); ++i) {
			int graphDist = allDistances[i][node];

			if (graphDist < 0) {
				return 1e100;
			}

			double toNode[3] = {
				p[0] - ports[i].center[0],
				p[1] - ports[i].center[1],
				p[2] - ports[i].center[2]
			};

			double len = vtkMath::Norm(toNode);
			double frontDot = 0.0;

			if (len > 1e-8) {
				toNode[0] /= len;
				toNode[1] /= len;
				toNode[2] /= len;
				frontDot = vtkMath::Dot(toNode, ports[i].normal);
			}

			if (requireFront && frontDot <= -0.05) {
				return 1e100;
			}

			score += static_cast<double>(graphDist) * 1000.0;
			score += len;
			score += std::max(0.0, 0.20 - frontDot) * 250.0;
		}

		if (node < static_cast<int>(adjList.size())) {
			score -= static_cast<double>(adjList[node].size()) * 10.0;
		}

		return score;
	};

	int hubNode = -1;
	double bestNodeScore = 1e100;

	for (int pass = 0; pass < 4; ++pass) {
		bool requireJunction = pass < 2;
		bool requireFront = pass == 0 || pass == 2;

		bestNodeScore = 1e100;
		hubNode = -1;

		for (int node = 0; node < nodeCount; ++node) {
			double score = ScoreSkeletonNode(node, requireJunction, requireFront);

			if (score < bestNodeScore) {
				bestNodeScore = score;
				hubNode = node;
			}
		}

		if (hubNode >= 0) {
			break;
		}
	}

	if (hubNode < 0) {
		if (false) qDebug() << "[3D branch] cannot find hub node.";
		GenerateNBranch();
		return;
	}

	double hubCenter[3];
	GetNodePoint(hubNode, hubCenter);

	double maxPortDim = 0.0;

	for (BranchPort3D& port : ports) {
		double dir[3] = {
			port.center[0] - hubCenter[0],
			port.center[1] - hubCenter[1],
			port.center[2] - hubCenter[2]
		};

		port.distToHub = vtkMath::Norm(dir);

		if (port.distToHub > 1e-8) {
			vtkMath::Normalize(dir);
		}
		else {
			dir[0] = -port.normal[0];
			dir[1] = -port.normal[1];
			dir[2] = -port.normal[2];
			NormalizeSafe(dir, port.u);
			port.distToHub = std::max(port.width, port.height) * 3.0;
		}

		for (int d = 0; d < 3; ++d) {
			port.dirFromHub[d] = dir[d];
		}

		maxPortDim = std::max(maxPortDim, std::max(port.width, port.height));
	}

	// ============================================================
	// 1. 三维分支统一分组：所有进入 3D 的端面都先分组。
	//    当前模板只处理已经遇到的两类：
	//    - 主干 1+1：直接生成主干，再把其它 1/2 组作为侧枝插入；
	//    - 主干 1+2：先用二维中心 hub 生成主扇区，再复用侧枝插入流程。
	// ============================================================

	struct PortGroup3D
	{
		std::vector<int> members;
		double dir[3] = { 0.0, 0.0, 0.0 };
		double center[3] = { 0.0, 0.0, 0.0 };
	};

	std::vector<PortGroup3D> portGroups;
	const bool enablePortGrouping = true;

	auto FinalizePortGroup = [&](PortGroup3D& group)
	{
		if (group.members.empty()) {
			return;
		}

		group.center[0] = group.center[1] = group.center[2] = 0.0;
		group.dir[0] = group.dir[1] = group.dir[2] = 0.0;

		for (int p : group.members) {
			for (int d = 0; d < 3; ++d) {
				group.center[d] += ports[p].center[d];
				group.dir[d] += ports[p].dirFromHub[d];
			}
		}

		double invCount = 1.0 / std::max(1, static_cast<int>(group.members.size()));
		for (int d = 0; d < 3; ++d) {
			group.center[d] *= invCount;
		}

		if (vtkMath::Norm(group.dir) < 1e-6) {
			for (int d = 0; d < 3; ++d) {
				group.dir[d] = group.center[d] - hubCenter[d];
			}
		}

		NormalizeSafe(group.dir, ports[group.members.front()].dirFromHub);
	};

	auto PairProjectedDirDot = [&](int a, int b) -> double
	{
		double sep[3] = {
			ports[b].center[0] - ports[a].center[0],
			ports[b].center[1] - ports[a].center[1],
			ports[b].center[2] - ports[a].center[2]
		};

		if (vtkMath::Norm(sep) > 1e-8) {
			vtkMath::Normalize(sep);
		}
		else {
			sep[0] = sep[1] = sep[2] = 0.0;
		}

		double da[3] = { ports[a].dirFromHub[0], ports[a].dirFromHub[1], ports[a].dirFromHub[2] };
		double db[3] = { ports[b].dirFromHub[0], ports[b].dirFromHub[1], ports[b].dirFromHub[2] };

		double ca = vtkMath::Dot(da, sep);
		double cb = vtkMath::Dot(db, sep);

		for (int d = 0; d < 3; ++d) {
			da[d] -= ca * sep[d];
			db[d] -= cb * sep[d];
		}

		if (vtkMath::Norm(da) < 1e-8 || vtkMath::Norm(db) < 1e-8) {
			return vtkMath::Dot(ports[a].dirFromHub, ports[b].dirFromHub);
		}

		vtkMath::Normalize(da);
		vtkMath::Normalize(db);
		return vtkMath::Dot(da, db);
	};

	auto NormalDotForGrouping = [&](int a, int b) -> double
	{
		double na[3] = { ports[a].normal[0], ports[a].normal[1], ports[a].normal[2] };
		double nb[3] = { ports[b].normal[0], ports[b].normal[1], ports[b].normal[2] };

		if (vtkMath::Norm(na) < 1e-8 || vtkMath::Norm(nb) < 1e-8) {
			return -1.0;
		}

		vtkMath::Normalize(na);
		vtkMath::Normalize(nb);
		return vtkMath::Dot(na, nb);
	};

	// 分组的大前提：两个端面法向必须近似同向。
	// 反向或接近垂直的端面即使距离很近，也不能组成一个 pair 组。
	const double pairNormalDotLimit = 0.50;

	auto PairNormalCompatible = [&](int a, int b) -> bool
	{
		return NormalDotForGrouping(a, b) >= std::max(pairNormalDotLimit, 0.82);
	};

	auto FirstStepFromHubToPort = [&](int portIndex) -> int
	{
		if (portIndex < 0 ||
			portIndex >= static_cast<int>(ports.size()) ||
			hubNode < 0 ||
			hubNode >= static_cast<int>(adjList.size()) ||
			portIndex >= static_cast<int>(allDistances.size()) ||
			hubNode >= static_cast<int>(allDistances[portIndex].size()))
		{
			return -1;
		}

		const int distHub = allDistances[portIndex][hubNode];
		if (distHub < 0) {
			return -1;
		}

		if (distHub == 0) {
			return ports[portIndex].anchorNode;
		}

		for (int neighbor : adjList[hubNode]) {
			if (neighbor < 0 ||
				neighbor >= static_cast<int>(allDistances[portIndex].size()))
			{
				continue;
			}

			if (allDistances[portIndex][neighbor] == distHub - 1) {
				return neighbor;
			}
		}

		return -1;
	};

	auto PairHardCompatible = [&](int a, int b) -> bool
	{
		if (!PairNormalCompatible(a, b)) {
			return false;
		}

		const double dirDot =
			vtkMath::Dot(ports[a].dirFromHub, ports[b].dirFromHub);
		if (dirDot < 0.30) {
			return false;
		}

		const double centerDist =
			std::sqrt(vtkMath::Distance2BetweenPoints(ports[a].center, ports[b].center));
		const double distanceLimit =
			std::max(
				maxPortDim * 2.25,
				0.55 * (
					std::max(ports[a].width, ports[a].height) +
					std::max(ports[b].width, ports[b].height)));

		if (centerDist > distanceLimit) {
			return false;
		}

		double normalMean[3] = {
			ports[a].normal[0] + ports[b].normal[0],
			ports[a].normal[1] + ports[b].normal[1],
			ports[a].normal[2] + ports[b].normal[2]
		};

		if (vtkMath::Norm(normalMean) > 1e-8) {
			vtkMath::Normalize(normalMean);

			double delta[3] = {
				ports[b].center[0] - ports[a].center[0],
				ports[b].center[1] - ports[a].center[1],
				ports[b].center[2] - ports[a].center[2]
			};

			const double axialSep = std::abs(vtkMath::Dot(delta, normalMean));
			const double axialLimit =
				std::max(maxPortDim * 0.75, centerDist * 0.45);
			if (axialSep > axialLimit) {
				return false;
			}
		}

		return true;
	};

	auto PairScore = [&](int a, int b) -> double
	{
		if (!PairHardCompatible(a, b)) {
			return 1e100;
		}

		double dist = std::sqrt(vtkMath::Distance2BetweenPoints(ports[a].center, ports[b].center));
		double dotProjected = PairProjectedDirDot(a, b);
		return dist - dotProjected * maxPortDim * 3.0;
	};

	auto PairSoftScore = [&](int a, int b) -> double
	{
		const double centerDist =
			std::sqrt(vtkMath::Distance2BetweenPoints(ports[a].center, ports[b].center));
		const double dirDot =
			vtkMath::Dot(ports[a].dirFromHub, ports[b].dirFromHub);
		const double normalDot = NormalDotForGrouping(a, b);
		const double projectedDot = PairProjectedDirDot(a, b);

		double normalMean[3] = {
			ports[a].normal[0] + ports[b].normal[0],
			ports[a].normal[1] + ports[b].normal[1],
			ports[a].normal[2] + ports[b].normal[2]
		};

		double axialSep = 0.0;
		if (vtkMath::Norm(normalMean) > 1e-8) {
			vtkMath::Normalize(normalMean);
			double delta[3] = {
				ports[b].center[0] - ports[a].center[0],
				ports[b].center[1] - ports[a].center[1],
				ports[b].center[2] - ports[a].center[2]
			};
			axialSep = std::abs(vtkMath::Dot(delta, normalMean));
		}

		double score = centerDist;
		score -= projectedDot * maxPortDim * 2.15;
		score -= std::max(0.0, normalDot) * maxPortDim * 0.80;
		score += std::max(0.0, 0.20 - dirDot) * maxPortDim * 2.35;
		score += axialSep * 0.70;

		int stepA = FirstStepFromHubToPort(a);
		int stepB = FirstStepFromHubToPort(b);
		if (stepA >= 0 && stepB >= 0) {
			score += (stepA == stepB) ? -maxPortDim * 1.25 : maxPortDim * 0.95;
		}

		return score;
	};

	auto PairSoftCompatible = [&](int a, int b) -> bool
	{
		const double normalDot = NormalDotForGrouping(a, b);
		if (normalDot < 0.35) {
			return false;
		}

		const double dirDot =
			vtkMath::Dot(ports[a].dirFromHub, ports[b].dirFromHub);
		if (dirDot < -0.15) {
			return false;
		}

		const double centerDist =
			std::sqrt(vtkMath::Distance2BetweenPoints(ports[a].center, ports[b].center));
		const double distanceLimit =
			std::max(
				maxPortDim * 3.75,
				0.95 * (
					std::max(ports[a].width, ports[a].height) +
					std::max(ports[b].width, ports[b].height)));

		if (centerDist > distanceLimit) {
			return false;
		}

		double normalMean[3] = {
			ports[a].normal[0] + ports[b].normal[0],
			ports[a].normal[1] + ports[b].normal[1],
			ports[a].normal[2] + ports[b].normal[2]
		};

		if (vtkMath::Norm(normalMean) > 1e-8) {
			vtkMath::Normalize(normalMean);

			double delta[3] = {
				ports[b].center[0] - ports[a].center[0],
				ports[b].center[1] - ports[a].center[1],
				ports[b].center[2] - ports[a].center[2]
			};

			const double axialSep = std::abs(vtkMath::Dot(delta, normalMean));
			const double axialLimit =
				std::max(maxPortDim * 1.35, centerDist * 0.78);
			if (axialSep > axialLimit) {
				return false;
			}
		}

		return true;
	};

	std::vector<std::tuple<double, int, int>> pairCandidates;
	pairCandidates.reserve(ports.size() * ports.size());
	for (int i = 0; i < static_cast<int>(ports.size()); ++i) {
		for (int j = i + 1; j < static_cast<int>(ports.size()); ++j) {
			if (!PairHardCompatible(i, j)) {
				continue;
			}

			pairCandidates.push_back(std::make_tuple(PairScore(i, j), i, j));
		}
	}

	std::sort(pairCandidates.begin(), pairCandidates.end(),
		[](const std::tuple<double, int, int>& a, const std::tuple<double, int, int>& b)
	{
		return std::get<0>(a) < std::get<0>(b);
	});

	std::vector<bool> grouped(ports.size(), false);
	const double pairScoreLimit = std::max(maxPortDim * 2.40, 1e-6);
	std::vector<std::pair<int, int>> selectedPairs;

	{
		const int portCount = static_cast<int>(ports.size());
		const int maxPairCount = portCount / 2;
		int bestPairCount = -1;
		double bestTotalScore = std::numeric_limits<double>::max();
		std::vector<bool> used(portCount, false);
		std::vector<std::pair<int, int>> currentPairs;
		std::vector<std::pair<int, int>> bestPairs;

		std::function<void(size_t, double)> SearchPairs =
			[&](size_t startIndex, double totalScore)
		{
			const int currentPairCount = static_cast<int>(currentPairs.size());
			int unusedCount = 0;
			for (bool isUsed : used) {
				if (!isUsed) {
					++unusedCount;
				}
			}

			if (currentPairCount + unusedCount / 2 < bestPairCount) {
				return;
			}

			if (currentPairCount > bestPairCount ||
				(currentPairCount == bestPairCount && totalScore < bestTotalScore))
			{
				bestPairCount = currentPairCount;
				bestTotalScore = totalScore;
				bestPairs = currentPairs;
			}

			if (currentPairCount >= maxPairCount) {
				return;
			}

			for (size_t idx = startIndex; idx < pairCandidates.size(); ++idx) {
				double score = std::get<0>(pairCandidates[idx]);
				if (score > pairScoreLimit) {
					continue;
				}

				int a = std::get<1>(pairCandidates[idx]);
				int b = std::get<2>(pairCandidates[idx]);
				if (used[a] || used[b]) {
					continue;
				}

				used[a] = true;
				used[b] = true;
				currentPairs.push_back({ a, b });
				SearchPairs(idx + 1, totalScore + score);
				currentPairs.pop_back();
				used[a] = false;
				used[b] = false;
			}
		};

		SearchPairs(0, 0.0);
		selectedPairs = bestPairs;
	}

	if (N == 5 && selectedPairs.size() < 2) {
		std::vector<std::tuple<double, int, int>> softPairCandidates;
		softPairCandidates.reserve(ports.size() * ports.size());
		for (int i = 0; i < static_cast<int>(ports.size()); ++i) {
			for (int j = i + 1; j < static_cast<int>(ports.size()); ++j) {
				if (!PairSoftCompatible(i, j)) {
					continue;
				}

				softPairCandidates.push_back(std::make_tuple(PairSoftScore(i, j), i, j));
			}
		}

		std::sort(softPairCandidates.begin(), softPairCandidates.end(),
			[](const std::tuple<double, int, int>& a, const std::tuple<double, int, int>& b)
		{
			return std::get<0>(a) < std::get<0>(b);
		});

		double bestTwoPairScore = std::numeric_limits<double>::max();
		std::vector<std::pair<int, int>> bestTwoPairs;
		const double softPairScoreLimit = std::max(maxPortDim * 6.0, 1e-6);
		for (size_t ia = 0; ia < softPairCandidates.size(); ++ia) {
			double scoreA = std::get<0>(softPairCandidates[ia]);
			if (scoreA > softPairScoreLimit) {
				continue;
			}

			int a0 = std::get<1>(softPairCandidates[ia]);
			int a1 = std::get<2>(softPairCandidates[ia]);
			for (size_t ib = ia + 1; ib < softPairCandidates.size(); ++ib) {
				double scoreB = std::get<0>(softPairCandidates[ib]);
				if (scoreB > softPairScoreLimit) {
					continue;
				}

				int b0 = std::get<1>(softPairCandidates[ib]);
				int b1 = std::get<2>(softPairCandidates[ib]);
				if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) {
					continue;
				}

				double totalScore = scoreA + scoreB;
				if (totalScore < bestTwoPairScore) {
					bestTwoPairScore = totalScore;
					bestTwoPairs.clear();
					bestTwoPairs.push_back({ a0, a1 });
					bestTwoPairs.push_back({ b0, b1 });
				}
			}
		}

		if (bestTwoPairs.size() == 2) {
			selectedPairs = bestTwoPairs;
			qDebug() << "[3D branch][PairRecover] N=5 recovered two pair groups";
		}
	}

	for (const auto& pair : selectedPairs) {
		int a = pair.first;
		int b = pair.second;
		PortGroup3D group;
		group.members.push_back(a);
		group.members.push_back(b);
		FinalizePortGroup(group);
		portGroups.push_back(group);
		grouped[a] = true;
		grouped[b] = true;
	}

	for (int p = 0; p < static_cast<int>(ports.size()); ++p) {
		if (grouped[p]) {
			continue;
		}

		PortGroup3D group;
		group.members.push_back(p);
		FinalizePortGroup(group);
		portGroups.push_back(group);
	}

	QString groupText;
	for (int g = 0; g < static_cast<int>(portGroups.size()); ++g) {
		if (!groupText.isEmpty()) {
			groupText += " + ";
		}
		groupText += QString::number(static_cast<int>(portGroups[g].members.size()));
	}
	qDebug() << "[3D branch][GroupSummary] N =" << N << ", groups =" << groupText;

	int mainA = 0;
	int mainB = 1;
	int mainGroupA = -1;
	int mainGroupB = -1;
	bool useVirtualMainEnd = false;
	std::vector<int> mainEndGroupPorts;
	std::vector<std::pair<int, int>> symmetricSidePairs;

	// 优先识别 1+1 主干：两个 singleton 组若方向最相反，则直接走普通四输入式主干插入模板。
	// 如果没有 1+1 主干，再允许 1+2 主扇区模板。
	if (enablePortGrouping && portGroups.size() >= 2) {
		double bestSingletonScore = -1e100;
		int bestSingleA = -1;
		int bestSingleB = -1;

		for (int ga = 0; ga < static_cast<int>(portGroups.size()); ++ga) {
			if (portGroups[ga].members.size() != 1) continue;

			for (int gb = ga + 1; gb < static_cast<int>(portGroups.size()); ++gb) {
				if (portGroups[gb].members.size() != 1) continue;

				double dot = vtkMath::Dot(portGroups[ga].dir, portGroups[gb].dir);
				double centerDist = std::sqrt(vtkMath::Distance2BetweenPoints(
					portGroups[ga].center,
					portGroups[gb].center
				));
				double score = (-dot) * 100.0 + centerDist * 0.02;

				if (score > bestSingletonScore) {
					bestSingletonScore = score;
					bestSingleA = ga;
					bestSingleB = gb;
				}
			}
		}

		if (bestSingleA >= 0 && bestSingleB >= 0) {
			mainGroupA = bestSingleA;
			mainGroupB = bestSingleB;
		}
		else {
			double bestMainGroupScore = -1e100;

			for (int ga = 0; ga < static_cast<int>(portGroups.size()); ++ga) {
				for (int gb = ga + 1; gb < static_cast<int>(portGroups.size()); ++gb) {
					int sa = static_cast<int>(portGroups[ga].members.size());
					int sb = static_cast<int>(portGroups[gb].members.size());

					// 当前只让 1+2 作为中心 hub 主扇区；其它组合先不强行套模板。
					if (!((sa == 1 && sb == 2) || (sa == 2 && sb == 1))) {
						continue;
					}

					double dot = vtkMath::Dot(portGroups[ga].dir, portGroups[gb].dir);
					double centerDist = std::sqrt(vtkMath::Distance2BetweenPoints(
						portGroups[ga].center,
						portGroups[gb].center
					));
					double score = (-dot) * 100.0 + centerDist * 0.02;

					if (score > bestMainGroupScore) {
						bestMainGroupScore = score;
						mainGroupA = ga;
						mainGroupB = gb;
					}
				}
			}
		}

		if (mainGroupA >= 0 && mainGroupB >= 0) {
			// 对 1+2 模板，单端口作为 A，多端口作为 B 端主扇区。
			if (portGroups[mainGroupA].members.size() > portGroups[mainGroupB].members.size()) {
				std::swap(mainGroupA, mainGroupB);
			}

			mainA = portGroups[mainGroupA].members.front();
			mainB = portGroups[mainGroupB].members.front();

			if (portGroups[mainGroupB].members.size() == 2) {
				useVirtualMainEnd = true;
				mainEndGroupPorts = portGroups[mainGroupB].members;
			}

			for (int g = 0; g < static_cast<int>(portGroups.size()); ++g) {
				if (g == mainGroupA || g == mainGroupB) {
					continue;
				}

				if (portGroups[g].members.size() == 2) {
					symmetricSidePairs.push_back({ portGroups[g].members[0], portGroups[g].members[1] });
				}
			}

			qDebug() << "[3D branch][GroupSummary] main groups ="
				<< static_cast<int>(portGroups[mainGroupA].members.size())
				<< "+"
				<< static_cast<int>(portGroups[mainGroupB].members.size())
				<< ", main ports =" << mainA << mainB
				<< ", side pairs =" << static_cast<int>(symmetricSidePairs.size());
		}
	}

	if (!enablePortGrouping || mainGroupA < 0 || mainGroupB < 0 || mainA == mainB) {
		double minDot = 1.0;

		for (int i = 0; i < static_cast<int>(ports.size()); ++i) {
			for (int j = i + 1; j < static_cast<int>(ports.size()); ++j) {
				double dot = vtkMath::Dot(ports[i].dirFromHub, ports[j].dirFromHub);

				if (dot < minDot) {
					minDot = dot;
					mainA = i;
					mainB = j;
				}
			}
		}
	}

	auto IsMainEndGroupPort = [&](int p) -> bool
	{
		if (!useVirtualMainEnd) {
			return false;
		}

		return std::find(mainEndGroupPorts.begin(), mainEndGroupPorts.end(), p) != mainEndGroupPorts.end();
	};

	auto FindSymmetricPairMate = [&](int p) -> int
	{
		for (const auto& pr : symmetricSidePairs) {
			if (pr.first == p) return pr.second;
			if (pr.second == p) return pr.first;
		}

		return -1;
	};

	// ============================================================
	// 1.4 1+2 主扇区模板开关
	// ============================================================
	// 分组只用于决策；真正生成时：
	// - mainEndGroupPorts 作为主干端部分流，不能再作为侧孔插入；
	// - symmetricSidePairs 作为成对侧孔插入，不能逐个独立投影。
	bool useGrouped122Template =
		enablePortGrouping &&
		useVirtualMainEnd &&
		mainEndGroupPorts.size() == 2 &&
		symmetricSidePairs.size() == 1;

	// 1+2+2 主骨架生成后，函数会把“选中的 host bridge 两端 + side0/side1”
	// 重映射成普通四输入流程。这个标记只用于继续启用成对侧枝的同面插入逻辑，
	// 不能再触发旧的 skip side-pair 逻辑。
	bool grouped122SideInsertMode = false;

	if (useGrouped122Template) {
		if (false) qDebug() << "[3D branch][Grouped122] enabled: main fan + reusable existing side-insert flow"
			<< "mainEndPorts =" << mainEndGroupPorts[0] << mainEndGroupPorts[1]
			<< "sidePair =" << symmetricSidePairs[0].first << symmetricSidePairs[0].second;

		// ============================================================
		// 1+2+2 结构化处理：
		//   A + B0 + B1 先用二维 GenerateNBranchFromRings 生成主骨架；
		//   然后从这次新生成的三条 bridge 里，选择最接近 sidePair 的 host bridge；
		//   把这条 host bridge 从 m_allCages 中取出，用原来已经跑通的四输入插入流程
		//   继续处理 side0/side1。
		//
		// 这样做的目的：
		//   - 主骨架仍然完全复用二维 NBranch 逻辑；
		//   - 侧枝插入不再重新写一套，而是复用下面原本成功的 4 输入插入流程；
		//   - 原来的普通 4 输入路径不受影响。
		// ============================================================

		auto MakeRealNBranchPort = [&](int portIndex) -> NBranchInputPort
		{
			NBranchInputPort out;
			out.sourceCageIndex = ports[portIndex].sel.cageIndex;
			out.sourceIsStart = ports[portIndex].sel.isStart;
			out.isVirtual = false;
			CopyRing(out.ring, ports[portIndex].ring);

			for (int d = 0; d < 3; ++d) {
				out.preferredDir[d] = ports[portIndex].dirFromHub[d];
			}
			if (vtkMath::Norm(out.preferredDir) > 1e-8) {
				vtkMath::Normalize(out.preferredDir);
				out.hasPreferredDir = true;
			}
			return out;
		};

		const int side0 = symmetricSidePairs[0].first;
		const int side1 = symmetricSidePairs[0].second;
		const int b0 = mainEndGroupPorts[0];
		const int b1 = mainEndGroupPorts[1];

		double aCenter[3];
		double b0Center[3];
		double b1Center[3];
		double s0Center[3];
		double s1Center[3];
		ComputeRingCenter(ports[mainA].ring, aCenter);
		ComputeRingCenter(ports[b0].ring, b0Center);
		ComputeRingCenter(ports[b1].ring, b1Center);
		ComputeRingCenter(ports[side0].ring, s0Center);
		ComputeRingCenter(ports[side1].ring, s1Center);

		double bMid[3] = {
			0.5 * (b0Center[0] + b1Center[0]),
			0.5 * (b0Center[1] + b1Center[1]),
			0.5 * (b0Center[2] + b1Center[2])
		};

		double sidePairCenter[3] = {
			0.5 * (s0Center[0] + s1Center[0]),
			0.5 * (s0Center[1] + s1Center[1]),
			0.5 * (s0Center[2] + s1Center[2])
		};

		double axisAtoB[3] = {
			bMid[0] - aCenter[0],
			bMid[1] - aCenter[1],
			bMid[2] - aCenter[2]
		};
		NormalizeSafe(axisAtoB, ports[mainA].dirFromHub);

		double mainRadius = 0.5 * std::min(ports[mainA].width, ports[mainA].height);
		double avgBRadius = 0.25 * (
			std::min(ports[b0].width, ports[b0].height) +
			std::min(ports[b1].width, ports[b1].height)
			);
		double sideRadius = 0.25 * (
			std::min(ports[side0].width, ports[side0].height) +
			std::min(ports[side1].width, ports[side1].height)
			);

		double distAToBMid = Distance(aCenter, bMid);
		double requestHubConcede = std::max(mainRadius * 1.78, maxPortDim * 0.30);
		requestHubConcede = std::max(requestHubConcede, sideRadius * 1.18);
		if (distAToBMid > 1e-8) {
			requestHubConcede = std::min(requestHubConcede, distAToBMid * 0.58);
		}

		double bSafeGap = std::max(maxPortDim * 0.045, avgBRadius * 0.20);
		bSafeGap = std::min(bSafeGap, maxPortDim * 0.120);

		std::vector<NBranchInputPort> mainFanPorts;
		mainFanPorts.reserve(3);
		mainFanPorts.push_back(MakeRealNBranchPort(mainA));
		mainFanPorts.push_back(MakeRealNBranchPort(b0));
		mainFanPorts.push_back(MakeRealNBranchPort(b1));

		g_NBranchHubRingRelax.enabled = true;
		g_NBranchHubRingRelax.targetInputIndex = 0;
		for (int d = 0; d < 3; ++d) {
			g_NBranchHubRingRelax.axisAtoB[d] = axisAtoB[d];
			g_NBranchHubRingRelax.aCenter[d] = aCenter[d];
		}
		g_NBranchHubRingRelax.requestMove = requestHubConcede;
		g_NBranchHubRingRelax.safeGap = bSafeGap;
		g_NBranchHubRingRelax.maxPortDim = maxPortDim;
		g_NBranchHubRingRelax.boundaryPoints.clear();
		for (auto& group : g_NBranchHubRingRelax.boundaryPointGroups) {
			group.clear();
		}
		g_NBranchHubRingRelax.boundaryPoints.reserve(8);

		// v11：只记录 B0/B1 两个已选择输入端面的原始 8 个角点。
		// 后续 ApplyNBranchHubRingRelaxInPlace() 会让每个 hub 点直接和这 8 个点做主方向坐标硬比较。
		// 不做边界采样，不做编号匹配，不按上下边界猜对应关系。
		for (int ringId = 0; ringId < 4; ++ringId) {
			g_NBranchHubRingRelax.boundaryPoints.push_back({
				ports[b0].ring[ringId][0],
				ports[b0].ring[ringId][1],
				ports[b0].ring[ringId][2]
				});
			g_NBranchHubRingRelax.boundaryPoints.push_back({
				ports[b1].ring[ringId][0],
				ports[b1].ring[ringId][1],
				ports[b1].ring[ringId][2]
				});
		}

		int beforeCageCount = static_cast<int>(m_allCages.size());
		int beforeHubCount = static_cast<int>(m_hubCages.size());

		NBranchBuildOptions mainFanOptions;
		mainFanOptions.addSdf = true;
		mainFanOptions.rebuildNow = false;
		mainFanOptions.debugLog = true;
		mainFanOptions.forceHubCenter = false;
		mainFanOptions.forceNormalUp = false;
		mainFanOptions.radiusScale = 0.8;
		mainFanOptions.minRadiusScale = 0.045;

		bool mainOk = GenerateNBranchFromRings(mainFanPorts, mainFanOptions);

		g_NBranchHubRingRelax.enabled = false;
		g_NBranchHubRingRelax.boundaryPoints.clear();
		for (auto& group : g_NBranchHubRingRelax.boundaryPointGroups) {
			group.clear();
		}

		if (!mainOk) {
			if (false) qDebug() << "[3D branch][Grouped122] main fan failed, fallback GenerateNBranch().";
			GenerateNBranch();
			return;
		}

		int afterCageCount = static_cast<int>(m_allCages.size());
		std::vector<int> createdHostCandidates;
		for (int ci = beforeCageCount; ci < afterCageCount; ++ci) {
			createdHostCandidates.push_back(ci);
		}

		if (createdHostCandidates.empty()) {
			if (false) qDebug() << "[3D branch][Grouped122] no generated host bridge, main fan only.";
			RebuildGlobalMesh();
			return;
		}

		auto SegmentCenter = [&](const CageSegment& seg, double c[3])
		{
			c[0] = c[1] = c[2] = 0.0;
			for (int i = 0; i < 4; ++i) {
				for (int d = 0; d < 3; ++d) {
					c[d] += seg.pStart[i][d] / 8.0;
					c[d] += seg.pEnd[i][d] / 8.0;
				}
			}
		};

		auto SegmentAxis = [&](const CageSegment& seg, double axis[3])
		{
			double cs[3], ce[3];
			ComputeRingCenter(seg.pStart, cs);
			ComputeRingCenter(seg.pEnd, ce);
			for (int d = 0; d < 3; ++d) axis[d] = ce[d] - cs[d];
			NormalizeSafe(axis, axisAtoB);
		};

		// 选择最接近 sidePair 的 fan bridge 作为 host。这样后续复用普通 4 输入插入流程时，
		// 侧枝会落到最合理的主骨架表面，而不是固定插到某一条 branch。
		int chosenHostIndex = createdHostCandidates.front();
		double bestHostScore = 1e100;
		for (int ci : createdHostCandidates) {
			const CageSegment& seg = m_allCages[ci];
			double cs[3], ce[3], axis[3];
			ComputeRingCenter(seg.pStart, cs);
			ComputeRingCenter(seg.pEnd, ce);
			for (int d = 0; d < 3; ++d) axis[d] = ce[d] - cs[d];
			double len = vtkMath::Norm(axis);
			if (len < 1e-8) continue;
			for (int d = 0; d < 3; ++d) axis[d] /= len;

			double toSide[3] = {
				sidePairCenter[0] - cs[0],
				sidePairCenter[1] - cs[1],
				sidePairCenter[2] - cs[2]
			};
			double t = vtkMath::Dot(toSide, axis);
			t = std::max(0.0, std::min(len, t));
			double closest[3] = {
				cs[0] + axis[0] * t,
				cs[1] + axis[1] * t,
				cs[2] + axis[2] * t
			};
			double score = Distance(sidePairCenter, closest);

			// 轻微优先选择 mainA 那条 bridge，因为侧枝通常从主入口侧插入；
			// 但如果 B 支路确实更近，距离项仍然会胜出。
			double segMid[3];
			SegmentCenter(seg, segMid);
			double distToA = Distance(segMid, aCenter);
			score += distToA * 0.03;

			if (score < bestHostScore) {
				bestHostScore = score;
				chosenHostIndex = ci;
			}
		}

		CageSegment chosenHost = m_allCages[chosenHostIndex];

		// 删除被选中的原 host bridge，下面的普通插入流程会用同样的两端 ring
		// 生成带开孔/缓冲区/侧枝桥接的新 host，避免重叠双层主干。
		m_allCages.erase(m_allCages.begin() + chosenHostIndex);
		if (chosenHostIndex < afterCageCount) {
			--afterCageCount;
		}

		// v14：1+2+2 专用中心 hub 安全规则段。
		// 被选中的 fan bridge 一端直接接二维 NBranch 生成的中心 hub。
		// 后续侧孔插入会在这个 bridge 上重建多 cut 过渡面；如果一直重建到中心 hub 边界，
		// 就会出现“插入过渡段直接吸收中心 hub 拓扑”的情况：一侧是多节点 cut，另一侧是中心 hub 的少节点/三角边界，
		// 视觉上红线通常就落在这里，导出时也容易形成不规则自由 wire。
		// 因此这里从靠近中心 hub 的一端预留一个很短的、完全规则的四边形安全段：
		//   hubRing -> safeRing 保持原始 4 节点规则连接，不开孔、不吸收 cut；
		//   safeRing -> otherRing 才交给下面普通四端面插入流程。
		// 这样不改变原有侧孔投影/吸附逻辑，只是在中心 hub 之前多留一圈干净接口。
		{
			double cStart[3];
			double cEnd[3];
			ComputeRingCenter(chosenHost.pStart, cStart);
			ComputeRingCenter(chosenHost.pEnd, cEnd);

			const double dStartToHub = Distance(cStart, hubCenter);
			const double dEndToHub = Distance(cEnd, hubCenter);
			const bool hubAtStart = dStartToHub <= dEndToHub;

			double hostLen = Distance(cStart, cEnd);
			double safeLen = std::max(maxPortDim * 0.070, sideRadius * 0.36);
			safeLen = std::min(safeLen, maxPortDim * 0.105);
			if (hostLen > 1e-8) {
				safeLen = std::min(safeLen, hostLen * 0.22);
			}
			else {
				safeLen = 0.0;
			}

			if (safeLen > 1e-8 && hostLen > safeLen * 1.8) {
				double safeRing[4][3];
				const double alpha = safeLen / hostLen;

				if (hubAtStart) {
					for (int i = 0; i < 4; ++i) {
						for (int d = 0; d < 3; ++d) {
							safeRing[i][d] = chosenHost.pStart[i][d] * (1.0 - alpha) + chosenHost.pEnd[i][d] * alpha;
						}
					}

					AddSegmentBetweenRings(chosenHost.pStart, safeRing);
					CopyRing(chosenHost.pStart, safeRing);
				}
				else {
					for (int i = 0; i < 4; ++i) {
						for (int d = 0; d < 3; ++d) {
							safeRing[i][d] = chosenHost.pEnd[i][d] * (1.0 - alpha) + chosenHost.pStart[i][d] * alpha;
						}
					}

					AddSegmentBetweenRings(safeRing, chosenHost.pEnd);
					CopyRing(chosenHost.pEnd, safeRing);
				}

				if (false) qDebug() << "[3D branch][Grouped122][HubSafeTransition] inserted clean regular segment before center hub"
					<< "hubAtStart =" << hubAtStart
					<< "hostLen =" << hostLen
					<< "safeLen =" << safeLen
					<< "alpha =" << alpha;
			}
		}

		auto MakeHostPortFromRing = [&](const double ring[4][3], const double normalHint[3]) -> BranchPort3D
		{
			BranchPort3D hp;
			hp.sel.cageIndex = -1;
			hp.sel.isStart = false;
			hp.anchorNode = hubNode;
			CopyRing(hp.ring, ring);
			ComputeRingCenter(hp.ring, hp.center);

			double rawU[3] = {
				hp.ring[1][0] - hp.ring[0][0],
				hp.ring[1][1] - hp.ring[0][1],
				hp.ring[1][2] - hp.ring[0][2]
			};
			double rawV[3] = {
				hp.ring[3][0] - hp.ring[0][0],
				hp.ring[3][1] - hp.ring[0][1],
				hp.ring[3][2] - hp.ring[0][2]
			};

			hp.width = vtkMath::Norm(rawU);
			hp.height = vtkMath::Norm(rawV);
			NormalizeSafe(rawU, ports[mainA].u);
			NormalizeSafe(rawV, ports[mainA].v);
			for (int d = 0; d < 3; ++d) {
				hp.u[d] = rawU[d];
				hp.v[d] = rawV[d];
				hp.normal[d] = normalHint[d];
				hp.dirFromHub[d] = normalHint[d];
			}
			NormalizeSafe(hp.normal, axisAtoB);
			NormalizeSafe(hp.dirFromHub, hp.normal);
			hp.distToHub = Distance(hp.center, hubCenter);
			return hp;
		};

		double hostAxis[3];
		SegmentAxis(chosenHost, hostAxis);
		double negHostAxis[3] = { -hostAxis[0], -hostAxis[1], -hostAxis[2] };

		BranchPort3D hostStartPort = MakeHostPortFromRing(chosenHost.pStart, negHostAxis);
		BranchPort3D hostEndPort = MakeHostPortFromRing(chosenHost.pEnd, hostAxis);

		std::vector<BranchPort3D> remappedPorts;
		remappedPorts.reserve(4);
		remappedPorts.push_back(hostStartPort);
		remappedPorts.push_back(hostEndPort);
		remappedPorts.push_back(ports[side0]);
		remappedPorts.push_back(ports[side1]);

		ports = remappedPorts;
		mainA = 0;
		mainB = 1;
		mainGroupA = -1;
		mainGroupB = -1;
		useVirtualMainEnd = false;
		mainEndGroupPorts.clear();
		symmetricSidePairs.clear();

		// 1+2+2 主骨架生成后，剩余两个端口被重映射为 2/3。
		// 这里保留它们的 pair 关系，只用于后续普通四输入插入流程：
		// - 强制两个侧枝优先落到同一个 host side；
		// - 后面的 SameFaceBoundarySnap 会把它们分别吸附到该 side 的两条边界。
		// 注意 useGrouped122Template 已经被置为 false，所以不会重新进入 1+2+2 专用跳过逻辑。
		symmetricSidePairs.push_back({ 2, 3 });
		grouped122SideInsertMode = true;
		// 关键：主骨架已经生成完，后续只能走“普通四输入侧枝插入复用流程”。
		// 如果这里仍保持 true，下面旧逻辑会把 sidePair 再次跳过。
		useGrouped122Template = false;

		if (false) qDebug() << "[3D branch][Grouped122] main fan created, continue with existing 4-port insertion"
			<< "chosenHostIndex =" << chosenHostIndex
			<< "createdFanBridges =" << (afterCageCount - beforeCageCount + 1)
			<< "createdHubs =" << (static_cast<int>(m_hubCages.size()) - beforeHubCount)
			<< "requestHubConcede =" << requestHubConcede
			<< "bSafeGap =" << bSafeGap
			<< "bestHostScore =" << bestHostScore;
	}

	double mainStart[4][3];
	double mainEnd[4][3];

	CopyRing(mainStart, ports[mainA].ring);

	{
		double bestEnd[4][3];

		BridgeScore bestMainScore;
		bestMainScore.minEdgeDot = -1e100;
		bestMainScore.avgEdgeDot = -1e100;
		bestMainScore.totalLength = 1e100;
		bestMainScore.maxDiagDiff = 1e100;

		for (int family = 0; family < 2; ++family) {
			bool reverse = family == 1;

			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				BuildTubeCandidateRing(ports[mainB].ring, reverse, shift, candidate);

				BridgeScore score = EvaluateBridge(mainStart, candidate);

				if (IsBetterScore(score, bestMainScore)) {
					bestMainScore = score;
					CopyRing(bestEnd, candidate);
				}
			}
		}

		CopyRing(mainEnd, bestEnd);

		// N>=5 且主端是多端口组时，把主干终点平移到该端口组中心。
		// 这样不会把右侧双出口中的某一个当成唯一主干终点，能明显改善对称性。
		if (useVirtualMainEnd && mainGroupB >= 0) {
			double endCenter[3];
			ComputeRingCenter(mainEnd, endCenter);

			double delta[3] = {
				portGroups[mainGroupB].center[0] - endCenter[0],
				portGroups[mainGroupB].center[1] - endCenter[1],
				portGroups[mainGroupB].center[2] - endCenter[2]
			};

			for (int i = 0; i < 4; ++i) {
				for (int d = 0; d < 3; ++d) {
					mainEnd[i][d] += delta[d];
				}
			}
		}
	}

	// ============================================================
		// 2. 主干三段：左 transition / 中间开孔 / 右 transition
		//    三段侧面全部由统一 quad patch 重建，避免 T-junction。
		// ============================================================

	double mainStartCenter[3], mainEndCenter[3];
	ComputeRingCenter(mainStart, mainStartCenter);
	ComputeRingCenter(mainEnd, mainEndCenter);

	double trunkAxis[3] = {
		mainEndCenter[0] - mainStartCenter[0],
		mainEndCenter[1] - mainStartCenter[1],
		mainEndCenter[2] - mainStartCenter[2]
	};
	NormalizeSafe(trunkAxis, ports[mainA].normal);

	double trunkLen = Distance(mainStartCenter, mainEndCenter);

	double startInset = std::min(ports[mainA].width, ports[mainA].height) * 0.22;
	double endInset = std::min(ports[mainB].width, ports[mainB].height) * 0.22;

	startInset = std::max(startInset, maxPortDim * 0.035);
	endInset = std::max(endInset, maxPortDim * 0.035);

	startInset = std::min(startInset, maxPortDim * 0.14);
	endInset = std::min(endInset, maxPortDim * 0.14);

	double maxTotalInset = std::max(maxPortDim * 0.18, trunkLen * 0.32);
	if (startInset + endInset > maxTotalInset && (startInset + endInset) > 1e-8) {
		double scale = maxTotalInset / (startInset + endInset);
		startInset *= scale;
		endInset *= scale;
	}

	double hostStartInner[4][3];
	double hostEndInner[4][3];

	TranslateRing(mainStart, trunkAxis, startInset, hostStartInner);
	TranslateRing(mainEnd, trunkAxis, -endInset, hostEndInner);

	const int startBufferIndex = AddSegmentBetweenRings(mainStart, hostStartInner);
	const int trunkCageIndex = AddSegmentBetweenRings(hostStartInner, hostEndInner);
	const int endBufferIndex = AddSegmentBetweenRings(hostEndInner, mainEnd);

	for (int s = 0; s < 4; ++s) {
		m_suppressedCageSideFaces.push_back({ startBufferIndex, s });
		m_suppressedCageSideFaces.push_back({ trunkCageIndex, s });
		m_suppressedCageSideFaces.push_back({ endBufferIndex, s });
	}

	struct HostSide
	{
		int sideIndex = -1;

		double S0[3];
		double S1[3];
		double E1[3];
		double E0[3];

		double uAxis[3];
		double vAxis[3];
		double normal[3];

		double lenU = 1.0;
		double lenV = 1.0;
	};

	double trunkCenter[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < 4; ++i) {
		for (int d = 0; d < 3; ++d) {
			trunkCenter[d] += hostStartInner[i][d] / 8.0;
			trunkCenter[d] += hostEndInner[i][d] / 8.0;
		}
	}

	auto BuildHostSide = [&](const double ringA[4][3],
		const double ringB[4][3],
		int s,
		const double fallbackNormal[3],
		const double centerHint[3]) -> HostSide
	{
		int n = (s + 1) % 4;

		HostSide side;
		side.sideIndex = s;

		for (int d = 0; d < 3; ++d) {
			side.S0[d] = ringA[s][d];
			side.S1[d] = ringA[n][d];
			side.E1[d] = ringB[n][d];
			side.E0[d] = ringB[s][d];

			side.uAxis[d] = side.E0[d] - side.S0[d];
			side.vAxis[d] = side.S1[d] - side.S0[d];
		}

		side.lenU = vtkMath::Norm(side.uAxis);
		side.lenV = vtkMath::Norm(side.vAxis);

		if (side.lenU > 1e-8) {
			vtkMath::Normalize(side.uAxis);
		}
		else {
			side.uAxis[0] = trunkAxis[0];
			side.uAxis[1] = trunkAxis[1];
			side.uAxis[2] = trunkAxis[2];
			side.lenU = 1.0;
		}

		if (side.lenV > 1e-8) {
			vtkMath::Normalize(side.vAxis);
		}
		else {
			side.vAxis[0] = ports[mainA].u[0];
			side.vAxis[1] = ports[mainA].u[1];
			side.vAxis[2] = ports[mainA].u[2];
			side.lenV = 1.0;
		}

		vtkMath::Cross(side.uAxis, side.vAxis, side.normal);
		NormalizeSafe(side.normal, fallbackNormal);

		double faceCenter[3] = {
			(side.S0[0] + side.S1[0] + side.E1[0] + side.E0[0]) * 0.25,
			(side.S0[1] + side.S1[1] + side.E1[1] + side.E0[1]) * 0.25,
			(side.S0[2] + side.S1[2] + side.E1[2] + side.E0[2]) * 0.25
		};

		double outDir[3] = {
			faceCenter[0] - centerHint[0],
			faceCenter[1] - centerHint[1],
			faceCenter[2] - centerHint[2]
		};

		if (vtkMath::Dot(side.normal, outDir) < 0.0) {
			for (int d = 0; d < 3; ++d) {
				side.normal[d] = -side.normal[d];
			}
		}

		return side;
	};

	std::vector<HostSide> hostSides(4);
	std::vector<HostSide> startTransitionSides(4);
	std::vector<HostSide> endTransitionSides(4);

	for (int s = 0; s < 4; ++s) {
		hostSides[s] = BuildHostSide(hostStartInner, hostEndInner, s, ports[mainA].normal, trunkCenter);
		startTransitionSides[s] = BuildHostSide(mainStart, hostStartInner, s, hostSides[s].normal, trunkCenter);
		endTransitionSides[s] = BuildHostSide(hostEndInner, mainEnd, s, hostSides[s].normal, trunkCenter);
	}

	auto BilinearOnSide = [](const HostSide& side, double u, double v, double out[3])
	{
		for (int d = 0; d < 3; ++d) {
			double pStart = side.S0[d] * (1.0 - v) + side.S1[d] * v;
			double pEnd = side.E0[d] * (1.0 - v) + side.E1[d] * v;
			out[d] = pStart * (1.0 - u) + pEnd * u;
		}
	};

	auto SmoothStep01 = [](double x) -> double
	{
		x = std::max(0.0, std::min(1.0, x));
		return x * x * (3.0 - 2.0 * x);
	};

	auto Bell01 = [&](double x) -> double
	{
		double s = SmoothStep01(x);
		return 4.0 * s * (1.0 - s);
	};

	auto CurvedOnSide = [&](const HostSide& side, double u, double v, double out[3])
	{
		BilinearOnSide(side, u, v, out);

		double bu = Bell01(u);
		double bv = Bell01(v);

		double bulge = std::min(side.lenV, maxPortDim) * 0.060 * bu * bv;

		for (int d = 0; d < 3; ++d) {
			out[d] += side.normal[d] * bulge;
		}
	};

	auto TransitionOnSide = [&](const HostSide& side, double u, double v, double out[3])
	{
		BilinearOnSide(side, u, v, out);

		double bu = Bell01(u);
		double bv = Bell01(v);

		double bulge = std::min(side.lenV, maxPortDim) * 0.030 * bu * bv;

		for (int d = 0; d < 3; ++d) {
			out[d] += side.normal[d] * bulge;
		}
	};

	auto PointToSideUV = [&](const HostSide& side, const double p[3], double& outU, double& outV)
	{
		double local[3] = {
			p[0] - side.S0[0],
			p[1] - side.S0[1],
			p[2] - side.S0[2]
		};

		outU = vtkMath::Dot(local, side.uAxis) / std::max(side.lenU, 1e-8);
		outV = vtkMath::Dot(local, side.vAxis) / std::max(side.lenV, 1e-8);
	};

	auto FindPathToHub = [&](int startNode, int endNode) -> std::vector<int>
	{
		std::vector<int> parent(nodeCount, -1);
		std::vector<int> visited(nodeCount, 0);

		std::queue<int> q;
		q.push(startNode);
		visited[startNode] = 1;

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (curr == endNode) {
				break;
			}

			if (curr < 0 || curr >= static_cast<int>(adjList.size())) {
				continue;
			}

			for (int next : adjList[curr]) {
				if (next < 0 || next >= nodeCount) continue;
				if (visited[next]) continue;

				visited[next] = 1;
				parent[next] = curr;
				q.push(next);
			}
		}

		std::vector<int> path;

		if (!visited[endNode]) {
			path.push_back(startNode);
			path.push_back(endNode);
			return path;
		}

		int curr = endNode;

		while (curr >= 0) {
			path.push_back(curr);

			if (curr == startNode) {
				break;
			}

			curr = parent[curr];
		}

		std::reverse(path.begin(), path.end());
		return path;
	};

	auto IntersectSegmentWithSidePlane = [&](const HostSide& side,
		const double a[3],
		const double b[3],
		double hit[3],
		double& outU,
		double& outV) -> bool
	{
		double faceCenter[3] = {
			(side.S0[0] + side.S1[0] + side.E1[0] + side.E0[0]) * 0.25,
			(side.S0[1] + side.S1[1] + side.E1[1] + side.E0[1]) * 0.25,
			(side.S0[2] + side.S1[2] + side.E1[2] + side.E0[2]) * 0.25
		};

		double ab[3] = {
			b[0] - a[0],
			b[1] - a[1],
			b[2] - a[2]
		};

		double denom = vtkMath::Dot(ab, side.normal);

		if (std::abs(denom) < 1e-9) {
			return false;
		}

		double diff[3] = {
			faceCenter[0] - a[0],
			faceCenter[1] - a[1],
			faceCenter[2] - a[2]
		};

		double t = vtkMath::Dot(diff, side.normal) / denom;

		if (t < -0.05 || t > 1.05) {
			return false;
		}

		t = std::max(0.0, std::min(1.0, t));

		for (int d = 0; d < 3; ++d) {
			hit[d] = a[d] + ab[d] * t;
		}

		PointToSideUV(side, hit, outU, outV);
		return true;
	};

	struct ProjectedHole
	{
		int portIndex = -1;
		int sideIndex = -1;

		double uc = 0.5;
		double vc = 0.5;

		// 原始骨架投影中心。后续 HoleProtect 可能会临时同步切线，
		// 但边界吸附必须回到这个中心附近，不能被另一孔拖到远端。
		double baseUc = 0.5;
		double baseVc = 0.5;

		bool lockPairU = false;
		double lockedPairU0 = 0.0;
		double lockedPairU1 = 1.0;

		double halfU = 0.05;
		double halfV = 0.05;

		// 期望孔半长/半宽：用于同面双孔贴边后，重新恢复输入端面尺寸。
		// halfU/halfV 可能会在避让循环里被压缩，但 desiredHalfU 不应被压缩。
		double desiredHalfU = 0.05;
		double desiredHalfV = 0.05;

		double u0 = 0.0;
		double u1 = 0.0;
		double v0 = 0.0;
		double v1 = 0.0;

		double holeRing[4][3];
		double collarRing[4][3];
	};

	std::vector<std::vector<ProjectedHole>> holesOnSide(4);

	auto IntersectSkeletonPathWithSide = [&](const HostSide& side,
		const BranchPort3D& port,
		double& outU,
		double& outV,
		double& outScore,
		bool requireInsideFace = false,
		bool allowFallbackProjection = true) -> bool
	{
		std::vector<int> path = FindPathToHub(port.anchorNode, hubNode);

		bool found = false;
		double bestScore = -1e100;
		double bestU = 0.5;
		double bestV = 0.5;

		for (int k = 0; k + 1 < static_cast<int>(path.size()); ++k) {
			double a[3], b[3];
			GetNodePoint(path[k], a);
			GetNodePoint(path[k + 1], b);

			double hit[3];
			double u = 0.5;
			double v = 0.5;

			if (!IntersectSegmentWithSidePlane(side, a, b, hit, u, v)) {
				continue;
			}

			double cu = std::max(0.0, std::min(1.0, u));
			double cv = std::max(0.0, std::min(1.0, v));
			double outside = std::abs(u - cu) + std::abs(v - cv);

			// 严格模式用于侧枝插入：交点必须真的落在当前 host 面内部。
			// 不能把已经穿出当前面的交点 clamp 回 0/1，否则插入点会跨到对面。
			if (requireInsideFace && outside > 1e-4) {
				continue;
			}

			double dirScore = vtkMath::Dot(port.dirFromHub, side.normal);
			double score = dirScore * 20.0 - outside * 60.0;

			if (u >= -0.05 && u <= 1.05 && v >= -0.05 && v <= 1.05) {
				score += 10.0;
			}

			if (score > bestScore) {
				bestScore = score;
				bestU = cu;
				bestV = cv;
				found = true;
			}
		}

		if (!found && !allowFallbackProjection) {
			return false;
		}

		if (!found) {
			double a[3];
			GetNodePoint(port.anchorNode, a);

			double hit[3];
			double u = 0.5;
			double v = 0.5;

			if (IntersectSegmentWithSidePlane(side, a, hubCenter, hit, u, v)) {
				double cu = std::max(0.0, std::min(1.0, u));
				double cv = std::max(0.0, std::min(1.0, v));
				double outside = std::abs(u - cu) + std::abs(v - cv);
				double dirScore = vtkMath::Dot(port.dirFromHub, side.normal);

				bestScore = dirScore * 20.0 - outside * 60.0 + 2.0;
				bestU = cu;
				bestV = cv;
				found = true;
			}
		}

		if (!found) {
			double faceCenter[3] = {
				(side.S0[0] + side.S1[0] + side.E1[0] + side.E0[0]) * 0.25,
				(side.S0[1] + side.S1[1] + side.E1[1] + side.E0[1]) * 0.25,
				(side.S0[2] + side.S1[2] + side.E1[2] + side.E0[2]) * 0.25
			};

			double p[3] = {
				port.center[0],
				port.center[1],
				port.center[2]
			};

			double diff[3] = {
				p[0] - faceCenter[0],
				p[1] - faceCenter[1],
				p[2] - faceCenter[2]
			};

			double h = vtkMath::Dot(diff, side.normal);

			for (int d = 0; d < 3; ++d) {
				p[d] -= side.normal[d] * h;
			}

			PointToSideUV(side, p, bestU, bestV);

			bestU = std::max(0.0, std::min(1.0, bestU));
			bestV = std::max(0.0, std::min(1.0, bestV));

			bestScore = vtkMath::Dot(port.dirFromHub, side.normal) * 5.0;
			found = true;
		}

		outU = bestU;
		outV = bestV;
		outScore = bestScore;

		return found;
	};

	// ============================================================
// 3. 支管投影：允许多个孔落在同一侧面；N>=5 时保护端口组和对称 pair
// ============================================================

	for (int p = 0; p < static_cast<int>(ports.size()); ++p) {
		if (p == mainA) {
			continue;
		}

		// 五端 1+2+2 模板：主端 pair 不能再侧插，后面会走端部分流。
		if (useGrouped122Template && IsMainEndGroupPort(p)) {
			continue;
		}

		// 五端 1+2+2 模板：本轮只做主干 + 主端 pair 的骨架连接。
		// 另一组 side pair 暂时不插入主干侧面，避免提前破坏 host 表面。
		int preMate = FindSymmetricPairMate(p);
		if (useGrouped122Template && !grouped122SideInsertMode && preMate >= 0) {
			if (false) qDebug() << "[3D branch][Grouped122] skip side pair insertion before remap, port =" << p;
			continue;
		}

		// 普通模式下 mainB 仍然是主干终点；虚拟主端组模式下，mainB 也要作为真实端口插入。
		if (!useVirtualMainEnd && p == mainB) {
			continue;
		}

		// 1+2+2 remap 后，side pair 会在较小编号的端口处一次性生成两个孔。
		// 较大编号的 mate 必须跳过，否则它会再次按“普通单孔”流程投影，
		// 产生一个多余孔/多余控制笼，并可能从侧柱桥接到错误的 host 边界。
		if (grouped122SideInsertMode && preMate >= 0 && p > preMate) {
			if (false) qDebug() << "[3D branch][Grouped122SidePair] skip duplicated mate insertion, port ="
				<< p << "mate =" << preMate;
			continue;
		}

		int bestSide = -1;
		double bestU = 0.5;
		double bestV = 0.5;
		double bestScore = -1e100;

		for (int s = 0; s < 4; ++s) {
			double u = 0.5;
			double v = 0.5;
			double score = -1e100;

			IntersectSkeletonPathWithSide(hostSides[s], ports[p], u, v, score);

			if (score > bestScore) {
				bestScore = score;
				bestSide = s;
				bestU = u;
				bestV = v;
			}
		}

		if (bestSide < 0) {
			continue;
		}

		// 五端 1+2+2 模板：侧向 pair 成对生成孔。
		// v4 修正：孔中心必须来自“各自骨架路径与 host side 的交点”，
		// 不能再用 pairCenter 抽象投影。pair 只负责强制同 side、同边界、按 U 顺序防交叉。
		int pairMate = FindSymmetricPairMate(p);
		if ((useGrouped122Template || grouped122SideInsertMode) && pairMate >= 0 && p < pairMate) {
			int pairSide = -1;
			double bestPairScore = -1e100;
			double hitU0 = 0.5;
			double hitV0 = 0.5;
			double hitU1 = 0.5;
			double hitV1 = 0.5;

			for (int s = 0; s < 4; ++s) {
				double u0 = 0.5;
				double v0 = 0.5;
				double score0 = -1e100;
				double u1 = 0.5;
				double v1 = 0.5;
				double score1 = -1e100;

				bool ok0 = IntersectSkeletonPathWithSide(hostSides[s], ports[p], u0, v0, score0, true, false);
				bool ok1 = IntersectSkeletonPathWithSide(hostSides[s], ports[pairMate], u1, v1, score1, true, false);

				if (!ok0 || !ok1) {
					continue;
				}

				double edgeScore0 = std::max(0.0, 1.0 - std::min(v0, 1.0 - v0) * 2.0);
				double edgeScore1 = std::max(0.0, 1.0 - std::min(v1, 1.0 - v1) * 2.0);
				double oppositeEdgeBonus = ((v0 < 0.5) != (v1 < 0.5)) ? 18.0 : -6.0;
				double sepU = std::abs(u1 - u0);

				// 以真实骨架交点为主；同面双侧枝优先落到同一 host side 的两条相对边界。
				// 后续贴边时只扩展 V 到边界，U 长度仍以输入端面尺寸为主。
				double pairScore = score0 + score1 + oppositeEdgeBonus + (edgeScore0 + edgeScore1) * 8.0 + sepU * 8.0;

				if (pairScore > bestPairScore) {
					bestPairScore = pairScore;
					pairSide = s;
					hitU0 = std::max(0.0, std::min(1.0, u0));
					hitV0 = std::max(0.0, std::min(1.0, v0));
					hitU1 = std::max(0.0, std::min(1.0, u1));
					hitV1 = std::max(0.0, std::min(1.0, v1));
				}
			}

			if (pairSide < 0) {
				if (false) qDebug() << "[3D branch][Grouped122SidePair] cannot find skeleton-side hits, fallback skip pair"
					<< p << pairMate;
				continue;
			}

			const HostSide& side = hostSides[pairSide];

			int leftPort = p;
			int rightPort = pairMate;
			double leftU = hitU0;
			double rightU = hitU1;
			double leftHitV = hitV0;
			double rightHitV = hitV1;

			if (rightU < leftU) {
				std::swap(leftPort, rightPort);
				std::swap(leftU, rightU);
				std::swap(leftHitV, rightHitV);
			}

			leftU = std::max(0.08, std::min(0.92, leftU));
			rightU = std::max(0.08, std::min(0.92, rightU));

			const double minPairSep = 0.18;
			if (rightU - leftU < minPairSep) {
				double mid = std::max(0.18, std::min(0.82, 0.5 * (leftU + rightU)));
				leftU = std::max(0.06, mid - 0.5 * minPairSep);
				rightU = std::min(0.94, mid + 0.5 * minPairSep);
			}

			// 两个端口落到同一 host side 时，直接分别预分配到两条相对边界。
			// 之前先把二者放到同一条边，后续避让循环会把 U 长度压小，导致孔长度小于输入端面。
			auto PortVMetric = [&](int portIndex) -> double
			{
				double diff[3] = {
					ports[portIndex].center[0] - side.S0[0],
					ports[portIndex].center[1] - side.S0[1],
					ports[portIndex].center[2] - side.S0[2]
				};
				return vtkMath::Dot(diff, side.vAxis);
			};

			double leftMetric = PortVMetric(leftPort);
			double rightMetric = PortVMetric(rightPort);

			// 侧枝孔贴哪一条边，优先由真实骨架与当前 host 面的交点决定。
			// 如果两个真实交点本来就在 V 的两侧，就保持各自所在侧；
			// 只有交点退化到同侧时，才用端面中心在 V 方向的相对位置兜底分配。
			double leftV = (leftHitV < 0.5) ? 0.055 : 0.945;
			double rightV = (rightHitV < 0.5) ? 0.055 : 0.945;
			if ((leftV < 0.5) == (rightV < 0.5)) {
				leftV = (leftMetric <= rightMetric) ? 0.055 : 0.945;
				rightV = (leftMetric <= rightMetric) ? 0.945 : 0.055;
			}

			auto MakePairHole = [&](int portIndex, double uc, double vc) -> ProjectedHole
			{
				ProjectedHole hole;
				hole.portIndex = portIndex;
				hole.sideIndex = pairSide;
				hole.uc = uc;
				hole.vc = vc;
				hole.baseUc = uc;
				hole.baseVc = vc;

				const double portLong = std::max(ports[portIndex].width, ports[portIndex].height);
				const double portShort = std::min(ports[portIndex].width, ports[portIndex].height);

				// U 是插入孔的长度方向。优先保证它接近输入端面长边；
				// MeasureRingExtent 受投影方向影响可能偏小，不能只依赖它。
				double extentU = std::max(
					MeasureRingExtent(ports[portIndex].ring, side.uAxis),
					portLong * 0.92
				);
				double extentV = std::max(
					MeasureRingExtent(ports[portIndex].ring, side.vAxis),
					portShort * 0.72
				);

				hole.halfU = 0.5 * extentU / std::max(side.lenU, 1e-8);
				hole.halfV = 0.5 * extentV / std::max(side.lenV, 1e-8);

				// 长度允许更大，避免孔被压成短小矩形；宽度仍然相对保守，因为贴边会扩到边界。
				hole.halfU = std::max(0.040, std::min(0.46, hole.halfU));
				hole.halfV = std::max(0.030, std::min(0.42, hole.halfV));

				hole.desiredHalfU = hole.halfU;
				hole.desiredHalfV = hole.halfV;

				return hole;
			};

			ProjectedHole leftHole = MakePairHole(leftPort, leftU, leftV);
			ProjectedHole rightHole = MakePairHole(rightPort, rightU, rightV);

			{
				const double margin = 0.035;
				const double localMinSpan = 0.045;
				double lo = margin;
				double hi = 1.0 - margin;
				double centerU = 0.5 * (leftU + rightU);
				centerU = std::max(lo, std::min(hi, centerU));

				double spanU = std::max(2.0 * leftHole.desiredHalfU, 2.0 * rightHole.desiredHalfU);
				spanU = std::max(spanU, std::abs(rightU - leftU) + 2.0 * std::max(leftHole.desiredHalfU, rightHole.desiredHalfU));
				spanU = std::max(spanU, localMinSpan);
				spanU = std::min(spanU, hi - lo);

				double lockU0 = centerU - 0.5 * spanU;
				double lockU1 = centerU + 0.5 * spanU;

				if (lockU0 < lo) { lockU1 += (lo - lockU0); lockU0 = lo; }
				if (lockU1 > hi) { lockU0 -= (lockU1 - hi); lockU1 = hi; }

				lockU0 = std::max(lo, std::min(hi, lockU0));
				lockU1 = std::max(lo, std::min(hi, lockU1));
				if (lockU1 - lockU0 < localMinSpan) {
					lockU1 = std::min(hi, lockU0 + localMinSpan);
					lockU0 = std::max(lo, lockU1 - localMinSpan);
				}

				leftHole.lockPairU = true;
				leftHole.lockedPairU0 = lockU0;
				leftHole.lockedPairU1 = lockU1;
				rightHole.lockPairU = true;
				rightHole.lockedPairU0 = lockU0;
				rightHole.lockedPairU1 = lockU1;
			}

			holesOnSide[pairSide].push_back(leftHole);
			holesOnSide[pairSide].push_back(rightHole);

			if (false) qDebug() << "[3D branch][Grouped122SidePair] skeleton-hit driven"
				<< "ports =" << leftPort << rightPort
				<< "side =" << pairSide
				<< "u =" << leftU << rightU
				<< "hitV =" << leftHitV << rightHitV
				<< "edgeV =" << leftV << rightV
				<< "desiredHalfU =" << holesOnSide[pairSide][holesOnSide[pairSide].size() - 2].desiredHalfU
				<< holesOnSide[pairSide][holesOnSide[pairSide].size() - 1].desiredHalfU;

			continue;
		}

		// 右侧多端口主端组：所有组内真实端口都靠近 end 插入，不再让其中一个消失。
		if (IsMainEndGroupPort(p)) {
			bestU = 0.88;
			bestV = (bestV < 0.5) ? 0.08 : 0.92;

			if (false) qDebug() << "[3D branch][EndGroupPort]"
				<< "port =" << p
				<< "side =" << bestSide
				<< "u =" << bestU
				<< "v =" << bestV;
		}

		// 对称 pair 的第二个端口强制使用第一个端口的 side，避免一个被识别到侧面、另一个到底面。
		int mate = FindSymmetricPairMate(p);
		if (mate >= 0 && mate < p) {
			bool forced = false;

			for (int s = 0; s < 4 && !forced; ++s) {
				for (const ProjectedHole& existing : holesOnSide[s]) {
					if (existing.portIndex == mate) {
						bestSide = s;

						double forcedU = 0.5;
						double forcedV = 0.5;
						double forcedScore = -1e100;
						bool forcedHitOk = IntersectSkeletonPathWithSide(hostSides[bestSide], ports[p], forcedU, forcedV, forcedScore, true, false);

						bestU = std::max(0.0, std::min(1.0, forcedU));
						// 第二个对称侧枝必须落在 mate 的对侧边，不能穿过当前面跑到同一侧/对面。
						bestV = existing.vc < 0.5 ? 0.92 : 0.08;
						if (!forcedHitOk) {
							bestU = existing.uc;
						}

						if (false) qDebug() << "[3D branch][SymPairForce]"
							<< "port =" << p
							<< "mate =" << mate
							<< "side =" << bestSide
							<< "u =" << bestU
							<< "v =" << bestV;

						forced = true;
						break;
					}
				}
			}
		}

		bestU = std::max(0.08, std::min(0.92, bestU));
		bestV = std::max(0.04, std::min(0.96, bestV));

		const HostSide& side = hostSides[bestSide];

		const double portLong = std::max(ports[p].width, ports[p].height);
		const double portShort = std::min(ports[p].width, ports[p].height);

		double extentU = std::max(
			MeasureRingExtent(ports[p].ring, side.uAxis),
			portLong * 0.92
		);

		double extentV = std::max(
			MeasureRingExtent(ports[p].ring, side.vAxis),
			portShort * 0.78
		);

		double halfU = 0.5 * extentU / std::max(side.lenU, 1e-8);
		double halfV = 0.5 * extentV / std::max(side.lenV, 1e-8);

		halfU = std::max(halfU, 0.035);
		halfV = std::max(halfV, 0.035);

		halfU = std::min(halfU, 0.46);
		halfV = std::min(halfV, 0.42);

		ProjectedHole hole;
		hole.portIndex = p;
		hole.sideIndex = bestSide;
		hole.uc = bestU;
		hole.vc = bestV;
		hole.baseUc = bestU;
		hole.baseVc = bestV;
		hole.halfU = halfU;
		hole.halfV = halfV;
		hole.desiredHalfU = halfU;
		hole.desiredHalfV = halfV;

		holesOnSide[bestSide].push_back(hole);
	}

	// ============================================================
	// 4. 同侧多孔尺寸避让 + 孔边界保护
	// ============================================================

	const double minGap = 0.025;
	const double minHalf = 0.020;
	const double minSpan = 0.045;
	const double edgeSnap = enablePortGrouping ? 0.160 : 0.120;
	const double cutSnapTol = 1e-5;

	auto Clamp01 = [](double x) -> double
	{
		return std::max(0.0, std::min(1.0, x));
	};

	auto CommitHoleBox = [&](ProjectedHole& h,
		double u0, double u1,
		double v0, double v1)
	{
		u0 = Clamp01(u0);
		u1 = Clamp01(u1);
		v0 = Clamp01(v0);
		v1 = Clamp01(v1);

		if (u1 < u0) std::swap(u0, u1);
		if (v1 < v0) std::swap(v0, v1);

		if (u1 - u0 < minSpan) {
			double c = Clamp01(0.5 * (u0 + u1));
			u0 = c - 0.5 * minSpan;
			u1 = c + 0.5 * minSpan;

			if (u0 < 0.0) { u1 -= u0; u0 = 0.0; }
			if (u1 > 1.0) { u0 -= (u1 - 1.0); u1 = 1.0; }

			u0 = Clamp01(u0);
			u1 = Clamp01(u1);
		}

		if (v1 - v0 < minSpan) {
			double c = Clamp01(0.5 * (v0 + v1));
			v0 = c - 0.5 * minSpan;
			v1 = c + 0.5 * minSpan;

			if (v0 < 0.0) { v1 -= v0; v0 = 0.0; }
			if (v1 > 1.0) { v0 -= (v1 - 1.0); v1 = 1.0; }

			v0 = Clamp01(v0);
			v1 = Clamp01(v1);
		}

		// 距离 side 边界很近时直接贴边，避免生成极窄四边形条带。
		if (u0 < edgeSnap) u0 = 0.0;
		if (1.0 - u1 < edgeSnap) u1 = 1.0;
		if (v0 < edgeSnap) v0 = 0.0;
		if (1.0 - v1 < edgeSnap) v1 = 1.0;

		if (u1 - u0 < minSpan) {
			if (u0 <= 0.0) u1 = std::min(1.0, minSpan);
			else if (u1 >= 1.0) u0 = std::max(0.0, 1.0 - minSpan);
			else {
				double c = Clamp01(0.5 * (u0 + u1));
				u0 = Clamp01(c - 0.5 * minSpan);
				u1 = Clamp01(c + 0.5 * minSpan);
			}
		}

		if (v1 - v0 < minSpan) {
			if (v0 <= 0.0) v1 = std::min(1.0, minSpan);
			else if (v1 >= 1.0) v0 = std::max(0.0, 1.0 - minSpan);
			else {
				double c = Clamp01(0.5 * (v0 + v1));
				v0 = Clamp01(c - 0.5 * minSpan);
				v1 = Clamp01(c + 0.5 * minSpan);
			}
		}

		h.u0 = u0;
		h.u1 = u1;
		h.v0 = v0;
		h.v1 = v1;

		h.uc = 0.5 * (h.u0 + h.u1);
		h.vc = 0.5 * (h.v0 + h.v1);
		h.halfU = 0.5 * (h.u1 - h.u0);
		h.halfV = 0.5 * (h.v1 - h.v0);
	};

	auto CommitHoleBoxSelectiveSnap = [&](ProjectedHole& h,
		double u0, double u1,
		double v0, double v1,
		bool allowUSnap,
		bool allowVSnap)
	{
		u0 = Clamp01(u0);
		u1 = Clamp01(u1);
		v0 = Clamp01(v0);
		v1 = Clamp01(v1);

		if (u1 < u0) std::swap(u0, u1);
		if (v1 < v0) std::swap(v0, v1);

		if (u1 - u0 < minSpan) {
			double c = Clamp01(0.5 * (u0 + u1));
			u0 = c - 0.5 * minSpan;
			u1 = c + 0.5 * minSpan;

			if (u0 < 0.0) { u1 -= u0; u0 = 0.0; }
			if (u1 > 1.0) { u0 -= (u1 - 1.0); u1 = 1.0; }

			u0 = Clamp01(u0);
			u1 = Clamp01(u1);
		}

		if (v1 - v0 < minSpan) {
			double c = Clamp01(0.5 * (v0 + v1));
			v0 = c - 0.5 * minSpan;
			v1 = c + 0.5 * minSpan;

			if (v0 < 0.0) { v1 -= v0; v0 = 0.0; }
			if (v1 > 1.0) { v0 -= (v1 - 1.0); v1 = 1.0; }

			v0 = Clamp01(v0);
			v1 = Clamp01(v1);
		}

		// 关键修复：1+2+2 同面侧枝只允许 V 边贴 host 的两条边界，
		// U 边只负责保持孔中心和输入端面长度，不能因为靠近整块面的 U=0/U=1 被二次吸附。
		if (allowUSnap) {
			if (u0 < edgeSnap) u0 = 0.0;
			if (1.0 - u1 < edgeSnap) u1 = 1.0;
		}
		if (allowVSnap) {
			if (v0 < edgeSnap) v0 = 0.0;
			if (1.0 - v1 < edgeSnap) v1 = 1.0;
		}

		if (u1 - u0 < minSpan) {
			if (u0 <= 0.0) u1 = std::min(1.0, minSpan);
			else if (u1 >= 1.0) u0 = std::max(0.0, 1.0 - minSpan);
			else {
				double c = Clamp01(0.5 * (u0 + u1));
				u0 = Clamp01(c - 0.5 * minSpan);
				u1 = Clamp01(c + 0.5 * minSpan);
			}
		}

		if (v1 - v0 < minSpan) {
			if (v0 <= 0.0) v1 = std::min(1.0, minSpan);
			else if (v1 >= 1.0) v0 = std::max(0.0, 1.0 - minSpan);
			else {
				double c = Clamp01(0.5 * (v0 + v1));
				v0 = Clamp01(c - 0.5 * minSpan);
				v1 = Clamp01(c + 0.5 * minSpan);
			}
		}

		h.u0 = u0;
		h.u1 = u1;
		h.v0 = v0;
		h.v1 = v1;

		h.uc = 0.5 * (h.u0 + h.u1);
		h.vc = 0.5 * (h.v0 + h.v1);
		h.halfU = 0.5 * (h.u1 - h.u0);
		h.halfV = 0.5 * (h.v1 - h.v0);
	};

	auto CenteredIntervalPreservePairU = [&](const ProjectedHole& h,
		double center,
		double span)->std::array<double, 2>
	{
		center = Clamp01(center);
		span = std::max(minSpan, std::min(0.98, span));

		if (h.lockPairU) {
			const double margin = std::max(minGap, 0.020);
			double lo = margin;
			double hi = 1.0 - margin;

			if (span > hi - lo) {
				span = hi - lo;
			}

			double a = center - 0.5 * span;
			double b = center + 0.5 * span;

			if (a < lo) { b += (lo - a); a = lo; }
			if (b > hi) { a -= (b - hi); b = hi; }

			a = std::max(lo, a);
			b = std::min(hi, b);
			return { Clamp01(a), Clamp01(b) };
		}

		double a = center - 0.5 * span;
		double b = center + 0.5 * span;

		if (a < 0.0) { b -= a; a = 0.0; }
		if (b > 1.0) { a -= (b - 1.0); b = 1.0; }

		a = Clamp01(a);
		b = Clamp01(b);
		return { a, b };
	};

	auto CommitHoleBoxPreservePairU = [&](ProjectedHole& h,
		double u0, double u1,
		double v0, double v1)
	{
		if (h.lockPairU) {
			u0 = h.lockedPairU0;
			u1 = h.lockedPairU1;
		}

		CommitHoleBoxSelectiveSnap(h, u0, u1, v0, v1, !h.lockPairU, true);
	};

	auto FitHoleFixedCenter = [&](ProjectedHole& h)
	{
		h.uc = Clamp01(h.uc);
		h.vc = Clamp01(h.vc);

		h.halfU = std::max(minHalf, std::min(0.46, std::max(h.halfU, h.desiredHalfU)));
		h.halfV = std::max(minHalf, std::min(0.46, std::max(h.halfV, h.desiredHalfV)));

		CommitHoleBox(h,
			h.uc - h.halfU,
			h.uc + h.halfU,
			h.vc - h.halfV,
			h.vc + h.halfV);
	};

	auto FitHoleFixedCenterPreservePairU = [&](ProjectedHole& h)
	{
		h.uc = Clamp01(h.uc);
		h.vc = Clamp01(h.vc);

		h.halfU = std::max(minHalf, std::min(0.46, std::max(h.halfU, h.desiredHalfU)));
		h.halfV = std::max(minHalf, std::min(0.46, std::max(h.halfV, h.desiredHalfV)));

		auto uRange = CenteredIntervalPreservePairU(h, h.uc, 2.0 * h.halfU);
		CommitHoleBoxPreservePairU(h,
			uRange[0],
			uRange[1],
			h.vc - h.halfV,
			h.vc + h.halfV);
	};

	auto CutInsideInterval = [&](double cut, double a0, double a1) -> bool
	{
		return cut > a0 + edgeSnap * 0.15 && cut < a1 - edgeSnap * 0.15;
	};

	auto HoleCutsWouldSplit = [&](const ProjectedHole& a, const ProjectedHole& b) -> bool
	{
		return
			CutInsideInterval(b.u0, a.u0, a.u1) ||
			CutInsideInterval(b.u1, a.u0, a.u1) ||
			CutInsideInterval(a.u0, b.u0, b.u1) ||
			CutInsideInterval(a.u1, b.u0, b.u1) ||
			CutInsideInterval(b.v0, a.v0, a.v1) ||
			CutInsideInterval(b.v1, a.v0, a.v1) ||
			CutInsideInterval(a.v0, b.v0, b.v1) ||
			CutInsideInterval(a.v1, b.v0, b.v1);
	};

	for (int s = 0; s < 4; ++s) {
		auto& holes = holesOnSide[s];

		if (holes.empty()) {
			continue;
		}

		for (ProjectedHole& h : holes) {
			FitHoleFixedCenterPreservePairU(h);
		}

		for (int iter = 0; iter < 24; ++iter) {
			bool changed = false;

			for (int i = 0; i < static_cast<int>(holes.size()); ++i) {
				for (int j = i + 1; j < static_cast<int>(holes.size()); ++j) {
					ProjectedHole& a = holes[i];
					ProjectedHole& b = holes[j];

					double overlapU = std::min(a.u1, b.u1) - std::max(a.u0, b.u0);
					double overlapV = std::min(a.v1, b.v1) - std::max(a.v0, b.v0);

					if (overlapU <= -minGap || overlapV <= -minGap) {
						continue;
					}

					double du = std::abs(a.uc - b.uc);
					double dv = std::abs(a.vc - b.vc);

					double allowU = std::max(0.0, du - minGap);
					double allowV = std::max(0.0, dv - minGap);

					double sumHalfU = a.halfU + b.halfU;
					double sumHalfV = a.halfV + b.halfV;

					double scaleU = (sumHalfU > 1e-8) ? allowU / sumHalfU : 0.0;
					double scaleV = (sumHalfV > 1e-8) ? allowV / sumHalfV : 0.0;

					scaleU = std::max(0.0, std::min(1.0, scaleU));
					scaleV = std::max(0.0, std::min(1.0, scaleV));

					if (scaleU >= scaleV) {
						double scale = std::max(0.55, scaleU * 0.96);
						a.halfU *= scale;
						b.halfU *= scale;
					}
					else {
						double scale = std::max(0.55, scaleV * 0.96);
						a.halfV *= scale;
						b.halfV *= scale;
					}

					FitHoleFixedCenterPreservePairU(a);
					FitHoleFixedCenterPreservePairU(b);
					changed = true;
				}
			}

			if (!changed) {
				break;
			}
		}

		// 防止同一 side 上某个孔的 cut 穿过另一个孔内部。
		for (int iter = 0; iter < 8; ++iter) {
			bool changed = false;

			for (int i = 0; i < static_cast<int>(holes.size()); ++i) {
				for (int j = i + 1; j < static_cast<int>(holes.size()); ++j) {
					ProjectedHole& a = holes[i];
					ProjectedHole& b = holes[j];

					if (!HoleCutsWouldSplit(a, b)) {
						continue;
					}

					double du = std::abs(a.uc - b.uc);
					double dv = std::abs(a.vc - b.vc);

					int aMate = FindSymmetricPairMate(a.portIndex);
					int bMate = FindSymmetricPairMate(b.portIndex);
					bool isSymmetricPair = (aMate == b.portIndex) || (bMate == a.portIndex);

					bool aNearLowerV = a.v0 < edgeSnap || a.vc < 0.20;
					bool aNearUpperV = (1.0 - a.v1) < edgeSnap || a.vc > 0.80;
					bool bNearLowerV = b.v0 < edgeSnap || b.vc < 0.20;
					bool bNearUpperV = (1.0 - b.v1) < edgeSnap || b.vc > 0.80;

					bool bothBoundaryHoles = (aNearLowerV || aNearUpperV) && (bNearLowerV || bNearUpperV);

					if (isSymmetricPair || bothBoundaryHoles) {
						double u0 = std::min(a.u0, b.u0);
						double u1 = std::max(a.u1, b.u1);

						double av0 = a.v0;
						double av1 = a.v1;
						double bv0 = b.v0;
						double bv1 = b.v1;

						if (a.vc <= b.vc) {
							av0 = 0.0;
							av1 = std::max(av1, minSpan);
							bv1 = 1.0;
							bv0 = std::min(bv0, 1.0 - minSpan);
						}
						else {
							bv0 = 0.0;
							bv1 = std::max(bv1, minSpan);
							av1 = 1.0;
							av0 = std::min(av0, 1.0 - minSpan);
						}

						CommitHoleBoxPreservePairU(a, u0, u1, av0, av1);
						CommitHoleBoxPreservePairU(b, u0, u1, bv0, bv1);

						if (false) qDebug() << "[3D branch][HoleProtect] boundary pair align on side"
							<< s
							<< "ports =" << a.portIndex << b.portIndex
							<< "u =" << u0 << u1
							<< "a.v =" << a.v0 << a.v1
							<< "b.v =" << b.v0 << b.v1;
					}
					else if (dv >= du) {
						double u0 = std::min(a.u0, b.u0);
						double u1 = std::max(a.u1, b.u1);

						CommitHoleBoxPreservePairU(a, u0, u1, a.v0, a.v1);
						CommitHoleBoxPreservePairU(b, u0, u1, b.v0, b.v1);

						if (false) qDebug() << "[3D branch][HoleProtect] sync U cuts on side"
							<< s
							<< "ports =" << a.portIndex << b.portIndex
							<< "u0 u1 =" << u0 << u1;
					}
					else {
						double v0 = std::min(a.v0, b.v0);
						double v1 = std::max(a.v1, b.v1);

						CommitHoleBoxPreservePairU(a, a.u0, a.u1, v0, v1);
						CommitHoleBoxPreservePairU(b, b.u0, b.u1, v0, v1);

						if (false) qDebug() << "[3D branch][HoleProtect] sync V cuts on side"
							<< s
							<< "ports =" << a.portIndex << b.portIndex
							<< "v0 v1 =" << v0 << v1;
					}

					changed = true;
				}
			}

			if (!changed) {
				break;
			}
		}

		// ------------------------------------------------------------
		// Same-face side-branch boundary snap, v2
		//
		// 只把“每个孔离 host 边界最近的那一条边”吸附到边界；
		// 不再把两个孔强制共用同一个 U 范围，也不强制一个贴 v=0、另一个贴 v=1。
		//
		// 这样可以避免：
		// 1) 右侧孔的左边被错误拉到最左侧；
		// 2) 孔的 U 长度被压小，导致侧枝根部扭曲；
		// 3) 为了边界吸附产生过多无意义切线。
		//
		// 尺寸策略：
		// - U/V span 优先恢复到输入端面投影得到的 desiredHalfU / desiredHalfV；
		// - 在空间允许时保持输入端面尺寸；
		// - 只有被吸附的那条边移动到 0 或 1，另一侧按尺寸自然展开。
		// ------------------------------------------------------------
		if (!holes.empty()) {
			auto DesiredSpanU = [&](const ProjectedHole& h) -> double
			{
				double span = std::max(h.u1 - h.u0, 2.0 * h.desiredHalfU);
				span = std::max(span, minSpan);
				span = std::min(span, 0.86);
				return span;
			};

			auto DesiredSpanV = [&](const ProjectedHole& h) -> double
			{
				double span = std::max(h.v1 - h.v0, 2.0 * h.desiredHalfV);
				span = std::max(span, minSpan);
				span = std::min(span, 0.86);
				return span;
			};

			auto AttachHoleToNearestBoundary = [&](ProjectedHole& h)
			{
				// v3：吸附边界以后，重新按输入端面尺寸恢复孔大小。
				// grouped122 / symmetric pair 的支撑侧枝，只吸附 host 面的 V 两侧边界，
				// 不允许再把右侧孔的 U 边吸附到整块面的最左/最右边界。
				// U 方向作为孔的长度方向，优先保持输入端面长边对应长度；
				// V 方向作为孔的宽度方向，吸附到最近 V 边界后保持输入端面短边对应宽度。

				double spanU = DesiredSpanU(h);
				double spanV = DesiredSpanV(h);

				// 尺寸下限来自当前孔，尺寸上限来自可用空间。
				// 这样前面的避让不会把孔压短；只有真的空间不够时才缩小。
				spanU = std::max(spanU, 2.0 * h.desiredHalfU);
				spanV = std::max(spanV, 2.0 * h.desiredHalfV);
				spanU = std::max(spanU, minSpan);
				spanV = std::max(spanV, minSpan);
				spanU = std::min(spanU, 0.92);
				spanV = std::min(spanV, 0.92);

				double centerU = Clamp01(h.baseUc);
				double centerV = Clamp01(h.baseVc);

				int mate = FindSymmetricPairMate(h.portIndex);
				// HostSide 的 U 轴就是主骨架连接方向（ringA -> ringB）。
				// 因此“与主骨架方向平行的两条边”就是 V=0 / V=1 两条边。
				// 成对侧枝/1+2+2 remap 后的侧枝必须只贴这两条边；
				// 普通非成对孔仍然走最近边逻辑，避免误伤其它单侧孔。
				bool forceMainParallelBoundaryOnly = grouped122SideInsertMode || (mate >= 0) || (holes.size() == 2);

				// 先按照输入端面尺寸恢复原始孔框。
				auto BuildCenteredBox = [&]()
				{
					double u0 = centerU - 0.5 * spanU;
					double u1 = centerU + 0.5 * spanU;
					double v0 = centerV - 0.5 * spanV;
					double v1 = centerV + 0.5 * spanV;

					if (u0 < 0.0) { u1 -= u0; u0 = 0.0; }
					if (u1 > 1.0) { u0 -= (u1 - 1.0); u1 = 1.0; }
					if (v0 < 0.0) { v1 -= v0; v0 = 0.0; }
					if (v1 > 1.0) { v0 -= (v1 - 1.0); v1 = 1.0; }

					u0 = Clamp01(u0); u1 = Clamp01(u1);
					v0 = Clamp01(v0); v1 = Clamp01(v1);

					return std::array<double, 4>{ u0, u1, v0, v1 };
				};

				auto box = BuildCenteredBox();
				double u0 = box[0];
				double u1 = box[1];
				double v0 = box[2];
				double v1 = box[3];

				// 注意：1+2+2 的同面侧枝，本质是背景支撑贴 host 两侧边界，
				// 所以这里只在 V0/V1 里选最近边界。不能再选 U0/U1。
				int attachEdge = 2;
				double distU0 = u0;
				double distU1 = 1.0 - u1;
				double distV0 = v0;
				double distV1 = 1.0 - v1;
				double bestDist = 1e100;

				if (forceMainParallelBoundaryOnly) {
					if (distV0 <= distV1) {
						attachEdge = 2;
						bestDist = distV0;
					}
					else {
						attachEdge = 3;
						bestDist = distV1;
					}
				}
				else {
					attachEdge = 0;
					bestDist = distU0;
					if (distU1 < bestDist) { bestDist = distU1; attachEdge = 1; }
					if (distV0 < bestDist) { bestDist = distV0; attachEdge = 2; }
					if (distV1 < bestDist) { bestDist = distV1; attachEdge = 3; }
				}

				double snapThreshold = grouped122SideInsertMode ? 0.48 : 0.30;
				bool shouldSnap = forceMainParallelBoundaryOnly || bestDist < snapThreshold;

				// 吸附后仍然保留 spanU/spanV：
				// - 如果吸附 V 边，U 继续围绕骨架投影中心展开，V 从边界向内展开；
				// - 如果普通模式吸附 U 边，V 继续围绕骨架投影中心展开，U 从边界向内展开。
				if (shouldSnap) {
					if (attachEdge == 2) {
						v0 = 0.0;
						v1 = std::min(1.0, spanV);
						auto uRange = CenteredIntervalPreservePairU(h, centerU, spanU);
						u0 = uRange[0];
						u1 = uRange[1];
					}
					else if (attachEdge == 3) {
						v1 = 1.0;
						v0 = std::max(0.0, 1.0 - spanV);
						auto uRange = CenteredIntervalPreservePairU(h, centerU, spanU);
						u0 = uRange[0];
						u1 = uRange[1];
					}
					else if (attachEdge == 0) {
						u0 = 0.0;
						u1 = std::min(1.0, spanU);
						v0 = centerV - 0.5 * spanV;
						v1 = centerV + 0.5 * spanV;
						if (v0 < 0.0) { v1 -= v0; v0 = 0.0; }
						if (v1 > 1.0) { v0 -= (v1 - 1.0); v1 = 1.0; }
					}
					else {
						u1 = 1.0;
						u0 = std::max(0.0, 1.0 - spanU);
						v0 = centerV - 0.5 * spanV;
						v1 = centerV + 0.5 * spanV;
						if (v0 < 0.0) { v1 -= v0; v0 = 0.0; }
						if (v1 > 1.0) { v0 -= (v1 - 1.0); v1 = 1.0; }
					}
				}

				u0 = Clamp01(u0); u1 = Clamp01(u1);
				v0 = Clamp01(v0); v1 = Clamp01(v1);

				// 如果因为边界空间不足导致 span 被压缩，这里只接受真实空间限制；
				// 不再接受前面 HoleProtect 造成的无意义压缩。
				// grouped122 / symmetric pair 只允许 V 吸附，禁止 CommitHoleBox 的 U 方向二次贴边，
				// 否则左孔的右边会因为 u1 接近 1 被吸到整块 host 的右边界。
				CommitHoleBoxSelectiveSnap(h, u0, u1, v0, v1, !forceMainParallelBoundaryOnly, true);

				if (false) qDebug() << "[3D branch][SameFaceBoundarySnap] nearest-edge preserve input size"
					<< "side =" << s
					<< "port =" << h.portIndex
					<< "forceVOnly =" << forceMainParallelBoundaryOnly
					<< "attachEdge =" << attachEdge
					<< "bestDist =" << bestDist
					<< "snap =" << shouldSnap
					<< "base =" << centerU << centerV
					<< "u =" << h.u0 << h.u1
					<< "v =" << h.v0 << h.v1
					<< "desiredSpan =" << spanU << spanV
					<< "actualSpan =" << (h.u1 - h.u0) << (h.v1 - h.v0);
			};

			// 优先只处理明确的 symmetric pair / grouped122 remapped pair。
			bool snappedPair = false;
			for (int i = 0; i < static_cast<int>(holes.size()) && !snappedPair; ++i) {
				for (int j = i + 1; j < static_cast<int>(holes.size()) && !snappedPair; ++j) {
					int mateI = FindSymmetricPairMate(holes[i].portIndex);
					int mateJ = FindSymmetricPairMate(holes[j].portIndex);
					bool isPair = (mateI == holes[j].portIndex) || (mateJ == holes[i].portIndex);

					bool isRemappedSidePair =
						(static_cast<int>(ports.size()) == 4 &&
						((holes[i].portIndex == 2 && holes[j].portIndex == 3) ||
							(holes[i].portIndex == 3 && holes[j].portIndex == 2)));

					if (!isPair && !isRemappedSidePair) {
						continue;
					}

					AttachHoleToNearestBoundary(holes[i]);
					AttachHoleToNearestBoundary(holes[j]);
					double u0 = std::min(holes[i].u0, holes[j].u0);
					double u1 = std::max(holes[i].u1, holes[j].u1);
					CommitHoleBoxPreservePairU(holes[i], u0, u1, holes[i].v0, holes[i].v1);
					CommitHoleBoxPreservePairU(holes[j], u0, u1, holes[j].v0, holes[j].v1);
					snappedPair = true;
				}
			}

			// 没有显式 pair 但同一面刚好两个孔时，也按最近边吸附。
			if (!snappedPair && holes.size() == 2) {
				AttachHoleToNearestBoundary(holes[0]);
				AttachHoleToNearestBoundary(holes[1]);
			}
			else if (!snappedPair && holes.size() == 1) {
				AttachHoleToNearestBoundary(holes[0]);
			}
		}


		// ------------------------------------------------------------
		// Final hole-integrity normalize, v4
		// ------------------------------------------------------------
		// SameFaceBoundarySnap 会保留每个孔自己的最近边界吸附和输入尺寸，
		// 但同一 host side 上多个孔如果在一个参数方向分离、另一个参数方向重叠，
		// 后续公共 U/V 切线会把“一个孔的边界线”贯穿到“另一个孔的内部”。
		//
		// 这里在真正生成 holeRing 之前做最后一次只针对拓扑安全的修正：
		// - 上下分离的孔，如果 U 边界互相切入，则统一 U 区间；
		// - 左右分离的孔，如果 V 边界互相切入，则统一 V 区间；
		// - 已经吸附到 V=0 / V=1 的边界不会被破坏，因为只改 U；
		// - 已经吸附到 U=0 / U=1 的边界不会被破坏，因为只改 V。
		//
		// 目标不是重新排布孔，而是保证每个孔内部没有任何公共 cut，
		// 从而插入孔始终保持干净的四节点边界，避免 T-junction 和 OCC 补洞崩溃。
		if (holes.size() >= 2) {
			auto BoundaryInsideInterval = [&](double cut, double a0, double a1) -> bool
			{
				return cut > a0 + cutSnapTol && cut < a1 - cutSnapTol;
			};

			auto HasBoundaryCutInsideU = [&](const ProjectedHole& a, const ProjectedHole& b) -> bool
			{
				return BoundaryInsideInterval(a.u0, b.u0, b.u1) ||
					BoundaryInsideInterval(a.u1, b.u0, b.u1) ||
					BoundaryInsideInterval(b.u0, a.u0, a.u1) ||
					BoundaryInsideInterval(b.u1, a.u0, a.u1);
			};

			auto HasBoundaryCutInsideV = [&](const ProjectedHole& a, const ProjectedHole& b) -> bool
			{
				return BoundaryInsideInterval(a.v0, b.v0, b.v1) ||
					BoundaryInsideInterval(a.v1, b.v0, b.v1) ||
					BoundaryInsideInterval(b.v0, a.v0, a.v1) ||
					BoundaryInsideInterval(b.v1, a.v0, a.v1);
			};

			auto SeparatedOrTouching = [&](double a0, double a1, double b0, double b1) -> bool
			{
				return a1 <= b0 + minGap * 0.35 || b1 <= a0 + minGap * 0.35;
			};

			auto ClampUnionInterval = [&](double& a0, double& a1, double b0, double b1)
			{
				double n0 = std::min(a0, b0);
				double n1 = std::max(a1, b1);

				// 极端情况下 union 过大时仍然保留一点外部面片，避免整张 side 被孔吃掉。
				const double maxSpan = 0.98;
				if (n1 - n0 > maxSpan) {
					double c = 0.5 * (n0 + n1);
					n0 = c - 0.5 * maxSpan;
					n1 = c + 0.5 * maxSpan;
				}

				if (n0 < 0.0) { n1 -= n0; n0 = 0.0; }
				if (n1 > 1.0) { n0 -= (n1 - 1.0); n1 = 1.0; }

				a0 = Clamp01(n0);
				a1 = Clamp01(n1);
			};

			for (int iter = 0; iter < 12; ++iter) {
				bool changed = false;

				for (int i = 0; i < static_cast<int>(holes.size()); ++i) {
					for (int j = i + 1; j < static_cast<int>(holes.size()); ++j) {
						ProjectedHole& a = holes[i];
						ProjectedHole& b = holes[j];

						bool uConflict = HasBoundaryCutInsideU(a, b);
						bool vConflict = HasBoundaryCutInsideV(a, b);

						bool separatedV = SeparatedOrTouching(a.v0, a.v1, b.v0, b.v1);
						bool separatedU = SeparatedOrTouching(a.u0, a.u1, b.u0, b.u1);

						// 上下分离：只统一 U，保留各自 V 边界吸附。
						if (uConflict && separatedV) {
							double u0 = a.u0;
							double u1 = a.u1;
							ClampUnionInterval(u0, u1, b.u0, b.u1);

							if (std::abs(a.u0 - u0) > cutSnapTol || std::abs(a.u1 - u1) > cutSnapTol ||
								std::abs(b.u0 - u0) > cutSnapTol || std::abs(b.u1 - u1) > cutSnapTol)
							{
								CommitHoleBoxPreservePairU(a, u0, u1, a.v0, a.v1);
								CommitHoleBoxPreservePairU(b, u0, u1, b.v0, b.v1);
								changed = true;

								if (false) qDebug() << "[3D branch][HoleIntegrity] sync U after snap"
									<< "side =" << s
									<< "ports =" << a.portIndex << b.portIndex
									<< "u =" << u0 << u1
									<< "a.v =" << a.v0 << a.v1
									<< "b.v =" << b.v0 << b.v1;
							}
						}

						// 左右分离：只统一 V，保留各自 U 边界吸附。
						if (vConflict && separatedU) {
							double v0 = a.v0;
							double v1 = a.v1;
							ClampUnionInterval(v0, v1, b.v0, b.v1);

							if (std::abs(a.v0 - v0) > cutSnapTol || std::abs(a.v1 - v1) > cutSnapTol ||
								std::abs(b.v0 - v0) > cutSnapTol || std::abs(b.v1 - v1) > cutSnapTol)
							{
								CommitHoleBoxPreservePairU(a, a.u0, a.u1, v0, v1);
								CommitHoleBoxPreservePairU(b, b.u0, b.u1, v0, v1);
								changed = true;

								if (false) qDebug() << "[3D branch][HoleIntegrity] sync V after snap"
									<< "side =" << s
									<< "ports =" << a.portIndex << b.portIndex
									<< "v =" << v0 << v1
									<< "a.u =" << a.u0 << a.u1
									<< "b.u =" << b.u0 << b.u1;
							}
						}
					}
				}

				if (!changed) {
					break;
				}
			}
		}

	}


	// ------------------------------------------------------------
	// Global U integrity pass, v4
	// ------------------------------------------------------------
	// commonUCuts 会贯穿四个 host side。即使两个孔不在同一 side，
	// A 孔的 U 边界也可能成为 B 孔内部的贯穿 cut。
	// 所以在所有 side 的吸附完成后、真正生成 holeRing 之前，
	// 再做一次跨 side 的 U 区间同步，确保 common U cuts 不切穿任何孔。
	{
		auto BoundaryInsideInterval = [&](double cut, double a0, double a1) -> bool
		{
			return cut > a0 + cutSnapTol && cut < a1 - cutSnapTol;
		};

		auto HasBoundaryCutInsideU = [&](const ProjectedHole& a, const ProjectedHole& b) -> bool
		{
			return BoundaryInsideInterval(a.u0, b.u0, b.u1) ||
				BoundaryInsideInterval(a.u1, b.u0, b.u1) ||
				BoundaryInsideInterval(b.u0, a.u0, a.u1) ||
				BoundaryInsideInterval(b.u1, a.u0, a.u1);
		};

		auto ClampUnionInterval = [&](double& a0, double& a1, double b0, double b1)
		{
			double n0 = std::min(a0, b0);
			double n1 = std::max(a1, b1);

			const double maxSpan = 0.98;
			if (n1 - n0 > maxSpan) {
				double c = 0.5 * (n0 + n1);
				n0 = c - 0.5 * maxSpan;
				n1 = c + 0.5 * maxSpan;
			}

			if (n0 < 0.0) { n1 -= n0; n0 = 0.0; }
			if (n1 > 1.0) { n0 -= (n1 - 1.0); n1 = 1.0; }

			a0 = Clamp01(n0);
			a1 = Clamp01(n1);
		};

		for (int iter = 0; iter < 12; ++iter) {
			bool changed = false;

			for (int sa = 0; sa < 4; ++sa) {
				for (int ia = 0; ia < static_cast<int>(holesOnSide[sa].size()); ++ia) {
					for (int sb = sa; sb < 4; ++sb) {
						int jbStart = (sb == sa) ? ia + 1 : 0;
						for (int jb = jbStart; jb < static_cast<int>(holesOnSide[sb].size()); ++jb) {
							ProjectedHole& a = holesOnSide[sa][ia];
							ProjectedHole& b = holesOnSide[sb][jb];

							if (!HasBoundaryCutInsideU(a, b)) {
								continue;
							}

							double u0 = a.u0;
							double u1 = a.u1;
							ClampUnionInterval(u0, u1, b.u0, b.u1);

							if (std::abs(a.u0 - u0) > cutSnapTol || std::abs(a.u1 - u1) > cutSnapTol ||
								std::abs(b.u0 - u0) > cutSnapTol || std::abs(b.u1 - u1) > cutSnapTol)
							{
								CommitHoleBoxPreservePairU(a, u0, u1, a.v0, a.v1);
								CommitHoleBoxPreservePairU(b, u0, u1, b.v0, b.v1);
								changed = true;

								if (false) qDebug() << "[3D branch][HoleIntegrity] global sync U"
									<< "sideA =" << sa
									<< "portA =" << a.portIndex
									<< "sideB =" << sb
									<< "portB =" << b.portIndex
									<< "u =" << u0 << u1;
							}
						}
					}
				}
			}

			if (!changed) {
				break;
			}
		}
	}

	// 所有孔的最终 UV 已确定后，再生成真实 holeRing / collarRing。
	for (int s = 0; s < 4; ++s) {
		for (ProjectedHole& h : holesOnSide[s]) {
			const HostSide& side = hostSides[s];

			CurvedOnSide(side, h.u0, h.v0, h.holeRing[0]);
			CurvedOnSide(side, h.u1, h.v0, h.holeRing[1]);
			CurvedOnSide(side, h.u1, h.v1, h.holeRing[2]);
			CurvedOnSide(side, h.u0, h.v1, h.holeRing[3]);

			BranchPort3D& port = ports[h.portIndex];

			double holeCenter[3];
			ComputeRingCenter(h.holeRing, holeCenter);

			double toPort[3] = {
				port.center[0] - holeCenter[0],
				port.center[1] - holeCenter[1],
				port.center[2] - holeCenter[2]
			};
			NormalizeSafe(toPort, port.dirFromHub);

			double collarDir[3] = {
				toPort[0],
				toPort[1],
				toPort[2]
			};
			NormalizeSafe(collarDir, side.normal);

			double outDot = vtkMath::Dot(collarDir, side.normal);
			const double minOutDot = 0.35;

			if (outDot < minOutDot) {
				double add = minOutDot - outDot;
				collarDir[0] += side.normal[0] * add;
				collarDir[1] += side.normal[1] * add;
				collarDir[2] += side.normal[2] * add;
				NormalizeSafe(collarDir, side.normal);
			}

			double collarLift = std::min(port.width, port.height) * 0.10;
			collarLift = std::max(collarLift, maxPortDim * 0.020);
			collarLift = std::min(collarLift, maxPortDim * 0.075);

			for (int k = 0; k < 4; ++k) {
				for (int d = 0; d < 3; ++d) {
					h.collarRing[k][d] = h.holeRing[k][d] + collarDir[d] * collarLift;
				}
			}
		}
	}

	// ============================================================
	// 5. 构造受保护的 commonUCuts / sideVCuts
	// ============================================================

	auto AddUniqueCut = [&](std::vector<double>& cuts, double value)
	{
		value = Clamp01(value);

		for (double& c : cuts) {
			if (std::abs(c - value) < cutSnapTol) {
				c = 0.5 * (c + value);
				return;
			}
		}

		cuts.push_back(value);
	};

	auto IsHoleBoundaryCut = [&](const ProjectedHole& h, double cut, bool isU) -> bool
	{
		if (isU) {
			return std::abs(cut - h.u0) < cutSnapTol ||
				std::abs(cut - h.u1) < cutSnapTol;
		}

		return std::abs(cut - h.v0) < cutSnapTol ||
			std::abs(cut - h.v1) < cutSnapTol;
	};

	auto CutInsideHoleInterval = [&](const ProjectedHole& h, double cut, bool isU) -> bool
	{
		if (isU) {
			return cut > h.u0 + cutSnapTol && cut < h.u1 - cutSnapTol;
		}

		return cut > h.v0 + cutSnapTol && cut < h.v1 - cutSnapTol;
	};

	std::vector<ProjectedHole> allHolesForCommonU;

	for (int s = 0; s < 4; ++s) {
		for (const ProjectedHole& h : holesOnSide[s]) {
			allHolesForCommonU.push_back(h);
		}
	}

	std::vector<double> commonUCuts;
	std::vector<std::vector<double>> sideVCuts(4);

	AddUniqueCut(commonUCuts, 0.0);
	AddUniqueCut(commonUCuts, 1.0);

	for (const ProjectedHole& h : allHolesForCommonU) {
		AddUniqueCut(commonUCuts, h.u0);
		AddUniqueCut(commonUCuts, h.u1);
	}

	{
		std::vector<double> clean;
		clean.reserve(commonUCuts.size());

		for (double cut : commonUCuts) {
			cut = Clamp01(cut);

			bool isBoundaryOfSomeHole = false;
			bool isInsideAnyOtherHole = false;

			for (const ProjectedHole& h : allHolesForCommonU) {
				bool isBoundaryOfThisHole = IsHoleBoundaryCut(h, cut, true);
				if (isBoundaryOfThisHole) {
					isBoundaryOfSomeHole = true;
				}

				// 关键修复：不能因为 cut 是 A 孔边界，就允许它穿过 B 孔内部。
				// 这会把 B 孔边界切成多段，和四节点 h.holeRing 产生 T-junction。
				if (!isBoundaryOfThisHole && CutInsideHoleInterval(h, cut, true)) {
					isInsideAnyOtherHole = true;
					if (false) qDebug() << "[3D branch][HoleProtect] remove common U cut crossing another hole"
						<< "port =" << h.portIndex
						<< "cut =" << cut
						<< "hole u =" << h.u0 << h.u1
						<< "hole v =" << h.v0 << h.v1
						<< "boundaryOfSomeHole =" << isBoundaryOfSomeHole;
					break;
				}
			}

			bool keep = !isInsideAnyOtherHole;

			// 0/1 必须保留，用于整块 host side 外边界。
			if (std::abs(cut) < cutSnapTol || std::abs(cut - 1.0) < cutSnapTol) {
				keep = true;
			}

			if (keep) {
				bool exists = false;
				for (double c : clean) {
					if (std::abs(c - cut) < cutSnapTol) {
						exists = true;
						break;
					}
				}

				if (!exists) {
					clean.push_back(cut);
				}
			}
		}

		std::sort(clean.begin(), clean.end());
		commonUCuts.swap(clean);
	}

	for (int s = 0; s < 4; ++s) {
		AddUniqueCut(sideVCuts[s], 0.0);
		AddUniqueCut(sideVCuts[s], 1.0);

		for (const ProjectedHole& h : holesOnSide[s]) {
			AddUniqueCut(sideVCuts[s], h.v0);
			AddUniqueCut(sideVCuts[s], h.v1);
		}

		std::vector<double> cleanV;
		cleanV.reserve(sideVCuts[s].size());

		for (double cut : sideVCuts[s]) {
			cut = Clamp01(cut);

			bool isBoundaryOfSomeHole = false;
			bool isInsideAnyOtherHole = false;

			for (const ProjectedHole& h : holesOnSide[s]) {
				bool isBoundaryOfThisHole = IsHoleBoundaryCut(h, cut, false);
				if (isBoundaryOfThisHole) {
					isBoundaryOfSomeHole = true;
				}

				// 同理，V cut 也不能穿过另一个孔内部。
				if (!isBoundaryOfThisHole && CutInsideHoleInterval(h, cut, false)) {
					isInsideAnyOtherHole = true;
					if (false) qDebug() << "[3D branch][HoleProtect] remove side V cut crossing another hole"
						<< "side =" << s
						<< "port =" << h.portIndex
						<< "cut =" << cut
						<< "hole u =" << h.u0 << h.u1
						<< "hole v =" << h.v0 << h.v1
						<< "boundaryOfSomeHole =" << isBoundaryOfSomeHole;
					break;
				}
			}

			bool keep = !isInsideAnyOtherHole;

			// 0/1 必须保留，用于整块 host side 外边界。
			if (std::abs(cut) < cutSnapTol || std::abs(cut - 1.0) < cutSnapTol) {
				keep = true;
			}

			if (keep) {
				bool exists = false;
				for (double c : cleanV) {
					if (std::abs(c - cut) < cutSnapTol) {
						exists = true;
						break;
					}
				}

				if (!exists) {
					cleanV.push_back(cut);
				}
			}
		}

		std::sort(cleanV.begin(), cleanV.end());
		sideVCuts[s].swap(cleanV);
	}

	for (int s = 0; s < 4; ++s) {
		for (const ProjectedHole& h : holesOnSide[s]) {
			int uBoundaryCount = 0;
			int vBoundaryCount = 0;
			int uInsideCount = 0;
			int vInsideCount = 0;

			for (double u : commonUCuts) {
				if (std::abs(u - h.u0) < cutSnapTol || std::abs(u - h.u1) < cutSnapTol) {
					uBoundaryCount++;
				}
				else if (u > h.u0 + cutSnapTol && u < h.u1 - cutSnapTol) {
					uInsideCount++;
				}
			}

			for (double v : sideVCuts[s]) {
				if (std::abs(v - h.v0) < cutSnapTol || std::abs(v - h.v1) < cutSnapTol) {
					vBoundaryCount++;
				}
				else if (v > h.v0 + cutSnapTol && v < h.v1 - cutSnapTol) {
					vInsideCount++;
				}
			}

			if (false) qDebug() << "[3D branch][HoleDebug]"
				<< "side =" << s
				<< "port =" << h.portIndex
				<< "u =" << h.u0 << h.u1
				<< "v =" << h.v0 << h.v1
				<< "uBoundary =" << uBoundaryCount
				<< "vBoundary =" << vBoundaryCount
				<< "uInside =" << uInsideCount
				<< "vInside =" << vInsideCount;
		}
	}

	// ============================================================
		// 6. 生成 hubPoly：中间开孔段 + 左右 transition 全 quad
		// ============================================================

	vtkSmartPointer<vtkPolyData> hubPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> hubPts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> hubPolys = vtkSmartPointer<vtkCellArray>::New();

	auto InsertHubPoint = [&](const double p[3]) -> vtkIdType
	{
		for (vtkIdType id = 0; id < hubPts->GetNumberOfPoints(); ++id) {
			double q[3];
			hubPts->GetPoint(id, q);

			double dx = p[0] - q[0];
			double dy = p[1] - q[1];
			double dz = p[2] - q[2];

			if (dx * dx + dy * dy + dz * dz < 1e-12) {
				return id;
			}
		}

		return hubPts->InsertNextPoint(p);
	};

	auto QuadNormal = [](const double p0[3], const double p1[3], const double p2[3], double n[3])
	{
		double a[3] = {
			p1[0] - p0[0],
			p1[1] - p0[1],
			p1[2] - p0[2]
		};

		double b[3] = {
			p2[0] - p1[0],
			p2[1] - p1[1],
			p2[2] - p1[2]
		};

		vtkMath::Cross(a, b, n);

		if (vtkMath::Norm(n) > 1e-8) {
			vtkMath::Normalize(n);
		}
	};

	auto AddHubQuadOriented = [&](const HostSide& side,
		const double p0[3],
		const double p1[3],
		const double p2[3],
		const double p3[3])
	{
		double n[3];
		QuadNormal(p0, p1, p2, n);

		bool flip = vtkMath::Dot(n, side.normal) < 0.0;

		vtkIdType id0 = InsertHubPoint(p0);
		vtkIdType id1 = InsertHubPoint(p1);
		vtkIdType id2 = InsertHubPoint(p2);
		vtkIdType id3 = InsertHubPoint(p3);

		hubPolys->InsertNextCell(4);

		if (!flip) {
			hubPolys->InsertCellPoint(id0);
			hubPolys->InsertCellPoint(id1);
			hubPolys->InsertCellPoint(id2);
			hubPolys->InsertCellPoint(id3);
		}
		else {
			hubPolys->InsertCellPoint(id0);
			hubPolys->InsertCellPoint(id3);
			hubPolys->InsertCellPoint(id2);
			hubPolys->InsertCellPoint(id1);
		}
	};

	auto CellInsideHole = [](const std::vector<ProjectedHole>& holes,
		double u0, double u1, double v0, double v1) -> bool
	{
		double uc = 0.5 * (u0 + u1);
		double vc = 0.5 * (v0 + v1);

		for (const ProjectedHole& h : holes) {
			if (uc > h.u0 + 1e-7 && uc < h.u1 - 1e-7 &&
				vc > h.v0 + 1e-7 && vc < h.v1 - 1e-7)
			{
				return true;
			}
		}

		return false;
	};

	for (int s = 0; s < 4; ++s) {
		const HostSide& side = hostSides[s];

		const std::vector<double>& uCuts = commonUCuts;
		const std::vector<double>& vCuts = sideVCuts[s];

		for (int i = 0; i + 1 < static_cast<int>(uCuts.size()); ++i) {
			for (int j = 0; j + 1 < static_cast<int>(vCuts.size()); ++j) {
				double u0 = uCuts[i];
				double u1 = uCuts[i + 1];
				double v0 = vCuts[j];
				double v1 = vCuts[j + 1];

				if (u1 - u0 < 1e-6 || v1 - v0 < 1e-6) {
					continue;
				}

				if (CellInsideHole(holesOnSide[s], u0, u1, v0, v1)) {
					continue;
				}

				double p0[3], p1[3], p2[3], p3[3];

				CurvedOnSide(side, u0, v0, p0);
				CurvedOnSide(side, u1, v0, p1);
				CurvedOnSide(side, u1, v1, p2);
				CurvedOnSide(side, u0, v1, p3);

				AddHubQuadOriented(side, p0, p1, p2, p3);
			}
		}
	}

	auto AddHubTriOriented = [&](const HostSide& side,
		const double p0[3],
		const double p1[3],
		const double p2[3])
	{
		double a[3] = {
			p1[0] - p0[0],
			p1[1] - p0[1],
			p1[2] - p0[2]
		};

		double b[3] = {
			p2[0] - p1[0],
			p2[1] - p1[1],
			p2[2] - p1[2]
		};

		double n[3];
		vtkMath::Cross(a, b, n);

		if (vtkMath::Norm(n) > 1e-8) {
			vtkMath::Normalize(n);
		}

		bool flip = vtkMath::Dot(n, side.normal) < 0.0;

		vtkIdType id0 = InsertHubPoint(p0);
		vtkIdType id1 = InsertHubPoint(p1);
		vtkIdType id2 = InsertHubPoint(p2);

		hubPolys->InsertNextCell(3);

		if (!flip) {
			hubPolys->InsertCellPoint(id0);
			hubPolys->InsertCellPoint(id1);
			hubPolys->InsertCellPoint(id2);
		}
		else {
			hubPolys->InsertCellPoint(id0);
			hubPolys->InsertCellPoint(id2);
			hubPolys->InsertCellPoint(id1);
		}
	};

	auto BuildRealTransitionPatch = [&](const HostSide& side,
		const std::vector<double>& vCuts,
		bool innerAtEnd)
	{
		if (vCuts.size() < 2) {
			return;
		}

		// inner：靠中间开孔段的一侧，允许多节点。
		// outer：靠输入端的一侧，只允许 outer0 / outer1 两个边界点。
		double innerU = innerAtEnd ? 1.0 : 0.0;
		double midU = innerAtEnd ? 0.58 : 0.42;
		double outerU = innerAtEnd ? 0.0 : 1.0;

		double outer0[3];
		double outer1[3];

		// 注意：outer 边界只生成两个点。
		// 这里绝对不能按 vCuts 在 outerU 上生成多个点。
		TransitionOnSide(side, outerU, 0.0, outer0);
		TransitionOnSide(side, outerU, 1.0, outer1);

		std::vector<std::array<double, 3>> innerPts;
		std::vector<std::array<double, 3>> midPts;

		innerPts.reserve(vCuts.size());
		midPts.reserve(vCuts.size());

		for (double v : vCuts) {
			double pInner[3];
			double pMid[3];

			// innerPts 在多节点边界上。
			TransitionOnSide(side, innerU, v, pInner);

			// midPts 是过渡段内部点。
			// 这些点不能在 outer 边界上。
			TransitionOnSide(side, midU, v, pMid);

			innerPts.push_back({ pInner[0], pInner[1], pInner[2] });
			midPts.push_back({ pMid[0], pMid[1], pMid[2] });
		}

		const int M = static_cast<int>(innerPts.size());

		// 情况 1：inner 本来也是 2 点，直接一个 quad 接 outer0-outer1。
		if (M == 2) {
			double i0[3] = { innerPts[0][0], innerPts[0][1], innerPts[0][2] };
			double i1[3] = { innerPts[1][0], innerPts[1][1], innerPts[1][2] };

			// 外边界 edge 是 outer0-outer1，不插任何点。
			AddHubQuadOriented(side, i0, i1, outer1, outer0);
			return;
		}

		// 情况 2：inner 有 3 个及以上点。
		// 令内部过渡点 P1...P(K)，K = M - 2。
		// 这些 P 全部在 midU 内部线上，不在 outer 边界上。
		//
		// 例如 M = 4：
		// inner: I0 -- I1 -- I2 -- I3
		// mid:        P1 -- P2
		// outer: O0 ---------------- O1
		//
		// 生成：
		// Q0 = I0, I1, P1, O0
		// Q1 = I1, I2, P2, P1
		// Q2 = I2, I3, O1, P2
		//
		// 这正是你画的四点到两点过渡结构。
		for (int k = 0; k < M - 1; ++k) {
			double i0[3] = { innerPts[k][0], innerPts[k][1], innerPts[k][2] };
			double i1[3] = { innerPts[k + 1][0], innerPts[k + 1][1], innerPts[k + 1][2] };

			if (k == 0) {
				// 第一段：I0-I1-P1-O0
				double p1[3] = {
					midPts[1][0],
					midPts[1][1],
					midPts[1][2]
				};

				AddHubQuadOriented(side, i0, i1, p1, outer0);
			}
			else if (k == M - 2) {
				// 最后一段：I(M-2)-I(M-1)-O1-P(K)
				double pk[3] = {
					midPts[M - 2][0],
					midPts[M - 2][1],
					midPts[M - 2][2]
				};

				AddHubQuadOriented(side, i0, i1, outer1, pk);
			}
			else {
				// 中间段：Ik-I(k+1)-P(k+1)-Pk
				double pk[3] = {
					midPts[k][0],
					midPts[k][1],
					midPts[k][2]
				};

				double pk1[3] = {
					midPts[k + 1][0],
					midPts[k + 1][1],
					midPts[k + 1][2]
				};

				AddHubQuadOriented(side, i0, i1, pk1, pk);
			}
		}

		// 下面处理 inner 多点边界和 outer 单边界之间剩余的“右侧收口区域”。
		// 对于最常见的 M = 4，内部只有 P1/P2，直接一个 quad：
		// O0-O1-P2-P1
		//
		// 注意：outer 边界仍然只有 outer0-outer1，没有任何新增点落在这个边上。
		const int K = M - 2;

		if (K == 1) {
			// M = 3，五边区域无法全部四边形，保留一个小三角。
			double p1[3] = {
				midPts[1][0],
				midPts[1][1],
				midPts[1][2]
			};

			AddHubTriOriented(side, outer0, outer1, p1);
		}
		else if (K == 2) {
			// M = 4，完全四边形收口。
			double p1[3] = {
				midPts[1][0],
				midPts[1][1],
				midPts[1][2]
			};

			double p2[3] = {
				midPts[2][0],
				midPts[2][1],
				midPts[2][2]
			};

			AddHubQuadOriented(side, outer0, outer1, p2, p1);
		}
		else {
			// M > 4 的复杂多孔同侧情况。
			// 这里用一个内部 center 做局部收口，不在 outer 边界插点。
			// 这不是最优 quad 网格，但比从 outer0 扇形连所有点正确：
			// - outer0/outer1 边界保持完整；
			// - 新增点 center 在过渡段内部；
			// - 只在复杂情况使用少量 triangle。
			double center[3] = { 0.0, 0.0, 0.0 };

			center[0] += outer0[0] + outer1[0];
			center[1] += outer0[1] + outer1[1];
			center[2] += outer0[2] + outer1[2];

			for (int k = 1; k <= M - 2; ++k) {
				center[0] += midPts[k][0];
				center[1] += midPts[k][1];
				center[2] += midPts[k][2];
			}

			double denom = static_cast<double>(K + 2);
			center[0] /= denom;
			center[1] /= denom;
			center[2] /= denom;

			// outer 边界三角：outer0-outer1-center。
			AddHubTriOriented(side, outer0, outer1, center);

			// outer1 -> P(K)
			{
				double pk[3] = {
					midPts[M - 2][0],
					midPts[M - 2][1],
					midPts[M - 2][2]
				};

				AddHubTriOriented(side, outer1, pk, center);
			}

			// P(k+1) -> P(k)
			for (int k = M - 2; k > 1; --k) {
				double pa[3] = {
					midPts[k][0],
					midPts[k][1],
					midPts[k][2]
				};

				double pb[3] = {
					midPts[k - 1][0],
					midPts[k - 1][1],
					midPts[k - 1][2]
				};

				AddHubTriOriented(side, pa, pb, center);
			}

			// P1 -> outer0
			{
				double p1[3] = {
					midPts[1][0],
					midPts[1][1],
					midPts[1][2]
				};

				AddHubTriOriented(side, p1, outer0, center);
			}
		}
	};
	for (int s = 0; s < 4; ++s) {
		BuildRealTransitionPatch(startTransitionSides[s], sideVCuts[s], true);
		BuildRealTransitionPatch(endTransitionSides[s], sideVCuts[s], false);
	}

	if (hubPolys->GetNumberOfCells() > 0) {
		hubPoly->SetPoints(hubPts);
		hubPoly->SetPolys(hubPolys);
		m_hubCages.push_back(hubPoly);
	}

	// ============================================================
	// 7. 支管插入：所有孔全部连接
	// ============================================================

	for (int s = 0; s < 4; ++s) {
		for (const ProjectedHole& h : holesOnSide[s]) {
			BranchPort3D& port = ports[h.portIndex];

			int uInsideCount = 0;
			int vInsideCount = 0;

			for (double u : commonUCuts) {
				if (u > h.u0 + cutSnapTol && u < h.u1 - cutSnapTol) {
					uInsideCount++;
				}
			}

			for (double v : sideVCuts[s]) {
				if (v > h.v0 + cutSnapTol && v < h.v1 - cutSnapTol) {
					vInsideCount++;
				}
			}

			if (uInsideCount > 0 || vInsideCount > 0) {
				if (false) qDebug() << "[3D branch][HoleError] unsafe branch insertion, keep it for diagnosis"
					<< "side =" << s
					<< "port =" << h.portIndex
					<< "uInside =" << uInsideCount
					<< "vInside =" << vInsideCount;
			}

			AddSegmentBetweenRings(h.holeRing, h.collarRing);

			double bestTubeRing[4][3];

			BridgeScore bestScoreLocal;
			bestScoreLocal.minEdgeDot = -1e100;
			bestScoreLocal.avgEdgeDot = -1e100;
			bestScoreLocal.totalLength = 1e100;
			bestScoreLocal.maxDiagDiff = 1e100;

			for (int family = 0; family < 2; ++family) {
				bool reverse = family == 1;

				for (int shift = 0; shift < 4; ++shift) {
					double candidate[4][3];
					BuildTubeCandidateRing(port.ring, reverse, shift, candidate);

					BridgeScore score = EvaluateBridge(h.collarRing, candidate);

					if (IsBetterScore(score, bestScoreLocal)) {
						bestScoreLocal = score;
						CopyRing(bestTubeRing, candidate);
					}
				}
			}

			// 侧枝输入端不再额外添加 support 过渡段。
			// 保留 h.holeRing -> h.collarRing 作为插入孔缓冲，
			// 然后直接从 collarRing 连接到输入端面，避免把侧枝根部切得过碎。
			AddSegmentBetweenRings(h.collarRing, bestTubeRing);

			if (bestScoreLocal.minEdgeDot < 0.05) {
				if (false) qDebug() << "[3D branch] branch bridge twist risk, minEdgeDot ="
					<< bestScoreLocal.minEdgeDot
					<< "port =" << h.portIndex
					<< "side =" << s;
			}
		}
	}

	// ============================================================
	// 7.1 五端 1+2+2：主端 pair 的主段连接
	// ============================================================
	// 当前版本只实现：
	//   mainStart -> mainEnd 的主干控制笼已经在前面生成；
	//   mainEnd + 主端 pair 两个真实端面 调用二维 NBranch 核心生成局部分流。
	// 另一组 side pair 暂时不插入主干侧面，上面的投影循环已经跳过。
	if (useGrouped122Template && mainEndGroupPorts.size() == 2) {
		std::vector<NBranchInputPort> fanPorts;
		fanPorts.reserve(3);

		NBranchInputPort virtualMainEnd;
		virtualMainEnd.sourceCageIndex = -1;
		virtualMainEnd.sourceIsStart = false;
		virtualMainEnd.isVirtual = true;
		CopyRing(virtualMainEnd.ring, mainEnd);
		for (int d = 0; d < 3; ++d) {
			virtualMainEnd.preferredDir[d] = -trunkAxis[d];
		}
		virtualMainEnd.hasPreferredDir = true;
		fanPorts.push_back(virtualMainEnd);

		for (int idx = 0; idx < 2; ++idx) {
			int pidx = mainEndGroupPorts[idx];

			NBranchInputPort branchPort;
			branchPort.sourceCageIndex = ports[pidx].sel.cageIndex;
			branchPort.sourceIsStart = ports[pidx].sel.isStart;
			branchPort.isVirtual = false;
			CopyRing(branchPort.ring, ports[pidx].ring);

			double branchCenter[3];
			ComputeRingCenter(branchPort.ring, branchCenter);
			double virtualCenter[3];
			ComputeRingCenter(virtualMainEnd.ring, virtualCenter);

			for (int d = 0; d < 3; ++d) {
				branchPort.preferredDir[d] = branchCenter[d] - virtualCenter[d];
			}
			if (vtkMath::Norm(branchPort.preferredDir) > 1e-8) {
				vtkMath::Normalize(branchPort.preferredDir);
				branchPort.hasPreferredDir = true;
			}

			fanPorts.push_back(branchPort);
		}

		double virtualCenter[3];
		double b0Center[3];
		double b1Center[3];
		ComputeRingCenter(fanPorts[0].ring, virtualCenter);
		ComputeRingCenter(fanPorts[1].ring, b0Center);
		ComputeRingCenter(fanPorts[2].ring, b1Center);

		double fanHubCenter[3] = {
			virtualCenter[0] * 0.34 + b0Center[0] * 0.33 + b1Center[0] * 0.33,
			virtualCenter[1] * 0.34 + b0Center[1] * 0.33 + b1Center[1] * 0.33,
			virtualCenter[2] * 0.34 + b0Center[2] * 0.33 + b1Center[2] * 0.33
		};

		double v0[3] = {
			b0Center[0] - virtualCenter[0],
			b0Center[1] - virtualCenter[1],
			b0Center[2] - virtualCenter[2]
		};
		double v1[3] = {
			b1Center[0] - virtualCenter[0],
			b1Center[1] - virtualCenter[1],
			b1Center[2] - virtualCenter[2]
		};
		double fanNormalUp[3];
		vtkMath::Cross(v0, v1, fanNormalUp);
		NormalizeSafe(fanNormalUp, ports[mainA].normal);

		NBranchBuildOptions fanOptions;
		fanOptions.addSdf = true;
		fanOptions.rebuildNow = false;
		fanOptions.debugLog = true;
		fanOptions.forceHubCenter = true;
		fanOptions.hubCenter[0] = fanHubCenter[0];
		fanOptions.hubCenter[1] = fanHubCenter[1];
		fanOptions.hubCenter[2] = fanHubCenter[2];
		fanOptions.forceNormalUp = true;
		fanOptions.normalUp[0] = fanNormalUp[0];
		fanOptions.normalUp[1] = fanNormalUp[1];
		fanOptions.normalUp[2] = fanNormalUp[2];
		fanOptions.radiusScale = 0.75;
		fanOptions.minRadiusScale = 0.045;

		bool ok = GenerateNBranchFromRings(fanPorts, fanOptions);
		if (false) qDebug() << "[3D branch][Grouped122MainFanOut] reused GenerateNBranchFromRings"
			<< "ok =" << ok
			<< "virtual = mainEnd"
			<< "ports =" << mainEndGroupPorts[0] << mainEndGroupPorts[1]
			<< "sidePairSkipped =" << symmetricSidePairs[0].first << symmetricSidePairs[0].second
			<< "hubCenter =" << fanHubCenter[0] << fanHubCenter[1] << fanHubCenter[2];
	}

	if (false) qDebug() << "[3D branch] generated full quad transition multi-hole topology"
		<< "N =" << N
		<< "main ports =" << mainA << mainB
		<< "startBuffer =" << startBufferIndex
		<< "trunk =" << trunkCageIndex
		<< "endBuffer =" << endBufferIndex
		<< "hubNode =" << hubNode
		<< "center =" << hubCenter[0] << hubCenter[1] << hubCenter[2];

	RebuildGlobalMesh();
}

// ============================================================
// N=3 二维局部三角 hub
// 目的：替换旧的 N=3 二维中心 hub / 强制复制 A 面逻辑。
// 新规则：
//   1) A 端只决定三角 hub 的底边方向，不再复制 A 输入端面；
//   2) B0/B1 决定三角尖端方向；
//   3) 骨架交点 J 只作为局部锚点，不再强制等于 hub 六点几何中心；
//   4) 三角尺寸只由局部端口尺寸和 J 到端口距离限制，避免被外轮廓射线拉成长条。
// ============================================================
bool CageInteractorStyle::GenerateNBranchTriHub2D()
{
	const int N = static_cast<int>(m_selectedFacesForBranch.size());
	if (N != 3) {
		return false;
	}

	struct TriPort
	{
		SelectedFace sel;
		double ring[4][3];
		double oppositeRing[4][3];
		double center[3] = { 0.0, 0.0, 0.0 };
		double oppositeCenter[3] = { 0.0, 0.0, 0.0 };
		double normal[3] = { 1.0, 0.0, 0.0 };
		double branchDir[3] = { 1.0, 0.0, 0.0 };
		double width = 1.0;
		double depth = 1.0;
		int anchorNode = -1;
	};

	auto CopyRingLocal = [](double dst[4][3], const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	};

	auto RingCenterLocal = [](const double ring[4][3], double c[3])
	{
		c[0] = 0.0;
		c[1] = 0.0;
		c[2] = 0.0;
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				c[d] += ring[i][d] * 0.25;
			}
		}
	};

	auto DistanceLocal = [](const double a[3], const double b[3]) -> double
	{
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto ClampLocal = [](double v, double lo, double hi) -> double
	{
		return std::max(lo, std::min(v, hi));
	};

	auto NormalizeSafeLocal = [](double v[3], const double fallback[3])
	{
		if (vtkMath::Norm(v) < 1e-8) {
			if (fallback && vtkMath::Norm(fallback) > 1e-8) {
				v[0] = fallback[0];
				v[1] = fallback[1];
				v[2] = fallback[2];
			}
			else {
				v[0] = 1.0;
				v[1] = 0.0;
				v[2] = 0.0;
			}
		}
		vtkMath::Normalize(v);
	};

	auto SignedAreaLocal = [](const double ring[4][3], const double normal[3]) -> double
	{
		double area[3] = { 0.0, 0.0, 0.0 };
		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;
			double c[3];
			vtkMath::Cross(ring[i], ring[j], c);
			area[0] += c[0];
			area[1] += c[1];
			area[2] += c[2];
		}
		return vtkMath::Dot(area, normal);
	};

	auto ReverseRingLocal = [](double ring[4][3])
	{
		double tmp[4][3];
		for (int i = 0; i < 4; ++i) {
			const int j = (4 - i) % 4;
			for (int d = 0; d < 3; ++d) {
				tmp[i][d] = ring[j][d];
			}
		}
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				ring[i][d] = tmp[i][d];
			}
		}
	};

	auto FindNearestIndexLocal = [](const double p[3], const double ring[4][3]) -> int
	{
		int bestId = 0;
		double bestD2 = std::numeric_limits<double>::max();
		for (int i = 0; i < 4; ++i) {
			double d2 = vtkMath::Distance2BetweenPoints(p, ring[i]);
			if (d2 < bestD2) {
				bestD2 = d2;
				bestId = i;
			}
		}
		return bestId;
	};

	struct BridgeScoreLocal
	{
		double totalLength = 0.0;
		double maxDiagDiff = 0.0;
		double minEdgeDot = 1.0;
		double avgEdgeDot = 0.0;
	};

	auto DotEdgeNormalizedLocal = [](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double
	{
		double va[3] = { a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2] };
		double vb[3] = { b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2] };
		if (vtkMath::Norm(va) < 1e-10 || vtkMath::Norm(vb) < 1e-10) {
			return 0.0;
		}
		vtkMath::Normalize(va);
		vtkMath::Normalize(vb);
		return vtkMath::Dot(va, vb);
	};

	auto EvaluateBridgeLocal = [&](const double startRing[4][3], const double endRing[4][3]) -> BridgeScoreLocal
	{
		BridgeScoreLocal score;
		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;
			double dot = DotEdgeNormalizedLocal(startRing[i], startRing[j], endRing[i], endRing[j]);
			score.minEdgeDot = std::min(score.minEdgeDot, dot);
			score.avgEdgeDot += dot * 0.25;
			score.totalLength += DistanceLocal(startRing[i], endRing[i]);
			double diagA = DistanceLocal(startRing[i], endRing[j]);
			double diagB = DistanceLocal(startRing[j], endRing[i]);
			score.maxDiagDiff = std::max(score.maxDiagDiff, std::abs(diagA - diagB));
		}
		return score;
	};

	auto IsBetterBridgeScoreLocal = [](const BridgeScoreLocal& a, const BridgeScoreLocal& b) -> bool
	{
		if (a.minEdgeDot > b.minEdgeDot + 1e-6) return true;
		if (a.minEdgeDot < b.minEdgeDot - 1e-6) return false;
		if (a.avgEdgeDot > b.avgEdgeDot + 1e-6) return true;
		if (a.avgEdgeDot < b.avgEdgeDot - 1e-6) return false;
		if (a.maxDiagDiff < b.maxDiagDiff - 1e-6) return true;
		if (a.maxDiagDiff > b.maxDiagDiff + 1e-6) return false;
		return a.totalLength < b.totalLength;
	};

	auto BuildTubeCandidateRingLocal = [](const double inputRing[4][3], bool reverse, int shift, double outRing[4][3])
	{
		static const int forwardBase[4] = { 0, 1, 2, 3 };
		static const int reverseBase[4] = { 0, 3, 2, 1 };
		const int* base = reverse ? reverseBase : forwardBase;
		for (int i = 0; i < 4; ++i) {
			int id = base[(i + shift) % 4];
			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = inputRing[id][d];
			}
		}
	};

	// ------------------------------------------------------------
	// 1. 读取三个选中端面
	// ------------------------------------------------------------
	TriPort ports[3];

	for (int k = 0; k < 3; ++k) {
		const SelectedFace& sel = m_selectedFacesForBranch[k];
		if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
			if (false) qDebug() << "[TriHub2D] invalid selected cage index" << sel.cageIndex;
			return false;
		}

		const CageSegment& cage = m_allCages[sel.cageIndex];
		const double(*ring)[3] = sel.isStart ? cage.pStart : cage.pEnd;
		const double(*oppositeRing)[3] = sel.isStart ? cage.pEnd : cage.pStart;

		ports[k].sel = sel;
		CopyRingLocal(ports[k].ring, ring);
		CopyRingLocal(ports[k].oppositeRing, oppositeRing);
		RingCenterLocal(ports[k].ring, ports[k].center);
		RingCenterLocal(ports[k].oppositeRing, ports[k].oppositeCenter);

		for (int d = 0; d < 3; ++d) {
			ports[k].normal[d] = ports[k].center[d] - ports[k].oppositeCenter[d];
		}
		NormalizeSafeLocal(ports[k].normal, nullptr);

		const double e01 = DistanceLocal(ports[k].ring[0], ports[k].ring[1]);
		const double e12 = DistanceLocal(ports[k].ring[1], ports[k].ring[2]);
		const double e23 = DistanceLocal(ports[k].ring[2], ports[k].ring[3]);
		const double e30 = DistanceLocal(ports[k].ring[3], ports[k].ring[0]);
		const double lenA = 0.5 * (e01 + e23);
		const double lenB = 0.5 * (e12 + e30);

		ports[k].width = std::max(lenA, lenB);
		ports[k].depth = std::max(std::min(lenA, lenB), ports[k].width * 0.20);

		if (!m_skelNodes.empty()) {
			double tmpCenter[3] = { 0.0, 0.0, 0.0 };
			double tmpNormal[3] = { 1.0, 0.0, 0.0 };
			ports[k].anchorNode = FindNearestSkeletonNodeToSelectedFace(sel, tmpCenter, tmpNormal);
		}
	}

	// ------------------------------------------------------------
	// 2. 找真实骨架交点 J；如果失败，用三端面中心的几何中心兜底。
	// ------------------------------------------------------------
	double J[3] = { 0.0, 0.0, 0.0 };
	int junctionNode = -1;
	bool hasJunction = false;
	std::vector<std::vector<int>> adjList;

	if (!m_skelNodes.empty() && !m_skelEdges.empty()) {
		adjList = BuildSkeletonAdjacency();
		const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);

		auto GraphDistances = [&](int startNode) -> std::vector<int>
		{
			std::vector<int> dist(nodeCount, -1);
			if (startNode < 0 || startNode >= nodeCount) {
				return dist;
			}

			std::queue<int> q;
			dist[startNode] = 0;
			q.push(startNode);

			while (!q.empty()) {
				int u = q.front();
				q.pop();
				if (u < 0 || u >= static_cast<int>(adjList.size())) {
					continue;
				}

				for (int v : adjList[u]) {
					if (v < 0 || v >= nodeCount) continue;
					if (dist[v] >= 0) continue;
					dist[v] = dist[u] + 1;
					q.push(v);
				}
			}

			return dist;
		};

		bool anchorsOK = true;
		std::vector<int> distFromAnchor[3];
		for (int k = 0; k < 3; ++k) {
			if (ports[k].anchorNode < 0 || ports[k].anchorNode >= nodeCount) {
				anchorsOK = false;
				break;
			}
			distFromAnchor[k] = GraphDistances(ports[k].anchorNode);
		}

		if (anchorsOK) {
			double geomCenter[3] = { 0.0, 0.0, 0.0 };
			for (int k = 0; k < 3; ++k) {
				for (int d = 0; d < 3; ++d) {
					geomCenter[d] += ports[k].center[d] / 3.0;
				}
			}

			double bestScore = std::numeric_limits<double>::max();
			for (int node = 0; node < nodeCount; ++node) {
				if (node >= static_cast<int>(adjList.size())) continue;

				bool reachable = true;
				int sumDist = 0;
				for (int k = 0; k < 3; ++k) {
					if (distFromAnchor[k][node] < 0) {
						reachable = false;
						break;
					}
					sumDist += distFromAnchor[k][node];
				}
				if (!reachable) continue;

				double np[3] = {
					m_skelNodes[3 * node],
					m_skelNodes[3 * node + 1],
					m_skelNodes[3 * node + 2]
				};

				double degreePenalty = (adjList[node].size() >= 3) ? 0.0 : 1000000.0;
				double score = degreePenalty + static_cast<double>(sumDist) * 1000.0
					+ std::sqrt(vtkMath::Distance2BetweenPoints(np, geomCenter));

				if (score < bestScore) {
					bestScore = score;
					junctionNode = node;
				}
			}
		}

		if (junctionNode >= 0) {
			J[0] = m_skelNodes[3 * junctionNode];
			J[1] = m_skelNodes[3 * junctionNode + 1];
			J[2] = m_skelNodes[3 * junctionNode + 2];
			hasJunction = true;
		}
	}

	if (!hasJunction) {
		for (int k = 0; k < 3; ++k) {
			for (int d = 0; d < 3; ++d) {
				J[d] += ports[k].center[d] / 3.0;
			}
		}
	}

	// ------------------------------------------------------------
	// 3. 从 J 指向端口，判断 A/B0/B1。
	// 若有骨架交点，则优先用 “J -> 分支第一个相邻骨架节点” 的真实方向。
	// ------------------------------------------------------------
	for (int k = 0; k < 3; ++k) {
		for (int d = 0; d < 3; ++d) {
			ports[k].branchDir[d] = ports[k].center[d] - J[d];
		}
		NormalizeSafeLocal(ports[k].branchDir, ports[k].normal);
	}

	if (hasJunction && junctionNode >= 0 && !adjList.empty()) {
		const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);
		for (int k = 0; k < 3; ++k) {
			if (ports[k].anchorNode < 0 || ports[k].anchorNode >= nodeCount) continue;

			std::vector<int> path = FindPathBFS(junctionNode, ports[k].anchorNode, adjList);
			int firstNode = -1;
			if (path.size() >= 2) {
				firstNode = path[1];
			}
			else if (ports[k].anchorNode == junctionNode && junctionNode < static_cast<int>(adjList.size())) {
				double toCap[3] = {
					ports[k].center[0] - J[0],
					ports[k].center[1] - J[1],
					ports[k].center[2] - J[2]
				};
				NormalizeSafeLocal(toCap, ports[k].normal);

				double bestDot = -1e100;
				for (int nb : adjList[junctionNode]) {
					if (nb < 0 || nb >= nodeCount) continue;
					double v[3] = {
						m_skelNodes[3 * nb] - J[0],
						m_skelNodes[3 * nb + 1] - J[1],
						m_skelNodes[3 * nb + 2] - J[2]
					};
					NormalizeSafeLocal(v, toCap);
					double dot = vtkMath::Dot(v, toCap);
					if (dot > bestDot) {
						bestDot = dot;
						firstNode = nb;
					}
				}
			}

			if (firstNode >= 0 && firstNode < nodeCount) {
				double v[3] = {
					m_skelNodes[3 * firstNode] - J[0],
					m_skelNodes[3 * firstNode + 1] - J[1],
					m_skelNodes[3 * firstNode + 2] - J[2]
				};
				NormalizeSafeLocal(v, ports[k].branchDir);

				double toCap[3] = {
					ports[k].center[0] - J[0],
					ports[k].center[1] - J[1],
					ports[k].center[2] - J[2]
				};
				NormalizeSafeLocal(toCap, ports[k].normal);

				if (vtkMath::Dot(v, toCap) < 0.0) {
					for (int d = 0; d < 3; ++d) {
						v[d] = -v[d];
					}
				}

				for (int d = 0; d < 3; ++d) {
					ports[k].branchDir[d] = v[d];
				}
			}
		}
	}

	auto AngleBetweenPorts = [&](int i, int j) -> double
	{
		double dot = vtkMath::Dot(ports[i].branchDir, ports[j].branchDir);
		dot = ClampLocal(dot, -1.0, 1.0);
		return std::acos(dot);
	};

	const double a01 = AngleBetweenPorts(0, 1);
	const double a12 = AngleBetweenPorts(1, 2);
	const double a20 = AngleBetweenPorts(2, 0);

	int b0 = 0;
	int b1 = 1;
	if (a12 < a01 && a12 <= a20) {
		b0 = 1;
		b1 = 2;
	}
	else if (a20 < a01 && a20 < a12) {
		b0 = 2;
		b1 = 0;
	}

	int A = -1;
	for (int k = 0; k < 3; ++k) {
		if (k != b0 && k != b1) {
			A = k;
			break;
		}
	}
	if (A < 0) {
		return false;
	}

	// ------------------------------------------------------------
	// 4. 建立局部二维标架。
	// x：J -> A；y：局部平面内与 x 垂直；up：三分支所在平面法向。
	// ------------------------------------------------------------
	double up[3];
	vtkMath::Cross(ports[b0].branchDir, ports[b1].branchDir, up);
	if (vtkMath::Norm(up) < 1e-8) {
		vtkMath::Cross(ports[A].branchDir, ports[b0].branchDir, up);
	}
	NormalizeSafeLocal(up, ports[A].normal);

	double x[3] = { ports[A].branchDir[0], ports[A].branchDir[1], ports[A].branchDir[2] };
	const double xUp = vtkMath::Dot(x, up);
	for (int d = 0; d < 3; ++d) {
		x[d] -= xUp * up[d];
	}
	NormalizeSafeLocal(x, ports[A].branchDir);

	double y[3];
	vtkMath::Cross(up, x, y);
	NormalizeSafeLocal(y, ports[b0].branchDir);

	// 让 B0 位于 y 正侧，B1 位于 y 负侧，保证 tip-baseL / tip-baseR 分配稳定。
	double toB0[3] = {
		ports[b0].center[0] - J[0],
		ports[b0].center[1] - J[1],
		ports[b0].center[2] - J[2]
	};
	if (vtkMath::Dot(toB0, y) < 0.0) {
		for (int d = 0; d < 3; ++d) {
			y[d] = -y[d];
		}
	}

	// ------------------------------------------------------------
	// 5. 生成局部三角形。
	// 关键：A 只给底边方向，不复制 A 面、不射线贴外轮廓。
	// ------------------------------------------------------------
	const double avgBWidth = 0.5 * (ports[b0].width + ports[b1].width);
	double localW = std::min(ports[A].width, avgBWidth);
	localW = std::max(localW, 1e-3);

	double localDepth = (ports[A].depth + ports[b0].depth + ports[b1].depth) / 3.0;
	localDepth = ClampLocal(localDepth, localW * 0.20, localW * 1.20);

	const double distA = std::max(DistanceLocal(J, ports[A].center), localW);
	const double distB0 = std::max(DistanceLocal(J, ports[b0].center), localW);
	const double distB1 = std::max(DistanceLocal(J, ports[b1].center), localW);
	const double distB = 0.5 * (distB0 + distB1);

	// A 侧底边离 J 的小偏移。这里是局部参数，不再等于 A 面位置。
	const double baseOffset = ClampLocal(0.22 * distA, 0.25 * localW, 0.75 * localW);

	// 尖端方向由 B0/B1 的合方向决定，并强制放在 A 的反方向一侧。
	double tipDir[3] = {
		ports[b0].branchDir[0] + ports[b1].branchDir[0],
		ports[b0].branchDir[1] + ports[b1].branchDir[1],
		ports[b0].branchDir[2] + ports[b1].branchDir[2]
	};
	const double tipUp = vtkMath::Dot(tipDir, up);
	for (int d = 0; d < 3; ++d) {
		tipDir[d] -= tipUp * up[d];
	}
	NormalizeSafeLocal(tipDir, ports[b0].branchDir);
	if (vtkMath::Dot(tipDir, x) > 0.0) {
		for (int d = 0; d < 3; ++d) {
			tipDir[d] = -tipDir[d];
		}
	}

	const double tipLen = ClampLocal(0.28 * distB, 0.45 * localW, 1.20 * localW);
	const double baseHalf = ClampLocal(0.50 * localW, 0.35 * localW, 0.85 * localW);

	double baseMid[3];
	double baseL[3];
	double baseR[3];
	double tip[3];

	for (int d = 0; d < 3; ++d) {
		baseMid[d] = J[d] + x[d] * baseOffset;
		baseL[d] = baseMid[d] + y[d] * baseHalf;
		baseR[d] = baseMid[d] - y[d] * baseHalf;
		tip[d] = J[d] + tipDir[d] * tipLen;
	}

	// ------------------------------------------------------------
	// 6. 生成三角 hub poly：6 点 + 上下两个三角 cap。
	// 点序：0 tipTop, 1 tipBot, 2 baseLTop, 3 baseLBot, 4 baseRTop, 5 baseRBot。
	// ------------------------------------------------------------
	const double halfDepth = 0.5 * localDepth;

	double triPts[6][3];
	auto MakeTopLocal = [&](const double p[3], double out[3])
	{
		for (int d = 0; d < 3; ++d) {
			out[d] = p[d] + up[d] * halfDepth;
		}
	};
	auto MakeBotLocal = [&](const double p[3], double out[3])
	{
		for (int d = 0; d < 3; ++d) {
			out[d] = p[d] - up[d] * halfDepth;
		}
	};

	MakeTopLocal(tip, triPts[0]);
	MakeBotLocal(tip, triPts[1]);
	MakeTopLocal(baseL, triPts[2]);
	MakeBotLocal(baseL, triPts[3]);
	MakeTopLocal(baseR, triPts[4]);
	MakeBotLocal(baseR, triPts[5]);

	vtkSmartPointer<vtkPolyData> hubPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

	for (int i = 0; i < 6; ++i) {
		pts->InsertNextPoint(triPts[i]);
	}

	// 上表面。
	polys->InsertNextCell(3);
	polys->InsertCellPoint(0);
	polys->InsertCellPoint(4);
	polys->InsertCellPoint(2);

	// 下表面，反向。
	polys->InsertNextCell(3);
	polys->InsertCellPoint(1);
	polys->InsertCellPoint(3);
	polys->InsertCellPoint(5);

	hubPoly->SetPoints(pts);
	hubPoly->SetPolys(polys);
	m_hubCages.push_back(hubPoly);

	// ------------------------------------------------------------
	// 7. 三条三角边分别生成 bridge cage，桥接到 A/B0/B1 输入端面。
	// ------------------------------------------------------------
	auto MakeSideRingLocal = [&](int top0, int top1, int bot1, int bot0, const double outwardNormal[3], double ring[4][3])
	{
		for (int d = 0; d < 3; ++d) {
			ring[0][d] = triPts[top0][d];
			ring[1][d] = triPts[top1][d];
			ring[2][d] = triPts[bot1][d];
			ring[3][d] = triPts[bot0][d];
		}

		if (SignedAreaLocal(ring, outwardNormal) < 0.0) {
			ReverseRingLocal(ring);
		}
	};

	auto AddBridgeLocal = [&](int portIndex, const double hubRing[4][3], const double outwardNormal[3])
	{
		double bestTubeRing[4][3];
		BridgeScoreLocal bestScore;
		bestScore.minEdgeDot = -1e100;
		bestScore.avgEdgeDot = -1e100;
		bestScore.totalLength = 1e100;
		bestScore.maxDiagDiff = 1e100;

		for (int family = 0; family < 2; ++family) {
			const bool reverse = (family == 1);
			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				BuildTubeCandidateRingLocal(ports[portIndex].ring, reverse, shift, candidate);

				BridgeScoreLocal score = EvaluateBridgeLocal(hubRing, candidate);
				if (IsBetterBridgeScoreLocal(score, bestScore)) {
					bestScore = score;
					CopyRingLocal(bestTubeRing, candidate);
				}
			}
		}

		// 如果桥接端面方向明显反了，则整体反向一次作为兜底。
		double areaHub = SignedAreaLocal(hubRing, outwardNormal);
		double areaEnd = SignedAreaLocal(bestTubeRing, outwardNormal);
		if (areaHub * areaEnd < 0.0) {
			ReverseRingLocal(bestTubeRing);
		}

		CageSegment bridge;
		CopyRingLocal(bridge.pStart, hubRing);
		CopyRingLocal(bridge.pEnd, bestTubeRing);
		m_allCages.push_back(bridge);

		if (false) qDebug() << "[TriHub2D][Bridge]"
			<< "port =" << portIndex
			<< "nearest0 =" << FindNearestIndexLocal(bestTubeRing[0], ports[portIndex].ring)
			<< "minEdgeDot =" << bestScore.minEdgeDot
			<< "avgEdgeDot =" << bestScore.avgEdgeDot
			<< "totalLength =" << bestScore.totalLength;
	};

	double ringA[4][3];
	double ringB0[4][3];
	double ringB1[4][3];

	// A 连接底边 baseL-baseR。
	MakeSideRingLocal(2, 4, 5, 3, ports[A].branchDir, ringA);

	// B0 位于 y 正侧，连接 tip-baseL。
	MakeSideRingLocal(0, 2, 3, 1, ports[b0].branchDir, ringB0);

	// B1 位于 y 负侧，连接 baseR-tip。
	MakeSideRingLocal(4, 0, 1, 5, ports[b1].branchDir, ringB1);

	AddBridgeLocal(A, ringA, ports[A].branchDir);
	AddBridgeLocal(b0, ringB0, ports[b0].branchDir);
	AddBridgeLocal(b1, ringB1, ports[b1].branchDir);

	// ------------------------------------------------------------
	// 8. SDF 信息保持原流程兼容。
	// ------------------------------------------------------------
	HubSDFInfo sdfInfo;
	for (int d = 0; d < 3; ++d) {
		sdfInfo.hubCenter[d] = J[d];
	}
	const double coreRadius = std::max(0.25 * localW, 1e-3);
	sdfInfo.coreRadius = coreRadius;
	sdfInfo.maxRadius = std::max(tipLen + baseOffset + baseHalf, localW);

	const int portOrder[3] = { A, b0, b1 };
	for (int ii = 0; ii < 3; ++ii) {
		int k = portOrder[ii];
		SDFTube tube;
		for (int d = 0; d < 3; ++d) {
			tube.center[d] = J[d];
			tube.axis[d] = ports[k].branchDir[d];
		}
		NormalizeSafeLocal(tube.axis, ports[k].normal);
		tube.radius = std::max(std::min(ports[k].width, ports[k].depth) * 0.5, coreRadius);
		sdfInfo.tubes.push_back(tube);
		sdfInfo.maxRadius = std::max(sdfInfo.maxRadius,
			DistanceLocal(J, ports[k].center) + std::max(ports[k].width, ports[k].depth));
	}
	m_sdfHubs.push_back(sdfInfo);

	if (false) qDebug() << "[TriHub2D] generated local triangle hub"
		<< "A =" << A
		<< "B0 =" << b0
		<< "B1 =" << b1
		<< "junctionNode =" << junctionNode
		<< "J =" << J[0] << J[1] << J[2]
		<< "baseOffset =" << baseOffset
		<< "baseHalf =" << baseHalf
		<< "tipLen =" << tipLen
		<< "localW =" << localW
		<< "localDepth =" << localDepth;

	return true;
}

void CageInteractorStyle::GenerateNBranch()
{
	m_use3DHubSdfProjection = false;

	int N = static_cast<int>(m_selectedFacesForBranch.size());
	if (N < 3) return;

	if (N == 3) {
		if (GenerateNBranchTriHub2D()) {
			RebuildGlobalMesh();
			return;
		}
		if (false) qDebug() << "[TriHub2D] failed, fallback to old GenerateNBranch().";
	}

	if (false) qDebug() << "\n==================================================";
	if (false) qDebug() << "[复杂分支生成] 开始，端面数量 N =" << N;

	// ============================================================
	// 0. 
	// ============================================================

	auto CopyRing = [](double dst[4][3], const double src[4][3])
	{
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				dst[i][d] = src[i][d];
			}
		}
	};

	auto RingCenter = [](const double ring[4][3], double c[3])
	{
		c[0] = 0.0;
		c[1] = 0.0;
		c[2] = 0.0;

		for (int i = 0; i < 4; ++i) {
			c[0] += ring[i][0] * 0.25;
			c[1] += ring[i][1] * 0.25;
			c[2] += ring[i][2] * 0.25;
		}
	};

	auto Distance = [](const double a[3], const double b[3]) -> double
	{
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto DotEdgeNormalized = [](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double
	{
		double va[3] = {
			a1[0] - a0[0],
			a1[1] - a0[1],
			a1[2] - a0[2]
		};

		double vb[3] = {
			b1[0] - b0[0],
			b1[1] - b0[1],
			b1[2] - b0[2]
		};

		double na = vtkMath::Norm(va);
		double nb = vtkMath::Norm(vb);

		if (na < 1e-10 || nb < 1e-10) {
			return 0.0;
		}

		va[0] /= na;
		va[1] /= na;
		va[2] /= na;

		vb[0] /= nb;
		vb[1] /= nb;
		vb[2] /= nb;

		return vtkMath::Dot(va, vb);
	};

	auto SignedAreaOnRing = [](const double ring[4][3], const double normal[3]) -> double
	{
		double area[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double c[3];
			vtkMath::Cross(ring[i], ring[j], c);

			area[0] += c[0];
			area[1] += c[1];
			area[2] += c[2];
		}

		return vtkMath::Dot(area, normal);
	};

	auto FindNearestIndex = [](const double p[3], const double ring[4][3]) -> int
	{
		int bestId = -1;
		double bestD2 = 1e100;

		for (int i = 0; i < 4; ++i) {
			double d2 = vtkMath::Distance2BetweenPoints(p, ring[i]);

			if (d2 < bestD2) {
				bestD2 = d2;
				bestId = i;
			}
		}

		return bestId;
	};

	struct BridgeScore
	{
		double minEdgeDot;
		double avgEdgeDot;
		double totalLength;
		double maxDiagDiff;
	};

	auto EvaluateBridge = [&DotEdgeNormalized, &Distance](
		const double startRing[4][3],
		const double endRing[4][3]) -> BridgeScore
	{
		BridgeScore score;
		score.minEdgeDot = 1.0;
		score.avgEdgeDot = 0.0;
		score.totalLength = 0.0;
		score.maxDiagDiff = 0.0;

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double dot = DotEdgeNormalized(
				startRing[i],
				startRing[j],
				endRing[i],
				endRing[j]
			);

			score.minEdgeDot = std::min(score.minEdgeDot, dot);
			score.avgEdgeDot += dot * 0.25;

			score.totalLength += Distance(startRing[i], endRing[i]);

			double diagA = Distance(startRing[i], endRing[j]);
			double diagB = Distance(startRing[j], endRing[i]);

			score.maxDiagDiff = std::max(score.maxDiagDiff, std::abs(diagA - diagB));
		}

		return score;
	};

	auto IsBetterScore = [](const BridgeScore& a, const BridgeScore& b) -> bool
	{
		// 
		if (a.minEdgeDot > b.minEdgeDot + 1e-6) return true;
		if (a.minEdgeDot < b.minEdgeDot - 1e-6) return false;

		// 
		if (a.avgEdgeDot > b.avgEdgeDot + 1e-6) return true;
		if (a.avgEdgeDot < b.avgEdgeDot - 1e-6) return false;

		// 
		if (a.maxDiagDiff < b.maxDiagDiff - 1e-6) return true;
		if (a.maxDiagDiff > b.maxDiagDiff + 1e-6) return false;

		//
		return a.totalLength < b.totalLength;
	};

	// ============================================================
	// 
	// ============================================================

	auto BuildTubeCandidateRing = [](
		const double tubeCanonical[4][3],
		bool reverse,
		int shift,
		double outRing[4][3],
		int outOrder[4])
	{
		static const int forwardBase[4] = { 0, 1, 2, 3 };
		static const int reverseBase[4] = { 0, 3, 2, 1 };

		const int* base = reverse ? reverseBase : forwardBase;

		for (int i = 0; i < 4; ++i) {
			int id = base[(i + shift) % 4];
			outOrder[i] = id;

			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = tubeCanonical[id][d];
			}
		}
	};

	auto DebugChosenBridge = [&](
		int sortedIndex,
		int originalTubeIndex,
		bool tubeIsStart,
		const double hubRing[4][3],
		const double tubeCanonical[4][3],
		const double bestTubeRing[4][3],
		const int bestOrder[4],
		bool bestReverse,
		const BridgeScore& bestScore,
		const double normalUp[3])
	{
		if (false) qDebug() << "\n==================== [桥接调试] ====================";
		if (false) qDebug() << "排序后序号 =" << sortedIndex
			<< "原始管口序号 =" << originalTubeIndex
			<< "选中端是否起始端 =" << tubeIsStart;
		if (false) qDebug() << "选中的点序类型 =" << (bestReverse ? "反向/循环平移" : "正向/循环平移");
		if (false) qDebug() << "选中的点序 ="
			<< bestOrder[0]
			<< bestOrder[1]
			<< bestOrder[2]
			<< bestOrder[3];

		if (false) qDebug() << "上方向 normalUp =" << normalUp[0] << normalUp[1] << normalUp[2];

		if (false) qDebug() << "--- 接头侧环 / 桥接起点环 ---";
		if (false) qDebug() << "hubRing[0] = 上侧(前一个)";
		if (false) qDebug() << "hubRing[1] = 下侧(前一个)";
		if (false) qDebug() << "hubRing[2] = 下侧(当前)";
		if (false) qDebug() << "hubRing[3] = 上侧(当前)";

		for (int i = 0; i < 4; ++i) {
			if (false) qDebug().nospace()
				<< "  hubRing[" << i << "] = ("
				<< hubRing[i][0] << ", "
				<< hubRing[i][1] << ", "
				<< hubRing[i][2] << ")";
		}

		if (false) qDebug() << "--- tubeCanonical ?---";

		for (int i = 0; i < 4; ++i) {
			if (false) qDebug().nospace()
				<< "  tubeCanonical[" << i << "] = ("
				<< tubeCanonical[i][0] << ", "
				<< tubeCanonical[i][1] << ", "
				<< tubeCanonical[i][2] << ")";
		}

		if (false) qDebug() << "--- bridge.pEnd ---";

		for (int i = 0; i < 4; ++i) {
			int nearId = FindNearestIndex(bestTubeRing[i], tubeCanonical);

			if (false) qDebug().nospace()
				<< "  bridge.pEnd[" << i << "] = tubeCanonical["
				<< nearId << "]  ("
				<< bestTubeRing[i][0] << ", "
				<< bestTubeRing[i][1] << ", "
				<< bestTubeRing[i][2] << ")";
		}

		if (false) qDebug() << "--- bridge ?hubRing[i] -> bridge.pEnd[i] ---";

		for (int i = 0; i < 4; ++i) {
			int nearId = FindNearestIndex(bestTubeRing[i], tubeCanonical);
			double len = Distance(hubRing[i], bestTubeRing[i]);

			if (false) qDebug().nospace()
				<< "  i=" << i
				<< " : hub[" << i << "] -> tube[" << nearId << "]"
				<< " length=" << len;
		}

		if (false) qDebug() << "--- ?---";

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double dot = DotEdgeNormalized(
				hubRing[i],
				hubRing[j],
				bestTubeRing[i],
				bestTubeRing[j]
			);

			double diagA = Distance(hubRing[i], bestTubeRing[j]);
			double diagB = Distance(hubRing[j], bestTubeRing[i]);

			if (false) qDebug().nospace()
				<< "  F" << i
				<< " edgeDot=" << dot
				<< " diagA=" << diagA
				<< " diagB=" << diagB;

			if (dot < 0.20) {
				if (false) qDebug() << "    !!! WARNING: ";
			}
		}

		double areaTube = SignedAreaOnRing(tubeCanonical, normalUp);
		double areaEnd = SignedAreaOnRing(bestTubeRing, normalUp);
		double areaHub = SignedAreaOnRing(hubRing, normalUp);

		if (false) qDebug() << "--- ---";
		if (false) qDebug() << "管口标准环面积 =" << areaTube;
		if (false) qDebug() << "最佳管口环面积 =" << areaEnd;
		if (false) qDebug() << "接头侧环面积 =" << areaHub;
		if (false) qDebug() << "最佳管口环与原管口环方向是否一致？"
			<< ((areaTube * areaEnd >= 0.0) ? "YES" : "NO");

		if (false) qDebug() << "--- 鏈€缁堣瘎鍒?---";
		if (false) qDebug() << "最小边方向点积 minEdgeDot =" << bestScore.minEdgeDot;
		if (false) qDebug() << "平均边方向点积 avgEdgeDot =" << bestScore.avgEdgeDot;
		if (false) qDebug() << "最大对角线差 maxDiagDiff =" << bestScore.maxDiagDiff;
		if (false) qDebug() << "总连接长度 totalLength =" << bestScore.totalLength;

		if (false) qDebug() << "==================== [桥接调试结束] ====================\n";
	};

	struct TubeData
	{
		int index;
		double center[3];

		// dir_out：hubCenter -> 端面中心的几何方向，用于原二维hub排序/摆放。
		double dir_out[3];

		// skeletonDir：当前选中端面所在控制笼的真实中心线方向，
		// 统一成“控制笼内部 -> 当前端面”的方向。它只作为备用方向。
		double skeletonDir[3];

		// branchDir：骨架交点 -> 当前分支相邻第一个骨架节点的真实骨架方向。
		// N=3 的 A/B 判断、主轴投影、整体避让移动都优先使用它。
		double branchDir[3];
		int anchorNode;

		double angle;
		double w;
		double h;
		double distToHub;
		double capRing[4][3];
		double cagePts8[8][3];
		bool isStart;
	};

	std::vector<TubeData> tubes(N);
	double geomCenter[3] = { 0.0, 0.0, 0.0 };

	// ============================================================
	// 1. 
	// ============================================================

	for (int k = 0; k < N; ++k) {
		auto sel = m_selectedFacesForBranch[k];
		if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
			if (false) qDebug() << "[复杂分支生成] 选中的控制笼序号无效 =" << sel.cageIndex;
			return;
		}

		tubes[k].index = k;
		tubes[k].isStart = sel.isStart;
		tubes[k].anchorNode = -1;

		for (int i = 0; i < 4; ++i) {
			const double* p = sel.isStart
				? m_allCages[sel.cageIndex].pStart[i]
				: m_allCages[sel.cageIndex].pEnd[i];

			for (int d = 0; d < 3; ++d) {
				tubes[k].capRing[i][d] = p[d];
			}
		}

		// N=3二维hub避让要用“B端面所在控制笼的8个节点”在主骨架方向上的投影范围，
		// 所以这里把该端面所属控制笼的起点/终点共8点都缓存下来。
		double startCenterForAxis[3] = { 0.0, 0.0, 0.0 };
		double endCenterForAxis[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				tubes[k].cagePts8[i][d] = m_allCages[sel.cageIndex].pStart[i][d];
				tubes[k].cagePts8[i + 4][d] = m_allCages[sel.cageIndex].pEnd[i][d];

				startCenterForAxis[d] += m_allCages[sel.cageIndex].pStart[i][d] * 0.25;
				endCenterForAxis[d] += m_allCages[sel.cageIndex].pEnd[i][d] * 0.25;
			}
		}

		// 真实骨架方向取“当前选中端面所在控制笼的中心线方向”，
		// 并统一成从控制笼内部指向当前端面的方向。
		// 选中 start 端面：End -> Start；选中 end 端面：Start -> End。
		double skelDir[3];
		if (sel.isStart) {
			skelDir[0] = startCenterForAxis[0] - endCenterForAxis[0];
			skelDir[1] = startCenterForAxis[1] - endCenterForAxis[1];
			skelDir[2] = startCenterForAxis[2] - endCenterForAxis[2];
		}
		else {
			skelDir[0] = endCenterForAxis[0] - startCenterForAxis[0];
			skelDir[1] = endCenterForAxis[1] - startCenterForAxis[1];
			skelDir[2] = endCenterForAxis[2] - startCenterForAxis[2];
		}

		if (vtkMath::Norm(skelDir) > 1e-8) {
			vtkMath::Normalize(skelDir);
		}
		else {
			skelDir[0] = 1.0;
			skelDir[1] = 0.0;
			skelDir[2] = 0.0;
		}

		for (int d = 0; d < 3; ++d) {
			tubes[k].skeletonDir[d] = skelDir[d];
			tubes[k].branchDir[d] = skelDir[d];
		}

		// 先绑定到最近骨架节点。真正用于 N=3 A/B 判断的 branchDir 会在 hubCenter
		// 初步确定后，由“骨架交点 -> 分支第一个相邻节点”重新计算；如果失败则保留 skeletonDir 兜底。
		if (!m_skelNodes.empty()) {
			double tmpCenter[3] = { 0.0, 0.0, 0.0 };
			double tmpNormal[3] = { 1.0, 0.0, 0.0 };
			tubes[k].anchorNode = FindNearestSkeletonNodeToSelectedFace(sel, tmpCenter, tmpNormal);
		}

		double* p0 = tubes[k].capRing[0];
		double* p1 = tubes[k].capRing[1];
		double* p2 = tubes[k].capRing[2];

		tubes[k].w = std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		tubes[k].h = std::sqrt(vtkMath::Distance2BetweenPoints(p1, p2));

		double pCenter[3] = { 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				pCenter[d] += tubes[k].capRing[i][d] * 0.25;
			}
		}

		for (int d = 0; d < 3; ++d) {
			tubes[k].center[d] = pCenter[d];
			geomCenter[d] += pCenter[d] / static_cast<double>(N);
		}
	}

	double hubCenter[3];
	std::memcpy(hubCenter, geomCenter, sizeof(double) * 3);

	// ============================================================
	// 2. 
	// ============================================================

	if (false) qDebug() << "--- 开始微调接头中心位置 ---";

	int maxIter = 50;

	for (int iter = 0; iter < maxIter; ++iter) {
		double shift[3] = { 0.0, 0.0, 0.0 };
		bool collision = false;

		for (int k = 0; k < N; ++k) {
			double vec[3] = {
				hubCenter[0] - tubes[k].center[0],
				hubCenter[1] - tubes[k].center[1],
				hubCenter[2] - tubes[k].center[2]
			};

			double dist = vtkMath::Norm(vec);
			double safeDist = tubes[k].w * 0.8;

			if (dist < safeDist && dist > 1e-5) {
				collision = true;

				double penetration = safeDist - dist;
				double step = std::min(penetration * 0.2, tubes[k].w * 0.05);

				for (int d = 0; d < 3; ++d) {
					shift[d] += (vec[d] / dist) * step;
				}
			}
		}

		if (!collision) break;

		double nextPos[3] = {
			hubCenter[0] + shift[0],
			hubCenter[1] + shift[1],
			hubCenter[2] + shift[2]
		};

		if (m_RefLocator) {
			double closestPoint[3];
			double dist2;
			vtkIdType cellId;
			int subId;

			m_RefLocator->FindClosestPoint(nextPos, closestPoint, cellId, subId, dist2);

			if (std::sqrt(dist2) > tubes[0].w * 0.5) {
				for (int d = 0; d < 3; ++d) {
					nextPos[d] = hubCenter[d];
				}
			}
		}

		std::memcpy(hubCenter, nextPos, sizeof(double) * 3);
	}

	// ============================================================
	// 3. 
	// ============================================================

	double minCenterDist = 1e18;

	for (int k = 0; k < N; ++k) {
		double geom_dir[3] = {
			tubes[k].center[0] - hubCenter[0],
			tubes[k].center[1] - hubCenter[1],
			tubes[k].center[2] - hubCenter[2]
		};

		tubes[k].distToHub = vtkMath::Norm(geom_dir);

		if (tubes[k].distToHub > 1e-6) {
			vtkMath::Normalize(geom_dir);
		}
		else {
			geom_dir[0] = 1.0;
			geom_dir[1] = 0.0;
			geom_dir[2] = 0.0;
		}

		for (int d = 0; d < 3; ++d) {
			tubes[k].dir_out[d] = geom_dir[d];
		}

		// skeletonDir 的长度来自真实控制笼中心线，但方向有时会因为 start/end 或输入顺序与
		// hubCenter->端面中心相反。这里只做符号统一，让它也指向当前端面一侧；
		// 后续夹角和 mainAxis 仍然使用真实控制笼中心线，而不是使用 geom_dir。
		if (vtkMath::Dot(tubes[k].skeletonDir, tubes[k].dir_out) < 0.0) {
			for (int d = 0; d < 3; ++d) {
				tubes[k].skeletonDir[d] = -tubes[k].skeletonDir[d];
			}
		}

		for (int d = 0; d < 3; ++d) {
			tubes[k].branchDir[d] = tubes[k].skeletonDir[d];
		}

		minCenterDist = std::min(minCenterDist, tubes[k].distToHub);
	}

	// N=3 的 hub 中心最终必须锁在真实骨架交点；后续几何约束不能再移动它。
	double n3LockedHubCenter[3] = { hubCenter[0], hubCenter[1], hubCenter[2] };
	bool n3HasLockedSkeletonJunction = false;

	// N=3 的 A/B 判定不能使用端面控制笼中线，也不能使用 hubCenter->端面中心。
	// 这里改为使用真实骨架图：先找三条端面 anchor 共同汇入的骨架交点，
	// 再取“交点 -> 各分支相邻第一个骨架节点”的三条向量作为夹角判断轴。
	if (N == 3 && !m_skelNodes.empty() && !m_skelEdges.empty()) {
		std::vector<std::vector<int>> adjList = BuildSkeletonAdjacency();
		const int nodeCount = static_cast<int>(m_skelNodes.size() / 3);

		auto GetNodePoint = [&](int node, double p[3]) -> bool
		{
			if (node < 0 || node >= nodeCount) return false;
			p[0] = m_skelNodes[3 * node];
			p[1] = m_skelNodes[3 * node + 1];
			p[2] = m_skelNodes[3 * node + 2];
			return true;
		};

		auto GraphDistances = [&](int startNode) -> std::vector<int>
		{
			std::vector<int> dist(nodeCount, -1);
			if (startNode < 0 || startNode >= nodeCount) return dist;

			std::queue<int> q;
			dist[startNode] = 0;
			q.push(startNode);

			while (!q.empty()) {
				int u = q.front();
				q.pop();
				if (u < 0 || u >= static_cast<int>(adjList.size())) continue;
				for (int v : adjList[u]) {
					if (v < 0 || v >= nodeCount) continue;
					if (dist[v] >= 0) continue;
					dist[v] = dist[u] + 1;
					q.push(v);
				}
			}
			return dist;
		};

		bool anchorsOK = (static_cast<int>(adjList.size()) == nodeCount);
		std::vector<std::vector<int>> distFromAnchor(3);
		for (int k = 0; k < 3; ++k) {
			if (tubes[k].anchorNode < 0 || tubes[k].anchorNode >= nodeCount) anchorsOK = false;
			else distFromAnchor[k] = GraphDistances(tubes[k].anchorNode);
		}

		int junctionNode = -1;
		double bestScore = std::numeric_limits<double>::max();
		if (anchorsOK) {
			for (int node = 0; node < nodeCount; ++node) {
				if (node >= static_cast<int>(adjList.size())) continue;

				bool reachable = true;
				int sumDist = 0;
				for (int k = 0; k < 3; ++k) {
					if (distFromAnchor[k][node] < 0) {
						reachable = false;
						break;
					}
					sumDist += distFromAnchor[k][node];
				}
				if (!reachable) continue;

				double jp[3];
				if (!GetNodePoint(node, jp)) continue;
				double hubD = std::sqrt(vtkMath::Distance2BetweenPoints(jp, hubCenter));

				// 强优先选择真实分叉节点；没有度数>=3节点时才允许普通节点兜底。
				double degreePenalty = (adjList[node].size() >= 3) ? 0.0 : 1000000.0;
				double score = degreePenalty + static_cast<double>(sumDist) * 1000.0 + hubD;
				if (score < bestScore) {
					bestScore = score;
					junctionNode = node;
				}
			}
		}

		if (junctionNode >= 0) {
			double jp[3];
			GetNodePoint(junctionNode, jp);
			for (int d = 0; d < 3; ++d) {
				n3LockedHubCenter[d] = jp[d];
				hubCenter[d] = jp[d];
			}
			n3HasLockedSkeletonJunction = true;

			for (int k = 0; k < 3; ++k) {
				std::vector<int> path = FindPathBFS(junctionNode, tubes[k].anchorNode, adjList);

				int firstNode = -1;
				if (path.size() >= 2) {
					firstNode = path[1];
				}
				else if (tubes[k].anchorNode == junctionNode && junctionNode < static_cast<int>(adjList.size())) {
					// 极少数情况下端面 anchor 就落在交点上，此时选最朝向端面中心的邻居。
					double toCap[3] = { tubes[k].center[0] - jp[0], tubes[k].center[1] - jp[1], tubes[k].center[2] - jp[2] };
					if (vtkMath::Norm(toCap) > 1e-8) vtkMath::Normalize(toCap);
					double bestDot = -1e100;
					for (int nb : adjList[junctionNode]) {
						double np[3];
						if (!GetNodePoint(nb, np)) continue;
						double v[3] = { np[0] - jp[0], np[1] - jp[1], np[2] - jp[2] };
						if (vtkMath::Norm(v) <= 1e-8) continue;
						vtkMath::Normalize(v);
						double dot = vtkMath::Dot(v, toCap);
						if (dot > bestDot) {
							bestDot = dot;
							firstNode = nb;
						}
					}
				}

				if (firstNode >= 0) {
					double fp[3];
					GetNodePoint(firstNode, fp);
					double branchDir[3] = { fp[0] - jp[0], fp[1] - jp[1], fp[2] - jp[2] };
					if (vtkMath::Norm(branchDir) > 1e-8) {
						vtkMath::Normalize(branchDir);
						double toCap[3] = { tubes[k].center[0] - jp[0], tubes[k].center[1] - jp[1], tubes[k].center[2] - jp[2] };
						if (vtkMath::Norm(toCap) > 1e-8) {
							vtkMath::Normalize(toCap);
							if (vtkMath::Dot(branchDir, toCap) < 0.0) {
								branchDir[0] = -branchDir[0];
								branchDir[1] = -branchDir[1];
								branchDir[2] = -branchDir[2];
							}
						}

						for (int d = 0; d < 3; ++d) {
							tubes[k].branchDir[d] = branchDir[d];
						}
					}
				}
			}

			if (false) qDebug() << "[2DHubSafety][N3BranchDir] use skeleton junction"
				<< "junctionNode =" << junctionNode
				<< "anchorNodes =" << tubes[0].anchorNode << tubes[1].anchorNode << tubes[2].anchorNode;
		}
		else {
			if (false) qDebug() << "[2DHubSafety][N3BranchDir] cannot find common skeleton junction; fallback to cage centerline direction.";
		}
	}

	// ============================================================
	// 4.
	// ============================================================

	double normalUp[3] = { 0.0, 0.0, 0.0 };
	double bestCrossNorm = 0.0;

	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			double c[3];
			vtkMath::Cross(tubes[i].dir_out, tubes[j].dir_out, c);

			double n = vtkMath::Norm(c);

			if (n > bestCrossNorm) {
				bestCrossNorm = n;
				normalUp[0] = c[0];
				normalUp[1] = c[1];
				normalUp[2] = c[2];
			}
		}
	}

	if (bestCrossNorm < 1e-8) {
		normalUp[0] = 0.0;
		normalUp[1] = 0.0;
		normalUp[2] = 1.0;
	}
	else {
		vtkMath::Normalize(normalUp);
	}

	double refAxis[3];
	double refSide[3];

	for (int d = 0; d < 3; ++d) {
		refAxis[d] = tubes[0].dir_out[d];
	}

	vtkMath::Cross(normalUp, refAxis, refSide);

	if (vtkMath::Norm(refSide) < 1e-8) {
		refSide[0] = 1.0;
		refSide[1] = 0.0;
		refSide[2] = 0.0;
	}
	else {
		vtkMath::Normalize(refSide);
	}

	vtkMath::Cross(refSide, normalUp, refAxis);
	vtkMath::Normalize(refAxis);

	for (int k = 0; k < N; ++k) {
		tubes[k].angle = std::atan2(
			vtkMath::Dot(tubes[k].dir_out, refSide),
			vtkMath::Dot(tubes[k].dir_out, refAxis)
		);
	}

	std::sort(tubes.begin(), tubes.end(), [](const TubeData& a, const TubeData& b) {
		return a.angle < b.angle;
	});

	// ============================================================
	// 5. 
	// ============================================================

	vtkSmartPointer<vtkPolyData> hubPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

	for (int k = 0; k < N; ++k) {
		int next_k = (k + 1) % N;

		double a1 = tubes[k].angle;
		double a2 = tubes[next_k].angle;

		if (a2 <= a1) {
			a2 += 2.0 * vtkMath::Pi();
		}

		double midAngle = (a1 + a2) * 0.5;

		double local_W = (tubes[k].w + tubes[next_k].w) * 0.5;
		double local_H = (tubes[k].h + tubes[next_k].h) * 0.5;

		double local_R = std::min(local_W * 0.8, minCenterDist * 0.9);

		double px = local_R * std::cos(midAngle);
		double py = local_R * std::sin(midAngle);

		double pt[3] = {
			hubCenter[0] + px * refAxis[0] + py * refSide[0],
			hubCenter[1] + px * refAxis[1] + py * refSide[1],
			hubCenter[2] + px * refAxis[2] + py * refSide[2]
		};

		// 
		pts->InsertNextPoint(
			pt[0] + normalUp[0] * (local_H * 0.5),
			pt[1] + normalUp[1] * (local_H * 0.5),
			pt[2] + normalUp[2] * (local_H * 0.5)
		);

		// 
		pts->InsertNextPoint(
			pt[0] - normalUp[0] * (local_H * 0.5),
			pt[1] - normalUp[1] * (local_H * 0.5),
			pt[2] - normalUp[2] * (local_H * 0.5)
		);
	}


	// ============================================================
	// 5.5 2D中心hub安全修正：尺寸匹配 + 端面越界避让
	//
	// 目的：
	// 1) 不改变二维中心hub的基本生成方式；
	// 2) 在背景空间允许时，让每个hub连接面尽量接近输入端面尺寸；
	// 3) 如果hub节点越过或过近输入端面，优先整体朝反方向避让；
	// 4) 如果整体避让会导致其它端面越界，再退化为局部单节点避让。
	// ============================================================
	{
		auto GetHubPoint = [&](int id, double p[3])
		{
			pts->GetPoint(id, p);
		};

		auto SetHubPoint = [&](int id, const double p[3])
		{
			pts->SetPoint(id, p);
		};

		auto GetHubRingForTube = [&](int k, double ring[4][3])
		{
			int prev_k = (k - 1 + N) % N;
			pts->GetPoint(2 * prev_k, ring[0]);
			pts->GetPoint(2 * prev_k + 1, ring[1]);
			pts->GetPoint(2 * k + 1, ring[2]);
			pts->GetPoint(2 * k, ring[3]);
		};

		auto PutHubRingForTube = [&](int k, const double ring[4][3])
		{
			int prev_k = (k - 1 + N) % N;
			pts->SetPoint(2 * prev_k, ring[0]);
			pts->SetPoint(2 * prev_k + 1, ring[1]);
			pts->SetPoint(2 * k + 1, ring[2]);
			pts->SetPoint(2 * k, ring[3]);
		};

		auto PointPlaneT = [&](const double p[3], int k) -> double
		{
			// 保留旧的中心平面诊断值：dir_out 指向输入端面；t >= 0 表示越过端面中心平面。
			double v[3] = {
				p[0] - tubes[k].center[0],
				p[1] - tubes[k].center[1],
				p[2] - tubes[k].center[2]
			};
			return vtkMath::Dot(v, tubes[k].dir_out);
		};

		auto SCoordForTube = [&](const double p[3], int k) -> double
		{
			return vtkMath::Dot(p, tubes[k].dir_out);
		};

		auto CapEntrySForTube = [&](int k) -> double
		{
			// 三维避让同款思想：不用端面中心平面做软判断，
			// 而是使用输入端面四个原始角点在 dir_out 方向上的最近入口。
			// hub 必须留在该入口之前，否则视觉上就会穿过端面。
			double entryS = std::numeric_limits<double>::max();
			for (int i = 0; i < 4; ++i) {
				entryS = std::min(entryS, SCoordForTube(tubes[k].capRing[i], k));
			}
			return entryS;
		};

		auto SafeGapForTube = [&](int k) -> double
		{
			double s = std::min(tubes[k].w, tubes[k].h);
			return std::max(s * 0.10, 1e-4);
		};

		auto CheckAllHubPointsSafe = [&](double* worstViolationOut, int* worstTubeOut) -> bool
		{
			double worst = -1e100;
			int worstTube = -1;

			for (int k = 0; k < N; ++k) {
				double gap = SafeGapForTube(k);
				double limitS = CapEntrySForTube(k) - gap;
				double ring[4][3];
				GetHubRingForTube(k, ring);

				for (int i = 0; i < 4; ++i) {
					// 精确入口限制：S(hub) 必须 <= minS(input cap vertices) - safeGap。
					double violation = SCoordForTube(ring[i], k) - limitS;
					if (violation > worst) {
						worst = violation;
						worstTube = k;
					}
				}
			}

			if (worstViolationOut) *worstViolationOut = worst;
			if (worstTubeOut) *worstTubeOut = worstTube;
			return worst <= 1e-7;
		};

		auto SnapshotPoints = [&]() -> std::vector<std::array<double, 3>>
		{
			std::vector<std::array<double, 3>> snap;
			snap.resize(static_cast<size_t>(pts->GetNumberOfPoints()));
			for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i) {
				double p[3];
				pts->GetPoint(i, p);
				snap[static_cast<size_t>(i)] = { p[0], p[1], p[2] };
			}
			return snap;
		};

		auto RestorePoints = [&](const std::vector<std::array<double, 3>>& snap)
		{
			for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i) {
				const auto& p = snap[static_cast<size_t>(i)];
				pts->SetPoint(i, p.data());
			}
		};

		auto RefSpaceAllowed = [&](double marginScale) -> bool
		{
			if (!m_RefLocator) return true;

			double refTol = 0.0;
			for (int k = 0; k < N; ++k) {
				refTol += std::min(tubes[k].w, tubes[k].h) / static_cast<double>(N);
			}
			refTol = std::max(refTol * marginScale, 1e-4);

			for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i) {
				double p[3];
				pts->GetPoint(i, p);

				double closestPoint[3];
				double dist2 = 0.0;
				vtkIdType cellId = -1;
				int subId = 0;
				m_RefLocator->FindClosestPoint(p, closestPoint, cellId, subId, dist2);

				if (std::sqrt(dist2) > refTol) {
					return false;
				}
			}
			return true;
		};

		auto MoveAllHubPoints = [&](const double delta[3])
		{
			for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i) {
				double p[3];
				pts->GetPoint(i, p);
				p[0] += delta[0];
				p[1] += delta[1];
				p[2] += delta[2];
				pts->SetPoint(i, p);
			}
		};

		auto MovePointRadiallyInPlane = [&](int id, double amount)
		{
			double p[3];
			pts->GetPoint(id, p);

			double v[3] = {
				p[0] - hubCenter[0],
				p[1] - hubCenter[1],
				p[2] - hubCenter[2]
			};

			// 只取二维hub平面内的径向分量，避免改变上下厚度方向。
			double up = vtkMath::Dot(v, normalUp);
			v[0] -= normalUp[0] * up;
			v[1] -= normalUp[1] * up;
			v[2] -= normalUp[2] * up;

			double len = vtkMath::Norm(v);
			if (len < 1e-8) return;

			v[0] /= len;
			v[1] /= len;
			v[2] /= len;

			p[0] += v[0] * amount;
			p[1] += v[1] * amount;
			p[2] += v[2] * amount;
			pts->SetPoint(id, p);
		};

		if (false) qDebug() << "[2DHubSafety] start size-match and endpoint-boundary protection, N =" << N;

		// ------------------------------------------------------------
		// A. 尺寸匹配：在不造成端面越界和不明显离开参考背景的前提下，
		//    逐步放大每个hub连接面宽度，使它尽量接近对应输入端面宽度。
		// ------------------------------------------------------------
		for (int iter = 0; iter < 12; ++iter) {
			bool changed = false;

			for (int k = 0; k < N; ++k) {
				double ring[4][3];
				GetHubRingForTube(k, ring);

				double topW = std::sqrt(vtkMath::Distance2BetweenPoints(ring[0], ring[3]));
				double botW = std::sqrt(vtkMath::Distance2BetweenPoints(ring[1], ring[2]));
				double curW = 0.5 * (topW + botW);

				double desiredW = tubes[k].w * 0.92;
				if (curW >= desiredW || curW < 1e-8) continue;

				int prev_k = (k - 1 + N) % N;
				int ids[4] = { 2 * prev_k, 2 * prev_k + 1, 2 * k + 1, 2 * k };

				double expand = std::min((desiredW - curW) * 0.28, std::min(tubes[k].w, tubes[k].h) * 0.08);
				if (expand <= 1e-8) continue;

				auto snap = SnapshotPoints();
				for (int ii = 0; ii < 4; ++ii) {
					MovePointRadiallyInPlane(ids[ii], expand);
				}

				double worst = 0.0;
				int worstTube = -1;
				if (!CheckAllHubPointsSafe(&worst, &worstTube) || !RefSpaceAllowed(0.85)) {
					RestorePoints(snap);
					if (false) qDebug() << "[2DHubSafety][SizeMatch] stop expand port =" << tubes[k].index
						<< "sorted =" << k
						<< "curW =" << curW
						<< "desiredW =" << desiredW
						<< "reason = safety/background limit"
						<< "worstViolation =" << worst
						<< "worstTube =" << worstTube;
				}
				else {
					changed = true;
					double ring2[4][3];
					GetHubRingForTube(k, ring2);
					double newW = 0.5 * (
						std::sqrt(vtkMath::Distance2BetweenPoints(ring2[0], ring2[3])) +
						std::sqrt(vtkMath::Distance2BetweenPoints(ring2[1], ring2[2]))
						);
					if (false) qDebug() << "[2DHubSafety][SizeMatch] expand port =" << tubes[k].index
						<< "sorted =" << k
						<< "width" << curW << "->" << newW
						<< "target =" << desiredW;
				}
			}

			if (!changed) break;
		}

		// N=3 专用避让需要保持“靠近主端面的四个 hub 节点”不被挪动，
// 因此这里先准备 N=3 的固定节点/待避让节点集合，后面的最终硬校验也复用它。
		std::vector<int> n3FixedAIds;
		std::vector<int> n3MoveBIds;
		std::vector<int> n3CheckBPorts;
		int n3TargetA = -1;
		bool n3ClampReady = false;

		// ------------------------------------------------------------
		// B. 越界避让：优先整体避让；整体避让不安全时，单节点避让。
		// ------------------------------------------------------------
		if (N != 3) {
			for (int iter = 0; iter < 20; ++iter) {
				double worst = 0.0;
				int worstTube = -1;
				if (CheckAllHubPointsSafe(&worst, &worstTube)) {
					if (false) qDebug() << "[2DHubSafety][Avoid] safe after iter =" << iter;
					break;
				}

				if (worstTube < 0) break;

				double moveLen = std::max(worst * 1.08, SafeGapForTube(worstTube) * 0.25);
				double delta[3] = {
					-tubes[worstTube].dir_out[0] * moveLen,
					-tubes[worstTube].dir_out[1] * moveLen,
					-tubes[worstTube].dir_out[2] * moveLen
				};

				auto snap = SnapshotPoints();
				double oldHubCenter[3] = { hubCenter[0], hubCenter[1], hubCenter[2] };
				MoveAllHubPoints(delta);
				hubCenter[0] += delta[0];
				hubCenter[1] += delta[1];
				hubCenter[2] += delta[2];

				double afterWorst = 0.0;
				int afterTube = -1;
				bool wholeMoveOK = CheckAllHubPointsSafe(&afterWorst, &afterTube) && RefSpaceAllowed(1.10);

				if (wholeMoveOK) {
					if (false) qDebug() << "[2DHubSafety][Avoid] whole-hub retreat applied"
						<< "iter =" << iter
						<< "triggerPort =" << tubes[worstTube].index
						<< "sorted =" << worstTube
						<< "moveLen =" << moveLen
						<< "beforeViolation =" << worst
						<< "afterWorst =" << afterWorst;
					continue;
				}

				RestorePoints(snap);
				hubCenter[0] = oldHubCenter[0];
				hubCenter[1] = oldHubCenter[1];
				hubCenter[2] = oldHubCenter[2];

				// 单节点避让：参考三维分支的硬限制方式，直接按输入端面四角点入口限制回退共享节点。
				// 这里不再用端面中心平面，也不通过 PutHubRing 覆盖整圈，避免相邻端口共享点互相写回。
				int localAdjusted = 0;
				for (int pass = 0; pass < 4; ++pass) {
					int passAdjusted = 0;
					for (int k = 0; k < N; ++k) {
						double limitS = CapEntrySForTube(k) - SafeGapForTube(k);
						int prev_k = (k - 1 + N) % N;
						int ids[4] = { 2 * prev_k, 2 * prev_k + 1, 2 * k + 1, 2 * k };

						for (int ii = 0; ii < 4; ++ii) {
							int id = ids[ii];
							double p[3];
							pts->GetPoint(id, p);

							double violation = SCoordForTube(p, k) - limitS;
							if (violation <= 1e-7) continue;

							p[0] -= tubes[k].dir_out[0] * violation;
							p[1] -= tubes[k].dir_out[1] * violation;
							p[2] -= tubes[k].dir_out[2] * violation;
							pts->SetPoint(id, p);
							++passAdjusted;
						}
					}

					localAdjusted += passAdjusted;
					if (passAdjusted == 0) break;
				}

				if (false) qDebug() << "[2DHubSafety][Avoid] local exact-cap-node retreat applied"
					<< "iter =" << iter
					<< "triggerPort =" << tubes[worstTube].index
					<< "sorted =" << worstTube
					<< "beforeViolation =" << worst
					<< "adjustedNodes =" << localAdjusted;
			}



		}
		else {
			if (false) qDebug() << "[2DHubSafety][Avoid] skip generic whole/local retreat for N=3; use two-node exact clamp only.";
		}


		// ------------------------------------------------------------
		// C/D. N=3 专用主骨架投影避让：
		//   1) A端面不再按“危险程度”选，而是按骨架夹角选：
		//      三根骨架里夹角最小的一对 = 两个 B；剩下那根 = A 主骨架。
		//   2) B端面的限制不再沿各自 B.dir_out 检查，而是统一投影到 A 主骨架方向 mainAxis 上。
		//   3) 对每个 B，取其所属控制笼 8 个节点在 mainAxis 上的投影范围 [bMin,bMax]。
		//      P 是 hub 中靠近 B 端面的两个节点，也就是 A 环之外的两个节点。
		//      如果 P 的 mainAxis 投影落入/没有离开 B 的投影范围，则整体 hub 沿 mainAxis 朝 A 端面移动。
		//   4) 整体移动必须同时满足：
		//      - P 离开 B 的投影范围，至少到 bMax + gap；
		//      - hub 的两个 Q 面/所有 hub 点不能越过 A 端面的入口，最多到 aMin - gap；
		//      - 参考背景空间允许。
		// ------------------------------------------------------------
		int finalClampCount = 0;
		if (N == 3) {
			auto ClampValue = [](double v, double lo, double hi) -> double
			{
				return std::max(lo, std::min(v, hi));
			};

			auto AngleBetweenDirs = [](const double a[3], const double b[3]) -> double
			{
				double c = vtkMath::Dot(a, b);
				c = std::max(-1.0, std::min(1.0, c));
				return std::acos(c);
			};

			auto AxisS = [&](const double p[3], const double axis[3]) -> double
			{
				return vtkMath::Dot(p, axis);
			};

			auto GetEightPointRangeOnAxis = [&](int k, const double axis[3], double& sMin, double& sMax)
			{
				sMin = std::numeric_limits<double>::max();
				sMax = -std::numeric_limits<double>::max();
				for (int i = 0; i < 8; ++i) {
					double sVal = AxisS(tubes[k].cagePts8[i], axis);
					sMin = std::min(sMin, sVal);
					sMax = std::max(sMax, sVal);
				}
			};

			auto MaxHubSOnAxis = [&](const double axis[3]) -> double
			{
				double maxS = -std::numeric_limits<double>::max();
				for (vtkIdType id = 0; id < pts->GetNumberOfPoints(); ++id) {
					double pnt[3];
					pts->GetPoint(id, pnt);
					maxS = std::max(maxS, AxisS(pnt, axis));
				}
				return maxS;
			};

			auto MinPositiveGap = [&]() -> double
			{
				double g = std::numeric_limits<double>::max();
				for (int k = 0; k < N; ++k) {
					g = std::min(g, std::min(tubes[k].w, tubes[k].h));
				}
				if (!std::isfinite(g) || g < 1e-8) g = 1.0;
				return std::max(g * 0.08, 1e-4);
			};

			// A/B 判定：三条真实骨架分叉线里，夹角最小的那一对是两个 B 端面；
			// 剩下那一条就是 A 主骨架。
			// 这里使用 branchDir：骨架交点 -> 该分支第一个相邻骨架节点。
			// 如果骨架图绑定失败，branchDir 已经在前面回退为控制笼中心线 skeletonDir。
			double angle01 = AngleBetweenDirs(tubes[0].branchDir, tubes[1].branchDir);
			double angle12 = AngleBetweenDirs(tubes[1].branchDir, tubes[2].branchDir);
			double angle20 = AngleBetweenDirs(tubes[2].branchDir, tubes[0].branchDir);

			double bestAngleScore = angle01; // 这里保留变量名，只用于日志输出；含义改为“最小 B-B 夹角”。
			int b0 = 0;
			int b1 = 1;
			if (angle12 < bestAngleScore) {
				bestAngleScore = angle12;
				b0 = 1;
				b1 = 2;
			}
			if (angle20 < bestAngleScore) {
				bestAngleScore = angle20;
				b0 = 2;
				b1 = 0;
			}

			for (int k = 0; k < 3; ++k) {
				if (k != b0 && k != b1) {
					n3TargetA = k;
					break;
				}
			}

			if (n3TargetA >= 0) {
				int prevA = (n3TargetA - 1 + N) % N;
				std::array<int, 4> aIds = { 2 * prevA, 2 * prevA + 1, 2 * n3TargetA + 1, 2 * n3TargetA };
				std::array<bool, 6> isA = { false, false, false, false, false, false };

				for (int id : aIds) {
					if (id >= 0 && id < 6) {
						isA[static_cast<size_t>(id)] = true;
						n3FixedAIds.push_back(id);
					}
				}
				for (int id = 0; id < 6; ++id) {
					if (!isA[static_cast<size_t>(id)]) n3MoveBIds.push_back(id);
				}
				for (int k = 0; k < N; ++k) {
					if (k != n3TargetA) n3CheckBPorts.push_back(k);
				}
				n3ClampReady = (n3MoveBIds.size() == 2 && n3CheckBPorts.size() == 2);

				// P / B范围 / Q限制 / hub整体移动的投影轴，也必须使用 A 的真实骨架分叉线。
				double mainAxis[3] = {
					tubes[n3TargetA].branchDir[0],
					tubes[n3TargetA].branchDir[1],
					tubes[n3TargetA].branchDir[2]
				};
				vtkMath::Normalize(mainAxis);

				double mainGap = MinPositiveGap();
				double aMin = 0.0, aMax = 0.0;
				GetEightPointRangeOnAxis(n3TargetA, mainAxis, aMin, aMax);

				auto ComputeNeedMoveForP = [&]() -> double
				{
					double need = 0.0;
					for (int bk : n3CheckBPorts) {
						double bMin = 0.0, bMax = 0.0;
						GetEightPointRangeOnAxis(bk, mainAxis, bMin, bMax);
						double requiredP = bMax + mainGap;

						for (int id : n3MoveBIds) {
							double pnt[3];
							pts->GetPoint(id, pnt);
							double pS = AxisS(pnt, mainAxis);
							// pS 落在 [bMin,bMax] 内，或者还没越过 bMax+gap，都要整体朝 A 推。
							need = std::max(need, requiredP - pS);
						}
					}
					return std::max(0.0, need);
				};

				auto ComputeAllowedMoveToA = [&]() -> double
				{
					double qLimit = aMin - mainGap;
					double maxHubS = MaxHubSOnAxis(mainAxis);
					return std::max(0.0, qLimit - maxHubS);
				};

				auto TryWholeMoveOnMainAxis = [&](double moveLen) -> bool
				{
					if (moveLen <= 1e-9) return true;
					double delta[3] = {
						mainAxis[0] * moveLen,
						mainAxis[1] * moveLen,
						mainAxis[2] * moveLen
					};
					auto snap = SnapshotPoints();
					double oldHubCenter[3] = { hubCenter[0], hubCenter[1], hubCenter[2] };
					MoveAllHubPoints(delta);
					hubCenter[0] += delta[0];
					hubCenter[1] += delta[1];
					hubCenter[2] += delta[2];

					bool ok = (ComputeAllowedMoveToA() >= -1e-7) && RefSpaceAllowed(1.20);
					if (!ok) {
						RestorePoints(snap);
						hubCenter[0] = oldHubCenter[0];
						hubCenter[1] = oldHubCenter[1];
						hubCenter[2] = oldHubCenter[2];
						return false;
					}
					return true;
				};

				if (n3ClampReady) {
					double needMove = ComputeNeedMoveForP();
					double allowMove = ComputeAllowedMoveToA();
					double appliedMove = std::min(needMove, allowMove);

					if (appliedMove > 1e-8) {
						if (TryWholeMoveOnMainAxis(appliedMove)) {
							++finalClampCount;
						}
					}

					// 保守巡边：上一版会一直朝 A 方向推到 A 限制或 RefSpaceAllowed 失败，
					// 但最近点距离不能判断“内外侧”，容易让 Q 点贴着背景外侧跑出模型。
					// 这里只允许在 P 已安全后做少量额外整体探索，避免 Q 过冲。
					int exploreSteps = 0;
					double exploredLen = 0.0;
					double maxExploreLen = mainGap * 0.45;
					for (int step = 0; step < 6; ++step) {
						double residualNeed = ComputeNeedMoveForP();
						if (residualNeed > 1e-6) break;

						double allow = ComputeAllowedMoveToA();
						if (allow <= 1e-8) break;
						if (exploredLen >= maxExploreLen) break;

						double stepLen = std::min(std::min(allow, mainGap * 0.08), maxExploreLen - exploredLen);
						if (stepLen <= 1e-9) break;
						if (!TryWholeMoveOnMainAxis(stepLen)) break;
						exploredLen += stepLen;
						++exploreSteps;
					}

					// 最终硬校验：只做主骨架投影校验，不再沿 B.dir_out 单点拉扯，避免把 P 拉回错误方向。
					double finalNeed = ComputeNeedMoveForP();
					double finalAllow = ComputeAllowedMoveToA();
					if (finalNeed > 1e-7 && finalAllow > 1e-8) {
						double addMove = std::min(finalNeed, finalAllow);
						if (TryWholeMoveOnMainAxis(addMove)) ++finalClampCount;
					}

					double remainNeed = ComputeNeedMoveForP();
					double remainAllow = ComputeAllowedMoveToA();

					pts->Modified();
					if (false) qDebug() << "[2DHubSafety][N3MainAxisAvoid] applied"
						<< "A(sorted/input) =" << n3TargetA << tubes[n3TargetA].index
						<< "B(sorted) =" << (n3CheckBPorts.size() > 0 ? n3CheckBPorts[0] : -1)
						<< (n3CheckBPorts.size() > 1 ? n3CheckBPorts[1] : -1)
						<< "Pnodes =" << (n3MoveBIds.size() > 0 ? n3MoveBIds[0] : -1)
						<< (n3MoveBIds.size() > 1 ? n3MoveBIds[1] : -1)
						<< "Q/A nodes =" << (n3FixedAIds.size() > 0 ? n3FixedAIds[0] : -1)
						<< (n3FixedAIds.size() > 1 ? n3FixedAIds[1] : -1)
						<< (n3FixedAIds.size() > 2 ? n3FixedAIds[2] : -1)
						<< (n3FixedAIds.size() > 3 ? n3FixedAIds[3] : -1)
						<< "angleScore =" << bestAngleScore
						<< "initialNeed =" << needMove
						<< "initialAllow =" << allowMove
						<< "applied =" << appliedMove
						<< "exploreSteps =" << exploreSteps
						<< "remainNeed =" << remainNeed
						<< "remainAllowToA =" << remainAllow;

					// ------------------------------------------------------------
					// N=3 紫色边界保护：A端与两个B端分别构成一条侧边界。
					// Q/A侧 hub 节点不允许越过这两条 A-B 边界。
					// 边界构造方式：
					//   对每个B端，取 A端四角中距离B中心最近的角点，
					//   再取 B端四角中距离A中心最近的角点，两点连成 A-B 边界线。
					//   以当前hubCenter所在侧作为安全侧；Q点若越线，则拉回到边界安全侧。
					// ------------------------------------------------------------
					int purpleClampCount = 0;
					auto ClosestCapPointTo = [&](int portId, const double target[3], double out[3])
					{
						double bestD2 = std::numeric_limits<double>::max();
						int bestI = 0;
						for (int ii = 0; ii < 4; ++ii) {
							double d2 = vtkMath::Distance2BetweenPoints(tubes[portId].capRing[ii], target);
							if (d2 < bestD2) {
								bestD2 = d2;
								bestI = ii;
							}
						}
						for (int d = 0; d < 3; ++d) out[d] = tubes[portId].capRing[bestI][d];
					};

					auto ClampPointToBoundarySide = [&](int pointId, const double a[3], const double b[3]) -> bool
					{
						if (pointId < 0 || pointId >= pts->GetNumberOfPoints()) return false;

						double lineDir[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
						if (vtkMath::Normalize(lineDir) < 1e-10) return false;

						double sideNormal[3];
						vtkMath::Cross(normalUp, lineDir, sideNormal);
						if (vtkMath::Normalize(sideNormal) < 1e-10) return false;

						auto SideValue = [&](const double p[3]) -> double
						{
							double v[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
							return vtkMath::Dot(v, sideNormal);
						};

						double safeSide = SideValue(hubCenter);
						if (std::abs(safeSide) < 1e-9) {
							safeSide = SideValue(tubes[n3TargetA].center);
						}
						const double sign = (safeSide >= 0.0) ? 1.0 : -1.0;
						const double minSide = sign * mainGap * 0.10;

						double pnt[3];
						pts->GetPoint(pointId, pnt);
						double curSide = SideValue(pnt);
						if (sign * curSide >= sign * minSide - 1e-8) return false;

						double delta = minSide - curSide;
						pnt[0] += sideNormal[0] * delta;
						pnt[1] += sideNormal[1] * delta;
						pnt[2] += sideNormal[2] * delta;
						pts->SetPoint(pointId, pnt);
						return true;
					};

					for (int bk : n3CheckBPorts) {
						double aBoundary[3], bBoundary[3];
						ClosestCapPointTo(n3TargetA, tubes[bk].center, aBoundary);
						ClosestCapPointTo(bk, tubes[n3TargetA].center, bBoundary);

						for (int qid : n3FixedAIds) {
							if (ClampPointToBoundarySide(qid, aBoundary, bBoundary)) {
								++purpleClampCount;
							}
						}
					}

					if (purpleClampCount > 0) {
						pts->Modified();
						if (false) qDebug() << "[2DHubSafety][N3PurpleBoundaryClamp] Q nodes clamped =" << purpleClampCount;
					}


					// ------------------------------------------------------------
					// N=3 二次巡边：在“不越过 A-B 紫色边界”的前提下，继续沿 A 主骨架方向
					// 探索背景尺寸边界。上一阶段先解决 P/B 避让，这里只做少量试探，
					// 每一步都用 RefSpaceAllowed 与 A-B 边界同时验算，失败立即回退。
					// ------------------------------------------------------------
					auto PointOnSafeSideOfABBoundary = [&](int pointId, const double a[3], const double b[3]) -> bool
					{
						if (pointId < 0 || pointId >= pts->GetNumberOfPoints()) return true;
						double lineDir[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
						if (vtkMath::Normalize(lineDir) < 1e-10) return true;

						double sideNormal[3];
						vtkMath::Cross(normalUp, lineDir, sideNormal);
						if (vtkMath::Normalize(sideNormal) < 1e-10) return true;

						auto SideValue = [&](const double p[3]) -> double
						{
							double v[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
							return vtkMath::Dot(v, sideNormal);
						};

						double safeSide = SideValue(hubCenter);
						if (std::abs(safeSide) < 1e-9) {
							safeSide = SideValue(tubes[n3TargetA].center);
						}
						const double sign = (safeSide >= 0.0) ? 1.0 : -1.0;
						const double minSide = sign * mainGap * 0.06;

						double pnt[3];
						pts->GetPoint(pointId, pnt);
						return sign * SideValue(pnt) >= sign * minSide - 1e-8;
					};

					auto AllQNodesInsidePurpleBoundary = [&]() -> bool
					{
						for (int bk : n3CheckBPorts) {
							double aBoundary[3], bBoundary[3];
							ClosestCapPointTo(n3TargetA, tubes[bk].center, aBoundary);
							ClosestCapPointTo(bk, tubes[n3TargetA].center, bBoundary);
							for (int qid : n3FixedAIds) {
								if (!PointOnSafeSideOfABBoundary(qid, aBoundary, bBoundary)) {
									return false;
								}
							}
						}
						return true;
					};

					auto TryMoveABSafe = [&](double moveLen) -> bool
					{
						if (moveLen <= 1e-9) return true;
						double delta[3] = {
							mainAxis[0] * moveLen,
							mainAxis[1] * moveLen,
							mainAxis[2] * moveLen
						};

						auto snap = SnapshotPoints();
						double oldHubCenter2[3] = { hubCenter[0], hubCenter[1], hubCenter[2] };
						MoveAllHubPoints(delta);
						hubCenter[0] += delta[0];
						hubCenter[1] += delta[1];
						hubCenter[2] += delta[2];

						const bool ok =
							RefSpaceAllowed(1.20) &&
							AllQNodesInsidePurpleBoundary() &&
							ComputeAllowedMoveToA() >= -1e-7;

						if (!ok) {
							RestorePoints(snap);
							hubCenter[0] = oldHubCenter2[0];
							hubCenter[1] = oldHubCenter2[1];
							hubCenter[2] = oldHubCenter2[2];
							return false;
						}
						return true;
					};

					int abSafeExploreSteps = 0;
					double abSafeExploreLen = 0.0;
					if (ComputeNeedMoveForP() <= 1e-6 && AllQNodesInsidePurpleBoundary()) {
						const double maxABSafeExplore = mainGap * 0.65;
						for (int step = 0; step < 10; ++step) {
							double allow = ComputeAllowedMoveToA();
							if (allow <= 1e-8) break;
							if (abSafeExploreLen >= maxABSafeExplore) break;
							double stepLen = std::min(std::min(allow, mainGap * 0.07), maxABSafeExplore - abSafeExploreLen);
							if (stepLen <= 1e-9) break;
							if (!TryMoveABSafe(stepLen)) break;
							abSafeExploreLen += stepLen;
							++abSafeExploreSteps;
						}
					}

					if (abSafeExploreSteps > 0) {
						pts->Modified();
						if (false) qDebug() << "[2DHubSafety][N3ABSafeBoundarySearch] steps =" << abSafeExploreSteps
							<< "length =" << abSafeExploreLen;
					}

					// 最终几何约束只修正 N=3 hub 的 6 个顶点，不改变骨架交点、
					// 端口排序、三角柱拓扑或后续 bridge 构造。
					// A 接口作为刚性基准；剩余两个点由 B0/B1 共同拟合，并限制在
					// 两个 B 端面相距最远的边所张成的平面以内。
					if (n3ClampReady) {
						auto HubIdsForTube = [&](int portId)->std::array<int, 4>
						{
							int prev = (portId - 1 + N) % N;
							return { 2 * prev, 2 * prev + 1, 2 * portId + 1, 2 * portId };
						};

						auto IsBsidePoint = [&](int pointId) -> bool
						{
							return std::find(n3MoveBIds.begin(), n3MoveBIds.end(), pointId) != n3MoveBIds.end();
						};

						// 1. A 接口严格复制 A 端面，只允许整体平移到当前 hub A 接口中心。
						double oldARing[4][3];
						GetHubRingForTube(n3TargetA, oldARing);
						double oldACenter[3] = { 0.0, 0.0, 0.0 };
						double capACenter[3] = { 0.0, 0.0, 0.0 };
						for (int i = 0; i < 4; ++i) {
							for (int d = 0; d < 3; ++d) {
								oldACenter[d] += oldARing[i][d] * 0.25;
								capACenter[d] += tubes[n3TargetA].capRing[i][d] * 0.25;
							}
						}

						double bestARing[4][3] = {};
						double bestAScore = std::numeric_limits<double>::max();
						for (int reverse = 0; reverse < 2; ++reverse) {
							for (int shift = 0; shift < 4; ++shift) {
								double candidate[4][3];
								int order[4];
								BuildTubeCandidateRing(tubes[n3TargetA].capRing, reverse != 0, shift, candidate, order);
								double score = 0.0;
								for (int i = 0; i < 4; ++i) {
									for (int d = 0; d < 3; ++d) candidate[i][d] += oldACenter[d] - capACenter[d];
									score += vtkMath::Distance2BetweenPoints(candidate[i], oldARing[i]);
								}
								if (score < bestAScore) {
									bestAScore = score;
									CopyRing(bestARing, candidate);
								}
							}
						}
						PutHubRingForTube(n3TargetA, bestARing);

						// 2. 每个 B 端面先用它与 A 共享的那条边求刚性平移，
						// 再给两个公共 B 侧点各产生一个目标；两个 B 的目标取平均。
						std::array<std::array<double, 3>, 2> bTargets = {
							std::array<double, 3>{ 0.0, 0.0, 0.0 },
							std::array<double, 3>{ 0.0, 0.0, 0.0 }
						};
						std::array<int, 2> bTargetCounts = { 0, 0 };

						for (int bk : n3CheckBPorts) {
							auto ids = HubIdsForTube(bk);
							double currentRing[4][3];
							GetHubRingForTube(bk, currentRing);
							double bestCandidate[4][3] = {};
							double bestScore = std::numeric_limits<double>::max();

							for (int reverse = 0; reverse < 2; ++reverse) {
								for (int shift = 0; shift < 4; ++shift) {
									double candidate[4][3];
									int order[4];
									BuildTubeCandidateRing(tubes[bk].capRing, reverse != 0, shift, candidate, order);
									double translation[3] = { 0.0, 0.0, 0.0 };
									int fixedCount = 0;
									for (int i = 0; i < 4; ++i) {
										if (IsBsidePoint(ids[i])) continue;
										for (int d = 0; d < 3; ++d) translation[d] += currentRing[i][d] - candidate[i][d];
										++fixedCount;
									}
									if (fixedCount == 0) continue;
									for (int d = 0; d < 3; ++d) translation[d] /= static_cast<double>(fixedCount);

									double score = 0.0;
									for (int i = 0; i < 4; ++i) {
										for (int d = 0; d < 3; ++d) candidate[i][d] += translation[d];
										const double weight = IsBsidePoint(ids[i]) ? 0.05 : 1.0;
										score += weight * vtkMath::Distance2BetweenPoints(candidate[i], currentRing[i]);
									}
									if (score < bestScore) {
										bestScore = score;
										CopyRing(bestCandidate, candidate);
									}
								}
							}

							for (int i = 0; i < 4; ++i) {
								for (int pi = 0; pi < 2; ++pi) {
									if (ids[i] != n3MoveBIds[pi]) continue;
									for (int d = 0; d < 3; ++d) bTargets[pi][d] += bestCandidate[i][d];
									++bTargetCounts[pi];
								}
							}
						}

						for (int pi = 0; pi < 2; ++pi) {
							if (bTargetCounts[pi] <= 0) continue;
							double target[3];
							for (int d = 0; d < 3; ++d) target[d] = bTargets[pi][d] / static_cast<double>(bTargetCounts[pi]);
							pts->SetPoint(n3MoveBIds[pi], target);
						}

						// A 环给出唯一的厚度向量。三个上点和三个下点必须使用同一个向量，
						// 因而上下三角形严格全等、平行，不允许 B0/B1 分别把上下点拉到不同平面。
						double prismThickness[3] = {
							0.5 * ((bestARing[0][0] - bestARing[1][0]) + (bestARing[3][0] - bestARing[2][0])),
							0.5 * ((bestARing[0][1] - bestARing[1][1]) + (bestARing[3][1] - bestARing[2][1])),
							0.5 * ((bestARing[0][2] - bestARing[1][2]) + (bestARing[3][2] - bestARing[2][2]))
						};
						double thicknessLength = vtkMath::Norm(prismThickness);
						if (thicknessLength < 1e-8) {
							thicknessLength = std::max(tubes[n3TargetA].h, 1e-4);
							for (int d = 0; d < 3; ++d) prismThickness[d] = normalUp[d] * thicknessLength;
						}
						double thicknessAxis[3] = {
							prismThickness[0] / thicknessLength,
							prismThickness[1] / thicknessLength,
							prismThickness[2] / thicknessLength
						};

						double bBase[3] = { 0.0, 0.0, 0.0 };
						for (int pointId : n3MoveBIds) {
							double p[3];
							pts->GetPoint(pointId, p);
							for (int d = 0; d < 3; ++d) bBase[d] += p[d] * 0.5;
						}
						for (int pointId : n3MoveBIds) {
							double p[3];
							const double side = (pointId % 2 == 0) ? 0.5 : -0.5;
							for (int d = 0; d < 3; ++d) p[d] = bBase[d] + side * prismThickness[d];
							pts->SetPoint(pointId, p);
						}

						// 六个顶点的几何中心严格回到骨架交点。只做统一平移，不改变 A 面尺寸/法向，
						// 也不改变上下三角形的相对关系。
						double prismCenter[3] = { 0.0, 0.0, 0.0 };
						for (int pointId = 0; pointId < 6; ++pointId) {
							double p[3];
							pts->GetPoint(pointId, p);
							for (int d = 0; d < 3; ++d) prismCenter[d] += p[d] / 6.0;
						}
						double centerCorrection[3] = {
							n3LockedHubCenter[0] - prismCenter[0],
							n3LockedHubCenter[1] - prismCenter[1],
							n3LockedHubCenter[2] - prismCenter[2]
						};
						MoveAllHubPoints(centerCorrection);
						for (int d = 0; d < 3; ++d) hubCenter[d] = n3LockedHubCenter[d];

						// 3. 从两个 B 端面各取一条边，选择边中点平面内距离最大的组合。
						// 两条边与其中点连线张成限制平面；B 侧上下点只能位于 A 中心同侧。
						struct EdgeInfo {
							double p0[3];
							double p1[3];
							double mid[3];
							double dir[3];
						};
						auto MakeEdge = [&](int portId, int edgeId) -> EdgeInfo
						{
							EdgeInfo edge = {};
							int next = (edgeId + 1) % 4;
							for (int d = 0; d < 3; ++d) {
								edge.p0[d] = tubes[portId].capRing[edgeId][d];
								edge.p1[d] = tubes[portId].capRing[next][d];
								edge.mid[d] = 0.5 * (edge.p0[d] + edge.p1[d]);
								edge.dir[d] = edge.p1[d] - edge.p0[d];
							}
							vtkMath::Normalize(edge.dir);
							return edge;
						};

						EdgeInfo farEdge0 = {};
						EdgeInfo farEdge1 = {};
						double farthestD2 = -1.0;
						for (int e0 = 0; e0 < 4; ++e0) {
							EdgeInfo edge0 = MakeEdge(n3CheckBPorts[0], e0);
							for (int e1 = 0; e1 < 4; ++e1) {
								EdgeInfo edge1 = MakeEdge(n3CheckBPorts[1], e1);
								double connector[3] = {
									edge1.mid[0] - edge0.mid[0],
									edge1.mid[1] - edge0.mid[1],
									edge1.mid[2] - edge0.mid[2]
								};
								double upPart = vtkMath::Dot(connector, thicknessAxis);
								for (int d = 0; d < 3; ++d) connector[d] -= thicknessAxis[d] * upPart;
								double d2 = vtkMath::Dot(connector, connector);
								if (d2 > farthestD2) {
									farthestD2 = d2;
									farEdge0 = edge0;
									farEdge1 = edge1;
								}
							}
						}

						if (farthestD2 > 1e-12) {
							double across[3] = {
								farEdge1.mid[0] - farEdge0.mid[0],
								farEdge1.mid[1] - farEdge0.mid[1],
								farEdge1.mid[2] - farEdge0.mid[2]
							};
							double boundaryNormal[3];
							// 限制平面沿统一厚度方向张开，因此它的法向严格位于 hub 中面内；
							// 回退不会再产生向上/向下分量。
							vtkMath::Cross(thicknessAxis, across, boundaryNormal);
							if (vtkMath::Normalize(boundaryNormal) > 1e-10) {
								auto PlaneSide = [&](const double p[3]) -> double
								{
									double v[3] = { p[0] - farEdge0.mid[0], p[1] - farEdge0.mid[1], p[2] - farEdge0.mid[2] };
									return vtkMath::Dot(v, boundaryNormal);
								};
								double currentACenter[3] = { 0.0, 0.0, 0.0 };
								for (int pointId : n3FixedAIds) {
									double p[3];
									pts->GetPoint(pointId, p);
									for (int d = 0; d < 3; ++d) currentACenter[d] += p[d] * 0.25;
								}
								double safeSide = PlaneSide(currentACenter);
								if (std::abs(safeSide) < 1e-9) safeSide = PlaneSide(tubes[n3TargetA].center);
								double sign = safeSide >= 0.0 ? 1.0 : -1.0;
								// 橙色平面是可用空间的极限，不是需要大幅退让的安全带。
								// 只保留数值容差，目标是让 B 侧上下节点主动贴近该平面。
								double boundaryGap = std::max(mainGap * 0.005, 1e-5);

								// 两个 B 点共用同一中面坐标，所以只计算一次平面内修正。
								double currentBBase[3] = { 0.0, 0.0, 0.0 };
								for (int pointId : n3MoveBIds) {
									double p[3];
									pts->GetPoint(pointId, p);
									for (int d = 0; d < 3; ++d) currentBBase[d] += p[d] * 0.5;
								}
								double signedSide = sign * PlaneSide(currentBBase);
								// 旧逻辑只处理越界(signedSide < gap)，因此合法但过于靠内的点不会
								// 向橙色极限平面展开。这里直接对齐目标距离，正负修正均允许。
								double correction = boundaryGap - signedSide;
								if (std::abs(correction) > 1e-8) {
									double bDelta[3] = {
										sign * boundaryNormal[0] * correction,
										sign * boundaryNormal[1] * correction,
										sign * boundaryNormal[2] * correction
									};
									for (int pointId : n3MoveBIds) {
										double p[3];
										pts->GetPoint(pointId, p);
										for (int d = 0; d < 3; ++d) p[d] += bDelta[d];
										pts->SetPoint(pointId, p);
									}
								}
							}
						}

						// 4. A 强制面在其对应截面内直接贴合背景外轮廓。
						// 初始面仍是 A 输入面的平移副本；这里只让左右两条边各自沿宽度轴
						// 对齐本侧最外截面交点，长了收回、短了伸出，厚度和法向不变。
						if (m_RefData && m_RefData->GetNumberOfPoints() > 0 &&
							m_RefData->GetNumberOfCells() > 0)
						{
							std::vector<int> aPairIds;
							for (int pointId : n3FixedAIds) {
								int pairId = pointId / 2;
								if (std::find(aPairIds.begin(), aPairIds.end(), pairId) == aPairIds.end()) {
									aPairIds.push_back(pairId);
								}
							}

							if (aPairIds.size() == 2) {
								auto PairBase = [&](int pairId, double base[3])
								{
									double top[3], bottom[3];
									pts->GetPoint(2 * pairId, top);
									pts->GetPoint(2 * pairId + 1, bottom);
									for (int d = 0; d < 3; ++d) base[d] = 0.5 * (top[d] + bottom[d]);
								};
								auto MovePairBase = [&](int pairId, const double newBase[3])
								{
									double top[3], bottom[3], oldBase[3];
									pts->GetPoint(2 * pairId, top);
									pts->GetPoint(2 * pairId + 1, bottom);
									for (int d = 0; d < 3; ++d) {
										oldBase[d] = 0.5 * (top[d] + bottom[d]);
										double delta = newBase[d] - oldBase[d];
										top[d] += delta;
										bottom[d] += delta;
									}
									pts->SetPoint(2 * pairId, top);
									pts->SetPoint(2 * pairId + 1, bottom);
								};

								double aBase0[3], aBase1[3];
								PairBase(aPairIds[0], aBase0);
								PairBase(aPairIds[1], aBase1);
								double aFaceCenter[3] = {
									0.5 * (aBase0[0] + aBase1[0]),
									0.5 * (aBase0[1] + aBase1[1]),
									0.5 * (aBase0[2] + aBase1[2])
								};
								double widthAxis[3] = {
									aBase1[0] - aBase0[0],
									aBase1[1] - aBase0[1],
									aBase1[2] - aBase0[2]
								};
								double halfWidth = 0.5 * vtkMath::Norm(widthAxis);
								vtkMath::Normalize(widthAxis);

								if (halfWidth > 1e-8) {
									vtkSmartPointer<vtkOBBTree> refTree = vtkSmartPointer<vtkOBBTree>::New();
									refTree->SetDataSet(m_RefData);
									refTree->BuildLocator();

									double bounds[6];
									m_RefData->GetBounds(bounds);
									double refDiagonal = std::sqrt(
										(bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
										(bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
										(bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));
									// 只查询 A 面附近的对应截面，避免射线穿到模型另一处无关分支。
									double maxRayLength = std::min(refDiagonal * 2.0, std::max(halfWidth * 4.0, refDiagonal * 0.05));
									double widthNeg[3] = { -widthAxis[0], -widthAxis[1], -widthAxis[2] };
									double wPos = 0.0, wNeg = 0.0;
									bool hasWPos = AutoStraightRayFarthestDistanceToSurface(refTree, aFaceCenter, widthAxis, maxRayLength, wPos);
									bool hasWNeg = AutoStraightRayFarthestDistanceToSurface(refTree, aFaceCenter, widthNeg, maxRayLength, wNeg);

									const double tangentTolerance = std::max(refDiagonal * 1e-5, 1e-6);
									bool movedNegativeEdge = false;
									bool movedPositiveEdge = false;

									// aPairIds[0] 位于 widthAxis 负侧，直接贴到负侧截面外轮廓。
									if (hasWNeg) {
										double allowed = std::max(wNeg - tangentTolerance, 1e-6);
										if (std::abs(halfWidth - allowed) > 1e-8) {
											double newBase[3];
											for (int d = 0; d < 3; ++d) newBase[d] = aFaceCenter[d] - widthAxis[d] * allowed;
											MovePairBase(aPairIds[0], newBase);
											movedNegativeEdge = true;
										}
									}

									// aPairIds[1] 位于 widthAxis 正侧，独立贴到正侧截面外轮廓。
									if (hasWPos) {
										double allowed = std::max(wPos - tangentTolerance, 1e-6);
										if (std::abs(halfWidth - allowed) > 1e-8) {
											double newBase[3];
											for (int d = 0; d < 3; ++d) newBase[d] = aFaceCenter[d] + widthAxis[d] * allowed;
											MovePairBase(aPairIds[1], newBase);
											movedPositiveEdge = true;
										}
									}

									if (false && (movedNegativeEdge || movedPositiveEdge)) {
										qDebug() << "[2DHubSafety][N3AdaptiveAEdgeClamp]"
											<< "negativeEdge =" << movedNegativeEdge
											<< "positiveEdge =" << movedPositiveEdge
											<< "halfWidth =" << halfWidth
											<< "available(-,+) =" << wNeg << wPos;
									}
								}
							}
						}

						pts->Modified();
						if (false) qDebug() << "[2DHubSafety][N3RigidAAndBFarEdgePlane] applied"
							<< "A(sorted/input) =" << n3TargetA << tubes[n3TargetA].index
							<< "B(sorted) =" << n3CheckBPorts[0] << n3CheckBPorts[1]
							<< "A rigid-fit error =" << bestAScore
							<< "far-edge distance2 =" << farthestD2
							<< "locked to skeleton junction =" << n3HasLockedSkeletonJunction;
					}
				}
			}
		}
		else {
			// N != 3 保持原来的所有 ring 节点硬限制。
			for (int pass = 0; pass < 4; ++pass) {
				int passClamp = 0;
				for (int k = 0; k < N; ++k) {
					double limitS = CapEntrySForTube(k) - SafeGapForTube(k);
					int prev_k = (k - 1 + N) % N;
					int ids[4] = { 2 * prev_k, 2 * prev_k + 1, 2 * k + 1, 2 * k };
					for (int ii = 0; ii < 4; ++ii) {
						int id = ids[ii];
						double pnt[3];
						pts->GetPoint(id, pnt);
						double violation = SCoordForTube(pnt, k) - limitS;
						if (violation <= 1e-7) continue;
						pnt[0] -= tubes[k].dir_out[0] * violation;
						pnt[1] -= tubes[k].dir_out[1] * violation;
						pnt[2] -= tubes[k].dir_out[2] * violation;
						pts->SetPoint(id, pnt);
						++passClamp;
					}
				}
				finalClampCount += passClamp;
				if (passClamp == 0) break;
			}
			if (finalClampCount > 0) {
				pts->Modified();
				if (false) qDebug() << "[2DHubSafety][FinalExactCapClamp] applied"
					<< "clampedNodes =" << finalClampCount;
			}
		}

		for (int k = 0; k < N; ++k) {
			double ring[4][3];
			GetHubRingForTube(k, ring);
			double topW = std::sqrt(vtkMath::Distance2BetweenPoints(ring[0], ring[3]));
			double botW = std::sqrt(vtkMath::Distance2BetweenPoints(ring[1], ring[2]));
			double leftH = std::sqrt(vtkMath::Distance2BetweenPoints(ring[0], ring[1]));
			double rightH = std::sqrt(vtkMath::Distance2BetweenPoints(ring[3], ring[2]));

			double maxT = -1e100;
			double maxExactViolation = -1e100;
			double limitS = CapEntrySForTube(k) - SafeGapForTube(k);
			for (int i = 0; i < 4; ++i) {
				maxT = std::max(maxT, PointPlaneT(ring[i], k));
				maxExactViolation = std::max(maxExactViolation, SCoordForTube(ring[i], k) - limitS);
			}

			if (false) qDebug() << "[2DHubSafety][Final] port =" << tubes[k].index
				<< "sorted =" << k
				<< "hubW(avg) =" << 0.5 * (topW + botW)
				<< "inputW =" << tubes[k].w
				<< "hubH(avg) =" << 0.5 * (leftH + rightH)
				<< "inputH =" << tubes[k].h
				<< "maxPlaneT(centerDiag) =" << maxT
				<< "centerRequired <=" << -SafeGapForTube(k)
				<< "exactCapViolation <=0 =" << maxExactViolation
				<< "capEntryS =" << CapEntrySForTube(k)
				<< "capLimitS =" << limitS;
		}
	}

	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

	// 
	polys->InsertNextCell(N);
	for (int k = 0; k < N; ++k) {
		polys->InsertCellPoint(k * 2);
	}

	// 
	polys->InsertNextCell(N);
	for (int k = N - 1; k >= 0; --k) {
		polys->InsertCellPoint(k * 2 + 1);
	}

	hubPoly->SetPoints(pts);
	hubPoly->SetPolys(polys);
	m_hubCages.push_back(hubPoly);

	if (false) qDebug() << "--- 原始接头已生成 ---";
	if (false) qDebug() << "接头中心 hubCenter =" << hubCenter[0] << hubCenter[1] << hubCenter[2];
	if (false) qDebug() << "上方向 normalUp =" << normalUp[0] << normalUp[1] << normalUp[2];
	if (false) qDebug() << "接头点数量 =" << pts->GetNumberOfPoints();
	if (false) qDebug() << "接头面 = 上表面/下表面";

	// ============================================================
	// 6. 
	// ============================================================

	for (int k = 0; k < N; ++k) {
		double tubeCanonical[4][3];
		double hubRing[4][3];

		// 6.1 
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				tubeCanonical[i][d] = tubes[k].capRing[i][d];
			}
		}

		// 6.2 
		//
		// hubRing[0] = U(prev/top)
		// hubRing[1] = D(prev/bottom)
		// hubRing[2] = D(k/bottom)
		// hubRing[3] = U(k/top)
		int prev_k = (k - 1 + N) % N;

		pts->GetPoint(2 * prev_k, hubRing[0]);
		pts->GetPoint(2 * prev_k + 1, hubRing[1]);
		pts->GetPoint(2 * k + 1, hubRing[2]);
		pts->GetPoint(2 * k, hubRing[3]);

		if (false) qDebug() << "\n--- [桥接候选搜索] 排序后序号 =" << k
			<< "original tube index =" << tubes[k].index
			<< "selected tubeIsStart =" << tubes[k].isStart
			<< "candidate family = forward + reverse, total 8";

		double bestTubeRing[4][3];
		int bestOrder[4] = { 0, 1, 2, 3 };
		bool bestReverse = false;

		BridgeScore bestScore;
		bestScore.minEdgeDot = -1e100;
		bestScore.avgEdgeDot = -1e100;
		bestScore.totalLength = 1e100;
		bestScore.maxDiagDiff = 1e100;

		// 6.3 
		for (int family = 0; family < 2; ++family) {
			bool reverse = (family == 1);

			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				int order[4];

				BuildTubeCandidateRing(
					tubeCanonical,
					reverse,
					shift,
					candidate,
					order
				);

				BridgeScore score = EvaluateBridge(hubRing, candidate);

				if (false) qDebug().nospace()
					<< "  family=" << (reverse ? "reverse" : "forward")
					<< " shift=" << shift
					<< " order=["
					<< order[0] << ","
					<< order[1] << ","
					<< order[2] << ","
					<< order[3] << "]"
					<< " minEdgeDot=" << score.minEdgeDot
					<< " avgEdgeDot=" << score.avgEdgeDot
					<< " maxDiagDiff=" << score.maxDiagDiff
					<< " totalLength=" << score.totalLength;

				if (IsBetterScore(score, bestScore)) {
					bestScore = score;
					bestReverse = reverse;
					CopyRing(bestTubeRing, candidate);

					for (int i = 0; i < 4; ++i) {
						bestOrder[i] = order[i];
					}
				}
			}
		}

		if (false) qDebug() << "  >>> 最佳点序类型 =" << (bestReverse ? "反向" : "正向")
			<< "最佳点序 =" << bestOrder[0] << bestOrder[1] << bestOrder[2] << bestOrder[3]
			<< "最佳最小边方向点积 =" << bestScore.minEdgeDot
			<< "平均边方向点积 =" << bestScore.avgEdgeDot
			<< "最大对角线差 =" << bestScore.maxDiagDiff
			<< "总连接长度 =" << bestScore.totalLength;

		DebugChosenBridge(
			k,
			tubes[k].index,
			tubes[k].isStart,
			hubRing,
			tubeCanonical,
			bestTubeRing,
			bestOrder,
			bestReverse,
			bestScore,
			normalUp
		);

		if (bestScore.minEdgeDot < 0.20)
		{
			if (false) qDebug() << "!!! minEdgeDot =" << bestScore.minEdgeDot;
			if (false) qDebug() << "[桥接警告] 最佳候选仍有偏扭风险";
		}

		// 6.4 
		CageSegment bridgeMain;

		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				bridgeMain.pStart[i][d] = hubRing[i][d];
				bridgeMain.pEnd[i][d] = bestTubeRing[i][d];
			}
		}

		m_allCages.push_back(bridgeMain);
	}

	// ============================================================
	// 7.
	// ============================================================

	HubSDFInfo sdfInfo;
	for (int d = 0; d < 3; ++d) {
		sdfInfo.hubCenter[d] = hubCenter[d];
	}

	sdfInfo.maxRadius = 0.0;

	for (int k = 0; k < N; ++k) {
		SDFTube tube;
		for (int d = 0; d < 3; ++d) {
			tube.center[d] = hubCenter[d];
			tube.axis[d] = tubes[k].dir_out[d];
		}

		if (vtkMath::Norm(tube.axis) < 1e-8) {
			tube.axis[0] = 1.0;
			tube.axis[1] = 0.0;
			tube.axis[2] = 0.0;
		}
		else {
			vtkMath::Normalize(tube.axis);
		}

		tube.radius = std::min(tubes[k].w, tubes[k].h) * 0.5;

		sdfInfo.tubes.push_back(tube);

		sdfInfo.maxRadius = std::max(
			sdfInfo.maxRadius,
			tubes[k].distToHub + std::max(tubes[k].w, tubes[k].h)
		);
	}

	m_sdfHubs.push_back(sdfInfo);

	RebuildGlobalMesh();
}

void CageInteractorStyle::ProjectToSmoothSDF(vtkSmartPointer<vtkPolyData> mesh)
{
	if (m_sdfHubs.empty() || !mesh) return;

	vtkPoints* pts = mesh->GetPoints();
	int numPts = pts->GetNumberOfPoints();
	double k_smooth = 4.0;

	for (int iter = 0; iter < 5; ++iter)
	{
		for (int i = 0; i < numPts; ++i) {
			double p[3];
			pts->GetPoint(i, p);

			int bestHub = -1;
			double minDistToHub = 1e9;
			for (size_t h = 0; h < m_sdfHubs.size(); ++h) {
				double d = sqrt(vtkMath::Distance2BetweenPoints(p, m_sdfHubs[h].hubCenter));
				if (d < minDistToHub) { minDistToHub = d; bestHub = h; }
			}

			// 
			if (bestHub == -1 || minDistToHub > m_sdfHubs[bestHub].maxRadius * 2.5) continue;

			auto& hub = m_sdfHubs[bestHub];

			auto calcF = [&](double pos[3]) {
				double sum = 0.0;
				double coreRadius = hub.coreRadius;
				if (coreRadius <= 1e-8) {
					for (const auto& tube : hub.tubes) {
						coreRadius = std::max(coreRadius, tube.radius * 1.35);
					}
				}

				if (coreRadius > 1e-8) {
					double hubDist = std::sqrt(vtkMath::Distance2BetweenPoints(pos, hub.hubCenter));
					double sphereSdf = hubDist - coreRadius;
					sum += std::exp(-k_smooth * sphereSdf);
				}

				for (const auto& tube : hub.tubes) {
					double ap[3] = { pos[0] - tube.center[0], pos[1] - tube.center[1], pos[2] - tube.center[2] };
					double dot = vtkMath::Dot(ap, tube.axis);
					// 
					if (dot < 0) dot = 0;

					double proj[3] = { tube.center[0] + dot * tube.axis[0], tube.center[1] + dot * tube.axis[1], tube.center[2] + dot * tube.axis[2] };
					double dist = sqrt(vtkMath::Distance2BetweenPoints(pos, proj));

					double sdf = dist - tube.radius;
					sum += std::exp(-k_smooth * sdf);
				}
				return -(1.0 / k_smooth) * std::log(sum);
			};

			double F_val = calcF(p);

			double grad[3], eps = 1e-4;
			for (int d = 0; d < 3; d++) {
				double pp[3] = { p[0], p[1], p[2] }; pp[d] += eps;
				grad[d] = (calcF(pp) - F_val) / eps;
			}

			double gradSq = vtkMath::Dot(grad, grad);
			if (gradSq > 1e-8) {
				double weight = std::exp(-minDistToHub / (hub.maxRadius * 1.5));

				// 
				double maxDisp = hub.maxRadius * 0.05;
				for (int d = 0; d < 3; d++) {
					double disp = weight * F_val * (grad[d] / gradSq) * 0.3; // 
					disp = std::max(-maxDisp, std::min(maxDisp, disp));
					p[d] -= disp;
				}
			}
			pts->SetPoint(i, p);
		}
	}
	pts->Modified();
}
void CageInteractorStyle::ComputeRMF(const double x0[3], const double t0[3], const double u0[3],
	const double x1[3], const double t1[3], double u1[3], double v1[3])
{
	double v_1[3] = { x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2] };
	double c1 = vtkMath::Dot(v_1, v_1);

	// 
	if (c1 < 1e-8) {
		double dot = vtkMath::Dot(u0, t1);
		for (int i = 0; i < 3; i++) u1[i] = u0[i] - dot * t1[i];
		vtkMath::Normalize(u1);
		vtkMath::Cross(t1, u1, v1);
		return;
	}

	double rL[3], tL[3];
	double dot_u0_v1 = vtkMath::Dot(u0, v_1);
	double dot_t0_v1 = vtkMath::Dot(t0, v_1);

	for (int i = 0; i < 3; i++) {
		rL[i] = u0[i] - (2.0 / c1) * dot_u0_v1 * v_1[i];
		tL[i] = t0[i] - (2.0 / c1) * dot_t0_v1 * v_1[i];
	}

	double v_2[3] = { t1[0] - tL[0], t1[1] - tL[1], t1[2] - tL[2] };
	double c2 = vtkMath::Dot(v_2, v_2);

	if (c2 < 1e-8) {
		for (int i = 0; i < 3; i++) u1[i] = rL[i];
		vtkMath::Normalize(u1);
		vtkMath::Cross(t1, u1, v1);
		return;
	}

	double dot_rL_v2 = vtkMath::Dot(rL, v_2);
	for (int i = 0; i < 3; i++) {
		u1[i] = rL[i] - (2.0 / c2) * dot_rL_v2 * v_2[i];
	}
	vtkMath::Normalize(u1);
	vtkMath::Cross(t1, u1, v1);
}


bool CageInteractorStyle::GenerateSingleBranchLoopTopology(vtkRenderer* ren)
{
	int selectedCount = static_cast<int>(m_selectedFacesForBranch.size());

	if (selectedCount != 1 && selectedCount != 2) {
		if (false) qDebug() << "[圆环生成] 请选择 1 个或 2 个控制笼端面。";
		return false;
	}

	if (m_skelNodes.empty() || m_skelEdges.empty()) {
		if (false) qDebug() << "[圆环生成] 骨架数据为空。";
		return false;
	}

	std::vector<LoopAnchorInfo> anchors;
	std::vector<double> inputWidths;
	std::vector<double> inputHeights;
	std::vector<std::array<std::array<double, 3>, 4>> inputRings;

	for (const SelectedFace& sel : m_selectedFacesForBranch) {
		if (sel.cageIndex < 0 || sel.cageIndex >= static_cast<int>(m_allCages.size())) {
			if (false) qDebug() << "[圆环生成] 选中的控制笼端面无效。";
			return false;
		}

		LoopAnchorInfo anchor;
		anchor.node = FindNearestSkeletonNodeToSelectedFace(sel, anchor.center, anchor.normal);

		if (anchor.node < 0) {
			if (false) qDebug() << "[圆环生成] 选中端面绑定到骨架失败。";
			return false;
		}

		double capRing[4][3];
		if (!GetCanonicalCageCapRing(m_allCages, sel.cageIndex, sel.isStart, capRing)) {
			if (false) qDebug() << "[圆环生成] 读取选中端面环失败。";
			return false;
		}

		std::array<std::array<double, 3>, 4> storedRing;
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				storedRing[i][d] = capRing[i][d];
			}
		}

		inputRings.push_back(storedRing);

		double w = std::sqrt(vtkMath::Distance2BetweenPoints(capRing[0], capRing[1]));
		double h = std::sqrt(vtkMath::Distance2BetweenPoints(capRing[0], capRing[3]));

		if (w < 1e-8 || h < 1e-8) {
			if (false) qDebug() << "[圆环生成] 选中端面退化，无法生成。";
			return false;
		}

		inputWidths.push_back(w);
		inputHeights.push_back(h);
		anchors.push_back(anchor);
	}

	const LoopAnchorInfo& primaryAnchor = anchors.front();

	std::vector<std::vector<int>> adjList = BuildSkeletonAdjacency();

	bool isMultiBranchLoop = false;
	std::vector<std::vector<int>> connectorPaths;
	std::vector<int> cycleNodes;
	int attachNode = -1;

	auto nodePoint = [&](int node, double p[3]) {
		p[0] = m_skelNodes[3 * node];
		p[1] = m_skelNodes[3 * node + 1];
		p[2] = m_skelNodes[3 * node + 2];
	};

	auto findPathAvoidingNode = [&](int startNode, int endNode, int blockedNode) {
		std::vector<int> emptyPath;

		if (startNode < 0 || endNode < 0 ||
			startNode >= static_cast<int>(adjList.size()) ||
			endNode >= static_cast<int>(adjList.size()) ||
			startNode == blockedNode || endNode == blockedNode) {
			return emptyPath;
		}

		std::vector<int> parent(adjList.size(), -1);
		std::vector<bool> visited(adjList.size(), false);
		std::queue<int> q;

		visited[startNode] = true;
		q.push(startNode);

		while (!q.empty()) {
			int curr = q.front();
			q.pop();

			if (curr == endNode) {
				break;
			}

			for (int next : adjList[curr]) {
				if (next == blockedNode || visited[next]) {
					continue;
				}

				visited[next] = true;
				parent[next] = curr;
				q.push(next);
			}
		}

		if (!visited[endNode]) {
			return emptyPath;
		}

		std::vector<int> path;
		for (int curr = endNode; curr != -1; curr = parent[curr]) {
			path.push_back(curr);
		}

		std::reverse(path.begin(), path.end());
		return path;
	};

	if (m_useForcedLoopTopology) {
		cycleNodes = m_forcedLoopCycleNodes;
		connectorPaths = m_forcedLoopConnectorPaths;
		isMultiBranchLoop = (selectedCount == 2);

		if (cycleNodes.empty() || connectorPaths.empty() || connectorPaths[0].empty()) {
			if (false) qDebug() << "[圆环生成] 强制使用的圆环轨迹为空。";
			return false;
		}

		attachNode = connectorPaths[0].back();
		if (false) qDebug() << "[圆环生成] 使用自动识别得到的圆环轨迹，选中端面数量 =" << selectedCount
			<< "圆环骨架点数量 =" << static_cast<int>(cycleNodes.size())
			<< "连接路径数量 =" << static_cast<int>(connectorPaths.size())
			<< "连接交点 =" << attachNode;
	}
	else if (selectedCount == 1) {
		cycleNodes = FindBestSkeletonCycle(
			anchors,
			adjList,
			connectorPaths,
			isMultiBranchLoop
		);

		if (cycleNodes.empty() || connectorPaths.empty() || connectorPaths[0].empty()) {
			if (false) qDebug() << "[圆环生成] 未找到端面前方的闭合圆环。";
			return false;
		}

		if (isMultiBranchLoop) {
			ShowSkeletonCycle(ren, cycleNodes, connectorPaths, true);
			if (false) qDebug() << "[圆环生成] 检测到的圆环更像多分支圆环，此处只预览不生成。";
			return false;
		}

		attachNode = connectorPaths[0].back();
	}
	else {
		std::vector<int> inputPath = FindPathBFS(anchors[0].node, anchors[1].node, adjList);

		if (inputPath.size() < 3) {
			if (false) qDebug() << "[双端面圆环] 两个选中端面之间找不到骨架路径。";
			return false;
		}

		int junctionIndexInPath = -1;
		int bestDegree = -1;

		for (int i = 1; i + 1 < static_cast<int>(inputPath.size()); ++i) {
			int node = inputPath[i];
			int degree = static_cast<int>(adjList[node].size());

			if (degree >= 4 && degree > bestDegree) {
				bestDegree = degree;
				junctionIndexInPath = i;
			}
		}

		if (junctionIndexInPath < 0) {
			for (int i = 1; i + 1 < static_cast<int>(inputPath.size()); ++i) {
				int node = inputPath[i];
				int degree = static_cast<int>(adjList[node].size());

				if (degree >= 3 && degree > bestDegree) {
					bestDegree = degree;
					junctionIndexInPath = i;
				}
			}
		}

		if (junctionIndexInPath < 0) {
			if (false) qDebug() << "[双端面圆环] 两个选中端面之间没有找到交点。";
			return false;
		}

		attachNode = inputPath[junctionIndexInPath];

		std::set<int> inputNeighbors;
		inputNeighbors.insert(inputPath[junctionIndexInPath - 1]);
		inputNeighbors.insert(inputPath[junctionIndexInPath + 1]);

		std::vector<int> loopNeighbors;
		for (int n : adjList[attachNode]) {
			if (inputNeighbors.find(n) == inputNeighbors.end()) {
				loopNeighbors.push_back(n);
			}
		}

		std::vector<int> bestLoopPath;
		int bestA = -1;
		int bestB = -1;

		for (size_t i = 0; i < loopNeighbors.size(); ++i) {
			for (size_t j = i + 1; j < loopNeighbors.size(); ++j) {
				std::vector<int> candidatePath =
					findPathAvoidingNode(loopNeighbors[i], loopNeighbors[j], attachNode);

				if (candidatePath.empty()) {
					continue;
				}

				if (bestLoopPath.empty() || candidatePath.size() < bestLoopPath.size()) {
					bestLoopPath = candidatePath;
					bestA = loopNeighbors[i];
					bestB = loopNeighbors[j];
				}
			}
		}

		if (bestLoopPath.empty()) {
			if (false) qDebug() << "[双端面圆环] 无法从两个非输入分支追踪出圆环。";
			return false;
		}

		cycleNodes.clear();
		cycleNodes.push_back(attachNode);
		cycleNodes.insert(cycleNodes.end(), bestLoopPath.begin(), bestLoopPath.end());

		connectorPaths.clear();
		connectorPaths.push_back(FindPathBFS(anchors[0].node, attachNode, adjList));
		connectorPaths.push_back(FindPathBFS(anchors[1].node, attachNode, adjList));
		isMultiBranchLoop = true;

		if (false) qDebug() << "[双端面圆环] 已追踪到交点圆环:"
			<< "交点 =" << attachNode
			<< "输入邻接点 =" << *inputNeighbors.begin() << *inputNeighbors.rbegin()
			<< "圆环邻接点 =" << bestA << bestB
			<< "圆环骨架点数量 =" << static_cast<int>(cycleNodes.size());
	}

	auto attachIt = std::find(cycleNodes.begin(), cycleNodes.end(), attachNode);

	if (attachNode < 0 || attachIt == cycleNodes.end()) {
		if (false) qDebug() << "[圆环生成] 连接路径没有接触到检测出的圆环。";
		return false;
	}

	auto normalizeOrFallback = [](double v[3], const double fallback[3]) {
		if (vtkMath::Norm(v) < 1e-8) {
			v[0] = fallback[0];
			v[1] = fallback[1];
			v[2] = fallback[2];
		}

		if (vtkMath::Norm(v) < 1e-8) {
			v[0] = 1.0;
			v[1] = 0.0;
			v[2] = 0.0;
		}

		vtkMath::Normalize(v);
	};

	double baseWidth = inputWidths[0];
	double baseHeight = inputHeights[0];

	for (size_t i = 1; i < inputWidths.size(); ++i) {
		baseWidth = 0.5 * (baseWidth + inputWidths[i]);
		baseHeight = 0.5 * (baseHeight + inputHeights[i]);
	}

	// ============================================================
	// 1.
	// ============================================================

	std::vector<int> orderedCycle;

	int attachIdx = static_cast<int>(std::distance(cycleNodes.begin(), attachIt));

	for (int i = 0; i < static_cast<int>(cycleNodes.size()); ++i) {
		int id = cycleNodes[(attachIdx + i) % static_cast<int>(cycleNodes.size())];
		orderedCycle.push_back(id);
	}

	orderedCycle.push_back(attachNode);

	std::vector<double> cumulative(orderedCycle.size(), 0.0);

	for (size_t i = 1; i < orderedCycle.size(); ++i) {
		double p0[3], p1[3];
		nodePoint(orderedCycle[i - 1], p0);
		nodePoint(orderedCycle[i], p1);

		cumulative[i] =
			cumulative[i - 1] +
			std::sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
	}

	double loopLength = cumulative.back();

	if (loopLength < 1e-8) {
		if (false) qDebug() << "[圆环生成] 检测到的圆环长度为零。";
		return false;
	}

	auto sampleCyclePoint = [&](double targetLen, double out[3]) {
		while (targetLen < 0.0) {
			targetLen += loopLength;
		}

		while (targetLen >= loopLength) {
			targetLen -= loopLength;
		}

		if (targetLen <= 1e-8) {
			nodePoint(orderedCycle.front(), out);
			return;
		}

		for (size_t i = 1; i < orderedCycle.size(); ++i) {
			if (cumulative[i] < targetLen) {
				continue;
			}

			double segLen = cumulative[i] - cumulative[i - 1];

			double alpha = 0.0;

			if (segLen > 1e-8) {
				alpha = (targetLen - cumulative[i - 1]) / segLen;
			}

			double p0[3], p1[3];
			nodePoint(orderedCycle[i - 1], p0);
			nodePoint(orderedCycle[i], p1);

			for (int d = 0; d < 3; ++d) {
				out[d] = p0[d] * (1.0 - alpha) + p1[d] * alpha;
			}

			return;
		}

		nodePoint(orderedCycle.front(), out);
	};

	// ============================================================
	// 2. 
	// ============================================================

	double loopCenter[3] = { 0.0, 0.0, 0.0 };

	for (int node : cycleNodes) {
		double p[3];
		nodePoint(node, p);

		for (int d = 0; d < 3; ++d) {
			loopCenter[d] += p[d] / static_cast<double>(cycleNodes.size());
		}
	}

	double loopNormal[3] = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < static_cast<int>(cycleNodes.size()); ++i) {
		int j = (i + 1) % static_cast<int>(cycleNodes.size());

		double pi[3], pj[3];
		nodePoint(cycleNodes[i], pi);
		nodePoint(cycleNodes[j], pj);

		double a[3] = {
			pi[0] - loopCenter[0],
			pi[1] - loopCenter[1],
			pi[2] - loopCenter[2]
		};

		double b[3] = {
			pj[0] - loopCenter[0],
			pj[1] - loopCenter[1],
			pj[2] - loopCenter[2]
		};

		double c[3];
		vtkMath::Cross(a, b, c);

		for (int d = 0; d < 3; ++d) {
			loopNormal[d] += c[d];
		}
	}

	normalizeOrFallback(loopNormal, primaryAnchor.normal);

	if (vtkMath::Dot(loopNormal, primaryAnchor.normal) < 0.0) {
		for (int d = 0; d < 3; ++d) {
			loopNormal[d] = -loopNormal[d];
		}
	}

	auto computeRadial = [&](const double p[3], double radial[3]) {
		for (int d = 0; d < 3; ++d) {
			radial[d] = p[d] - loopCenter[d];
		}

		double dotN = vtkMath::Dot(radial, loopNormal);

		for (int d = 0; d < 3; ++d) {
			radial[d] -= dotN * loopNormal[d];
		}

		normalizeOrFallback(radial, primaryAnchor.normal);
	};

	auto radialDistanceToLoopCenter = [&](const double p[3]) -> double {
		double v[3] = {
			p[0] - loopCenter[0],
			p[1] - loopCenter[1],
			p[2] - loopCenter[2]
		};

		double dotN = vtkMath::Dot(v, loopNormal);

		for (int d = 0; d < 3; ++d) {
			v[d] -= dotN * loopNormal[d];
		}

		return vtkMath::Norm(v);
	};

	std::vector<double> loopRadii;

	for (int node : cycleNodes) {
		double p[3];
		nodePoint(node, p);

		double r = radialDistanceToLoopCenter(p);

		if (r > 1e-8) {
			loopRadii.push_back(r);
		}
	}

	if (loopRadii.empty()) {
		if (false) qDebug() << "[圆环生成] 估算圆环半径失败。";
		return false;
	}

	std::sort(loopRadii.begin(), loopRadii.end());

	double loopRadius = loopRadii[loopRadii.size() / 2];

	double minRingRadialWidth = std::max(baseWidth * 0.6, loopRadius * 0.12);

	// ============================================================
	// 3. 

	struct RayHitItem {
		double t;
		double p[3];
	};

	auto collectIntersectionsNoFail = [&](const double p0[3], const double p1[3], std::vector<RayHitItem>& hits) -> bool {
		hits.clear();

		if (!m_RefLocator) {
			return false;
		}

		vtkNew<vtkPoints> hitPts;
		vtkNew<vtkIdList> cellIds;

		int ok = m_RefLocator->IntersectWithLine(
			p0,
			p1,
			hitPts.GetPointer(),
			cellIds.GetPointer()
		);

		if (!ok || hitPts->GetNumberOfPoints() <= 0) {
			return false;
		}

		double dir[3] = {
			p1[0] - p0[0],
			p1[1] - p0[1],
			p1[2] - p0[2]
		};

		normalizeOrFallback(dir, primaryAnchor.normal);

		for (vtkIdType i = 0; i < hitPts->GetNumberOfPoints(); ++i) {
			double hp[3];
			hitPts->GetPoint(i, hp);

			double v[3] = {
				hp[0] - p0[0],
				hp[1] - p0[1],
				hp[2] - p0[2]
			};

			double t = vtkMath::Dot(v, dir);

			if (t < 1e-6) {
				continue;
			}

			RayHitItem item;
			item.t = t;
			item.p[0] = hp[0];
			item.p[1] = hp[1];
			item.p[2] = hp[2];

			hits.push_back(item);
		}

		if (hits.empty()) {
			return false;
		}

		std::sort(
			hits.begin(),
			hits.end(),
			[](const RayHitItem& a, const RayHitItem& b) {
			return a.t < b.t;
		}
		);

		// 
		std::vector<RayHitItem> uniqueHits;

		for (const auto& h : hits) {
			if (!uniqueHits.empty()) {
				double d2 = vtkMath::Distance2BetweenPoints(uniqueHits.back().p, h.p);

				if (d2 < 1e-6) {
					continue;
				}
			}

			uniqueHits.push_back(h);
		}

		hits.swap(uniqueHits);

		return !hits.empty();
	};

	struct SectionInfo {
		double inner[3];
		double outer[3];
		double center[3];
		double radial[3];
		double innerRadius;
		double outerRadius;
		double width;
		bool valid;
	};

	auto sampleBackgroundSection = [&](const double skelPoint[3], SectionInfo& sec) -> bool {
		sec.valid = false;
		sec.width = 0.0;
		sec.innerRadius = 0.0;
		sec.outerRadius = 0.0;

		double radial[3];
		computeRadial(skelPoint, radial);

		for (int d = 0; d < 3; ++d) {
			sec.radial[d] = radial[d];
		}

		double rayLength = std::max(loopLength * 2.0, baseWidth * 30.0);
		rayLength = std::max(rayLength, 80.0);

		double rayStart[3];
		double rayEnd[3];

		for (int d = 0; d < 3; ++d) {
			rayStart[d] = loopCenter[d] - radial[d] * rayLength * 0.05;
			rayEnd[d] = loopCenter[d] + radial[d] * rayLength;
		}

		std::vector<RayHitItem> hits;

		if (!collectIntersectionsNoFail(rayStart, rayEnd, hits)) {
			return false;
		}

		std::vector<RayHitItem> positiveHits;

		for (const auto& h : hits) {
			double v[3] = {
				h.p[0] - loopCenter[0],
				h.p[1] - loopCenter[1],
				h.p[2] - loopCenter[2]
			};

			double t = vtkMath::Dot(v, radial);

			if (t > 1e-5) {
				positiveHits.push_back(h);
			}
		}

		if (positiveHits.size() < 2) {
			return false;
		}

		std::sort(
			positiveHits.begin(),
			positiveHits.end(),
			[&](const RayHitItem& a, const RayHitItem& b) {
			double va[3] = {
				a.p[0] - loopCenter[0],
				a.p[1] - loopCenter[1],
				a.p[2] - loopCenter[2]
			};
			double vb[3] = {
				b.p[0] - loopCenter[0],
				b.p[1] - loopCenter[1],
				b.p[2] - loopCenter[2]
			};
			return vtkMath::Dot(va, radial) < vtkMath::Dot(vb, radial);
		}
		);

		// 以骨架圆环半径作为中心线，分别寻找中心线内外最近的背景边界。
		double centerGap = std::max(baseWidth * 0.05, loopRadius * 0.01);
		const RayHitItem* innerHit = nullptr;
		const RayHitItem* outerHit = nullptr;

		for (const auto& h : positiveHits) {
			double vh[3] = {
				h.p[0] - loopCenter[0],
				h.p[1] - loopCenter[1],
				h.p[2] - loopCenter[2]
			};
			double r = vtkMath::Dot(vh, radial);

			if (r < loopRadius - centerGap) {
				innerHit = &h;
			}
			else if (r > loopRadius + centerGap) {
				outerHit = &h;
				break;
			}
		}

		if (!innerHit || !outerHit) {
			return false;
		}

		double vi[3] = {
			innerHit->p[0] - loopCenter[0],
			innerHit->p[1] - loopCenter[1],
			innerHit->p[2] - loopCenter[2]
		};
		double vo[3] = {
			outerHit->p[0] - loopCenter[0],
			outerHit->p[1] - loopCenter[1],
			outerHit->p[2] - loopCenter[2]
		};

		sec.innerRadius = vtkMath::Dot(vi, radial);
		sec.outerRadius = vtkMath::Dot(vo, radial);
		sec.width = sec.outerRadius - sec.innerRadius;

		if (sec.innerRadius < 1e-5 || sec.width < minRingRadialWidth) {
			return false;
		}

		for (int d = 0; d < 3; ++d) {
			sec.inner[d] = innerHit->p[d];
			sec.outer[d] = outerHit->p[d];
			sec.center[d] = (sec.inner[d] + sec.outer[d]) * 0.5;
		}

		sec.valid = true;
		return true;
	};

	// ============================================================
	// 4. 
	// ============================================================

	const int segmentCount = 6;

	double firstBoundaryLen =
		(selectedCount == 1)
		? -loopLength / static_cast<double>(segmentCount * 2)
		: 0.0;

	struct RingFrame {
		double p[4][3];
		double center[3];
		double inner[3];
		double outer[3];
		double radial[3];
		double up[3];
		double width;
	};

	std::vector<RingFrame> rings(segmentCount + 1);

	std::vector<double> validWidths;
	std::vector<double> validInnerRadii;
	std::vector<double> validOuterRadii;

	const int estimateSectionCount = 48;

	for (int i = 0; i < estimateSectionCount; ++i) {
		double sampleLen =
			firstBoundaryLen +
			static_cast<double>(i) *
			loopLength / static_cast<double>(estimateSectionCount);

		double skelP[3];
		sampleCyclePoint(sampleLen, skelP);

		SectionInfo sec;

		if (sampleBackgroundSection(skelP, sec)) {
			validWidths.push_back(sec.width);
			validInnerRadii.push_back(sec.innerRadius);
			validOuterRadii.push_back(sec.outerRadius);
		}
	}

	double fallbackWidth = loopRadius * 0.42;
	double fallbackInnerRadius = std::max(loopRadius - fallbackWidth * 0.5, baseWidth * 0.5);
	double fallbackOuterRadius = loopRadius + fallbackWidth * 0.5;

	if (!validWidths.empty()) {
		std::sort(validWidths.begin(), validWidths.end());
		fallbackWidth = validWidths[validWidths.size() / 2];
	}

	if (!validInnerRadii.empty()) {
		std::sort(validInnerRadii.begin(), validInnerRadii.end());
		fallbackInnerRadius = validInnerRadii[validInnerRadii.size() / 2];
	}

	if (!validOuterRadii.empty()) {
		std::sort(validOuterRadii.begin(), validOuterRadii.end());
		size_t outerIndex = static_cast<size_t>(
			static_cast<double>(validOuterRadii.size() - 1) * 0.85
			);
		fallbackOuterRadius = validOuterRadii[outerIndex];
	}

	if (fallbackOuterRadius - fallbackInnerRadius < minRingRadialWidth) {
		double midRadius = (fallbackOuterRadius + fallbackInnerRadius) * 0.5;
		fallbackInnerRadius = std::max(midRadius - minRingRadialWidth * 0.5, baseWidth * 0.5);
		fallbackOuterRadius = fallbackInnerRadius + minRingRadialWidth;
	}

	double radiusPadding = std::max(baseWidth * 0.08, (fallbackOuterRadius - fallbackInnerRadius) * 0.04);
	fallbackInnerRadius = std::max(fallbackInnerRadius - radiusPadding * 0.35, baseWidth * 0.3);
	fallbackOuterRadius += radiusPadding;

	// 外圈采样点在背景边界上时，相邻外点连成的边会缩在背景内侧。
	// 按正多边形外接半径外扩，使外侧边线近似与背景最外侧相切。
	const double piValue = 3.14159265358979323846;
	double outerTangentScale = 1.0 / std::cos(piValue / static_cast<double>(segmentCount));
	if (selectedCount == 2) {
		outerTangentScale *= 1.03;
	}
	fallbackOuterRadius *= outerTangentScale;
	fallbackWidth = fallbackOuterRadius - fallbackInnerRadius;

	for (int i = 0; i < segmentCount; ++i) {
		RingFrame& rf = rings[i];

		double skelP[3];
		double sampleLen =
			firstBoundaryLen +
			static_cast<double>(i) *
			loopLength / static_cast<double>(segmentCount);

		sampleCyclePoint(sampleLen, skelP);

		double radial[3];
		computeRadial(skelP, radial);

		double innerRadius = fallbackInnerRadius;
		double outerRadius = fallbackOuterRadius;

		if (outerRadius - innerRadius < minRingRadialWidth) {
			outerRadius = innerRadius + minRingRadialWidth;
		}

		for (int d = 0; d < 3; ++d) {
			rf.radial[d] = radial[d];
			rf.inner[d] = loopCenter[d] + radial[d] * innerRadius;
			rf.outer[d] = loopCenter[d] + radial[d] * outerRadius;
			rf.center[d] = (rf.inner[d] + rf.outer[d]) * 0.5;
			rf.up[d] = loopNormal[d];
		}

		rf.width = outerRadius - innerRadius;
	}

	double ringHeight = baseHeight;

	// ============================================================
	// 5. 
	// ============================================================

	auto makeRing = [&](RingFrame& rf) {
		for (int d = 0; d < 3; ++d) {
			rf.p[0][d] = rf.inner[d] - rf.up[d] * ringHeight * 0.5;
			rf.p[1][d] = rf.outer[d] - rf.up[d] * ringHeight * 0.5;
			rf.p[2][d] = rf.outer[d] + rf.up[d] * ringHeight * 0.5;
			rf.p[3][d] = rf.inner[d] + rf.up[d] * ringHeight * 0.5;
		}
	};

	for (int i = 0; i < segmentCount; ++i) {
		makeRing(rings[i]);
	}

	// 
	for (int k = 0; k < 4; ++k) {
		for (int d = 0; d < 3; ++d) {
			rings[segmentCount].p[k][d] = rings[0].p[k][d];
		}
	}

	for (int d = 0; d < 3; ++d) {
		rings[segmentCount].center[d] = rings[0].center[d];
		rings[segmentCount].inner[d] = rings[0].inner[d];
		rings[segmentCount].outer[d] = rings[0].outer[d];
		rings[segmentCount].radial[d] = rings[0].radial[d];
		rings[segmentCount].up[d] = rings[0].up[d];
	}

	rings[segmentCount].width = rings[0].width;

	// ============================================================
	// 6. 
	// ============================================================

	int addedSegments = 0;
	int loopStartCageIndex = static_cast<int>(m_allCages.size());

	for (int i = 0; i < segmentCount; ++i) {
		CageSegment seg;

		CopyRing4(seg.pStart, rings[i].p);
		CopyRing4(seg.pEnd, rings[i + 1].p);

		m_allCages.push_back(seg);
		++addedSegments;
	}

	// ============================================================
	// 7. 

	auto bridgeDistance = [](const double a[3], const double b[3]) -> double {
		return std::sqrt(vtkMath::Distance2BetweenPoints(a, b));
	};

	auto bridgeEdgeDot = [](const double a0[3], const double a1[3],
		const double b0[3], const double b1[3]) -> double {
		double va[3] = {
			a1[0] - a0[0],
			a1[1] - a0[1],
			a1[2] - a0[2]
		};
		double vb[3] = {
			b1[0] - b0[0],
			b1[1] - b0[1],
			b1[2] - b0[2]
		};

		double na = vtkMath::Norm(va);
		double nb = vtkMath::Norm(vb);

		if (na < 1e-10 || nb < 1e-10) {
			return 0.0;
		}

		for (int d = 0; d < 3; ++d) {
			va[d] /= na;
			vb[d] /= nb;
		}

		return vtkMath::Dot(va, vb);
	};

	struct LoopBridgeScore {
		double minEdgeDot;
		double avgEdgeDot;
		double totalLength;
		double maxDiagDiff;
	};

	auto evaluateLoopBridge = [&](const double startRing[4][3],
		const double endRing[4][3]) -> LoopBridgeScore {
		LoopBridgeScore score;
		score.minEdgeDot = 1.0;
		score.avgEdgeDot = 0.0;
		score.totalLength = 0.0;
		score.maxDiagDiff = 0.0;

		for (int i = 0; i < 4; ++i) {
			int j = (i + 1) % 4;

			double dot = bridgeEdgeDot(startRing[i], startRing[j], endRing[i], endRing[j]);
			score.minEdgeDot = std::min(score.minEdgeDot, dot);
			score.avgEdgeDot += dot * 0.25;
			score.totalLength += bridgeDistance(startRing[i], endRing[i]);

			double diagA = bridgeDistance(startRing[i], endRing[j]);
			double diagB = bridgeDistance(startRing[j], endRing[i]);
			score.maxDiagDiff = std::max(score.maxDiagDiff, std::abs(diagA - diagB));
		}

		return score;
	};

	auto isBetterLoopBridge = [](const LoopBridgeScore& a, const LoopBridgeScore& b) -> bool {
		if (a.minEdgeDot > b.minEdgeDot + 1e-6) return true;
		if (a.minEdgeDot < b.minEdgeDot - 1e-6) return false;

		if (a.avgEdgeDot > b.avgEdgeDot + 1e-6) return true;
		if (a.avgEdgeDot < b.avgEdgeDot - 1e-6) return false;

		if (a.maxDiagDiff < b.maxDiagDiff - 1e-6) return true;
		if (a.maxDiagDiff > b.maxDiagDiff + 1e-6) return false;

		return a.totalLength < b.totalLength;
	};

	auto buildLoopTubeCandidate = [](const double tubeRing[4][3],
		bool reverse,
		int shift,
		double outRing[4][3],
		int outOrder[4]) {
		static const int forwardBase[4] = { 0, 1, 2, 3 };
		static const int reverseBase[4] = { 0, 3, 2, 1 };

		const int* base = reverse ? reverseBase : forwardBase;

		for (int i = 0; i < 4; ++i) {
			int id = base[(i + shift) % 4];
			outOrder[i] = id;

			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = tubeRing[id][d];
			}
		}
	};

	struct LoopBridgeCandidate {
		int inputIndex;
		int segmentIndex;
		bool reverse;
		int order[4];
		double face[4][3];
		double inputRing[4][3];
		LoopBridgeScore score;
	};

	auto makeLoopOuterFace = [&](int segmentIndex, double outFace[4][3]) {
		int nextIndex = segmentIndex + 1;

		for (int d = 0; d < 3; ++d) {
			outFace[0][d] = rings[segmentIndex].p[1][d];
			outFace[1][d] = rings[segmentIndex].p[2][d];
			outFace[2][d] = rings[nextIndex].p[2][d];
			outFace[3][d] = rings[nextIndex].p[1][d];
		}
	};

	auto copyInputRing = [&](int inputIndex, double outRing[4][3]) {
		for (int i = 0; i < 4; ++i) {
			for (int d = 0; d < 3; ++d) {
				outRing[i][d] = inputRings[inputIndex][i][d];
			}
		}
	};

	auto buildBestLoopBridge = [&](int inputIndex, int segmentIndex, LoopBridgeCandidate& outCandidate) {
		double inputRing[4][3];
		copyInputRing(inputIndex, inputRing);

		makeLoopOuterFace(segmentIndex, outCandidate.face);

		outCandidate.inputIndex = inputIndex;
		outCandidate.segmentIndex = segmentIndex;
		outCandidate.reverse = false;
		outCandidate.score.minEdgeDot = -1e100;
		outCandidate.score.avgEdgeDot = -1e100;
		outCandidate.score.totalLength = 1e100;
		outCandidate.score.maxDiagDiff = 1e100;

		for (int family = 0; family < 2; ++family) {
			bool reverse = (family == 1);

			for (int shift = 0; shift < 4; ++shift) {
				double candidate[4][3];
				int order[4];

				buildLoopTubeCandidate(inputRing, reverse, shift, candidate, order);

				LoopBridgeScore score = evaluateLoopBridge(outCandidate.face, candidate);

				if (isBetterLoopBridge(score, outCandidate.score)) {
					outCandidate.score = score;
					outCandidate.reverse = reverse;
					CopyRing4(outCandidate.inputRing, candidate);

					for (int i = 0; i < 4; ++i) {
						outCandidate.order[i] = order[i];
					}
				}
			}
		}
	};

	auto bridgeAssignmentScore = [](const LoopBridgeCandidate& a, const LoopBridgeCandidate& b) -> double {
		return a.score.totalLength + b.score.totalLength
			+ a.score.maxDiagDiff * 8.0 + b.score.maxDiagDiff * 8.0
			- a.score.minEdgeDot * 20.0 - b.score.minEdgeDot * 20.0;
	};

	auto addLoopBridge = [&](const LoopBridgeCandidate& candidate) {
		CageSegment bridgeToInput;
		CopyRing4(bridgeToInput.pStart, candidate.face);
		CopyRing4(bridgeToInput.pEnd, candidate.inputRing);
		m_allCages.push_back(bridgeToInput);

		m_suppressedCageSideFaces.push_back({
			loopStartCageIndex + candidate.segmentIndex,
			1
			});

		if (false) qDebug() << "[圆环生成] 已添加圆环外侧面到输入端面的桥接段:"
			<< "输入端序号 =" << candidate.inputIndex
			<< "圆环段序号 =" << candidate.segmentIndex
			<< "点序 =" << candidate.order[0] << candidate.order[1] << candidate.order[2] << candidate.order[3]
			<< "是否反向 =" << candidate.reverse
			<< "最小边方向点积 =" << candidate.score.minEdgeDot
			<< "平均边方向点积 =" << candidate.score.avgEdgeDot
			<< "最大对角线差 =" << candidate.score.maxDiagDiff
			<< "总连接长度 =" << candidate.score.totalLength;
	};

	if (selectedCount == 1) {
		LoopBridgeCandidate candidate;
		buildBestLoopBridge(0, 0, candidate);
		addLoopBridge(candidate);
	}
	else {
		LoopBridgeCandidate a0, a1, b0, b1;
		buildBestLoopBridge(0, 0, a0);
		buildBestLoopBridge(1, segmentCount - 1, a1);
		buildBestLoopBridge(0, segmentCount - 1, b0);
		buildBestLoopBridge(1, 0, b1);

		double scoreA = bridgeAssignmentScore(a0, a1);
		double scoreB = bridgeAssignmentScore(b0, b1);

		if (scoreA <= scoreB) {
			addLoopBridge(a0);
			addLoopBridge(a1);
		}
		else {
			addLoopBridge(b0);
			addLoopBridge(b1);
		}
	}

	if (addedSegments == 0) {
		if (false) qDebug() << "[圆环生成] 没有生成任何控制笼段。";
		return false;
	}

	if (!m_useForcedLoopTopology) {
		ShowSkeletonCycle(ren, cycleNodes, connectorPaths, false);
	}

	RebuildGlobalMesh();

	if (false) qDebug() << "[圆环生成] 已生成闭合 6 段圆环控制笼:"
		<< "圆环段数 =" << addedSegments
		<< "连接交点 =" << attachNode
		<< "圆环长度 =" << loopLength
		<< "圆环半径 =" << loopRadius
		<< "基础宽度 =" << baseWidth
		<< "备用宽度 =" << fallbackWidth
		<< "最小径向宽度 =" << minRingRadialWidth
		<< "备用外半径 =" << fallbackOuterRadius
		<< "外扩比例 =" << outerTangentScale
		<< "基础高度 =" << baseHeight
		<< "有效背景截面数量 =" << static_cast<int>(validWidths.size());

	return true;
}

void CageInteractorStyle::RecognizeAutoTopology(vtkRenderer* ren)
{
	if (!ren) {
		if (false) qDebug() << "[路线识别] renderer 为空。";
		return;
	}

	// false：只识别并显示轨迹，不生成 hub / bridge / loop
	AutoDetectAndShowTopology(ren, false);
}
void CageInteractorStyle::BuildAutoTopology(vtkRenderer* ren)
{
	if (!ren) {
		if (false) qDebug() << "[自动构建] renderer 为空。";
		return;
	}

	// true：识别并直接生成 hub / bridge / loop
	AutoDetectAndShowTopology(ren, true);
}
