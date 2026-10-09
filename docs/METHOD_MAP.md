# 论文方法与实现对应关系

| 论文内容 | 源文件 | 关键函数 |
| --- | --- | --- |
| 第3节 平均曲率流骨架提取 | SkeletonExtractor.cpp | extractSkeleton |
| 第3节 非设计域骨架剔除 | SkeletonExtractor.cpp | filterSkeletonInsideClosedSurface |
| 第3节 节点合并 | CageInteractorStyle.cpp | MergeSkeletonNodes |
| 第4.1节 分支区域与连续段识别 | StraightSkeletonSegmenter.cpp | ExtractStraightSegments、TrimPolyline |
| 第4.2节 弧弦比与直线段判断 | StraightSkeletonSegmenter.cpp | IsPolylineStraight、AppendStraightRuns |
| 第4.2节 连续段控制网格 | CageInteractorStyle.cpp | BuildStraightCagesFromSegments、ComputeRMF |
| 第4.3.2节 共面分支连接 | CageInteractorStyle.cpp | GenerateNBranchFromRings、GenerateNBranchTriHub2D、GenerateNBranch |
| 第4.3.3节 三维分支连接 | CageInteractorStyle.cpp | GenerateNBranch3DSkeleton |
| 细分曲面与全局控制拓扑 | CageInteractorStyle.cpp | RebuildGlobalMesh |
| 曲面实体化和 STEP 导出 | CageInteractorStyle.cpp | outfile |

表中对应关系用于定位实现，不声称论文的每个公式都在当前代码中逐项以同名参数出现。

## 已核实的实现细节

`extractSkeleton` 直接读入三角网格，调用 CGAL `contract_until_convergence` 和 `convert_to_skeleton`，没有在此入口先执行各向同性重网格化。局部骨架半径取关联表面点距离的75%分位数。

原始直线判断同时要求弧弦比不超过1.04、最大偏离不超过 max(弦长×0.035, 模型对角线×0.003)、局部方向偏角不超过12°。分支端部裁剪距离为局部表面距离的2.5倍，并受到段长比例上限约束。缺少表面定位器时采用模型对角线×0.035的回退尺度。

原始分段入口跳过拓扑闭环或近似几何闭环；闭环需要后续专门的拓扑处理，不能用分段示例独立完成。

## 自动计算和人工交互

骨架提取、非设计域过滤、连续段识别与控制网格构建分别具有计算入口。源程序还提供节点合并、删边、骨架延伸、控制网格调整和旋转等交互模式。分支与导出函数依赖已选择端面、场景状态或用户对话框，因此本材料不声称原始软件全流程无需人工操作。

## 简化版本改动

只改写分段模块的表面距离接口并移除显示层。传入 `SurfaceDistance` 回调时应返回输入点到设计域表面的非负最近距离，且与坐标使用相同长度单位。默认不传回调时，示例结果采用回退尺度；与使用真实表面距离的原程序裁剪结果可能不同。其余分段数值常量和算法逻辑保留。
