# Skeleton-driven subdivision reconstruction for topology optimization

论文《基于骨架驱动的细分网格模型的拓扑优化模型重建方法》的关键步骤源代码与简化运行示例。此材料用于回应审稿意见第5条。

## 内容和范围

- `original/`：从作者现有 C++ 程序直接摘出的6个文件，逐字节保留，包含平均曲率流骨架提取、非设计域骨架过滤、骨架分段、控制网格生成、分支连接、OpenSubdiv 细分及 OpenCASCADE STEP 导出。
- `portable/`：从 `StraightSkeletonSegmenter.cpp` 提取的可独立编译分段示例，保留原始拓扑遍历、分支端部裁剪和直线判断算法；删除 VTK 显示函数，将表面距离查询改为可选回调。默认示例采用原程序已有的尺度回退值。
- `examples/`：人工构造的直线、Y 型和闭环骨架，用于演示算法；不是论文算例，不作为几何误差或力学性能证据。
- `docs/`：论文与函数对应关系、运行环境、人工交互说明、修订回复建议。

## 运行简化程序

需要 CMake 3.16或更高版本，以及支持 C++17 的编译器；此示例不依赖 Qt、VTK、CGAL、OpenSubdiv 或 OpenCASCADE。

```sh
cmake -S . -B build
cmake --build build --config Release
```

Windows 多配置生成器：

```sh
build/Release/segment_demo.exe examples/y_skeleton.txt segments.obj
```

Linux/macOS 或单配置生成器：

```sh
./build/segment_demo examples/y_skeleton.txt segments.obj
```

输入首行是节点数和边数，随后依次为节点的 x y z 坐标及边的两个0起始节点索引。输出 OBJ 中的 `v` 和 `l` 表示分段后的骨架折线，`source_chain` 表示来源骨架链。它不是实体 CAD 或控制网格输出。

## 原始核心代码的使用

原始文件用于公开关键算法实现，依赖作者主程序提供的 VTK 场景和 Qt 交互环境，不能仅凭本仓库的 CMake 构建完整 GUI。具体入口和依赖见 `docs/METHOD_MAP.md` 与 `docs/ENVIRONMENT.md`。

简化程序演示论文第4.1节及第4.2节的骨架段识别部分，不涵盖 STL 输入、表面截面采样、三维分支连接、细分曲面实体化或有限元分析。完整关键实现公开于 `original/CageInteractorStyle.cpp`，不得将此分段示例描述为论文全流程自动复现。

## 溯源和验证

`source_manifest.json` 记录原始文件的 SHA-256；`docs/VALIDATION.md` 记录实际执行的验证结果。公开代码对应2026年10月9日的程序快照，不能据此保证它与论文此前各算例运行时的版本完全相同。

第三方库不随代码分发，使用者应分别遵守其许可。本材料尚未指定通用开源许可证；公开可读不等于授予任意再分发许可。
