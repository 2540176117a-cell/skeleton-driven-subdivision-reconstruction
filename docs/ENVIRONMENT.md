# 运行环境与集成说明

简化分段程序只依赖 C++17 标准库。原始核心文件使用 CGAL、Eigen、Boost、Qt、VTK、OpenSubdiv 和 OpenCASCADE。

作者工程配置中可核实的版本包括 Qt 5.9.1、VTK 8.2、Eigen 3.4.0、MSVC v141；`CageInteractorStyle.h` 显式使用 OpenSubdiv `v3_7_0` 命名空间。CGAL、Boost 和 OpenCASCADE 的精确版本未从本次工程配置确认，不应将未核实的版本写成论文实验环境。

原始源代码中的第三方接口和交互层保留，不附带第三方头文件、库、DLL 或原始工程的本机绝对路径构建配置。`outfile` 中的 `C:/` 仅是文件保存对话框的初始位置。

集成时：先将封闭三角网格交给 `SkeletonExtractor`；将非设计域封闭表面交给过滤接口；用参考表面构建 VTK `vtkCellLocator`；调用 `ExtractStraightSegments`；为 `CageInteractorStyle` 设置参考表面与骨架数据；调用建笼入口；按场景状态和端面选择调用分支接口；最后调用导出接口。涉及交互的接口仍需 Qt 应用和 VTK 渲染器，仓库未提供完整主程序。

原始论文、返稿图片、有限元数据库、个人目录、第三方库、完整 GUI 与二进制程序均不在此公开包中。
