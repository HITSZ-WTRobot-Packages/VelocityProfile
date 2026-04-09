# VelocityProfile

面向嵌入式运动控制的速度规划模块集合。

## 模块概览

- `Core/`：公共接口定义，目前包含 `IVelocityProfile`；
- `SCurve/`：带初末状态约束的 S 曲线速度规划实现与仿真工具。

## 构建

- `cmake -S . -B build`
- `cmake --build build`
- `cmake --build build --target VelocityProfileSCurve`

## 文档与验证

- `SCurve/README.md`：模块介绍、轨迹组成和使用方式；
- `SCurve/src/带初末状态约束的 S 形曲线路径规划.md`：算法推导说明；
- `SCurve/simulation/plot_display.py`：单轴曲线验证；
- `SCurve/simulation/plot_chassis_display.py`：底盘运动验证。
