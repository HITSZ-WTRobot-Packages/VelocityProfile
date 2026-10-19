# VelocityProfile

面向嵌入式运动控制的单轴速度曲线包集合，作为宿主 STM32 项目的 Git SubModule 使用，不提供独立固件构建入口。

## 包与导航

| 包 | 职责 | 文档 |
| --- | --- | --- |
| `VelocityProfile::Core` | 与具体算法解耦的位置、速度、加速度采样接口 | [Core/README.md](Core/README.md) |
| `VelocityProfile::SCurve` | 依赖 Core，提供带初末状态和速度、加速度、加加速度约束的 S 曲线规划 | [SCurve/README.md](SCurve/README.md) |

包边界和依赖以各目录的 `cpkg.toml` 为准。Core 定义抽象，SCurve 实现规划；控制周期调度、控制器跟踪及执行器驱动由消费方负责。

## 阅读入口

- 接口的职责与单位边界：[Core](Core/README.md)。
- 算法使用、可行性限制及完整数学推导：[SCurve](SCurve/README.md)。
- Python 数学仿真命令：[SCurve 仿真与部署注意](SCurve/README.md#仿真与部署注意)。仿真工具不是仓库级固件构建或硬件兼容性验证。
- 维护规则：[AGENTS.md](AGENTS.md)。
