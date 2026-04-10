/**
 * @file    s_curve.hpp
 * @author  syhanjin LIJunHong659
 * @date    2026-01-28
 */
#ifndef S_CURVE_HPP
#define S_CURVE_HPP
#include <cstdint>
#include "IVelocityProfile.hpp"

// 最大二分查找误差
#ifndef S_CURVE_MAX_BS_ERROR
#    define S_CURVE_MAX_BS_ERROR (0.001f)
#endif

namespace velocity_profile
{

/**
 * @brief 带初末速度/加速度约束的单轴 S 曲线速度规划器。
 *
 * 对外职责包括：
 * - 根据起点/终点的位置、速度、加速度和全局约束构造一条可采样轨迹；
 * - 通过 CalcX/CalcV/CalcA 提供任意时刻的期望位置、速度、加速度；
 * - 通过 getTotalTime() 与 success() 暴露轨迹总时长和求解结果。
 *
 * 若边界条件与给定限幅不兼容，则构造后 success() 返回 false。
 */
class SCurveProfile final : public IVelocityProfile
{
public:
    /**
     * @brief 轨迹全局约束。
     *
     * 三个参数均按绝对值理解；构造函数内部会先取其绝对值，再根据起终点坐标推导运动方向。
     */
    struct Config
    {
        float max_spd;  ///< 速度上限。
        float max_acc;  ///< 加速度上限。
        float max_jerk; ///< 加加速度上限。
    };

    /**
     * @brief 构造一条满足初末状态约束的 S 曲线。
     *
     * @param cfg  全局约束配置。
     * @param xs   起点位置。
     * @param vs   起点速度。
     * @param as   起点加速度。
     * @param xe   终点位置。
     * @param ve   终点速度，默认为 0。
     * @param ae   终点加速度，默认为 0。
     */
    SCurveProfile(
            const Config& cfg, float xs, float vs, float as, float xe, float ve = 0, float ae = 0);

    /**
     * @brief 计算时刻 @p t 的位置。
     *
     * 当 @p t 超出轨迹末端时返回终点位置；当轨迹构造失败时返回 0。
     */
    [[nodiscard]] float CalcX(float t) const override;

    /**
     * @brief 计算时刻 @p t 的速度。
     *
     * 返回值已恢复到原始坐标方向；当轨迹构造失败时返回 0。
     */
    [[nodiscard]] float CalcV(float t) const override;

    /**
     * @brief 计算时刻 @p t 的加速度。
     *
     * 返回值已恢复到原始坐标方向；当轨迹构造失败时返回 0。
     */
    [[nodiscard]] float CalcA(float t) const override;

    /// @brief 获取整条轨迹的总时长。
    [[nodiscard]] float getTotalTime() const override { return total_time_; }

    /// @brief 判断当前边界条件下是否成功构造出可行轨迹。
    [[nodiscard]] bool success() const override { return success_; }

private:
    /**
     * @brief 内部辅助类型，表示一段单侧单调加速过程。
     */
    class SCurveAccel
    {
    public:
        SCurveAccel();

        /// @brief 初始化一段内部单侧加速过程。
        void init(float vs, float vp, float am, float jm);

        /// @brief 查询单边过程在时刻 @p t 已累计的位移。
        [[nodiscard]] float getDistance(float t) const;

        /// @brief 查询单边过程在时刻 @p t 的速度。
        [[nodiscard]] float getVelocity(float t) const;

        /// @brief 查询单边过程在时刻 @p t 的加速度。
        [[nodiscard]] float getAcceleration(float t) const;

        /// @brief 获取单边过程的总位移。
        [[nodiscard]] float getTotalDistance() const { return total_distance_; }

        /// @brief 获取单边过程的总时长。
        [[nodiscard]] float getTotalTime() const { return total_time_; }

    private:
        bool  has_uniform_; ///< 是否存在匀加速平台。
        float vs_;          ///< 该标准过程的起始速度。
        float jm_;          ///< 该标准过程使用的固定加加速度。

        float total_time_;     ///< 单边过程总时长。
        float total_distance_; ///< 单边过程总位移。

        float t1_; ///< 加加速段结束时刻。
        float x1_; ///< 加加速段结束时累计位移。
        float v1_; ///< 加加速段结束时速度。
        float t2_; ///< 匀加速段结束时刻；无匀加速段时与 t1_ 相等。

        float ap_; ///< 峰值加速度。
        float vp_; ///< 过程目标峰值速度。
    };

    /**
     * @brief 内部边界规整结果。
     */
    struct SidePrepare
    {
        float t_pre;   ///< 预处理时长。
        float x_pre;   ///< 预处理位移。
        float v_base;  ///< 规整后的基线速度。
        float t_shift; ///< 对内部过程的时间偏移。
        float vp_min;  ///< 峰值速度下界。
        bool  valid;   ///< 规整结果是否有效。
    };

    bool success_; ///< 轨迹是否成功构造。

    bool  has_const_; ///< 是否存在匀速段。
    float direction_; ///< 轨迹在原始坐标系下的方向，取值为 +1 或 -1。
    float jm_;        ///< 当前轨迹使用的加加速度上限。
    float vp_;        ///< 当前轨迹的实际峰值速度。

    float xs_; ///< 原始坐标系下的起点位置。
    float xe_; ///< 原始坐标系下的终点位置。
    float ve_; ///< 归一化后的终点速度。
    float ae_; ///< 归一化后的终点加速度。

    float vs_;     ///< 归一化后的起点速度。
    float as_;     ///< 归一化后的起点加速度。
    float t1_pre_; ///< 起点显式预处理时长。
    float x1_pre_; ///< 起点预处理结束后的绝对位置。
    float ts1_;    ///< 起点过程的时间偏移。
    float xs1_;    ///< 起点过程的位移基准。

    float t3_pre_; ///< 末端过程预处理时长。
    float x3_pre_; ///< 末端过程预处理位移。
    float ts3_;    ///< 末端过程的时间偏移。
    float xs3_;    ///< 末端过程的位移基准。

    float t1_; ///< 第一阶段结束时刻。
    float t2_; ///< 第二阶段结束时刻。
    float x1_; ///< 阶段切换位置。

    float total_time_; ///< 整条轨迹总时长。

    SCurveAccel process1_{}; ///< 起点侧内部过程。
    SCurveAccel process3_{}; ///< 终点侧内部过程。

    /// @brief 查询末端过程在反向时间 @p tau 下的累计位移。
    [[nodiscard]] float getReverseDistance(float tau) const;

    /// @brief 查询末端过程在反向时间 @p tau 下的速度。
    [[nodiscard]] float getReverseVelocity(float tau) const;

    /// @brief 查询末端过程在反向时间 @p tau 下的加速度。
    [[nodiscard]] float getReverseAcceleration(float tau) const;

#ifdef DEBUG
    uint32_t binary_search_count_;
#endif
};

} // namespace velocity_profile
#endif // S_CURVE_HPP
