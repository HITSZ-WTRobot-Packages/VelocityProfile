/**
 * @file    s_curve.cpp
 * @author  syhanjin LIJunHong659
 * @date    2026-04-09
 */
#include "s_curve.hpp"

#include <cmath>

namespace velocity_profile
{
namespace
{
/**
 * 构造 SCurveProfile 时，需要反复判断：
 * “若峰值速度取 vp，则起点侧和终点侧一共要消耗多少位移？”
 *
 * 为了让这一判断可以在二分搜索中高频执行，这里提供一套只关注位移估算的轻量结构，
 * 避免在每次试探峰值速度时都完整初始化 SCurveAccel 对象。
 */
constexpr float kHalf     = 0.5f;
constexpr float kOneSixth = 1.0f / 6.0f;
constexpr float kZero     = 0.0f;

struct FastEvalConfig
{
    float am;             ///< 最大加速度。
    float jm;             ///< 最大加加速度。
    float am_square;      ///< am^2，避免在循环内重复乘法。
    float jerk_ramp_time; ///< 从 0 提升到最大加速度所需时长 am / jm。
    float inv_double_am;  ///< 1 / (2 * am)，用于距离公式复用。
};

struct FastEvalSide
{
    float v_base;  ///< 标准单边过程的等效起始速度。
    float x_pre;   ///< 本侧显式预处理已经占用的位移。
    float t_shift; ///< 在标准单边过程上需要跳过的前置时间。
};

struct FastEvalProfile
{
    bool  has_uniform;    ///< 是否存在匀加速段。
    float v_base;         ///< 过程起始速度。
    float vp;             ///< 候选峰值速度。
    float t1;             ///< 加加速段结束时刻。
    float t2;             ///< 匀加速段结束时刻；三角曲线时与 t1 相同。
    float x1;             ///< 加加速段结束位移。
    float v1;             ///< 加加速段结束速度。
    float total_time;     ///< 单边过程总时长。
    float total_distance; ///< 单边过程总位移。
};

struct FastEvalResult
{
    float dx1;   ///< 起点侧消耗位移。
    float dx3;   ///< 终点侧消耗位移。
    float delta; ///< dx1 + dx3 - len，用于判定候选峰速偏大还是偏小。
};

/**
 * @brief 为指定起始速度和候选峰值速度构造单边标准加速过程的轻量描述。
 */
[[nodiscard]] FastEvalProfile BuildFastEvalProfile(const FastEvalConfig& cfg,
                                                   const float           v_base,
                                                   const float           vp)
{
    FastEvalProfile profile{};
    const float     delta_v = vp - v_base;

    profile.v_base = v_base;
    profile.vp     = vp;

    if (cfg.jm * delta_v > cfg.am_square)
    {
        profile.has_uniform    = true;
        profile.t1             = cfg.jerk_ramp_time;
        profile.t2             = delta_v / cfg.am;
        profile.v1             = v_base + kHalf * cfg.am * profile.t1;
        profile.x1             = v_base * profile.t1 + kOneSixth * cfg.am * profile.t1 * profile.t1;
        profile.total_time     = profile.t2 + profile.t1;
        profile.total_distance = (vp * vp - v_base * v_base) * cfg.inv_double_am +
                                 kHalf * (v_base + vp) * profile.t1;
        return profile;
    }

    const float peak_acc = sqrtf(cfg.jm * delta_v);

    profile.has_uniform    = false;
    profile.t1             = peak_acc / cfg.jm;
    profile.t2             = profile.t1;
    profile.v1             = v_base + kHalf * peak_acc * peak_acc / cfg.jm;
    profile.x1             = v_base * profile.t1 + kOneSixth * peak_acc * profile.t1 * profile.t1;
    profile.total_time     = 2.0f * profile.t1;
    profile.total_distance = (v_base + vp) * profile.t1;
    return profile;
}

/**
 * @brief 计算轻量单边过程在时刻 @p t 内已经走过的位移。
 */
[[nodiscard]] float EvaluateFastEvalDistance(const FastEvalConfig&  cfg,
                                             const FastEvalProfile& profile,
                                             const float            t)
{
    if (t <= kZero)
        return 0;

    if (t < profile.t1)
        return profile.v_base * t + kOneSixth * cfg.jm * t * t * t;

    if (profile.has_uniform && t < profile.t2)
    {
        const float dt = t - profile.t1;
        return profile.x1 + profile.v1 * dt + kHalf * cfg.am * dt * dt;
    }

    if (t < profile.total_time)
    {
        const float dt = profile.total_time - t;
        return profile.total_distance - profile.vp * dt + kOneSixth * cfg.jm * dt * dt * dt;
    }

    return profile.total_distance;
}

/**
 * @brief 计算某一侧在候选峰值速度 @p vp 下总共会占用的位移。
 *
 * 位移由两部分组成：
 * - 显式预处理段 `x_pre`；
 * - 标准单边过程的总位移减去时间平移前被“吸收”的那一小段。
 */
[[nodiscard]] float EvaluateSideDistance(const FastEvalConfig& cfg,
                                         const FastEvalSide&   side,
                                         const float           vp)
{
    if (vp <= side.v_base)
        return side.x_pre;

    const FastEvalProfile profile        = BuildFastEvalProfile(cfg, side.v_base, vp);
    const float           shift_distance = EvaluateFastEvalDistance(cfg, profile, side.t_shift);

    return side.x_pre + profile.total_distance - shift_distance;
}

/**
 * @brief 计算当前候选峰值速度相对目标总位移的偏差。
 *
 * delta > 0 表示两侧过程消耗的位移过长，峰值速度应该降低；
 * delta < 0 表示仍有剩余路程，可继续提高峰值速度或引入匀速段。
 */
[[nodiscard]] FastEvalResult EvaluateDistanceDelta(const FastEvalConfig& cfg,
                                                   const FastEvalSide&   start,
                                                   const FastEvalSide&   end,
                                                   const float           len,
                                                   const float           vp)
{
    FastEvalResult result{};
    result.dx1   = EvaluateSideDistance(cfg, start, vp);
    result.dx3   = EvaluateSideDistance(cfg, end, vp);
    result.delta = result.dx1 + result.dx3 - len;
    return result;
}
} // namespace

SCurveProfile::SCurveAccel::SCurveAccel() :
    has_uniform_(false), vs_(0), jm_(0), total_time_(0), total_distance_(0), t1_(0), x1_(0), v1_(0),
    t2_(0), ap_(0), vp_(0)
{
}

/**
 * 内部单边加速器只处理“速度单调上升”的标准情形。
 *
 * 初始化时会根据 `vp - vs` 是否足以碰到最大加速度 `am`，分成两种实现：
 * - 梯形加速度曲线：`加加速 -> 匀加速 -> 减加速`
 * - 三角加速度曲线：`加加速 -> 减加速`
 *
 * 后续 `getDistance/getVelocity/getAcceleration` 只需要按这些阶段分段求值即可。
 */
void SCurveProfile::SCurveAccel::init(const float vs,
                                      const float vp,
                                      const float am,
                                      const float jm)
{
    // 先判断是否能触及最大加速度：
    // - 能触及：加速度曲线为梯形（存在匀加速段）；
    // - 不能触及：加速度曲线为三角形（只经历加加速和减加速）。
    has_uniform_ = jm * (vp - vs) > am * am;
    vs_          = vs;
    vp_          = vp;
    jm_          = jm;
    if (has_uniform_) // 梯形加速曲线（存在匀加速段）
    {
        ap_ = am;

        t1_ = am / jm;        // 加加速段时长
        t2_ = (vp - vs) / am; // 加加速度段与减加速段时刻分界（加加速度+匀加速度一共的时间）

        v1_ = vs + 0.5f * am * t1_; // 加加速段与匀加速段速度分界

        x1_ = vs * t1_ + 1 / 6.0f * am * t1_ * t1_;

        total_time_     = t2_ + t1_;
        total_distance_ = (vp * vp - vs * vs) / (2.0f * am) +
                          0.5f * (vs + vp) * t1_ /*(vs + vp) * am / (2.0f * jm)*/;
    }
    else // 三角加速曲线（无匀加速段）
    {
        ap_ = sqrtf(jm * (vp - vs));

        t1_ = ap_ / jm;
        t2_ = t1_; // 没有匀加速段

        v1_ = vs + 0.5f * ap_ * ap_ / jm;

        x1_ = vs * t1_ + 1 / 6.0f * ap_ * t1_ * t1_;

        total_time_     = 2.0f * sqrtf((vp - vs) / jm);
        total_distance_ = (vs + vp) * sqrtf((vp - vs) / jm);
    }
}
float SCurveProfile::SCurveAccel::getDistance(const float t) const
{
    if (t <= 0)
    {
        return 0;
    }
    if (t < t1_)
    {
        return vs_ * t + 1 / 6.0f * jm_ * t * t * t;
    }
    // 只有存在匀加速平台时才会进入这一段；三角曲线会直接落到最后一段。
    if (has_uniform_ && t < t2_)
    {
        const float _t = t - t1_;
        return x1_ + v1_ * _t + 0.5f * ap_ * _t * _t;
    }
    if (t < total_time_)
    {
        const float _t = total_time_ - t;
        return total_distance_ - vp_ * _t + 1 / 6.0f * jm_ * _t * _t * _t;
    }
    return total_distance_;
}
float SCurveProfile::SCurveAccel::getVelocity(const float t) const
{
    if (t <= 0)
    {
        return vs_;
    }
    if (t < t1_)
    {
        return vs_ + 0.5f * jm_ * t * t;
    }
    if (has_uniform_ && t < t2_)
    {
        return v1_ + ap_ * (t - t1_);
    }
    if (t < total_time_)
    {
        const float _t = total_time_ - t;
        return vp_ - 0.5f * jm_ * _t * _t;
    }
    return vp_;
}
float SCurveProfile::SCurveAccel::getAcceleration(const float t) const
{
    if (t <= 0)
        return 0;
    if (t < t1_)
        return jm_ * t;
    if (has_uniform_ && t < t2_)
        return ap_;
    if (t < total_time_)
        return jm_ * (total_time_ - t);
    return 0;
}

/**
 * SCurveProfile 的求解过程分为 4 步：
 * 1. 方向归一化：把位移统一折算为正向问题，减少分支；
 * 2. 单侧边界规整：把起点/终点的速度与加速度约束转成内部过程可处理的形式；
 * 3. 峰值速度判定：先尝试直接到 `vm` 并插入匀速段；
 * 4. 无匀速时二分：若路程不足以容纳匀速段，则二分峰值速度直到刚好满足位移约束。
 *
 * 求解完成后，轨迹采样阶段再根据时刻落在哪个区间，分别调用起点过程、匀速段或终点逆过程。
 */
SCurveProfile::SCurveProfile(
        const Config& cfg, float xs, float vs, float as, float xe, float ve, float ae)
{
    // 约束全部按绝对值解释，方向统一由起终点位置决定。
    auto vm = fabsf(cfg.max_spd);
    auto am = fabsf(cfg.max_acc);
    auto jm = fabsf(cfg.max_jerk);

    // 全部折算到“正向移动”求解，最后再通过 direction_ 恢复原始方向。
    const float dir = xe > xs ? 1.0f : -1.0f;
    direction_      = dir;
    xs_             = xs;
    xe_             = xe;
    const float len = fabsf(xe - xs);
    vs              = vs * dir;
    as              = as * dir;
    ve              = ve * dir;
    ae              = ae * dir;

    vs_ = vs;
    as_ = as;
    ve_ = ve;
    ae_ = ae;
    jm_ = jm;

    // 如果距离很小，且初末速度/加速度都比较相近，直接返回
    if (len < 1e-3f && fabsf(ve - vs) < 1e-3f && fabsf(ae - as) < 1e-3f)
    {
        t1_pre_     = 0;
        t1_         = 0;
        has_const_  = false;
        t2_         = 0;
        t3_pre_     = 0;
        x3_pre_     = 0;
        ts1_        = 0;
        ts3_        = 0;
        xs1_        = 0;
        xs3_        = 0;
        total_time_ = 0;
        success_    = true;
        return;
    }

    if (fabsf(vs) > vm || fabsf(as) > am || fabsf(ve) > vm || fabsf(ae) > am)
    {
        // 初末边界本身已经违反约束时，不存在可行解。
        success_ = false;
        return;
    }

    /**
     * 把单侧边界条件规整成“标准单边加速器可消费”的形式。
     *
     * 实现分两类：
     * - `a0 < 0`：当前加速度与运动方向相反，必须先显式把加速度抬回 0；
     * - `a0 >= 0`：当前加速度与运动方向同向，可通过时间平移吸收进标准过程。
     *
     * 规整结果会给出：
     * - 显式预处理的时间/位移；
     * - 进入内部单边加速器时的等效基线速度；
     * - 该侧允许的峰值速度下界。
     */
    const auto prepareSide = [vm, jm](const float v0, const float a0)
    {
        SidePrepare ret{};
        ret.valid = true;
        // 两条路径：
        // - a0 < 0：当前加速度方向与运动方向相反，需要先“抬加速度”到 0（显式预处理）
        //   这会产生 t_pre、x_pre，并得到一个等效的起始速度 v_base（抬加速后的基线速度）
        // - a0 >= 0：通过时间平移将非零加速度并入标准加速器（不显式改变位移），
        //   这会产生 t_shift 和 v_base，以及更高的峰速下界 vp_min
        if (a0 < 0)
        {
            // 逆向加速度必须先抬到 0。等效起始速度：
            // v_base = v0 - a0^2 / (2 * jm)
            ret.v_base = v0 - 0.5f * a0 * a0 / jm;
            // 若抬加速度后的速度超过最大速度约束，则本侧无效
            if (fabsf(ret.v_base) > vm)
            {
                ret.valid = false;
                return ret;
            }
            // 抬加速后该侧对 vp 的下界就是 v_base
            ret.vp_min = ret.v_base;
            // 预处理时长 t_pre = -a0 / jm（把加速度升到 0 所需时间）
            ret.t_pre = -a0 / jm;
            // 预处理位移 x_pre = v0 * t_pre + (1/3) * a0 * t_pre^2
            ret.x_pre = v0 * ret.t_pre + 1 / 3.0f * a0 * ret.t_pre * ret.t_pre;
            // 时间平移为 0（无时间平移）
            ret.t_shift = 0;
        }
        else
        {
            // 同向加速度通过时间平移并入基础加速过程
            // 需要的峰速下界：vp_min = v0 + a0^2 / (2 * jm)
            ret.vp_min = v0 + 0.5f * a0 * a0 / jm;
            // 若最大速度不足以满足该下界，则本侧无效
            if (vm < ret.vp_min)
            {
                ret.valid = false;
                return ret;
            }
            // 无显式预处理
            ret.t_pre = 0;
            ret.x_pre = 0;
            // 时间平移 t_shift = a0 / jm（在 SCurve 内跳过相应的前段）
            ret.t_shift = a0 / jm;
            // 等效起始速度 v_base = v0 - 0.5 * a0 * t_shift （与抬加速表达等价）
            ret.v_base = v0 - 0.5f * a0 * ret.t_shift;
        }
        // 下界不能为负（将其截为 0）
        if (ret.vp_min < 0)
            ret.vp_min = 0;
        return ret;
    };

    const SidePrepare start = prepareSide(vs, as);
    if (!start.valid)
    {
        success_ = false;
        return;
    }

    // 终点边界通过时间反演变成“逆过程起点”约束：v_rs = v_e, a_rs = -a_e
    const SidePrepare endr = prepareSide(ve, -ae);
    if (!endr.valid)
    {
        success_ = false;
        return;
    }

    t1_pre_ = start.t_pre;
    x1_pre_ = xs + dir * start.x_pre;
    ts1_    = start.t_shift;
    t3_pre_ = endr.t_pre;
    x3_pre_ = endr.x_pre;
    ts3_    = endr.t_shift;

    // 峰值速度必须同时满足起点侧和终点侧的最低要求。
    float vp_min = fmaxf(start.vp_min, endr.vp_min);
    if (vp_min < 0)
        vp_min = 0;
    if (vm < vp_min)
    {
        success_ = false;
        return;
    }

    // 去掉两端显式预处理后，中间仍需由主过程覆盖的净位移。
    const float len0 = len - start.x_pre - endr.x_pre;
    if (len0 < -S_CURVE_MAX_BS_ERROR)
    {
        // 两端固定预处理已经把路走“超了”，因此一定无解。
        success_ = false;
        return;
    }

    const FastEvalConfig fast_eval_cfg{
        am, jm, am * am, am / jm, 0.5f / am,
    };
    const FastEvalSide start_eval{ start.v_base, start.x_pre, start.t_shift };
    const FastEvalSide end_eval{ endr.v_base, endr.x_pre, endr.t_shift };

    // 先尝试“峰值速度能到 vm，并且中间还留有匀速段”的情况。
    // 若此时两侧过程消耗的位移仍小于总位移，则剩余部分就是匀速段。
    const FastEvalResult vm_eval =
            EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, vm);
    const float x_const = -vm_eval.delta;
    if (x_const > 0)
    {
        process1_.init(start.v_base, vm, am, jm);
        process3_.init(endr.v_base, vm, am, jm);
        xs1_ = process1_.getDistance(ts1_);
        xs3_ = process3_.getDistance(ts3_);

        // 存在匀速段
        has_const_ = true;

        // 这里是时刻分界点，所以需要加上开头部分
        t1_ = t1_pre_ + process1_.getTotalTime() - ts1_;

        const float t_const = x_const / vm;
        t2_                 = t1_ + t_const;
        total_time_         = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;

        x1_ = xs_ + direction_ * vm_eval.dx1;

        vp_ = vm;

        success_ = true;
        return;
    }

    // 若在 vm 下已经没有余量，则轨迹退化为“加速后立即减速”。
    // 由于总位移随峰值速度单调变化，这里可以直接对峰值速度做二分搜索。
    float l = vp_min, r = vm;
    float delta_d = len0, dx1 = 0, dx3 = 0;

#ifdef DEBUG
    binary_search_count_ = 0;
#endif
    while (r - l > S_CURVE_MAX_BS_ERROR)
    {
#ifdef DEBUG
        binary_search_count_++;
#endif
        const float          mid = 0.5f * (l + r);
        const FastEvalResult mid_eval =
                EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, mid);
        delta_d = mid_eval.delta;
        if (fabsf(delta_d) <= S_CURVE_MAX_BS_ERROR)
        {
            r = l = mid;
            break;
        }
        if (delta_d > 0)
            r = mid;
        else
            l = mid;
    }

    // 用最终二分结果重建精确的单边过程，并刷新分界时刻/位置。
    // 轻量估算只用于搜索；真正落盘到成员变量时仍以完整过程为准。
    const float vp = 0.5f * (l + r);
    process1_.init(start.v_base, vp, am, jm);
    process3_.init(endr.v_base, vp, am, jm);
    xs1_    = process1_.getDistance(ts1_);
    xs3_    = process3_.getDistance(ts3_);
    dx1     = start.x_pre + process1_.getTotalDistance() - xs1_;
    dx3     = endr.x_pre + process3_.getTotalDistance() - xs3_;
    delta_d = dx1 + dx3 - len;

    if (delta_d > S_CURVE_MAX_BS_ERROR)
    {
        // 即使 vp 降到最低也无法找到解
        success_ = false;
        return;
    }

    has_const_ = false;
    // 这里是时刻分界点，所以需要加上开头部分
    t1_         = t1_pre_ + process1_.getTotalTime() - ts1_;
    t2_         = t1_;
    total_time_ = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;

    x1_ = xs_ + direction_ * dx1;

    vp_ = vp;

    success_ = true;
}

float SCurveProfile::getReverseDistance(const float tau) const
{
    if (tau <= 0)
        return 0;
    if (tau < t3_pre_)
    {
        // 末端逆过程的显式预处理段。
        // 这里使用的是“从终点往回看”的时间变量 tau，因此符号与正向过程不同。
        const float tau2 = tau * tau;
        const float tau3 = tau2 * tau;
        return ve_ * tau - 0.5f * ae_ * tau2 + 1 / 6.0f * jm_ * tau3;
    }

    // 进入标准单边加速器后，减去时间平移前已经被吸收的位移基准 xs3_。
    return x3_pre_ + process3_.getDistance(tau - t3_pre_ + ts3_) - xs3_;
}

float SCurveProfile::getReverseVelocity(const float tau) const
{
    if (tau <= 0)
        return ve_;
    if (tau < t3_pre_)
        return ve_ - ae_ * tau + 0.5f * jm_ * tau * tau;

    return process3_.getVelocity(tau - t3_pre_ + ts3_);
}

float SCurveProfile::getReverseAcceleration(const float tau) const
{
    if (tau <= 0)
        return -ae_;
    if (tau < t3_pre_)
        return -ae_ + jm_ * tau;

    return process3_.getAcceleration(tau - t3_pre_ + ts3_);
}

float SCurveProfile::CalcX(const float t) const
{
    if (!success_)
        return 0;
    // 外部采样统一走“按时间落区间分段求值”的实现：
    // 起点预处理 -> 主加速 -> 匀速 -> 末端逆过程。
    // 起始之前
    if (t <= 0)
        return xs_;
    // 起点预处理：先把与运动方向相反的初始加速度抬回 0。
    if (t < t1_pre_)
    {
        const float t2 = t * t;
        const float t3 = t2 * t;
        return xs_ + direction_ * (vs_ * t + 0.5f * as_ * t2 + 1 / 6.0f * jm_ * t3);
    }
    // 主加速段：调用标准单边加速器，并扣掉时间平移前已经吸收的那部分位移。
    if (t < t1_)
        return x1_pre_ + direction_ * (process1_.getDistance(t - t1_pre_ + ts1_) - xs1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return x1_ + direction_ * vp_ * (t - t1_);
    // 减速段通过“从终点反推”计算，能自然兼容终点速度/加速度约束。
    if (t < total_time_)
        return xe_ - direction_ * getReverseDistance(total_time_ - t);
    return xe_;
}

float SCurveProfile::CalcV(const float t) const
{
    if (!success_)
        return 0;
    // 起始之前
    if (t <= 0)
        return direction_ * vs_;
    // 起点预处理段。
    if (t < t1_pre_)
        return direction_ * (vs_ + as_ * t + 0.5f * jm_ * t * t);
    // 主加速段。
    if (t < t1_)
        return direction_ * process1_.getVelocity(t - t1_pre_ + ts1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return direction_ * vp_;
    // 逆过程映射到正向时间后的减速段。
    if (t < total_time_)
        return direction_ * getReverseVelocity(total_time_ - t);
    return direction_ * ve_;
}

float SCurveProfile::CalcA(const float t) const
{
    if (!success_)
        return 0;
    // 起始之前
    if (t <= 0)
        return direction_ * as_;
    // 起点预处理段。
    if (t < t1_pre_)
        return direction_ * (as_ + jm_ * t);
    // 主加速段。
    if (t < t1_)
        return direction_ * process1_.getAcceleration(t - t1_pre_ + ts1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return 0;
    // 逆过程映射到正向时间后的减速段。
    if (t < total_time_)
        return -direction_ * getReverseAcceleration(total_time_ - t);
    return direction_ * ae_;
}
} // namespace velocity_profile
