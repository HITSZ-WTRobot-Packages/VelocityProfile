/**
 * @file    s_curve.cpp
 * @author  syhanjin
 * @date    2026-01-28
 */
#include "s_curve.hpp"

#include <cmath>

namespace velocity_profile
{
SCurveProfile::SCurveAccel::SCurveAccel() :
    has_uniform_(false), vs_(0), jm_(0), total_time_(0), total_distance_(0), t1_(0), x1_(0), v1_(0),
    t2_(0), x2_(0), v2_(0), ap_(0), vp_(0)
{
}

void SCurveProfile::SCurveAccel::init(const float vs,
                                      const float vp,
                                      const float am,
                                      const float jm)
{
    has_uniform_ = jm * (vp - vs) > am * am;
    vs_          = vs;
    vp_          = vp;
    jm_          = jm;
    if (has_uniform_)                       // 梯形加速曲线（存在匀加速段）
    {
        ap_ = am;

        t1_             = am / jm;          // 加加速段时长
        t2_             = (vp - vs) / am;   // 加加速度段与减加速段时刻分界（加加速度+匀加速度一共的时间）
        const float dt2 = t2_ - t1_;        // 匀加速段时长

        v1_ = vs + 0.5f * am * t1_;         // 加加速段与匀加速段速度分界
        v2_ = vp - 0.5f * am * t1_;         // 匀加速段与减加速段速度分界

        x1_ = vs * t1_ + 1 / 6.0f * am * t1_ * t1_;
        x2_ = x1_ + v1_ * dt2 + 0.5f * am * dt2 * dt2;

        total_time_     = t2_ + t1_;
        total_distance_ = (vp * vp - vs * vs) / (2.0f * am) +
                          0.5f * (vs + vp) * t1_ /*(vs + vp) * am / (2.0f * jm)*/;
    }
    else                                    // 三角加速曲线（无匀加速段）   
    {
        ap_ =  sqrtf(jm * (vp - vs));

        t1_ = ap_ / jm;
        t2_ = t1_;                          // 没有匀加速段

        v1_ = vs + 0.5f * ap_ * ap_ / jm;
        v2_ = v1_;

        x1_ = vs * t1_ + 1 / 6.0f * ap_ * t1_ * t1_;
        x2_ = x1_;

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
    // 除了这段以外两种情况都一样
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

SCurveProfile::SCurveProfile(const Config& cfg,
                             float         xs,
                             float         vs,
                             float         as,
                             float         xe,
                             float         ve,
                             float         ae)
{
    auto vm =  fabsf(cfg.max_spd);
    auto am =  fabsf(cfg.max_acc);
    auto jm =  fabsf(cfg.max_jerk);

    // 全部归到正向移动判断
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

    // 如果目标基本没有变化，则不动
    if (len < 1e-6f)
    {
        if ( fabsf(vs - ve) > 1e-3f ||  fabsf(as - ae) > 1e-3f)
        {
            success_ = false;
            return;
        }
        t1_pre_         = 0;
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

    if ( fabsf(vs) > vm ||  fabsf(as) > am ||  fabsf(ve) > vm ||  fabsf(ae) > am)
    {
        // 初末边界超限则无法构建曲线
        success_ = false;
        return;
    }

    

    const auto prepareSide = [vm, jm](const float v0, const float a0) {
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
            if ( fabsf(ret.v_base) > vm)
            {
                ret.valid = false;
                return ret;
            }
            // 抬加速后该侧对 vp 的下界就是 v_base
            ret.vp_min  = ret.v_base;
            // 预处理时长 t_pre = -a0 / jm（把加速度升到 0 所需时间）
            ret.t_pre   = -a0 / jm;
            // 预处理位移 x_pre = v0 * t_pre + (1/3) * a0 * t_pre^2
            ret.x_pre   = v0 * ret.t_pre + 1 / 3.0f * a0 * ret.t_pre * ret.t_pre;
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
            ret.t_pre   = 0;
            ret.x_pre   = 0;
            // 时间平移 t_shift = a0 / jm（在 SCurve 内跳过相应的前段）
            ret.t_shift = a0 / jm;
            // 等效起始速度 v_base = v0 - 0.5 * a0 * t_shift （与抬加速表达等价）
            ret.v_base  = v0 - 0.5f * a0 * ret.t_shift;
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
    ars_                   = -ae;
    vrs_                   = ve;
    const SidePrepare endr = prepareSide(vrs_, ars_);
    if (!endr.valid)
    {
        success_ = false;
        return;
    }

    t1_pre_     = start.t_pre;
    x1_pre_     = xs + dir * start.x_pre;
    ts1_    = start.t_shift;
    t3_pre_ = endr.t_pre;
    x3_pre_ = endr.x_pre;
    ts3_    = endr.t_shift;

    float vp_min =  fmaxf(start.vp_min, endr.vp_min);
    if (vp_min < 0)
        vp_min = 0;
    if (vm < vp_min)
    {
        success_ = false;
        return;
    }

    const float len0 = len - start.x_pre - endr.x_pre;
    if (len0 < -S_CURVE_MAX_BS_ERROR)
    {
        // 固定位移已经超过目标距离，无法构造满足边界的轨迹
        success_ = false;
        return;
    }

    // 首先假定存在匀速段
    {
        const float vp = vm;
        process1_.init(start.v_base, vp, am, jm);
        process3_.init(endr.v_base, vp, am, jm);

        xs1_            = process1_.getDistance(ts1_);
        xs3_            = process3_.getDistance(ts3_);
        const float dx1 = start.x_pre + process1_.getTotalDistance() - xs1_;
        const float dx3 = endr.x_pre + process3_.getTotalDistance() - xs3_;
        const float x_const = len - dx1 - dx3;
        if (x_const > 0)
        {
            // 存在匀速段
            has_const_ = true;

            // 这里是时刻分界点，所以需要加上开头部分
            t1_ = t1_pre_ + process1_.getTotalTime() - ts1_;

            const float t_const = x_const / vm;
            t2_                 = t1_ + t_const;
            total_time_         = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;

            x1_ = xs_ + direction_ * dx1;
            x2_ = x1_ + direction_ * x_const;

            vp_ = vm;

            success_ = true;
            return;
        }
    }
    // 不存在匀速段，求最大速度
    // 由于代数表达式太过复杂，这里直接上二分查找
    float l = vp_min, r = vm;
    float delta_d = len0, dx1 = 0, dx3 = 0;

#ifdef DEBUG
    binary_search_count_ = 0;
#endif
    // 最大迭代次数约 13 次
    while (r - l > S_CURVE_MAX_BS_ERROR)
    {
#ifdef DEBUG
        binary_search_count_++;
#endif
        const float mid = 0.5f * (l + r);
        process1_.init(start.v_base, mid, am, jm);
        process3_.init(endr.v_base, mid, am, jm);
        // 更新，用于判断是否找到解
        xs1_    = process1_.getDistance(ts1_);
        xs3_    = process3_.getDistance(ts3_);
        dx1     = start.x_pre + process1_.getTotalDistance() - xs1_;
        dx3     = endr.x_pre + process3_.getTotalDistance() - xs3_;
        delta_d = dx1 + dx3 - len;
        if (delta_d < S_CURVE_MAX_BS_ERROR && delta_d > -S_CURVE_MAX_BS_ERROR)
        {
            r = l = mid;
            break;
        }
        if (delta_d > 0)
            r = mid;
        else
            l = mid;
    }

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
    x2_ = x1_;

    vp_ = vp;

    success_ = true;
}

float SCurveProfile::getReverseDistance(const float tau) const
{
    if (tau <= 0)
        return 0;
    if (tau < t3_pre_)
    {
        const float tau2 = tau * tau;
        const float tau3 = tau2 * tau;
        return vrs_ * tau + 0.5f * ars_ * tau2 + 1 / 6.0f * jm_ * tau3;
    }

    return x3_pre_ + process3_.getDistance(tau - t3_pre_ + ts3_) - xs3_;
}

float SCurveProfile::getReverseVelocity(const float tau) const
{
    if (tau <= 0)
        return vrs_;
    if (tau < t3_pre_)
        return vrs_ + ars_ * tau + 0.5f * jm_ * tau * tau;

    return process3_.getVelocity(tau - t3_pre_ + ts3_);
}

float SCurveProfile::getReverseAcceleration(const float tau) const
{
    if (tau <= 0)
        return ars_;
    if (tau < t3_pre_)
        return ars_ + jm_ * tau;

    return process3_.getAcceleration(tau - t3_pre_ + ts3_);
}

float SCurveProfile::CalcX(const float t) const
{
    if (!success_)
        return 0;
    // 起始之前
    if (t <= 0)
        return xs_;
    // 反向加速度刹车
    if (t < t1_pre_)
    {
        const float t2 = t * t;
        const float t3 = t2 * t;
        return xs_ + direction_ * (vs_ * t + 0.5f * as_ * t2 + 1 / 6.0f * jm_ * t3);
    }
    // 加速过程
    if (t < t1_)
        return x1_pre_ + direction_ * (process1_.getDistance(t - t1_pre_ + ts1_) - xs1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return x1_ + direction_ * vp_ * (t - t1_);
    // 减速过程
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
    // 反向加速度刹车
    if (t < t1_pre_)
        return direction_ * (vs_ + as_ * t + 0.5f * jm_ * t * t);
    // 加速过程
    if (t < t1_)
        return direction_ * process1_.getVelocity(t - t1_pre_ + ts1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return direction_ * vp_;
    // 减速过程
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
    // 反向加速度刹车
    if (t < t1_pre_)
        return direction_ * (as_ + jm_ * t);
    // 加速过程
    if (t < t1_)
        return direction_ * process1_.getAcceleration(t - t1_pre_ + ts1_);
    // 匀速过程
    if (has_const_ && t < t2_)
        return 0;
    // 减速过程
    if (t < total_time_)
        return -direction_ * getReverseAcceleration(total_time_ - t);
    return direction_ * ae_;
}
} // namespace velocity_profile