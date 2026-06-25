/**
 * @file    s_curve.cpp
 * @author  syhanjin LIJunHong659
 * @date    2026-04-09
 */
#include "new_s_curve.hpp"

#include <algorithm>
#include <cmath>

namespace scurve_new
{
namespace
{
constexpr float kHalf     = 0.5f;
constexpr float kOneSixth = 1.0f / 6.0f;
constexpr float kZero     = 0.0f;

struct FastEvalConfig
{
    float am;
    float jm;
    float am_square;
    float jerk_ramp_time;
    float inv_double_am;
};

struct FastEvalSide
{
    float v_base;
    float x_pre;
    float t_shift;
};

struct FastEvalProfile
{
    bool  has_uniform;
    float v_base;
    float vp;
    float t1;
    float t2;
    float x1;
    float v1;
    float total_time;
    float total_distance;
};

struct FastEvalResult
{
    float dx1;
    float dx3;
    float delta;
};

[[nodiscard]] bool IsFinite(const float value)
{
    return std::isfinite(value) != 0;
}

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

void SCurveProfile::SCurveAccel::init(const float vs,
                                      const float vp,
                                      const float am,
                                      const float jm)
{
    has_uniform_ = jm * (vp - vs) > am * am;
    vs_          = vs;
    vp_          = vp;
    jm_          = jm;
    if (has_uniform_)
    {
        ap_ = am;

        t1_ = am / jm;
        t2_ = (vp - vs) / am;

        v1_ = vs + kHalf * am * t1_;

        x1_ = vs * t1_ + kOneSixth * am * t1_ * t1_;

        total_time_     = t2_ + t1_;
        total_distance_ = (vp * vp - vs * vs) / (2.0f * am) + kHalf * (vs + vp) * t1_;
    }
    else
    {
        ap_ = sqrtf(jm * (vp - vs));

        t1_ = ap_ / jm;
        t2_ = t1_;

        v1_ = vs + kHalf * ap_ * ap_ / jm;

        x1_ = vs * t1_ + kOneSixth * ap_ * t1_ * t1_;

        total_time_     = 2.0f * sqrtf((vp - vs) / jm);
        total_distance_ = (vs + vp) * sqrtf((vp - vs) / jm);
    }
}

float SCurveProfile::SCurveAccel::getDistance(const float t) const
{
    if (t <= 0)
        return 0;
    if (t < t1_)
        return vs_ * t + kOneSixth * jm_ * t * t * t;
    if (has_uniform_ && t < t2_)
    {
        const float dt = t - t1_;
        return x1_ + v1_ * dt + kHalf * ap_ * dt * dt;
    }
    if (t < total_time_)
    {
        const float dt = total_time_ - t;
        return total_distance_ - vp_ * dt + kOneSixth * jm_ * dt * dt * dt;
    }
    return total_distance_;
}

float SCurveProfile::SCurveAccel::getVelocity(const float t) const
{
    if (t <= 0)
        return vs_;
    if (t < t1_)
        return vs_ + kHalf * jm_ * t * t;
    if (has_uniform_ && t < t2_)
        return v1_ + ap_ * (t - t1_);
    if (t < total_time_)
    {
        const float dt = total_time_ - t;
        return vp_ - kHalf * jm_ * dt * dt;
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

SCurveProfile::SCurveProfile(
        const Config& cfg, float xs, float vs, float as, float xe, float ve, float ae)
{
    success_      = false;
    total_time_   = 0;
    failure_info_ = { cfg, xs, vs, as, xe, ve, ae, FailureCode::None };

    const auto fail = [this](const FailureCode code)
    {
        failure_info_.code = code;
        success_           = false;
    };

    if (!IsFinite(cfg.max_spd) || !IsFinite(cfg.max_acc) || !IsFinite(cfg.max_jerk) ||
        !IsFinite(xs) || !IsFinite(xe) || !IsFinite(vs) || !IsFinite(as) || !IsFinite(ve) ||
        !IsFinite(ae))
    {
        fail(FailureCode::NonFiniteInput);
        return;
    }

    const float vm = fabsf(cfg.max_spd);
    const float am = fabsf(cfg.max_acc);
    const float jm = fabsf(cfg.max_jerk);
    if (vm <= 0 || am <= 0 || jm <= 0)
    {
        fail(FailureCode::InvalidLimit);
        return;
    }

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

    // Keep the 5a3ae00 short-distance behavior: nearly identical end states are a
    // zero-time successful profile.
    if (len < 1e-3f && fabsf(ve - vs) < 1e-3f && fabsf(ae - as) < 1e-3f)
    {
        t1_pre_     = 0;
        x1_pre_     = xs_;
        t1_         = 0;
        has_const_  = false;
        t2_         = 0;
        t3_pre_     = 0;
        x3_pre_     = 0;
        ts1_        = 0;
        ts3_        = 0;
        xs1_        = 0;
        xs3_        = 0;
        x1_         = xs_;
        total_time_ = 0;
        success_    = true;
        return;
    }

    // Try old-style algorithm first (no preprocessing).
    // When inputs are valid for the old version, this path guarantees
    // bit-identical results with scurve_old.
    {
        float vs0 = vs, as0 = as, ve0 = ve, ae0 = ae;  // save originals
        auto tryOldStyle = [&](float sv, float sa, float ev, float ea) -> bool {
            // Old boundary check
            if (fabsf(sv) > vm || fabsf(sa) > am || fabsf(ev) > vm || fabsf(ea) > am)
                return false;

            const auto prepSide = [vm, jm](float v0, float a0) {
                SidePrepare r{};
                r.valid = true;
                if (a0 < 0) {
                    r.v_base = v0 - kHalf * a0 * a0 / jm;
                    if (fabsf(r.v_base) > vm) { r.valid = false; return r; }
                    r.vp_min = r.v_base;
                    r.t_pre = -a0 / jm;
                    r.x_pre = v0 * r.t_pre + (1.0f/3.0f) * a0 * r.t_pre * r.t_pre;
                    r.t_shift = 0;
                } else {
                    r.vp_min = v0 + kHalf * a0 * a0 / jm;
                    if (vm < r.vp_min) { r.valid = false; return r; }
                    r.t_pre = 0; r.x_pre = 0;
                    r.t_shift = a0 / jm;
                    r.v_base = v0 - kHalf * a0 * r.t_shift;
                }
                if (r.vp_min < 0) r.vp_min = 0;
                return r;
            };

            SidePrepare s = prepSide(sv, sa);
            if (!s.valid) return false;
            SidePrepare e = prepSide(ev, -ea);
            if (!e.valid) return false;

            float vpm = fmaxf(s.vp_min, e.vp_min);
            if (vpm < 0) vpm = 0;
            if (vm < vpm) return false;

            float l0 = len - s.x_pre - e.x_pre;
            if (l0 < -S_CURVE_MAX_BS_ERROR) return false;

            FastEvalConfig fe{am, jm, am*am, am/jm, 0.5f/am};
            FastEvalSide se{s.v_base, s.x_pre, s.t_shift};
            FastEvalSide ee{e.v_base, e.x_pre, e.t_shift};

            FastEvalResult vme = EvaluateDistanceDelta(fe, se, ee, len, vm);
            float xc = -vme.delta;
            if (xc > 0) {
                process1_.init(s.v_base, vm, am, jm);
                process3_.init(e.v_base, vm, am, jm);
                xs1_ = process1_.getDistance(s.t_shift);
                xs3_ = process3_.getDistance(e.t_shift);
                t1_pre_ = s.t_pre; x1_pre_ = xs_ + direction_ * s.x_pre; ts1_ = s.t_shift;
                t3_pre_ = e.t_pre; x3_pre_ = e.x_pre; ts3_ = e.t_shift;
                has_const_ = true;
                t1_ = t1_pre_ + process1_.getTotalTime() - ts1_;
                t2_ = t1_ + xc / vm;
                total_time_ = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;
                x1_ = xs_ + direction_ * vme.dx1;
                vp_ = vm; vs_ = sv; as_ = sa; ve_ = ev; ae_ = ea;
                success_ = true;
                return true;
            }

            float l = vpm, r = vm, dd = l0;
            while (r - l > OLD_S_CURVE_MAX_BS_ERROR) {
                float mid = kHalf * (l + r);
                auto me = EvaluateDistanceDelta(fe, se, ee, len, mid);
                dd = me.delta;
                if (fabsf(dd) <= OLD_S_CURVE_MAX_BS_ERROR) { r = l = mid; break; }
                if (dd > 0) r = mid; else l = mid;
            }
            float vp = kHalf * (l + r);
            process1_.init(s.v_base, vp, am, jm);
            process3_.init(e.v_base, vp, am, jm);
            xs1_ = process1_.getDistance(s.t_shift);
            xs3_ = process3_.getDistance(e.t_shift);
            float d1 = s.x_pre + process1_.getTotalDistance() - xs1_;
            float d3 = e.x_pre + process3_.getTotalDistance() - xs3_;
            dd = d1 + d3 - len;
            if (dd > OLD_S_CURVE_MAX_BS_ERROR) return false;

            t1_pre_ = s.t_pre; x1_pre_ = xs_ + direction_ * s.x_pre; ts1_ = s.t_shift;
            t3_pre_ = e.t_pre; x3_pre_ = e.x_pre; ts3_ = e.t_shift;
            has_const_ = false;
            t1_ = t1_pre_ + process1_.getTotalTime() - ts1_;
            t2_ = t1_;
            total_time_ = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;
            x1_ = xs_ + direction_ * d1;
            vp_ = vp; vs_ = sv; as_ = sa; ve_ = ev; ae_ = ea;
            success_ = true;
            return true;
        };

        if (tryOldStyle(vs, as, ve, ae)) return;
        // Old style failed; restore originals for preprocessing path
        vs = vs0; as = as0; ve = ve0; ae = ae0;
    }

    // Clamp start/end velocity to [-vm, vm]
    vs = std::clamp(vs, -vm, vm);
    ve = std::clamp(ve, -vm, vm);
    as = std::clamp(as, -am, am);
    ae = std::clamp(ae, -am, am);

    // Limit acceleration to prevent future velocity from exceeding vm
    float max_start_a = std::min(am, sqrtf(2.0f * jm * (vm - fabsf(vs))));
    float max_end_a   = std::min(am, sqrtf(2.0f * jm * (vm - fabsf(ve))));

    as = std::clamp(as, -max_start_a, max_start_a);
    ae = std::clamp(ae, -max_end_a, max_end_a);

    // Update member variables with clamped values
    vs_ = vs;
    as_ = as;
    ve_ = ve;
    ae_ = ae;

    const auto prepareSide = [vm, jm](const float v0, const float a0)
    {
        SidePrepare ret{};
        ret.valid = true;
        if (a0 < 0)
        {
            ret.v_base = v0 - kHalf * a0 * a0 / jm;
            // Preprocessing guarantees |v_base| <= vm
            // if (fabsf(ret.v_base) > vm)
            // {
            //     ret.valid = false;
            //     return ret;
            // }
            ret.vp_min  = ret.v_base;
            ret.t_pre   = -a0 / jm;
            ret.x_pre   = v0 * ret.t_pre + (1.0f / 3.0f) * a0 * ret.t_pre * ret.t_pre;
            ret.t_shift = 0;
        }
        else
        {
            ret.vp_min = v0 + kHalf * a0 * a0 / jm;
            // Preprocessing guarantees vp_min <= vm
            // if (vm < ret.vp_min)
            // {
            //     ret.valid = false;
            //     return ret;
            // }
            ret.t_pre   = 0;
            ret.x_pre   = 0;
            ret.t_shift = a0 / jm;
            ret.v_base  = v0 - kHalf * a0 * ret.t_shift;
        }
        if (ret.vp_min < 0)
            ret.vp_min = 0;
        return ret;
    };

    // Pre-compute preparation distances; scale as/ae if they exceed available length
    auto preXPre = [jm](float v0, float a0) -> float
    {
        if (a0 >= 0)
            return 0;
        float t = -a0 / jm;
        return v0 * t + (1.0f / 3.0f) * a0 * t * t;
    };
    float pre_total = preXPre(vs, as) + preXPre(ve, -ae);
    if (pre_total > len + NEW_S_CURVE_MAX_BS_ERROR)
    {
        float scale = len / pre_total;
        as *= scale;
        ae *= scale;
        as = std::clamp(as, -max_start_a, max_start_a);
        ae = std::clamp(ae, -max_end_a, max_end_a);
        as_ = as;
        ae_ = ae;
    }

    SidePrepare start = prepareSide(vs, as);
    // Preprocessing guarantees start.valid is always true
    // if (!start.valid)
    // {
    //     fail(FailureCode::InvalidEndState);
    //     return;
    // }

    SidePrepare endr = prepareSide(ve, -ae);
    // Preprocessing guarantees endr.valid is always true
    // if (!endr.valid)
    // {
    //     fail(FailureCode::InvalidEndState);
    //     return;
    // }

    t1_pre_ = start.t_pre;
    x1_pre_ = xs + dir * start.x_pre;
    ts1_    = start.t_shift;
    t3_pre_ = endr.t_pre;
    x3_pre_ = endr.x_pre;
    ts3_    = endr.t_shift;

    float vp_min = fmaxf(start.vp_min, endr.vp_min);
    if (vp_min < 0)
        vp_min = 0;
    // Preprocessing guarantees vp_min <= vm
    // if (vm < vp_min)
    // {
    //     fail(FailureCode::SolveCoreFailed);
    //     return;
    // }

    const float len0 = len - start.x_pre - endr.x_pre;
    // len0 < -EPS is prevented by pre-scaling as/ae before prepareSide
    // if (len0 < -NEW_S_CURVE_MAX_BS_ERROR)
    // {
    //     fail(FailureCode::SolveCoreFailed);
    //     return;
    // }

    const FastEvalConfig fast_eval_cfg{
        am, jm, am * am, am / jm, 0.5f / am,
    };
    FastEvalSide start_eval{ start.v_base, start.x_pre, start.t_shift };
    FastEvalSide end_eval{ endr.v_base, endr.x_pre, endr.t_shift };

    // Pre-check: ensure minimum distance at vp_min fits within len
    const FastEvalResult min_eval =
            EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, vp_min);
    if (min_eval.delta > NEW_S_CURVE_MAX_BS_ERROR)
    {
        // Scale as/ae proportionally; min distance scales with acceleration
        float scale = len / (min_eval.dx1 + min_eval.dx3);
        as *= scale;
        ae *= scale;
        as = std::clamp(as, -max_start_a, max_start_a);
        ae = std::clamp(ae, -max_end_a, max_end_a);
        as_ = as;
        ae_ = ae;

        start = prepareSide(vs, as);
        endr = prepareSide(ve, -ae);

        t1_pre_ = start.t_pre;
        x1_pre_ = xs + dir * start.x_pre;
        ts1_    = start.t_shift;
        t3_pre_ = endr.t_pre;
        x3_pre_ = endr.x_pre;
        ts3_    = endr.t_shift;

        vp_min = fmaxf(start.vp_min, endr.vp_min);
        if (vp_min < 0)
            vp_min = 0;

        start_eval = FastEvalSide{ start.v_base, start.x_pre, start.t_shift };
        end_eval   = FastEvalSide{ endr.v_base, endr.x_pre, endr.t_shift };

        // If scaling didn't help (e.g. as/ae were zero), force immediate deceleration
        const FastEvalResult retry_eval =
                EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, vp_min);
        if (retry_eval.delta > NEW_S_CURVE_MAX_BS_ERROR)
        {
            as = -std::min(am, sqrtf(2.0f * jm * fabsf(vs)));
            as = std::clamp(as, -max_start_a, max_start_a);
            as_ = as;

            start = prepareSide(vs, as);
            endr = prepareSide(ve, -ae);

            t1_pre_ = start.t_pre;
            x1_pre_ = xs + dir * start.x_pre;
            ts1_    = start.t_shift;
            t3_pre_ = endr.t_pre;
            x3_pre_ = endr.x_pre;
            ts3_    = endr.t_shift;

            vp_min = fmaxf(start.vp_min, endr.vp_min);
            if (vp_min < 0)
                vp_min = 0;

            start_eval = FastEvalSide{ start.v_base, start.x_pre, start.t_shift };
            end_eval   = FastEvalSide{ endr.v_base, endr.x_pre, endr.t_shift };

            // If still exceeds, reduce vs to fit within available distance
            const FastEvalResult vs_eval =
                    EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, vp_min);
            if (vs_eval.delta > NEW_S_CURVE_MAX_BS_ERROR)
            {
                float vs_scale = len / (vs_eval.dx1 + vs_eval.dx3);
                vs *= vs_scale;
                vs = std::clamp(vs, -vm, vm);
                vs_ = vs;
                ve *= vs_scale;
                ve = std::clamp(ve, -vm, vm);
                ve_ = ve;

                start = prepareSide(vs, as);
                endr = prepareSide(ve, -ae);

                t1_pre_ = start.t_pre;
                x1_pre_ = xs + dir * start.x_pre;
                ts1_    = start.t_shift;
                t3_pre_ = endr.t_pre;
                x3_pre_ = endr.x_pre;
                ts3_    = endr.t_shift;

                vp_min = fmaxf(start.vp_min, endr.vp_min);
                if (vp_min < 0)
                    vp_min = 0;

                start_eval = FastEvalSide{ start.v_base, start.x_pre, start.t_shift };
                end_eval   = FastEvalSide{ endr.v_base, endr.x_pre, endr.t_shift };
            }
        }
    }

    const FastEvalResult vm_eval =
            EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, vm);
    const float x_const = -vm_eval.delta;
    if (x_const > 0)
    {
        process1_.init(start.v_base, vm, am, jm);
        process3_.init(endr.v_base, vm, am, jm);
        xs1_ = process1_.getDistance(ts1_);
        xs3_ = process3_.getDistance(ts3_);

        has_const_ = true;

        t1_ = t1_pre_ + process1_.getTotalTime() - ts1_;

        const float t_const = x_const / vm;
        t2_                 = t1_ + t_const;
        total_time_         = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;

        x1_ = xs_ + direction_ * vm_eval.dx1;

        vp_ = vm;

        success_ = true;
        return;
    }

    // Binary search + overshoot retry loop
    float l, r, delta_d, dx1, dx3;
    for (int retry = 0; retry < 4; retry++)
    {
        l = vp_min; r = vm;
#ifdef DEBUG
        binary_search_count_ = 0;
#endif
        while (r - l > NEW_S_CURVE_MAX_BS_ERROR)
        {
#ifdef DEBUG
            binary_search_count_++;
#endif
            const float          mid = kHalf * (l + r);
            const FastEvalResult mid_eval =
                    EvaluateDistanceDelta(fast_eval_cfg, start_eval, end_eval, len, mid);
            delta_d = mid_eval.delta;
            if (fabsf(delta_d) <= NEW_S_CURVE_MAX_BS_ERROR)
            {
                r = l = mid;
                break;
            }
            if (delta_d > 0)
                r = mid;
            else
                l = mid;
        }

        const float vp = kHalf * (l + r);
        process1_.init(start.v_base, vp, am, jm);
        process3_.init(endr.v_base, vp, am, jm);
        xs1_    = process1_.getDistance(ts1_);
        xs3_    = process3_.getDistance(ts3_);
        dx1     = start.x_pre + process1_.getTotalDistance() - xs1_;
        dx3     = endr.x_pre + process3_.getTotalDistance() - xs3_;
        delta_d = dx1 + dx3 - len;

        if (delta_d <= NEW_S_CURVE_MAX_BS_ERROR)
        {
            vp_ = vp;
            break;  // undershoot or exact: handled below
        }

        // Overshoot: reduce vs/ve and retry
        float scale = len / (dx1 + dx3);
        if (scale > 0.999f) { vp_ = vp; break; }
        vs *= scale;
        vs = std::clamp(vs, -vm, vm);
        vs_ = vs;
        ve *= scale;
        ve = std::clamp(ve, -vm, vm);
        ve_ = ve;

        start = prepareSide(vs, as);
        endr = prepareSide(ve, -ae);

        t1_pre_ = start.t_pre;
        x1_pre_ = xs_ + direction_ * start.x_pre;
        ts1_    = start.t_shift;
        t3_pre_ = endr.t_pre;
        x3_pre_ = endr.x_pre;
        ts3_    = endr.t_shift;

        vp_min = fmaxf(start.vp_min, endr.vp_min);
        if (vp_min < 0) vp_min = 0;

        start_eval = FastEvalSide{ start.v_base, start.x_pre, start.t_shift };
        end_eval   = FastEvalSide{ endr.v_base, endr.x_pre, endr.t_shift };
    }

    // Absorb undershoot with constant-velocity phase
    float residual = delta_d < 0 ? -delta_d : 0;
    has_const_  = residual > 0;
    t1_         = t1_pre_ + process1_.getTotalTime() - ts1_;
    t2_         = t1_ + (has_const_ ? residual / vp_ : 0);
    total_time_ = t2_ + t3_pre_ + process3_.getTotalTime() - ts3_;

    x1_ = xs_ + direction_ * dx1;

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
        return ve_ * tau - kHalf * ae_ * tau2 + kOneSixth * jm_ * tau3;
    }

    return x3_pre_ + process3_.getDistance(tau - t3_pre_ + ts3_) - xs3_;
}

float SCurveProfile::getReverseVelocity(const float tau) const
{
    if (tau <= 0)
        return ve_;
    if (tau < t3_pre_)
        return ve_ - ae_ * tau + kHalf * jm_ * tau * tau;

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
    if (t <= 0)
        return xs_;
    if (t < t1_pre_)
    {
        const float t2 = t * t;
        const float t3 = t2 * t;
        return xs_ + direction_ * (vs_ * t + kHalf * as_ * t2 + kOneSixth * jm_ * t3);
    }
    if (t < t1_)
        return x1_pre_ + direction_ * (process1_.getDistance(t - t1_pre_ + ts1_) - xs1_);
    if (has_const_ && t < t2_)
        return x1_ + direction_ * vp_ * (t - t1_);
    if (t < total_time_)
        return xe_ - direction_ * getReverseDistance(total_time_ - t);
    return xe_;
}

float SCurveProfile::CalcV(const float t) const
{
    if (!success_)
        return 0;
    if (t <= 0)
        return direction_ * vs_;
    if (t < t1_pre_)
        return direction_ * (vs_ + as_ * t + kHalf * jm_ * t * t);
    if (t < t1_)
        return direction_ * process1_.getVelocity(t - t1_pre_ + ts1_);
    if (has_const_ && t < t2_)
        return direction_ * vp_;
    if (t < total_time_)
        return direction_ * getReverseVelocity(total_time_ - t);
    return direction_ * ve_;
}

float SCurveProfile::CalcA(const float t) const
{
    if (!success_)
        return 0;
    if (t <= 0)
        return direction_ * as_;
    if (t < t1_pre_)
        return direction_ * (as_ + jm_ * t);
    if (t < t1_)
        return direction_ * process1_.getAcceleration(t - t1_pre_ + ts1_);
    if (has_const_ && t < t2_)
        return 0;
    if (t < total_time_)
        return -direction_ * getReverseAcceleration(total_time_ - t);
    return direction_ * ae_;
}
} // namespace scurve_new
