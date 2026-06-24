/**
 * @file    s_curve.hpp
 * @author  syhanjin LIJunHong659
 * @date    2026-01-28
 */
#ifndef S_CURVE_HPP
#define S_CURVE_HPP

#include <cstdint>

#include "IVelocityProfile.hpp"

#ifndef S_CURVE_MAX_BS_ERROR
#    define S_CURVE_MAX_BS_ERROR (0.001f)
#endif

namespace velocity_profile
{

class SCurveProfile final : public IVelocityProfile
{
public:
    struct Config
    {
        float max_spd;
        float max_acc;
        float max_jerk;

        constexpr Config operator*(const float ratio) const
        {
            return { max_spd * ratio, max_acc * ratio, max_jerk * ratio };
        }
        constexpr Config operator/(const float factor) const
        {
            return { max_spd / factor, max_acc / factor, max_jerk / factor };
        }
    };

    enum class FailureCode : uint8_t
    {
        None                        = 0,
        NonFiniteInput              = 1,
        InvalidLimit                = 2,
        InvalidEndState             = 3,
        SolveCoreFailed             = 4,
        VelocityClampRecoveryFailed = 5,
        StopPrefixFailed            = 6,
        PrefixHandoffFailed         = 7,
        InternalError               = 8,
    };

    struct FailureInfo
    {
        Config      cfg{};
        float       xs{};
        float       vs{};
        float       as{};
        float       xe{};
        float       ve{};
        float       ae{};
        FailureCode code{ FailureCode::None };
    };

    SCurveProfile(
            const Config& cfg, float xs, float vs, float as, float xe, float ve = 0, float ae = 0);

    [[nodiscard]] float CalcX(float t) const override;
    [[nodiscard]] float CalcV(float t) const override;
    [[nodiscard]] float CalcA(float t) const override;
    [[nodiscard]] float getTotalTime() const override { return total_time_; }
    [[nodiscard]] bool success() const override { return success_; }
    [[nodiscard]] const FailureInfo& failureInfo() const { return failure_info_; }

private:
    class SCurveAccel
    {
    public:
        SCurveAccel();
        void init(float vs, float vp, float am, float jm);

        [[nodiscard]] float getDistance(float t) const;
        [[nodiscard]] float getVelocity(float t) const;
        [[nodiscard]] float getAcceleration(float t) const;
        [[nodiscard]] float getTotalDistance() const { return total_distance_; }
        [[nodiscard]] float getTotalTime() const { return total_time_; }

    private:
        bool  has_uniform_;
        float vs_;
        float jm_;

        float total_time_;
        float total_distance_;

        float t1_;
        float x1_;
        float v1_;
        float t2_;

        float ap_;
        float vp_;
    };

    struct SidePrepare
    {
        float t_pre;
        float x_pre;
        float v_base;
        float t_shift;
        float vp_min;
        bool  valid;
    };

    bool success_{ false };

    bool  has_const_{ false };
    float direction_{ 1.0f };
    float jm_{ 0.0f };
    float vp_{ 0.0f };

    float xs_{ 0.0f };
    float xe_{ 0.0f };
    float ve_{ 0.0f };
    float ae_{ 0.0f };

    float vs_{ 0.0f };
    float as_{ 0.0f };
    float t1_pre_{ 0.0f };
    float x1_pre_{ 0.0f };
    float ts1_{ 0.0f };
    float xs1_{ 0.0f };

    float t3_pre_{ 0.0f };
    float x3_pre_{ 0.0f };
    float ts3_{ 0.0f };
    float xs3_{ 0.0f };

    float t1_{ 0.0f };
    float t2_{ 0.0f };
    float x1_{ 0.0f };

    float total_time_{ 0.0f };

    SCurveAccel process1_{};
    SCurveAccel process3_{};

    FailureInfo failure_info_{};

    [[nodiscard]] float getReverseDistance(float tau) const;
    [[nodiscard]] float getReverseVelocity(float tau) const;
    [[nodiscard]] float getReverseAcceleration(float tau) const;

#ifdef DEBUG
    uint32_t binary_search_count_{ 0 };
#endif
};

} // namespace velocity_profile

#endif // S_CURVE_HPP
