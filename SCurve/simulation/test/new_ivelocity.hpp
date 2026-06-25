/**
 * @file    IVelocityProfile.hpp
 * @author  syhanjin
 * @date    2026-01-29
 * @brief   速度曲线规划基类
 */
#ifndef NEW_IVELOCITY_PROFILE_HPP
#define NEW_IVELOCITY_PROFILE_HPP

namespace scurve_new
{

class IVelocityProfile
{
public:
    virtual ~IVelocityProfile() = default;

    virtual bool success() const { return false; }

    // calculate position at time t
    virtual float CalcX(float t) const = 0;

    // calculate velocity at time t
    virtual float CalcV(float t) const = 0;

    // calculate acceleration at time t
    virtual float CalcA(float t) const = 0;

    // get total duration of the motion
    virtual float getTotalTime() const = 0;
};

} // namespace scurve_new

#endif // NEW_IVELOCITY_PROFILE_HPP
