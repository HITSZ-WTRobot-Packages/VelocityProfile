/**
 * C++ SCurve evaluation module — mirrors s_curve_eval.py
 * Samples at 0.001s, checks continuity, limits, boundary conditions.
 */
#pragma once
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

struct EvalResult
{
    bool   planned;
    bool   boundary_ok;
    bool   continuity_ok;
    bool   limits_ok;
    bool   all_pass;
    float  max_v;
    float  max_a;
    float  max_j;
    std::vector<std::string> failures;
};

template <typename SCurve>
EvalResult evaluate(const SCurve& curve, float xs, float xe,
                    float vm, float am, float jm)
{
    EvalResult r{};
    r.planned     = curve.success();
    r.boundary_ok = true;
    r.continuity_ok = true;
    r.limits_ok   = true;

    if (!curve.success())
    {
        r.failures.push_back("planning failed");
        r.all_pass = false;
        if (curve.getTotalTime() <= 0.0f) return r;
    }

    float T = curve.getTotalTime();
    if (T < 1e-6f)
    {
        r.boundary_ok   = true;
        r.continuity_ok = true;
        r.limits_ok     = true;
        r.all_pass      = r.planned;
        return r;
    }

    // Sample at 0.001s
    constexpr float dt    = 0.001f;
    int             n     = (int)(T / dt) + 1;
    if (n < 2) { n = 2; }

    std::vector<float> xv(n + 1), vv(n + 1), av(n + 1);
    for (int i = 0; i <= n; i++)
    {
        float t = (i < n) ? (float)i * dt : T;
        xv[i] = curve.CalcX(t);
        vv[i] = curve.CalcV(t);
        av[i] = curve.CalcA(t);
    }

    // --- Boundary check (relative tolerance) ---
    float scale = fmaxf(fabsf(xs), fabsf(xe));
    float bound_eps = fmaxf(1e-5f, 5e-6f * scale);
    if (fabsf(xv[0] - xs) >= bound_eps)
    {
        r.failures.push_back("x(0) != xs");
        r.boundary_ok = false;
    }
    if (fabsf(xv[n] - xe) >= bound_eps)
    {
        r.failures.push_back("x(T) != xe");
        r.boundary_ok = false;
    }

    // --- Continuity ---
    float max_dx = vm * dt;
    float max_dv = am * dt;
    float max_da = jm * dt;
    float thresh_x = fmaxf(1e-3f, 5.0f * max_dx);
    float thresh_v = fmaxf(0.1f, fmaxf(5.0f * max_dv, 0.05f * vm));
    float thresh_a = fmaxf(0.5f, 5.0f * max_da);

    for (int i = 0; i < n; i++)
    {
        float dx = fabsf(xv[i + 1] - xv[i]);
        float dv = fabsf(vv[i + 1] - vv[i]);
        float da = fabsf(av[i + 1] - av[i]);

        if (dx > thresh_x)
        {
            r.failures.push_back("x discontinuity");
            r.continuity_ok = false;
            break;
        }
        if (dv > thresh_v)
        {
            r.failures.push_back("v discontinuity");
            r.continuity_ok = false;
        }
        if (da > thresh_a)
        {
            r.failures.push_back("a discontinuity");
            r.continuity_ok = false;
        }
    }

    // --- Limits ---
    r.max_v = 0; r.max_a = 0; r.max_j = 0;
    for (int i = 0; i <= n; i++)
    {
        if (fabsf(vv[i]) > r.max_v) r.max_v = fabsf(vv[i]);
        if (fabsf(av[i]) > r.max_a) r.max_a = fabsf(av[i]);
    }
    for (int i = 0; i < n; i++)
    {
        float j = (av[i + 1] - av[i]) / dt;
        if (fabsf(j) > r.max_j) r.max_j = fabsf(j);
    }

    float vm_tol = fmaxf(0.01f, 1e-4f * vm);
    float am_tol = fmaxf(0.1f, 1e-3f * am);
    if (r.max_v > vm + vm_tol)
    {
        r.failures.push_back("v exceeds limit");
        r.limits_ok = false;
    }
    if (r.max_a > am + am_tol)
    {
        r.failures.push_back("a exceeds limit");
        r.limits_ok = false;
    }
    // Jerk from numerical diff can spike at phase boundaries
    float jm_tol = fmaxf(1000.0f, jm * 0.5f);
    if (r.max_j > jm + jm_tol)
    {
        r.failures.push_back("j exceeds limit");
        r.limits_ok = false;
    }

    r.all_pass = r.planned && r.boundary_ok && r.continuity_ok && r.limits_ok;
    return r;
}
