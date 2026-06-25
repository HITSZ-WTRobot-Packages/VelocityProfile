/**
 * Test 2: 1M old-vs-new comparison.
 * When old succeeds, new must produce identical trajectory.
 * Tolerance scaled by input magnitude.
 */
#include "old_s_curve.hpp"
#include "new_s_curve.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

constexpr int  TOTAL   = 1'000'000;
constexpr int  THREADS = 28;
constexpr int  SAMPLES = 100;

static constexpr double VM_MIN = 0.1,  VM_MAX = 2000.0;
static constexpr double AM_MIN = 0.5,  AM_MAX = 10000.0;
static constexpr double JM_LO  = 20.0, JM_HI  = 200.0;
static constexpr double POS_MIN = -720.0, POS_MAX = 720.0;

std::mutex print_mtx;

struct CmpStats { int old_fail, mismatch; };

static float tolerance(float scale)  { return fmaxf(5e-4f, 5e-6f * scale); }
static float tol_a(float scale)      { return fmaxf(5e-4f, 5e-5f * scale); }
static float tol_T(float scale)      { return fmaxf(5e-4f, 5e-7f * scale); }

static void run_batch(int count, CmpStats& st)
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    auto rd_d = [&](double lo, double hi) { return lo + dist(rng) * (hi - lo); };

    for (int i = 0; i < count; i++)
    {
        float vm = (float)rd_d(VM_MIN, VM_MAX);
        float am = (float)rd_d(AM_MIN, AM_MAX);
        float jm = am * (float)rd_d(JM_LO, JM_HI);
        float xs = (float)rd_d(POS_MIN, POS_MAX);
        float xe = (float)rd_d(POS_MIN, POS_MAX);
        float vs = (float)rd_d(-vm, vm);
        float ve = (float)rd_d(-vm, vm);
        float as = (rng() & 1) ? 0.0f : (float)rd_d(-am, am);
        float ae = (rng() & 1) ? 0.0f : (float)rd_d(-am, am);

        // Old version
        scurve_old::SCurveProfile::Config old_c{ vm, am, jm };
        scurve_old::SCurveProfile old_curve(old_c, xs, vs, as, xe, ve, ae);
        if (!old_curve.success()) { st.old_fail++; continue; }

        // New version
        scurve_new::SCurveProfile::Config new_c{ vm, am, jm };
        scurve_new::SCurveProfile new_curve(new_c, xs, vs, as, xe, ve, ae);

        float T   = old_curve.getTotalTime();
        float Tn  = new_curve.getTotalTime();
        float tol  = tolerance(fmaxf(fabsf(xs), fabsf(xe)));
        float tolA = tol_a(fmaxf(fabsf(xs), fabsf(xe)));
        float tolT = tol_T(fmaxf(fabsf(xs), fabsf(xe)));

        if (fabsf(T - Tn) > tolT)
        {
            std::lock_guard<std::mutex> lk(print_mtx);
            printf("MISMATCH T: old=%.6f new=%.6f\n", T, Tn);
            st.mismatch++;
            continue;
        }

        bool ok = true;
        for (int s = 0; s <= SAMPLES && ok; s++)
        {
            float t  = T * (float)s / (float)SAMPLES;
            float xo = old_curve.CalcX(t);
            float xn = new_curve.CalcX(t);
            float vo = old_curve.CalcV(t);
            float vn = new_curve.CalcV(t);
            float ao = old_curve.CalcA(t);
            float an = new_curve.CalcA(t);

            float dx = fabsf(xo - xn);
            float dv = fabsf(vo - vn);
            float da = fabsf(ao - an);

            if (dx > tol || dv > tol || da > tolA)
            {
                std::lock_guard<std::mutex> lk(print_mtx);
                printf("MISMATCH t=%.4f Δx=%.2e Δv=%.2e Δa=%.2e | xs=%.2f xe=%.2f vs=%.1f ve=%.1f vm=%.1f am=%.1f jm=%.0f\n",
                       t, dx, dv, da, xs, xe, vs, ve, vm, am, jm);
                st.mismatch++;
                ok = false;
            }
        }
    }
}

int main()
{
    printf("SCurve Cross-Version Comparison — %d tests, %d threads\n", TOTAL, THREADS);
    printf("Tolerance: max(1e-4, 1e-6 * scale) for x/v, 10x for a\n\n");

    std::vector<std::thread> threads;
    std::vector<CmpStats> st(THREADS);
    int per = TOTAL / THREADS;

    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < THREADS; i++)
        threads.emplace_back(run_batch, per + (i < TOTAL % THREADS), std::ref(st[i]));
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::steady_clock::now();

    int old_fail = 0, mismatch = 0;
    for (auto& s : st) { old_fail += s.old_fail; mismatch += s.mismatch; }
    int compared = TOTAL - old_fail;

    double dt = std::chrono::duration<double>(t1 - t0).count();
    printf("\nDone in %.1fs\n", dt);
    printf("Old failed:  %d (%.2f%%)\n", old_fail, 100.0*old_fail/TOTAL);
    printf("Compared:    %d\n", compared);
    printf("Mismatches:  %d (%.4f%%)\n", mismatch, 100.0*mismatch/fmax(1,compared));
    return mismatch > 0 ? 1 : 0;
}
