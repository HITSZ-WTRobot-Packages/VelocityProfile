/**
 * Performance benchmark: measures average planning time.
 * Compile twice: once with -DO2_MODE (old + O2), once with -DOFAST_MODE (new + Ofast).
 */
#ifdef O2_MODE
#  include "old_s_curve.hpp"
   namespace sc = scurve_old;
#else
#  include "new_s_curve.hpp"
   namespace sc = scurve_new;
#endif

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>

constexpr int  TOTAL   = 1'000'000;
constexpr int  THREADS = 28;

static constexpr double VM_MIN = 0.1,  VM_MAX = 2000.0;
static constexpr double AM_MIN = 0.5,  AM_MAX = 10000.0;
static constexpr double JM_LO  = 20.0, JM_HI  = 200.0;
static constexpr double POS_MIN = -720.0, POS_MAX = 720.0;

struct BenchStats { double time_ms; int total, success; };

static void run_batch(int count, BenchStats& st)
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

        sc::SCurveProfile::Config cfg{ vm, am, jm };

        auto t0 = std::chrono::steady_clock::now();
        sc::SCurveProfile curve(cfg, xs, vs, as, xe, ve, ae);
        auto t1 = std::chrono::steady_clock::now();

        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        st.time_ms += us * 0.001;
        st.total++;
        if (curve.success()) st.success++;
    }
}

int main()
{
#ifdef O2_MODE
    printf("SCurve Benchmark: OLD + O2\n");
#else
    printf("SCurve Benchmark: NEW + Ofast\n");
#endif
    printf("%d tests, %d threads\n", TOTAL, THREADS);

    std::vector<std::thread> threads;
    std::vector<BenchStats> st(THREADS);
    int per = TOTAL / THREADS;

    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < THREADS; i++)
        threads.emplace_back(run_batch, per + (i < TOTAL % THREADS), std::ref(st[i]));
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::steady_clock::now();

    double total_ms = 0, wall_s = std::chrono::duration<double>(t1 - t0).count();
    int total = 0, success = 0;
    for (auto& s : st) { total += s.total; total_ms += s.time_ms; success += s.success; }

    printf("Total:        %d tests\n", total);
    printf("Success:      %d (%.2f%%)\n", success, 100.0*success/total);
    printf("Avg time:     %.3f us/test\n", total_ms * 1000.0 / total);
    printf("Wall time:    %.1f s (%.0f tests/s)\n", wall_s, total / wall_s);
    return 0;
}
