/**
 * Test 1: 1M random SCurve evaluations on the NEW version.
 * All tests must pass. No fixed seeds.
 */
#include "new_s_curve.hpp"
#include "eval.hpp"
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

static constexpr double VM_MIN = 0.1,  VM_MAX = 2000.0;
static constexpr double AM_MIN = 0.5,  AM_MAX = 10000.0;
static constexpr double JM_LO  = 20.0, JM_HI  = 200.0;
static constexpr double POS_MIN = -720.0, POS_MAX = 720.0;

std::mutex print_mtx;

struct Stats { int total, fails; };

static void run_batch(int count, Stats& st)
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

        scurve_new::SCurveProfile::Config cfg{ vm, am, jm };
        scurve_new::SCurveProfile curve(cfg, xs, vs, as, xe, ve, ae);

        auto r = evaluate(curve, xs, xe, vm, am, jm);
        if (!r.all_pass)
        {
            std::lock_guard<std::mutex> lk(print_mtx);
            printf("FAIL: ");
            for (auto& f : r.failures) printf("%s ", f.c_str());
            printf("| xs=%.2f xe=%.2f vs=%.2f ve=%.2f vm=%.1f am=%.1f jm=%.0f\n",
                   xs, xe, vs, ve, vm, am, jm);
            st.fails++;
        }
        st.total++;
    }
}

int main()
{
    printf("SCurve Random Validation (NEW) — %d tests, %d threads\n", TOTAL, THREADS);

    std::vector<std::thread> threads;
    std::vector<Stats> st(THREADS);
    int per = TOTAL / THREADS;

    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < THREADS; i++)
        threads.emplace_back(run_batch, per + (i < TOTAL % THREADS), std::ref(st[i]));
    for (auto& t : threads) t.join();
    auto t1 = std::chrono::steady_clock::now();

    int tot = 0, fails = 0;
    for (auto& s : st) { tot += s.total; fails += s.fails; }

    double dt = std::chrono::duration<double>(t1 - t0).count();
    printf("Done in %.1fs  %d/%d tests  failures: %d\n", dt, tot, TOTAL, fails);
    return fails > 0 ? 1 : 0;
}
