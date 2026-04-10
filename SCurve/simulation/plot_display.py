# -*- coding: utf-8 -*-
from s_curve import SCurve

import numpy as np
import matplotlib.pyplot as plt


def plot_scurve(xs, xe, vs, as_, vm, am, jm, ve=0.0, ae=0.0, dt=0.001):
    s = SCurve()
    if s.init(xs, xe, vs, as_, vm, am, jm, ve, ae) == s.S_CURVE_FAILED:
        print("S curve init failed")
        return

    total_time = s.total_time

    if total_time <= 0.0:
        t = np.array([0.0])
    else:
        t = np.arange(0, total_time + dt, dt)

    x = np.array([s.calc_x(ti) for ti in t])
    v = np.array([s.calc_v(ti) for ti in t])
    a = np.array([s.calc_a(ti) for ti in t])
    if len(t) >= 2:
        j = np.gradient(a, t, edge_order=1)
    else:
        j = np.zeros_like(a)

    print("===== SCurve quick check =====")
    print(f"total_time: {total_time:.6f}s")
    print(f"peak_v: {s.vp:.6f}, has_const: {s.has_const}")
    print(f"x(0)/x(T): {x[0]:.6f} / {x[-1]:.6f} (target {xs:.6f} / {xe:.6f})")
    print(f"v(0)/v(T): {v[0]:.6f} / {v[-1]:.6f} (target {vs:.6f} / {ve:.6f})")
    print(f"a(0)/a(T): {a[0]:.6f} / {a[-1]:.6f} (target {as_:.6f} / {ae:.6f})")
    print(f"max|v|={np.max(np.abs(v)):.6f} <= vm({vm:.6f})")
    print(f"max|a|={np.max(np.abs(a)):.6f} <= am({am:.6f})")
    print(f"max|j|≈{np.max(np.abs(j)):.6f} <= jm({jm:.6f}) [numerical]")

    # x(t)
    plt.figure()
    plt.plot(t, x)
    plt.title("x(t)")
    plt.xlabel("t")
    plt.ylabel("x")
    plt.grid()

    # v(t)
    plt.figure()
    plt.plot(t, v)
    plt.title("v(t)")
    plt.xlabel("t")
    plt.ylabel("v")
    plt.grid()

    # a(t)
    plt.figure()
    plt.plot(t, a)
    plt.title("a(t)")
    plt.xlabel("t")
    plt.ylabel("a")
    plt.grid()

    plt.figure()
    plt.plot(t, j)
    plt.title("j(t)")
    plt.xlabel("t")
    plt.ylabel("j")
    plt.grid()

    plt.show()

if __name__ == "__main__":
    plot_scurve(
        xs=0.0,
        xe=5.0,
        vs=1.0,
        as_=-1.0,
        vm=10.0,
        am=1.2,
        jm=2.4,
        ve=1.0,
        ae=-1.0,
    )
