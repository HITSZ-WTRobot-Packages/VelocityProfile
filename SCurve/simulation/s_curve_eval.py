# -*- coding: utf-8 -*-
"""SCurve trajectory evaluation module.

Evaluates a planned S-curve against its input constraints and checks
correctness: planning success, boundary conditions, continuity, and limits.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

EPS = 1.0e-6
BOUNDARY_EPS = 1.0e-5  # relaxed tolerance for x(0)/x(T) boundary checks

# Multiplier for discontinuity detection: a jump larger than MULT * max_expected_change
# is considered a discontinuity.  max_expected_change = limit * dt.
DISCONTINUITY_MULTIPLIER = 5.0

# Jerk limit tolerance: numerical differentiation amplifies noise.
JERK_LIMIT_TOLERANCE = 0.01


@dataclass
class EvalReport:
    """Result of evaluating an S-curve trajectory."""

    planned: bool = False
    boundary_ok: bool = False
    continuity_ok: bool = False
    limits_ok: bool = False
    all_pass: bool = False

    failures: list[str] = field(default_factory=list)

    max_v: float = 0.0
    max_a: float = 0.0
    max_j: float = 0.0

    def __repr__(self) -> str:
        lines = [
            "=== SCurve Eval Report ===",
            f"  planned:       {self.planned}",
            f"  boundary_ok:   {self.boundary_ok}",
            f"  continuity_ok: {self.continuity_ok}",
            f"  limits_ok:     {self.limits_ok}",
            f"  => all_pass:   {self.all_pass}",
            f"  max_v: {self.max_v:.4f}  max_a: {self.max_a:.4f}  max_j: {self.max_j:.4f}",
        ]
        if self.failures:
            lines.append("  failures:")
            for f in self.failures:
                lines.append(f"    - {f}")
        return "\n".join(lines)


def evaluate(
    curve,       # SCurve instance
    xs: float,
    xe: float,
    vm: float,
    am: float,
    jm: float,
    sample_dt: float = 0.001,
) -> EvalReport:
    """Evaluate a planned S-curve trajectory.

    Args:
        curve: An SCurve instance (must have calc_x, calc_v, calc_a, success, total_time).
        xs, xe: Expected start/end positions.
        vm, am, jm: Limits (max speed, max acceleration, max jerk).
        sample_dt: Sampling interval in seconds (default 0.001).

    Returns:
        EvalReport with all check results.
    """
    report = EvalReport()

    # --- 1. Planning success ---
    report.planned = curve.success
    if not curve.success:
        report.failures.append("planning failed (curve.success is False)")
        report.all_pass = False
        # Continue with remaining checks on whatever data is available.
        # If total_time is 0 or negative, we can't sample — stop here.
        if curve.total_time <= 0.0:
            return report

    T = curve.total_time

    # Zero-time trajectory: only check planning success and boundary
    if T < EPS:
        report.boundary_ok = True
        report.continuity_ok = True
        report.limits_ok = True
        report.all_pass = report.planned
        return report

    # --- Sample trajectory at fixed interval ---
    dt = sample_dt
    n_samples = int(T / dt) + 1
    if n_samples < 2:
        n_samples = 2
        dt = T
    times = [i * dt for i in range(n_samples)]
    # Ensure last sample is exactly at T
    if times[-1] < T - EPS:
        times.append(T)

    xs_vals = [curve.calc_x(t) for t in times]
    vs_vals = [curve.calc_v(t) for t in times]
    as_vals = [curve.calc_a(t) for t in times]

    # Jerk from numerical differentiation of acceleration
    js_vals: list[float] = []
    for i in range(len(times) - 1):
        da = as_vals[i + 1] - as_vals[i]
        js_vals.append(da / dt)
    js_vals.append(0.0)  # last point, no next sample

    # --- 2. Boundary conditions ---
    bound_eps = max(BOUNDARY_EPS, 3e-5 * max(abs(xs), abs(xe)))
    x0_ok = abs(xs_vals[0] - xs) < bound_eps
    xT_ok = abs(xs_vals[-1] - xe) < bound_eps
    report.boundary_ok = x0_ok and xT_ok
    if not x0_ok:
        report.failures.append(f"x(0)={xs_vals[0]:.6f} != xs={xs}")
    if not xT_ok:
        report.failures.append(f"x(T)={xs_vals[-1]:.6f} != xe={xe}")

    # --- 3. Continuity ---
    # Expected maximum change per sample under the given limits
    max_dx = vm * dt
    max_dv = am * dt
    max_da = jm * dt
    thresh_x = max(1e-3, DISCONTINUITY_MULTIPLIER * max_dx)
    thresh_v = max(0.1, DISCONTINUITY_MULTIPLIER * max_dv, 0.02 * vm)
    thresh_a = max(0.5, DISCONTINUITY_MULTIPLIER * max_da)

    continuity_failures: list[str] = []
    for i in range(len(times) - 1):
        dx = abs(xs_vals[i + 1] - xs_vals[i])
        dv = abs(vs_vals[i + 1] - vs_vals[i])
        da = abs(as_vals[i + 1] - as_vals[i])

        if dx > thresh_x:
            continuity_failures.append(f"x discontinuity at t≈{times[i]:.6f}: Δx={dx:.6f}")
            break
        if dv > thresh_v:
            continuity_failures.append(f"v discontinuity at t≈{times[i]:.6f}: Δv={dv:.4f}")
        if da > thresh_a:
            continuity_failures.append(f"a discontinuity at t≈{times[i]:.6f}: Δa={da:.4f}")

    report.continuity_ok = len(continuity_failures) == 0
    report.failures.extend(continuity_failures)

    # --- 4. Limits ---
    max_v = max(abs(v) for v in vs_vals)
    max_a = max(abs(a) for a in as_vals)
    max_j = max(abs(j) for j in js_vals)

    report.max_v = max_v
    report.max_a = max_a
    report.max_j = max_j

    limit_failures: list[str] = []
    if max_v > vm + EPS:
        limit_failures.append(f"v exceeds limit: max|v|={max_v:.4f} > vm={vm}")
    if max_a > am + EPS:
        limit_failures.append(f"a exceeds limit: max|a|={max_a:.4f} > am={am}")
    if max_j > jm + JERK_LIMIT_TOLERANCE:
        limit_failures.append(f"j exceeds limit: max|j|={max_j:.4f} > jm={jm}")

    report.limits_ok = len(limit_failures) == 0
    report.failures.extend(limit_failures)

    # --- Final ---
    report.all_pass = (
        report.planned
        and report.boundary_ok
        and report.continuity_ok
        and report.limits_ok
    )

    return report
