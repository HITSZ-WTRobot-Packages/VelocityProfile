# -*- coding: utf-8 -*-
"""Randomized SCurve validation with multiprocessing.

Generates random parameters, plans S-curves, evaluates them, and reports failures.
"""

from __future__ import annotations

import math
import json
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

from s_curve import SCurve
from s_curve_eval import evaluate, EvalReport

# --- Parameter ranges ---
# 覆盖线位移 (m, m/s) 和角度 (deg, deg/s) 两种量纲

VM_RANGE = (0.1, 2000.0)      # 最大速度
AM_RANGE = (0.5, 10000.0)     # 最大加速度
JM_RANGE = (5.0, 100000.0)    # 最大加加速度
POS_RANGE = (-720.0, 720.0)   # ±两圈

TOTAL_TESTS = 1_000_000
NUM_WORKERS = 28
BATCH_SIZE = 10_000  # tests per worker batch


@dataclass
class TestCase:
    xs: float
    xe: float
    vs: float
    as_: float
    vm: float
    am: float
    jm: float
    ve: float = 0.0
    ae: float = 0.0


@dataclass
class FailureRecord:
    params: dict
    failures: list[str]


def rand_float(lo: float, hi: float) -> float:
    return lo + random.random() * (hi - lo)


def generate_params() -> TestCase:
    """Generate one random set of S-curve parameters."""
    vm = rand_float(*VM_RANGE)
    am = rand_float(*AM_RANGE)
    jm = am * rand_float(20.0, 200.0)

    xs = rand_float(*POS_RANGE)
    xe = rand_float(*POS_RANGE)
    vs = rand_float(-vm, vm)
    ve = rand_float(-vm, vm)
    as_ = rand_float(-am, am) if random.random() < 0.5 else 0.0
    ae = rand_float(-am, am) if random.random() < 0.5 else 0.0

    return TestCase(xs=xs, xe=xe, vs=vs, as_=as_, vm=vm, am=am, jm=jm, ve=ve, ae=ae)


def run_batch(count: int) -> list[FailureRecord]:
    """Run a batch of random tests. Returns list of failures."""
    failures: list[FailureRecord] = []

    for _ in range(count):
        tc = generate_params()

        sc = SCurve()
        ret = sc.init(
            xs=tc.xs, xe=tc.xe,
            vs=tc.vs, as_=tc.as_,
            vm=tc.vm, am=tc.am, jm=tc.jm,
            ve=tc.ve, ae=tc.ae,
        )

        report = evaluate(sc, xs=tc.xs, xe=tc.xe, vm=tc.vm, am=tc.am, jm=tc.jm)

        if not report.all_pass:
            failures.append(FailureRecord(
                params={
                    "xs": tc.xs, "xe": tc.xe,
                    "vs": tc.vs, "as": tc.as_,
                    "vm": tc.vm, "am": tc.am, "jm": tc.jm,
                    "ve": tc.ve, "ae": tc.ae,
                    "success": sc.success,
                    "total_time": sc.total_time,
                },
                failures=report.failures,
            ))

    return failures


def _worker_init():
    """Suppress broken-pipe warnings in worker processes."""
    pass


def main():
    print(f"Starting {TOTAL_TESTS:,} tests with {NUM_WORKERS} workers...")
    t0 = time.time()

    # Distribute work: each worker gets a unique seed and runs BATCH_SIZE tests
    total_batches = TOTAL_TESTS // BATCH_SIZE
    batches_per_worker = total_batches // NUM_WORKERS
    remainder = total_batches % NUM_WORKERS

    tasks: list[int] = []
    for _ in range(total_batches):
        tasks.append(BATCH_SIZE)

    all_failures: list[FailureRecord] = []
    completed = 0

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        for failures in executor.map(run_batch, tasks):
            all_failures.extend(failures)
            completed += BATCH_SIZE
            if completed % (BATCH_SIZE * 10) == 0:
                elapsed = time.time() - t0
                rate = completed / elapsed if elapsed > 0 else 0
                sys.stderr.write(
                    f"\r  {completed:>10,} / {TOTAL_TESTS:,}  "
                    f"({100*completed/TOTAL_TESTS:.1f}%)  "
                    f"failures: {len(all_failures)}  "
                    f"rate: {rate:,.0f}/s"
                )
                sys.stderr.flush()

    elapsed = time.time() - t0
    sys.stderr.write("\n")
    sys.stderr.flush()

    print(f"\nDone in {elapsed:.1f}s ({TOTAL_TESTS/elapsed:,.0f} tests/s)")
    print(f"Total failures: {len(all_failures)} / {TOTAL_TESTS:,} "
          f"({100*len(all_failures)/TOTAL_TESTS:.4f}%)")

    # Print failures
    if all_failures:
        print(f"\n{'='*60}")
        print(f"FAILURES ({len(all_failures)}):")
        print(f"{'='*60}")
        for i, f in enumerate(all_failures):
            print(f"\n--- Failure #{i+1} ---")
            for key, val in f.params.items():
                print(f"  {key}: {val}")
            print(f"  failures: {f.failures}")

    # Save to file
    output_path = "scurve_failures.json"
    if all_failures:
        with open(output_path, "w") as fp:
            json.dump(
                [{"params": f.params, "failures": f.failures} for f in all_failures],
                fp, indent=2,
            )
        print(f"\nFailures saved to {output_path}")
    else:
        print(f"\nNo failures to save.")


if __name__ == "__main__":
    main()
