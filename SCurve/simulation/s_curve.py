# -*- coding: utf-8 -*-

"""
本文件是通过 AI 将 s_curve.c 转换得来
--------------------------------------------------------------------------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


S_CURVE_MAX_BS_ERROR = 0.001


@dataclass(slots=True)
class SidePrepare:
    t_pre: float = 0.0
    x_pre: float = 0.0
    v_base: float = 0.0
    t_shift: float = 0.0
    vp_min: float = 0.0
    valid: bool = True


@dataclass(slots=True)
class FastEvalConfig:
    am: float
    jm: float
    am_square: float
    jerk_ramp_time: float
    inv_double_am: float


@dataclass(slots=True)
class FastEvalSide:
    v_base: float
    x_pre: float
    t_shift: float


@dataclass(slots=True)
class FastEvalProfile:
    has_uniform: bool = False
    v_base: float = 0.0
    vp: float = 0.0
    t1: float = 0.0
    t2: float = 0.0
    x1: float = 0.0
    v1: float = 0.0
    total_time: float = 0.0
    total_distance: float = 0.0


@dataclass(slots=True)
class FastEvalResult:
    dx1: float = 0.0
    dx3: float = 0.0
    delta: float = 0.0


def _build_fast_eval_profile(cfg: FastEvalConfig, v_base: float, vp: float) -> FastEvalProfile:
    profile = FastEvalProfile(v_base=v_base, vp=vp)
    delta_v = vp - v_base

    if cfg.jm * delta_v > cfg.am_square:
        profile.has_uniform = True
        profile.t1 = cfg.jerk_ramp_time
        profile.t2 = delta_v / cfg.am
        profile.v1 = v_base + 0.5 * cfg.am * profile.t1
        profile.x1 = v_base * profile.t1 + (1.0 / 6.0) * cfg.am * profile.t1 * profile.t1
        profile.total_time = profile.t2 + profile.t1
        profile.total_distance = (
            (vp * vp - v_base * v_base) * cfg.inv_double_am
            + 0.5 * (v_base + vp) * profile.t1
        )
        return profile

    peak_acc = math.sqrt(cfg.jm * delta_v)
    profile.has_uniform = False
    profile.t1 = peak_acc / cfg.jm
    profile.t2 = profile.t1
    profile.v1 = v_base + 0.5 * peak_acc * peak_acc / cfg.jm
    profile.x1 = v_base * profile.t1 + (1.0 / 6.0) * peak_acc * profile.t1 * profile.t1
    profile.total_time = 2.0 * profile.t1
    profile.total_distance = (v_base + vp) * profile.t1
    return profile


def _evaluate_fast_eval_distance(
    cfg: FastEvalConfig,
    profile: FastEvalProfile,
    t: float,
) -> float:
    if t <= 0.0:
        return 0.0

    if t < profile.t1:
        return profile.v_base * t + (1.0 / 6.0) * cfg.jm * t * t * t

    if profile.has_uniform and t < profile.t2:
        dt = t - profile.t1
        return profile.x1 + profile.v1 * dt + 0.5 * cfg.am * dt * dt

    if t < profile.total_time:
        dt = profile.total_time - t
        return profile.total_distance - profile.vp * dt + (1.0 / 6.0) * cfg.jm * dt * dt * dt

    return profile.total_distance


def _evaluate_side_distance(cfg: FastEvalConfig, side: FastEvalSide, vp: float) -> float:
    if vp <= side.v_base:
        return side.x_pre

    profile = _build_fast_eval_profile(cfg, side.v_base, vp)
    shift_distance = _evaluate_fast_eval_distance(cfg, profile, side.t_shift)
    return side.x_pre + profile.total_distance - shift_distance


def _evaluate_distance_delta(
    cfg: FastEvalConfig,
    start: FastEvalSide,
    end: FastEvalSide,
    length: float,
    vp: float,
) -> FastEvalResult:
    result = FastEvalResult()
    result.dx1 = _evaluate_side_distance(cfg, start, vp)
    result.dx3 = _evaluate_side_distance(cfg, end, vp)
    result.delta = result.dx1 + result.dx3 - length
    return result


class SCurveAccel:
    def __init__(self):
        self.has_uniform = False
        self.vs = 0.0
        self.jm = 0.0
        self.total_time = 0.0
        self.total_distance = 0.0
        self.t1 = 0.0
        self.x1 = 0.0
        self.v1 = 0.0
        self.t2 = 0.0
        self.ap = 0.0
        self.vp = 0.0

    def init(self, vs: float, vp: float, am: float, jm: float) -> None:
        self.has_uniform = jm * (vp - vs) > am * am
        self.vs = vs
        self.vp = vp
        self.jm = jm

        if self.has_uniform:
            self.ap = am
            self.t1 = am / jm
            self.t2 = (vp - vs) / am
            self.v1 = vs + 0.5 * am * self.t1
            self.x1 = vs * self.t1 + (1.0 / 6.0) * am * self.t1 * self.t1
            self.total_time = self.t2 + self.t1
            self.total_distance = (
                (vp * vp - vs * vs) / (2.0 * am)
                + 0.5 * (vs + vp) * self.t1
            )
            return

        self.ap = math.sqrt(jm * (vp - vs))
        self.t1 = self.ap / jm
        self.t2 = self.t1
        self.v1 = vs + 0.5 * self.ap * self.ap / jm
        self.x1 = vs * self.t1 + (1.0 / 6.0) * self.ap * self.t1 * self.t1
        self.total_time = 2.0 * math.sqrt((vp - vs) / jm)
        self.total_distance = (vs + vp) * math.sqrt((vp - vs) / jm)

    def get_distance(self, t: float) -> float:
        if t <= 0.0:
            return 0.0
        if t < self.t1:
            return self.vs * t + (1.0 / 6.0) * self.jm * t * t * t
        if self.has_uniform and t < self.t2:
            dt = t - self.t1
            return self.x1 + self.v1 * dt + 0.5 * self.ap * dt * dt
        if t < self.total_time:
            dt = self.total_time - t
            return self.total_distance - self.vp * dt + (1.0 / 6.0) * self.jm * dt * dt * dt
        return self.total_distance

    def get_velocity(self, t: float) -> float:
        if t <= 0.0:
            return self.vs
        if t < self.t1:
            return self.vs + 0.5 * self.jm * t * t
        if self.has_uniform and t < self.t2:
            return self.v1 + self.ap * (t - self.t1)
        if t < self.total_time:
            dt = self.total_time - t
            return self.vp - 0.5 * self.jm * dt * dt
        return self.vp

    def get_acceleration(self, t: float) -> float:
        if t <= 0.0:
            return 0.0
        if t < self.t1:
            return self.jm * t
        if self.has_uniform and t < self.t2:
            return self.ap
        if t < self.total_time:
            return self.jm * (self.total_time - t)
        return 0.0

    def __repr__(self) -> str:
        return (
            "SCurveAccel(\n"
            f"    t1: {self.t1}, t2: {self.t2}, total_time: {self.total_time},\n"
            f"    x1: {self.x1}, total_distance: {self.total_distance},\n"
            f"    vs: {self.vs}, v1: {self.v1}, vp: {self.vp},\n"
            f"    ap: {self.ap}, has_uniform: {self.has_uniform}\n"
            ")"
        )


class SCurve:
    S_CURVE_FAILED = 0
    S_CURVE_SUCCESS = 1

    def __init__(self):
        self.process1 = SCurveAccel()
        self.process3 = SCurveAccel()
        self._reset_state()

    def _reset_state(self) -> None:
        self.success = False
        self.has_const = False
        self.direction = 1.0
        self.vp = 0.0
        self.vs = 0.0
        self.as_ = 0.0
        self.ve = 0.0
        self.ae = 0.0
        self.jm = 0.0

        self.xs = 0.0
        self.xe = 0.0
        self.t1_pre = 0.0
        self.x1_pre = 0.0
        self.ts1 = 0.0
        self.xs1 = 0.0

        self.vrs = 0.0
        self.ars = 0.0
        self.t3_pre = 0.0
        self.x3_pre = 0.0
        self.ts3 = 0.0
        self.xs3 = 0.0

        self.x1 = 0.0
        self.x2 = 0.0

        self.t1 = 0.0
        self.t2 = 0.0
        self.total_time = 0.0

    def __repr__(self) -> str:
        return (
            "SCurve(\n"
            f"    success: {self.success}, has_const: {self.has_const}, direction: {self.direction},\n"
            f"    vp: {self.vp}, xs: {self.xs}, x1: {self.x1}, x2: {self.x2}, xe: {self.xe},\n"
            f"    ts1: {self.ts1}, ts3: {self.ts3}, t1: {self.t1}, t2: {self.t2}, total_time: {self.total_time},\n"
            f"    process1: {self.process1},\n"
            f"    process3: {self.process3},\n"
            ")"
        )

    def init(
        self,
        xs: float,
        xe: float,
        vs: float,
        as_: float,
        vm: float,
        am: float,
        jm: float,
        ve: float = 0.0,
        ae: float = 0.0,
    ) -> int:
        self._reset_state()
        self.process1 = SCurveAccel()
        self.process3 = SCurveAccel()

        vm = abs(vm)
        am = abs(am)
        jm = abs(jm)

        direction = 1.0 if xe > xs else -1.0
        length = abs(xe - xs)

        vs *= direction
        as_ *= direction
        ve *= direction
        ae *= direction

        self.direction = direction
        self.xs = xs
        self.xe = xe
        self.vs = vs
        self.as_ = as_
        self.ve = ve
        self.ae = ae
        self.jm = jm
        self.vrs = ve
        self.ars = -ae

        if length < 1e-6:
            if abs(vs - ve) > 1e-3 or abs(as_ - ae) > 1e-3:
                return self.S_CURVE_FAILED

            self.x1_pre = xs
            self.x1 = xs
            self.x2 = xs
            self.success = True
            return self.S_CURVE_SUCCESS

        if abs(vs) > vm or abs(as_) > am or abs(ve) > vm or abs(ae) > am:
            return self.S_CURVE_FAILED

        def prepare_side(v0: float, a0: float) -> SidePrepare:
            result = SidePrepare()

            if a0 < 0.0:
                result.v_base = v0 - 0.5 * a0 * a0 / jm
                if abs(result.v_base) > vm:
                    result.valid = False
                    return result

                result.vp_min = result.v_base
                result.t_pre = -a0 / jm
                result.x_pre = v0 * result.t_pre + (1.0 / 3.0) * a0 * result.t_pre * result.t_pre
                result.t_shift = 0.0
            else:
                result.vp_min = v0 + 0.5 * a0 * a0 / jm
                if vm < result.vp_min:
                    result.valid = False
                    return result

                result.t_pre = 0.0
                result.x_pre = 0.0
                result.t_shift = a0 / jm
                result.v_base = v0 - 0.5 * a0 * result.t_shift

            if result.vp_min < 0.0:
                result.vp_min = 0.0
            return result

        start = prepare_side(vs, as_)
        if not start.valid:
            return self.S_CURVE_FAILED

        endr = prepare_side(ve, -ae)
        if not endr.valid:
            return self.S_CURVE_FAILED

        self.t1_pre = start.t_pre
        self.x1_pre = xs + direction * start.x_pre
        self.ts1 = start.t_shift
        self.t3_pre = endr.t_pre
        self.x3_pre = endr.x_pre
        self.ts3 = endr.t_shift

        vp_min = max(start.vp_min, endr.vp_min)
        if vp_min < 0.0:
            vp_min = 0.0
        if vm < vp_min:
            return self.S_CURVE_FAILED

        len0 = length - start.x_pre - endr.x_pre
        if len0 < -S_CURVE_MAX_BS_ERROR:
            return self.S_CURVE_FAILED

        fast_eval_cfg = FastEvalConfig(
            am=am,
            jm=jm,
            am_square=am * am,
            jerk_ramp_time=am / jm,
            inv_double_am=0.5 / am,
        )
        start_eval = FastEvalSide(start.v_base, start.x_pre, start.t_shift)
        end_eval = FastEvalSide(endr.v_base, endr.x_pre, endr.t_shift)

        vm_eval = _evaluate_distance_delta(fast_eval_cfg, start_eval, end_eval, length, vm)
        x_const = -vm_eval.delta

        if x_const > 0.0:
            self.process1.init(start.v_base, vm, am, jm)
            self.process3.init(endr.v_base, vm, am, jm)
            self.xs1 = self.process1.get_distance(self.ts1)
            self.xs3 = self.process3.get_distance(self.ts3)

            self.has_const = True
            self.t1 = self.t1_pre + self.process1.total_time - self.ts1

            t_const = x_const / vm
            self.t2 = self.t1 + t_const
            self.total_time = self.t2 + self.t3_pre + self.process3.total_time - self.ts3

            self.x1 = self.xs + self.direction * vm_eval.dx1
            self.x2 = self.x1 + self.direction * x_const
            self.vp = vm
            self.success = True
            return self.S_CURVE_SUCCESS

        left = vp_min
        right = vm
        delta_d = len0
        dx1 = 0.0
        dx3 = 0.0

        while right - left > S_CURVE_MAX_BS_ERROR:
            mid = 0.5 * (left + right)
            mid_eval = _evaluate_distance_delta(fast_eval_cfg, start_eval, end_eval, length, mid)
            delta_d = mid_eval.delta

            if abs(delta_d) <= S_CURVE_MAX_BS_ERROR:
                left = mid
                right = mid
                break

            if delta_d > 0.0:
                right = mid
            else:
                left = mid

        vp = 0.5 * (left + right)
        self.process1.init(start.v_base, vp, am, jm)
        self.process3.init(endr.v_base, vp, am, jm)
        self.xs1 = self.process1.get_distance(self.ts1)
        self.xs3 = self.process3.get_distance(self.ts3)
        dx1 = start.x_pre + self.process1.total_distance - self.xs1
        dx3 = endr.x_pre + self.process3.total_distance - self.xs3
        delta_d = dx1 + dx3 - length

        if delta_d > S_CURVE_MAX_BS_ERROR:
            return self.S_CURVE_FAILED

        self.has_const = False
        self.t1 = self.t1_pre + self.process1.total_time - self.ts1
        self.t2 = self.t1
        self.total_time = self.t2 + self.t3_pre + self.process3.total_time - self.ts3

        self.x1 = self.xs + self.direction * dx1
        self.x2 = self.x1
        self.vp = vp
        self.success = True
        return self.S_CURVE_SUCCESS

    def _get_reverse_distance(self, tau: float) -> float:
        if tau <= 0.0:
            return 0.0
        if tau < self.t3_pre:
            tau2 = tau * tau
            tau3 = tau2 * tau
            return self.ve * tau - 0.5 * self.ae * tau2 + (1.0 / 6.0) * self.jm * tau3

        return self.x3_pre + self.process3.get_distance(tau - self.t3_pre + self.ts3) - self.xs3

    def _get_reverse_velocity(self, tau: float) -> float:
        if tau <= 0.0:
            return self.ve
        if tau < self.t3_pre:
            return self.ve - self.ae * tau + 0.5 * self.jm * tau * tau

        return self.process3.get_velocity(tau - self.t3_pre + self.ts3)

    def _get_reverse_acceleration(self, tau: float) -> float:
        if tau <= 0.0:
            return -self.ae
        if tau < self.t3_pre:
            return -self.ae + self.jm * tau

        return self.process3.get_acceleration(tau - self.t3_pre + self.ts3)

    def calc_x(self, t: float) -> float:
        if not self.success:
            return 0.0
        if t <= 0.0:
            return self.xs
        if t < self.t1_pre:
            t2 = t * t
            t3 = t2 * t
            return self.xs + self.direction * (
                self.vs * t + 0.5 * self.as_ * t2 + (1.0 / 6.0) * self.jm * t3
            )
        if t < self.t1:
            return self.x1_pre + self.direction * (
                self.process1.get_distance(t - self.t1_pre + self.ts1) - self.xs1
            )
        if self.has_const and t < self.t2:
            return self.x1 + self.direction * self.vp * (t - self.t1)
        if t < self.total_time:
            return self.xe - self.direction * self._get_reverse_distance(self.total_time - t)
        return self.xe

    def calc_v(self, t: float) -> float:
        if not self.success:
            return 0.0
        if t <= 0.0:
            return self.direction * self.vs
        if t < self.t1_pre:
            return self.direction * (self.vs + self.as_ * t + 0.5 * self.jm * t * t)
        if t < self.t1:
            return self.direction * self.process1.get_velocity(t - self.t1_pre + self.ts1)
        if self.has_const and t < self.t2:
            return self.direction * self.vp
        if t < self.total_time:
            return self.direction * self._get_reverse_velocity(self.total_time - t)
        return self.direction * self.ve

    def calc_a(self, t: float) -> float:
        if not self.success:
            return 0.0
        if t <= 0.0:
            return self.direction * self.as_
        if t < self.t1_pre:
            return self.direction * (self.as_ + self.jm * t)
        if t < self.t1:
            return self.direction * self.process1.get_acceleration(t - self.t1_pre + self.ts1)
        if self.has_const and t < self.t2:
            return 0.0
        if t < self.total_time:
            return -self.direction * self._get_reverse_acceleration(self.total_time - t)
        return self.direction * self.ae
