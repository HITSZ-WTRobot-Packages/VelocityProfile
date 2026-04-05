# 带初末状态约束的 S 形曲线路径规划

## 问题定义

在区间 $[x_s,x_e]$ 上规划最短时间的运动轨迹 $x(t)$，要求：位置、速度、加速度连续；并受约束：

- 最大速度 $|v|\le v_m$；最大加速度 $|a|\le a_m$；最大加加速度（jerk）$|j|\le j_m$。
- 边界条件：初始 $v(0)=v_s,\ a(0)=a_s$；终止 $v(T)=v_e,\ a(T)=a_e$。

本文在原有只支持 $v_e=0,a_e=0$ 的实现上做扩展，使其兼容任意满足上界约束的终端速度与加速度。

## 1 原始算法回顾（简要）

核心思想：把运动方向统一为正向（记 $\mathrm{dir}=\mathrm{sign}(x_e-x_s)$），在归一化坐标系内求解，再乘回方向。核心模块为 `SCurveAccel`，它给出从“加速度为 0 的初速 $v_{in}$”到“加速度为 0 的峰值 $v_p$”的最短时间分段解。

`SCurveAccel` 的两类解：

- 有匀加速段（3 段）：当 $j_m(v_p-v_{in})>a_m^2$，先 jerk 增加加速度到 $a_m$，再匀加速，最后减小加速度。
- 无匀加速段（2 段）：当增速不足以达到 $a_m$，最大加速度受 jerk 限制，前后对称。

旧实现处理初始加速度 $a_s$：若 $a_s<0$（与运动方向相反），先做一段把加速度抬到 0 的预处理；若 $a_s\ge0$，通过时间平移把带初加速度的情况嵌入 `SCurveAccel`。

旧实现把减速段当作“从 0 加速到 $v_p$ 的逆过程”，因此只能满足终态 $v_e=0,a_e=0$。

## 2 扩展思路 — 用时间反演处理终端约束

关键：把终端约束通过时间反演（mirror）转换成“反向起点约束”，然后用与起点相同的单侧处理逻辑。

令 $$\tau=T-t$$ 为反演时间，则反演域的初始条件为
$$
v_{rs}=v_e,
\qquad a_{rs}=-a_e
$$
（注意 $a$ 在时间反演下符号改变），此时反演域内的“加速到 $v_p$”过程与起点侧用同一套 `SCurveAccel` 逻辑求解。

## 3 单侧预处理（起点/反演端通用）

对任一侧边界 $(v_0,a_0)$，定义 `prepareSide(v_0,a_0)`：

- 若 $a_0<0$（与该侧正向过程相反），必须先把加速度抬到 0：
	$$
	t_{pre}=-\frac{a_0}{j_m},\qquad v_{base}=v_0-\frac{a_0^2}{2j_m}
	$$
	预处理位移
	$$
	x_{pre}=v_0 t_{pre}+\frac{1}{3}a_0 t_{pre}^2
	$$

- 若 $a_0\ge0$（同向），可以通过时间平移嵌入标准解：
	$$
	t_{shift}=\frac{a_0}{j_m},\qquad v_{base}=v_0-\frac{a_0^2}{2j_m}
	$$
	此外该侧对峰值速度的下界为
	$$
	v_{p,min}=v_0+\frac{a_0^2}{2j_m}.
	$$

处理后得到：该侧等效起始速度 `v_base`、预处理位移 `x_pre`、时间平移 `t_shift` 与最低可行峰值 `v_{p,min}`。

## 4 长度匹配与二分求解

设总路程 $L=|x_e-x_s|$，起点侧贡献为
$$
\Delta x_s(v_p)=x_{pre,s}+\big(S_s(v_p)-S_s(t_{shift,s})\big)
$$
其中 $S_s(\cdot)$ 表示 `SCurveAccel` 的位移（`getDistance`），类似地终点反演侧贡献为
$$
\Delta x_e(v_p)=x_{pre,e}+\big(S_e(v_p)-S_e(t_{shift,e})\big)
$$

若存在匀速段则取 $v_p=v_m$，并令
$$
L=\Delta x_s(v_m)+\Delta x_e(v_m)+x_{const},\qquad x_{const}>0
$$
无匀速段时，需求根
$$
F(v_p)=\Delta x_s(v_p)+\Delta x_e(v_p)-L=0
$$
在区间
$$
v_p\in[\max(v_{p,min,s},v_{p,min,e}),\ v_m]
$$
用二分（或其它单峰算法）求解；代码中以二分为准并用容差 $\varepsilon$ 截止。

## 5 时间拼接与总时长

对每侧定义加速器输出总时长 $T_{acc}$ 与时间偏移后有效时长
$$
T_s=t_{pre,s}+T_{acc,s}-t_{shift,s},\qquad T_e=t_{pre,e}+T_{acc,e}-t_{shift,e}
$$
若存在匀速段额外时长 $T_c=x_{const}/v_p$，则总时长
$$
T=T_s+T_c+T_e
$$

## 6 轨迹求值（在代码中的实现要点）

- 对 $t\le0$ 返回初始 $x_s,v_s,a_s$（已乘方向）。
- 起点预处理段按多项式直接计算（jerk 常数段）。
- 加速段使用 `SCurveAccel.getDistance/Velocity/Acceleration`，并考虑时间平移 `t_shift`。
- 减速段通过反演函数计算（代码中的 `getReverseDistance/Velocity/Acceleration`）：令 $\tau=T-t$，
	$$
	x(t)=x_e-\mathrm{dir}\cdot X_r(\tau),\quad v(t)=\mathrm{dir}\cdot V_r(\tau),\quad a(t)=-\mathrm{dir}\cdot A_r(\tau)
	$$
	如此保证在 $t=T$ 时满足 $v(T)=v_e,\ a(T)=a_e$。

## 7 可行性与失败条件

实现中应返回失败（不可构造）的情形包括：

- 初末速度或加速度超过各自上限。
- 单侧预处理后速度仍超限。
- 两侧预处理位移之和已超过总路程 $L$。
- 二分或求解结束后仍无法在容差内匹配位移。

## 8 与旧版的兼容性

当 $v_e=0,a_e=0$ 时，反演侧退化为原先固定的“从 0 到 $v_p$ 的逆过程”，新版严格向后兼容旧行为，旧的调用无需改动（构造函数保持默认参数）。

---



