# LQMDini.for — 兰州量子分子动力学（LQMD）初始化程序文档

## 1. 程序概述

`LQMDini.for`（5837 行）是 LQMD 输运模型的**初始化与基态制备程序**。其核心功能是：

- 为孤立原子核生成稳定的基态相空间配置（位置 + 动量）
- 通过费米关联动量 `pfc` 的迭代调整，使结合能与实验值匹配
- 通过 RMS 半径漂移检测剔除不稳定配置
- 输出稳定配置供动力学主程序 (`main.for`) 使用
- 摩擦冷却方法采用的是线性摩擦冷却

## 2. 输入文件与参数

主输入文件为 `LQMDin.in`，读取于第 65-82 行：

| 行 | 变量 | 含义 | 典型值 (Th-232) |
|----|------|------|-----------------|
| 1 | `npr, nzpr` | 总核子数、质子数 | 232, 90 |
| 2 | `progEb, beta2` | 目标结合能 (MeV/u)、四极形变 | -7.615, 0.0 |
| 3 | `alpha, beta, gamma` | Skyrme 两体/三体强度、密度幂次 | -0.2825, 0.220, 1.1667 |
| 4 | `gsur, gtau, eta, gsym` | 表面项、有效质量项、eta、表面对称项 | 0.0120, 0.0097, 1.6666, -0.0027 |
| 5 | `symm, rho0, c1, c2` | 对称能强度、饱和密度(fm⁻³)、波包宽度系数 | 0.032, 0.16, 0.08, 0.92 |
| 6 | `asym, bsym, csym, gamsy` | MDI 对称能参数 (nsym=3 时启用) | 0.043, -0.01675, 0., 2. |
| 7 | `vpa, vpb, vpc, q01, p01, e01` | Pauli 势参数 (通常全 0 = 关闭) | 0., 0., 0., 1.6, 0.12, 30 |
| 8 | `ish, rmax, r0, ds` | 壳修正开关(0=液滴)、RMS 容差(fm)、半径参数(fm)、弥散度(fm) | 0, 1.0, 1.16, 0.65 |
| 9 | `initi, nruns` | 初始化模式(0=新生成)、独立运行次数 | 0, 500 |
| 10 | `nt, dt` | 时间步数、步长 (fm/c) | 200, 0.5 |
| 11 | `Ri0, rdt` | 熔合势垒初始间距(fm)、递减步长(仅 initi=1 使用) | 25., 0.05 |
| 12 | `cdelt, cmeps, cmom, xsym, momp` | 动量相关势参数 | 0.00176, 500., -0.0848, 0., 1 |
| 13 | `nsym, xmu` | 对称能类型(1/2/3)、nsym=2 时的密度幂次 | 1, 0.5 |
| 14 | `gamma_fric, tcool` | 摩擦冷却强度 (c/fm)、冷却时长 (fm/c) | 0.02, 2000.0 |
| 14a | `gamma_fric_R` | 位置阻尼强度 (c/fm)，0=关闭 | 0.0 |
| 15 | `iniout` | 稳定 run 相空间输出文件名 | tini.dat |
| 16 | `rmsout` | RMS 半径演化输出文件名 | rms.dat |
| 17 | `aveout` | 平均结合能演化输出文件名 | ave.dat |
| 18 | `dtiniout` | 每时间步核子坐标输出文件名 | dtini.dat |

### 重要说明

- **`rmax`**（第8行）是**核心稳定性判据**：当核的 RMS 半径偏离初始值超过 `rmax` fm 时，该 run 被标记为不稳定（`lswitch=0`）并中止演化。
- **`progEb`**（第2行）是结合能目标值，`init0` 会通过调节 `pfc` 迭代匹配到此值。
- **Skyrme 参数**（第3-6行）必须与动力学程序 `LQMD.in` 的对应参数完全一致，否则初始化时找到的基态在动力学演化中会成为激发态 → 碎裂。

## 3. 主程序流程

```
program LQMDini
│
├─ 读取 LQMDin.in 参数 (行 65-82)
├─ 初始化常数、波包宽度、库仑表 (行 84-358)
│
└─ do nru = 1, nruns                    ← 外层：独立运行统计循环
   │
   ├─ 设置随机数种子
   ├─ call init0                        ← 生成并冷却基态组态
   │
   ├─ do k = 1, nt                      ← 内层：时间演化循环 (0 → nt*dt fm/c)
   │   │
   │   ├─ 写入 dtini.dat (每核子一行)    ← rtime, x, y, z, iso
   │   ├─ call updtim / updtsh / constr  ← 传播位置与动量（含位置阻尼）
   │   ├─ call prot                      ← 计算 RMS → 实时写入 rms.dat
   │   │   └─ 稳定性判据: |rms - radj| > rmax ? → lswitch=0 (不稳定)
   │   ├─ call avene                     ← 计算结合能 → 实时写入 ave.dat
   │   │
   │   ├─ if(lswitch==0) 中止演化       ← 不稳定：跳出时间循环
   │   │
   │   └─ if(k==nt)                     ← 最终时间步：稳定！
   │       ├─ 写入 0result.inf:  "Running NO: N is stable !!"
   │       ├─ 写入 tini.dat:  相空间坐标 (x,y,z,px,py,pz,em,iso,isy,id)
   │       ├─ 写入 rms.dat:   RMS 半径时间演化
   │       └─ 写入 ave.dat:   平均结合能时间演化
   │
   └─ 关闭文件，结束
```

## 4. 核心子程序详解

### 4.1 `init0`（行 2107-2138）— 结合能调谐主控

```
1. 设置初始 pfc = 0.92 * sqrt((ρ₀/2)^(2/3) / (2m₀))
2. call init1   → 生成核子位置 + 费米动量
3. call avene  → 计算结合能 Eb
4. if Eb < progEb - 0.2 → pfc -= 0.002  (降低截断 = 更多高动量 = 更大动能 = 更浅结合)
5. if Eb > progEb + 0.2 → pfc += 0.002  (提高截断 = 更少高动量 = 更小动能 = 更深结合)
6. 重复直至 |Eb - progEb| < 0.2 MeV/u 或超过 140 次迭代
7. 最终以 laesw=1 调用 avene，输出能量分解到 007E.out
```

**物理本质**：`pfc` 是费米面以下的动量截断。调节 `pfc` 等价于调节核子的平均动能，从而调节总结合能。这是 QMD 中初始化稳定基态的核心技巧。

### 4.2 `init1`（行 1560-2104）— 完整基态组态生成

流程：

| 步骤 | 行号 | 操作 |
|------|------|------|
| 1 | 1707-1763 | **位置采样**：均匀球内随机采样，半径 `r1 = 1.28*A^(1/3) - 0.76 + 0.8*A^(-1/3)`，拒绝距离太近的核子对 (nn<1.5fm, np<0.5fm) |
| 2 | 1715 | **形变应用**：`r *= (1 + β₂(3cos²θ-1)/2)`，四极形变 |
| 3 | 1769-1790 | **RMS 检验**：RMS 半径与目标 RMS0N/RMS0P 比较，偏差>0.3fm 则重采 |
| 4 | 1848 | `call densitn` → 计算每个核子处的局域密度 `denstl(i)` |
| 5 | 1861-1909 | **费米动量采样**：`p_F(r) = ħc(3ρ(r)/(8π))^(1/3)`（局域 Thomas-Fermi），动量大小 `pr = (p_F - pfc) * ran^(1/3)`，方向均匀随机 |
| 6 | 1902-1947 | **Pauli 检验**：检查 6D 相空间占用，违反则重采 |
| 7 | 1948-1950 | `call delrot` → 减去集体转动角动量 |
| 8 | 1956-2044 | **CM 归零**：质心平移到原点，存储到输出数组 `xi/yi/zi/pxi/pyi/pzi` |

**关键的硬编码参数**：
- `DISTNN = 1.5` fm：两个中子间的最小距离
- `DISTNP = 0.5` fm：中子-质子间的最小距离
- `PHMANN = 89.5`：中子-中子 Pauli 阻塞的相空间阈值（6D）
- `PHMANP = 44.8`：中子-质子 Pauli 阻塞的相空间阈值（6D）
- `RMS0N/RMS0P`：变形校正后的中子/质子目标 RMS 半径

### 4.3 `prot`（行 2361-2612）— 诊断与稳定性判据

**稳定性判据**（行 2490-2496）：
```fortran
if(iflag1 < 0 .and. epp == 0 .and. abs(rms - radj) > rmax) then
    lswitch = 0   ! 不稳定：RMS 半径漂移超出容差
else
    lswitch = 1   ! 稳定
end if
```

- `rms` = 当前 RMS 半径
- `radj` = 初始目标 RMS 半径（由 `RMS0N`/`RMS0P` 计算）
- `rmax` = 容差（从 `LQMDin.in` 第8行读取，当前为 1.0 fm）

当 `lswitch=0` 且无碰撞动力学（`epp=0`）时，时间演化立即中止，该 run 被跳过，不写入 `tini.dat`。

### 4.4 `avene`（行 4480-4793）— 结合能/势能完整计算

计算孤立核子的所有能量分量，返回 `aene`（MeV/u）和 `tene`（GeV）：

| 能量分量 | 公式 | 代码变量 |
|----------|------|----------|
| 动能 | Σ p²/(2m) + Thomas-Fermi 修正 | `tkine` |
| Skyrme 两体 | (α/2) Σ ρᵢ | `utwo` |
| Skyrme 三体 | β/(γ+1) Σ ρᵢ^γ | `uthr` |
| 库仑直接 | 误差函数核电荷分布积分 | `ucou` |
| 库仑交换 | Cₑₓ·Z^(4/3) | `couex` |
| 对称能 | nsym=1/2/3 三种选项 | `usym` |
| 表面梯度 | gsur · Σ (6/λₗ - 4r²ᵢⱼ/λₗ²) | `ugradi` |
| 表面对称 | gsym · Σ(iso₁-iso₂)² · f(rᵢⱼ) | `ussr` |
| 有效质量 | gtau · Σ ρᵢ^η | `utau` |
| 动量相关 | momp=0/1/2 (Aichelin/Yukawa) | `umomd` |
| Pauli 势 | 高斯 + Fermi-Dirac 形式 | `upaul` |
| 壳修正 | Woods-Saxon 型壳修正势 | `Esh` |

**总能量**：`tene = tkine + usky + ucou + usym + ussr + ugradi + utau + couex + umomd + upaul + Esh`

输出写入 `007E.out`（当 `laesw=1` 时），包含每一步的各能量分量。

### 4.5 `updtim`（行 2995-3509）— 完整 LQMD 时间演化

单步传播核子相空间坐标：

1. 计算速度（含动量相关修正）：`vel = p/e + δUmomd/δp`
2. 传播位置：`r(t+dt) = r(t) + v·dt`
3. 计算力（所有势能项的梯度）：
   - Skyrme 局域 + 密度依赖
   - 库仑（查表法 `cy`/`cyu`）
   - 对称能（nsym=1/2/3）
   - 表面梯度
   - 有效质量
   - 动量相关（momp=0/1/2）
   - Pauli 势
4. 更新动量：`p(t+dt) = p(t) + F·dt`

### 4.6 `densitn`（行 2247-2270）— 坐标空间密度

```fortran
do i = 1, np
   do j = 1, np
      width = sl(i) + sl(j)                        ! 两核子波包宽度之和
      dd = exp(-dr² / width) / (π·width)^(3/2)      ! 高斯波包重叠
      denstl(i) += dd                               ! 总密度
      dens1(i)  += dd              (仅同种同位旋)    ! 质子/中子密度
   end do
end do
```

每个核子有独立的高斯宽度 `sl(i) = (c1·A^(1/3) + c2)²`（参见 `width0`，行 5157-5173）。

### 4.7 `init2`（行 1313-1556）— 从文件读取初始配置（碰撞模式）

当 `initi ≠ 0` 时使用，从 `aini(nru)` 和 `bini(nru)` 文件读取弹核和靶核，并执行：
- 质心归零
- Euler 角随机旋转
- 速度叠加（相对论 boost）
- 时间步长减半（预测-校正格式初始化）

## 5. 稳定性调优指南

### 5.1 稳定性判据链

```
init1 生成组态 → avene 计算结合能 → |Eb - progEb| > 0.2 → 调节 pfc → 重试
                                    ↓
                              结合能匹配 ✓
                                    ↓
             时间演化 (updtim, 0 → nt*dt fm/c)
                                    ↓
             每步检查: |rms(t) - rms(0)| > rmax ?
                ↓ YES                    ↓ NO (全程)
           lswitch=0                   lswitch=1
           不稳定，中止               k=nt 时输出
           (不写入 tini.dat)          "stable !!"
                                    写入 tini.dat
```

### 5.2 为什么初试核不稳定（常见根因）

| 症状 | 可能原因 | 修复方向 |
|------|----------|----------|
| RMS 迅速膨胀 (>rmax) | Skyrme 参数与动力学不一致 | 统一 LQMDin.in 与 LQMD.in 的 Skyrme 参数 |
| 结合能无法收敛 | `progEb` 不在当前 EOS 能给出的范围内 | 放宽 `progEb`（如 -7.615 → -7.0） |
| 大部分 run 飘移 | `rmax` 过严 | 增大 `rmax`（如 0.6 → 1.0） |
| 演化前期稳定，后期飘移 | `nt` 不够长，未检测到慢不稳定性 | 增大 `nt`（如 100 → 200） |
| 核子初始间距不合理 | `DISTNN`/`DISTNP`/`PHMANN`/`PHMANP` 不匹配 | 调整 `init1` 中硬编码的阈值（需改源码） |
| 密度分布不符合 Woods-Saxon | 均匀球采样不够精确 | 增大拒绝采样迭代次数，或启用 `RADRI`（Woods-Saxon 径向采样） |

### 5.3 推荐的调优步骤

1. **统一哈密顿量**：确保 `LQMDin.in` 第 3-6 行的 Skyrme 参数与 `LQMD.in` 第 5-7 行**完全一致**
2. **检查结合能匹配**：运行后查看 `0result.inf` 和 `ave.dat`，确认结合能收敛到 `progEb ± 0.2 MeV/u`
3. **监控 RMS 演化**：查看 `rms.dat`，RMS 应在约 10-20 fm/c 后稳定在小振幅振荡（<0.3 fm）
4. **检查能量分解**：启用 `laesw=1` 查看 `007E.out`，确认各能量分量物理合理
5. **统计成功率**：800+ runs 中应有 >80% 标记为 "stable !!"

## 6. 输出文件清单

| 文件 | 单元 | 说明 |
|------|------|------|
| `0result.inf` | 4 | 每个 run 的稳定性结果 + CM 误差 |
| `tini.dat` | 50 | 稳定 run 的最终相空间坐标（替换 coll.out 用作动力学初始条件） |
| `rms.dat` | 60 | 每时间步实时写入：`rtime, rms`（无数组缓存） |
| `ave.dat` | 70 | 每时间步实时写入：`rtime, aene`（无数组缓存） |
| `dtini.dat` | 80 | 每时间步每核子一行：`rtime, x, y, z, iso` |
| `007E.out` | 3 | 一行 16 列数值：`time, tkine, ucou, couex, ...`（每 100 步，laesw=1） |
| `Uint-R.dat` | 400 | 静态熔合势垒 vs. 间距（仅 initi=1） |
| `rmsnp.dat` | 55 | 中子/质子 RMS 半径预计算值 |

## 7. 关键 COMMON 块速查

| COMMON 块 | 核心变量 |
|-----------|----------|
| `(未命名)` | `x,y,z,px,py,pz(500)` — 核子相空间；`iso(500)` — 同位旋 (1=p,0=n)；`id(500)` — 粒子类型 |
| `/dd/` | `dt, ntstep, lswitch, rmax, radj` — 时间步与稳定性控制 |
| `/pot/` | Skyrme 参数、库仑查询表 `ih(100000)` |
| `/skyrme/` | `al, sal, width` — 波包宽度与 Skyrme 固定量 |
| `/grad/` | `gsur, rho0, cdelt, cmeps, cmnn, cmnp, momp` — 表面与动量相关参数 |
| `/bind/` | `progEb, aene, pfc, laesw` — 结合能调谐 |
| `/veloc/` | `velx, vely, velz(500)` — 速度分量 |
| `/INICLW/` | `DISTNN, DISTNP, PHMANN, PHMANP` — 初始化距离/相空间阈值 |
| `/frictn/` | `gamma_fric, tcool, gamma_fric_R` — 摩擦冷却参数 |
| `/rp0/` | `xi, yi, zi, pxi, pyi, pzi, emi(5000)` — 输出相空间数组 |
| `/gg/` | `nzta, nzpr, nsym` — 质子数与对称能类型 |
| `/wavel/` | `sl(500)` — 每个核子的高斯波包宽度 |
| `/deform/` | `beta2` — 四极形变参数 |
| `/shell/` | `ish, r0, ds, Ri0` — 壳修正与核半径参数 |

## 8. 版本与编译信息

- **文件**: `LQMDini_v4_0.for`，新增位置阻尼 + 输出格式统一
- **当前版本**: v4_0（详见 `LQMDini_changelog.md`）
- **编译器**: Intel Fortran `ifort`/`ifx`
- **关联输入**: `LQMDin.in`（与程序放在同一目录）
- **产物**: `tini.dat` → 可重命名为 `coll.out` 用于动力学主程序
