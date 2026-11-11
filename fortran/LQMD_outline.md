# LQMD_v5_0 — 兰州量子分子动力学（LQMD）主程序文档

## 1. 程序概述

`LQMD_v5_0` 是 LQMD（Lanzhou Quantum Molecular Dynamics）输运模型的**动力学主程序**，约 19,000 行 Fortran 代码，由 5 个源文件组成。其核心功能是：

- 模拟重离子碰撞、质子/介子诱导散裂反应的全过程
- 通过高斯波包 + 平均场（Skyrme 型）+ 两体碰撞（含 Pauli 阻塞）传播核子相空间
- 支持表面并合模型（Surface Coalescence）模拟预平衡发射
- 支持共振态产生/衰变、奇异粒子产生、反核子反应等多道过程

## 2. 文件架构

```
LQMD_v5_0.for                  ← 驱动器: INCLUDE 4 个子模块
  main_v5_0.for  (2,006 行)    ← 主程序 PROGRAM LQMD: 输入/输出/主循环
  evolt_v5_0.for (9,745 行)    ← 演化+碰撞+力+表面并合 (117 个子程序)
  nncoll_v5_0.for(3,833 行)    ← 强子-强子碰撞动力学 (12 个子程序)
  cross_v5_0.for (3,014 行)    ← 截面数据库 (~40 个函数)

ctrlL_v5_0.for (597 行)        ← 独立批处理控制程序 (PROGRAM CONTROL)
```

## 3. main_v5_0.for — 主程序

### 3.1 输入文件

主输入文件为 `LQMD.in`（单元 2003），辅助文件包括 `massfrdm.dat`（结合能表）和 `arndt.in`（π-N 截面）。

### 3.2 输入参数（LQMD.in，26 行）

| 行 | 参数 | 含义 |
|----|------|------|
| 1 | `nrmin, nrmax, nth_run, nth_noh` | 事件范围、并行标识 |
| 2 | `npr, nzpr, nta, nzta, nanti` | 弹核 A/Z、靶核 A/Z、反核子数 |
| 3 | `x00min, x00max, x00stp, x00cut, nx00` | 碰撞参数范围与扫描 |
| 4 | `epp, nt, dt, ITMIN, ITMAX` | 束流能量(GeV/u)、时间步数、步长(fm/c)、输出窗口 |
| 5 | `alpha, beta, gamma` | Skyrme 两体/三体/密度幂次 |
| 6 | `gsur, gtau, eta, gsym, symm, rho0` | 表面/有效质量/对称能/饱和密度 |
| 7 | `asym, bsym, csym, gamsy` | 对称能密度依赖 (nsym=3) |
| 8 | `vpa, vpb, vpc, q01, p01, e01` | Pauli 势参数 |
| 9 | `r0, ds, c1, c2` | 半径参数、弥散度、波包宽度系数 |
| 10 | `rsur, rmin` | 附加表面半径、最小半径截断 |
| 11 | `betap, betat, betac` | 弹核/靶核/复合核四极形变 β₂ |
| 12 | `ish, ipot, iwrite` | 壳修正(0=液滴)、势类型、输出间隔 |
| 13 | `cdelt, cmeps, cmom, xsym, momp` | 动量相关势参数 |
| 14 | `iflag1~iflag7` | 主控标志（见 3.3 节） |
| 15 | `mdel, mdin, npop, nkaon` | 介质修正标志 |
| 16 | `nsym, xmu` | 对称能类型(1/2/3)、密度幂次 |
| 17 | `kfile, v01, initi, ibin, ird, crsc, ifile, gradi, nmint` | I/O 与初始化控制 |
| 18 | `iprin, ifdelt, iflagcro, iflagpau, ipoli` | 打印/Delta/截面/Pauli/极化 |
| 19 | `OUTBARYON, OUTMESON` | 重子/介子输出开关 |
| 20 | `GETX00CUT, GETX00CUT_IN, GETX00CUT_ST, G_ROT` | 自动碰撞参数扫描 |
| 22 | `SURF_COL, PRINT_CLUST, PRINT_POT, NDATA` | 表面并合开关 |
| 23 | `K_START_SURF, K_END_SURF` | 表面并合时间窗口 |
| 24 | `RR0, D_SURF, COL_PATH, CLUST_LEAST_DIST_COEF` | 表面并合几何参数 |
| 25 | `H0_SURF(1..5)` | p/n/d/t/³He 的相空间约束参数 |
| 26 | `H1_SURF, INCL4_6` | 势垒缩放、A=4/6 团簇开关 |

### 3.3 关键控制标志

| 标志 | 值 | 含义 |
|------|-----|------|
| `iflag1` | <0 | 跳过碰撞（纯平均场） |
| | ≥0 | 启用碰撞 |
| `iflag2` | 0 | 标准欧拉更新 (updt) |
| | 1 | Yukawa/Coulomb/Viola 更新 |
| `iflag3` | 3 | 隐式重子平均场 (updtim) + 介子场 (udtmsn) |
| `iflag4` | 1 | 相对论碰撞 (rcollbb+resdecay+rcollmeson+mesonabs) |
| | 2 | 非相对论碰撞 (nrcoll) |
| `iflag5` | 0 | Yukawa/Coulomb 势 |
| | 1 | Viola 势 + 动量相关 |
| `iflag6` | 1 | 氘核形成开启 (ideut=2) |
| `iflag7` | 1 | 壳演化 (updtsh) + 颈动力学 (neck) |
| `iflagcro` | 0 | 旧截面 (rcollbb) |
| | 1 | 新截面 (crossnew + rcollnew) |
| `iflagpau` | 1/2/3 | Pauli 阻塞方法 (R34/重叠/密度) |
| `SURF_COL` | 1 | 表面并合开启 |
| `momp` | 0 | 无动量相关 |
| | 1 | Aichelin 对数型动量相关 |
| | 2 | Yukawa 型动量相关 |
| `nsym` | 1 | 线性对称能 |
| | 2 | Csym*(ρ/ρ₀)^μ |
| | 3 | Skyrme 型 + MDI 对称能 |

### 3.4 主循环流程

```
PROGRAM LQMD
│
├─ 读取 LQMD.in、massfrdm.dat、arndt.in
├─ 初始化: 截面表、库仑表、Pauli 参数、波包宽度
│
├─ do 碰撞参数扫描 (x00min → x00max, step x00stp)
│   │
│   └─ do 事件循环 (nru = nrmin → nrmax)
│       │
│       ├─ 初始化: 随机种子、碰撞参数采样、init/init2
│       │
│       ├─ do ideut 循环 (1 或 2 次，iflag6=1 时开氘核)
│       │   │
│       │   └─ do 时间步循环 (k = 1 → nt)
│       │       │
│       │       ├─ [表面并合] CONSTR_AND_TUNNEL → EMISSION_TEST
│       │       ├─ [碰撞] rcollbb/resdecay/rcollmeson/mesonabs
│       │       ├─ [平均场] updtim + udtmsn (+ updtsh + neck)
│       │       └─ [输出] 重子/介子/表面并合数据
│       │
│       └─ 统计累加
│
└─ 最终输出、关闭文件
```

## 4. evolt_v5_0.for — 演化/碰撞/力/表面并合模块

这是最大的文件（9,745 行），包含 117 个子程序，是 LQMD 的物理核心。

### 4.1 子程序分类索引

#### 时间演化器（4 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `updt` | 3573-3631 | 最简欧拉传播：自由流 + Skyrme 力 |
| `updtyc` | 3633-3758 | 含 Yukawa/Coulomb + 对称能力 |
| `updtv` | 3760-3958 | 速度相关传播：含动量相关 + Pauli 势 |
| `updtim` | 3960-4510 | **主力演化器**：完整 Skyrme+Coulomb+对称+表面+有效质量+动量相关+Pauli |

#### 介子演化（1 个）
| `udtmsn` | 4512-4856 | 介子传播：Coulomb + in-medium 色散 + 洛伦兹力(kaon) |

#### 碰撞处理（8 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `rcollbb` | 56-324 | 重子-重子碰撞（含共振产生/Pauli 阻塞） |
| `rcollmeson` | 326-534 | 介子-介子 + 介子-重子碰撞 |
| `resdecay` | 536-730 | 共振态衰变 (Δ→Nπ, N*→Nπ/Nη) |
| `mesonabs` | 733-910 | 介子重吸收 (πN→Δ, ηN→N*) |
| `multcoll` | 913-1090 | 多体碰撞 (NNπ→NN 吸收) |
| `rcollv` | 1092-1374 | 速度相关 NN 碰撞 |
| `rcollnew` | 1376-1589 | 参数化截面 NN 碰撞 |
| `nrcoll` | 1820-1926 | 非相对论 NN 碰撞 |

#### Pauli 阻塞（6 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `paulin` | 1928-2055 | 原始 Pauli: 6D 相空间重叠积分 |
| `pauli12` | 2057-2080 | 简化 Pauli: Fn 占领概率 |
| `paur34` | 2082-2241 | 改进 Pauli (iflagpau=1): 双粒子重叠法 |
| `pauover` | 2243-2380 | 直接高斯重叠 Pauli (iflagpau=2) |
| `paurho` | 2384-2557 | 密度相关 Pauli (iflagpau=3): 局域费米球 |
| `constr` | 6584-6717 | 约束动力学 Pauli: 高占领粒子散射 |

#### 初始化（2 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `init` | 2585-2876 | 碰撞初始化：Woods-Saxon 采样 + Fermi 动量 |
| `init2` | 2878-3161 | 散裂初始化：从文件读取弹核/靶核 |

#### 密度计算（5 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `densit` | 3163-3185 | 总重子密度函数 |
| `densitn` | 3235-3262 | 每核子局域密度计算 |
| `densim` | 3264-3288 | 动量空间密度 |
| `dendgn` | 3291-3345 | 3D 格点密度分布 (50³) |
| `denscal` | 4860-4871 | 标量密度（相对论修正） |

#### 表面并合模块（14 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `CONSTR_AND_TUNNEL` | 8596-8692 | **递归团簇构建器**：从领头核子出发，逐层构建所有可能的团簇 |
| `CONSTR_LOOP` | 8466-8591 | 迭代团簇构建（非递归版本） |
| `PICK_NUCLEON` | 8698-8903 | **核心筛选**：对给定团簇，筛选符合并合条件的候选核子 |
| `EMISSION_TEST` | 8907-8955 | 发射能量检验：总能 + 结合能 + 隧穿概率 |
| `TUNNEL1` | 8958-9237 | **量子隧穿计算**：转移矩阵法计算 Coulomb+WS 势垒穿透 |
| `SORT_CANDID` | 9242-9304 | 候选团簇去重：通过核子 ID 的幂和比对 |
| `GET_CLUST_KIN` | 9306-9377 | 团簇运动学：动量、角动量计算 |
| `CLUSTER_STRONG_DECAY` | 9380-9487 | 不稳定团簇级联强衰变 (p/n/α 发射) |
| `STRONG_DECAY` | 9490-9728 | 衰变模式判定：Q 值 + 相空间 |
| `GET_TOTAL_ENERGY` | 8344-8437 | 单核子总能量（动能+势能各分量） |

#### 团簇识别（2 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `cluster` | 6987-7107 | 并合模型团簇识别: Δr<3.5fm 且 Δp<0.2GeV/c |
| `minspa` | 5135-5196 | 最小生成树: Δr<3.0fm |

#### 能量计算（2 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `avene` | 5985-6266 | 完整能量分解：动能+Skyrme+Coulomb+对称+表面+动量+Pauli+壳修正 |
| `total_energy` | 6378-6396 | 总能量求和 |

#### 诊断/分析（8 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `prot` | 3348-3468 | RMS半径/动量/共振态计数诊断 |
| `collnum` | 3470-3571 | 碰撞/阻塞统计累加 |
| `rdensity` | 5462-5559 | 径向密度分布 |
| `axyzdensity` | 5749-5847 | x/y/z 轴向密度投影 |
| `maxdensit` | 5850-5881 | 最大密度追踪 |
| `neck` | 6719-6903 | 颈区分析：密度极小点、粒子交换 |
| `center` | 6906-6964 | 弹核/靶核/体系质心计算 |
| `search` | 6269-6305 | 海森堡不确定性检验 |

#### 数学工具（10 个）
| 子程序 | 行号 | 说明 |
|--------|------|------|
| `coulgn` | 5884-5944 | 库仑势查表（双精度 erfc） |
| `coulg1` | 5043-5100 | 库仑势查表（旧版） |
| `erf` | 5102-5114 | 误差函数 |
| `erfc` | 5946-5982 | 双精度补误差函数 |
| `xexp` / `xdexp` | 6398-6436 | 安全指数函数 |
| `ranf` / `RANFL` | 5244-5298 | 随机数生成器 |
| `LORENTZ` | 7498-7514 | 洛伦兹变换 |
| `TWOPARTDECAY` | 7543-7568 | 两体衰变运动学 |

### 4.2 updtim — 主力时间演化器详解

```
1. 密度计算: call densitn → 每核子局域密度 dens1(i)
2. 动量相关速度: dvdpx/dvdpy/dvdpz = δUmomd/δp
   - momp=0: 无
   - momp=1: Aichelin 对数型
   - momp=2: Gauss 积分型
3. 速度赋值: vel = p/E + dvdp（含 Pauli 势修正）
4. 位置传播: r(t+dt) = r(t) + vel*dt
5. 力累加（双循环）:
   - Coulomb (fc): 查表法 cy(j)
   - Skyrme (f1): α·2ρ + β·2γ/(γ+1)·ρ^γ + gtau·2η·ρ^η
   - 密度依赖 (f3): -2βγ/(γ+1)·ρ^(γ-1) 非局域力
   - 有效质量 (ft): gtau·η·ρ^(η-1) 非局域力
   - 表面对称 (fsg): gsym·梯度(isospin 密度差)
   - 表面 (fgr): gsur·(grad ρ)²/ρ
   - 对称能 (fs): nsym=1/2/3 三种形式
   - 动量相关 (fv): cfm·log(p²·cmeps+1)² 或 Gauss 积分
   - Pauli 势 (fh): vpa·exp(-dr²/wq²) + vpb·exp(-p²/wp²) + vpc/Fermi
6. 动量更新: p(t+dt) = p(t) + dt*ΣF
```

## 5. nncoll_v5_0.for — 强子碰撞动力学

12 个子程序，处理所有强子-强子碰撞道。

### 5.1 子程序列表

| 子程序 | 行号 | 说明 |
|--------|------|------|
| `prodrsn` | 3-456 | 共振态产生: NN→NΔ, NN→NN*, NN→NΔK 等 |
| `collres` | 458-1100 | 含共振态碰撞: NΔ, NN*, NΛ, NΣ, ΔΔ |
| `prkbyby` | 1102-1592 | 重子-重子碰撞产生 Kaon: NN→NYK, NΔ→NYK |
| `prkpnby` | 1595-1963 | π/η/ω/ρ-重子碰撞产生 Kaon |
| `prkknhy` | 1965-2133 | Kaon-Hyperon 逆过程: KN→πN |
| `prbbkkb` | 2135-2737 | 重子-重子碰撞产生 KKbar |
| `prpnkkb` | 2739-3008 | π-重子碰撞产生 KKbar |
| `pnpnkkb` | 3010-3131 | ππ→KKbar 对产生 |
| `stranex` | 3133-3232 | 奇异交换: πY→KbarN |
| `elaskaon` | 3234-3291 | Kaon 弹性散射/电荷交换 |
| `absantikn` | 3293-3385 | 反 Kaon 吸收 |
| `proantinn` | 3388-3833 | 反重子-重子碰撞: 弹性/电荷交换/湮灭 |

### 5.2 碰撞判据

两体碰撞采用几何判据：在质心系中计算最近距离与碰撞参数 `b`，若 `b < √(σ_total/π)` 且轨迹在时间步内交叉 → 碰撞发生。终态动量按各向同性或角分布采样，终态相空间受 Pauli 阻塞检验。

### 5.3 Kaon 介质修正

- Kaon 在介质中获得排斥势 (~+28 MeV at ρ₀)
- 反 Kaon 获得吸引势 (~-100 MeV at ρ₀)
- 有效质量: `m* = √(m² - a·ρ + (bρ)²) ± bρ`，影响产生阈值

## 6. cross_v5_0.for — 截面数据库

约 40 个函数，提供所有反应道的参数化截面（单位 mb）。

### 6.1 截面分类

| 类别 | 代表函数 | 说明 |
|------|----------|------|
| NN 基本 | `CROSStot`, `CROSSels`, `CROSSins` | pp/np 总/弹性/非弹截面 |
| 共振产生 | `sigmand`, `sigmans`, `sigmadd` | NN→NΔ, NN→NN*, NN→ΔΔ |
| 共振吸收 | `sigddnn` | ΔΔ→NN 细致平衡 |
| Kaon 产生 | `signlk`, `signsk`, `sigdlk` 等 16 个 | NN/NΔ/πN/πΔ→NYK/ΔYK |
| KKbar 对 | `cnnkkb`, `cpbkkb` 等 6 个 | NN→NNKKbar, πN→NKKbar 等 |
| 奇异反应 | `spnlk`, `spnsk` 等 4 个 | πN→YK 参数化 |
| Kaon 散射 | `CKELS` | KN 弹性 |
| 反 Kaon | `CROKB` | KbarN 总/弹性/非弹/吸收 |
| Hyperon | `CHYPN` | YN 弹性+转换 (ΛN↔ΣN) |
| 反重子 | `SIGANTIBB` | 反核子-核子全道 (弹性/电荷交换/湮灭/超子产生) |
| 反 Delta | `NbarN_to_NbarDelta_dm` | pbar n ↔ nbar Δ⁻ (OPEM 模型) |

### 6.2 参考来源

- NN 截面: Cugnon 参数化 + PDG 高能外推
- Δ/N* 产生: Huber & Aichelin (Nucl. Phys. A 573, 1994)
- Kaon 产生: Tsushima, Sibirtsev, Thomas, Li (PRC 59, 1999)
- 反核子: Cugnon-Vandermeulen + PDG

## 7. ctrlL_v5_0.for — 批处理控制程序（独立可执行程序）

`PROGRAM CONTROL`（597 行），是独立的批处理调度器，**不被 LQMD_v5_0.for 包含**。

### 功能

- 读取 `INPUT.in`（任务数、每任务事件数、可执行文件名、超时）
- 读取 `LQMD0.in`（LQMD 物理参数模板，自动为每个子任务生成定制 LQMD.in）
- 为 `NFILE` 个子任务创建编号目录 (001, 002, ...)
- 编译 LQMD.for → 复制可执行文件 → 写入定制 LQMD.in → 后台提交 (nohup)
- 监控所有子任务的 `monitor*.dat` 文件，检测卡死/崩溃
- 崩溃恢复：kill 进程 → 从最后完成事件重新提交
- 全部完成后合并输出文件 (initangm.dat, SURF_COAL.dat, coll.out 等) 到 `results/`
- 可选调用后处理程序 `ctrlA.for`

### 子程序

| 子程序 | 说明 |
|--------|------|
| `Folder001` | 整数→3 位目录名 (1→"001") |
| `Folder0001` | 整数→4 位 nohup ID (1→"0001") |

## 8. 输出文件

| 文件 | 条件 | 内容 |
|------|------|------|
| `coll.out` | 始终 | 每时间步核子相空间 (x,y,z,px,py,pz,iso,isy,id) |
| `SURF_COAL.dat` | SURF_COL=1 | 表面并合发射团簇 (Z,A,Ek,px,py,pz,Lx,Ly,Lz) |
| `COUNTSUR.out` | SURF_COL=1 | 表面并合团簇计数统计 |
| `initangm.dat` | 始终 | 初始角动量分布 |
| `monitor*.dat` | ctrlL 批处理 | 各子任务进度 (当前事件号/碰撞参数) |
| `maxdensit.out` | 始终 | 最大密度时间演化 |

## 9. 关键 COMMON 块

| COMMON | 核心变量 |
|--------|----------|
|（无名）| `x,y,z,px,py,pz(500)` — 核子相空间; `iso(500)` — 同位旋; `id(500)` — 粒子类型; `idr(500)` — 弹/靶标签 |
| `/aa/` | `iflag1~iflag4, ncoll, nblock` — 控制标志与碰撞计数 |
| `/dd/` | `dt, ntstep, nru, lswitch, initi` — 时间步与运行控制 |
| `/pot/` | `alpha, beta, gamma, gsym, gtau, eta, Cex, asym, bsym, csym` — Skyrme 参数 |
| `/grad/` | `gsur, rho0, cdelt, cmeps, cmnn, cmnp, xsym, momp` — 表面/动量相关参数 |
| `/skyrme/` | `al, sal, width` — 高斯波包宽度 |
| `/coul/` | `cy(5000), cyu(5000)` — 库仑势查表 |
| `/wavel/` | `sl(500)` — 每核子波包宽度 |
| `/densi1/` | `denstl(500), dens1(500)` — 每核子局域密度 |
| `/medium/` | `mdel, mdin, nsym, npop, nkaon` — 介质修正标志 |
| `/SURF/` | `VANISH(500), H0_SURF(5), TREE, CLUST_STOR, TOT_CLUST_NUM` — 表面并合核心数据 |
| `/SURF1/` | `CLUST_LEAST_DIST, R0_SURF, D_SURF, COL_PATH` — 表面并合几何参数 |
| `/SURF_EXT/` | `H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)` — 结合能数据库 |
| `/TUNNEL/` | `Z_RESID, A_RESID` — 剩余核（用于隧穿计算） |
| `/rp0/` | `xi,yi,zi,pxi,pyi,pzi(5000)` — 初始相空间存档 |

## 10. 版本与编译信息

- **版本**: v5_0（含表面并合 + 结团时刻点输出）
- **文件**: `LQMD_v5_0.for` + 4 个 INCLUDE 文件，共 ~19,000 行
- **编译器**: Intel Fortran `ifort`/`ifx`
- **关联输入**: `LQMD.in`（主控参数）、`LQMD0.in`（批处理模板）、`INPUT.in`（ctrlL 批处理用）
- **产物**: `LQMD.exe`（动力学主程序）、`CONTROL.exe`（批处理控制器）
- **更新日志**: `LQMD_changelog.md`
