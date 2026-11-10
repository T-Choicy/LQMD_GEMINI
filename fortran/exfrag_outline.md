# exfrag 程序框架大纲

> 此文件为历史概览，保留原文供对照。当前版本请以 [exfrag 工程代码说明书](docs/exfrag_%E5%B7%A5%E7%A8%8B%E4%BB%A3%E7%A0%81%E8%AF%B4%E6%98%8E%E4%B9%A6.md)及其参数、接口、验证附录为准。下文中的输入名、全程退激、参考系和开关作用等概括未逐项反映当前实现，例如实际输入是 `coll.dat`，COLD 也包含未调用 GEMINI 的对象。

## 物理链总览

```
LQMD输运处理级联输出 (coll.out)
       │
       ▼
    MST模型  ── 相空间聚类，识别热碎片 (Z, A, E*, J)
       │
       ▼
  GEMINI 统计衰变  ── 热碎片逐级退激 (蒸发/裂变/IMF/γ)，生成冷碎片
       │
       ▼
  截面统计输出  ── d²σ/dEdΩ, 同位素截面, 动能分布, 激发能演化 …
```

即: 程序读入 LQMD 输运模型产生的核子相空间，用MST法识别末态热碎片，对每个热碎片调用 GEMINI 做统计衰变，最后在实验室系统计冷碎片的各类截面和分布。

---

## 文件结构

```
exfrag.f                         ← 主驱动 (2 行)
  ├── INCLUDE 'GEMMOD.f'         ← 17 个 MODULE，GEMINI 基础设施 (10784 行)
  └── INCLUDE 'anahyp_v3.0.f'    ← 主程序体 MAIN PROGRAM (4050 行)
```

| 文件 | 行数 | 角色 |
|------|------|------|
| `exfrag.f` | 2 | 驱动文件，INCLUDE 其余两个 |
| `GEMMOD.f` | 10784 | 模块库：类型定义、物理函数、GEMINI 衰变引擎 |
| `anahyp.f` | 4050 | 主程序：输入→MST→GEMINI→统计输出的完整工作流 |

---

## 一、GEMMOD.f — 模块速查表

GEMMOD.f 包含 17 个 MODULE，按功能分层：

| # | 模块 | 行号 | 功能 | 核心入口 |
|---|------|------|------|----------|
| M1 | `gem_kind` | 4–50 | 全局类型与常量 | `info`(控制参数), `identity`(碎片身份), `vell`/`spinn` |
| M2 | `trig` | 54–226 | 角度制三角函数(单/双精度重载) | `SIND, COSD, TAND, ASIND, ACOSD, ATAND, ATAN2D` |
| M3 | `sierk_stuff` | 231–1082 | Sierk 裂变势垒 | `bfis_sierk, yrast_sierk, fisrot` |
| M4 | `sortsub` | 1086–1519 | GEMINI 事件文件读写 | `open_event_file, get_event, start_sort, name` |
| M5 | `spectra` | 1523–1672 | 能谱直方图 | `hist_ek, hist_write, hist_normalize` |
| M6 | `light_stuff` | 1676–2150 | 轻粒子多重数输出 | `write_mult` |
| M7 | `mass_stuff` | 2154–2577 | 质量表与结合能 | `getmass, prepare_mass, prepare_mass_tf` |
| M8 | `functs` | 2581–3046 | 通用数学函数 | `gauss, exp_decay, ran_j, get_little_a` |
| M9 | `tl_stuff` | 3052–3975 | 穿透系数(光学模型) | `get_tl, read_tl, prepare_tl, prepare_deformed_emission` |
| M10 | `stat_mod` | 3979–4110 | 能级密度 | `get_level_density, collective_enhancement, alloc_storage` |
| M11 | `yrast_mod` | 4116–4756 | Yrast 线与鞍点 | `yrast, yrastld, saddle, sadfit, read_sad` |
| M12 | `angle` | 4762–5953 | 粒子发射角度生成 | `sphere_lp, sphere_imf, quantum_lp, fast_fission_angle` |
| M13 | `evap` | 5959–6475 | 简单蒸发模型 | `evaporation_simple, light_particle, prepare_evap` |
| M14 | `evaporation` | 6481–7262 | Hauser-Feshbach 蒸发 | `Hauser_Feshbach, prepare_Hauser_Feshbach` |
| M15 | `nonstatistical_mod` | 7268–7973 | 非统计发射 | `non_statistical, prepare_nonstatistical` |
| M16 | `statistical_mod` | 7978–10129 | **GEMINI 衰变核心** | `decay, statistical, decay_widths, gamma_ray, IMF, fission_width` |
| M17 | `gemsub` | 10134–10784 | **顶层接口** | `event, input, fusion, dmp, INPUT_ANA` |

**模块调用层级**: `gemsub`(顶层) → `statistical_mod`(衰变核心) → `evaporation`/`evap`(蒸发) + `angle`(角度) + `tl_stuff`(穿透) + `yrast_mod`(鞍点) + `functs`/`stat_mod`(能级密度) + `sierk_stuff`(裂变势垒)

---

## 二、anaHyp_v3.0.f — 主程序流程

### 总流程

```
[区块0] 声明与常量定义
    ↓
[区块1] 读取输入文件 + 初始化数组 + 打开文件
    ↓
[区块2] 参考系运动学 (Lab→CM→NN-CM)
    ↓
[区块3] 事件循环 (三层: 事件→时间步→核子)
    │      ├── 读 QMD 核子相空间
    │      ├── COLM=1: 并合模型 → 热碎片
    │      └── GEMINI=1: 统计衰变 → 冷碎片
    ↓
[区块8] 截面归一化 + 最终输出
```

### 区块 0: 声明 (第 44–262 行)

- **物理常量**: `RMASS=0.938918` GeV (核子质量), `PMASS=0.138039` GeV (π介子质量), `rho0=0.16` fm⁻³ (饱和密度)
- **Bin 参数**: `dkin=1` MeV (动能bin), `dthe=0.1`° (角bin), `MG=1000` (最大核子数)
- **核心 COMMON 块**: `XYZ`(位置), `MOM`(动量+同位旋+电荷), 势能分量(`fr, tkin, utwo, uthr, ucou, umod, usym, ussr, ugrd, couex, utau, aene, excen, bdgd`)
- **碎片数组**: `IC`(核子→簇归属), `XJR(簇,1:3)`(质子/中子/超子数), `ICL`(簇核子数)
- **谱分布**: `RAPIPI, PTRSPI, PTPION`(π介子), `drflow, epflow`(流), 碎片谱(`RPPRLG, PTRSCL, XKINSPC`), 激发能/速度/同位素分布(`ENERGY, velocity, Aisot`)

### 区块 1: 初始化 (第 263–900 行)

读取 `anahyp.in` → 设置控制开关和物理参数 → 读取 `massfrdm.dat` 结合能表 → 打开输入/输出文件 → 数组清零 → 读 QMD 数据头 (事件数、碰撞参数、入射能、波包宽度)

### 区块 2: 参考系运动学 (第 901–1034 行)

建立 Lab→CM→N-N CM 三套参考系, 计算 Mandelstam-s、弹靶快度、几何截面 `csf = 10×2π×x00×dx0` (mb)

### 区块 3: 事件循环 — 核心 (第 1035–3101 行)

```
DO L1 = 1, npar                     ← 事件循环
  ├── 读取 SURF_COAL 表面碎片 (可选)
  ├── DO it = 1, kktot             ← 时间步循环
  │     ├── DO i = 1, mass         ← 核子循环: 读入 (x,y,z,px,py,pz,em,lb,isy,idr,idat)
  │     │     ├── 累积质心坐标, 统计共振态/介子/超子计数
  │     │     └── 计算中心密度 rhoc
  │     ├── 质心归一化 + 介子/核子分析 (FSVT π谱填充)
  │     ├── [COLM] 并合模型: MST 聚类 → 热碎片 (Z,A)
  │     └── [GEMINI] 统计衰变:
  │           ├── 计算每个碎片的激发能 E* (势能求和)
  │           ├── 筛选合格碎片 (A>0, Z>MCFR, E*∈[DECNEXC1,DECNEXC2])
  │           ├── 初始化 GEMINI: input()→prepare_mass→prepare_HF→alloc_storage
  │           ├── fusion(L,cn) + event(cn,weight)  ← GEMINI 核心
  │           └── 衰变产物回传 + 洛伦兹变换到母核系
  ├── CM→Lab 洛伦兹变换 (所有碎片)
  └── 冷碎片在线统计: CNPCS, IPCS, DDCSS, ELEKED, TEELECS
```

**并合模型(MST)** (第 1463–1729 行): 两个核子属于同一碎片的判据为 **Δr² ≤ rr0² 且 Δp² ≤ pp0²**, 辅以簇间距检查确保空间局域性。

**激发能计算** (第 2075–2512 行): 基于 QMD 高斯波包, 对簇内所有核子对求和: Skyrme 二体/三体势 + 库仑 + 对称能 + 动量相关(Pauli)势 + 梯度势 + 有效质量项, 总能量减基态结合能得 E*。

**GEMINI 衰变** (第 2075–2830 行): 对每个合格热碎片调用 `gemsub:event()`, 级联衰变至所有子核 E* 低于发射阈值。支持 n/p/d/t/³He/α 蒸发、裂变、IMF 发射、γ 退激。Hauser-Feshbach 模式需光学势穿透系数, 简单蒸发模式快速近似。

### 区块 8: 最终输出 (第 3119–3598 行)

| 输出量 | 公式 | 单位 |
|--------|------|------|
| DDCSS 双微分截面 | `d²σ/dEdΩ = count × π·X00CUT² / (2π·sinθ·dkin·dthe·N_eff)` | mb/(sr·MeV) |
| LNKED 动能分布 | `dσ/dE = count × π·X00CUT² / (N_eff × dkin)` | mb/MeV |
| CNPCS 截面 | `σ = count × π·X00CUT² × 10 / N_eff` | mb |
| FSVT π⁻/π⁺比 | `ratio ± ratio × sqrt((σ_num/N_num)² + (σ_den/N_den)²)` | — |

冷碎片与热碎片分别输出 (CNPCS_HOT / CNPCS_COLD)。N_eff 为有效事件数。

---

## 三、输入输出文件速查

### 输入

| 文件 | 句柄 | 内容 |
|------|------|------|
| `anahyp.in` | 5 | 运行控制、模型开关、物理参数 |
| `massfrdm.dat` | 1412 | 结合能表 (IZ, IA, ELS, I1, Be) |
| `coll.out` | INFILE+2 | LQMD 级联输出 (核子+介子相空间) |
| `initangm.dat` | INFILE+1 | 初始角动量 (ANGM_DEC=1) |
| `SURF_COAL.dat` | INFILE+3 | 表面并合碎片 (SURF_COAL=1) |
| `DDCSS_LNKED_in.in` | — | 待分析核素+角度列表 |
| `IPCS_in.in` | — | 待分析同位素列表 |
| `ELEKED_in.in` | — | 待分析元素列表 |
| `TEELECS_in.in` | — | 时间演化元素列表 |

### 输出 (目录 `anaXXX/`)

| 文件 | 句柄 | 开关 | 内容 |
|------|------|------|------|
| `MONITOR` | OUTFILE+5 | 始终 | 有效事件数 |
| `DDCSS_LNKED_out` | OUTFILE+1 | DDCSS_LNKED | d²σ/dEdΩ + 轻核动能分布 |
| `ELEKED_out` | OUTFILE+2 | ELEKED | 元素动能分布 |
| `CNPCS_HOT` | OUTFILE+3 | CNPCS | 热碎片 Z/N/A 截面 |
| `CNPCS_COLD` | OUTFILE+4 | CNPCS | 冷碎片 Z/N/A 截面 |
| `TEELECS_out` | OUTFILE+6 | TEELECS | 元素截面时间演化 |
| `FSVT.dat` | — | FSVT | π⁻/π⁺ 谱比 |
| `IPCS_out.dat` | — | IPCS | 同位素产生截面 |
| `TEAMUL.dat` | — | TEAE | 最重碎片 E* 时间演化 |
| `EXCCROSS.dat` | — | TEAE | 激发能截面分布 |

---

## 四、控制开关速查

| 开关 | 含义 | 默认 |
|------|------|------|
| `COLM` | 并合模型 (Coalescence) | 1 |
| `GEMINI` | GEMINI 统计衰变 | 1 |
| `SURF_COAL` | 表面并合模型 | 0 |
| `ANGM_DEC` | 角动量关联衰变 | 0 |
| `DDCSS_LNKED` | 双微分截面 + 轻核动能 | 0 |
| `CNPCS` | 核素产生截面 | 0 |
| `IPCS` | 同位素产生截面 | 0 |
| `TEELECS` | 元素截面时间演化 | 0 |
| `ELEKED` | 元素动能分布 | 0 |
| `LNVICS` | 速度分布 | 0 |
| `FSVT` | π⁻/π⁺ 谱比 (对称能探针) | 0 |
| `TEAE` | 最重碎片 E* 时间演化 | 0 |
| `MST` | 最小生成树聚类模式 (=1) | 1 |
| `READMESON` | 读取介子相空间 | 1 |
| `OUTENE` | 输出激发能各分量 | 0 |
| `TESTEX` | 测试激发能计算 | 0 |
| `DIFF_T` | 扩散温度截断 | 0 |
| `PRINTDEC` | 打印衰变产物 | 0 |

---

## 五、关键子程序调用链

```
MAIN (anahyp_v3.0.f)
  ├── gemsub:INPUT_ANA()           ← 读取 GEMINI 控制参数
  ├── [每事件, 每时间步]
  │     ├── [COLM] MST 聚类         ← 并合模型碎片识别
  │     ├── [GEMINI] 激发能计算     ← QMD 波包势能求和
  │     └── [GEMINI] 对每个热碎片:
  │           ├── gemsub:input(IZ, IN, nexc)
  │           ├── mass_stuff:prepare_mass
  │           ├── evaporation:prepare_Hauser_Feshbach() 或 evap:prepare_evap()
  │           ├── tl_stuff:read_tl + prepare_deformed_emission
  │           ├── yrast_mod:read_sad (若 IMF 选项≥2)
  │           ├── angle:quantum_ang_generator (若角度生成)
  │           ├── stat_mod:alloc_storage
  │           ├── gemsub:fusion(L, cn)        ← 构造复合核
  │           ├── gemsub:event(cn, weight)    ← ** GEMINI 级联衰变核心 **
  │           │     └── statistical_mod:decay → evaporation / evap / angle / …
  │           └── yrast_mod:de_c              ← 释放内存
  └── 洛伦兹变换 + 截面统计输出
```

---

## 六、适用场景与输出范围

- **入射道**: 重离子碰撞 (HICs)、质子/反质子/π介子/K介子诱导反应
- **能量范围**: 中能 (~100 MeV/n) 至高能 (GeV/n)
- **物理模型链**: LQMD(级联) → Coalescence(碎片化) → GEMINI(退激)
- **输出产物**: n, p, d, t, ³He, α, IMF, 裂变碎片, γ, π, K, Λ/Σ/Ξ 超核
- **输出物理量**: 双微分截面 d²σ/dEdΩ, 同位素截面 σ(Z,N), 动能分布 dσ/dE, 激发能演化, π⁻/π⁺ 谱比
