# exfrag 工程代码说明书

版本：2026-09-05；对象：[exfrag.for](../exfrag.for)、[anahyp.for](../anahyp.for)、[GEMMOD.for](../GEMMOD.for)。本文解释当前源码，不以历史注释代替实现，也不代表已经完成物理模型标定。

配套查阅：[参数手册](exfrag_%E5%8F%82%E6%95%B0%E6%89%8B%E5%86%8C.md)｜[全部过程与接口索引](exfrag_%E6%BA%90%E7%A0%81%E7%B4%A2%E5%BC%95.md)｜[验证记录与已知问题](exfrag_%E9%AA%8C%E8%AF%81%E4%B8%8E%E5%B7%B2%E7%9F%A5%E9%97%AE%E9%A2%98.md)。

**公式入门：** 新增[对照 Charity（2008）的公式精讲](exfrag_%E5%85%AC%E5%BC%8F%E7%B2%BE%E8%AE%B2_%E5%AF%B9%E7%85%A7Charity2008.md)，从能量账本、态密度和概率逐步推导，并逐项映射文献式（1）–（23）与当前实现。配合[交互教学计算器](exfrag_%E5%85%AC%E5%BC%8F%E4%BA%A4%E4%BA%92%E6%95%99%E5%AD%A6.html)可直接观察参数变化的结果；其中的简化计算不替代 GEMINI 运行。

## 1. 先理解这个程序负责哪一段计算

`exfrag` 是输运计算之后的碎片分析和统计退激程序。它读取已经产生的核子相空间记录，识别初级碎片，估算每个碎片的内部能量与激发能，再对符合条件的碎片调用 GEMINI，最后统计产物。它不自行演化核子运动方程，也不重新计算输运阶段的核子碰撞。

入口的全部源码是：

```fortran
      INCLUDE 'GEMMOD.for'
      INCLUDE 'anahyp.for'
```

因此，用户所说的“两个程序”在工程上分别是 **17 个 MODULE 组成的库**和**一个无名主程序及其 6 个辅助过程**。文本包含后形成一个编译单元、一个可执行程序。不要同时把 `exfrag.for`、`GEMMOD.for`、`anahyp.for` 当作三个独立源文件一起链接，否则会重复定义；也不要把另一个主程序 `LQMDini.for` 加进同一目标。

当前文件长度为：`exfrag.for` 2 行、`GEMMOD.for` 10786 行、`anahyp.for` 4009 行。附录列出 123 个过程，其中 GEMMOD 117 个、anahyp 6 个。所谓“模块”既可能指真正的 Fortran MODULE，也可能指 anahyp 中的逻辑区段；本文明确区分两者。

```text
LQMD 初始化与输运（本手册只核对数据接口）
          │ 核子/可选介子记录 coll.dat；可选 initangm.dat
          ▼
anahyp 无名主程序
  输入 → 时间片读取 → 相空间并合 → 热碎片与 E*、J
                                          │
                                          ▼
                gemsub::input → fusion → event
                                          │
                      statistical_mod::decay
                      ├─ 已知轻核能级衰变
                      ├─ Hauser–Feshbach/简化蒸发
                      ├─ IMF/裂变
                      └─ γ 通道与停止判据
                                          │
              frag(1:n_frag) 末端产物 ← 递归衰变
                                          ▼
anahyp 回填 → 参考系处理 → HOT/COLD、同位素、能谱、激发能统计
```

GEMINI 的顺序二体衰变背景可参见 Charity 收录于 IAEA 裂变/散裂模型讲习班的说明；这里的具体分支和修改仍以本地代码为准，不能直接套用其它 GEMINI/GEMINI++ 版本的默认参数或行为。[IAEA 原始讲习班资料](https://nds.iaea.org/records/zbw8e-5gd76)

### 1.1 适用边界

主路径的并合只接收 `idr<=0`、`idat<=0` 的核子，`XJR(:,3)` 在当前路径没有建立超子簇的有效流程。源码虽然保留 π/K/超子相关数组与旧注释，不能据此宣称当前 exfrag 已完整支持超核统计衰变。`idatpr` 被写死为 0，外层 `npr>0` 又使部分介子入射运动学分支不可达。

最可信的阅读与验证起点是普通核子入射、普通核子碎片、一个冻结时间片。多时间片、特殊粒子、表面并合、量子角分布和批次合并都需要单独验证。源码本身没有给出一个可据以宣称精度的统一束流能量范围。

### 1.2 “热”“冷”在这里的实际含义

HOT 是 GEMINI 调用之前的并合簇。COLD 是调用之后还保留在碎片列表中的核素，包括：GEMINI 末端产物、低电荷而直接按基态处理的轻簇，以及因激发能、电荷或动量门槛而未退激的原始簇。因此 **COLD 文件并不保证其中每个对象都经过统计退激**。关闭 GEMINI 时，HOT/COLD 的差别也不再代表物理退激。

GEMINI 遇到无可用粒子通道、或 γ 宽度占总宽度至少 99% 时，会把剩余 `Ex` 存入 `exx`，把 `Ex` 设为 0 后返回。不是逐个追踪所有低能 γ 到严格基态；主程序还会过滤 `Z<0` 的光子伪粒子。[停止逻辑](../GEMMOD.for#L8166)

## 2. 新使用者的阅读路线

| 目标 | 阅读位置 | 应当能回答的问题 |
|---|---|---|
| 跑对程序 | 本文 §3、参数手册 §1 | 工作目录是什么？究竟读哪个文件？ |
| 跟踪一个事件 | 本文 §4、§5 | 核子如何变成簇？数组什么时候改变含义？ |
| 理解激发能 | 本文 §6 | E* 是总量还是每核子量？基态质量来自哪里？ |
| 理解退激 | 本文 §7–§9 | 哪些核进入 GEMINI？通道怎么竞争？ |
| 使用输出 | 本文 §10 | 是截面、产额还是加权和？分母是什么？ |
| 做参数扫描 | 参数手册 §5–§7 | 哪些参数应与输运一致？哪些只是数值收敛参数？ |
| 修改或扩展 | 本文 §11、源码索引、问题附录 | 哪些隐式状态、分配和接口必须同步？ |

建议先用纸笔追踪一个只有几个核子的人工事件，再看一个重核事件的 `IC → XJR/ICL → EXCEN → cn → frag`。只看最终质量谱，很难发现输入错位或碎片身份不一致。

## 3. 构建、运行目录和输入协议

### 3.1 已核实的构建关系

当前 [Console1.vfproj](../Console1.vfproj) 的 Source Files 只有 `LQMDini.for`；已有 `Console1.exe` 不能当作 exfrag 使用。历史 BuildLog 也显示其链接的是 LQMDini 对象。

本次用本机 Intel Fortran 2025.3.2 对原始 exfrag 完成独立编译和链接。先打开已初始化 oneAPI 与 Visual Studio x64 环境的命令提示符，再在独立构建目录执行：

```bat
ifx /nologo /c /fixed /extend-source:132 /warn:all /check:all /traceback /Od /I"D:\Spallation\fortran" "D:\Spallation\fortran\exfrag.for"
ifx /nologo exfrag.obj /exe:exfrag.exe
```

第一行生成 `.obj/.mod`，第二行生成可执行文件。构建目录应可写，输出不应混入源目录。`/fixed` 指定固定格式，`/extend-source:132` 指定扩展有效列；不要把 `.for` 随意当作自由格式。`/check:all` 和 `/traceback` 适合诊断，不是物理正确性的证明。源码包含 Intel 扩展，例如整数逻辑表达式、实数下标及旧式 PAUSE，GNU 编译器不是已经验证的等价替代。Intel 官方文档明确区分编译与链接，并要求模块定义先于引用。[Intel 编译调用说明](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-1/invoke-the-compiler.html)

直接从普通 PowerShell 调用 ifx，编译可以成功而链接报 `could not find 'link'`。需要让环境初始化与编译链接处在同一进程环境中。仓库根目录的 `compile_fortran.ps1` 在子 `cmd` 中运行 setvars，环境不会自动回传到父 PowerShell；不能把打印“初始化成功”当作链接环境已准备好。

### 3.2 工作目录就是全部相对路径的基准

运行时不是相对于 exe 所在目录寻找输入，而是相对于当前工作目录。单次模式的目录应为：

```text
case/
├─ Input file/
│  ├─ anahyp.in                 # 主控制文件，含 GEMINI 参数尾部
│  ├─ massfrdm.dat              # QMD 侧的原子质量过剩表
│  ├─ IPCS_in.in                # 仅 IPCS=1 需要
│  ├─ DDCSS_LNKED_in.in         # 仅相应开关开启需要
│  ├─ ELEKED_in.in              # 当前示例与读取列数有冲突，见问题附录
│  └─ TEELECS_in.in             # 仅 TEELECS=1 需要
├─ coll.dat                    # NTH_RUN=0
├─ initangm.dat                # ANGM_DEC=1
├─ SURF_COAL.dat               # SURF_COAL=1
├─ tbl/                        # GEMINI 文本表
└─ gdat/                       # GEMINI 二进制表/缓存
```

在 PowerShell 中运行某案例：

```powershell
Set-Location -LiteralPath 'D:\Spallation\tmp\my-exfrag-case'
& 'D:\Spallation\tmp\my-exfrag-build\exfrag.exe' *> run.log
```

上面是操作方式示例，目录须由使用者为实际案例建立。**当前自带参数组合存在已复现的 unit 18 冲突，不能复制后直接视为生产基线。** 单独关闭冲突开关可以验证较窄路径，但不能证明其余输出都可用。

`NTH_RUN≠0` 时输入切换到 `results/coll.dat` 等；主要输出变成 `anaNNN/名称NNN_MMMM.dat`。`NNN` 来自 NTH_RUN 的三位字符，`MMMM` 来自 NTH_NOH 的四位字符。程序不创建 `anaNNN` 目录，也不真正恢复随机状态或中间统计量。

### 3.3 anahyp.in 的读取规则

它不是按键值解析的配置文件。`read(2003,*)` 逐行消费，每个节标题也消耗一个物理记录，后续数值位置固定。不能在参数段中插入空行、删掉标题，或按照注释自行添加一个参数。当前实际读 67 行：anahyp 读到第 45 行，`INPUT_ANA()` 接着从同一个 unit 2003 读第 46–67 行。

常见错觉有三种：

1. `Input file/INPUT_ANA.in`、`anahyp0.in`、`test.inp` 名称相近，但当前 exfrag 不读取它们。
2. 示例第 34 行注释写了五个名字，源码只读 `CNPCS, IPCS, TEELECS, TECOLD` 四个值，没有 `LCPs` 这一项。
3. `test.inp` 中的独立 GEMINI 初始核、`ratio` 等项目不属于当前输入协议；当前 `ratio_ana` 由调用处固定传入 `-1.`。

所有 67 行的逐项作用与样例值见参数手册。样例值是仓库当前配置，不等于推荐或自动默认值。

### 3.4 coll.dat 的精确格式

[读取位置](../anahyp.for#L920)与[粒子记录位置](../anahyp.for#L1225)：

```fortran
965  format(/1x,I4,5(F8.3,1x))
966  format(1x,4(I3,1x))
9669 format(1x,I4,7f9.4,5I4)
9670 format(1x,I4,7f9.4,6I4)
```

| 记录 | 依序字段 | 解释 |
|---|---|---|
| 头部第一物理行 | 跳过一行 | format 前导 `/` 的效果；通常为说明文字 |
| 头部数值行 | `nt,x00min,x00max,epp,dt,al` | 步数、冲击参数范围、束流动能 GeV/u、时间步 fm/c、宽度相关值 |
| 核组成行 | `nta,npr,nzta,nzpr` | 靶质量、弹质量、靶电荷、弹电荷；顺序不是弹在前 |
| 普通核子行 | `k,x,y,z,xk,yk,zk,em,LB,isy,id,idr` | 步号；fm；GeV/c；静质量 GeV；电荷、同位旋、身份、共振类别 |
| 可选介子行 | `k,xpion,ypion,zpion,xkpion,ykpion,zkpion,empion,idms,icp` | 紧跟每一条核子行，不是事件末尾独立的介子列表 |
| 表面并合核子行 | 普通行后追加 `VANISH,EMD_NUCL` | 已发射/删除标志、关联发射状态，使用 9670 |

普通路径虽使用带 `5I4` 的 format，实际只给 4 个整数接收项；不能仅凭 format 的长度推断实际变量数。数据记录是定宽而非 CSV，自由空格分隔不一定能读对。每事件固定 `mass=npr+nta` 条核子记录，每事件固定 `kktot` 个时间片；`READMESON=1` 时记录数再翻倍。主程序没有事件头自动寻找功能，也没有可靠 EOF 恢复。

`idatpr>0` 的额外身份分支当前不可达，而普通路径没有给 `idat` 赋值；`READMESON=0` 时一些介子数组仍被后续语句引用。这些依赖编译器/存储初值的代码必须作为缺陷处理，不能以“通常为零”代替显式初始化。

头部 `epp` 才是这里用于参考系变换的束流能量。分析配置文件中的 alpha/beta 等不会改变原来输运生成的束流。`X00CUT` 只参与面积归一化，没有在事件读取处真正剔除 `b>X00CUT` 的事件。

### 3.5 GEMINI 表文件与缓存

| 用途 | 优先路径 | 文本回退/注意 |
|---|---|---|
| 基础质量表 | `gdat/mass_tf.gdat` | `tbl/mass_tf.tbl`，随后 `prepare_mass` 还读 `tbl/mass.tbl` |
| 球形穿透系数 | `gdat/tl.gdat` | 源码回退写的是根目录 `tl.tbl`，而仓库文件在 `tbl/tl.tbl` |
| 变形穿透系数 | `gdat/tl_def.gdat` | `tbl/tl_def.tbl` |
| 极化穿透系数 | `gdat/tl_pol.gdat` | `tbl/tl_pol.tbl`；只在 polarization 开启时需要 |
| 鞍点插值表 | `gdat/sad.gdat` | `tbl/sad.tbl`；IMF option≥2 |
| 角分布 CDF | `gdat/q_ang.gdat` | 重建后却写到根目录 `q_ang.gdat`，读写位置不一致 |
| K 求和辅助表 | `gdat/plm2.gdat` | 可按 `prepare_plm2` 重建 |

`.gdat` 是顺序非格式化二进制记录，不是平台无关序列化；编译器的记录标记、字节序和数组尺寸都要一致。修改文本表后，程序可能继续使用旧缓存；不能只替换 `.tbl` 就声称新模型生效。应在独立案例副本中明确选择缓存版本或重建路径，保留哈希，避免并行案例共享可写缓存。

## 4. anahyp 的执行骨架与状态机

| 源码行 | 逻辑区段 | 产物/状态 |
|---|---|---|
| 28–268 | USE、COMMON、参数、数组 | 建立主程序共享状态 |
| 283–342 | 参数与 XBIN 质量过剩表 | 读取后无统一参数合法性检查 |
| 347–487 | 输入输出路径与 unit | 命名规则、开关、文件句柄 |
| 498–916 | 直方图清零和常量 | 非所有变量都清零 |
| 920–1069 | 输运头、运动学、并合阈值 | `mass,RCR,PCR,VPR1` |
| 1091–1335 | 事件/时间片/核子读取 | 原始核子与可选介子记录 |
| 1342–1589 | 选时间窗、密度、介子统计、总角动量 | 临时质心平移后又还原 |
| 1592–1993 | 并合、编号、热簇统计 | `IC,XJR,ICL,CL*,CLK*` |
| 1996–2047 | 总角动量减去外部轨道角动量 | `TOTANGM3` |
| 2090–2380 | 簇内相空间与激发能 | `AENE,BDGD,EXCEN,CLKIN` |
| 2381–2502 | 能量诊断、最大残余核统计 | `E_SH,COMPEXC,EXCCROSS*` |
| 2508–2858 | 热簇逐一调用 GEMINI | 追加新槽位、清除母簇槽位 |
| 2861–3103 | 表面产物回填、CM→Lab、冷簇统计 | 各类计数器 |
| 3109–3115 | 时间片/事件结束 | MONITOR、批次退出 |
| 3140–3673 | 归一化和输出 | 输出文件，主程序 END |
| 3680–4009 | ISOTAB/sl/xexp/scou/LORENTZ 等 | 外部辅助过程 |

等价的阅读用伪代码如下；它保留了当前不太直观的门控关系：

```text
读配置，打开文件，初始化
读 coll 头，mass = npr + nta
for L1 = 1..npar:
    若表面并合开启，先读本事件表面产物
    for it = 1..kktot:
        读 mass 个核子及可选介子；必要时读初始 J
        若为批次模式且 L1 不在区间，跳到事件尾（有多时间片错位风险）
        NNEVENT += 1                    # 是读取的时间片数，不是唯一事件数
        if kkmin <= it <= kkmax:
            临时扣除全系统均值，分析介子，再还原
            if COLM == 1:
                并合；累计热碎片；角动量；计算 E*
                TEAE/ECHECK 等统计
                if GEMINI == 1:
                    遍历原热碎片，过滤、准备 GEMINI、生成衰变事件
                    追加 daughter，删除原热碎片
                对现存碎片做参考系处理
                填 COLD/IPCS/DDCSS 等
    输出 MONITOR 或按 NRMAX 退出
按各输出分支自己的分母写文件
```

`GEMINI` 的 if 嵌套在 `COLM` 之内；`COLM=0,GEMINI=1` 不会凭空获得碎片执行退激。`kkmin/kkmax` 是文件中的时间片序号，不是 fm/c；时间由 `k*dt` 得到。GEMINI 在选中的每个时间片都可能运行，时间窗并不自动表示“只分析最后一个冻结态”。

## 5. 数据结构与并合算法

### 5.1 索引是最容易读错的部分

| 名称 | 索引域与含义 | 生命周期/单位 |
|---|---|---|
| `X/Y/Z(i)` | 核子索引 i | fm；有暂时减去全系统中心的阶段 |
| `XK/YK/ZK(i)` | 核子索引 i | GeV/c；有暂时减均值的阶段 |
| `em(i),LB(i),isy(i),idr(i)` | 核子属性 | GeV、电荷、同位旋编码、共振标志 |
| `IC(i)` | 核子 → 簇编号 | 0 表示未纳入；退激后只是重新分配标签 |
| `XJR(c,1:3)` | 簇 c 的 Z、N、超子占位计数 | 实数数组保存本应整数的数量 |
| `ICL(c)` | 簇 c 的质量数 A | 应与 XJR(c,1)+XJR(c,2) 一致 |
| `CLX/CLY/CLZ(c)` | 初期是核子坐标和 | 末尾除以 A 后成为中心；新 daughter 却提前写入中心，存在语义不一致 |
| `CLKX/CLKY/CLKZ(c)` | 初期是总动量 | 2907 后经 boost 再除 A，成为 Lab 每核子动量 |
| `XFR/YFR/ZFR(i)` | 核子在所属簇的相对坐标 | 通过算术均值平移得到 |
| `XKFR/YKFR/ZKFR(i)` | 核子相对所属簇的平均动量 | 不是严格相对论静止系 boost |
| `fr(i),frp(i)` | 核子局域密度相关量 | 前者无量纲，后者 fm⁻³ |
| `tkin/utwo/uthr/ucou/…(c)` | 簇的各项总能量 | GeV；不是每核子值 |
| `AENE(c)` | 簇的总内部能量 | GeV，可以为负 |
| `BDGD(c),EXCEN(c)` | 总结合能、总激发能 | MeV |
| `DIFFUSION(c)` | 退激前簇动能/A | MeV/u；不是扩散系数，daughter 未继承 |
| `IFSTAT(c),IFSURF(c)` | 统计/表面产物标签 | 决定 DDCSS 是否统计 |
| `para` | GEMINI 控制对象 | MODULE 全局状态，非每事件私有 |
| `cn,frag(:)` | GEMINI 母核/末端产物 | 结构体字段见 §8.1 |

Fortran 不区分大小写，`NEXC` 因 I–N 隐式整型规则是整数；`EXCEN` 是实数。主程序没有 `IMPLICIT NONE`，阅读时必须先查声明再凭名称判断类型。源码索引中 GEMINI 过程的 INTENT 只覆盖显式实参，不覆盖 MODULE 的读写副作用。

### 5.2 并合不是只有一个距离判据

[主路径](../anahyp.for#L1620)内部设置 `MST=1`，并不是从 anahyp.in 读取。对候选核子 i、j，先要求：

\[
|\mathbf r_i-\mathbf r_j|^2\le r_{r0}^{\,2},\qquad
|\mathbf p_i-\mathbf p_j|^2\le p_{p0}^{\,2}.
\]

`rr0` 单位 fm，`pp0` 单位 GeV/c。代码先平方生成 `RCR/PCR`，比较时不再开方。它在当前输入坐标/动量系中做差，没有逐对洛伦兹变换。

若 j 尚无簇，代码计算 i 所属簇的中心和 RMS 半径：

\[
\bar{\mathbf r}_C=\frac1{A_C}\sum_{k\in C}\mathbf r_k,\qquad
r_{\rm rms,C}=\sqrt{\frac1{A_C}\sum_{k\in C}|\mathbf r_k-\bar{\mathbf r}_C|^2}.
\]

当该簇已有至少两个核子，还要求新核子满足：

\[
|\mathbf r_j-\bar{\mathbf r}_C|\le r_{\rm rms,C}+rr0.
\]

若两核子已经分别有簇，则计算两簇并集的 RMS，并要求：

\[
r_{\rm rms,C_1\cup C_2}\le1.14\,A_{C_1\cup C_2}^{1/3}\quad\mathrm{fm}.
\]

通过后把较大编号改为较小编号。因此它是带额外紧致度约束的逐对合并算法，并不等同于只按相空间邻接关系求连通分量的标准 MST。改变核子记录顺序可能改变增量合并过程，应做顺序扰动检查。外层两重核子循环里还有全核子扫描，最坏代价可接近 O(A³)。

`ISOTAB(Z,N)` 有一长串允许核素区间，但在并合中的调用被注释掉；当前没有通过该函数进行稳定性筛选。`Iso/Iram`、`rhy/RSH` 等遗留声明也不能当作已启用的物理分支。

### 5.3 聚类之后应成立的不变量

\[
A_c=\#\{i:IC(i)=c\},\quad
Z_c=\sum_{i:IC(i)=c}LB_i,\quad
N_c=A_c-Z_c.
\]

总核子数应等于所有合法簇质量和；加上被排除的共振/特殊粒子及表面已发射成分后，才能与完整事件的守恒量比较。

当前源码先在 1792–1812 累加 `XJR/ICL/CL*/CLK*`，后在 1815–1836 只压缩 `IC` 编号，没有重排前述簇数组。这会破坏上述不变量，是参数扫描前必须处理的问题。示例：旧标签 `[1,1,1,3]` 变为 `[1,1,1,2]`，但质量表仍在槽位 1 与 3 中；后续只遍历到新 `ICMAX=2`，槽位 3 会被漏掉。不是增大 `MG` 或改变 `rr0` 就能可靠修复。

## 6. 激发能重建：逐项公式与代码映射

### 6.1 符号和宽度

先区分总激发能、热激发能和发射相对动能，可参见[公式精讲第3节](exfrag_%E5%85%AC%E5%BC%8F%E7%B2%BE%E8%AE%B2_%E5%AF%B9%E7%85%A7Charity2008.md#L58)。Charity文献从热核已形成之后开始，不直接推导本节的QMD能量重建公式。

本节所有能量和势项先以 GeV 计算，最后转 MeV。令 C 为一个簇，`q_i=LB(i)`，`s_i=isy(i)`。通常核子同位旋取 ±1，但应与输运文件编码核对；若输入改用 ±1/2，对称项会缩小四倍。

波包宽度来自 [sl](../anahyp.for#L3926)：

\[
L_i=\left(0.08 A_{\rm origin}^{1/3}+0.92\right)^2\;\mathrm{fm}^2,
\quad A_{\rm origin}=\begin{cases}npr&i\le npr\\nta&i>npr.\end{cases}
\]

它依赖核子原始顺序与弹靶质量，**不按当前碎片 A 重新计算**。输入头的 `al` 后来被局部赋值覆盖，不能通过改头部 al 控制这里的能量宽度。`TESTEX=1` 另用 `MASS1` 计算宽度，但其准备代码已注释，见问题附录。

定义簇内相对变量：

\[
\mathbf r'_i=\mathbf r_i-\frac1{A_C}\sum_{j\in C}\mathbf r_j,
\quad \mathbf p'_i=\mathbf p_i-\frac1{A_C}\sum_{j\in C}\mathbf p_j.
\]

能量循环中的宽度和重叠核为：

\[
a_{ij}=2(L_i+L_j),\quad
t_{ij}=\frac{\exp[-|\mathbf r'_i-\mathbf r'_j|^2/a_{ij}]}
 {\rho_0(\pi a_{ij})^{3/2}},\qquad\rho_0=0.16\;\mathrm{fm}^{-3}.
\]

\[
f_i=\sum_{j\in C,j\ne i}t_{ij},\qquad
f_i^{(p)}=\rho_0\sum_{j\in C,j\ne i}q_iq_jt_{ij}.
\]

`fr=f_i`、`frp=f_i^(p)`；只有同簇核子对贡献能量。前面的“高密度核子计数”使用 `L_i+L_j` 且包括自项，中心密度用 `2L_i`；不能将这些诊断密度与这里的 `fr` 直接混用。[能量重叠核](../anahyp.for#L2151)

### 6.2 内部动能与 Skyrme 项

\[
T_C=\sum_{i\in C}\frac{|\mathbf p'_i|^2}{2m_i},\quad
U_2=\frac{\alpha}{2}\sum_i f_i,\quad
U_3=\frac{\beta}{1+\gamma}\sum_i f_i^\gamma.
\]

变量映射为 `tkin,utwo,uthr`，`beta=batat`，`gamma=gmmat`。[对应源码](../anahyp.for#L2197)

注意这里内部动能是非相对论表达式，而平动动能另用平方根公式。alpha<0 通常提供吸引，beta>0 提供随密度增长的排斥；它们与 gamma 是一个成套参数化，不能把三者独立调整而仍声称保持原输运模型的 EOS。

### 6.3 表面、体对称、表面对称与有效质量项

记

\[
D_{ij}=\frac6{a_{ij}}-\frac{4r_{ij}'^{\,2}}{a_{ij}^2}.
\]

则逐对累加的最终表达式为：

\[
U_{\rm surf}=g_{\rm sur}\sum_{i<j\in C}D_{ij}t_{ij},
\]

\[
U_{\rm sym}=2\,symm\sum_{i<j\in C}s_is_jt_{ij},\qquad
U_{\rm ssym}=2\,g_{\rm sym}\sum_{i<j\in C}s_is_jD_{ij}t_{ij},
\]

\[
U_\tau=g_{\tau11}\sum_{i\in C} f_i^{\eta11}.
\]

对应 `ugrd,usym,ussr,utau`。代码对同一个簇槽位分别通过 IS1、IS2 加半项，合起来得到上式。`gsur/gsym` 因 D 带 fm⁻²，量纲为 GeV·fm²；`symm/gtau11` 为 GeV。不要把 `gsym` 与 `symm` 混为同一对称能系数。[逐对势项](../anahyp.for#L2172)

### 6.4 Coulomb 直接项和交换项

理想的平滑直接势结构为：

\[
U_{\rm C}=\sum_{i<j\in C}q_iq_j\,
0.00144\frac{\operatorname{erf}(r_{ij}/\sqrt{a_{ij}})}{r_{ij}}
\quad\mathrm{GeV}.
\]

实际实现并不直接使用连续的 r：`idr22=int(40*r²+1.5)` 经 `ih` 查表，`nccc` 转换为约 0.01 fm 精度的距离，`scou` 使用 `r=(nccc-1)/100`。若 `nccc==1` 或 `5000`，它直接返回 0；因此 r→0 的正确有限极限并未按解析式实现。`idr2` 上限截到 100000，超出距离网格不能认为继续遵循 1/r。

交换项是：

\[
U_{\rm C,ex}=c_{\rm ex}\sum_i[f_i^{(p)}]^{1/3},\quad
c_{\rm ex}=-\frac34(0.00144)\left(\frac3\pi\right)^{1/3}.
\]

`erf` 用多项式近似：设 `t=1/(1+0.3275911*x)`，

\[
\operatorname{erf}(x)\approx1-e^{-x^2}
(0.254829592t-0.284496736t^2+1.421413741t^3
-1.453152027t^4+1.061405429t^5).
\]

此处调用 x≥0。`xexp` 则把指数参数夹到 [-87,87]，是避免溢出的数值保护，并非物理截断。[scou/erf/xexp](../anahyp.for#L3941)

### 6.5 动量相关项 MDI

\[
U_{\rm MDI}=\sum_{i<j\in C}C_{ij}\,
\ln^2\!\left(1+500|\mathbf p'_i-\mathbf p'_j|^2\right)t_{ij}.
\]

\[
C_{ij}=\begin{cases}
C_0(1+x_s),&q_i=q_j,\\
C_0(1-x_s),&q_i\ne q_j,
\end{cases}\quad x_s=-0.65,
\quad C_0=\begin{cases}0,&momp=0\\0.00176,&momp\ne0.\end{cases}
\]

动量单位 GeV/c，因此 500 的量纲对应 (GeV/c)⁻²。当前 `momp` 只控制这部分强度是否置零，没有从输入选择多套 MDI 模型。该项的对数平方是动量依赖势，不能因为旧打印文字出现 Pauli 就称作实际 Pauli 阻塞算法。

### 6.6 总内部能量、结合能和激发能

对 `Z>MCFR` 的簇：

\[
E_{\rm int}=T+U_2+U_3+U_C+U_{C,ex}+U_{\rm surf}
+U_{\rm sym}+U_{\rm ssym}+U_{\rm MDI}+U_\tau.
\]

`massfrdm.dat` 第五列的 `Be` 实际作为**质量过剩**存入 `XBIN(Z,A)`，不是 B/A。可以用氘验证：`XBIN(1,2)=13.136`，

\[
B(Z,A)=7.289Z+8.071N-XBIN(Z,A),\qquad
E^*=\max(0,1000E_{\rm int}+B).
\]

氘给出 B≈2.224 MeV，说明这个解释与代码吻合。`AENE` 是 E_int，`BDGD` 是 B，`EXCEN` 是总 E*。若需每核子激发能，要显式除 A。

对 `Z<=MCFR` 的簇，直接设 `AENE=-B/1000`，所以 E*=0；即使该簇的内部相对动量很大，也不会按完整势能重建。负 E* 同样被截为零，原始负偏差不再保留。输出中零激发能比例过高应先查此门槛与基态一致性，不能直接解释为“产生了大量冷核”。

缺失质量表项的 `XBIN` 被初始化为 0，不报错即可生成错误 B。修改质量库前应分别核对 QMD 侧 `massfrdm.dat` 和 GEMINI 侧 `tbl/gdat`，两套库不自动保持一致。

完整 E_int 不含独立的 Yukawa、Pauli 或额外壳修正能项；`OUTENE/TESTEX/PRINTDEC` 中若仍打印 `uyuk11/upaul11/Esh11`，不能据此认为这些项已实际计算。[总和与截断](../anahyp.for#L2206)

### 6.7 平动动能与 ECHECK 的解释

热簇有效质量使用：

\[
M_C=Z m_p+N m_n+E_{\rm int},\quad
K_C=\sqrt{M_C^2+|\mathbf P_C|^2}-M_C.
\]

这里 `m_p=0.938271998`、`m_n=0.939565330` GeV，`CLKIN` 仍为 GeV。

ECHECK 的 before 是被筛选核子的相对论自由动能和；after 是各簇 `sqrt((M_free+AENE)^2+P²)-M_free` 的和，已经包含内部能量。它们不采用同一完整哈密顿量，所以 `DeltaE≠0` 不是独立的严格能量不守恒证明，`DeltaE≈0` 也不是整个输运加退激链守恒的验收。ECHECK 的代码仅检查并合与能量重建之前/之后，不检查 GEMINI 全级联终态。[诊断公式](../anahyp.for#L2381)

## 7. 角动量与 GEMINI 接口

### 7.1 初始角动量的来源

`ANGM_DEC=1` 时每事件首个时间片从 `initangm.dat` 读 `L12,Jx,Jy,Jz`。随后估计：

\[
\mathbf J_{\rm remain}=\mathbf J_{\rm initial}
-\frac1{\hbar c}\sum_C\mathbf R_C\times\mathbf P_C
-\frac1{\hbar c}\sum_m\mathbf r_m\times\mathbf p_m
-\sum_{\rm surf}\mathbf J_{\rm surf},
\qquad(\hbar c)^{-1}\approx5.067731.
\]

再把模平方对应到 J(J+1)：

\[
J_{\rm eff}=\frac{\sqrt{1+4|\mathbf J_{\rm remain}|^2}-1}{2}.
\]

`ANGM_LIGHT` 比较的是 `ICL`，即**质量数 A**，不是注释所写的电荷 Z。多于一个 `A>ANGM_LIGHT` 的簇时跳过本次计算；代码没有先清空 TOTANGM3。初始向量本身又被逐时间片原地减去轨道项，第二个时间片未恢复原初始值。`ANGM_DEC=0` 时，随后传给 `fusion` 的 L 也没有显式赋值。这些情况不能当作“自动设置 J=0”。

通过门槛的重簇使用 `L=int(TOTANGM3)`；小簇设 L=0。这个接口只传入自旋大小，`fusion` 另设固定自旋方向；并未传递主程序计算的完整 J 向量。[角动量计算](../anahyp.for#L1996)｜[L 赋值](../anahyp.for#L2634)

### 7.2 一个热簇真正进入退激的全部条件

依序检查：

1. `COLM=1` 且 `GEMINI=1`，当前 it 在分析窗口内。
2. `A=Z+N>0`，`Z>MCFR`。
3. `|P_C|² <= (xmom*ICL(C))²`，其中 `xmom=sqrt(epp*(epp+2*0.938))`。
4. `NEXC=int(EXCEN(C)) > DECNEXC1`。
5. 若 NEXC>DECNEXC2，将 NEXC 截为 DECNEXC2 的整型值。
6. `Z<=DECIZ`。

因此 `DECNEXC2` 是**裁剪上限**，不是把超限碎片拒绝，也不是每核子上限。E*=2.9 MeV 先变成整数 2，当门槛为 2.0 时不会进入；E*=3.0 才通过。每个通过筛选的热簇只调用一次 `event`，不是 `R_num=100` 次。[门槛与输入转换](../anahyp.for#L2519)

### 7.3 调用顺序与隐式状态

```fortran
call input(IZ, IN, NEXC, -1.)          ! IZ=Z，IN=N，不是 A
call prepare_mass
if (para%tl_iwbc) then
  call prepare_Hauser_Feshbach()
  call read_tl()
else
  call prepare_evap()
endif
if (para%I_angle .or. para%k_sum) call factorial()
if (para%ratio /= 1.0 .and. .not. para%k_sum) &
  call prepare_deformed_emission(para%ratio,1.)
if (para%Z_imf_min > 2) call prepare_5li_5he()
if (para%imf_option >= 2) call read_sad()
if (para%I_angle) call quantum_ang_generator()
call alloc_storage()
! L 必须由调用者明确赋值；当前主程序部分分支未保证这一点
call fusion(L,cn)
call event(cn,weight)
call de_c
! 读取 frag 后清理
call DE_ARRAY_tl
call DE_ARRAY_angle
```

这是阅读用的自由格式接口骨架，并非完整可独立编译文件。完整调用还需要先从 unit 2003 执行 `INPUT_ANA()`，定义 USE/变量并准备表文件。`DE_ARRAY_angle` 同时释放 `stat_mod::storage`，名称并不能体现全部副作用。不能把清理顺序缩减成“只释放角度数组”。

上面的原样调用还包含一项已知状态问题：`prepare_5li_5he` 对 SAVE 数组原地归一化并累加成 CDF，每个热簇再调用一次会把旧 CDF 当成密度再次累计。外部复用前应确认这种重复准备已被正式修正；手册没有替用户修改该实现。

## 8. GEMMOD 的 17 个模块：职责、公式和使用条件

### 8.1 gem_kind：共享协议，而不是物理计算器

`gem_kind` 定义 `r4/r8`、π、半径常数以及 `info/vell/spinn/identity`。`para` 是 info 的全局实例。它没有独立动力学公式；任务是保证所有模块对量的意义一致。

| identity 字段 | 意义 | 注意 |
|---|---|---|
| `Z,A` | 电荷、质量数 | γ 的负 Z 是内部标志，不是普通核素 |
| `Id` | 衰变分支身份 | 用移位/置位标记二体分支，非原 QMD 核子编号 |
| `level` | 已知离散能级编号 | >0 时转非统计已知能级衰变 |
| `Ex` | 激发能 MeV | 终止时常置 0，旧注释的“终态为 -1”不适用当前递归返回 |
| `exx` | 停止时保留的剩余激发能 | 回填主程序没有完整使用它 |
| `time` | GEMINI 内部时间 | 单位由 0.65824/Γ 确定为 10⁻²¹ s |
| `J,K,M` | 自旋、体轴投影、实验轴投影 | ℏ 单位，允许半整数 |
| `spin%theta/phi` | 自旋方向 | 度 |
| `vel%v/theta/phi` | 速度大小和方向 | cm/ns、度、度 |

重要常数为 `r0=1.16 fm`、`sep=2 fm`、`r0ld=1.2249 fm`、`l_max_quantum=50`。这几个半径参数属于 GEMINI 几何，不是并合 rr0。[类型定义](../GEMMOD.for#L4)

### 8.2 trig：角度制函数

提供 SIND/COSD/TAND 及反函数的 r4/r8 重载：`SIND(θ)=sin(πθ/180)`，逆函数把弧度乘 180/π。`ATAN2D(y,x)` 用两个分量保留象限。模块没有拟合参数；接口名里的 D 表示 degree，不能当作 double precision。主程序的 `ACOS` 与这些角度制函数同时存在，是角度单位错误的重要来源。[trig](../GEMMOD.for#L54)

### 8.3 sierk_stuff：有限程旋转核势垒拟合

主要入口是 `jmax_sierk`、`bfis_sierk`、`yrast_sierk`、`M_inertia_sierk`、`bs_sierk`。先调用 `jmax_sierk(Z,A)` 更新模块内 Z、A 与多项式缓存，再调用其它只带 J 的入口；后者不是无状态函数。

勒让德多项式满足：

\[
P_0(x)=1,\quad P_1(x)=x,\quad
P_n(x)=\frac{(2n-1)xP_{n-1}(x)-(n-1)P_{n-2}(x)}n.
\]

令 z=0.01Z、a=0.0025A，代码用常数系数矩阵计算：

\[
J_{\max}=\sum_{k=1}^{7}\sum_{i=1}^{5}emxcof_{ki}P_{k-1}(z)P_{i-1}(a),
\quad B_f(0)=\mathbf P_z^T\,elzcof\,\mathbf P_a.
\]

J<J20 段用 `B_f(J)=B_f(0)(1+qa J²+qb J³)`；较高 J 用 J/Jmax 的多项式，负势垒及 J>Jmax 置 0。Yrast 与鞍点表面积同样是 Z、A、J 的多项式展开；本模块前面大量 DATA 常数是拟合表，不是建议日常逐个调整的自由参数。

本地边界检查为 `19<=Z<=111`、`1.2Z+0.01Z²<=A<=5.8Z-0.024Z²`；另在非零 J 时对 Z>102 返回哨兵值。界外不能把 999 或 -1 当物理势垒/角动量。模型理论来源是有限程表面能和弥散 Coulomb 修正的宏观旋转核模型。[Sierk 原始论文](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.33.2039)｜[本地入口](../GEMMOD.for#L638)

### 8.4 sortsub：历史 GEMINI 事件文件接口

负责文件头、事件读取、元素名与排序入口；没有新的核物理公式。`get_event` 按 GEMINI 自己的事件文件格式重建事件，与 `coll.dat` 完全不同。`open_event_file/start_sort` 包含交互式文件选择流程。当前 anahyp 并未调用这些入口；不要寻找一个自动产生的 `.evt` 再交给 exfrag，当前接口通过 `frag` 内存数组返回产物。[sortsub](../GEMMOD.for#L1086)

### 8.5 spectra：库内的能谱直方图

`hist_ek(ek,de,mode,yield)` 将中心 ek、宽度 de 的小区间分摊到能谱格点；`hist_normalize(gamma_total,j)` 将当前谱除以总宽度后累加到自旋 j 对应的数组切片，没有在此额外乘自旋权重。`hist_write/read/print` 处理 `.spe` 数据。期望的分摊为 `H_b += yield × 重叠长度/de`；当前源码把 `r_min/r_max` 先除以 `deltaE`，边缘格又混用格号与能量宽度，不能保证任意 `deltaE` 下等价于标准重叠积分。同一格分支直接加 yield，且没有先检查格号上界。若启用此历史功能，需要补积分守恒与边界测试。

当前 `ek_loop` 的 `hist_ek`、`decay_widths` 的 `hist_zero/hist_normalize` 调用被注释，所以主程序的 DDCSS 不是从这个模块直接输出。主程序有自己的 `cross_count/LNKEB`。[spectra](../GEMMOD.for#L1523)

### 8.6 light_stuff：轻核通道与已知能级表

保存 `Z1,A1,S,DEF,width,limit,level_no,mode_gs,exotic` 等并行表，每一个 mode 不一定对应一个不同核素，也可能是同一核素的不同能级。范围 -2:51，负/零 mode 包含多中子占位；`prepare_Hauser_Feshbach` 实际从 mode=1 扫描，不能说全部声明通道都自动启用。

候选通道要求 `Z1(mode)<Z_imf_min` 且 `exotic(mode)<=exotic_index`。粒子自旋用显式耦合，或在关闭 `ssum` 时乘 `(2s+1)` 简并因子。`DEF` 是质量过剩/能级相关能量，参与分离能 `Ba=M_d+DEF-M_parent+ex1`。`write_mult` 写概率累计多重数；当前主程序虽更新 mult0，不调用最终输出入口。[通道表](../GEMMOD.for#L1676)

### 8.7 mass_stuff：质量与配对修正

`prepare_mass_tf/prepare_mass` 填表，`getmass` 返回实验质量过剩、液滴基准质量过剩及可选 beta2。表索引用 `ind=chz(Z)-cha1(Z)+A`；cha1/cha2 列出 Z=0..100 的质量范围，界外调用 `massfr`，而不是简单返回“实验质量不存在”。

`massfr` 的质量过剩由核子质量过剩、体积、有限程表面、Coulomb/交换、Wigner、Coulomb 弥散、质子有限尺寸、常数和配对基准项组成。令 I=(N−Z)/A，最主要的部分是：

\[
\Delta M=m_n^\Delta N+m_H^\Delta Z
-a_v(1-k_vI^2)A+a_s(1-k_sI^2)A^{2/3}F(x)
+\frac{0.6e^2}{r_0}\frac{Z^2-0.76361Z^{4/3}}{A^{1/3}}
+E_{wig}+E_{cd}+E_{cpf}+E_{a0}+E_{ca}+E_p,
\]

\[
x=0.68/(1.16A^{1/3}),\quad
F(x)=1-3x^2+(1+1/x)(2+3x+3x^2)e^{-2/x}.
\]

系数 `av=15.9937,kv=1.927,as=21.13,ks=2.3`。尾项如下，能量单位 MeV，长度单位 fm：

\[
E_{wig}=36|I|-1.433\times10^{-5}Z^{2.39},\quad
a_d=0.7071\times0.99,\quad y=a_d/(r_0A^{1/3}),
\]

\[
E_{cd}=-\frac{3Z^2e^2a_d^2}{(r_0A^{1/3})^3}
\left[1-1.875y+2.625y^3
-0.75e^{-2/y}(1+4.5y+7y^2+3.5y^3)\right],
\]

\[
k_F=\frac{(7.06858Z/A)^{0.33333}}{r_0},\quad q=(r_pk_F)^2,
\]

\[
E_{cpf}=-\frac{0.125r_p^2e^2}{r_0^3}
(3.0208-0.113541667q+0.0012624q^2)\frac{Z^2}{A},
\]

\[
E_{a0}=4.4,\quad E_{ca}=0.212\frac{Z^2-N^2}{A},\quad
E_p=12/\sqrt A-10/A.
\]

这里 `r0=1.16,rp=0.8,e²=1.4399764,mnΔ=8.071431,mHΔ=7.289034`；末尾 Ep 在 massfr 中无奇偶 IF，虽然注释写 odd-odd，实际作为统一基准项相加。其它路径另用 get_pairing 进行奇偶修正。不要给 `massfr` 的返回值再次加 A×938 MeV 后再当质量过剩使用。[massfr](../GEMMOD.for#L2376)

`get_pairing` 的配对基准为：偶偶=0，奇 Z 偶 N 为 `4.8/Z^(1/3)`，偶 Z 奇 N 为 `4.8/N^(1/3)`，奇奇为两者之和再减 `6.6/A^(2/3)`；若奇奇且 N=Z 还加 `30/A`。这不是最常见的 ±12/√A 教科书记号，不应擅自替换。[配对修正](../GEMMOD.for#L2539)

`droplet_mass` 是另一组液滴参数化，当前主要 getmass 回退使用 massfr。两者不要同时叠加。该备用过程的公式为：

\[
\delta=\frac{(N-Z)/A+(3c_1/16Q)Z/A^{2/3}}
 {1+(9J/4Q)A^{-1/3}},\quad
\epsilon=\frac{-2a_2/A^{1/3}+L\delta^2+c_1Z^2/A^{4/3}}K,
\]

\[
E_v=(-a_1+J\delta^2-K\epsilon^2/2)A,\quad
E_s=(a_2+9J^2\delta^2/4Q)A^{2/3},
\]

\[
E_c=c_1Z^2/A^{1/3}-c_2Z^2A^{1/3}-c_3Z^2/A
-c_4Z/2^{1/3}-c_5Z^2.
\]

`a1=15.96,a2=20.69,J=36.8,Q=17,K=240,L=100,c1=0.73219,c2=0.00016302,c3=1.28846,c4=0.55911,c5=0.00049274`；此处 J 是对称能系数，不是自旋。再加 Wigner 项：N=Z且N为奇数时30/A，其余30|N−Z|/A；最后加核子质量过剩 `8.07169N+7.28922Z`。该过程中的额外奇偶项和原子结合能扣除被注释，不能补入当前公式。[droplet_mass](../GEMMOD.for#L2480)

### 8.8 functs：概率、角动量代数与能级密度参数 a

关于“为什么态密度是指数”“a与温度怎样联系”和壳效应消退，参见[公式精讲第4–5节](exfrag_%E5%85%AC%E5%BC%8F%E7%B2%BE%E8%AE%B2_%E5%AF%B9%E7%85%A7Charity2008.md#L99)，其中单独说明a随U变化时温度导数公式的适用条件。

`gauss` 用 Box–Muller：`r=sqrt(-2 ln u)`、随机 φ，交替返回 `r cosφ`、`r sinφ`。`exp_decay(Γ)` 返回：

\[
\Delta t=-0.65824\ln(u+10^{-37})/\Gamma.
\]

Γ 单位 MeV，对应时间单位 10⁻²¹ s；配置注释的 `10E-21` 不应被误读为 10⁻²⁰ s。`cl_go` 通过缩放阶乘计算 Clebsch–Gordan 系数，接口整数传入的是角动量/投影的两倍；`factorial` 必须先准备。`transform` 旋转方向，`addvectors` 是速度的笛卡尔相加，属于非相对论运动学。

`get_little_a` 先原地修正传入的 U，再计算 a、熵：

\[
U_1=U+E_{pair},\quad
U_2=\begin{cases}
U_1+\delta W,&mass\_option=0,\\
U_1,&mass\_option=1,\\
U_1+\delta W[1-e^{-(U_1+E_{rot})/18.52}],&mass\_option=2,
\end{cases}
\]

\[
a=\frac A{aden\_0}\times aa\_scale\quad(aden\_type=0),\quad
S=2\sqrt{aU_2},\quad T\simeq\sqrt{U_2/a}.
\]

U≤0 时熵返回 0。`aa_scale` 是每次调用的实参，常为 1；`para%a_scale` 只在特定裂变路径传入。所有 22 种 aden_type 的代码公式列在参数手册，不能把 aden_0 在所有模式中都解释成 A/a。[get_little_a](../GEMMOD.for#L2873)

`ran_j` 在允许区间均匀试抽 x，然后按 `exp(2 sqrt(a(U−ac*x²/2))−rua)` 接受，用于热转动模态。`zbeta(A)=A/(1.98+0.0155 A^(2/3))` 是某些 a 模型的参考稳定线。

### 8.9 tl_stuff：穿透系数、变形平均与插值

先 `read_tl` 载入表，再对某发射 mode、残余核 `(Z2,A2)` 调用 `prepare_tl`；随后 `get_tl(l,ek,…)` 获取通道透射概率。它没有实时求解光学模型 Schrödinger 方程，使用预计算表的参数化。

对每个 l，表系数生成 `c1,c2,c3`，普通分支：

\[
q_l(E)=c_1E+c_2+c_3/E,\qquad T_l(E)=\frac1{1+e^{q_l(E)}}.
\]

`prepare_tl` 的 c2、c3 是 `(l+1)` 的二次多项式。`get_tll` 在非物理驻点附近改用 100 或驻点值，指数>65 返回 T=0。mode≤1 的能量先加 1，某些质子/α mode 还有随 A2−A_EAL 的能标修正。

变形/膨胀模式以不同半径表做二次插值，再按表面积或 `|Y_lm|²` 权重平均 T；`ratio` 为轴比，`expand` 为半径放大率。`look_up` 缓存不能跨不同几何/残余核不清空地复用。l 上限 50 与角分布表共同约束；改上限必须同步重建缓存和数组。[get_tl](../GEMMOD.for#L3124)｜[get_tll](../GEMMOD.for#L3756)

### 8.10 stat_mod：缩放能级密度与采样存储

模块保存所有候选通道的 `Za/Aa/gamma/time/U_therm`，以及候选轻粒子的 `dau` 和积分 CDF `storage`。真正用于宽度比的密度不是一个未缩放的绝对 ρ：

\[
T=\sqrt{U/a},\qquad \sigma^2=\mathcal I T/40.848,
\]

\[
\widetilde\rho(U,J)=\frac{(2J+1)e^{S-ee}}
 {pe\;a^{1/4}[1+U^{5/4}(\sigma^2)^{3/2}]}.
\]

先以 `ee=S_parent,pe=1` 计算母核密度，再设 `pe=2π×母核缩放密度`，以后用差熵避免指数过大。不要把 `get_level_density` 的结果单独拿去当绝对能级密度而忽略 `pe/ee`。

`collective_enhancement` 在第一条有效执行语句处就 RETURN，下方 beta2/费米淡出公式没有生效。另一个名叫 `collective` 的标量在对称裂变宽度中生效，二者不是同一个开关。[stat_mod](../GEMMOD.for#L3981)

### 8.11 yrast_mod：旋转基线、非对称鞍点与形变

`yrast` 给指定 Z、A、J、K、ratio 的最低旋转能；ratio≤0 时 Z<19 走旋转液滴近似，否则走 Sierk；高于拟合自旋区间时以转动惯量外延。ratio>0 时：

\[
E_{rot}=\frac{34.540}{A^{5/3}}
[B_{r\perp}(J(J+1)-K^2)+B_{r\parallel}K^2]+E_{def},
\]

\[
E_{def}=(B_s-1)17.9439(1-1.7826I^2)A^{2/3}
+(B_c-1)0.7053Z^2/A^{1/3}.
\]

`def_rot` 按旋转椭球长短轴比计算 Bs/Bc 及惯量比；球形时全为 1。主要惯量比为 `Br_perp=2/(ratio^(-2/3)+ratio^(4/3))`、`Br_para=ratio^(2/3)`。

`sadfit` 对裂变性 x、旋转量 y、质量不对称度 α 的表做插值；x 超过约 0.61743 被夹断。`saddle` 用插值修正、液滴质量和接触 Coulomb 能得到实际分裂通道的鞍点：

\[
E_{sad}=sad\_array(\min(A_1,A_2))+M_{ld,1}+M_{ld,2}
+\frac{1.44Z_1Z_2}{R_1+R_2+sep},
\]

\[
\mathcal I_{orb}=\frac{A_1A_2}{A_1+A_2}(R_1+R_2+sep)^2,\quad
\mathcal I_{eff}=\left[(\mathcal I_1+\mathcal I_2)^{-1}
-(\mathcal I_1+\mathcal I_2+\mathcal I_{orb})^{-1}\right]^{-1}.
\]

`sad_array` 已扣除了参考质量和 Coulomb 项，不能再把 E_sad 简化成 B_f+Q 后重复加质量。先 `read_sad` 分配系数，再 `sadfit`，再 `saddle`；最后 `de_c` 释放系数。[yrast](../GEMMOD.for#L4544)｜[saddle](../GEMMOD.for#L4662)

### 8.12 angle：方向、自旋耦合与能量分配

`sphere_lp` 用经典角动量三角形构造发射方向，例如：

\[
\cos\theta=\frac{J_0^2+J_2^2-(l+s)^2}{2J_0J_2},
\qquad v_{rel}=\sqrt{2E_{rel}/(931.5 A_{red})}\,c.
\]

两体反冲速度按 `v1=A2/(A1+A2)*vrel`、`v2=A1/(A1+A2)*vrel` 分配，并与母核速度相加；实际代码有自己的质量/速度常数。退激角度内部以度计算。

`sphere_imf` 为两碎片分配弯曲、摇摆、倾斜、扭转热转动分量。各分量以 `ran_j` 的熵权重抽样，总耗能 `ΔE=1/2 Σ a_mode J_mode²`；剩余热能分配满足概率与 `exp[2√(a1 U1)+2√(a2 U2)]` 成正比。由 U1/U2 加各自旋转能构造 daughter Ex。

`quantum_lp/quantum_imf` 用 CG 系数平方、球谐/旋转函数构造角概率；`quantum_ang_generator` 将其预积分为 CDF，`find_theta` 反查。当前 `para%quantum=.false.` 被 input 写死，但仍可能准备 q_ang 表；“调用 quantum_ang_generator”不等于启用了量子角分布。

库内方向计算并不弥补主程序只传 J 大小、未传方向的接口限制。对实验室角关联的定量应用必须另外验证母核方向和 boost。[angle](../GEMMOD.for#L4764)

### 8.13 evap：简化蒸发

先准备最多 15 个简化 mode；分离能 `Ba=M_d+DEF−M_parent`，带电粒子的 Coulomb 势垒 `Ec=1.44 Z1 Z2/separation_c`，另有 `l(l+1)*const` 离心势垒。

\[
U_2=E_0^*-Ba-E_c-E_{rot,2}-l(l+1)const.
\]

低 U2 时显式求和，宽度增量为 `ρ̃(U2,J2) T (2s+1)`；U2>40 MeV 用 `light_particle_quick` 的积分近似。选定通道后，在 `0<=ek<=min(4T,U2)` 内按

\[
p(ek)\propto e^{2\sqrt{a(U_2-ek)}-S(U_2)}
\frac{1+U_2^2}{1+(U_2-ek)^2}
\]

接受/拒绝采样。这里的 ek 是超出 Coulomb 和离心势垒的部分，回传总相对动能会再加回势垒。

`tl_iwbc=.false.` 开启此近似；即便 `.true.`，衰变过程中 `Z0<=10` 也会落入此分支。因此“全程 Hauser–Feshbach”不是当前实现的准确说法。[evap](../GEMMOD.for#L5961)

### 8.14 evaporation：Hauser–Feshbach 蒸发

积分、三角耦合条件和宽度的数值计算示例，见[公式精讲第6节](exfrag_%E5%85%AC%E5%BC%8F%E7%B2%BE%E8%AE%B2_%E5%AF%B9%E7%85%A7Charity2008.md#L230)。特别注意文献式（1）的J表示发射粒子的总角动量。

物理骨架是对残余自旋、发射轨道角动量和能量求和/积分：

\[
\Gamma_c=\frac1{2\pi\rho_p(E_p^*,J_p)}
\sum_{J_d,l,j}\int T_{cl}(\epsilon)\,
\rho_d(E_p^*-Ba-\epsilon,J_d)\,W_{\rm coupling}\,d\epsilon.
\]

本地通过前述 `pe/ee` 缩放密度实现分母；K 求和/极化分支使用更细的耦合权重，不能把一个额外 `(2l+1)` 普遍叠加到所有分支。具体流程：

```text
Hauser_Feshbach
  → 分离能、离散态宽度采样、残余惯量
  → s2_loop（残余自旋）
      → 可选 k2_loop（残余 K）
      → ek_loop（动能）
          → ell_loop（l 及耦合）→ get_tl
          → get_little_a → get_level_density
          → gam_jj = rho_scaled × de × sum_tl
  → 从 storage CDF 选残余自旋/能量
  → 从 sub_store CDF 选 l、l+s、K
  → 填 dau 和 gamma
```

一般 `de=1 MeV`；尾部停止比例 `lower_ek/lower_ell=0.001`、`lower_s2/lower_k2=0.01`，与主配置的 `threshold` 是两组不同控制。`n_light_normal=24000`、K 求和时 250000；`n_sub=800` 约束内部候选表。

共振能级的额外激发能可按截断 Breit–Wigner 型 `width/2*tan(179.998u−89.999)` 抽样；5He/5Li/8Be 还用特殊线形表。`exotic_index` 决定可选能级集合，不是一个可连续缩放的发射概率。

源码 6602 写的是 `r1=r0*a1(mode)*.33333333`，不是 `r0*a1(mode)**.33333333`。它影响起始惯量与能量扫描位置；应作为实现疑点核实，不能在公式推导中悄悄改成预期立方根。[HF 主入口](../GEMMOD.for#L6532)

### 8.15 nonstatistical_mod：已知轻核离散态衰变

按母核 `(Z,A,level)` 的 SELECT CASE，选择预置分支、线宽、粒子组合和可用能量。最多返回 3 个 daughter，例如多中子或不稳定轻核的逐级解体。衰变时间仍用 `exp_decay(Γ)`，角度和反冲由 angle 辅助。它不是高激发连续能级密度模型，作用是处理已知轻核态，避免对离散结构直接使用统计连续近似。

能量分配基本受 `E_parent + Q = ΣE_daughter + E_relative` 约束，但特定支路使用内置能谱表和近似质量常数。`prepare_nonstatistical` 通过 SAVE first 只准备一次，与重复调用 `prepare_5li_5he` 的行为不同。[non_statistical](../GEMMOD.for#L7431)

### 8.16 statistical_mod：通道竞争、IMF、裂变与 γ

宽度转概率/寿命、裂变势垒、Moretto鞍点、热能分配及E1/E2的教学推导，见[公式精讲第7–11节](exfrag_%E5%85%AC%E5%BC%8F%E7%B2%BE%E8%AE%B2_%E5%AF%B9%E7%85%A7Charity2008.md#L298)。其中明确区分本地Fortran与文献所述GEMINI++处方。

`decay` 根据 `level` 分配至非统计或统计过程。`decay_widths` 的核心状态是：

\[
U_0=E_0^*-E_{yrast},\quad
\Gamma_{tot}=\sum_c\Gamma_c.
\]

先蒸发，再 E1/E2 γ，再根据 imf_option 添加对称裂变或所有质量不对称度；最多 2000 个候选通道。`threshold>0` 会按当前通道相对宽度裁掉后续搜索列表里的小概率蒸发模式，有可能影响罕见产物。

**无延迟的基本通道竞争**是 `P(c)=Γc/Γtot`，先累计宽度再用均匀随机数选择。概率分支有时间延迟修正和可选偏置权重，当前 `N_weight=0`，主分析直方图也没有使用 GEMINI 的 weight。

**按时间竞争**对每个通道抽等待时间，取最早者。IMF 延迟依赖质量不对称度：

\[
\eta=(A_0-2A_1)/A_0,\quad
t_d(\eta)=t_{delay}\exp[-\eta^2/(2\,sig_{delay}^2)].
\]

sharp_delay 时先等待剩余延迟再加指数寿命；非 sharp 分支模拟宽度随时间爬升。当前 input 只有在 `imf_option!=1,t_delay!=0,sig_delay!=0` 时设置 time_wise=true。[时间竞争](../GEMMOD.for#L8207)

**IMF 模式 2**遍历 Z1 从 Z_imf_min 到 Z0/2，以质量表范围决定 A1，再调用 `complex_frag`。其鞍点热能为 `U_s=Ex−E_sad`，加上对应母核壳/配对修正后，用近似过渡态宽度：

\[
\Gamma_{IMF}\approx T_s\,
\frac{(2J_0+1)e^{S_s-ee}}
 {pe\;a_s^{1/4}[1+U_s^{5/4}(\sigma_s^2)^{3/2}]}.
\]

开启 Lestone 时对 K 额外求和，再除 `(2J+1)`；额外旋转能是 `20.424K²/I_eff`。注意这个路径 `get_little_a` 最后传入 1.0，未使用 `para%a_scale`，也未直接使用 b_scale。

**对称裂变模式 1**用 `E_s=E_rot+b_scale*B_Sierk`（仅当 b_scale>0 才缩放），`a_s` 乘 a_scale，按 Δε=0.5 MeV 积分鞍点密度。宽度还乘 `collective=hω/T0`，其中 `hω=36.34 sqrt((1−x)/A)`、`x=Ec/(2Es)` 上限夹为 1；Kramer≠0 时再乘 Kramer。因此 Kramer=0 是“不施加修正”，不是“关闭裂变”。

**IMF 模式 3**先生成非对称通道，再用对称裂变宽度重新归一化靠近对称峰的那部分；`a_scale/b_scale` 在此可以通过 fission_width 间接起作用。

对称裂变的 Viola 相对动能为：

\[
K_{Viola}=0.1189Z^2/A^{1/3}+7.3\;\mathrm{MeV},\quad
E_{scission}=K_{Viola}-(M_0-M_1-M_2).
\]

**γ 发射**仅显式构建 E1/E2 两种统计通道：

\[
d\Gamma_{E\lambda}\propto E_\gamma^{2\lambda+1}
\,A^{2\lambda/3}\,(2\pi w_\lambda strength_\lambda)
\,\widetilde\rho_d\,dE_\gamma,
\]

`w1=6.8e-8,w2=4.9e-14`，步长一般 1 MeV、U0≤15 时 0.2 MeV。GDR 分支被常量 `.FALSE.` 禁用。E1/E2 强度≤0 会被恢复为 0.1/50，不能用 0 禁用 γ。

**形变涨落** `fluctuate` 在 K_sum 时搜索形变与 K，权重包含 `sqrt((5/4)*B/er0*T) × ρ̃ × sqrt(Br_parallel) / dqdx × Δq`，K≠0 加倍；`Δq=0.2`、`|q|<=5`、候选上限 10000。它改写父核 K 和当前 ratio，是更昂贵且需要独立验证的分支。[衰变核心](../GEMMOD.for#L7980)

### 8.17 gemsub：公开事件入口和递归级联

`INPUT_ANA` 只读共享 unit 2003 的参数尾部。`input(Z,N,NEXC,ratio)` 把热簇与配置转换到 para，同时强制设置 R_num=100、N_weight=0、I_angle=true、quantum=false、dExcn=0 等。`ratio=-1` 转为 `para%expand=true,para%ratio=0`，衰变中球形半径膨胀率为 1.1，A>170 时为 1.15。

`fusion(L,cn)` 建立母核，速度为 0、时间为 0；奇质量且 ssum=true 时，在 `L+1/2` 与 `|L−1/2|` 之间随机选择。不要把输入整数 L 直接等同于最后 cn%J。

`event` 清空 n_frag、n_binary、weight=1，调用 `init_decay`，然后递归：

```text
decay_event(parent):
    decay(parent) → daughter(1:n_extra)
    对每一个 daughter:
        若 level>0 且 Ex>=0，继续递归
        或 A>4 且 Ex>0，继续递归
        否则追加到 frag
```

当前 `frag` 存的是递归末端产物，不是整个历史树。`max_frag=100`；插入前没有全面的上限保护，重系统高多重数不能只靠外层 MG=1000 保证安全。

`dmp` 是可另行调用的事件文件写出入口，当前主程序未调用；`crash` 写 `gem.crash`，STATUS=NEW 意味着旧文件存在也可能影响异常记录。[gemsub](../GEMMOD.for#L10137)

## 9. 退激回填和参考系：怎样追踪同一个量

### 9.1 产物回填

对每个有效 daughter：`A>0,Z>=0` 才写回。新槽位 `IS=ICMAX1+nf`，设置 XJR、ICL、动量、CLKIN 与 IFSTAT=1。主程序从母簇核子列表里顺序挑 NFRA 个，把其 IC 改到新槽位；这是 bookkeeping，并不是 GEMINI 决定了哪一个具体 QMD 核子进入哪个冷核。

随后旧母核 XJR/ICL/坐标/动量置零。不能用退激后的 IC 与原始核子动量重新计算冷碎片内部能量；二者不再是一套自洽微观态。新槽位增长基于 n_frag，包括被过滤的伪粒子可能留下的空槽，需检查容量。

### 9.2 GEMINI 速度到 QMD 动量

GEMINI 返回 `v` 单位 cm/ns、θ/φ 单位度。主程序采用近似：

\[
\mathbf p'_f=A_f(0.938)\frac{v_f}{30}
(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta),
\]

再用 QMD 侧质量表算产物基态质量及能量。这个 p=mv 转换是非相对论近似，不应被描述为完整相对论 p=γmv。

母簇速度则由 `M=0.938*ICL`、`E=sqrt(M²+P²)` 得到 β=P/E，忽略了结合能和 E* 对母质量的修正，然后调用 `LORENTZ(-β)` 把 daughter 静止母核系推回原系。[回填](../anahyp.for#L2667)

### 9.3 LORENTZ 的公式和实参别名问题

该辅助过程实现：

\[
\gamma=(1-\beta^2)^{-1/2},\quad
\mathbf p'=\mathbf p+\frac{\gamma-1}{\beta^2}(\boldsymbol\beta\cdot\mathbf p)\boldsymbol\beta
-\gamma E\boldsymbol\beta,
\]

\[
E'=\gamma(E-\boldsymbol\beta\cdot\mathbf p).
\]

负 β 表示反向 boost。小 β²<10⁻⁸ 时令除法项为 0；超光速/接近光速时把 `1−β²` 截到 10⁻¹⁰，不能把这个保护当作输入物理合法。

主程序 2915 的调用把 `E_THIS,CLKX,CLKY,CLKZ` 同时传作输入与输出，而子程序先写 APZ 后才用 AZ 算 AP0，会因别名让能量读到更新后的动量；编译优化下也存在未定义别名行为风险。应以不同临时变量做输入/输出才能验证洛伦兹不变量，不能因函数公式本身正确就忽略调用方式。[LORENTZ](../anahyp.for#L3991)

### 9.4 “Lab 谱”需审查的具体差异

主程序用 VPR1 做 CM→Lab boost，VPR1 的分子写成：

\[
P_{code}=\sqrt{(npr\,epp)^2+2(npr\,0.938)epp}.
\]

若 epp 为 GeV/u，物理上弹核总动量应为 `npr*sqrt(epp*(epp+2m))`，第二项应随 npr² 缩放；npr=1 时二者一致，重离子入射时不同。后面的 BETA_T 另写了 npr²，但并未替代实际使用的 VPR1。

boost 后 `CLK*` 被改成 Lab 每核子动量，但 `CLKIN` 未同步重新计算。DDCSS 用新的动量方向求 θ，却用旧 CLKIN 分 bin，可能混合两个参考系。CNPCS 的 Z/N/A 计数不直接依赖这些动能，但动能/角谱不能据文件标题直接当作经过验证的 Lab 双微分截面。

## 10. 输出解释与归一化

定义面积 `S=10π X00CUT²` mb。因 1 fm²=10 mb，这个换算是正确的。它假定事件样本对应面积元 `2πb db` 下的均匀圆盘采样；若输入是均匀 b、固定 b、非零 bmin 环带或带权样本，需要使用匹配的权重，不能直接乘同一 S。代码不检验这个前提。

### 10.1 各文件到底包含什么

| 文件 | 当前内容 | 单次模式的计算方式/注意 |
|---|---|---|
| `CNPCS_HOT.dat` | 热簇按 Z、N、A 的三个边缘分布 | `S*count/NEVENT`，单位 mb |
| `CNPCS_COLD.dat` | 最终列表同样三个分布及 10 槽核素计数 | 同上；不是完整二维 σ(Z,N) |
| `IPCS_out.dat` | 所选元素 Z 的中子数分布 | `S*count(Z,N)/NEVENT`，N=0..200 |
| `DDCSS_LNKED_out.dat` | 所选核素、所选角度的能量谱 | 标准分支还输出轻核 dσ/dE；仅 IFSTAT/IFSURF 对象 |
| `ELEKED_out…dat` | 计划中的指定核素动能分布 | 实际最终写 unit 21，且仅 NTH_RUN=1 执行；已知不一致 |
| `TEELECS_out…dat` | 指定 `(Z,N)` 的 HOT/COLD 时间计数 | 实际输出 count/NEVENT，无面积因子，属于平均多重数 |
| `TEAMUL.dat` | 按事件类别分组的最大残余核时间演化 | 时间、E*/A、Amax、Z>MCFR 的簇数；不是总 E* |
| `EXCCROSS.dat` | 末时间片最大簇 E* 的分箱截面；按 A 累加的激发能 | 第一块 mb/MeV；第二块实际 MeV/事件，不是表头写的 mb |
| `GEMINI_INPUT.dat` | 每时间片最大 ICL 簇的槽号、A、Z、J、E* | 未经过真正退激筛选；不一定等于所有 GEMINI 输入列表 |
| `E_SH.out` | ECHECK 的 before/after/差/相对差 | MeV 与无量纲；不是壳修正表 |
| `FSVT.dat` | π−/π+ 与 π 的动能密度等 | 条件性路径，快度二阶矩有错误，见问题附录 |
| `MONITOR…dat` | NNEVENT 的进度记录 | 仅 NTH_RUN≠0；多时间片时不是唯一事件数 |
| `LNVICS.dat` | 速度分布占位 | 当前未实现有效填充/输出 |

CNPCS 的 Z 与 N 输出只到 150，并从 1 开始，省略 Z=0 自由中子/N=0 项；A 输出到 300。数组声明较大不等于文件会完整写出所有值。各边缘分布之和代表带多重数的 inclusive 截面，通常不应直接与单一反应总截面相等。

### 10.2 双微分截面与分箱

标准 DDCSS 输出：

\[
\left.\frac{d^2\sigma}{dE\,d\Omega}\right|_{bin}
\approx\frac{S\,C_{E,\theta}}
 {NEVENT\,\Delta E\,[2\pi\sin\theta\,\Delta\theta]}
\times SPACFAC.
\]

`dkin=1 MeV`，**`dthe=0.1 rad≈5.73°`**，输出角度表用度，先乘 π/180。`SPACFAC` 是展示缩放因子，不是物理效应；比较实验前要除回。

填充用 `ceiling(θ/dthe)`，输出却用 `nint(θ/dthe)`，会造成某些角度对应到前一个 bin。例如 θ=30°，填充 bin=6，输出读 bin=5；不是同一角窗。接近 0° 时 nint 可为 0，而数组下标从 1 开始；sinθ 分母也会退化。接近 π 的最后一格又被 `nthe<mnthe` 排除。

推荐修正后的输出应使用明确区间 `[θlo,θhi]` 和精确 `ΔΩ=2π(cosθlo−cosθhi)`，并让填充/读取使用同一 bin 定义；这是改进方向，当前源码并未实现。

动能 `ceiling(E/ΔE)` 对正 E 对应 `((i−1)ΔE,iΔE]`，输出的是上边界 `iΔE`，不是 bin 中心。大尺度分箱 `DDCSS_LNKED_S=1` 内部按 1、10、100… MeV 重分箱，并不读取样例文件后面的 Split/BinN 列表；该分支也不输出标准分支的 LNKED 部分。[输出分箱](../anahyp.for#L3153)

### 10.3 DIFF_T/B_DIFF 的真正作用

对满足 `0<=DIFFUSION<=DIFF_T` 的候选，抽一个均匀随机数，保留概率为：

\[
P_{keep}=1-\exp\left[\frac{B_{DIFF}}
 {DIFFUSION-DIFF_T-10^{-6}}\right].
\]

DIFFUSION 是退激前的 K/A，不是 E*/A，也不是介质温度。B_DIFF 应为非负且与 K/A 同能量量纲；为负会产生不合理概率。daughter 新槽的 DIFFUSION 保持事件初值 −1，通常不会受到该门控。因此它不是所有冷产物的统一温度修正，也不作用于 CNPCS/IPCS 截面。

### 10.4 EXCCROSS 和时间演化

每个所选时间片寻找最大质量簇；同质量时 `>=` 选最后一个，而 GEMINI_INPUT 的 maxloc 选第一个。TEAE 的 `MAXEXC` 取 `k*dt>=EXPNT` 后、且尚未标记多重碎裂的最大 E*/A，并按 `[EX1DOWN,EX1UP)`、`[EX2DOWN,EX2UP)`、`[EX3DOWN,∞)` 给整个事件分类。

EXCCROSS 第一块只在 `it==kkmax` 填入 `ceil(E*/10)`，每 bin 宽 10 MeV，上限 1000 MeV，零能量有独立的 0 槽。它不是把所有热碎片的 E* 都计进去。第二块 `EXCCROSS2(A)=Σ E*_heaviest`，单次模式除总事件数；不是同一 A 内的条件平均，也不是截面。

TEELeCS 与 XT/YT/ZT 用原始步号 k 直接作下标，上限 300；输出时间却用 k×常量 DELT=0.5，不是输入 dt。若 dt 改变或原始 k>300，输出时间与数组范围都要重新核对。

### 10.5 事件分母与批次合并

单次模式 `NEVENT` 是外层已进入的事件数；选时间窗内累计 HOT/COLD 时会累计多个时间片，但多数归一化只除 NEVENT，没有除时间片数。虽然计算过 `factor=NEVENT*(kkmax−kkmin+1)`，相关输出没用它。

批次 `NTH_RUN≠0` 时 CNPCS、标准 DDCSS/LNKED、EXCCROSS 有意取消 `/NEVENT`，输出的是面积加权计数；而 IPCS、FSVT、粗分箱 DDCSS 仍使用各自原有分母。不能把全部批次输出按一种算法直接平均。

对统一面积且单时间片的 CNPCS 批次，应按：

\[
\sigma=\frac{\sum_r S\,C_r}{\sum_r N_{r,\mathrm{actual}}}.
\]

N_actual 必须是确实读入并参与该统计的独立事件数；多时间片 MONITOR 不满足这个定义。NRMIN 跳过事件时还会跳出时间片循环，若 kktot>1，剩余记录不被消费，下一事件会错位。建议当前只在一个时间片、事件范围和分母已经独立核对时讨论批次运行。

### 10.6 π 谱与误差

FSVT 使用 `pT=sqrt(px²+py²)`、`y=0.5 ln[(E+pz)/(E−pz)]`，π 动能格宽 5 MeV，pT 格宽 0.02 GeV/c，|y/yproj| 格宽 0.02。比值 R=N−/N+，源码用事件多重数的一、二阶矩估计均值误差并乘 1.96，再按独立分子分母传播。

这不是 Poisson bin 误差的简单代用，且源码将 `SNNYPPION` 累加为 `NTRPION²` 而非 `NYPPION²`。此外分子与分母相关性未进入误差公式。使用快度比误差之前必须修正并独立检验；不应把当前输出直接称为经验证的 95% 置信区间。

## 11. 工程使用、调试和扩展时的操作原则

### 11.1 建立基线的顺序

先检查输入协议与 unit 冲突，再验证 `IC/ICL/XJR` 对应关系和 A/Z 守恒，然后验证单核基态的原始内部能量、质量表与激发能；之后才考察 GEMINI 和最终谱。基线至少记录：源文件哈希、编译器版本与选项、全部输入和表文件哈希、事件数、时间片范围、束流系、采样 b 分布、随机种子方案、输出单位和每个输出的实际分母。

每个扫描点独立目录和独立进程。MODULE 全局变量、SAVE 数组、COMMON 和随机流使当前库不是可重入线程安全接口；不能直接 OpenMP 并行调用 event。独立进程还要防止多个进程重建同一 gdat 缓存。

### 11.2 最有价值的断点与监测量

| 断点 | 看什么 | 正常关系/异常征兆 |
|---|---|---|
| 921 后 | 弹靶、束流、dt | 与生成 coll 的输入一致 |
| 1329 后 | 每条粒子记录 | k 一致，静质量/电荷合理，无整行错位 |
| 1812、1836 后 | IC、ICL、XJR、动量和 | 编号变了后统计量必须仍对应 |
| 2219 前后 | AENE、B、未截断 E* | 基态偏差不应靠零截断掩盖 |
| 2560 前 | NEXC 与跳过原因 | 与 real EXCEN 的截断差别清楚 |
| 2644 后 | cn%A/Z/Ex/J、para | 不是只看输入文件注释 |
| 2645 后 | frag、n_frag、weight | ΣA、ΣZ 与母核核子数/电荷一致，光子单独处理 |
| 2915 前后 | 输入/输出四矢量和 CLKIN | 质量壳不变量、动能与角度参考系一致 |
| 最终输出前 | 原始 counts、NEVENT、NNEVENT | 单次/批次/多时间片分母明确 |

### 11.3 参数扫描不应替代缺陷排查

质量谱明显损失一个重残余核，先检查编号重排；程序在启动时 EOF，先检查 unit 18/19；退激尾部对事件顺序异常敏感，先检查重复线形准备；重离子动能谱偏移，先检查 VPR1 与 CLKIN。只有上述条件得到控制后，才把差异归因于 rr0/pp0、E*、J、势垒或能级密度。

详细扫描矩阵、各参数的方向性影响、数值收敛与物理标定分离方法见参数手册。那里列出的扫描幅度是工程诊断示例，未经本项目实验数据拟合，不能作为最优参数引用。

### 11.4 如果需要扩展，哪些改动必须成组进行

扩充粒子数要同时检查 MG、核素直方图 300/200 上限、GEMINI max_frag、候选通道与插值表范围；扩充 l_max_quantum 要同步重建 q_ang/plm2/tl 相关维度；增加轻核 mode 要同步 Z1/A1/S/DEF/width/level/exotic/光学表/非统计分支；新增输出要先定义对象集合、参考系、bin 边界和事件权重，再接入计数器。

这些是维护时应遵循的一致性约束，并非本次已实施的修改。原程序源文件和物理数据在本次说明书任务中保持原样。

## 12. 本说明书的验证边界

本次完成了当前源码的入口、模块、参数、公式和文件流核对，独立编译链接，以及受控人工数据/接口诊断。项目的 check/build/test 顶层脚本仍是占位。完整真实输运事件 `coll.dat` 未随当前工作区数据集提供，现有历史截面文件不能证明它们由当前这份源码生成，也不能反推所有输入状态。

请把[验证与已知问题](exfrag_%E9%AA%8C%E8%AF%81%E4%B8%8E%E5%B7%B2%E7%9F%A5%E9%97%AE%E9%A2%98.md)作为本手册的一部分阅读：其中逐项标明已复现事实、源码确定行为和需要进一步验证的风险。本文公式中简化的理论骨架都给出实现差别；精确接口和未展开的大型拟合系数可从[源码索引](exfrag_%E6%BA%90%E7%A0%81%E7%B4%A2%E5%BC%95.md)跳转。
