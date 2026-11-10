# Visual Studio 编译与运行

适用：Visual Studio 2022、Intel oneAPI Fortran，选择 **x64**。2026-10-21 修复了模块包含顺序、重复编译、定格式行宽及两个程序的输出目录。

## 1. 打开解决方案并选择启动程序

打开 `D:\Spallation\fortran\Console1.sln`。若 VS 已经打开此解决方案，接受工程重载，或关闭后重新打开。

| VS 中的项目 | 编译入口 | 用途 |
|---|---|---|
| Console1 | LQMD.for | 运行输运，生成初级事件数据 |
| exfrag | exfrag.for | 读取已有 coll.dat、initangm.dat，重建激发能并运行 GEMINI 退激 |

在解决方案资源管理器中，右键需要运行的项目，选择“设为启动项目”。已有初级事件、只想比较质量表时，选择 **exfrag**。

工具栏选择 `Debug` 或 `Release`，平台选择 `x64`，再执行“生成 → 重新生成解决方案”。Debug 用于查错，Release 用于正式计算。Win32 配置保留用于旧环境；本轮仅验证 x64 的 ifx 编译。

两个项目只编译各自的入口。`mass_select.for`、`main.for`、`GEMMOD.for`、`anahyp.for` 等由入口通过 INCLUDE 编入，不能再作为独立编译项重复添加。`mass_select.for` 是模块，不是第三个可执行程序。

## 2. 检查工作目录并运行

右键启动项目 → 属性 → 调试，将“工作目录”设为：

```text
$(ProjectDir)
```

在本机就是 `D:\Spallation\fortran\`。程序在此目录下寻找 `Input file`、`gdat`、`tbl`、`coll.dat` 等相对路径。不要把工作目录设为 exe 所在的输出目录。

按 `Ctrl+F5` 直接运行，或按 `F5` 调试。运行前确认输入对应你要计算的事件数；程序会写出结果文件。已有结果需要保留时，应先复制运行目录并在副本中计算；六表比较可直接使用下面的隔离目录脚本。

生成后的路径为：

| 配置 | LQMD | exfrag |
|---|---|---|
| Debug x64 | x64\Debug\LQMD\LQMD.exe | x64\Debug\exfrag\exfrag.exe |
| Release x64 | x64\Release\LQMD\LQMD.exe | x64\Release\exfrag\exfrag.exe |

这些路径相对于 `fortran`。旧的 `x64\Debug\Console1.exe` 不是本轮新生成的程序。

## 3. 设置质量表

修改 `Input file\m_select.in`，例如使用 HFB02 的统一对比模式：

```text
2 ! m_select
1 ! mass_framework
314159 ! comparison_seed
```

编号为 `0=原默认、1=FRDM12、2=HFB02、3=KTUY05、4=WS4+RBF、5=本地LDMM候选质量面`。

- `mass_framework=0` 保留原计算路径，支持 0、1；2–5 明确拒绝。
- `mass_framework=1` 是用户选择的共同退激处方，支持 0–5。固定参考结构，仅比较基态质量带来的阈值和激发能变化。5 的文件来源仍未确认，运行日志会标记。
- 统一模式中 LQMD 使用共同的 M0 产生初级事件；将同一份 `coll.dat`、`initangm.dat` 交给 exfrag，分别选择质量表重放。改变 `m_select.in` 后重新运行即可，无须重新编译。

完整物理边界、缺项回退和比较方法见[质量表统一对比](质量表统一对比.md)与[质量表选择](质量表选择.md)。

在已加载 Intel oneAPI 环境的 PowerShell 中，进入仓库根目录，使用已经生成的 Release exfrag 比较六组。可先打开 Intel oneAPI 的 Intel 64 命令提示符，再输入 `pwsh`；普通 cmd 也可先执行：

```bat
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
pwsh
```

然后执行：

```powershell
./scripts/compare_mass_tables.ps1 `
  -SourceDirectory 'D:\Spallation\fortran' `
  -Executable 'D:\Spallation\fortran\x64\Release\exfrag\exfrag.exe' `
  -OutputDirectory 'D:\Spallation\tmp\mass-compare-vs-run1' `
  -Models 0,1,2,3,4,5 -Seed 314159
```

输出目录必须尚不存在。脚本复制初级事件及所需输入到每组目录，统一设置 `mass_framework=1` 和种子，运行 exfrag，并核对输入哈希与正常结束标记。源目录要有完整、相互匹配的 `coll.dat`、`initangm.dat` 和输入文件；批量运行时间取决于事件数。

## 4. 本次编译错误的原因

- `#7013 MASS_SELECTION`：源码目录残留 gfortran 生成的 `.mod`，Intel ifx 无法读取。已将确认的这类文件归档到 `tmp\task012\round5-vs\`，保留原件可追溯。以后跨编译器检查须使用独立的模块输出目录。
- `#6580`、`#6406`：在模块无效、定义顺序错误时可产生连锁错误；将 INCLUDE 片段重复列为编译项也会重复定义。LQMD 现由 wrapper 先包含模块，再包含主程序，两个项目各只有一个编译入口。
- `#5082`、`LINE_`、`LINE_N` 等：定格式源被默认 72 列截断。项目已设置 `/fixed /extend-source:132`。

两项目的 `.obj`、`.mod` 分别输出到各自配置目录的 `obj` 子目录，避免互相覆盖。若自行修改了 VS 属性，请检查以上设置后重新生成。
