# 怎样查看和解释结果

## 先找到文件

`run_example.ps1` 最后打印Results与run.json路径。典型布局：

```text
runs/<体系能量模型>/<运行标识>/
  run.json                         版本、提交、脏状态、能量、事件、种子、程序SHA及退出码
  controller.stdout.log            总控日志
  controller.stderr.log            总控错误
  source/results/calculation/
    manifest.json                  分片方案、输入哈希、模型和阶段状态
    kinematics/shard-N/attempt-N/   原始coll、initangm及输运日志
    analysis/model-M/shard-N/attempt-N/  退激输出与日志
    merged/
      coll.dat / initangm.dat       合并输运轨迹
      model-M/CNPCS_HOT.dat         MST初级热碎片分布
      model-M/CNPCS_COLD.dat        GEMINI后冷产物分布
      distributions.csv            完整A/Z/N分布，含零bin
      target_isotopes.csv          五个目标核素截面
```

五表批跑在外面另有model-1至model-5工作区及batch-manifest/status，汇总目录在批跑根下。不要混用单表运行和五表父目录。

## 查看225Ac及其它核素

```powershell
$merged='你的运行目录/source/results/calculation/merged'
Import-Csv "$merged/target_isotopes.csv" | Format-Table model_name,isotope,cross_section_mb,events
Import-Csv "$merged/target_isotopes.csv" | Where-Object isotope -eq '225Ac'
```

列为model、model_name、isotope、Z、N、A、cross_section_mb、events。这里直接给出程序归一化后的毫巴截面，不要再除一次事件数。默认五目标为225Ac、227Ac、227Th、225Ra、223Ra；必须在CNPCS_in.in中恰好配置一次。不是按文本行号猜核素。

`distributions.csv` 的stage是HOT/COLD，variable是A/Z/N，bin为对应整数，cross_section_mb为该bin总截面。A分布同一个A下会包含多个元素，不能把A=225当作225Ac。HOT与COLD描述两个阶段，不能相加当作最终产额。

## 画分布图

安装Python绘图依赖后：

```powershell
python figure/scripts/plot_cross_section.py "$merged/model-1/CNPCS_COLD.dat" --yscale log --output ./figure/output/cold.png
```

对数轴无法显示零值；图中没画出的零bin不等于可以无证据外推阈值。CSV保留零值。针对225Ac跨能量比较，应从每次运行的target_isotopes.csv取同一模型、同一口径的数值并保留各自events和版本。

## 归一化和统计

控制器按分片事件数加权：σ总=Σ(σ片×N片)/ΣN片，不能简单累加各片截面。当前原生文本有打印精度限制，合并不会恢复被四舍五入掉的有效位数。

在本例固定截止bcut=8.25fm且使用既定几何采样归一化时，几何截面为π×bcut²×10mb，约2138mb；若每个事件至多有一个指定重核，计数近似为σ/σ几何×N。多个碎片/其它采样不能直接套这个二项解释，应使用事件记录及实际权重。2事件仅测试链路，不能用它判断225Ac产额。零计数也不证明真实截面为零。

CNPCS_COLD是程序退激后的直接残核输出；没有在本仓库自动计入所有后续放射性母体衰变、照射/冷却时间或实验化学收集效应。因此与实验累积截面比较必须先统一口径。改变质量表会改变分离能/Q和重建激发能，不等于验证了整套模型。

## 判断运行是否完整

检查run.json退出0、manifest的L/A/merge均completed、各片state退出0，以及程序结束标志和所需文件。单纯“有CNPCS文件”不足以认定本次成功。`MST_ACCOUNTING.dat`记录逐事件预期A/Z及HOT/COLD A/Z；出现不闭合会使当前程序非零退出。账目闭合证明记账一致，不能单独证明所有物理假设正确。

coll.dat是输运相空间轨迹，不是质量表；initangm.dat应与同一批事件配套保存。不要从不同运行中拼接二者，也不要通过修改头部伪造事件数。重新分析需原始输入、质量表和种子来源，不应仅保存最后CSV。
