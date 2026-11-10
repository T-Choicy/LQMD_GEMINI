# 运行手册

## 1. 首先跑一个小样

完成编译后从仓库根目录运行：

```powershell
./scripts/run_example.ps1
```

默认质子+232Th、80MeV实验室入射动能，模型1(FRDM12)、统一质量框架1、SO(3)取向，1片×2事件。脚本复制全部必要输入和数值表到新运行目录的 `source/`，再用 `examples/p232Th-80MeV/` 中的物理模板覆盖工作副本。原始 `fortran/` 和例子不会被运行改写。

只准备输入、不启动物理计算：

```powershell
./scripts/run_example.ps1 -PrepareOnly
```

每次使用新目录；准备检查和正式运行是两次独立调用，会生成不同目录。支持 `-RunDirectory <不存在的路径>`。已存在目录会被拒绝，失败运行保留日志，不重用它掩盖错误。

## 2. 五个并行程序、10000事件

```powershell
./scripts/run_example.ps1 -Shards 5 -EventsPerShard 2000 -Model 1
```

这是单个质量表5片，共10000个输运事件，之后每片退激并按事件数加权汇总。Model可选1–5；每用户共享上限25个计算槽，不能把5片当作50000事件。默认种子基值314159，分片种子由控制器派生，保存于manifest。

## 3. 改物理参数

物理模板在 `examples/p232Th-80MeV/`：`LQMD0.in`、`anahyp0.in`、`CNPCS_in.in`。LQMD能量在第4数据行的 `epp`，单位 **GeV/u**，质子40/50/60/80MeV分别为0.040/0.050/0.060/0.080。这里是入射动能，不是人为指定复合核激发能。模板位置式读取，不要删标题、随意换行或改变字段数。

改变例子后提交版本再正式运行。质子能量可用 `./scripts/run_example.ps1 -EnergyMeV 40 -Shards 5 -EventsPerShard 2000` 指定，脚本将输入转换为0.040GeV/u，并同步目录名与run.json；不要手改物理能量却留下错误的元数据。当前仅参数化p+232Th，换反应体系应建立独立例子。根fortran模板保留800MeV参考参数，但控制输入默认仅2事件，不代表80MeV例子。

通用控制器读取某个工作区的 `Input file/INPUT.in`。七个有效数据行依次为NFILE、NRUN_PER、NAME、T_1ST、T_PER、CALL_ANA、FAIL；CALL_ANA=1自动分析，FAIL须为0。通过控制器运行时物理参数来自 `LQMD0.in` / `anahyp0.in`，直接物理程序才读取 `LQMD.in` / `anahyp.in`。`m_select.in` 三个有效值依次为质量表、质量框架、种子。

## 4. 五张质量表批跑

先在仓库 `fortran/Input file/` 设置物理模板和控制文件：NFILE=5、NRUN_PER=2000、唯一NAME、CALL_ANA=1、FAIL=0，质量框架=1。然后：

```powershell
./scripts/run_windows_mass_models.ps1 -RunDirectory ./runs/five-models-80MeV-10k
```

先检查可加 `-PrepareOnly`。这会为1–5各建一组独立工作区、每表5个LQMD进程，总计最多25个，各表10000事件；输运基态均冻结M0，分析按所选质量表。各表配对种子，不能把它们当成50000个独立物理初态样本。

后端JSON五表路径则是一组共用轨迹、五张表分析，和五组输运批跑组织方式不同。已有完成的五表运行，可使用 `import_windows_kinematics.ps1 -PreviousRunDirectory <旧批次> -RunDirectory <已准备的新批次>`：脚本核对布局和输入，不接受任意coll文件拼接。外部coll至少还需要匹配的initangm、完整事件范围、头部、原物理输入和归一化信息，不能凭文件名认定一致。

## 5. 中断与失败

看总控日志、各片 `state.json`、`program.stderr.log`、结束标志与输出完整性。控制器拒绝输入/程序哈希变更的原运行目录；修改参数后开新目录。自动衔接分析失败时保留输运结果，在相同工作区和未改输入的前提下使用ctrlA重试；成功分片会被验证后复用。禁止手动将failed改为completed。

每次正式实验至少保存：Git标签及提交号、dirty状态、编译器版本、exe SHA、输入和质量表SHA、事件数、分片、随机种子、运行日志、轨迹和汇总结果。`run.json`与控制器manifest提供其中大部分；GitHub不替代课题组数据备份。
