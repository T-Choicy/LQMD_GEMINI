# 四控制器的共同流程（TASK-029）

`ctrlL.for` / `ctrlL_win.f` 调用同一 L 阶段，`ctrlA.for` / `ctrlA_win.for` 调用同一 A 阶段。四文件均包含 `controller_dispatch.for`，调用现有 `scripts/windows_event_controller.ps1`，不再各自维护模板解析、分片和汇总代码。物理主程序、分析程序和输入模板没有改变。

原 `.in` 入口读取 `INPUT.in` / `INPUT_ANA.in`，按 `m_select.in` 指定的 1–5 中一张表运行；保持 Windows 版既有范围，不能用这个入口选择 0。JSON 配置的 `Models=[1,2,3,4,5]` 在共用轨迹上逐表分析。Windows 的五组输运批跑仍使用 `run_windows_mass_models.ps1`，详见 [原使用说明](RUNNING.md)。这两种五表流程不能当作同一种统计抽样。

每片使用私有输入和输出，完整复制模板后只调整事件、种子及选表等控制字段。分析 Models 行按当前程序的两字段格式生成，保留 `%%Gemini settings` 后的完整内容，包括后加的 GEMINI 参数。`CALL_ANA=1` 在成功输运后自动进入分析。结果仍在 `results/<NAME>/`，包括 `manifest.json`、分片日志及 `merged/` 中的截面与来源记录。旧的共享 `LQMD.out` / `anahyp.out` 不再用于判定完成。

Windows 使用系统 PowerShell 5.1；Linux 使用 **PowerShell 7 (`pwsh`)**，需要先安装到 PATH，另需预编译的 LQMD 和 exfrag。普通版因此不再是仅依赖 Fortran 的独立调度器。Linux `.in` 入口优先使用输入源目录的 `LQMD` / `exfrag`，否则使用项目 `fortran/bin/` 中的同名程序；JSON 可明确提供路径。不能把 Windows `.exe` 当作 Linux 程序使用。

Linux 编译示例（仅编译控制器，不编译物理程序）：

```sh
cd /path/to/Spallation
mkdir -p fortran/bin
gfortran -ffixed-line-length-132 -I fortran fortran/ctrlL.for -o fortran/bin/ctrlL
gfortran -ffixed-line-length-132 -I fortran fortran/ctrlA.for -o fortran/bin/ctrlA
./fortran/bin/ctrlL --prepare-only
# 确认准备的私有输入后，去掉 --prepare-only 运行；CALL_ANA=1 自动分析。
```

Intel 编译器也可使用相应固定格式及 132 列选项。Windows 按原构建脚本或项目编译，公共 include 路径为 `fortran/`。入口应从项目树内运行，或通过含路径的程序名启动，以便找到共用脚本。支持显式配置路径、`--help`、`--prepare-only`、`--merge-only`；错误返回非零。含空格路径受支持；Linux 入口路径不支持单引号，Windows 入口路径不支持双引号或百分号，会明确拒绝。

Windows 保留原有共享命名信号量；Linux 用本用户 `~/.local/state/spallation/compute-slots/` 下的 25 个内核文件锁，进程退出后释放。限制针对本用户使用这些控制器启动的计算进程，不限制手工启动的程序。运行期间不要删除或替换该锁目录。

公开版测试入口和覆盖边界见 [测试说明](TESTING.md)。模拟测试只验证控制流程与数据组织，不能认证核反应产额。
