# 安装与编译

## Windows 标准环境

- Git；Windows PowerShell 5.1（系统自带）。运行脚本建议使用上层示例给出的 `powershell -NoProfile -ExecutionPolicy Bypass -File ...`，不要修改全局执行策略。
- Intel oneAPI Fortran `ifx`，以及可被其识别的 Visual Studio C++ 链接环境。该版本在 Intel 2025.3 / VS2022 上实际构建验证。
- Python只用于查看、绘图和Python测试，不参与Fortran物理计算。绘图支持Python3.9+，依赖见 `figure/requirements.txt`。

从仓库根目录执行：

```powershell
./scripts/build_windows_controllers.ps1
```

生成 `fortran/x64/Release/win_batch/` 内的 `LQMD.exe`、`exfrag.exe`、`ctrlL_win.exe`、`ctrlA_win.exe`。不要只编译 `anahyp.f` 或 `GEMMOD.f`：退激的入口是 `exfrag.f`，模块顺序由include保证。固定格式需要132列，构建器已设置相应选项。

已有四个程序时构建器跳过；**修改源码或切换版本后必须重建**：

```powershell
./scripts/build_windows_controllers.ps1 -Rebuild
```

重建先将旧产物归档至 `tmp/build/windows-controllers/`，检测到活动物理程序时拒绝重建。不要在运行中替换程序。构建命令保存在这个目录的 `compile_windows_controllers.cmd`。构建出错会返回非零；查看终端编译错误及批处理，不能只看是否存在旧exe。

若VS装在非标准位置，先在当前PowerShell会话设置实际目录，例如：

```powershell
$env:VS2022INSTALLDIR='F:\Program Files\Microsoft Visual Studio\2022\Professional'
./scripts/build_windows_controllers.ps1 -Rebuild
```

也可通过 `-Compiler <实际ifx.exe> -SetVars <实际setvars.bat>` 指定编译器；该路径属于你的机器，不要提交到共享脚本。

## 绘图与软件测试

```powershell
python -m venv .venv
./.venv/Scripts/python.exe -m pip install -r figure/requirements.txt
./.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
./.venv/Scripts/python.exe scripts/check_repository.py
```

这组Python测试仅验证解析和绘图。Fortran调度的模拟测试见 [测试说明](TESTING.md)，物理小样见运行手册。没有把占位的check/build/test脚本作为测试通过凭据。

## Linux 范围

Linux需要PowerShell7 `pwsh`；普通控制器可用gfortran固定格式132列编译。物理程序仍需另配兼容编译器，不能直接运行Windows exe。本版本没有承诺gfortran能完整编译既有Intel扩展物理代码。详细控制器机制见 [四入口细节](controller_details.md)。Windows入门脚本 `run_example.ps1` 使用Windows程序和隐藏窗口，不作为Linux运行入口。
