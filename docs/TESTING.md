# 测试与验证覆盖

## 无需历史数据的Python测试

```powershell
python -m pip install -r figure/requirements.txt
python -m unittest discover -s tests -p 'test_*.py' -v
python scripts/check_repository.py
```

解析测试使用明确标记的合成CNPCS夹具，绘图测试使用内存数组，不依赖data/或tmp/taskNNN。仓库检查核验版本、载荷清单、禁止提交的路径及大小、include和本地文档链接。

## Windows原生控制器模拟测试

需要Intel ifx及链接环境：

```powershell
./tests/windows_controller_smoke.ps1
./scripts/build_windows_controllers.ps1
./tests/windows_legacy_controller_smoke.ps1
```

第一套验证分片事件、种子、五表加权汇总、损坏/失败拒绝及恢复；第二套验证五表各5个模拟程序和自动分析。日志保存在tmp/tests下，模拟器没有核反应物理含义。非标准编译器路径按脚本参数指定；VS路径按安装文档在当前会话配置。

编译后可再运行 `./tests/test_example.ps1`，验证例子准备、能量与元数据同步、输入保护和错误参数拒绝，不额外跑物理事件。

Linux `pwsh -File tests/controller_slots.ps1` 检查25个计算槽的本进程及跨进程互斥与释放，测试前不得有其它同用户的控制器计算在运行。它不是Linux物理计算验证。

## 真实程序小样

编译后执行 `./scripts/run_example.ps1 -Shards 1 -EventsPerShard 2`。检查run.json、manifest、各阶段退出码、CNPCS及MST_ACCOUNTING；再决定正式统计量。发布v0.1.0的验证摘要放在VALIDATION.md。

改动影响物理公式、激发能、随机数或守恒时，仅跑上述软件测试不够，应另行设计对应轨迹回归并归档数据。没有把仓库CI里的Python绿色状态称为Fortran物理认证。
