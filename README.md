# LQMD_GEMINI

课题组使用的 **LQMD 输运 → MST 碎片识别 → GEMINI 统计退激 → 截面汇总与绘图** 代码。当前版本 **v0.1.1**，标准可复现入门流程为 Windows / Intel oneAPI Fortran。Linux 控制器已有独立测试，但本版本未认证 Linux 上完整物理链。

## 从这里开始

1. [安装与编译](docs/INSTALL.md)：需要什么软件、怎样编译、编译失败看哪里。
2. [运行手册](docs/RUNNING.md)：先跑2事件，再跑5片×2000事件，单表、五表和已有轨迹分析的区别。
3. [数据说明](docs/DATA.md)：结果位置、文件和列的含义、查看225Ac、画对数截面图、归一化与统计边界。
4. [版本与组内协作](CONTRIBUTING.md)：分支、提交、审查、标签、实验记录和历史版本复现。
5. [上传与归档范围](docs/REPOSITORY_POLICY.md)：什么进Git、什么放计算目录、什么应单独归档。
6. [模型与已知限制](docs/PHYSICS.md)、[来源与署名](NOTICE.md)。

```powershell
git clone https://github.com/T-Choicy/LQMD_GEMINI.git
cd LQMD_GEMINI
git switch --detach v0.1.1
powershell -NoProfile -ExecutionPolicy Bypass -File ./scripts/build_windows_controllers.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File ./scripts/run_example.ps1
```

最后一条运行质子入射232Th、80 MeV、1片×2事件、FRDM12/统一框架，自动接续退激。它验证安装与文件链，不给出可靠产额。每次生成独立 `runs/p232Th80-model1/<运行标识>/`，不会覆盖旧运行。若要改束流能量，按运行手册修改例子的模板，并保持单位一致。

## 仓库结构

```text
fortran/                  物理源码、工程、必要输入和数值表
mass_tables/              保留原格式的质量表来源副本
examples/p232Th-80MeV/    可直接使用的小样物理模板
scripts/                  通用构建、调度、运行与分析入口
tests/                    可独立运行的软件测试与合成夹具
figure/                   NumPy/Matplotlib绘图工具
docs/                     使用、数据与模型说明
VERSION / CHANGELOG.md    版本号及变更记录
PROVENANCE.json           发布文件来源和字节哈希
runs/                     本机运行与结果，自动生成，不进Git
tmp/                      构建/测试中间文件，不进Git
```

历史工作区中的大规模轨迹、运行结果、临时环境、旧聊天记录、依赖私有归档的实验脚本未上传。测试里的临时路径是运行时工作目录，不是需要上传的数据。测试目录不包含实际临时产物。
