# 课题组版本与协作规则

## 使用一个确定版本

VERSION当前为0.1.1，对应不可移动的注释标签v0.1.1。正式计算优先checkout标签并重新编译，保留Git提交号和程序SHA。只有版本号而没有提交号，不足以识别未发布改动。

```powershell
git fetch --tags origin
git switch --detach v0.1.1
./scripts/build_windows_controllers.ps1 -Rebuild
```

不要把main上的后续修复自动替换到正在运行或已发表的实验。旧标签永远保留，不以force更新；旧计算结果不覆盖，必要时以新版本新目录重算并说明差异。

## 修改代码

```powershell
git switch main
git pull --ff-only
git switch -c feature/简短主题
# 修改、运行相关测试与小样、记录物理影响
git add 具体文件
git diff --cached
git commit -m "fix: 具体问题及修复"
git push -u origin HEAD
```

发起PR，经组内另一位成员审查后合并main。避免直接在main上多人改写；不提交exe、轨迹、账户凭据或本机路径。物理公式、激发能、MST判据、质量表、随机数或归一化变更须明确写出旧/新行为与验证范围，不能以“重构”隐藏物理变化。维护者应在GitHub规则设置中为main启用所需审查；本文件是流程约定，只有服务端实际启用后才是技术强制门禁。

## 发布新版本

采用MAJOR.MINOR.PATCH。0.x阶段仍在完善研究接口：破坏输入/输出或更改重要物理默认值提升MINOR；兼容修复提升PATCH，并在CHANGELOG标明它是否可能改变产额。版本号不能代替物理描述。稳定公共契约后再考虑1.0.0。

发布前更新VERSION、CHANGELOG、PROVENANCE，跑仓库检查、Python测试、真实编译、控制器模拟测试及受影响物理回归。提交后在干净工作区生成不可移动标签，再推送：

```powershell
git tag -a v0.1.2 -m 'v0.1.2: 发布内容与验证范围'
git push origin main
git push origin v0.1.2
```

这里的0.1.2只是后续发布示例，不能在未完成变更时照抄创建。发现已发布标签有问题时发布新标签，不删改旧标签。更新PROVENANCE时按其生成脚本重算载荷SHA；不要将提交自身哈希嵌入自身文件形成循环。

每次正式计算保存版本、提交号、dirty状态、输入/质量表/程序SHA、编译器、事件数、种子、分片和归一化。脏工作区只用于调试；若确需研究试算，归档补丁并标明未发布，不能称它等同某个标签。
