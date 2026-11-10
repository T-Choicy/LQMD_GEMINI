# figtools

`figtools` 是面向 CASCADE/anahyp 多块 ASCII `.dat` 输出的轻量绘图工具。
运行时只依赖 Python 3.9+、NumPy 1.24+ 和 Matplotlib 3.7+。

## 安装

```powershell
python -m pip install -r figure/requirements.txt
```

也可以在项目根目录执行 `python -m pip install -e figure`。

## 生成第一张图

```powershell
python figure/scripts/plot_cross_section.py 你的运行目录/source/results/calculation/merged/model-1/CNPCS_COLD.dat
```

默认输出到 `figure/output/CNPCS_COLD_cross_section.png`，生成一行三列的 Z、N、A
截面分布图。可用 `--yscale log` 使用对数纵轴，用 `--quantities Z A` 选择面板，
用 `--output path/to/result.png` 指定输出路径；未指定 `--title` 时会从输入路径推导标题。
脚本会打印实际选用的字体：优先 Arial，依次降级到 Helvetica、Liberation Sans、
DejaVu Sans。

## 多能量对比图（PPT）

为每个物理量分别生成一张会议报告用的多能量对比图：

```powershell
python figure/scripts/plot_energy_comparison.py --root data/p+232Th --mode COLD
```

默认生成 `figure/output/p+232Th_Z_energy.png`、`p+232Th_N_energy.png` 和
`figure/output/p+232Th_A_energy.png`。每张图叠加按能量升序排列的 80、100、200、
300、600、1000 MeV 曲线，默认使用 log 纵轴、`cividis` 顺序色带、纯线和 PPT 字号。
可用 `--quantities Z A`、`--output-dir path/to/output`、`--yscale linear` 和
`--title` 调整输出。扫描器允许发现其它合法能量目录，但本 PPT 入口只绘制上述确认的六个
能量。

## 扩展一种图

在 `figtools/io.py` 中为新数据块增加规范化键和必要的数据转换，在 `figtools/plots.py`
中新增一个接收 `Result`、自建 `matplotlib.figure.Figure` 并返回 `Figure` 的函数，
再在 `figure/scripts/` 增加薄的 argparse 入口。不要使用 `matplotlib.pyplot` 的全局
绘图 API；这样同一进程中连续生成的图可以互不干扰。

## 测试

```powershell
python -m unittest tests.test_io -v
python -m unittest tests.test_plots -v
```

顶层 `scripts/check.ps1`、`build.ps1`、`test.ps1` 仍是项目骨架中的占位入口，不能将其
成功退出视为真实检查、构建或测试通过。
