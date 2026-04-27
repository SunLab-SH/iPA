## iPA — Copilot / AI agent quick instructions

目标：帮助 AI 代码代理快速理解、修改并贡献 iPA（Integrated Processing and Analysis Toolkit）项目。

简短“契约”（Inputs / Outputs / 成功标准）
- 输入：修改点通常是 `ipa/` 下的模块、`examples/` 脚本或 `docs/` 文档；数据位于项目根的 `data/` 目录。
- 输出：功能变更应包含可运行的例子或小测试（参见 `examples/`），并保持包可通过 `pip install -e .` 本地安装。
- 成功：所有新增 Python 文件通过基本导入（`python -c "import ipa"`）且文档或示例能在本地环境复现。

大局架构（快速概览，需读多处文件以掌握细节）
- 包根：`ipa/`，主要子包为 `processing/`（图像处理管线）、`analysis/`（定量分析）和 `data_loader/`（不同显微镜文件格式的读取器）。
- 示例与入口：`examples/` 下有多种演示脚本（例如可视化和 RDF 分析），这些脚本是快速验证改动的首选场景。
- 配置与资源：模型权重和大文件使用 Git LFS；示例数据放在项目根的 `data/`（见 `README.md` 的数据说明）。

关键文件与目录（常引用）
- `README.md`（`iPA/README.md`）：安装、数据组织和总体流程说明，优先阅读。
- `setup.py` 与 `requirements.txt`：安装和 extras（如 `ml`、`fileformats`）定义了可选依赖。
- `ipa/processing/`、`ipa/analysis/`、`ipa/data_loader/`：主要代码位置，修改这些目录内逻辑时请添加或更新对应示例。
- `examples/`：运行示例以验证改动。
- `workflow_images/figure_1_v25.jpg`：项目流程图，帮助理解模块边界。

开发 / 构建 / 测试 快速命令（以项目根 `iPA/` 为工作目录）
- 推荐环境：Conda, Python 3.7 或 3.8（见 `README.md`）。
- 安装：
```powershell
conda create -n ipa_env python=3.7
conda activate ipa_env
pip install -r requirements.txt
pip install -e .
``` 
- 若需全部可选依赖： `pip install -e .[all]` 或安装 `setup.py` 中的 extras（如 `ml` 包含 PyTorch）。
- Git LFS：项目依赖大文件（`.pth` 等），执行 `git lfs pull` 拉取模型权重。
- 文档：使用 Sphinx（docs/）构建，通常 `make -C docs html` 或 `sphinx-build -b html docs/ docs/_build`。
- 测试：项目使用 `pytest`（在 `extras_require['dev']` 中）。新增逻辑时至少提供一个快速 pytest 用例。

项目约定与习惯用法（可被 AI 代码代理依赖）
- 数据位置固定：脚本假定 `data/` 在项目根并包含 `cryoET/`, `sim_images/`, `sxt_images/`, `wfm_images/`。不要在示例中硬编码绝对路径。
- 轻量可复现：大改动应配合 `examples/` 中的最小运行脚本或 notebook 来证明行为。
- 依赖分组：将重型 ML 依赖放在 `ml` extras，许多工具仅在需要时安装（避免在 CI 中一键安装所有 heavy deps）。

集成点与外部依赖
- 模型权重：通过 Git LFS 托管（注意在 CI 或新环境中执行 `git lfs pull`）。
- 文件格式支持：项目通过 `readlif`, `czifile`, `nd2reader`, `aicsimageio` 等处理专有显微镜格式；对这些模块的更改会影响 `ipa/data_loader/`。
- 可视化：使用 `plotly`, `matplotlib`, `seaborn`；交互式示例通常在 `examples/` 可见。

具体示例（说明常见修改场景）
- 新增数据加载器：在 `ipa/data_loader/` 下添加模块，导出统一接口（函数名/类名应与现有 loader 保持一致），并在 `examples/` 添加小脚本演示读取和可视化一帧。
- 改进处理算法：在 `ipa/processing/` 中修改函数后，运行 `examples/` 中对应脚本或添加 pytest 短测试来验证输出形状与类型。
- 文档与图像：若更改了 pipeline，更新 `workflow_images/`（若需要）并在 docs 中更新相应章节。

常见边界/注意点（edge cases）
- 大文件/模型未拉取：在 CI 或新环境执行前先 `git lfs pull`，否则示例/测试会因缺少 `.pth` 文件失败。
- 平台差异：开发推荐 Linux/macOS；Windows 下某些 `make` 或 shell 命令需要用 PowerShell 或等效命令替换。
- Python 版本：代码以 3.7-3.8 兼容为主，避免使用 3.10+ 的新语法，除非同时更新 `setup.py`。

验证与交付（AI 代理应完成的最小交付）
1. 变更代码后，确保 `python -c "import ipa"` 不报错。
2. 提供或更新至少一个 `examples/` 脚本来复现改动（或一个 pytest 测试）。
3. 若涉及模型/大文件，说明需要执行 `git lfs pull` 并在 PR 描述中标注。

如果本文件中有遗漏或不准确的地方，请指出要补充的区域（例如：常用脚本名、CI 细节或特定模块的约定）。

----
最后一步：请告诉我是否要把此文件提交到仓库的其他分支，或在合并时保留现有项目级说明（如果有）。
