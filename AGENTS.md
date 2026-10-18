## ATST-Tools developing guide

### 简介

ATST-Tools (ASE Transition State Tools for ABACUS and ML potentials)：建立用 ASE 等科学计算 Python Package 调用 ABACUS / DeePMD-kit 作为第一性原理/机器学习势计算后端，完成高阶科学计算工作流的封装 Python Code。

### 开发要求

- 代码设计需遵循Zen of Python原则。
- 代码库尽可能集成和封装，CLI设计尽可能兼顾易用和可扩展，仓库核心代码架构需要具备足够可扩展性。
- Unit Test覆盖度足够且粒度合适，Example中需要覆盖项目各方面功能并作为用户快速上手入口。
- 核心代码库各个函数需要具有精练且完整的，Google Style的docstring。
- 程序性输出（stdout/stderr、warning、logging、异常消息、CLI help/错误）统一英文；文档保持现有语言，新增或修改的注释与 docstring 使用英文。
- 项目一段开发任务结束后，在 `docs/` 的合适位置完成更新；维护规则见
  `docs/developer/DOCUMENTATION.md`。

### 仓库设计
- 优化 user interface，使用 CLI + YAML inputs 进行工作流交互，并在examples/目录下留下示例交互模式。
- 对于 ABACUS 的 ASE interface，ATST-Tools 以独立安装的 abacuslite 为唯一 ASE calculator backend。ATST-Tools 的基本定位是 abacuslite CLI wrapper，通过 abacuslite 完成 ABACUS 相关设置，并配合 ASE 开展计算任务；尽可能不重写 ASE 的已有实现，做到 ASE-native。
- temp_repos 不同步到 git 仓库，其中放置开发时可参考代码仓库（不保证存在于任意 checkout；依赖它的快照 drift-check 需先准备对应上游树）。

### 版本号语义

- minor（2.x.0）：保留给阶段性/重大功能发布——新增独立功能点，或仓库渐进开发一段时间后的整体发布。
- patch（2.2.x）：承担小功能加入与新增优化，以及 bug 修复。
- 版本号在 `pyproject.toml` 单一版本源维护；发版时同步新增
  `docs/releases/RELEASE_NOTES_<版本>.md`，并更新功能矩阵
  `docs/reports/FEATURE_STATUS_MATRIX.md`。
- 每个发布须在 release notes 中说明向后兼容性；仅纯增量（不破坏既有
  YAML/CLI/API 契约）的功能与修复可进入当前 minor 或 patch 序列。

### 维护验证环境（SAI）

- SAI 相关 module、队列和软件版本仅用于维护者验证证据，不构成终端用户的运行前提。服务器上有可直接使用的ABACUS LTS 3.10.1 和 DeePMD-kit 3.1.3。
- image-level NEB/AutoNEB 已受支持；真实并行验证需在具备站点兼容 MPI 的 Python 环境中进行，并使用新创建的 conda 环境 atst-dev。
- 开展调用ABACUS和DeePMD-kit的测试需要将任务通过slurm脚本交到4V100节点上，使用GPU节点计算。对于ABACUS，你需要在INPUT中设置ks_solver cusolver (在默认的basis lcao下)。

### 基本边界
- abacuslite 是项目的 ABACUS-ASE backend，`pyproject.toml` 声明
  `abacuslite>=1.1.0,<2`。独立包由 `QuantumMisaka/abacuslite` 维护；缺包或已安装包
  损坏均应失败，不回退到旧实现。`src/atst_tools/external/ASE_interface/abacuslite`
  仅保留薄转发模块供 Paimon 等既有调用方使用，计算器和 parser 类对象仍来自
  独立包。原 vendored 基线、补丁登记和双环境验收记录是历史审计证据，
  `ABACUSLITE_SNAPSHOT.md`/`PATCHES.md` 不再驱动在线 drift CI 或运行时选择。
- ase-abacus 是 legacy ABACUS-ASE backend 参考基线，位于 temp_repos/ase-abacus
  （若本地存在）；main 分支运行时不再使用，仅作功能对照，不能参与项目开发。
- temp_repo下存放有本项目的可参考代码库，该目录下内容不进入git仓库。
- 可拓展基于ase的分子动力学计算功能，并为其他的基于ase的模拟功能提供可扩展设计。

### 用户快速使用
- `examples` 目录下包含项目的快速上手案例。
- `docs/skills/atst-cli/SKILL.md` 中包含 atst-cli 的快速使用说明

### 文档治理入口

开发任务结束前必须按文档治理机制更新对应入口。先判断本次变更影响的是用户路径、开发者路径还是项目管理路径，再同步长期文档、功能矩阵和发布说明。

#### 用户入口
- `README.md`：项目目标、支持 workflow、快速开始、参数入口和状态入口。
- `docs/index.md`：用户、开发者、项目管理者三条阅读路径。
- `docs/user/USER_GUIDE_CN.md`：中文 10 分钟快速上手、ABACUS/DP 注意事项。
- `examples/README.md`：示例学习路径和可运行配置说明。
- `docs/user/CONFIG_REFERENCE.md`：手写 YAML 语义参考。
- `docs/user/YAML_INPUT_VARIABLES.md`：由 schema 生成的 YAML 参数总表。

#### 开发者入口
- `docs/developer/HANDOVER.md`：维护者日常 checklist，新增 workflow、YAML 字段、CLI、backend、example、report、release 时先查这里。
- `docs/developer/DOCUMENTATION.md`：文档入口、目录职责、归档规则和检查命令。
- `docs/developer/YAML_INPUT_GOVERNANCE.md`：YAML schema、生成参数文档和测试治理。
- `docs/reports/FEATURE_STATUS_MATRIX.md`：当前功能支持状态。
- `docs/superpowers/specs/`、`docs/superpowers/plans/`：设计 spec 与待执行计划；
  完成后吸收结论并移入 `docs/archive/`。

#### 变更后检查
- 文档-only 变更至少运行：
  ```bash
  git diff --check -- README.md docs examples/README.md AGENTS.md
  rg -n "^<<<<<<<|^=======|^>>>>>>>" README.md docs examples/README.md AGENTS.md
  ```
- 运行 `conda run -n atst-dev python scripts/check_docs_governance.py`
  （冲突标记、相对链接、入口文件与 HTML 基础解析）。
- 修改 YAML schema 时，重新生成 `docs/user/YAML_INPUT_VARIABLES.md`，并运行 `tests/unit/test_config.py`。
- 新增 report 时写清日期与证据；被取代的材料移入 `docs/archive/reports/`，
  状态变化同步到 `FEATURE_STATUS_MATRIX.md`。
