# 文档系统交接与维护

**版本**: 2026-10-17
**日期**: 2026-10-17
**状态**: 维护
**责任人**: ATST-Tools maintainers

本文档是维护者日常 checklist。任何功能、YAML、CLI、backend、example、report 或
release 变更，都先从对应小节确认需要同步的文档。

## 1. 日常入口

- 用户入口：[README.md](../../README.md)、[docs/index.md](../index.md)、
  [USER_GUIDE_CN.md](../user/USER_GUIDE_CN.md)。
- 参数入口：[CONFIG_REFERENCE.md](../user/CONFIG_REFERENCE.md) 和
  [YAML_INPUT_VARIABLES.md](../user/YAML_INPUT_VARIABLES.md)。
- 源码入口：[源码导览](../../src/atst_tools/README.md)，用于理解模块目标、主要内容、
  调用关系与修改位置；新增目录或调整职责时同步更新。
- 开发入口：[YAML_INPUT_GOVERNANCE.md](YAML_INPUT_GOVERNANCE.md)、
  [DOCUMENTATION.md](DOCUMENTATION.md)、
  [EXAMPLE_VALIDATION_OPERATIONS.md](EXAMPLE_VALIDATION_OPERATIONS.md) 和
  [GOVERNANCE_AND_RELEASE_GATES.md](GOVERNANCE_AND_RELEASE_GATES.md)。
- 项目状态入口：[FEATURE_STATUS_MATRIX.md](../reports/FEATURE_STATUS_MATRIX.md)。

README 提供安装和首次使用路径；参数、后端与 API 细节由各自用户参考维护。
功能矩阵展示当前支持状态并链接验收报告，下游平台的发布与验收事实留在对应报告中。
GitHub 为开发规范源，Gitee 由维护者手动同步；PyPI 只分发 Python 包，示例从仓库取得。

## 1.1 新 checkout 与 CPU 开发准备

从 GitHub checkout 开发时，先创建 Python 3.10+ 环境，再以 editable 模式安装
项目及测试依赖。下面命令不需要 ABACUS、GPU、SAI 或 `mpi4py`，适合跑普通
CPU-only unit tests：

若当前环境已从其他 checkout 安装了 `atst_tools`，在 worktree 中运行测试前使用
`PYTHONPATH=src` 或 `pip install -e .`，并确认 `atst_tools.__file__` 指向当前
worktree。

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -c "import atst_tools; print(atst_tools.__file__)"
python -m pytest tests/unit -q
```

`atst_tools.__file__` 应指向当前 checkout 的 `src/atst_tools`。此环境仅能证明
不需计算后端的普通测试；ABACUS/DP、image-level MPI 和站点运行时验证另按
[EXAMPLE_VALIDATION_OPERATIONS.md](EXAMPLE_VALIDATION_OPERATIONS.md) 的证据边界执行。
不要把未执行的站点或计算验证写成通过。

CLI 和 Python API 共用内部执行实现 `atst_tools.workflows.execution`；服务层在
确定调用方的 MPI world 后调度执行 helper，CLI 层负责参数、模板与输出。独立
`validate_ts`/`atst vibration validate` 只读已有结果，不是 workflow 的自动后处理。

## 1.2 已交付计划的外部交接

2026-10-16 将两份本仓已完成实施计划归档；以下交接仍开放，归档不代表外部验收完成。

| 交接 | 责任方与关闭条件 | 本仓承接 |
| --- | --- | --- |
| abacuslite 实现维护 | calculator/parser 修复在独立 `QuantumMisaka/abacuslite` 仓实施；本仓仅维护 ATST adapter 与旧路径转发。 | 本页 backend checklist；历史快照账本只作审计。 |
| abacuslite 后续 SIF/harness | 父仓 ABACUS owner 负责交付与其验收；历史直接 SAI 力读取检验已通过，不替代后续平台交付。 | 后端快照和版本说明；不在本仓重开已完成实现。 |
| convergence 消费 | Paimon owner 完成 fixture 映射与联合验收；父仓 companion plan 为 `toolbox/ABACUS/docs/superpowers/plans/2026-09-17-matesbench-toolbox-hardening-plan.md`。 | [API 阶段事实](../user/PYTHON_API_REFERENCE.md)、`docs/reports/data/convergence_fixtures_20260917/`；真实运行时和 curated 产物更新须有相应证据。 |
| GPU P6 平台迁移 | ABACUS/MLIP-Agent owner 完成共享 image 参数/runner、guide、交接测试与 E2E/平台验收。 | [唯一活跃计划](../superpowers/plans/2026-09-20-atst-gpu-node-tuning-plan.md)，保留 DP 容器/FT²DP 条件项。 |

当前设计来源继续保留；是否归档取决于其决策是否已被长期参考吸收，不以代码已交付
为唯一依据。逐项处置见 [本地集成与收口](../reports/MAIN_REVIEW_FIXES_2026-10-16.md#6-本地合入与在途收口2026-10-16)。

## 2. 新增或修改 workflow

- 接入 `utils/config_schema.py` 的判别联合与类型表，更新 CLI 模板和
  `api/services.py` 的调度；按能力核对返回结构适配、合成产物、profiles/plots。
  新工作流应自己产出阶段事实和 manifest，不依靠 API 从轨迹推断收敛。
- 明确输入/重启身份、冻结结果和产物路径契约；复合工作流由外层拥有顶层 manifest，
  子任务交回阶段事实。同步验证 CLI 与 API 的科学执行一致性，保留各自既有输出契约。
- 更新 `README.md` 的支持列表和 quick start 说明。
- 更新 `docs/user/USER_GUIDE_CN.md`、`CLI_REFERENCE.md`、`CONFIG_REFERENCE.md`。
- 添加或更新 `examples/<case>/config*.yaml` 和 `examples/README.md`。
- 更新 `docs/reports/FEATURE_STATUS_MATRIX.md`。
- 添加或更新 workflow 测试、example 测试和必要验证报告。
- 运行期消息（print/warning/logging/异常）统一英文；改动收敛或诊断文案时同步更新
  锁定语义 token 的测试。
- 固定晶胞准备复用 `utils/geometry.py`：输入容差 `1e-5 Å、rtol=0`，容差内保真；
  NEB 所需适配只改内部运行副本的必要周期向量。必须在端点校验/计算前完成，
  cell/PBC 变化后丢弃旧 SPC、端点标记和恒电势事实。调用层传已有 world 协调 root
  warning 和已有失败同步，不在几何原语创建 communicator；长期说明见
  [CONFIG_REFERENCE.md](../user/CONFIG_REFERENCE.md#22-nudged-elastic-band-neb)。

## 3. 新增或修改 YAML 字段

- 修改 `src/atst_tools/utils/config_schema.py`。
- 重新生成 `docs/user/YAML_INPUT_VARIABLES.md`。
- 在 `docs/user/CONFIG_REFERENCE.md` 解释语义、默认值影响和常见配置组合。
- 更新示例 YAML。
- 运行 `tests/unit/test_config.py` 和相关 workflow/example tests。

## 4. 新增或修改 CLI

- 更新 `docs/user/CLI_REFERENCE.md`。
- 更新 `docs/skills/atst-cli/SKILL.md` 中的操作片段。
- 必要时更新 README 的轻量命令列表。
- 添加 CLI 测试或更新现有命令测试。

## 4.1 新增或修改稳定 Python API / process runner

- 保持 `atst_tools.api.__all__` 的既有稳定 root imports；runner 不得成为新的
  root import 或改变 `atst` CLI 行为。
- 更新 `docs/user/PYTHON_API_REFERENCE.md`、README、中文用户指南、release notes
  和功能矩阵，说明 JSON schema、退出码、manifest 权威性与 MPI 责任边界。
- 对 `python -m atst_tools.api.runner` 运行 API/document、clean-wheel 和实际 MPI
  门禁；runner 只能消费 API，绝不启动 Slurm、`mpirun` 或 `srun`。

## 5. 新增或修改 calculator backend

2.3.1 的发布与安装验证证据见
[发布说明](../releases/RELEASE_NOTES_2.3.1.md)；历史工作树处置见
[集成收尾记录](../reports/ABACUSLITE_230_CLOSEOUT_2026-10-19.md)。
旧 `integrate/abacuslite-atst229-20261017` 已被主线方案取代，其科学修复已经吸收；
后续工作从当前主线开展，不恢复该分支的 vendored 回退。

- 更新 README backend section 和项目边界说明。
- 更新 `docs/user/CONFIG_REFERENCE.md` 的 calculator 配置说明。
- 更新 `docs/user/USER_GUIDE_CN.md` 的运行时依赖和环境边界。
- 更新 feature/status reports 和 calculator factory tests。
- 若 backend 涉及 ABACUS/DP 环境，补充验证报告或在现有报告中记录边界。
- 2.3.1 声明 `abacuslite>=1.2.1,<2` 为运行时依赖，ABACUS 只使用
  独立安装的包。缺包或包内部导入损坏应失败，不回退。旧
  `atst_tools.external.ASE_interface.abacuslite` 顶层、`core`、
  `io.generalio`/`latestio`/`legacyio` 和 `utils.ksampling` 模块路径为薄转发层；
  Paimon 的直接 `generalio` 导入继续可用且与独立包模块同一对象。
- `.github/workflows/tests.yml` 的 `installed-abacuslite-consumer` job 构建 ATST
  wheel，正常 PR/push 使用已发布 abacuslite 1.2.1，schedule 跟踪兼容范围内
  最新 PyPI 版本，手动 `abacuslite_ref` 可从协调提交构建候选 wheel。只有一套
  external 安装环境；从 checkout 外调用
  `scripts/run_abacuslite_consumer.py --atst-root ROOT --backend external --report JSON`，
  再执行 `--verify-report JSON`。无 live vendored lane 或双报告比较。
- runner 核验已安装 ATST/abacuslite 的 site-packages 来源、安装文件指纹、
  ATST wheel 哈希、首次 runtime 消费者无旧路径导入、原始 SCF/MD fixture 的
  独立数值常量、完整离线测试 case ID/结果及 pytest 起止来源。测试因文档/example
  相对路径从仓根运行；`ATST_TEST_INSTALLED_WHEEL=1` 与共享子进程环境 helper
  保留安装来源，惰性 `sitecustomize` 防源码遮蔽。新增 Python 子进程测试复用
  该 helper。显式旧路径转发测试允许兼容导入，但不将其算成第二套实现。
- `ABACUSLITE_SNAPSHOT.md` 与 `PATCHES.md` 保留冻结历史；配对验收的 wheel
  hash、原始日志和 2026-10-17 结果保留在 Git 历史。旧 drift checker、vendored
  实现 CI 和快照 mutation 测试已退役；新的 parser 修复由独立 abacuslite 仓测试。
- 维护 `abacuslite` property-derived keywords 时，显式用户输入与自动派生值需按
  ABACUS 开关语义比较，例如 `True`、`1`、`"1"` 等价；但 `False`/`0` 与请求
  `forces` 或 `stress` 自动需要的 `"1"` 仍应作为冲突报错。

## 6. 新增或修改依赖

- 默认 `dependencies` 只放普通安装必须具备的轻量运行时依赖；重型后端、
  MPI、绘图和发布工具优先进入 `[project.optional-dependencies]`。
- Sella 是一等 workflow backend，默认依赖需保持 `sella>=2.5,<3`，除非有
  明确兼容性回退证据。
- DP 工作流使用 `dp` extra，MPI image-level parallelism 使用 `parallel`
  extra，NEB plotting helper 使用 `plot` extra。
- 基础安装测试必须保持不依赖 `mpi4py`；使用 `parallel` extra 的用户在 image-parallel
  验证前必须确认站点兼容的 MPI import 与 launcher。
- 修改依赖范围时，更新 README、`docs/user/USER_GUIDE_CN.md` 和本文件。
- 添加或更新 package metadata 测试，并运行构建检查：
  `pytest tests/unit/test_package_metadata.py -q`、`python -m build` 和
  `python -m twine check --strict dist/*`。
- 对包含稳定 Python API 的 release，运行
  `python scripts/verify_wheel_api.py`；它在临时目录构建 wheel、在临时 venv
  安装 wheel（默认会从配置的包索引解析依赖），并验证公共 root imports、H2/Au 的受控
  public API example、installed API runner JSON handoff 和显式验证 CLI/API 一致性，不在工作树保留 build artifact。若 release 环境已
  准备 `mpi4py` 和 `mpiexec`，额外运行
  `python scripts/verify_wheel_api.py --mpi-smoke`；该 two-rank smoke 限时，
  验证 image-parallel NEB/AutoNEB 的 rank-0 endpoint failure 以及 rank-local
  optimizer construction failure 会使所有 rank 以 typed API failure 退出，且在
  找不到 launcher 时明确跳过。

## 7. 新增或修改 example

- 更新 `examples/README.md` 的学习路径、目录说明和 chemical systems 表。
- 涉及示例的本地维护检查、SAI/Slurm 执行、curated outputs 或运行溯源时，先查
  `docs/developer/EXAMPLE_VALIDATION_OPERATIONS.md`；示例 README 不承载这些维护者操作。
- 更新 `examples/reference_results.json`，或明确该 example 暂无 reference 结果。
- 保证新增输入在 `inputs/` 或受控路径下，生成输出不进入 git。
- 运行 example 解析、dry-run 或 reference-result 测试。

## 7.1 新增或修改基准工具链（bench）

术语（2026-09-21 定）：**bench** = 测量工具链（`atst_tools.bench`）；**batch** = 一次运行
（一个 manifest 在一个已有 allocation 内跑一遍）；**harness** 是这批模块与产物的历史名，
仅保留在下面列出的冻结字符串里，不再用于新文字或新代码。

三件套：

| 工具 | 命令 | 作用 |
| --- | --- | --- |
| 批次执行器 | `python -m atst_tools.bench.batch_runner --manifest cases.json --out runs --devices 0,1 --slots 1 --cpu-budget 8 --threads 2` | 在既有分配内执行一张用例清单：槽位调度、超时回收、每 case 报告 + 批次汇总 |
| 变体扫描 | `python -m atst_tools.bench.sweep --manifest cases.json --out runs/sweep --devices 0,1 --slots 1,2,3 --repeats 3` | 变体 × 重复矩阵，交替顺序，出 `sweep_summary.json` |
| 测量记录 | `python -m atst_tools.bench.record --manifest cases.json --out runs/bench_record.json --run-dir runs --job-id … --partition … --qos …` | 自描述记录：修订/分支/dirty、清单与夹具哈希、结果树哈希、operator 字段 |

站点一键入口：`scripts/sai_runtime_bench.sbatch`（计时 sweep → 证据 pass → 记录）。
通道夹具（2026-09-22）：`record` 的 `--fixture`（可重复）对文件夹具逐字段不变；目录夹具按 `file_count`/`tree_sha256` 记录，新增可重复的
`--fixture-labeled LABEL=PATH`（与前者按命令行顺序共同写入 `inputs.fixtures`，元素多一个
`label`；PATH 是目录时另记 `file_count`/`tree_sha256`，缺件 warning 措辞不变）。
`scripts/sai_runtime_bench.sbatch` 因此按通道传参：`MODEL` → `dp_model=`、`ABACUS_INPUTS`
（文件或目录）→ `abacus_inputs=`，label 集合写入 `--note`，ABACUS 通道不再记成 DP 模型哈希。
清单模板：`examples/runtime_batch_cases.example.json`。清单字段：`case_id`、`config`、
`slots`、`threads`、`timeout_s`、`workdir`、`launcher`、`ranks`、`env`；**每个 case 必须有
自己的目录**（未给 `workdir` 时按配置目录）。2026-10-16 未发布修复以解析后的
实际路径判重，包括显式/默认目录和符号链接别名；case_id 必须是安全的单路径段，
允许内部空格，进度协议须完整传递 ID。env 可传普通自定义值，受管设备、线程、缓存
和绑定事实由批次控制；YAML/runtime、显式 omp 和 runner 参数不能突破该预算或
改写 coordinator 管理的输入/输出位置。相对 config 在加载清单时按批次启动 cwd
解析一次，预检与两种 worker 使用同一绝对路径，避免 chdir 后读到同名异文件。
具体证据及兼容说明见 [九项问题整改验收](../reports/MAIN_REVIEW_FIXES_2026-10-16.md)。

线程记账为 `case.threads × ranks`；`ranks` 优先使用 case 显式声明，否则尝试从
launcher 的 rank 参数识别，未识别时为 1。任意 calculator command 启动的后端
子进程数不由 bench 自动审计；使用嵌套 ABACUS MPI 时需显式计入其资源需求。
这是调度预算，不是操作系统层的 CPU 硬隔离。

产物（新运行）：`<out>/<case_id>/case_report.json`（每 case）、`batch_summary.json`（批次）、
`worker.out` / `worker.err`（每 case 报告目录内）；用例配置带 `runtime.telemetry` 时另有
`runtime_evidence.json`；批量侧 `bench_record.json`。

**共享 worker 模式（`--share-worker`，2026-09-21 新增）**：整批 case 在同一个 worker 进程里顺序执行，
模型加载与首次调用预热只付一次（维护者实测：2 例 CCQN 批 68.2 s → 34.2 s，第 2 例只花 3.1 s）。语义与约束：
要求 `--slots 1`；拒绝 per-case `launcher`、非一致的 per-case `env`、超预算线程与空 `--devices`；
共享 worker 以「已绑定进程」身份启动（用 `runtime.launch.build_child_environment` 发布
`ATST_RUNTIME_BOUND`/`ATST_INHERITED_DEVICES`/`ATST_EFFECTIVE_DEVICES`/`ATST_BINDING`），所以带 `runtime:`
段的配置可以运行，**与绑定矛盾的请求仍逐 case fail-closed**（例如只持有 0 号卡时请求 `devices: [1]`）；
单进程没有 rank 轮转，配置里的 `round_robin` 不产生旋转；`--timeout` 与 `--stop-on-failure` 仍生效
（超时杀 worker → 该 case `timeout`、其余 `failed`/`shared worker exited`）。新增产物：
`<out>/batch_jobs.json`、`<out>/worker.out|err`；每 case 报告多一个 `"shared_worker": true` 键。

命名沿革（先读）：2026-09-21 之前上述产物叫 `harness_case.json` / `harness_summary.json` /
`harness_worker.{out,err}`；**归档切片 `docs/reports/data/**` 仍带旧名，照原样保留**，
`bench_record` 读取时同时接受新旧汇总名（新名优先）。两条历史拼写**冻结不动**：schema 字符串
`atst-bench-harness-v1`（版本化文档标识，不是模块名）与证据取值 `ATST_THREADS_SOURCE=harness`
（归档 sidecar 的 `threads_source` 就是它）。

边界：它不是生产队列/调度器（不申请资源、不跨作业排队），也不新增 `atst batch` 公共
接口；科学逻辑在 `workflows/`。回归：`tests/unit/test_bench_batch_runner.py`、
`tests/integration/test_bench_batch_runner_mpi.py`。

## 8. 新增 report 或移动旧 report

- 新增 report 时写清日期、环境和证据范围。
- 只有当 report 是核心入口或当前重点证据时，才从 `docs/index.md` 链接。
- 被取代的 report 先吸收结论，再移到 `docs/archive/reports/`。
- 新增 spec/plan 放入 `docs/superpowers/specs/` 或 `docs/superpowers/plans/`；
  完成后吸收结论并移入 `docs/archive/`。

## 9. 准备 release

- 包版本以 `pyproject.toml` 的 `[project].version` 为唯一来源；
  `atst_tools.package_version()` 在源码树读取该字段，安装后读取发行包元数据。
- 更新 `docs/releases/RELEASE_NOTES_<version>.md`。
- 更新 README badge、版本说明和 release scope。
- 运行离线门禁：`python scripts/check_release_readiness.py --tag v<version>`。
- 按 [GOVERNANCE_AND_RELEASE_GATES.md](GOVERNANCE_AND_RELEASE_GATES.md) 完成本地/CI
  机械门禁；当前模型家族的独立审查是正常终审。触及 AGENTS、SKILL 或治理契约时，
  维护者可按需要求跨家族审查作为 owner-judged advisory evidence；可移植的冻结
  diff 审查证据与 launcher 产物同样有效，缺少跨家族审查不阻断交付。
- GitHub 分支保护、`pypi` Environment 审批/标签限制、PyPI Trusted Publisher、
  push/tag/release、Gitee 手动同步及发布后渲染检查均由管理员/维护者负责，
  不是仓库文件修改的副作用。
- 若 release 包含示例或 curated-output 变更，按
  `docs/developer/EXAMPLE_VALIDATION_OPERATIONS.md` 复核运行证据与溯源。
- 更新 `FEATURE_STATUS_MATRIX.md`。
- 运行文档链接、格式和相关测试。

## 10. 最小检查

```bash
git diff --check -- README.md docs examples/README.md
rg -n "^<<<<<<<|^=======|^>>>>>>>" README.md docs examples/README.md
python scripts/check_docs_governance.py
```

若修改 YAML schema：

```bash
conda run -n atst-dev python -m atst_tools.utils.config_docs --output docs/user/YAML_INPUT_VARIABLES.md
conda run -n atst-dev pytest tests/unit/test_config.py -q
```

改动运行期消息（print/warning/logging/异常/CLI 文案）时：

```bash
conda run -n atst-dev pytest tests/unit/test_program_output_language.py -q
```
