# 九项主线审查问题整改与验收

**版本**: 2026-10-16（2.2.8 后未发布修复）
**日期**: 2026-10-16
**状态**: 九项修复已合入本地 main；未推送、未发布
**责任人**: ATST-Tools maintainers

生命周期：validation；报告级别：L2。需求来自用户逐项确认的九项缺陷；历史现象、
严重性及负向复现保留在 [主线审查 §6](MAIN_BRANCH_REVIEW_2026-10-15.md#6-当前实现独立系统复核2026-10-15)。
本报告记录修复后行为与验收，不改写已发布 2.2.8 的内容或历史测试结果。

## 1. 依据与实施边界

先检查已有 SPEC/PLAN，复用科学契约、runtime 设计及完成/收敛分离规则。原双入口
整改已完成，现有 GPU 计划剩余平台事项不覆盖本轮；因此新增九项整改计划，采用
writing-plans + subagent-driven-development，按 NEB、IRC/DMF、manifest、bench、
恒电势五个独立写入范围实施，由控制器统一集成和文档维护。

执行分支为 `fix/review-nine-findings-20261016`，worktree 为
`/home/james/.codex/worktrees/5738/atst-tools`，BASE 为
`5ad5de582313f853ae3d713df2fa5a6f2bd9cd5e`。原 main 的 27 项未提交文档/证据完整
带入并逐项校验；修复验收阶段原 checkout 保持不动，带入内容不计作本轮源码实现。
随后获授权执行本地合入，见 §6。
恢复快照位于 `/tmp/atst-nine-fixes-20261016/`。

环境使用 `atst-dev`，执行检查明确 `PYTHONPATH=<worktree>/src`，避免环境中的
editable 安装指向其他 checkout。本次不改变版本、不发布、不推送、不更新父仓
gitlink，也没有运行新的 ABACUS/DP GPU 或平台 P6 验收。

## 2. 修复与验收对应

| 原问题 | 修复后行为 | 主要回归 |
| --- | --- | --- |
| SCI-01 / P1：NEB 重启错带 | TS guess 分段路径按 `n_images + 3` 恢复，普通路径保持 `+2`；不改变已生成路径。 | 实际 ASE 生成/执行；6 帧与 30 帧恢复完整带及 x=0 初态，5/29 帧截断拒绝。 |
| SCI-02 / P2：IRC 默认输出 | 缺省/null normalized 输出使用有效默认文件名，自定义路径保持。 | 双方向完成后的读回及 manifest；受控 Sella optimizer。 |
| SCI-03 / P2：descent IRC 结果 | 主/normalized/返回帧保留独立冻结能量和力；有效 live cache 或匹配的已存分方向帧提供结果。 | 真实 ASE FIRE/受控 calculator；逐帧几何/能量/力一致，保存不增加求值；缓存回退拒绝错几何。 |
| SCI-04 / P2：DMF 元素顺序 | 优化器初始化前拒绝原子数或有序元素不符，不重排终态。 | HO→OH、原子数不同均拒绝；合法 HO→HO 继续到优化器加载边界。 |
| IO-01 / P2：manifest 路径 | writer 记录必要的 `path_base`，summary 按明确基准读取；API enrichment 保留路径语义。 | 根目录、嵌套、外置绝对 manifest、绝对 artifact、异 cwd、缺失产物、API round-trip。 |
| BENCH-01 / P2：工作目录 | 按解析后的实际目录判重，默认/显式路径、相对别名和符号链接同址均拒绝。 | 独立目录仍可运行；冲突在启动前发现。 |
| BENCH-02 / P2：环境越过预算 | 受管环境以批次分配为准；核验 YAML/runtime、显式 omp 和 runner 参数，普通自定义 env 保留。 | 无 runtime 配置的独立/共享子进程实际环境；额外线程表示和受管参数路径。 |
| BENCH-03 / P2：case ID | 拒绝路径分隔、`.`/`..` 和控制字符，安全的含空格 ID 完整传递。 | `../escape` 拒绝；`foo bar` 成功且报告归属正确。 |
| CP-01 / P2：facts-v1 互读 | 新写出补回旧版要求的历史 `profile_status`；新读者继续接受含字段旧数据及缺字段过渡数据。 | 实际 2.2.8 reader/publisher/restart 夹具，双向读取与电子数恢复；几何/身份不符仍拒绝。 |

回归集中于 `tests/unit/test_neb_restart_contract.py`、`test_irc_result_contract.py`、
`test_dmf_endpoint_contract.py`、`test_manifest_path_contract.py`、
`test_bench_boundary_contract.py`、`test_cp_facts_compatibility.py`。
旧 reader 夹具固定历史逻辑，测试时无需 Git、网络或外部 checkout。

## 3. 兼容性与迁移

### 科学行为

普通和 TS 分段 NEB 的生成图像数保持不变；修复重启的错误切带。恢复旧轨迹时保留
原 `make` 设置，不把已被错误切带并继续计算的轨迹当作自动修复的数据；该类历史
计算需核对端点和原始完整带后重新启动。

IRC 保存已有结果，不额外调用后端；缺证据时不虚构能量、力或收敛。DMF 新校验
拒绝先前会被错误解释的端点，用户应在输入准备时建立正确的原子对应。

### Manifest

artifact 和 plot 字符串仍以工作目录为基准；可选顶层 `path_base` 是 manifest
目录到该工作目录的相对路径。根目录布局省略 `.`，绝对 artifact 路径不变。
API JSON 继续提供绝对 workdir/manifest，保持 artifact 字符串及既有字段。

缺 `path_base` 的旧 manifest 继续按 manifest 目录解析。旧嵌套文件如果记录的
实际是 workdir-relative 路径，必须由知道原始工作目录的使用者补明基准；程序不
搜索同名文件。旧消费者需理解新增基准后再消费嵌套布局，默认根目录布局兼容。
直接使用内部 writer 时，`artifact_root` 可显式指定调用方的路径基准。

### Bench

批次拥有设备、线程、可写缓存和绑定事实，普通 env 仍传递；受管 env 不再由
case 覆盖。超出已记账线程的 YAML/runtime 或显式 omp 请求在启动前拒绝；runner
参数不能改写协调者管理的 config/workdir/result-json 或资源控制。
相对 config 在清单加载时按批次启动 cwd 固定为绝对路径，预检、独立 worker 和
共享 job 使用同一文件；切换 workdir 不再改变配置身份。合法含空格 ID 完整传递。

线程记账依旧使用 `case.threads × ranks`，rank 取显式声明或已识别的 launcher
参数。任意 calculator command 的嵌套进程不自动进入审计；后端 MPI 的资源需求
须在案例中显式计入。bench 是有限实验调度，不能据此声称实现操作系统硬隔离。

### 恒电势

facts-v1 写出恢复两种历史 label：`reference_only_unvalidated` 和
`compensated_gate_pending_scientific_acceptance`。它们只解决旧 reader 的必填字段
兼容，不作为科学验收状态，也不恢复 calculator 的历史状态推断。

过渡期缺字段 frames 在当前 reader 中继续可读；原 2.2.8 reader 不会因新代码
修复而改变。向旧版交接这类 frame 时需经过当前实现的有效性校验并重新发布 facts，
而非跳过几何、Hamiltonian identity 或 residual 校验。

## 4. 验证与审查记录

五个工作包均先观察缺陷回归 RED，再完成修复与 GREEN；控制器核对实际差异和日志，
独立 reviewer 检查科学、兼容和资源边界。bench 初审发现预检与隔离 worker 可能读取
不同目录下的同名配置，补充独立/共享回归后复审通过；整体冻结源码终审未发现剩余
实质问题。最终集成检查在该修复后重新执行；文档终审核对修复状态、迁移说明、
账本与冻结输入/日志哈希后通过。

| 验证范围 | 结果 | 原始证据 |
| --- | --- | --- |
| 完整默认套件 | 1347 passed、32 skipped、3 warnings | [完整日志](data/review_fixes_20261016/full-tests-final.log) |
| opt-in 本地真实 MPI | 27 passed | [MPI 日志](data/review_fixes_20261016/mpi-tests-final.log) |
| 临时构建 wheel、干净安装 API/CLI/runner | passed | [安装检查](data/review_fixes_20261016/wheel-api-final.log) |
| 文档治理与文档回归 | passed；26 passed | [治理检查](data/review_fixes_20261016/governance-final.log)、[文档测试](data/review_fixes_20261016/docs-tests-final.log) |
| 相关回归范围 | NEB 45、IRC/DMF 39、manifest 116、bench 96、恒电势 46 passed；范围有交集，不相加 | [命令、RED/GREEN 与哈希索引](data/review_fixes_20261016/verification.json) |

默认跳过项中的 27 项 MPI 已由单独 opt-in 执行；其余 5 项依赖分发制品、真实
ABACUS 示例、外部 Toolbox 或 DP 模型，跳过原因保留在日志。三个 warning 来自
既有 IDPP 未完全收敛负例和 ASE logfile flush 弃用提示。wheel 脚本在临时目录
构建与安装，未保留发行制品，也不替代完整发行检查。

验证索引记录冻结源码/测试/构建输入 SHA-256、命令及退出码；原 main 的 27 项
预存改动和状态逐项对照未变。文档治理、26 项文档测试及差异/冲突检查通过。完成计划
移入 pending-delete，结论由本报告和长期参考承接。

## 5. 尚未宣称的结果

本轮 CPU 与本地 MPI 证据只验证科学工作流契约、进程/文件交接及资源输入边界。
受控 Sella IRC 不能代替真实 Sella 路径科学验收，DMF 输入检查不等于 IPOPT 求解
或 TS 验证，恒电势互读不等于真实电化学有效性；平台 P6 和发行仍各自独立。


## 6. 本地合入与在途收口（2026-10-16）

用户随后明确授权核验、合入与归档。九项修复、回归和前轮架构/源码导览文档以
`9e5f3b8` 提交，并从 `5ad5de5` 快进进入本地 main；没有重写已有历史。冻结的
198 项源码/测试/构建输入与验收记录一致，复用完整套件、MPI 和 wheel 证据；合入
后在 main 实跑六个新增回归文件，**63 passed**。

原 main 的 27 项未提交内容先保存为 stash
`09bfc9a6eae61dde6dcba0ec3c293a970a235cbd`，逐项核验文件内容，并建立独立恢复引用
`refs/archive/closeout-20261016/pre-integration-main`。其内容已由合入提交吸收，
不重复应用该 stash。父仓 HEAD 与 index gitlink（`45d8c78`）保持不变；未推送或发布。

| 对象 | 处置 | 理由与承接 |
| --- | --- | --- |
| abacuslite 08-04 计划 | 本仓收口、归档 | Tasks 1–8 已交付；父仓 SIF/harness 与延后的上游同步由 HANDOVER §1.2 保留。 |
| convergence 09-17 计划 | 本仓收口、归档 | A0–A6 随 2.2.6 交付；API 参考/fixtures 承接事实，Paimon 联合消费仍由其 owner 验收。 |
| prepare/stress 08-12 设计 | 归档 | CLI/API/配置参考、2.2.3 发布记录及测试已吸收。 |
| 用户分发/MPI 安装 09-04 设计 | 归档 | 用户/开发者长期入口、2.2.4 发布验收及测试已吸收。 |
| abacuslite 08-04 设计审查 | 归档 | 结论已纳入仍活跃的跨仓设计，保留历史原文。 |
| GPU 09-20 计划 | 继续活跃 | P6 共享 image 参数/runner、owner guide、交接测试和平台验收未完成；DP 容器/FT²DP 条件项仍开放，TF 专项仍 not-planned。 |
| 其他现行 SPEC / GPU P0 review | 保留 | 仍指导科学、runtime、后端、治理与跨仓边界；GPU/runtime/runner 的过时状态已更新，历史基线保留。 |

工作区和恢复资料采用保留式归档：

- 本轮 Codex 管理的工作树 `/home/james/.codex/worktrees/5738/atst-tools` 保留路径、
  分支和注册；完整文件快照（含 ignored 产物）、独立 Git bundle 与恢复说明归档至
  `/home/james/scratch/atst-tools-closeout-20261016-124622`。
  当前工具只有创建/列出能力，没有解除挂载能力；不移动当前受管路径以免破坏会话
  附件。该分支已吸收，不再有独立未合入代码。
- 旧 `fix/ase-ts-dual-interface-20261013` 的 `3163c06` 已是 main 祖先；其工作树
  `/home/james/scratch/atst-tools-closeout-20261015-153534/atst-tools-worktree` 已在归档
  位置，tracked/untracked 干净，历史 ignored 产物保留，无需重复搬移或删除分支。

归档和文档状态经独立审查；最终文档治理、检查日志与恢复说明随收口保存。所有
移动均保留原文，不把跨仓待办或平台 P6 写成已验收。
