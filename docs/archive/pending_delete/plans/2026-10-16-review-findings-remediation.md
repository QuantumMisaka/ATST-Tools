# 九项主线审查问题整改计划

**版本**: 2026-10-16
**日期**: 2026-10-16
**状态**: 已完成，结论已吸收到整改报告与长期参考；待删除复核
**责任人**: ATST-Tools maintainers

## 目标、来源与已有设计

用户于 2026-10-16 明确授权以 writing-plans + subagent-driven-development
确认并逐步修复主线审查九项问题。问题与 CPU 复现见
[主线审查 §6](../../../reports/MAIN_BRANCH_REVIEW_2026-10-15.md#6-当前实现独立系统复核2026-10-15)，
架构约束见 [架构复核](../../../reports/ARCHITECTURE_REVIEW_2026-10-16.md)。

已检查现有 SPEC/PLAN：

- 复用 [科学契约与交付设计](../../../superpowers/specs/2026-10-13-scientific-contracts-and-delivery-design.md)
  的 ASE-native、几何保真、公共执行、manifest 所有权、CLI/API 兼容边界。
  原双入口整改计划已完成并移入 pending-delete；不重开已完成任务。
- 复用 [runtime 接口设计](../../../superpowers/specs/2026-09-21-atst-runtime-interface-design.md) 的
  早期绑定、线程语义、已有 allocation 与有限批次边界；GPU 计划 P0–P5 已完成，
  剩余 P6 平台工作不归入本轮。
- workflow convergence 计划本仓部分已交付；保留完成/收敛/TS 验证分离。
  abacuslite force-read 计划剩余跨仓工作不覆盖本轮九项问题。
- 没有覆盖本轮九项缺陷的待执行计划，因此新增本文件。需求直接来自用户问题表，
  不另建重复 SPEC；实现本计划不需再次批准。

## 工作区、所有权与验证基线

- 原 checkout 的 main 和 27 个 dirty 文件保持原处；基线完整带入 native worktree
  `/home/james/.codex/worktrees/5738/atst-tools`，分支
  `fix/review-nine-findings-20261016`，BASE `5ad5de582313f853ae3d713df2fa5a6f2bd9cd5e`。
- 原 tracked binary diff、完整 untracked tar、状态及文件哈希位于
  `/tmp/atst-nine-fixes-20261016/`；迁入后逐项匹配。已有文档不能算本轮新增实现。
  各实现所有者的生产源码/测试均相对 BASE 干净；新增测试和后续 diff 纳入评审。
- 所有检查显式设置 `PYTHONPATH=<worktree>/src`，使用 `conda run --no-capture-output
  -n atst-dev python -m pytest ... -o addopts=-ra`；解释器已核验为 Python 3.10 环境，
  实际包路径是本 worktree，ASE 3.28.0，包版本 2.2.8。
- 复用同源码的 1284 passed/32 skipped、真实本地 MPI 27 passed 作为整改前基线。
  新行为测试先观察正确原因的 RED，再实施并记录 GREEN；既有复现脚本断言缺陷存在，
  不能直接当修复后的成功门禁，也不改写历史失败证据。
- 子代理只写分配的文件，不操作共享 index、不提交、不再派发子代理；文档与计划由
  控制器统一更新。修改共享 helper/API 前协调所有权；最终冻结 diff 后独立终审。
- 当前授权覆盖本地修复和验证；不发版、不改版本、不推送、不改父仓 gitlink。
  CPU/受控后端测试不冒充 GPU 科学或平台 P6 验收。

## 设计裁定

- R1：NEB 保留现有普通/TS 分段生成的图像数语义和既有轨迹，统一重启的带宽判断；
  不能通过改变历史生成结果来掩盖错误切带。6 帧与可被错误带宽整除的 30 帧均须验证。
- R2：IRC 默认/显式 null 使用有效输出默认值，显式路径保留；冻结已有计算结果时
  不主动触发额外后端计算，不共享活 calculator，也不将完成冒充收敛。
- R3：DMF 在构造优化路径前验证端点有序元素和原子数；不静默重排目标原子。
- R4：manifest 与消费者采用明确一致的路径基准，保留根目录旧行为、绝对路径和
  API/runner 交接；不能使用“哪个文件存在就选哪个”的搜索掩盖来源歧义。具体兼容
  编码由 Task 3 核验现有消费者后在此记录。
  实施裁定：writer 默认产物相对调用 cwd，嵌套或外置 manifest 保存可选 `path_base`
  （manifest 目录到产物根的相对位置）；根目录布局省略。缺字段旧 manifest 仍按
  manifest-parent 解析。API 保留现有 workdir-relative artifact 字符串，重写保留
  基准；直接 helper 可显式指定 `artifact_root`，不按绝对路径的内外位置猜语义。
- R5：bench 默认/显式/别名目录按实际解析结果判重；禁止 case_id 形成路径越界，
  安全的含空格 ID 应能正确报告。最终 worker 受管设备/线程/缓存事实必须符合预算，
  普通自定义 env 保留；不能静默扩大资源或改变普通单次运行的 OMP 优先级。
  实施裁定：批次内检查最终 runtime/YAML/args 的预算一致性；显式 calculator/MD
  omp 高于已记账 case threads 时可在 bench 拒绝，较小值保持。case args 不得
  覆盖 coordinator 拥有的 config/workdir/result-json 隔离路径。
  终审补充：在 load_cases 一次固定 config 绝对路径，预检与两种 worker 读取同一文件。
- R6：恒电势保持 facts-v1，恢复旧读者需要的 `profile_status` 兼容字段及其历史
  语义；新读者继续接受已产生的缺字段新格式。字段仅作兼容，不恢复已移除的非事实
  科学认证/状态推断。使用真实旧读取逻辑验证新→旧、旧→新与缺字段过渡数据。

## 实施与验收

### Task 1: NEB 重启身份（SCI-01，P1）

- [x] 统一普通与 TS guess 分段路径的实际带宽，重启完整恢复末带端点/顺序；保留
  不完整轨迹拒绝与既有 MPI 失败同步。
- [x] RED/GREEN 覆盖 n_images=3、TS guess→6 帧、30 帧历史、普通路径、损坏/截断
  轨迹；从 owning/public 执行边界验证初态不从 x=0 漂成 x=0.15。
- 所有者范围：`workflows/execution.py`、`utils/restart_helpers.py`，必要的
  `utils/idpp.py` 共享计数 helper，以及独立 `tests/unit/test_neb_restart_contract.py`。
  不占用现有大型 `test_workflows.py`，以便 Task 2 并行。

### Task 2: IRC 结果与 DMF 端点（SCI-02/03/04，P2）

- [x] IRC 缺省/null 双方向成功写出默认 normalized 轨迹；自定义路径仍有效。
- [x] descent 主/normalized/分方向轨迹都保留对应帧有效能量和力；保持几何、方向
  与阶段事实，确认结果保存不增加后端求值。
- [x] DMF 有序元素/原子数不匹配在优化初始化前明确拒绝；合法输入保留。
- 所有者范围：`workflows/irc.py`、`workflows/dmf.py`、独立
  `tests/unit/test_irc_result_contract.py` 和 `tests/unit/test_dmf_endpoint_contract.py`。
  优先复用已有 helper，修改共享 geometry/neb_endpoints 前交控制器协调。

### Task 3: manifest 路径契约（IO-01，P2）

- [x] 核验 writer、summary、API/runner、plots/runtime 引用等消费者后选择最小兼容
  修复并记录 R4 补充；子目录 manifest 正确找到工作目录中的 neb.traj。
- [x] 覆盖根目录旧布局、嵌套 manifest、绝对路径、不同 cwd 消费、缺失文件及 API
  round-trip；API enrich/rewrite 不二次重定位或改坏已有路径。
- 所有者范围：`utils/artifacts.py`、`utils/summary.py`、必要的 `api/services.py` /
  `api/models.py`，独立 `tests/unit/test_manifest_path_contract.py`。不修改工作流
  生产者文件；跨范围需求先交控制器。

### Task 4: bench 隔离与资源边界（BENCH-01/02/03，P2）

- [x] 默认配置目录/显式 workdir/相对别名/符号链接同址拒绝；合法独立目录可运行。
- [x] `case.env` 无法覆盖预算拥有的设备、线程、缓存及绑定事实；无 runtime 段
  也满足约束。核对 runtime 请求与批次预算，保留普通 env，独立与共享 worker 均覆盖。
- [x] `../escape` 等路径越界 ID 拒绝；`foo bar` 等安全 ID 成功运行并正确归属报告。
- 所有者范围：`bench/batch_runner.py`、`bench/batch_worker.py`、必要的 `bench/sweep.py`，
  独立 `tests/unit/test_bench_boundary_contract.py`，必要 owning bench tests。
  runtime 共用实现不由本包自行改动；范围外的资源优先级问题交控制器裁定。

### Task 5: 恒电势 facts-v1 互读（CP-01，P2）

- [x] 写出向旧 2.2.8 读取器兼容的事实；新读取器保留兼容旧含字段和过渡缺字段数据。
- [x] 测试旧读取逻辑实际接受新产物并能提取合法 restart electron count；几何/身份
  不符仍拒绝。不能仅断言字典含字段；夹具不依赖联网、外部 checkout 或测试时 git。
- 所有者范围：`calculators/constant_potential.py`、独立
  `tests/unit/test_cp_facts_compatibility.py` 及必要的小型历史 reader 夹具。

### Task 6: 集成、终审与交付

- [x] 控制器检查各包实际 diff、RED/GREEN 原始日志和剩余风险；关键科学/兼容边界
  交独立 reviewer，冻结所有写入后完成整体审查并修复有效问题。
- [x] 所有修复合并后运行完整默认 tests、opt-in 本地真实 MPI；补合适的 installed
  wheel API/runner 验证。无需重新执行未变 vendored drift 或 GPU 性能实验。
- [x] 更新用户参考/源码导览/HANDOVER 的受影响契约，功能矩阵记录未发布修复，
  主线审查保留历史并追加关闭证据，不把旧发布包写成已修复。
- [x] 新结果写入活跃整改报告并登记账本；计划完成后吸收结论、移入 pending-delete
  复核并更新归档 inventory。交付保留 worktree 与恢复资料，不自行合入 main。

## 执行记录

- 2026-10-16：完成已有 SPEC/PLAN 检索、隔离建立、27 项 dirty 内容迁入校验与环境
  路径确认。Task 1–5 可按上述独立写入范围并行，Task 6 等待全部实现与验收。
- Task 1 首轮：先观察 3 failed/2 passed（6 帧误拒、30 帧错切、截断误收），修复后
  5 passed，相关 owning checks 45 passed；原始日志 `/tmp/atst-nine-fixes-20261016/neb-*.log`。
  实际只改 execution.py 和新增 regression；独立审查待完成，暂不标记整体验收。
- Task 1/3/5 已完成独立审查，无实质问题。Task 3 相关 116 passed；Task 5 相关
  46 passed，历史 reader 夹具直接执行旧验证/重启逻辑，不依赖测试时 Git。
  Task 2 已实现并通过 39 项相关检查，独立审查中。Task 4 仍在完成入口预算检查。

- Task 2 独立审查通过。Task 4 终审发现相对 config 在子进程 cwd 下可指向另一文件，
  新回归 2 failed 后修复为一次冻结绝对路径，相关 96 passed，复审和整体源码终审通过。
- Task 6 最终源码冻结后：默认套件 1347 passed/32 skipped；本地真实 MPI 27 passed；
  clean-wheel API/CLI/runner 门禁通过。文档验收及完整日志/哈希记录于
  [整改验收](../../../reports/MAIN_REVIEW_FIXES_2026-10-16.md) 和随附 verification.json。
  原 main 的 27 项预存文件及状态保持不变；工作树保留，没有提交、合入或推送。
