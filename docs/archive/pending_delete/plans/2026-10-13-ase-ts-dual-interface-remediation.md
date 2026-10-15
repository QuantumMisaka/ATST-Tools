# ASE 过渡态计算与 CLI/Agent 双入口整改计划

**版本**: 2026-10-13
**日期**: 2026-10-13
**状态**: 本地实施与验收完成，未发布；待删除复核
**责任人**: ATST-Tools maintainers

## 目标、依据与授权

基于 ASE 调用 ABACUS/DP 完成过渡态计算，CLI/YAML 与 Agent Python Interface
共同消费科学执行和结果层。修复 [审查报告 §8](../../../reports/ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md#8-面向-ase-过渡态计算及双入口的复核与优化)
列出的缺陷与交接债务；遵循 [设计 rev.8](../../../superpowers/specs/2026-10-13-scientific-contracts-and-delivery-design.md)。
用户明确要求 writing-plans + subagent-driven-development 并授权逐步修复，计划固化不
重新索取实施许可。当前不执行发布、推送、外部消息或 GPU 作业；真实后端科学验证
依项目 SAI/Slurm 约定，在有可用站点条件时另补，不以 CPU/受控 calculator 冒充。

用户进一步确认：TS 验证是独立功能，由用户/Agent 显式启动；计算交付优化结果及
已有收敛事实，不能自动启动振动/IRC或声称已完成 TS 验证。既有 D2S 振动默认 false，
显式启用的历史路径作为用户 opt-in 保持兼容。独立验证可消费现有振动结果；未知
驻点性保持 unknown，不静默启动 calculator。回调异常采用所有 rank 同步失败退出；
TS JSON 保留既有 pass/fail，增明确验证范围、驻点性与完整检查状态。

## 隔离、基线与分工

- 原 checkout：`/home/james/work/sidereus/app-tools/toolbox/ABACUS/deps/atst-tools`，main
  上已有已验收晶胞/文档 dirty 改动；保持原处不动。
- 执行 worktree：`/home/james/.codex/worktrees/3fbf/atst-tools`；分支
  `fix/ase-ts-dual-interface-20261013`，HEAD `45d8c78de8a605c127cd4c36742d0c9525a0a6b6`。
- 基线已完整迁入；tracked binary patch、untracked tar、状态和 pre-dispatch 内容快照
  位于 `/tmp/atst-remediation-20261013/`。已有 dirty 不能算子代理新增成果。
- `atst-dev` editable 指向另一 checkout。所有执行检查须明确
  `PYTHONPATH=/home/james/.codex/worktrees/3fbf/atst-tools/src`，并从 worktree 根调用
  `conda run --no-capture-output -n atst-dev python -m pytest ... -o addopts=-ra`。
- 子代理独立文件写入；不提交、不操作共享 index、不派发下级 reviewer。控制器集成
  与登记长期文档/计划，独立当前家族 reviewer 在重要工作包和最终集成边界审查。
- 原晶胞工作包默认全量 1230 passed/26 skipped、真实 MPI23 passed；新隔离源码迁入
  后先验证实际 module 路径，未改领域复用证据，实际疑点补针对性检查。

## 实施与验收

### Task 1: 完整周期 IDPP 的一般 MIC

- [x] 修复合法 unwrapped 多胞与高度斜胞 pair translation，保留 NEB-like IDPP
  更新、固定 pair translation、ASE 约束和输入保真。优先复用 ASE find_mic；保持
  确定的等距映像处理，不用与旧错误扫描逐位一致替代科学验收。
- [x] 回归 x=1/28 Å、cell10、三个相同 frame 的 objective/force 为零，以及整数
  晶格平移不变性、斜胞、partial PBC、mic=False、约束与既有示例。
- 预计所有者：`utils/idpp.py`、`tests/unit/test_idpp.py`、必要 `test_ccqn.py`。
  changed behavior 用 RED/GREEN；记录既有路径变化及原因，不能改历史输出伪造复现。

### Task 2: 分发许可与 Relax 输入回退

- [x] 显式打包第三方 LICENSE/UPSTREAM 与必要 abacuslite provenance；补 PyDMF
  公开来源 URL（仅经本地上游 git/原始材料证实者），不把私有临时目录当公开来源。
- [x] 修复 Relax 回退找到 init.traj 后未读取 atoms；正常/回退/缺失/restart 通过
  CLI/API或 owning workflow 边界测试。
- [x] 隔离 PEP517 sdist→wheel 构建，检查实际许可正文与 provenance 成员；twine
  与 clean-wheel API gate 各按其实际范围报告。
- 预计所有者：`pyproject.toml`/打包清单、`external/pydmf/UPSTREAM.md`、
  `workflows/relax.py`、专用 distribution/relax 测试。与 Task1/3/4 可并行；共享
  版本/index 操作由控制器串行，当前不 bump 版本或发布。

### Task 3: 进度失败的 MPI 与 Agent 交接

- [x] root 进度 stream/callback 失败同步到所有 rank，转 typed workflow failure；
  包括开始、image/profiles/plots 等尾部 root 操作，失败不让 peer 进入另一 collective。
- [x] 用户选定同步失败退出；保留进度开关及事件格式，不增加 advisory 默认或新 knob。
- [x] 起始与计算后回调失败、stream failure、成功事件路径、2-rank bounded exit、
  结构化 API/runner 错误（可调用边界）回归，清理测试进程组。
- 预计所有者：`api/services.py`、`test_api.py`、相关 MPI integration 测试。
  Task5 修改此文件须等待 Task3 写入和验收完成。

### Task 4: 优化事实与显式独立 TS 验证

- [x] Dimer 保存真实 optimizer return/nsteps，输出 StageRecord/manifest；D2S
  复用该事实，保持执行完成与收敛分离，不自动把未收敛改成异常。
- [x] 独立验证扩展现有 summary 的 scope、stationarity 与明确综合检查状态，保留
  旧字段和 pass/fail 兼容语义。缺少有效力/结构/振动范围事实时保持 unknown；不能
  从孤立虚频数量断言整个结构 TS 已确认。受约束/子集范围必须明确，不能默认为全体系。
- [x] 计算 workflow 默认不启动验证；已有显式 vibration/opt-in D2S 仍可交付筛查。
  提供可由 CLI 和 Python 调用的显式校验入口，读取已有输出，不暗中启动新 calculator。
- [x] 相邻 Dimer/D2S/vibration、unknown/malformed/finite data、子集/约束范围、默认
  无验证及显式验证两入口检查。新增结果为增量字段，不改旧退出码/公共 imports。
- 预计所有者：`mep/dimer.py`、`workflows/d2s.py`/`vibration.py`、
  `utils/ts_validation.py`、独立 validation API helper、专用 tests 和 test_workflows.py。
  CLI/API root export 与 docs/generated YAML若需要由控制器协调，避免 Task3 同文件写入。

### Task 5: CLI/Agent 共用执行层

- [x] 把 scripts.main 执行 helper/runner 下沉到职责清晰的内部模块；API services
  直接 dispatch 共用层，preflight 直接 utils.abacus_io，CLI 留模板/参数/展示。
- [x] 保留稳定 root API、YAML、CLI、import-light 绑定、外部 world、backend command
  行为、当前 CLI/API manifest/preflight 的历史差异。不建设插件或调度框架。
- [x] baseline/candidate owning checks、相同 fixture 的 CLI/API结果、已有 monkeypatch
  保护转向真正执行层；wheel安裝与真实MPI再验，确保API import不通过CLI执行层。
- 依赖 Task3 与 Task4 文件交接完成；Task1可以先完成。预计 files：新执行模块、
  `scripts/main.py`、`api/services.py` 和相邻 unit/integration tests。

### Task 6: 后端身份、双入口与兼容验证

- [x] 共用 manifest 明确实际 backend/source、可用版本/snapshot补丁事实、DP model
  与 head 或 supplied calculator 类型；复用已有实例/配置事实，不二次加载模型或
  强制逐次大文件hash。未知事实明确unknown，计算结果身份与software环境分工。
- [x] CLI/Python 同配置对比实际 cell/PBC、冻结能量/力、阶段事实、异常和 artifacts。
  结构 snapshot 与冻结结果消费说明保持清晰；process runner 每个case独立workdir。
- [x] 普通 CPU CI 保持最低 Python 3.10，增加较新受支持 Python/ASE环境；干净 wheel
  用公共入口验证，不把site环境/大模型作为贡献者前提。
- 依赖 Task5；预计 metadata/backend helpers、services与共享 manifest、专用 paired
  tests、CI与verify_wheel_api。程序性输出英文，测试按观察行为而非实现快照设计。

### Task 7: 文档收敛、集成与最终独立审查

- [x] 更新用户/API/CLI参考和示例：计算与独立验证、同步callback失败、Dimer收敛、
  运行副本/几何边界与backend identity；新增或修改 schema则生成YAML变量表并测config。
- [x] 补贡献者 editable 安装/普通CPU测试入口；发布事实引用唯一release证据，校正
  本地确证tag/SHA与应用不启动外层MPI但backend执行配置命令的表述。
- [x] 更新报告、FEATURE_STATUS_MATRIX、DOCUMENTATION_STATUS_REPORT，与spec当前
  状态；新执行 plan登记，完成结论吸收后移入pending_delete并维护治理表。
- [x] 集成全量 tests、实际MPI launcher、vendored tests/drift、隔离build/产物内容、
  clean-wheel API、docs script/26tests/diff/conflict检查；任何失败诊断后修复，未执行
  的GPU/站点/发布范围准确列出。
- [x] 独立 reviewer 对完整 attributable delta + 既有基线复核，修复重要问题后交付
  retained worktree。原 main dirty 内容、发布版本与remote保持原状。

## 进度与裁决

- 2026-10-13：计划建立；隔离worktree与dirty基线迁入已完成，先启动Tasks1–4，
  待交接后执行Task5–7。全部工作为本地未发布开发。
- Ruling P1：保留现有原始工作区，使用native managed worktree +完整patch/tar基线
  迁入 — 避免main上并行编辑和丢失已验收改动 — 新分支交付后整合应按原main最新
  状态处理，不覆盖用户之后的独立编辑。
- Ruling P2（用户确认）：进度失败所有rank同步中止；TS验证保留旧筛查状态并增
  scope/驻点性/综合状态，且只有显式验证请求才执行 — 保持科学计算与验证责任清楚，
  Agent不会挂住或把筛查pass当整个TS已验证 — 新事实与既有结果保持兼容。

### 2026-10-15 恢复执行

用户确认按已有顺序继续，以 writing-plans + subagent-driven-development 执行。
复用本 PLAN 与 SPEC rev.8；本轮按 L3 协调关联工作包，不新建重复计划。
原 main 工作区保持不动；既有 worktree 的完整 tracked/untracked 恢复基线位于
`/tmp/atst-remediation-20261015/resume-baseline`，差异和状态在同级目录。
Tasks 1–4 既有源码与原始日志已找到，当前进行冻结快照独立复核，不因已有实现直接勾选。
Task 5 公共执行层、Task 4 显式验证入口补齐及 Task 6 身份记录按不重叠文件派发。
本轮日志集中于 `/tmp/atst-remediation-20261015/`，最终必要证据将吸收入长期报告。
平台 P6 是独立下游工作：本计划交付本仓接口与证据，不把平台发布或 GPU 验收默认为完成。

- 2026-10-15 冻结快照独立审查：Tasks 1–3 可接受，复用原始 RED/GREEN、
  分发构建与双 rank MPI 日志。Task 4 四项 P2（缺失范围误报全体系、默认固定原子
  范围、验证覆盖结构输入、CLI 缺口）已补实现和回归；专项 19 项通过，控制器
  科学相关组合 105 项通过。全量集成与最终独立审查仍待完成。
- CI 兼容矩阵配置为 Python 3.10/3.12 × ASE 3.28.0/当前可解析版本；本地配置
  检查通过不代表远端四个环境已经执行。当前本地测试环境为 Python 3.10.19、ASE 3.28.0。

- 2026-10-15 最终收口：全量 1284 passed/32 skipped，真实 MPI 27 passed，vendored
  28 项/2 skipped，drift 0，隔离 sdist→wheel、实际 notice 两项、Twine strict、
  clean-wheel API/CLI/runner/显式验证及 MPI smoke 均通过。独立终审的 Dimer import
  与显式固定原子振动兼容发现已修复并复核通过。结果吸收至系统审查 §9 和长期文档。
  保留 worktree，不提交、不推送、不发布；平台 P6 和真实 GPU 科学验收独立跟踪。
