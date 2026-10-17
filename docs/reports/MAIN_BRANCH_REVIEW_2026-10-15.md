# main 分支更新与提交批次审查

**版本**: 2026-10-15
**日期**: 2026-10-15
**状态**: 维护
**责任人**: ATST-Tools maintainers

> §1–§5 保留首次历史整理审查；§6 是对 `5ad5de5` 的独立系统复核。
> 提交清单完整性得到确认，另确认 1 项 P1、8 项 P2（含原有恒电势风险）。§1–6 保留
> 修复前审查事实；后续源码整改及验收见 §7 与 [整改报告](MAIN_REVIEW_FIXES_2026-10-16.md)。

## 1. 范围与方法

冻结对象为 main `065da5a`（本报告提交自身不在统计内），共 **734 个可达提交**，
包括侧枝提交及 merge commit；不只统计 first-parent。远端主线基线为 `45d8c78`，
最新发行 tag 为 `v2.2.8` → `b6a68f9`。当前包版本仍为 2.2.8，不能据此把
main 等同于已发布 wheel。

本轮完成提交图、发布说明、分支恢复引用、当前相关实现与已有验证证据的交叉核对。
历史批次属于结构/交付审查，不代表重新逐行审查 734 个提交或重跑所有历史科学计算。
B11 复用已完成的整体代码独立审查、修复复审及合入后 1284 passed/32 skipped；
本轮另对集成与发布边界进行独立复核。审查发现与待办见 §4。

提交不携带原始分支名。merge 消息可提供分支线索；快进、rebase 和手动吸收只能
根据保存的 tip、祖先关系及补丁语义解释，不能凭 `branch --contains` 推断唯一作者分支。
以下按互斥版本区间分批，来源表另列，因此不会把共享祖先重复计为某分支贡献。

## 2. 不重叠的提交批次

| 批次 | Git 范围（左端不含、右端含） | 提交数 | 内容 |
| --- | --- | ---: | --- |
| B01 | `v1.5.1` | 291 | 早期脚本与 legacy 基线 |
| B02 | `v1.5.1..v2.0.2` | 103 | 包结构、CLI/YAML 与 ASE 工作流 |
| B03 | `v2.0.2..v2.1.4` | 42 | MD、依赖治理与实验性 DMF |
| B04 | `v2.1.4..v2.2.0` | 56 | 稳定 Python API 与进程 runner |
| B05 | `v2.2.0..v2.2.4` | 77 | Agent 可观测性与输入健壮性 |
| B06 | `v2.2.4..v2.2.5` | 2 | 周期 SCF 帧身份与非收敛提示 |
| B07 | `v2.2.5..v2.2.6` | 26 | 统一收敛事实与交接契约 |
| B08 | `v2.2.6..v2.2.7` | 110 | GPU runtime、bench 与恒电势候选 |
| B09 | `v2.2.7..v2.2.8` | 20 | DP PyTorch 与运行证据收口 |
| B10 | `v2.2.8..45d8c78` | 5 | 发布后的主线修订 |
| B11 | `45d8c78..065da5a` | 2 | 本轮整改及可恢复收尾 |

### B01：早期脚本与 legacy 基线

早期 NEB、AutoNEB、Dimer 脚本、轨迹保存与解耦；这是历史基线，不代表当前推荐入口。

### B02：包结构、CLI/YAML 与 ASE 工作流

统一 pip 包与 atst run；abacuslite/DP 后端、schema、CCQN/D2S/IRC、端点治理、两阶段 NEB、image-level MPI 和 STRU 输入。

### B03：MD、依赖治理与实验性 DMF

新增 ASE/native ABACUS MD 与后处理；Python 3.10 基线、optional extras、ABACUS preflight；实验性 DMF/D2S、vendored snapshot 和 CI。

### B04：稳定 Python API 与进程 runner

统一配置/运行/结果文档接口，CLI 经 API service 路由；嵌入式 CCQN、MPI 失败同步、wheel/runner 验证与用户/维护者文档分层。

### B05：Agent 可观测性与输入健壮性

进度事件、只读 profiles/plots；端点信任标记、AutoNEB CI 受控错误；atst prepare、轨迹 stress、安装 profile 与 MPI 诊断。

### B06：周期 SCF 帧身份与非收敛提示

vendored abacuslite 周期等价帧匹配；拒绝不匹配结构，并承接 fd7c8da 的提示意图。

### B07：统一收敛事实与交接契约

StageRecord 与英文 advisory；工作流分阶段/方向记录、D2S/DMF 事实、消费者 fixtures、语言与 wheel 门禁。完成状态不等同科学收敛。

### B08：GPU runtime、bench 与恒电势候选

设备与线程预算、隔离 worker、逐 rank 绑定、sidecar、有限 batch/sweep/record、SAI 验证与多轮修复；恒电势开发候选并入。

### B09：DP PyTorch 与运行证据收口

PT 制品/head 诊断、DP 线程键、attempt 归属、只读模型共享/重启回归、MPS 三态证据和通道制品记录。

### B10：发布后的主线修订

2.2.8 发布记录、MPI rebind 证据、恒电势 profile_status 删除及科学引用修正；这些源码修订不能反写成 v2.2.8 wheel 已包含。

### B11：本轮整改及可恢复收尾

周期几何/IDPP、Relax/许可打包、MPI 进度同步、Dimer 与显式 TS 验证、公共 execution 与 backend identity、入门文档，以及合入/备份/分支清理。

完整 [commits.tsv](data/main_history_20261015/commits.tsv) 含批次、完整 SHA、
父提交、作者日期和原始标题，按每批拓扑顺序排列。生成时已断言：无跨批重复、
无漏项，集合与 `git rev-list 065da5a` 完全一致。日期只辅助阅读，以 Git 图为准。

代表性实现提交：B04 `8908b61` / `446805f` / `c9b7908`；B07 `7fb89be` /
`99d8d47` / `c3df747`；B08 `b9fecb0` / `729cf33` / `186e2a3` / `33f7c1b`；
B09 `0fd205a` / `7f335aa` / `1ab1024`；B10 `715056b`；B11 `3163c06` / `065da5a`。
这些是阅读起点，不替代批次完整清单。

## 3. 分支来源与吸收情况

| 原分支 | 保存的 tip | 当前处置 |
| --- | --- | --- |
| `abacus-demo-hardening` | `fd7c8da` | 非祖先；独有 fd7c8da 的意图已吸收；远端已清理，恢复引用保留 |
| `abacuslite-upstream-issue-fixes` | `7f16fe2` | main 祖先；已清理远端，保留恢复引用 |
| `dev/gpu-tuning-closeout-20260922` | `6a7bf94` | main 祖先；已清理远端，保留恢复引用 |
| `feat/atst-prepare-and-trajectory-stress` | `fd214ba` | main 祖先；已清理远端，保留恢复引用 |
| `feat/dmf` | `7fa9bb7` | main 祖先；已清理远端，保留恢复引用 |
| `feat/neb-parallel` | `e42649c` | main 祖先；已清理远端，保留恢复引用 |
| `feature/gpu-node-tuning` | `f03b3e6` | main 祖先；已清理远端，保留恢复引用 |
| `fix/autoneb-ci-endpoint` | `93eb0a6` | main 祖先；已清理远端，保留恢复引用 |
| `fix/issue-31-neb-stru-mobility` | `dea1f4b` | main 祖先；已清理远端，保留恢复引用 |
| `fix/review-findings-20260921` | `0f5fb51` | main 祖先；已清理远端，保留恢复引用 |
| `refactor/unify-structure` | `ecb3161` | main 祖先；已清理远端，保留恢复引用 |
| `fix/ase-ts-dual-interface-20261013` | `3163c06` | 已快进进入 main；分支和归档 worktree 保留 |
| `develop` | `56d46da` | main 祖先；长期分支保留 |
| `v1.4.2` / `v1.5.0` | `70b22a5` / `47c513f` | main 祖先；版本分支保留，与同名/相近 tag 分开理解 |

来源主线：`refactor/unify-structure` 对应包化，`feat/neb-parallel` 对应并行工作流，
`feat/dmf` 对应实验性候选生成，`fix/issue-31-neb-stru-mobility` 对应 STRU/约束，
`fix/autoneb-ci-endpoint` 对应端点信任和 CI 边界，`feat/atst-prepare-and-trajectory-stress`
对应配置反向生成和应力保留；GPU 三条分支是先实施、后审查修复、再收尾的嵌套主线，
不是三份独立功能。abacuslite 分支承担 backend 快照与解析治理。

Git 明确保留的 merge 入口包括 `aa3607d` / `684ac8e`（结构统一）、`28bf91a`
（并行集成）、`234155a`（Issue 31）、`ea36eb8`（DMF）、`7738f2c`（稳定 API）、
`21c5571`（AutoNEB 端点）。其余分支按保存的 tip 和补丁解释，不虚构 merge 节点。

`abacus-demo-hardening` 的独有提交 `fd7c8da` 不在 main 的提交集合内，不能称为
Git 已合并。其兼容性提示意图在 2.2.5 发布说明有明确承接记录，2.2.6 又统一为
StageRecord/advisory；不应重新 cherry-pick 旧中文提示实现。

## 4. 审查结论与后续批次

- **历史不重写。** B01–B09 已形成版本边界；不 squash、不 rebase 已发布历史。
  B11 的 `3163c06` 涵盖 95 个文件，是一次完整集成提交，回退粒度较粗；本报告用
  科学算法、执行/API、打包、文档四条阅读路径解释，不将其伪装成多个既有 commit。
- **发布范围比本地 ahead 范围更大。** origin/main..main 只有两条新增提交，但下一次
  发版还须审阅整个 v2.2.8..main，尤其 B10 的恒电势字段变更。不能只总结 `3163c06`。
- **P2：恒电势跨版本读取不兼容且缺少迁移披露。** `715056b` 删除恒电势
  `profile_status` 输出与必填校验，但保持 `atst-constant-potential-facts-v1`。
  v2.2.8 的 `_CP_REQUIRED_FACT_KEYS` 包含此字段，旧 `_constant_potential_facts_error`
  会拒绝新版 facts。新读旧可行，旧读新会失败；当前用户/API 文档没有移除说明。
  这是已在 origin/main 的历史风险，并非 B11 引入。下一次发布前须决定消费侧迁移、
  兼容字段或协议演进策略，补跨版本回归和发布说明；不能自动归类为“纯增量 patch”。
  本轮不恢复字段或更改 schema。
- **验收分层。** B08/B09 有既有 SAI 证据，但本轮不是新增 GPU 科学验收。
  恒电势是开发候选，DMF 是实验性路径；P6 平台迁移、镜像与消费者联合验收仍开放。
- **可恢复性成立。** 11 条已清理分支均保留恢复引用和 bundle；隔离工作区可用，
  main 原脏树已实际恢复演练。具体路径与日志见
  [系统审查 §11](ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md#11-合入与工作区清理核对2026-10-15)。

建议后续提交按可独立验证的行为边界分批：①恒电势消费兼容裁定及回归；②明确下次
版本的 release notes / 功能矩阵 / 发行门禁；③另在父仓提交 gitlink、镜像和消费验收。
这三批尚未执行，不将本报告当作发布、平台部署或父仓修改的授权。

独立复核范围为 `origin/main..065da5a` 的完整性、版本/发布边界与恢复事实，并扩展
检查 `715056b`：未发现 B11 新的代码/集成问题；上述 P2 得到确认。复核检查了
归档 bundle、恢复副本、分支引用及父仓状态，复用既有测试证据，不新增 GPU 或跨仓
验收结论。历史 B01–B09 没有逐提交重审，不能据此声称没有潜在缺陷。

## 5. 证据入口

- [2.2.0 发布说明](../releases/RELEASE_NOTES_2.2.0.md)：API/runner。
- [2.2.5 发布说明](../releases/RELEASE_NOTES_2.2.5.md)：周期帧身份及旧分支语义承接。
- [2.2.6 发布说明](../releases/RELEASE_NOTES_2.2.6.md)：收敛契约。
- [2.2.7 发布说明](../releases/RELEASE_NOTES_2.2.7.md)：runtime/bench/恒电势候选。
- [2.2.8 发布说明](../releases/RELEASE_NOTES_2.2.8.md)：PT 与证据收口。
- [系统审查](ARCHITECTURE_DISTRIBUTION_REVIEW_2026-10-13.md)：本轮实现、文档与合入验证。

## 6. 当前实现独立系统复核（2026-10-15）

### 6.1 结论与审查边界

**历史批次整理准确，架构演进与项目定位一致；已有整改及本地测试证据成立，但当前
实现仍有可复现缺陷，不能据测试全绿判定可无条件发布。** 本轮确认 9 项问题：
1 项 P1、8 项 P2；其中 CP-01 是 §4 已知问题，其余 8 项为本轮新确认的历史缺陷。
未确认 `3163c06` 新引入的代码回归，不把这句话扩大为对该提交全部路径的证明。

审查冻结在 `5ad5de582313f853ae3d713df2fa5a6f2bd9cd5e`，开始时工作树干净；
该树的 `src/`、`tests/`、`pyproject.toml`、`scripts/` 和 `.github/` 自 `3163c06`
未变。本轮仅更新审查、状态账本与证据，不实施运行代码修复，不提交、不推送、不发布。
工作流按 L2 有界审查执行；四个当前模型家族的独立 reviewer 分别审查科学工作流、
API/结果契约、runtime/bench 和状态证据，主代理核对实现、复现输出并运行本地检查。

- 独立重算全部批次集合与 734 条提交的 SHA、父提交、作者日期、标题，完全匹配。
  当前 HEAD 共 735 个可达提交；较本地 `origin/main=45d8c78` 领先 3 条，较
  `v2.2.8=b6a68f925c8c5c9e8a5afc13890d16f853956e9f` 多 8 条。未刷新远端引用。
- B01–B03 的主要成果是从脚本转向可安装包、统一配置与工作流；B04–B07 建立稳定
  API、进程交接和收敛事实；B08–B09 加入资源约束、运行证据和 DP 制品诊断；
  B10–B11 处理发布后修订及科学输入、执行层、分发问题。此演进没有偏离
  ASE/abacuslite wrapper 定位，也没有把 bench 扩张为集群调度器。
- 当前科学执行由共享 `workflows.execution` 承担，CLI/API 共同消费；后端优先导入、
  vendored 来源、optional extras 和结果事实分层合理。历史 MIC/IDPP、Relax 回退、
  MPI 进度失败同步、Dimer 终止事实及 notice 打包的整改仍有对应回归与产物证据。
- 历史提交没有全部逐行重审。代码审查优先当前可达路径及近几批高风险合约；没有
  新跑 ABACUS/DP GPU、真实 Sella IRC 积分或 IPOPT 求解，没有完成 P6 平台验收。

### 6.2 已确认问题

以下行号绑定冻结源码；优先级表示修复建议，不代表本轮已经修复。

| ID / 等级 | 触发与实际影响 | 源码位置 | 版本归属 |
| --- | --- | --- | --- |
| SCI-01 / **P1** | `make.ts_guess` 生成 `n_images+3` 帧，重启按 `n_images+2` 切带；通常报错，累计帧数恰好整除时会把中间像当作固定初态继续计算。 | [execution.py](../../src/atst_tools/workflows/execution.py) L302–307；[idpp.py](../../src/atst_tools/utils/idpp.py) L745–749 | v2.2.8 已存在；B11 迁移执行位置时保留。 |
| SCI-02 / P2 | 默认 Sella IRC、`direction=both` 且省略 `normalized_trajectory`：schema 注入 `None`，两方向结束后 `write(None, ...)` 失败，manifest 尚未写出。 | [irc.py](../../src/atst_tools/workflows/irc.py) L79–81、121；[config_schema.py](../../src/atst_tools/utils/config_schema.py) L812 | v2.2.8 已存在。 |
| SCI-03 / P2 | descent IRC 对分支执行 `Atoms.copy()` 后丢失 calculator，主轨迹/normalized 轨迹及返回帧无能量和力；分方向轨迹虽有结果，却未列入 manifest。 | [irc.py](../../src/atst_tools/workflows/irc.py) L306–318 | v2.2.8 已存在。 |
| SCI-04 / P2 | DMF 不核对端点有序元素，`HO` → `OH` 被接受，实际末像变成首帧元素 `HO` 配末帧坐标，候选对应错误终态。 | [dmf.py](../../src/atst_tools/workflows/dmf.py) L163–169；[PyDMF](../../src/atst_tools/external/pydmf/dmf/dmf.py) L257–259 | v2.2.8 已存在；实验性路径。 |
| IO-01 / P2 | `artifact_manifest: reports/artifacts.json` 时，工作流记录相对运行目录的 `neb.traj`，summary 相对 manifest 目录读取，成功产物被报告为缺失。 | [execution.py](../../src/atst_tools/workflows/execution.py) L530；[summary.py](../../src/atst_tools/utils/summary.py) L289–296 | 早于 v2.2.5；v2.2.8 已存在。 |
| BENCH-01 / P2 | 显式 `workdir` 与默认配置目录即使解析为同一实际目录，也因检查键前缀不同而通过；路径别名同样可绕过，运行产物可能互相覆盖。 | [batch_runner.py](../../src/atst_tools/bench/batch_runner.py) L233–242、316–330 | v2.2.8 已存在，之后文件未改。 |
| BENCH-02 / P2 | `case.env` 在绑定后覆盖设备、线程、缓存变量；分配卡 0 / 线程 1 的子进程可收到卡 7 / 线程 99，调度预算与实际运行不符。 | [batch_runner.py](../../src/atst_tools/bench/batch_runner.py) L555–571 | v2.2.8 已存在，之后文件未改。 |
| BENCH-03 / P2 | `case_id` 仅检查非空：`../escape` 使报告路径越出输出目录；含空格的 ID 被进度解析截断，成功 worker 被记为失败。 | [batch_runner.py](../../src/atst_tools/bench/batch_runner.py) L133–136、979–989、1284；[batch_worker.py](../../src/atst_tools/bench/batch_worker.py) L116–146 | v2.2.8 已存在，之后文件未改。 |
| CP-01 / P2 | 沿用 facts-v1 却删除 `profile_status`，旧版 restart 拒绝新版 facts；允许重算时丢弃 checkpoint 电子数，新读旧仍可行。 | [constant_potential.py](../../src/atst_tools/calculators/constant_potential.py) L57、80–107、572 起 | `715056b`，发布后 B10 已知风险。 |

**科学路径复现**：SCI-01 使用真实 ASE EMT/FIRE，`n_images=3` 首次运行得到 6 帧，
重启却要求 5 帧；5 个 6 帧带共 30 帧通过整除检查，初态 x 从 0 变为 0.15。
SCI-02 仅用受控优化器替代 Sella 积分，两方向完成后由实际 workflow/schema/ASE IO
复现 `TypeError`，不将替身完成称为真实 IRC 科学验收。SCI-03 使用 EMT/FIRE，
分方向轨迹读回有 calculator，交付主轨迹读回为 `calc=None`。SCI-04 使用真实
PyDMF 初始化，仅替代缺失的 `cyipopt.Problem` 空壳，未执行 solve。

**执行与消费复现**：IO-01 经公共 `run_workflow` 和 ASE NEB/FIRE/LennardJones
完成，实际 `neb.traj` 存在而 summary 为 `n_missing=1, complete=False`；
有限步数未收敛已正确记录为 advisory，不是本项缺陷。BENCH-01 比较接受后的真实目录；
BENCH-02 通过真实 batch 子进程捕获环境，普通无 runtime 请求路径没有后续纠正；
BENCH-03 经实际 shared-worker 进度消费链确认 worker exit 0、有结果文件却报告失败。
没有实际使用 GPU 7，越界路径与环境覆盖仅在 CPU 临时目录中验证，不宣称权限提升。
CP-01 用当前 CPU fixture 生成真实 facts，再在内存加载 tag `v2.2.8` 读取器，旧
restart 抛缺字段错误，当前 restart 返回约 10 个电子；补回旧字段后新旧均接受。

### 6.3 处置顺序与验收边界

1. **先修 SCI-01**：从实际初始链或可信带宽元数据恢复图像数，统一 TS guess 计数。
   覆盖非整除报错与恰好整除的静默错带，验证重启前后固定端点不变。
2. **保护输入和交付**：SCI-02 明确 `None` 的禁用/默认语义并统一写出条件；
   SCI-03 冻结已有能量/力，避免额外后端计算；SCI-04 在构建路径及 calculator 前
   核对原子数和有序原子序数，必要时经显式对齐；IO-01 统一写读路径基准，兼顾
   默认布局和已有消费者，不能只改 summary 后破坏另一种历史布局。
3. **修 bench 隔离**：用实际规范化路径做冲突判断，限制资源/绑定/缓存环境覆盖，
   对 case ID 建立安全单路径段及进度编码规则；补真实子进程的预算与报告验证。
4. **再裁定发布兼容**：CP-01 需要兼容字段或明确的消费迁移方案。单独改 schema
   字符串不会让旧消费者自动接受；应补双向跨版本回归、发布说明与版本语义裁定。

另有两项维护观察，不计入上述 9 项代码问题：

- **发行自动检查缺口**：`tests/unit/test_distribution.py` 依赖 `ATST_DIST_SDIST` /
  `ATST_DIST_WHEEL`，默认跳过；`.github/workflows/publish-pypi.yml` 在构建前跑测试，
  构建后仅跑 Twine/API 检查。本次手工 notice 内容检查有效，但 CI 未自动防止同类
  漏打包复发。建议在上传分发物前显式运行这两项；这不是当前修复产物仍漏许可的证明。
- **状态文字漂移已在本轮校正**：feature matrix 保留了 9 月 20 日“未发布/不在 PyPI”
  的历史段落，但候选已随 2.2.7/2.2.8 包分发；站点重绑定三次复验与后续上游确认
  也需区分。当前矩阵明确发行、开发候选、未发布整改和 P6 平台验收四个边界。

Ruling: 本轮保留可复现问题及下一步验收，不顺带修改运行逻辑或选定 CP 协议策略。
原因是用户要求系统审查确认，输入计数/消费兼容的修复应有独立可审阅范围；代价是
这些缺陷仍存在，下一次发布评估必须明确处理，不能把审查完成写成整改完成。

### 6.4 当前验证与证据

| 检查 | 结果 | 结论范围 |
| --- | --- | --- |
| 全部历史分区及每条元数据重算 | 734/734 匹配，无重复/漏项 | 原清单止于 `065da5a`，不应伪增到 735。 |
| 当前 checkout 默认 tests | 1284 passed、32 skipped、3 warnings，exit 0 | Python 3.10.19 / ASE 3.28.0；不覆盖本轮全部新增边界。 |
| `ATST_RUN_MPI_TESTS=1` integration | 27 passed，exit 0 | 实际本地 MPI launcher、受控 calculator/进程 fixture，不是 GPU 科学作业。 |
| 三组额外 CPU 复现 | 均 exit 0，症状得到确认 | 复现脚本的成功表示缺陷存在，不表示缺陷已修复。 |
| 既有修复 wheel 与当前源码比较 | 111 个包内文件已对照，仅 execution.py 空白变化，AST 相同 | wheel SHA 与既有归档一致，可复用相应构建/notice/Twine/clean-install 证据；未新发布。 |
| 文档治理、文档测试与差异检查 | governance 通过，26 tests passed，diff check 通过 | 宽泛冲突扫描仅命中 pytest 日志分隔标题；精确冲突标记检查无命中。 |

32 项默认 skips 为 27 项 opt-in MPI、2 项需指定分发物的检查，以及完整 ABACUS
run-dir、toolbox 等价、真实 DP 模型各 1 项。MPI 27 项本轮另跑通过；其余 5 项没有
冒称本轮默认验收通过。没有读取实时 PyPI/Gitee/远端 CI 状态，发行包边界复用先前
归档的官方 wheel 检查；其 PyDMF LICENSE 缺失不能被本地修复产物的检查覆盖。

证据目录：[confirmation](data/main_history_20261015/confirmation/)：

- [命令与结果汇总](data/main_history_20261015/confirmation/verification.json)、
  [文档测试](data/main_history_20261015/confirmation/docs-tests.log)、
  [文档治理](data/main_history_20261015/confirmation/governance.log)。
- [历史重算](data/main_history_20261015/confirmation/history-check.json)、
  [当前默认测试](data/main_history_20261015/confirmation/default-tests.log)、
  [当前 MPI 测试](data/main_history_20261015/confirmation/mpi-tests.log)。
- [科学复现脚本](data/main_history_20261015/confirmation/science_findings.py)与
  [输出](data/main_history_20261015/confirmation/science-findings.log)；
  [API/CP 复现脚本](data/main_history_20261015/confirmation/contract_findings.py)与
  [输出](data/main_history_20261015/confirmation/contract-findings.log)；
  [bench 复现脚本](data/main_history_20261015/confirmation/runtime_findings.py)与
  [输出](data/main_history_20261015/confirmation/runtime-findings.log)。
- [分发内容对照](data/main_history_20261015/confirmation/distribution-comparison.json)；
  构建/notice/Twine/clean-install 原始日志仍在
  [整改证据](data/system_review_20261013/remediation_20261015/verification.json) 所列位置。

复现命令（从仓库根执行，产物进入系统临时目录）：

```bash
PYTHONPATH="$PWD/src" conda run -n atst-dev python docs/reports/data/main_history_20261015/confirmation/science_findings.py
PYTHONPATH="$PWD/src" conda run -n atst-dev python docs/reports/data/main_history_20261015/confirmation/contract_findings.py
PYTHONPATH="$PWD/src" conda run -n atst-dev python docs/reports/data/main_history_20261015/confirmation/runtime_findings.py
```

审查后的状态为：核心工作流及交付架构已有稳定基础；9 项缺陷待独立整改；
恒电势仍为开发候选、DMF 仍为实验性；P6 共享 image/runner 迁移、下游 E2E 与平台
验收仍开放。源码未发布整改、已发行 2.2.8 和站点历史验证须继续分别描述。


## 7. 九项问题后续整改（2026-10-16，未发布）

§6 的九项缺陷已在隔离分支 `fix/review-nine-findings-20261016` 修复，并通过逐包
回归、独立源码终审和集成检查。冻结源码默认套件 **1347 passed、32 skipped**，
本地真实 MPI **27 passed**，clean-wheel public API 门禁通过。bench 终审额外发现
预检与 worker 的相对配置路径不一致，已统一绝对路径并完成 RED/GREEN 与复审。

逐项行为、旧轨迹/manifest/facts 迁移、原始日志及源码哈希见
[九项问题整改验收](MAIN_REVIEW_FIXES_2026-10-16.md)。本节只更新后续处理状态，
不重写 §1–6 的历史证据。原 main 与发布的 2.2.8 均未改变，未执行新的 GPU 科学
验收或平台 P6。
