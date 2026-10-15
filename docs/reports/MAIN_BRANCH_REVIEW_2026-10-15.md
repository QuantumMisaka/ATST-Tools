# main 分支更新与提交批次审查

**版本**: 2026-10-15
**日期**: 2026-10-15
**状态**: 维护
**责任人**: ATST-Tools maintainers

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
