# ATST-Tools 源码架构复核

**版本**: 2026-10-16（源码基线 `5ad5de5`）
**日期**: 2026-10-16
**状态**: 审查完成；下文保留修复前基线，后续整改另见验收报告
**责任人**: ATST-Tools maintainers

生命周期：review；报告级别：L3。本报告判断现有结构是否适合继续迭代，长期模块说明
维护在 [源码导览](../../src/atst_tools/README.md)。逐项缺陷与复现仍由
[主线审查 §6](MAIN_BRANCH_REVIEW_2026-10-15.md#6-当前实现独立系统复核2026-10-15)
承载，本轮没有修复或关闭这些缺陷。

后续九项问题的实现与验收见 [整改报告](MAIN_REVIEW_FIXES_2026-10-16.md)。本报告的
导入快照、源码定位及测试数字继续对应修复前基线；§4 前两项的关闭状态由整改报告
维护，第三项仍是随扩展处理的设计债。

## 1. 总体判断与范围

**现有架构适合继续发展：保留单一 Python 包内的模块分工、CLI/API 共享科学执行、
ASE calculator 接口，以及 runtime 与 bench 的分工。当前改进重点是共性契约在
不同工作流中的一致落实。** 未发现需要整体重构、拆包或建立插件框架才能解决的问题。

ATST-Tools 的核心是围绕 ASE 组织 ABACUS/DP 的过渡态及相关科学工作流。
Relax、Vibration、IRC 为路径生成和候选验证提供配套能力，MD 提供相关模拟能力；
运行证据及 bench 支持维护者比较运行配置。模块存在不等于同等成熟：DMF 仍为实验性，
恒电势仍为开发候选，具体范围以 [功能矩阵](FEATURE_STATUS_MATRIX.md) 为准。

源码基线为 `5ad5de582313f853ae3d713df2fa5a6f2bd9cd5e`；工作树已有源码导览、
历史审查及其文档治理改动。本轮再次确认生产 Python 源码、tests、pyproject 和 CI
配置相对 HEAD 无差异，因此可复用同基线的测试与缺陷复现。没有重新审阅全部历史
提交，也没有执行新的 ABACUS/DP GPU、真实 Sella IRC、IPOPT 求解或平台验收。

## 2. 实际依赖与启动顺序

对非 vendored 的 72 个 Python 文件做 AST 显式导入扫描：包根 1，api 5，bench 5，
calculators 7，mep 6，runtime 7，scripts 4，utils 28，workflows 9。
逐文件哈希、依赖边和动态导入位置保存在
[依赖快照](data/architecture_review_20261016/import-audit.json)，扫描方法见
[只读脚本](data/architecture_review_20261016/import_audit.py)。

扫描没有发现显式模块导入构成的循环；没有 API → CLI 的显式依赖，也没有科学模块
导入 bench。此结论只覆盖 AST 中的 `import` / `from`，不将目录级的双向引用
误判为模块循环，也不把静态图当作动态加载与所有执行路径均无环的证明。
`api.__getattr__` 和轻 CLI 的动态加载另对照其目标映射及入口测试核验。

| 实际关系 | 判断与边界 |
| --- | --- |
| scripts → api → workflows；workflows → mep/calculators | 主执行方向清楚。`scripts.main` 对 execution helper 的直接导入兼顾历史兼容，API 不经 CLI 执行科学计算。 |
| api → mep/utils | 嵌入式 `run_ccqn` 和只读分析入口有直接调用，不要求所有能力都经过一个统一工作流类。 |
| bench → api/runtime | bench 消费既有计算与启动能力，科学模块不需要知道批次的存在。 |
| runtime.errors → api.models | 复用轻量错误基类；属于受控交叉依赖，不是 runtime 调用 API 服务。应保持 models 不加载科学库。 |
| utils → calculators/runtime | `abacus_io` 使用后端接口；geometry/neb_endpoints 清理、读写恒电势事实。utils 已包含领域逻辑，不能当作严格无依赖的底层库。 |

六个新进程分别导入包根、api、api.runner、scripts.cli、runtime.cli_dispatch、
runtime.launch，均未加载检查清单中的 numpy/ASE/DP/Torch/TensorFlow/JAX/Sella/mpi4py。
这是 bootstrap 轻入口的证据。访问 `api.run_workflow` 会加载服务与科学依赖，嵌入式
宿主仍须先建立环境；轻导入不意味着可以在科学库初始化后任意重绑设备。
代码依据：[api 惰性导出](../../src/atst_tools/api/__init__.py)、
[CLI 轻入口](../../src/atst_tools/scripts/cli.py)、
[runner](../../src/atst_tools/api/runner.py)、
[runtime 错误](../../src/atst_tools/runtime/errors.py)。

## 3. 各边界是否合适

### 3.1 科学步骤、算法与计算后端

`workflows` 组织输入、阶段和交付；`mep` 提供搜索组件。但当前 `mep` 也包含任务
适配器：例如 `CCQNOptimizer` 继承 ASE Optimizer，`AbacusCCQN` 负责装配及输出；
AutoNEB 也有 engine 与 runner。这种组合可以工作，代价是新增代码时需判断其属于
算法还是编排，不能仅按目录名机械决定依赖。

建议新增算法消费 ASE Atoms/calculator，新增复合任务进入 workflows；已有适配器
只在实际扩展需要时拆分并保留兼容入口。证据见
[ccqn.py](../../src/atst_tools/mep/ccqn.py)、
[autoneb.py](../../src/atst_tools/mep/autoneb.py) 和
[d2s.py](../../src/atst_tools/workflows/d2s.py)。

calculator 工厂负责命令、模型和配置装配，ABACUS 的 ASE 行为由 abacuslite 提供，
DP 使用 DeePMD-kit 的 calculator，恒电势以 calculator 装饰器加入固定几何电子循环。
独立 abacuslite 优先；仅包本身缺失时回退快照，不掩盖其内部依赖缺失。
这是适合项目定位的接口复用。`BaseCalculator` 当前没有生产代码引用或继承，不能把
它描述成所有后端实际遵循的基类；无需为形式统一再包一层 ASE calculator。

证据见 [factory.py](../../src/atst_tools/calculators/factory.py)、
[后端解析](../../src/atst_tools/calculators/abacuslite_backend.py)、
[dp.py](../../src/atst_tools/calculators/dp.py)、
[base.py](../../src/atst_tools/calculators/base.py)。external 保留上游来源、基线、补丁
与许可记录的方式应继续执行；本轮未重新进行快照 drift 或分发构建验收。

### 3.2 配置与扩展

Pydantic schema 负责类型、默认值和组合校验，归一化配置由 CLI/API 共同使用。
目前 12 类 calculation 采用显式判别联合、类型表和服务调度分支，扩展步骤可追踪。
同时，模板、返回结构适配、合成产物及可选 profiles/plots 存在相应分支：新增类型
不是只加一个类或 schema 即完成。

当前规模保留显式映射更容易审查。先在 [HANDOVER §2](../developer/HANDOVER.md#2-新增或修改-workflow)
补齐接入清单；只有出现实际第三方工作流接入或反复的注册遗漏，再评估集中描述表。
本次不新增 plugin registry。证据见
[schema](../../src/atst_tools/utils/config_schema.py) 的 `CalculationConfig` / 类型表、
[services](../../src/atst_tools/api/services.py) 的 dispatch/result helpers 和
[CLI 模板](../../src/atst_tools/scripts/main.py)。

### 3.3 科学事实与结果交接

现有所有权安排应保持：执行阶段记录优化器的收敛事实，复合任务聚合子阶段，
manifest 登记产物，API 将其转换为调用方结果，runner 提供 root-rank JSON。
API 缺少工作流 manifest 时合成的阶段保持 `converged: null`；它不凭轨迹或退出码
补出科学收敛。CLI 保留原有产物策略，API 要求可交接结果，这一兼容差异有明确目的。
执行完成、优化收敛及显式 TS 验证继续分开。

风险在落实方式：IO-01 的生产者与 summary 对相对路径基准理解不同；CP-01 的
事实内容变更破坏旧读者，却沿用 v1 标识；IRC 保存轨迹时遗漏有效冻结结果。
这些问题说明需要共同的路径、身份、缓存和版本规则。已有
[artifacts.py](../../src/atst_tools/utils/artifacts.py)、
[convergence.py](../../src/atst_tools/utils/convergence.py)、
[geometry.py](../../src/atst_tools/utils/geometry.py) 和
[neb_endpoints.py](../../src/atst_tools/utils/neb_endpoints.py) 可作为收敛位置。
不能直接把“复制有效 SPC 缓存”和“主动求能量/力”合为同一操作；后者可能启动计算。

API finalizer 目前可补全 manifest、追加 plots/runtime/backend identity，然后形成
结果；这是有意安排的多阶段写入。改动时应保持科学阶段事实由工作流拥有，避免
为了“单一写入者”大改所有工作流；先修复实际可见的路径与协议问题。

### 3.4 runtime、bench 与外层资源

职责链为：调用方/平台取得 allocation → bench 在其内安排有限 case → runtime
在科学库初始化前建立设备和线程环境 → API 执行科学任务。普通单次运行绕过 bench。
bench 默认调用进程 runner，共享 worker 则在同一子进程依次调用 API；两种模式
有不同隔离和模型复用成本，维持显式选择合理。

bench 随包安装有实际用途：在与被测代码相同版本的环境中复现有限性能实验。
它使用 API/runtime，没有资源申请、持久队列或平台控制面，不需要因目录解释困难
而移出包。examples 提供负载，tests 验证行为，workdir 承载科学计算读写，`--out`
收集实验报告，这四者各有用途。

三项 bench 缺陷属于入口约束不足：目录应按解析后的实际位置比较；case_id 应是
安全单路径段并能无歧义地进入进度协议；用户 env 与配置不应改变已经计入预算的
设备/线程/缓存事实。整改应核对最终 worker 环境与 YAML runtime 请求的一致性。
这需要修复 bench 输入与启动边界，不要求引入生产调度平台。

## 4. 迭代顺序与设计债

以下为后续实施建议，不表示本轮已经完成对应修改。

| 顺序 | 工作 | 完成判据 |
| --- | --- | --- |
| 1：科学正确性 | SCI-01 分段 NEB 重启；SCI-02/03 IRC；SCI-04 DMF 输入身份 | 使用既有复现建立修复前后对照；覆盖完整带宽/错误切带、缺省路径、冻结能量/力及元素顺序。 |
| 2：交接与隔离 | IO-01、CP-01、BENCH-01/02/03 | 生产者与消费者对路径及版本一致；兼容策略明确；别名目录、受管 env、case-ID 不能绕过约束。CP 兼容决策须在下一次发布前完成。 |
| 3：随实际扩展整理 | MEP 适配器与算法混合、utils 领域依赖、service 多处分支、未使用的 BaseCalculator | 在触及相应能力时收敛职责；不破坏已公开入口，新增能力有明确接入与验证位置。 |

第一、二行沿用主线审查的 **1 项 P1、8 项 P2**，不是本轮新增的九项问题。第三行
属于维护成本和一致性风险，不应直接标成已发生故障或发布阻断项。

Ruling：保留现有模块布局，以实际契约缺陷驱动局部整改。理由：主执行方向清楚，
已有共享组件可承载修复，目前无跨项目插件或生产调度需求。若未来大量工作流
扩展导致映射频繁遗漏，再评估集中注册与模块拆分；过早重构的代价是迁移现有
CLI/API、manifest 及科学行为，却不能替代对九项边界问题的验证。

## 5. 证据与交付边界

- 独立复核覆盖执行/API 契约、科学/后端和 runtime/bench 三个范围；结论均为保留
  主架构、局部收紧契约。静态检查与已复现缺陷相互印证，没有据此宣布所有路径正确。
- 复用 2026-10-15 同生产源码的默认套件 **1284 passed、32 skipped** 和本地真实
  MPI **27 passed**；原始命令和日志见
  [验证记录](data/main_history_20261015/confirmation/verification.json)。不是本轮重跑，
  也不是新的 GPU 或发布验收。
- 本轮入口与共享执行测试 **21 passed**、六个新进程导入检查通过，文档治理检查
  及文档测试 **26 passed**，差异与实际冲突标记检查通过；结果记录在
  [检查日志](data/architecture_review_20261016/checks.log)。
- 生产运行代码保持不变；交付为本报告、依赖快照及源码 README/HANDOVER/状态账本
  更新。既有九项问题、平台 P6 与发布兼容决策仍按各自入口跟踪。
