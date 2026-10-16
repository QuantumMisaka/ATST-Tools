# GPU P6 接手：默认入口与 image 内 MPI 多卡

**版本**: 2026-10-16
**日期**: 2026-10-16
**状态**: 可行性核验完成；平台修复、真实计算与发布验收尚未完成
**责任人**: ATST-Tools maintainers / ABACUS Toolbox maintainers

生命周期：review；报告级别：L3。执行状态由既有
[GPU 计划 P6](../superpowers/plans/2026-09-20-atst-gpu-node-tuning-plan.md#p6atst-交付及后续平台接入)
维护；本报告承载接手时的代码事实和离线复现，不另建实施计划。

## 1. 维护者决定与结论

2026-10-16 维护者要求接手 P6，并明确：“以默认行为为主，探索有没有办法启动
高阶的单 image 内部 mpirun 多卡运行行为”。因此不把原计划的多 image 共享少量 GPU
作为此次优先开发入口，也不重释 `n_gpu` 或修改默认并行行为。

**image 内 MPI 多卡已有基础入口；将它做成可靠的平台交接路径，优先于新增混合并行。**
以下三种拓扑必须分别描述，不能都叫“多卡 NEB”：

| 拓扑 | Python / image 层 | ABACUS 层 | 当前状态 |
| --- | --- | --- | --- |
| 既有 image 并行 | 一个 Python rank 对应一个活动 image | 每图单 rank、单卡 | 保留；平台启用 `parallel=True` 时沿用既有契约 |
| image 串行、每图多卡 | 单 Python runner，依次计算各 image | 每次 `mpirun -np G abacus` 使用同一组 G 卡 | 代码已有通路；跨任务资源和 mapping 交接尚有缺口 |
| 多 image 并发、每图多卡 | M 个并发 image，各持互斥 GPU 组 | 每组 G 个 ABACUS ranks，总计 M×G 卡 | 高阶探索；当前平台明确拒绝直接嵌套 MPI，尚无实现或验收 |

平台 prepare 的 `parallel` 参数默认 False，资源默认 1 GPU；“保留一图一卡”指
启用既有 image-parallel 模式后的映射规则，不表示本轮把默认改成 image 并行。
单个 image 内的多卡能力不等于同时运行多个 image，也不自动证明性能或显存收益。

## 2. 基线与可复现证据

父仓 `app-tools`：`45f36bd14bc28c2dd2c6602d0bcde8175563295a`；所读 ATST 工作树：
`767b0cd98244353660a6528147140d9f53b7086a`。父仓 index 的 ATST gitlink 仍为
`45d8c78de8a605c127cd4c36742d0c9525a0a6b6`，本次没有同步它。父仓已有他人在途文档，
本次未修改。设计与证据写在独立 ATST 分支 `docs/gpu-p6-inner-mpi-20261016`。

环境预检：`abacus-env` 可导入 `adam_community`、`ase`、`docstring_parser`，
`adam-cli --help` 成功；解释器为 `/home/james/apps/miniforge3/envs/abacus-env/bin/python`，
ASE 3.29.0。以下现有测试在父仓 `toolbox/ABACUS` 执行，**58 passed / 22 warnings**
（警告为 ASE/NumPy 弃用消息），不是 GPU 或平台端到端验收：

```bash
conda run -n abacus-env python -m pytest \
  tests/test_transition_runtime_policy.py \
  tests/test_transition_neb_parallel.py \
  tests/test_atst_runtime_pass_through.py -q
```

专项离线复现：[脚本](data/gpu_p6_20261016/probe_platform.py)、
[机器可读结果](data/gpu_p6_20261016/probe.json)。在 ATST 根执行，最后一个参数为
待审父仓的 ABACUS 路径：

```bash
conda run -n abacus-env python docs/reports/data/gpu_p6_20261016/probe_platform.py \
  /path/to/app-tools/toolbox/ABACUS
```

该脚本调用真实的资源 hook、配置构造、mapping 校验和 run 适配器，但替换平台框架、
SIF 刷新、check-input/runner 子进程及绘图，不执行 ABACUS、MPI、GPU 或调度操作。
`adapter_status_with_stub_runner=success` 只证明抵达模拟执行边界，不表示计算成功。
结果保存关键源文件哈希及双方 HEAD；它是当前基线的审计复现，不是期望缺陷永远存在的回归测试。

## 3. 当前入口实际做了什么

以下父仓路径均相对 `toolbox/ABACUS/`。

1. `abacus_transition.py` 的 prepare 工具继承
   `toolkits/transition_state.py::_TransitionRuntimeMixin.resAlloc`：显式 `n_gpu=4`
   实测申请 GPU=4 / CPU=16。只看 class 的 GPU=1 会误判；本轮初步猜测已据此纠正。
2. `toolkits/transition_atst.py::prepare_transition_workflow` 在 `parallel=False`、
   `n_gpu>1` 时，通过 `_build_abacus_calculator` 调用
   `utils/env_setting.py::setting_abacus_calc_env`。在模拟的合法 4 卡 mapping 下生成
   `command: mpirun -np 4 --map-by ppr:4:node:pe=4 abacus`、`mpi: 1`、`omp: 4`。
   此处 `mpi: 1` 用于避免二次包装；实际 ABACUS rank 数已包含在 command 中。
   GPU 分支由 `n_gpu` 决定 ranks，`mpi_mode` 并不是第二个多卡开关。
3. ATST `calculators/factory.py::_build_abacus_command` 保留已显式给出的 launcher；
   `workflows/execution.py` / `mep/neb.py` 在非 image 并行路径逐图调用 calculator。
   因而“单 Python runner → 每图 mpirun 多卡”已有基本组合，无需启动外层 Python MPI。
4. `parallel=True` 时，平台 `_validate_image_parallel_calculator` 明确要求
   `command='abacus'`、`mpi=1`。不能把 command 换成 mpirun 就绕过此契约。

现有参数组合示意为 prepare 显式 `parallel=False, n_gpu=4`，run 同样显式
`n_gpu=4`；它不是已通过本轮平台验收的操作配方，必须先补齐下列交接检查。

## 4. 已复现的交接缺口

| 项目 | 观察 | 影响及修复目标 |
| --- | --- | --- |
| run 的内层资源一致性 | 保存的 command 为 `mpirun -np 4 ...`，run 省略 `n_gpu` 时 resAlloc 返回 GPU=1，实际适配器仍抵达模拟正式 runner | 在 check-input 和正式运行前核对配置意图与本次显式资源；不得只校验 image-parallel 配置 |
| prepare mapping 跨任务固化 | prepare 将当前 `--map-by` 写入 YAML；非 image-parallel run 不调用 mapping 刷新。模拟更换本次 mapping 后，仍传出旧 command | 本次执行使用当前 allocation/mapping；不能把上游任务的站点拓扑当作可移植科学配置 |
| prepare 与未来资源意图耦合 | 正常 prepare hook 会申请未来意图对应的卡数；若只在 1 卡 allocation 中构造 4 卡配置，现有真实 gate 拒绝 | 不是 prepare 必然错误分配；是否将纯准备改为离线资源声明属于额外设计选择，不在本轮自动改变 |

前两项是优先修复边界。修复应复用现有资源 helper、API runner 和 prepared 输入链，
不在 `resAlloc` 内读 YAML，不在 run 静默改 solver，不把未知站点 wrapper 的任意 command
解析成可信资源事实。受管理的生成命令和用户自定义命令的兼容边界需在实施前收束。

另有既有 image 模式的计数回归需要纳入 P6 默认路径复核：父仓
`_neb_interior_image_count` 直接使用 `make.n_images`，未计入 TS guess 额外的内部图；
本次仅源码定位，未将其宣称为已复现或已修复。

## 5. 高阶混合并行的可行方向

“ABACUS 单次运行支持 MPI 多卡”是后端前提；在 image 层组合它还需要互斥 GPU 组、
CPU/OMP 预算、进程回收、结果顺序和端点归属。Python communicator 不能直接当作
外部 ABACUS 可执行文件的 communicator。直接把外层 rank 数乘以 G 也不符合当前
ATST `validate_image_parallel_world` 的一 rank 一内部图契约。

建议探索一个**不处在外层 MPI world 中的协调进程**，为每个活动 image 的独立 worker
分配互斥 GPU 组，再由每个 worker 启动本组的 ABACUS MPI 作业。科学计算继续通过
ASE calculator / 既有 NEB、AutoNEB 流程交接；具体 worker 机制尚未选定，不能据此
宣称现有 API 已提供进程池。CPU 线程、进程环境必须按 worker 隔离，不能通过并发修改
全局 `os.environ` 实现 GPU 分组。

对照候选是站点支持的独立 Slurm job steps，由平台 launcher 承担组隔离；是否适合
当前 SIF/SAI 链路需要真实 gate。直接递归 mpirun 不作为首选：Open MPI 曾明确报告
递归启动检测，移除一两个环境变量不能证明资源归属、PMIx 会话和清理安全。
参见 [Open MPI 的递归调用问题记录](https://github.com/open-mpi/ompi/issues/8171)
及 [Slurm 启动支持文档](https://docs.open-mpi.org/en/main/launching-apps/slurm.html)。
这些资料提供风险依据，不是本项目混合并行的验收证据。

ASE 的 [NEB 实现](https://docs.ase-lib.org/_modules/ase/mep/neb.html) 提供串行、线程与
MPI 的图像求值路径，但并不替 ATST 分配 GPU 或创建相互隔离的 ABACUS MPI 作业。
这里的方案属于从现有接口推导的候选设计，尚未实现。

每组 G 卡是否能解决显存不足取决于 solver、矩阵分布和体系；不能把多卡总显存直接
当作单图可用显存。性能验收要比较逐步能量/力、收敛事实、单 band 延迟、卡时与峰值
显存，不能只看 GPU 利用率。本轮保留所有科学参数和默认行为。

## 6. 交付边界

本次完成用户方向固化、既有入口追踪、58 项离线基线检查和专项复现；未改生产源码。
下一实施包先修复 §4 的受管理串行 image 内 MPI 交接，再做同一 allocation 下的真实
单图多卡验收；混合并行的入口与 worker 架构需经候选验证后另行确认。
原“多图共享少卡”平台扩展保留 deferred，ATST 已有能力不撤回；MLIP 的平台资源
所有权与 DP 条件项继续独立跟踪。

真实计算、本地 SAI E2E、父仓 gitlink 同步、平台发布各自按 owner 授权执行，不能以
离线测试代替。GPU P6 保持开放，不能归档为已完成。

独立 reviewer 已复核本轮文档、探针及相关调用路径，未发现 P1/P2 实质问题；
确认 prepare hook 的真实资源行为与模拟执行边界的表述。文档治理脚本与既有
`tests/unit/test_docs_governance.py` 的 26 项检查通过，差异空白检查通过。
