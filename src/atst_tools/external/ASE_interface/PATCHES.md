# vendored abacuslite 语义补丁清单

与上游逐字节差异的**语义补丁**（非结构适配：相对导入/删内嵌测试/注释 churn 由 checker 机械归一化）。每条登记后由 `scripts/check_abacuslite_snapshot.py` 归一化。

| 文件 | 位置 | 补丁 | 登记日期 | 上游状态 |
| --- | --- | --- | --- | --- |
| `abacuslite/io/legacyio.py` | SinglePointDFTCalculator 构造 | `efermi=ener['E_Fermi']` → `ener.get('E_Fermi')`（running log 缺 E_Fermi 容错） | 2026-08-04 | 未上游化（待 PR） |
| `abacuslite/io/latestio.py` | 同上 | 同上 | 2026-08-04 | 未上游化（待 PR） |
| `abacuslite/core.py` | `AbacusTemplate.read_results` | 将 `SinglePointDFTCalculator` 的 `fermi_level`/`eFermi` 显式保留为后端中立的 `efermi` 结果字段，避免 GenericFileIOCalculator 边界丢失恒电势求根所需的 Fermi 能级 | 2026-09-20 | 未上游化（待 PR） |
| `abacuslite/io/legacyio.py` | band parser | 容差块（`_legacy_band_parser_tolerant_block`） | 2026-05-10 | 已由 #7588 上游化（checker 归一化保留至基线推进后清理） |
| `abacuslite/core.py` | SCF 帧身份匹配与 STRU 文件名交接 | PBC-aware 周期性整数晶格平移判定；校验 cell/元素/映射/有限值，损坏帧 fail-closed；保留原生 relax/md | 2026-09-16 | 未上游化（待 PR） |
| `abacuslite/core.py` | `_select_scf_frame_for_structure` 默认值 | cell 分量与周期等价 Cartesian 残差的绝对容差从 `1e-4` 收紧到 `1e-5 Å`，`rtol=0`；原显式 `atol` 调用继续生效，原生 relax/md 不变 | 2026-10-13 | 未上游化（待 PR） |
| `abacuslite/io/generalio.py` | `_read_kpoint` | `0b01ed2`：坐标正则的小数部分改为非捕获组，并将权重读取从捕获组 6 修正为组 3（保证 point-KPT 小数坐标精确往返） | 2026-08-12 | 未上游化（待 PR） |
| `abacuslite/io/latestio.py` | `read_abacus_out` 的 `running_md.log`/`MD_dump` 分支 | 将 MD_dump 无量纲晶胞乘 `alat_in_angstrom`，MD Atoms 设 `pbc=True`；running-log SCF/relax 数据保持原有单位和 PBC 返回契约 | 2026-10-17 | 独立 abacuslite 候选 `dc020430bfa8c76edfc275b0f1eab2f12cc29d4d` 已修；vendored ABACUS 基线未推进 |

## 维护触点

- **core.py 帧选择（周期性误拒修复）**：接受每个原子独立的整数晶格平移，返回最新匹配
  的原始帧，不 wrap 或改写坐标。所有候选先检查数据完整性；有效的错位/错胞候选可跳过，
  损坏候选拒绝整个 SCF。`write_input` 的实际 STRU 文件名随 calculator 状态交接。
  原始 2026-08-04 多帧防御补丁与本次周期性误拒 bug fix 分开理解；本修复不改变
  RAW 力、约束投影或科学收敛判据。
- **2026-10-13 默认容差收紧**：以真实多帧 SCF 夹具、Cartesian/Direct/latestio
  回归及 `9e-6`/`11e-6`/`5e-5 Å` 的 cell/position 两侧测试核验。旧默认接受的
  `1e-5`–`1e-4 Å` 区间现在拒绝，不宣称所有历史输出格式均兼容；未对独立安装的
  外部 abacuslite 或环境 ASE 打补丁。行为验证后更新 checker 中此函数的 AST 登记，
  基线 SHA 保持不变。更完整的原始配对/GPU 兼容验证仍需后续站点计算证据。
- checker 先核对 `_FRAME_PATCH_AST_SHA256` 中登记的 helper 和三个 template 方法的
  AST 身份，再剥离登记补丁、恢复固定上游方法。不再按 helper 名字无条件忽略函数体。
  注释、格式和 docstring 变化不影响身份；可执行代码演进须先复核行为回归和实际 diff，
  再更新对应登记值及归一化。登记 hash 只证明补丁身份，不替代正确性测试或审查。
- latestio MD 晶胞补丁仅在两个精确代码块**同时**匹配时归一化；测试对晶胞
  缩放、PBC、分支条件和无关返回语句的 mutation 均要求 drift checker 拒绝。
