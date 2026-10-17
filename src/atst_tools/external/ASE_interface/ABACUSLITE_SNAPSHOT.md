# abacuslite vendored 快照基线

- 快照引入：2026-05-10（atst-tools `881a926`）
- drift-check CI 引入：2026-07-06（`7a04854`）钉扎 `762919f6`（= #7588 PR head）
- 当前基线：上游 `70f7ed69b5677c447afdc78e05240e93da660e66`（2026-07-22 commit；2026-07-23 `b8556c7` 写入 CI）
- 同步状态（2026-08-04 核对）：对上游 develop tip 零上游侧欠账；仅 fork 侧语义补丁（efermi 容错 + core.py 帧选择，见 PATCHES.md）
- 差异摘要：以 `scripts/check_abacuslite_snapshot.py --upstream <基线树> --vendored ...` 输出为准
- 本仓补丁更新（2026-10-17）：独立 abacuslite 候选 `dc020430bfa8c76edfc275b0f1eab2f12cc29d4d`
  的 latestio MD_dump 晶胞单位与 PBC 修复精确镜像到 vendored snapshot；补丁登记和
  drift checker 见 PATCHES.md。上述 ABACUS 上游快照基线保持不变。

- 本仓补丁更新（2026-09-16）：SCF 帧身份校验改为 PBC-aware，并传递实际 STRU 文件名；
  drift checker 对登记语义补丁先核验 AST 身份再归一化。该修复尚未上游化，不推进上述
  基线，也不将 2026-08-04 的同步状态外推为当前上游 tip 已核对。

- 本仓补丁更新（2026-10-13）：SCF 几何帧匹配的绝对默认容差由 `1e-4` 收紧至
  `1e-5 Å`，无相对项；经本仓帧选择回归与容差两侧用例确认后更新 AST 身份登记。
  见 `PATCHES.md` 的兼容边界；补丁尚未上游化，钉扎基线 SHA 不变。

## 同步流程

阶段拉取（发布 patch / 支持新 ABACUS / 上游新加固）后：更新本文件基线 SHA + 差异摘要；CI `ABACUS_DEVELOP_REF` 从本文件解析（单一事实源）。
