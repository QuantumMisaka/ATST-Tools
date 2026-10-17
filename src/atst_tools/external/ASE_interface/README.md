# ABACUS ASE interface compatibility paths

ATST-Tools uses the independently installed `abacuslite>=1.1.0,<2` package.
The `abacuslite/` Python modules here only forward historical import paths to
that installed package. They contain no calculator or parser implementation.
In particular, `atst_tools.external.ASE_interface.abacuslite.io.generalio`
continues to work for external callers such as Paimon and resolves to the
same module object as `abacuslite.io.generalio`.

`ABACUSLITE_SNAPSHOT.md` and `PATCHES.md` describe the retired vendored
snapshot. They are frozen audit records, not a live fallback or update process.
The old source and drift checker remain in Git history; current abacuslite
fixes belong in `QuantumMisaka/abacuslite`.
