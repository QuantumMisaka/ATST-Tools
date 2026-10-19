"""Tests for the FCP v2 reference producer helpers."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import types
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.io import write


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "run_fcp_v2_reference.py"


def _load_producer():
    spec = importlib.util.spec_from_file_location("run_fcp_v2_reference", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reads_last_evaluated_state_from_fcp_log(tmp_path):
    producer = _load_producer()
    log = tmp_path / "tmp-log-FCP.txt"
    log.write_text(
        "loop\tNELECT\tFermi(eV)\tFermishift(eV)\tmu(eV)\tUcal(V)\tconv(V)\t"
        "Ewithoutentropy(eV)\tEwithoutentropy_grand(eV)\tEtoten(eV)\t"
        "Etoten_grand(eV)\tCpersurf(e/V/A^2)\ttime(s)\n"
        "1\t216.5\t-3.2\t-1.0\t-4.2\t0.4\t0.0\t-1000\t-1001\t-999\t-1000\t0.02\t0\n"
        "2\t216.6\t-3.1\t-1.0\t-4.1\t0.4\t0.1\t-999\t-1000\t-998\t-999\t0.03\t1\n"
    )

    row = producer._read_last_fcp_log(log)

    assert row["iterations"] == 2
    assert row["nelec"] == pytest.approx(216.6)
    assert row["efermi_eV"] == pytest.approx(-3.1)
    assert row["mu_calc_eV"] == pytest.approx(-4.1)
    assert row["residual_mu_eV"] == pytest.approx(0.1)
    assert row["energy_eV"] == pytest.approx(-999.0)
    assert row["capacitance_per_area"] == pytest.approx(0.03)


def test_build_record_maps_fcp_state_and_forces():
    producer = _load_producer()
    atoms = Atoms("H", positions=[[0, 0, 0]], cell=[10, 10, 20], pbc=True)
    config = {
        "fcp": {
            "U": 0.4,
            "NELECT": 216.5,
            "NELECT0": 216.0,
            "FCPconv": 1e-5,
            "work_ref": 4.6,
        },
        "command": "abacus",
        "directory": "/tmp/fcp",
    }
    log = {
        "iterations": 3,
        "nelec": 216.5,
        "efermi_eV": -3.2,
        "fermishift_eV": -1.0,
        "mu_calc_eV": -4.2,
        "residual_mu_eV": 1e-6,
        "raw_energy_eV": -1000.0,
        "energy_eV": -1001.0,
        "capacitance_per_area": 0.02,
    }
    fcp = SimpleNamespace(C=0.4)
    source = Path("/tmp/v2_FCPelectrochem.py")

    record = producer._build_record(
        log=log,
        fcp=fcp,
        atoms=atoms,
        config=config,
        source=source,
        source_sha256="abc",
        energy=-1001.0,
        free_energy=-1000.5,
        forces=np.zeros((1, 3)),
    )

    assert record["nelec"] == pytest.approx(216.5)
    assert record["mu_target_eV"] == pytest.approx(-5.0)
    assert record["vacuum_level_eV"] == pytest.approx(1.0)
    assert record["forces_eV_per_A"] == [[0.0, 0.0, 0.0]]
    assert record["capacitance_e_per_V"] == pytest.approx(0.4)
    assert record["converged"] is True
    assert record["reference_sha256"] == "abc"


def test_requires_reference_fcp_keys():
    producer = _load_producer()
    with pytest.raises(ValueError, match="missing required keys"):
        producer._require_keys({"U": 0.4}, ("U", "NELECT"), "fcp")


def test_run_reference_with_fake_engine_and_fcp_source(tmp_path, monkeypatch):
    """Exercise the producer path without a real ABACUS installation."""
    producer = _load_producer()

    fake_abacus = types.ModuleType("ase.calculators.abacus")

    class FakeAbacusProfile:
        def __init__(self, command):
            self.command = command

    class FakeAbacus:
        def __init__(self, profile=None, directory=".", **kwargs):
            self.profile = profile
            self.directory = directory
            self.parameters = kwargs

    fake_abacus.AbacusProfile = FakeAbacusProfile
    fake_abacus.Abacus = FakeAbacus
    monkeypatch.setitem(sys.modules, "ase.calculators.abacus", fake_abacus)

    fake_source = tmp_path / "v2_FCPelectrochem.py"
    fake_source.write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "import numpy as np",
                "from ase.calculators.calculator import Calculator, all_changes",
                "",
                "class FCP(Calculator):",
                "    implemented_properties = ['energy', 'free_energy', 'forces']",
                "",
                "    def __init__(self, innercalc, fcptxt, U, NELECT, NELECT0,",
                "                 C=1/80, FCPconv=0.01, work_ref=4.6,",
                "                 max_FCP_iter=100, always_adjust=True,",
                "                 adaptive_lr=False, FCPmethod='Newton-fitting', **kwargs):",
                "        super().__init__()",
                "        self.innercalc = innercalc",
                "        self.fcptxt = fcptxt",
                "        self.C = 0.4",
                "        self.fermi = -3.2",
                "        self.fermishift = -1.0",
                "        self.FCPloop = 1",
                "",
                "    def calculate(self, atoms=None, properties=('energy',),",
                "                  system_changes=all_changes):",
                "        super().calculate(atoms, properties, system_changes)",
                "        out = Path(self.innercalc.directory)",
                "        out.mkdir(parents=True, exist_ok=True)",
                "        header = ('loop\\tNELECT\\tFermi(eV)\\tFermishift(eV)\\tmu(eV)\\t'",
                "                  'Ucal(V)\\tconv(V)\\tEwithoutentropy(eV)\\t'",
                "                  'Ewithoutentropy_grand(eV)\\tEtoten(eV)\\t'",
                "                  'Etoten_grand(eV)\\tCpersurf(e/V/A^2)\\ttime(s)\\n')",
                "        row = ('1\\t216.5\\t-3.2\\t-1.0\\t-4.2\\t0.4\\t1e-6\\t-1000\\t-1001\\t'",
                "               '-999\\t-1001\\t0.02\\t0\\n')",
                "        (out / 'tmp-log-FCP.txt').write_text(header + row)",
                "        self.results.update(energy=-1001.0, free_energy=-1000.5,",
                "                            forces=np.zeros((len(atoms), 3)))",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    structure = tmp_path / "input.extxyz"
    write(structure, Atoms("H", positions=[[0, 0, 0]], cell=[10, 10, 20], pbc=True))
    config = {
        "atoms": str(structure),
        "command": "abacus",
        "pseudo_dir": str(tmp_path),
        "basis_dir": str(tmp_path),
        "pseudopotentials": {"H": "H.upf"},
        "basissets": {"H": "H.orb"},
        "parameters": {"calculation": "scf"},
        "directory": str(tmp_path / "run"),
        "fcp": {
            "U": 0.4,
            "NELECT": 216.5,
            "NELECT0": 216.0,
            "FCPconv": 1e-5,
            "work_ref": 4.6,
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    record = producer.run_reference(config_path, fake_source)

    assert record["nelec"] == pytest.approx(216.5)
    assert record["efermi_eV"] == pytest.approx(-3.2)
    assert record["mu_calc_eV"] == pytest.approx(-4.2)
    assert record["mu_target_eV"] == pytest.approx(-5.0)
    assert record["vacuum_level_eV"] == pytest.approx(1.0)
    assert record["energy_eV"] == pytest.approx(-1001.0)
    assert record["free_energy_eV"] == pytest.approx(-1000.5)
    assert record["forces_eV_per_A"] == [[0.0, 0.0, 0.0]]
    assert record["converged"] is True
