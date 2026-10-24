import pytest

from atst_tools.calculators import factory
from atst_tools.utils.config_schema import ATSTConfig


class FakeProfile:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class FakeAbacus:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class FakeSocketIO:
    instances = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.closed = 0
        FakeSocketIO.instances.append(self)

    def close(self):
        self.closed += 1


@pytest.fixture
def fakes(monkeypatch):
    FakeSocketIO.instances = []
    monkeypatch.setattr(factory, "ATSTAbacusProfile", FakeProfile)
    monkeypatch.setattr(factory, "Abacus", FakeAbacus)
    monkeypatch.setattr(factory, "AbacusSocketIO", FakeSocketIO)
    monkeypatch.setattr(factory, "_SOCKET_CALCULATORS", [])


def _config(socket, workflow="relax", **calculation):
    return {
        "calculation": {"type": workflow, **calculation},
        "calculator": {
            "name": "abacus",
            "abacus": {
                "mpi": 2,
                "socket": socket,
                "parameters": {"calculation": "scf", "cal_force": 1, "ecutwfc": 50},
            },
        },
    }


def test_socket_disabled_builds_fileio_calculator(fakes):
    calc = factory.CalculatorFactory.get_calculator("abacus", _config(False), directory="run")

    assert isinstance(calc, FakeAbacus)
    assert "socket" not in calc.kwargs


def test_socket_true_builds_socketio_with_unique_unix_names(fakes):
    first = factory.CalculatorFactory.get_calculator("abacus", _config(True), directory="a")
    second = factory.CalculatorFactory.get_calculator("abacus", _config(True), directory="b")

    assert isinstance(first, FakeSocketIO)
    assert first.kwargs["inp"] == {"calculation": "scf", "cal_force": 1}
    assert first.kwargs["ecutwfc"] == 50
    assert "port" not in first.kwargs
    assert first.kwargs["unixsocket"].startswith("atst_")
    assert first.kwargs["unixsocket"] != second.kwargs["unixsocket"]
    assert first.kwargs["profile"].kwargs["command"] == "mpirun -np 2 abacus"


def test_socket_mapping_passes_port_and_timeout(fakes):
    calc = factory.CalculatorFactory.get_calculator(
        "abacus", _config({"port": 31500, "timeout": 30.0}), directory="run"
    )

    assert calc.kwargs["port"] == 31500
    assert calc.kwargs["timeout"] == 30.0
    assert "unixsocket" not in calc.kwargs


def test_socket_mapping_disabled_falls_back_to_fileio(fakes):
    calc = factory.CalculatorFactory.get_calculator(
        "abacus", _config({"enabled": False}), directory="run"
    )

    assert isinstance(calc, FakeAbacus)


def test_close_socket_calculators_closes_each_once(fakes):
    factory.CalculatorFactory.get_calculator("abacus", _config(True), directory="run")

    factory.close_socket_calculators()
    factory.close_socket_calculators()

    assert [calc.closed for calc in FakeSocketIO.instances] == [1]


@pytest.mark.parametrize("workflow", ["neb", "autoneb", "irc", "dmf"])
def test_socket_rejects_multi_calculator_workflows(fakes, workflow):
    with pytest.raises(ValueError, match="calculator.abacus.socket is supported only"):
        factory.CalculatorFactory.get_calculator("abacus", _config(True, workflow), directory="run")


def test_socket_rejects_native_md_driver(fakes):
    with pytest.raises(ValueError, match="driver=ase"):
        factory.CalculatorFactory.get_calculator(
            "abacus", _config(True, "md", driver="abacus_native"), directory="run"
        )


def test_socket_rejects_port_and_unixsocket_together(fakes):
    with pytest.raises(ValueError, match="either port or unixsocket"):
        factory.CalculatorFactory.get_calculator(
            "abacus", _config({"port": 1, "unixsocket": "x"}), directory="run"
        )


def _schema_config(socket, workflow="relax", **calculation):
    calc = {"type": workflow, "init_structure": "init.stru", **calculation}
    return {
        "calculation": calc,
        "calculator": {"name": "abacus", "abacus": {"socket": socket}},
    }


def test_schema_accepts_socket_bool_and_mapping():
    assert ATSTConfig.model_validate(_schema_config(True)).calculator.abacus.socket_enabled()
    mapped = ATSTConfig.model_validate(_schema_config({"unixsocket": "h2au", "timeout": 60}))
    assert mapped.calculator.abacus.socket.unixsocket == "h2au"
    assert not ATSTConfig.model_validate(_schema_config(False)).calculator.abacus.socket_enabled()


def test_schema_rejects_socket_for_unsupported_workflow():
    config = _schema_config(True, "vibration")
    ATSTConfig.model_validate(config)
    config["calculation"] = {"type": "irc", "init_structure": "ts.stru"}
    with pytest.raises(ValueError, match="calculator.abacus.socket is supported only"):
        ATSTConfig.model_validate(config)


def test_schema_rejects_unknown_socket_key():
    with pytest.raises(ValueError):
        ATSTConfig.model_validate(_schema_config({"host": "node1"}))
