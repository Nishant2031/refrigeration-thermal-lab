"""Headless smoke test: every tab renders and every calculate button runs without an exception."""
from pathlib import Path

import pytest

AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
APP = str(Path(__file__).resolve().parents[1] / "app_streamlit.py")


def _button(at, label):
    return next(b for b in at.button if b.label == label)


def test_all_tabs_run_without_exceptions():
    at = AppTest.from_file(APP, default_timeout=180).run()
    assert not at.exception
    for label in ("Run coupled CoolProp simulation", "Calculate required cooling", "Run UA screening",
                  "Check refrigerant energy balance", "Run plate-evaporator screening",
                  "Load bundled Copeland R407C example"):
        _button(at, label).click().run()
        assert not at.exception, label
    _button(at, "Fit map and score held-out rows").click().run()
    assert not at.exception
    assert not at.error, [e.value for e in at.error]
    _button(at, "Run datasheet check").click().run()
    assert not at.exception
    assert float(next(m.value for m in at.metric if m.label == "Mass flow: mean abs. error").split()[0]) < 0.5
    _button(at, "Run sweep").click().run()
    assert not at.exception
    assert any("diminishing returns" in i.value for i in at.info)
    values = {m.label: m.value for m in at.metric}
    assert values["Cooling capacity"] == "4.75 kW"
    assert values["Delivered duty"] in {"4.906 kW", "5.37 kW"}
