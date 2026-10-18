import json
import math

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString
from typer.testing import CliRunner

from rfi import cli
from rfi import constants as C
from rfi.metric import score_network
from rfi.result import Result
from tests.graphs import grid

runner = CliRunner()


@pytest.fixture
def small(monkeypatch):
    monkeypatch.setattr(C, "N_ORIGINS", 60)
    monkeypatch.setattr(C, "MIN_PAIRS_PER_BIN", 5)
    return grid(20, 100.0)


def _result(net) -> Result:
    res = score_network(net, 0.1, 0)
    geom = [LineString(net.xy[[a, b]]) for a, b in zip(net.u, net.v)]
    edges = gpd.GeoDataFrame({"length_m": net.length, "edge_freedom": res["edge_freedom"]},
                             geometry=geom, crs=32631)
    return Result(res["score"], res["bins"], res["circuity"], res["pairs"], edges,
                  [{"step": "graph", "nodes": len(net.xy)}], res["diagnostics"],
                  {"eps": 0.1, "seed": 0, "runtime_s": 0.1, "constants": {"DIST_BINS": C.DIST_BINS}})


def test_save_load_roundtrip(small, tmp_path):
    res = _result(small)
    res.save(tmp_path / "out")
    assert {p.name for p in (tmp_path / "out").iterdir()} == {
        "result.json", "pairs.parquet", "edges.parquet"}
    back = Result.load(tmp_path / "out")
    assert back.score == pytest.approx(res.score)
    assert back.pairs.equals(res.pairs)
    assert back.edges.crs == res.edges.crs
    np.testing.assert_allclose(back.edges["edge_freedom"], res.edges["edge_freedom"])
    assert back.bins[0]["range_m"] == res.bins[0]["range_m"]
    assert back.diagnostics["r_gt_1"] == 0


def test_nan_becomes_null_and_returns(small, tmp_path, monkeypatch):
    monkeypatch.setattr(C, "MIN_PAIRS_PER_BIN", 10**6)  # every bin dropped -> score NaN
    res = _result(small)
    assert math.isnan(res.score)
    res.save(tmp_path)
    assert json.loads((tmp_path / "result.json").read_text())["score"] is None
    assert math.isnan(Result.load(tmp_path).score)


def test_cli_area_forms():
    f = cli.area_from_args
    assert f("Le Marais, Paris", None, None, None, None) == "Le Marais, Paris"
    assert f(None, "2.1,41.3,2.2,41.4", None, None, None) == [2.1, 41.3, 2.2, 41.4]
    assert f(None, None, "41.39,2.165", 1500.0, None) == {"point": (41.39, 2.165),
                                                          "radius_m": 1500.0}


@pytest.mark.parametrize("args", [
    ("a", "1,2,3,4", None, None, None),      # two forms
    (None, None, None, None, None),          # none
    (None, None, "41.39,2.165", None, None),  # point without radius
    (None, "1,2,3", None, None, None),       # bad bbox
    ("a", None, None, 100.0, None),          # radius without point
])
def test_cli_rejects_bad_area_args(args):
    with pytest.raises(Exception, match="exactly one|needs|expects|only goes"):
        cli.area_from_args(*args)


def test_cli_score_invokes_pipeline(small, tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "rfi", lambda area, eps, seed: _result(small))
    out = tmp_path / "o"
    r = runner.invoke(cli.app, ["score", "--point", "41.39,2.165", "--radius", "500",
                                "--out", str(out)])
    assert r.exit_code == 0, r.output
    assert "RFI =" in r.output and (out / "result.json").exists()
    bad = runner.invoke(cli.app, ["score", "--out", str(out)])
    assert bad.exit_code != 0
