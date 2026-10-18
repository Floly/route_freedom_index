import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, box
from typer.testing import CliRunner

from rfi import cli
from rfi import constants as C
from rfi.metric import score_network
from rfi.result import Result
from rfi.viz.config import kepler_config
from rfi.viz.export import build_layers
from tests.graphs import grid

runner = CliRunner()
STATUSES = {"shortest", "usable", "unusable_in_ellipse", "ellipse"}


@pytest.fixture
def result(monkeypatch):
    monkeypatch.setattr(C, "N_ORIGINS", 60)
    monkeypatch.setattr(C, "MIN_PAIRS_PER_BIN", 5)
    net = grid(20, 100.0)
    res = score_network(net, 0.1, 0)
    geom = [LineString(net.xy[[a, b]]) for a, b in zip(net.u, net.v)]
    # a point near Barcelona in UTM 31N, so reprojection to WGS84 is meaningful
    off = np.array([430_000.0, 4_580_000.0])
    geom = [LineString(np.asarray(g.coords) + off) for g in geom]
    edges = gpd.GeoDataFrame(
        {"length_m": net.length, "edge_freedom": res["edge_freedom"], "u": net.u, "v": net.v,
         "ux": net.xy[net.u, 0] + off[0], "uy": net.xy[net.u, 1] + off[1],
         "vx": net.xy[net.v, 0] + off[0], "vy": net.xy[net.v, 1] + off[1]},
        geometry=geom, crs=32631)
    poly = gpd.GeoSeries([box(*edges.total_bounds)], crs=32631).to_crs(4326).iloc[0]
    meta = {"eps": 0.1, "seed": 0, "interior_wkt": poly.wkt, "download_wkt": poly.buffer(0.01).wkt}
    return Result(res["score"], res["bins"], res["circuity"], res["pairs"], edges, [],
                  res["diagnostics"], meta)


def test_network_roundtrip_from_edges(result):
    net = result.network()
    assert len(net.u) == len(result.edges)
    assert net.in_interior.any()
    # node coordinates stored on the edges reproduce edge lengths
    e = result.edges
    d = np.hypot(e["ux"] - e["vx"], e["uy"] - e["vy"])
    np.testing.assert_allclose(d, e["length_m"], rtol=1e-9)


def test_layers(result):
    layers = build_layers(result)
    assert set(layers) == {"area", "streets", "origins", "demo"}
    assert all(g.crs.to_epsg() == 4326 for g in layers.values())
    assert set(layers["area"]["role"]) == {"interior", "download"}
    assert {"length_m", "edge_freedom"} <= set(layers["streets"].columns)
    assert {"median_r", "n_pairs"} <= set(layers["origins"].columns)
    demo = layers["demo"]
    assert set(demo["status"]) <= STATUSES and "shortest" in set(demo["status"])
    assert demo["pair_id"].nunique() <= C.N_DEMO_PAIRS
    assert demo["pair_id"].map(type).eq(str).all()  # numeric ids turn into a kepler time filter
    ell = demo[demo["kind"] == "ellipse"]
    assert len(ell) == demo["pair_id"].nunique()


def test_every_demo_pair_has_shortest_edges(result):
    """Each demo pair shows at least one edge on a shortest path."""
    demo = build_layers(result)["demo"]
    assert (demo[demo["status"] == "shortest"].groupby("pair_id").size() >= 1).all()


def test_config_has_layers_and_pair_filter():
    cfg = kepler_config(41.39, 2.165, 3000.0)["config"]
    ids = [lay["config"]["dataId"] for lay in cfg["visState"]["layers"]]
    assert ids == ["area", "streets", "origins", "demo"]
    (flt,) = cfg["visState"]["filters"]
    assert flt["name"] == ["pair_id"] and flt["dataId"] == ["demo"]
    assert cfg["mapState"]["latitude"] == 41.39


def test_cli_map_reuses_matching_result(result, tmp_path, monkeypatch):
    area = {"point": (41.39, 2.165), "radius_m": 500.0}
    result.meta["area_input"] = {"point": [41.39, 2.165], "radius_m": 500.0}
    result.save(tmp_path)
    calls = []
    monkeypatch.setattr("rfi.viz.export.export_map", lambda res, path: calls.append(path))
    monkeypatch.setattr(cli, "rfi", lambda *a, **k: pytest.fail("rescored a matching result"))
    r = runner.invoke(cli.app, ["map", "--point", "41.39,2.165", "--radius", "500",
                                "--out", str(tmp_path / "map.html")])
    assert r.exit_code == 0, r.output
    assert calls == [tmp_path / "map.html"]
    assert cli._fresh_result(tmp_path, area, 0.2, 0) is None  # other eps -> rescore
