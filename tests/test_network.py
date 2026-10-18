import json
import os
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, box

from rfi import constants as C
from rfi.metric import score_network
from rfi.parser import network as nw
from rfi.parser.area import parse_area
from rfi.sampling import from_edges
from tests import graphs


def _segments(rows):
    return gpd.GeoDataFrame(
        pd.DataFrame(rows, columns=["id", "subtype", "class", "subclass"]).assign(
            geometry=[LineString([(0, 0), (1, 1)])] * len(rows)), crs=4326)


def test_walkable_filter_rules():
    segs = _segments([
        ("a", "road", "residential", None),
        ("b", "road", "footway", "sidewalk"),
        ("c", "road", "footway", "crosswalk"),
        ("d", "road", "motorway", None),
        ("e", "road", "trunk", "link"),
        ("f", "road", "primary", "link"),
        ("g", "rail", "standard_gauge", None),
        ("h", "road", "pedestrian", None),
    ])
    out, stats = nw._walkable(segs, drop_sidewalks=True)
    assert sorted(out["id"]) == ["a", "f", "h"]
    assert stats["before"]["residential"] == 1 and stats["after"].get("footway", 0) == 0
    out, _ = nw._walkable(segs, drop_sidewalks=False)
    assert sorted(out["id"]) == ["a", "b", "c", "f", "h"]


def test_ellipse_ok_inside_and_outside():
    xy = np.array([[500.0, 1500.0], [2500.0, 1500.0]])
    ok = nw.make_ellipse_ok(box(0, 0, 3000, 3000), xy)
    s, t = np.array([0, 0, 0]), np.array([1, 1, 1])
    # d = 2000. B=2100: half-width 320; B=2600: 830 (both inside); B=5000 sticks out of the box
    assert ok(s, t, np.array([2100.0, 2600.0, 5000.0])).tolist() == [True, True, False]


def test_ellipse_ok_checks_every_boundary_point(monkeypatch):
    xy = np.array([[500.0, 1500.0], [2500.0, 1500.0]])
    # polygon cuts off the top of the ellipse (y > 1700); the vertex at +b = 320 is outside
    ok = nw.make_ellipse_ok(box(0, 0, 3000, 1700), xy)
    assert not ok(np.array([0]), np.array([1]), np.array([2100.0]))[0]
    monkeypatch.setattr(C, "ELLIPSE_BOUNDARY_PTS", 4)  # vertices at +-a and +-b are enough here
    assert not nw.make_ellipse_ok(box(0, 0, 3000, 1700), xy)(
        np.array([0]), np.array([1]), np.array([2100.0]))[0]


def test_download_is_cached(tmp_path, monkeypatch):
    calls = []

    def fake_load(**kw):
        calls.append(kw)
        seg = gpd.GeoDataFrame({"id": ["s"], "subtype": ["road"], "class": ["residential"],
                                "subclass": [None]}, geometry=[LineString([(0, 0), (1, 1)])], crs=4326)
        con = gpd.GeoDataFrame({"id": ["c"]}, geometry=[Point(0, 0)], crs=4326)
        return {"segment": seg, "connector": con}

    monkeypatch.setattr(nw.c2g, "load_overture_data", fake_load)
    monkeypatch.setattr(C, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(nw, "_latest_release", lambda: "2099-01-01.0")
    poly = box(2.15, 41.38, 2.18, 41.40)
    seg1, _, meta1 = nw._download(poly)
    seg2, _, meta2 = nw._download(poly)
    assert len(calls) == 1 and calls[0]["keep_outer_neighbors"] is True
    assert calls[0]["types"] == ["segment", "connector"] and calls[0]["release"] == "2099-01-01.0"
    assert len(seg1) == len(seg2) == 1 and meta1 == meta2
    assert meta1["overture_release"] == "2099-01-01.0"
    folder = next(tmp_path.iterdir())
    assert (folder / "segments.parquet").exists() and (folder / "connectors.parquet").exists()
    nw._download(box(2.15, 41.38, 2.19, 41.40))  # different polygon, different cache entry
    assert len(calls) == 2


def _mini_city(n=5, step_m=120.0):
    """Tiny Overture-like street grid near Barcelona, WGS84."""
    lon0, lat0 = 2.165, 41.391
    dlon, dlat = step_m / 83000.0, step_m / 111000.0
    cid = {(i, j): f"c{i}_{j}" for i in range(n) for j in range(n)}
    xy = {k: (lon0 + k[0] * dlon, lat0 + k[1] * dlat) for k in cid}
    con = gpd.GeoDataFrame({"id": list(cid.values())}, geometry=[Point(xy[k]) for k in cid], crs=4326)
    rows = []
    for (i, j) in cid:
        for di, dj in ((1, 0), (0, 1)):
            if (i + di, j + dj) in cid:
                a, b = (i, j), (i + di, j + dj)
                rows.append({
                    "id": f"s{i}_{j}_{di}{dj}", "subtype": "road", "class": "residential",
                    "subclass": None, "level_rules": None,
                    "connectors": json.dumps([{"connector_id": cid[a], "at": 0.0},
                                              {"connector_id": cid[b], "at": 1.0}]),
                    "geometry": LineString([xy[a], xy[b]])})
    seg = gpd.GeoDataFrame(rows, crs=4326)
    return seg, con


def test_build_network_offline(monkeypatch):
    seg, con = _mini_city()
    monkeypatch.setattr(nw, "_download", lambda poly: (seg, con, {"overture_release": "test"}))
    interior = parse_area({"point": (41.391 + 0.002, 2.165 + 0.003), "radius_m": 300})
    built = nw.build_network(interior, 0.1)
    net = built.net
    assert len(net.xy) == 25 and len(net.u) == 40
    assert built.meta["overture_release"] == "test"
    assert net.in_interior.sum() > 0
    assert built.cleaning_stats[0]["step"] == "download"
    assert built.cleaning_stats[-1]["nodes"] == 25
    assert len(built.edges) == len(net.u) and np.isclose(built.edges["length_m"].sum(), net.length.sum())
    assert callable(built.ellipse_ok)


def test_from_edges_60k_nodes_under_one_second():
    n = 245  # 60k nodes, ~120k edges
    xy, u, v = graphs._grid_arrays(n, 100.0)
    t0 = time.perf_counter()
    net = from_edges(xy, u, v, np.full(len(u), 100.0))
    elapsed = time.perf_counter() - t0
    assert net.csr.nnz == 2 * len(u)
    assert elapsed < 1.0, elapsed


@pytest.mark.skipif(not os.environ.get("RFI_NETWORK_TESTS"), reason="needs internet")
def test_eixample_live(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "CACHE_DIR", str(tmp_path))
    interior = parse_area({"point": (41.3910, 2.1650), "radius_m": 1500})
    built = nw.build_network(interior, 0.1)
    assert len(built.net.xy) > 5000
    res = score_network(built.net, eps=0.1, seed=0, ellipse_ok=built.ellipse_ok)
    assert 0 < res["score"] <= 1
