import json
import logging

import pyproj
import pytest
from shapely import ops
from shapely.geometry import Polygon, box, mapping

from rfi import constants as C
from rfi.parser import area as area_mod
from rfi.parser.area import Interior, parse_area, validate_eps

POINT = {"point": (41.3910, 2.1650), "radius_m": 1500}


def test_point_radius():
    it = parse_area(POINT)
    assert isinstance(it, Interior)
    assert it.utm_crs.to_epsg() == 32631
    assert it.area_km2 == pytest.approx(3.14159 * 1.5**2, rel=0.01)
    assert it.diameter_m == pytest.approx(3000, abs=2)


def test_bbox():
    it = parse_area((2.15, 41.38, 2.18, 41.40))
    assert it.polygon.bounds == pytest.approx((2.15, 41.38, 2.18, 41.40))
    assert 4 < it.area_km2 < 6


def test_geojson_dict_file_and_shapely(tmp_path):
    poly = box(2.15, 41.38, 2.18, 41.40)
    a = parse_area(mapping(poly))
    b = parse_area({"type": "Feature", "properties": {}, "geometry": mapping(poly)})
    c = parse_area({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {}, "geometry": mapping(poly)}]})
    path = tmp_path / "area.geojson"
    path.write_text(json.dumps(mapping(poly)))
    d = parse_area(str(path))
    e = parse_area(poly)
    for it in (b, c, d, e):
        assert it.area_km2 == pytest.approx(a.area_km2)


def test_place_name(monkeypatch):
    monkeypatch.setattr(area_mod, "_geocode", lambda name: box(2.35, 48.85, 2.37, 48.87))
    it = parse_area("Le Marais, Paris")
    assert it.name == "Le Marais, Paris"
    assert it.utm_crs.to_epsg() == 32631


def test_swapped_bbox_detected():
    with pytest.raises(ValueError, match="swapped"):
        parse_area((35.60, 139.60, 35.70, 139.70))  # (S, W, N, E) order, Tokyo
    with pytest.raises(ValueError, match="swapped"):
        parse_area({"point": (120.0, 41.0), "radius_m": 100})  # (lon, lat) order


@pytest.mark.parametrize(
    "bad",
    [
        (2.18, 41.38, 2.15, 41.40),  # W > E
        (2.15, 41.40, 2.18, 41.38),  # S > N
        (2.15, 41.38, 2.18),  # wrong length
        (200.0, 41.38, 201.0, 41.40),  # lon out of range
        {"point": (41.39, 2.165), "radius_m": 0},
        {"point": (41.39, 2.165), "radius_m": -5},
        {"point": (41.39, 2.165)},  # radius missing
        {"type": "Point", "coordinates": [2.1, 41.3]},  # not polygonal
        Polygon(),  # empty
        "/no/such/file.geojson",
        "",
        42,
    ],
)
def test_invalid_inputs(bad):
    with pytest.raises((ValueError, TypeError)):
        parse_area(bad)


@pytest.mark.parametrize("eps", [0, -0.1, 0.51, 1.0, float("nan")])
def test_eps_out_of_range(eps):
    with pytest.raises(ValueError):
        validate_eps(eps)


@pytest.mark.parametrize("eps", [0.01, 0.1, 0.5])
def test_eps_ok(eps):
    assert validate_eps(eps) == eps


def test_download_geometry_buffer():
    it = parse_area(POINT)
    dl = it.download_geometry(0.1)
    expected = C.download_buffer_m(0.1)
    assert expected == pytest.approx(890, abs=5)  # SPEC: ≈ 890 m
    # buffered circle: radius grows by the buffer
    to_utm = pyproj.Transformer.from_crs(4326, it.utm_crs, always_xy=True).transform
    minx, _, maxx, _ = ops.transform(to_utm, dl).bounds
    assert (maxx - minx) / 2 == pytest.approx(1500 + expected, abs=3)
    assert dl.contains(it.polygon)


def test_warnings(caplog):
    with caplog.at_level(logging.WARNING):
        parse_area({"point": (41.39, 2.165), "radius_m": 500})
    assert any("diameter" in r.message for r in caplog.records)
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        parse_area({"point": (41.39, 2.165), "radius_m": 4000})  # ~50 km2
    assert any("area" in r.message for r in caplog.records)
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        parse_area(POINT)
    assert not caplog.records
