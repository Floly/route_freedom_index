"""Area input: place name, bbox, point + radius or polygon -> Interior. See parser/README.md."""

import json
from dataclasses import dataclass
from pathlib import Path

import pyproj
import shapely
from shapely import ops
from shapely.geometry import Point, box, shape
from shapely.geometry.base import BaseGeometry

from rfi import constants as C
from rfi.utils.log import RFILogger

log = RFILogger(__name__)


def validate_eps(eps: float) -> float:
    if not 0 < eps <= 0.5:  # also rejects NaN
        raise ValueError(f"eps must satisfy 0 < eps <= 0.5, got {eps}")
    return eps


def utm_crs_for(lon: float, lat: float) -> pyproj.CRS:
    zone = int((lon + 180) // 6) % 60 + 1
    return pyproj.CRS.from_epsg((32600 if lat >= 0 else 32700) + zone)


def _project(geom, src, dst):
    return ops.transform(pyproj.Transformer.from_crs(src, dst, always_xy=True).transform, geom)


@dataclass(frozen=True)
class Interior:
    polygon: BaseGeometry  # WGS84
    utm_crs: pyproj.CRS
    name: str
    area_km2: float
    diameter_m: float

    def polygon_utm(self) -> BaseGeometry:
        return _project(self.polygon, 4326, self.utm_crs)

    def download_geometry(self, eps: float) -> BaseGeometry:
        """WGS84 polygon to download: the interior plus the derived buffer."""
        buffered = self.polygon_utm().buffer(C.download_buffer_m(validate_eps(eps)))
        return _project(buffered, self.utm_crs, 4326)


def _geocode(name: str) -> BaseGeometry:
    import city2graph as c2g

    return c2g.get_boundaries(name).geometry.iloc[0]


def _check_lonlat(lon, lat, what):
    if not (-180 <= lon <= 180 and -90 <= lat <= 90):
        raise ValueError(f"{what}: lon/lat out of range (lon={lon}, lat={lat})")


def _from_bbox(bbox) -> BaseGeometry:
    if len(bbox) != 4:
        raise ValueError("bbox must be (W, S, E, N)")
    w, s, e, n = (float(x) for x in bbox)
    if (abs(s) > 90 or abs(n) > 90) and abs(w) <= 90 and abs(e) <= 90:
        raise ValueError("bbox looks like lat/lon swapped: expected (W, S, E, N) = (lon, lat, ...)")
    _check_lonlat(w, s, "bbox")
    _check_lonlat(e, n, "bbox")
    if not w < e:
        raise ValueError(f"bbox needs W < E, got W={w}, E={e}")
    if not s < n:
        raise ValueError(f"bbox needs S < N, got S={s}, N={n}")
    return box(w, s, e, n)


def _from_point(spec: dict) -> BaseGeometry:
    if "point" not in spec or "radius_m" not in spec:
        raise ValueError('point input needs {"point": (lat, lon), "radius_m": meters}')
    lat, lon = (float(x) for x in spec["point"])
    radius = float(spec["radius_m"])
    if abs(lat) > 90 and abs(lon) <= 90:
        raise ValueError("point looks like lat/lon swapped: expected (lat, lon)")
    _check_lonlat(lon, lat, "point")
    if not radius > 0:
        raise ValueError(f"radius_m must be > 0, got {radius}")
    utm = utm_crs_for(lon, lat)
    circle = _project(Point(lon, lat), 4326, utm).buffer(radius, quad_segs=32)
    return _project(circle, utm, 4326)


def _from_geojson(obj: dict) -> BaseGeometry:
    kind = obj.get("type")
    if kind == "FeatureCollection":
        geoms = [shape(f["geometry"]) for f in obj.get("features", [])]
        return shapely.union_all(geoms) if geoms else shapely.Polygon()
    if kind == "Feature":
        return shape(obj["geometry"])
    return shape(obj)


def _polygonal(geom: BaseGeometry) -> BaseGeometry:
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    if geom.geom_type == "GeometryCollection":
        geom = shapely.union_all([g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")])
    if geom.is_empty or geom.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("area must be a non-empty polygon")
    minx, miny, maxx, maxy = geom.bounds
    _check_lonlat(minx, miny, "polygon")
    _check_lonlat(maxx, maxy, "polygon")
    return geom


def parse_area(area) -> Interior:
    """Place name, (W,S,E,N), {"point": (lat,lon), "radius_m": m}, GeoJSON or shapely polygon."""
    name = "area"
    if isinstance(area, Path):
        area = str(area)
    if isinstance(area, str):
        if area.lower().endswith((".geojson", ".json")):
            if not Path(area).is_file():
                raise ValueError(f"GeoJSON file not found: {area}")
            geom, name = _from_geojson(json.loads(Path(area).read_text())), Path(area).stem
        elif area.strip():
            geom, name = _geocode(area.strip()), area.strip()
        else:
            raise ValueError("empty place name")
    elif isinstance(area, BaseGeometry):
        geom = area
    elif isinstance(area, dict):
        if "type" in area:
            geom = _from_geojson(area)
        else:
            geom, name = _from_point(area), "point"
    elif isinstance(area, (tuple, list)):
        geom, name = _from_bbox(area), "bbox"
    else:
        raise TypeError(f"unsupported area type: {type(area).__name__}")

    geom = _polygonal(geom)
    c = geom.centroid
    utm = utm_crs_for(c.x, c.y)
    minx, miny, maxx, maxy = _project(geom, 4326, utm).bounds
    interior = Interior(geom, utm, name, _project(geom, 4326, utm).area / 1e6,
                        max(maxx - minx, maxy - miny))
    if interior.area_km2 > C.MAX_AREA_KM2_WARN:
        log.warning("interior area is large for an 8 GB machine", area_km2=round(interior.area_km2, 1),
                    limit=C.MAX_AREA_KM2_WARN)
    if interior.diameter_m < C.MIN_DIAMETER_WARN_M - 1:
        log.warning("interior diameter is below the longest trip bin; long bins may stay empty",
                    diameter_m=round(interior.diameter_m), min_m=C.MIN_DIAMETER_WARN_M)
    return interior
