# Parser contract

`parse_area(area) -> Interior` turns any supported input into the interior polygon the score is
computed for. Everything downstream (download, graph, sampling) works from `Interior`.

## Accepted inputs

| form | Python value | CLI |
|---|---|---|
| place name | `"Le Marais, Paris"` (geocoded with `city2graph.get_boundaries`, needs internet) | `--place` |
| bbox | `(W, S, E, N)` in degrees, WGS84 | `--bbox W,S,E,N` |
| point + radius | `{"point": (lat, lon), "radius_m": 1500}` | `--point LAT,LON --radius M` |
| polygon | GeoJSON dict (geometry, Feature or FeatureCollection), path to a `.geojson` file, or a shapely `Polygon` / `MultiPolygon` in WGS84 | `--geojson file.geojson` |

A string ending in `.geojson` or `.json` that points to an existing file is read as a polygon;
any other string is a place name.

## Validation (raises `ValueError`)

- bbox: `W < E`, `S < N`, lon within ±180, lat within ±90.
- **Likely lat/lon swap**: a bbox whose S/N exceed ±90 while W/E fit within ±90 (or a point whose
  lat exceeds ±90 while lon fits within ±90) is rejected with a "looks like lat/lon swapped" message.
  Detection is only possible when a value is out of range for its role (e.g. lon > 90 given as lat);
  swapped values that are all within ±90 (Barcelona, Paris, New York) cannot be told apart.
- radius > 0; point lat in ±90, lon in ±180.
- polygon must be non-empty and polygonal; an invalid polygon is repaired with `make_valid`.
- `validate_eps(eps)`: `0 < eps <= 0.5`.

## Warnings (logged, not raised)

- interior area above `MAX_AREA_KM2_WARN` (30 km²)
- interior diameter (larger side of its bounding box, meters) below `MIN_DIAMETER_WARN_M` (3 km):
  the longest trip bins would be empty

## `Interior`

Frozen dataclass:

| field | meaning |
|---|---|
| `polygon` | interior, WGS84 shapely Polygon / MultiPolygon |
| `utm_crs` | `pyproj.CRS`, the UTM zone of the interior centroid; everything downstream is in meters |
| `name` | label for logs and results |
| `area_km2` | interior area, from the UTM projection |
| `diameter_m` | larger side of the interior bounding box in UTM |

Methods:

- `polygon_utm()`: interior projected to `utm_crs`.
- `download_geometry(eps)`: WGS84 polygon to download: the interior buffered by
  `constants.download_buffer_m(eps)` (≈ 890 m at eps = 0.1).
