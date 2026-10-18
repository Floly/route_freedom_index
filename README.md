# Route Freedom Index (RFI)

How freely can you wander a street network? RFI scores an area by how much of the nearby street network a pedestrian can actually use while staying within ε (default 10%) of the shortest route.

## Install
```bash
conda create -n rfi -c conda-forge python=3.12 geopandas scipy pyarrow
conda activate rfi
pip install -e ".[dev]"            # core + tests
pip install -e ".[dev,map]"        # + kepler.gl map export
pip install -e ".[dev,map,notebook]"
```

## Quick start
```bash
rfi score --point 41.3910,2.1650 --radius 1500 --out results/eixample/
rfi map   --point 41.3910,2.1650 --radius 1500 --out results/eixample/map.html
```
```python
from rfi import rfi
res = rfi({"point": (41.3910, 2.1650), "radius_m": 1500}, eps=0.1)
print(res.score)
```

## Defining an area
The tool accepts one of four inputs and turns it into a WGS84 polygon, the **interior**.

**Rule of thumb:** for comparing places, use **point + radius** with the same radius everywhere (1.5 km recommended).

### Size guidance
- **Diameter ≥ 3 km** is needed to fill the longest trip bin (1.5–3 km). Smaller areas still work, but the long bin is dropped and reported.
- **Area ≤ ~30 km²** keeps runtime and memory reasonable on a laptop. The tool warns above that (`MAX_AREA_KM2_WARN`).
- The tool downloads a little more than you ask for: the interior plus a buffer of about 900 m at ε = 0.1, so detours near the edge aren't cut off. This is automatic.

### 1. Point + radius (recommended)
**Get it:**
- Google Maps: right-click the spot. The first menu line is `lat, lon`. Click it to copy.
- OpenStreetMap: right-click → "Show address". The URL contains `#map=zoom/lat/lon`.

**Use it:**
```bash
rfi score --point 41.3910,2.1650 --radius 1500
```
```python
rfi(area={"point": (41.3910, 2.1650), "radius_m": 1500})
```
`--point` is **lat,lon**, because that's the order maps copy. The parser converts internally.

### 2. Place name
**Get it:** type a name as you would in OpenStreetMap search, e.g. `"Le Marais, Paris"`.

**Check it before scoring.** Nominatim sometimes returns a whole municipality or the wrong place:
```python
import city2graph as c2g
c2g.get_boundaries("Le Marais, Paris").explore()   # opens a map in Jupyter
```

**Use it:**
```bash
rfi score --place "Le Marais, Paris"
```

### 3. Bounding box
Format: `[min_lon, min_lat, max_lon, max_lat]` in WGS84, i.e. **west, south, east, north**.

> Common mistake: Google Maps gives lat,lon, but a bbox starts with **lon**. If your box ends up in the ocean or Antarctica, the order is swapped.

**Get it:**
- **bboxfinder.com:** draw a rectangle. The "Box" field shows `west,south,east,north`. Copy it.
- **openstreetmap.org → Export → "Manually select a different area":** drag the box. Read the four numbers as left = min_lon, bottom = min_lat, right = max_lon, top = max_lat.

**Use it:**
```bash
rfi score --bbox 2.155,41.385,2.175,41.400
```
```python
rfi(area=[2.155, 41.385, 2.175, 41.400])
```

### 4. Polygon (GeoJSON)
Use this when the area isn't a rectangle or circle, e.g. a neighbourhood bounded by a river or a ring road.

**Get it:** go to **geojson.io**, draw a polygon, then Save → GeoJSON. Exactly one Polygon or MultiPolygon is expected.

**Use it:**
```bash
rfi score --geojson my_area.geojson
```

### Quick sanity check (any input)
```python
from rfi.parser.area import parse_area
poly = parse_area(...)          # same input as rfi()
print(poly.area_km2, poly.diameter_km)
poly.explore()                  # confirm it's where you think it is
```

## Method in one paragraph
Sample pairs of intersections inside the area. For each pair with shortest walking distance d*, a street is **usable** if the best walk through it is at most (1+ε)·d*. That takes two shortest-path trees, with no path enumeration. Every usable street lies inside the ellipse with the two points as foci and focal-sum (1+ε)·d*. The per-pair score is the share of street length inside that ellipse that is usable, between 0 and 1. The area score is the median per distance bin (300–750, 750–1500, 1500–3000 m), averaged over bins. Data: Overture Maps via [city2graph](https://city2graph.net).

## Limitations
- Structural only: ignores land use, safety, sidewalk quality, slope, crossings and traffic.
- Depends on Overture/OSM completeness for footpaths.
- Edge usability means "some walk within budget uses it", not a count of distinct routes.
- Uniform walking speed (distance = time).
- Origins are sampled uniformly over intersections, which slightly oversamples dense junction clusters.
- Pairs whose detour ellipse leaves the downloaded area are dropped. They are mostly circuitous pairs, which can bias the score upward; the drop rate is reported per bin.

## Parser contract (for developers, `rfi/parser/area.py`)
Keep this section in sync with the code.

`parse_area(area) -> Interior` accepts:

| input | type | notes |
|---|---|---|
| place name | `str` | geocoded with `c2g.get_boundaries` |
| point + radius | `{"point": (lat, lon), "radius_m": float}` | circle buffered in local UTM, returned in WGS84 |
| bbox | `list`/`tuple` of 4 floats | `[min_lon, min_lat, max_lon, max_lat]` |
| polygon | `shapely` (Multi)Polygon, `GeoDataFrame`, or path to `.geojson` | GeoDataFrames are reprojected to WGS84 if needed |

It validates:
- lon ∈ [-180, 180] and lat ∈ [-90, 90]
- bbox has min < max on both axes
- **likely swapped order:** if a bbox's "lon" values fall within ±90 while its "lat" values exceed 90, raise an error that suggests swapping
- the polygon is valid and non-empty
- `0 < eps ≤ 0.5`

It warns if the area is > `MAX_AREA_KM2_WARN` (30 km²) or the diameter is < 3 km.

It returns `Interior` with:
- `geometry` (WGS84)
- `utm_crs`
- `area_km2`
- `diameter_km`
- `download_geometry(eps)`: the interior buffered by `½ · D_MAX_M · √((1+eps)² − 1) + BUFFER_MARGIN_M`, computed in `utm_crs` and returned in WGS84

## References
- Sevtsuk, A. & Mekonnen, M. (2012). Urban Network Analysis Toolbox. *International Journal of Geomatics and Spatial Analysis*, 22(2), 287–305.
- Alkhaja, N., Alawadi, K. & Ibrahim, H. M. (2023). Post-pandemic planning: Do we have enough and efficient access to parks? *Frontiers in Built Environment*, 9. doi:10.3389/fbuil.2023.1158430
- Tal, G. & Handy, S. (2012). Measuring nonmotorized accessibility and connectivity in a robust pedestrian network. *Transportation Research Record*, 2299, 48–56. doi:10.3141/2299-06
- Boeing, G. (2017). OSMnx: New methods for acquiring, constructing, analyzing, and visualizing complex street networks. *Computers, Environment and Urban Systems*.
- Boeing, G. (2019). Urban spatial order: street network orientation, configuration, and entropy. *Applied Network Science*, 4(1), 67. doi:10.1007/s41109-019-0189-1
- Sato, Y., Pietrostefani, E., Mahabir, R. & Arribas-Bel, D. (2026). City2Graph: A Python library for Heterogeneous Graph Neural Networks and spatial analysis in urban systems. *Computers, Environment and Urban Systems*, 130, 102492. doi:10.1016/j.compenvurbsys.2026.102492. Citation required when using city2graph in research.
