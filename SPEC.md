# Route Freedom Index (RFI): spec for Claude Code

---

## 0. Instructions for Claude Code

**How to work through this spec**
1. Read the whole file once. Then build in this order:
   - **Part 1 (Foundations)**
   - **Part 2 (MVP, milestones M1–M3)**
   - **Part 3 (Full solution, M4–M5)**, only after the user confirms the MVP
2. Do one milestone at a time. At the end of each:
   - run `pytest` and `ruff check .`
   - run the milestone's demo command
   - report back briefly: what was built, test results, demo output (score, runtime, peak RAM, diagnostics), deviations from the spec
   - **stop and wait for the user**
3. Before writing `rfi/parser/network.py`, check the city2graph API (https://city2graph.net/latest/api/) and the Overture transportation schema. If a function or field differs from this spec, adapt and report it. Don't invent APIs.

**Rules**
- Don't add user-facing parameters. Any new knob goes into `rfi/constants.py`.
- No `print` in library code. Use `RFILogger` (Part 1.6).
- Create only the files listed in the repo layout (Part 1.8). Ask before adding new files or dependencies outside `pyproject.toml`.
- Respect the performance budget (Part 1.7). If it's exceeded, follow its fix order.
- Keep the "Parser contract" section of `README.md` in sync with `rfi/parser/area.py`.

---

## Part 1. Foundations (shared by MVP and full solution)

### 1.1 Goal
A compact Python tool that takes an area and returns **one number**: how much freedom a pedestrian has to vary their route between two points without walking more than ε (default 10%) beyond the shortest path.

This is not a full walkability score. It measures one structural property: route-choice freedom.

Keep the project small:
- three user-facing parameters (`area`, `eps`, `seed`)
- every other setting is a named constant in `rfi/constants.py`

**Target machine: MacBook Air M1 (8 GB RAM).** Prefer simple, memory-light code over clever optimisation. No GPU, no torch.

### 1.2 Core metric

**Notation**
- G = (V, E): undirected walk network. Nodes are intersections, edge weight w = length in meters.
- s, t: sampled origin and destination nodes. d* = shortest-path distance s→t.
- ε: slack, validated to 0 < ε ≤ 0.5. The budget is B = (1+ε)·d*.
- d_s(·), d_t(·): shortest-path distances from s and from t. **Both s and t are sampled nodes**, so both rows already exist in the Dijkstra matrix (MVP step 3). No extra Dijkstra runs per pair.

**Usable edges.** Edge (u,v) with length w is usable if
`min(d_s(u) + w + d_t(v), d_s(v) + w + d_t(u)) ≤ B + ELLIPSE_TOL_M / 2`.
U = total length of usable edges.

**Reachable ellipse.** Every usable point x satisfies euclid(s,x) + euclid(x,t) ≤ B. So all usable edges lie inside the ellipse with foci s, t and focal-sum B.
- A = total length of edges whose **both endpoints** lie inside that ellipse (projected CRS, meters).
- The test is one numpy expression over node coordinates per pair: `|p − s| + |p − t| ≤ B + ELLIPSE_TOL_M`. No spatial index is needed.
- The ellipse tolerance is deliberately looser than the usable tolerance (`ELLIPSE_TOL_M` vs `ELLIPSE_TOL_M / 2`), so float rounding on exact ties (e.g. aligned s, t on a grid) cannot break U ⊆ A. Euclidean sums are computed in float64; float32 is only for storing the Dijkstra matrix.

**Per-pair score.** `r(s,t) = U / A`: of the streets within geometric reach under the budget, the share you can actually use.
- In theory U ⊆ A, so r ≤ 1. **Do not clip.**
- Check `r ≤ 1 + 1e-6`, count violations into `Result.diagnostics["r_gt_1"]`, and log a warning if any occur. A violation means a bug or a snapping artefact.

Why not U / d*: in a uniform grid, U grows roughly with d*² · density · sqrt(ε). U / A cancels trip length and density to first order.

**Area score**
1. Stratify pairs by **Euclidean** s–t distance into `DIST_BINS`.
2. Take the median r per bin.
3. **RFI** = mean of the bin medians over bins with ≥ `MIN_PAIRS_PER_BIN` pairs. Unfilled bins (small areas) are dropped and reported.

Also report, outside the score: pairs per bin, IQR per bin, median circuity d* / euclid.

### 1.3 Parameters

**User-facing (API and CLI)**
| name | default | meaning |
|---|---|---|
| `area` | required | place name, bbox, point + radius, or GeoJSON polygon. See "Defining an area" in `README.md`. |
| `eps` | 0.1 | allowed detour over the shortest route; 0 < ε ≤ 0.5 |
| `seed` | 0 | seed for all random sampling and subsampling |

**Constants (`rfi/constants.py`)**
| name | value | why it exists |
|---|---|---|
| `D_MAX_M` | 3000 | max Euclidean trip length |
| `DIST_BINS` | (300, 750), (750, 1500), (1500, 3000) m | trip-length strata |
| `N_ORIGINS` | 150 | sampled interior nodes. One Dijkstra each gives ~N²/2 candidate pairs. Controls runtime only. |
| `MAX_PAIRS_PER_BIN` | 300 | caps compute per bin |
| `MIN_PAIRS_PER_BIN` | 50 | a bin below this is dropped |
| `CIRCUITY_CAP` | 1.5 | pairs with d* > `CIRCUITY_CAP · euclid` are counted as "too circuitous" and excluded; sets the Dijkstra cutoff |
| `BUFFER_MARGIN_M` | 200 | extra margin on top of the ellipse half-width |
| `BOUNDARY_DROP_WARN` | 0.10 | warn if a bin loses more than this share of pairs to the boundary check |
| `SNAP_M` | 1.0 | merges line endpoints closer than this, so tiny drawing gaps don't create fake dead ends |
| `DROP_SIDEWALKS` | True | removes separately mapped sidewalks and crosswalks to get one edge per street |
| `N_DEMO_PAIRS` | 20 | example trips exported to the map (Part 3) |
| `MAX_AREA_KM2_WARN` | 30 | parser warns above this interior size (M1 Air budget) |
| `ELLIPSE_TOL_M` | 1e-3 | float tolerance on budget tests: ellipse uses `+ ELLIPSE_TOL_M`, usable uses `+ ELLIPSE_TOL_M / 2` |
| `R_TOL` | 1e-6 | tolerance for the `r ≤ 1 + R_TOL` check |

The validation notebook (Part 3) may override constants (e.g. `SNAP_M`, `DROP_SIDEWALKS`) for robustness checks. The API never exposes them.

**Derived values (computed, never configured)**
- **Download buffer** = `½ · D_MAX_M · √((1+ε)² − 1) + BUFFER_MARGIN_M`. This is the half-width of the ellipse for the longest straight trip, plus a margin: ≈ 890 m at ε = 0.1.
- **Dijkstra cutoff** = `(1+ε) · D_MAX_M · CIRCUITY_CAP`, ≈ 4950 m at ε = 0.1. This guarantees exact distances for every pair that isn't excluded as too circuitous.

### 1.4 Packaging: `pyproject.toml`
Create exactly this file at the repo root:

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "rfi"
version = "0.1.0"
description = "Route Freedom Index: how freely you can wander a street network"
readme = "README.md"
requires-python = ">=3.12"  # city2graph supports 3.12–3.14
dependencies = [
    "city2graph",      # base install only; no torch extras
    "geopandas",
    "shapely>=2",
    "pyproj",
    "numpy",
    "pandas",
    "scipy",
    "pyarrow",
    "typer",
    "psutil",
]

[project.optional-dependencies]
map = ["keplergl"]
dev = ["pytest", "ruff"]
notebook = ["jupyterlab", "matplotlib"]

[project.scripts]
rfi = "rfi.cli:app"

[tool.hatch.build.targets.wheel]
packages = ["rfi"]

[tool.ruff]
line-length = 100

[tool.pytest.ini_options]
testpaths = ["tests"]
```

Install with `pip install -e ".[dev]"`. Add `,map` for Part 3 M4 and `,notebook` for M5.

### 1.5 Output format
`rfi(area, eps=0.1, seed=0) -> Result`, a dataclass with:
- `score`
- per-bin stats
- circuity
- `pairs` (DataFrame)
- `edges` (GeoDataFrame with `edge_freedom`)
- `cleaning_stats`
- `diagnostics`
- `meta` (eps, seed, constants snapshot, data release/date, runtime)

`Result.save(out_dir)` writes:
- `result.json`: score, per-bin stats, circuity, cleaning_stats, diagnostics, meta
- `pairs.parquet`
- `edges.parquet` (GeoParquet)

`Result.load(out_dir)` reads them back.

### 1.6 Logging (`rfi/utils/log.py`)
One small logger class, reused everywhere. No `print` in library code.

```python
class RFILogger:
    """Thin wrapper over std logging. One instance per module: log = RFILogger(__name__)."""
    def __init__(self, name: str): ...
    def info(self, msg, **fields): ...      # appends fields as key=value
    def warning(self, msg, **fields): ...
    def debug(self, msg, **fields): ...
    @contextmanager
    def step(self, name: str, **fields):    # logs start, elapsed seconds, peak RSS (MB) on exit
        ...
```
- Configure the root handler once, in `rfi/__init__.py` or the CLI entry point.
- Format: `HH:MM:SS level module | message key=value ...`.
- Level comes from the env var `RFI_LOG_LEVEL` (default `INFO`). The CLI flag `-v` sets `DEBUG`.
- Every pipeline step (parse, download/cache, filter, split, graph, dijkstra, metric, save) runs inside `log.step(...)` and logs its main counts.

### 1.7 Performance budget (MacBook Air M1, 8 GB)
- 1.5 km radius interior, excluding the first download: scoring in **< 2 min**, peak RAM **< 2 GB**
- Dijkstra matrix: `N_ORIGINS × N` float32. For example, 150 × 60k nodes ≈ 36 MB.
- If a run exceeds the budget, the fix order is:
  1. lower `N_ORIGINS`
  2. lower `MAX_PAIRS_PER_BIN`
  3. only then optimise code

### 1.8 Repo layout
Files marked (P3) are created in Part 3. Everything else is created in Part 1 or 2.
```
pyproject.toml
README.md          # usage, defining an area, method, limitations, parser contract, references (provided)
rfi/
  __init__.py      # exposes rfi(), configures logging
  constants.py
  utils/
    log.py         # RFILogger
  parser/
    area.py
    network.py
  sampling.py
  metric.py
  result.py        # Result dataclass, save/load
  baselines.py     # (P3)
  viz/             # (P3)
    config.py
    export.py
  cli.py
cache/             # gitignored
results/           # gitignored
tests/
notebooks/validation.ipynb   # (P3)
```

---

## Part 2. MVP: validate that it works (M1–M3)

**Scope:** score one real area from the command line, backed by unit tests.
**Out of scope:** map, baselines, validation notebook (all in Part 3).

### 2.1 Area input (`rfi/parser/area.py`)
Implement `parse_area(area) -> Interior` exactly as specified in the "Parser contract" section of `README.md`:
- the four input forms
- validation, including detection of a likely lat/lon swap and the eps range
- warnings for areas above `MAX_AREA_KM2_WARN` or diameter below 3 km
- the `Interior` fields, including `download_geometry(eps)`, which uses the derived download buffer

### 2.2 Network (`rfi/parser/network.py`, city2graph + Overture Maps)
city2graph's data loader is Overture Maps. Overture segments come with **connectors**, so splitting at connectors rather than at every geometric crossing respects bridges and tunnels.

**Step 1. Download, cached**
```python
data = c2g.load_overture_data(
    area=download_polygon,           # WGS84 shapely polygon
    types=["segment", "connector"],
    save_to_file=False, return_data=True,
    keep_outer_neighbors=True,
)
```
- Cache to `cache/{area_hash}/{segments,connectors}.parquet` (GeoParquet). `area_hash` is the sha1 of the polygon WKT rounded to 6 decimals.
- Read from the cache on rerun. Deleting the folder refreshes the data.
- Record the Overture release (if the `overturemaps` library exposes it) or else the download date in `Result.meta`.

**Step 2. Walkable filter** (verify field names and values against the Overture transportation schema)
- keep `subtype == "road"`
- drop `class` in {`motorway`, `trunk`} and their link variants
- if `DROP_SIDEWALKS`, drop sidewalk and crosswalk subclasses
- log counts per class before and after

**Step 3. Project** to `Interior.utm_crs`. Everything downstream is in meters.

**Step 4. Split and snap**
```python
segs = c2g.process_overture_segments(
    segs, get_barriers=False, connectors_gdf=connectors, threshold=SNAP_M,
)
```

**Step 5. Build the graph**
```python
nodes, edges = c2g.segments_to_graph(segs, multigraph=True, directed=False)
nodes, edges = c2g.remove_isolated_components((nodes, edges))
edges["length_m"] = edges.geometry.length
```
- For routing, collapse parallel edges to the min `length_m` per unordered pair.
- For U and A, keep all edges.

**Step 6. Routing structure.** Map `node_id` to a contiguous 0..N-1 index (stable order, sorted by `node_id`). Build a symmetric `scipy.sparse.csr_matrix`.

**Step 7. Masks.** Compute `in_interior` per node. Sampling uses interior nodes only.

**Step 8. Stats.** Record node and edge counts and total length after each step in `Result.cleaning_stats`.

### 2.3 Sampling (`rfi/sampling.py`)
- Sample `N_ORIGINS` interior nodes (seeded), uniformly over nodes. This slightly oversamples dense junction clusters; accept it and note it in the README.
- One call: `scipy.sparse.csgraph.dijkstra(csr, indices=sampled, limit=cutoff)`. That gives a `N_ORIGINS × N` matrix; store it as float32.
- Form pairs i<j among the sampled nodes and assign each to a bin by Euclidean distance. Bins are half-open `[lo, hi)`; the last bin also includes `D_MAX_M`. Pairs with Euclidean distance outside `[DIST_BINS[0].lo, D_MAX_M]` (including < 300 m) are counted as `out_of_range`.
- `N_ORIGINS` is capped at the number of interior nodes.
- Drop and count pairs by exclusion reason in `Result.diagnostics`:
  - `unreachable`: d* is inf
  - `too_circuitous`: d* > `CIRCUITY_CAP · euclid`
  - `boundary`: the ellipse isn't inside the download area. Checked through an optional predicate passed by the caller; with no predicate (synthetic graphs) the count is 0.
- Subsample the remaining pairs to `MAX_PAIRS_PER_BIN` per bin (seeded).
- Report drop counts **per bin**. Warn if `boundary` drops exceed `BOUNDARY_DROP_WARN` in any bin, since dropped pairs are mostly circuitous ones and their loss biases RFI upward.

### 2.4 Compute (`rfi/metric.py`)
For each kept pair, vectorized with numpy over all edges and nodes:
- usable mask over edges, then U
- ellipse mask over nodes, then edges with both endpoints inside, then A
- r = U / A, with the violation check from Part 1.2

Accumulate `edge_freedom` counters (usable count and in-ellipse count per edge) in the same loop. `edge_freedom = usable count / in-ellipse count`; it is NaN where the in-ellipse count is 0. A plain Python loop over pairs is fine; at most ~900 pairs.

### 2.5 API and CLI
- `rfi(area, eps=0.1, seed=0) -> Result` in `rfi/__init__.py` (format in Part 1.5)
- CLI (`typer`): `rfi score <area args> [--eps 0.1] [--seed 0] [-v] --out results/<name>/`
- Area args (see "Defining an area" in `README.md`):
  - `--place "Le Marais, Paris"`
  - `--bbox W,S,E,N`
  - `--point LAT,LON --radius M`
  - `--geojson file.geojson`

### 2.6 Unit tests (pytest), asserting orderings rather than absolute values
- Functions read constants as `C.X` at call time (never as default argument values), so tests can override them with `monkeypatch`.
- r(full grid) > r(grid with 30% of edges removed) > r(tree)
- tree: U equals the shortest-path length for every pair, on a tree where every off-path edge has `2w > ε·d*` (long edges, ε = 0.05); the test asserts that condition. (a side spur of length w is usable iff `2w ≤ ε·d*`).
- reverse test: a spur with `2w ≤ ε·d*` MUST enter U (checks the usable formula itself)
- grid plus cul-de-sacs: r decreases as the share of branches grows
- duplicating every edge in parallel changes the score by less than 5%
- r ≤ 1 + 1e-6 on all synthetic graphs
- parser: each input form, swapped-bbox detection, invalid inputs, eps out of range

### 2.7 MVP milestones
1. **M1: Core on synthetic graphs.**
   - Implement `pyproject.toml`, `constants.py`, `utils/log.py`, `metric.py` and `sampling.py`.
   - Test on graphs built in tests (grid, tree, cul-de-sac grid).
   - *Done when* all metric unit tests pass.
2. **M2: Parser.**
   - Implement `parser/area.py` and `parser/network.py` with caching and cleaning stats.
   - *Done when* parser tests pass, and a short Python snippet builds the graph for a 1.5 km point in Eixample (41.3910, 2.1650) and logs its cleaning stats.
3. **M3: End to end + CLI.**
   - Implement `rfi()`, `result.py` (save/load) and `rfi score`.
   - *Done when* `rfi score --point 41.3910,2.1650 --radius 1500 --out results/eixample/` finishes within the performance budget and writes all three output files.

**→ Stop after M3 and report. Continue to Part 3 only when the user confirms.**

---

## Part 3. Full solution (M4–M5)

### 3.1 Interactive map (kepler.gl, no hosting; extra `rfi[map]`)
Produce a **standalone HTML file** with the `keplergl` Python package (`save_to_html`). It opens from disk, but needs internet for kepler's JS and basemap tiles. If the default basemap needs a Mapbox token, use a free style in the config. If that fails, fall back to `GeoDataFrame.explore()`.

**Datasets** (EPSG:4326, coordinates rounded to 6 decimals, geometry simplified to ~1 m):

| name | geometry | key fields |
|---|---|---|
| `area` | interior + download-area polygons | `role`, `rfi`, `eps`, `n_pairs` |
| `streets` | all edges | `length_m`, `class`, `edge_freedom` |
| `origins` | sampled interior nodes | `median_r`, `n_pairs` |
| `demo` | edges and ellipses for `N_DEMO_PAIRS` pairs | `pair_id`, `kind` (edge/ellipse), `status` (shortest/usable/unusable_in_ellipse), `r` |

**Config is generated in code.** `rfi/viz/config.py` builds the kepler config dict programmatically:
- layers and colors: streets by `edge_freedom`, demo edges by `status`, origins by `median_r`
- a filter on `pair_id`
- an initial map view centered on the interior

`rfi/viz/export.py` provides `export_map(result, path)`, which uses it. Restyling by hand in Jupyter is optional and not required for any milestone.

**CLI:** `rfi map <area args> --out map.html` reuses `results/<name>/` if present, otherwise runs scoring first.

### 3.2 Baselines (`rfi/baselines.py`, used by the validation notebook)
1. Intersection density (degree ≥ 3 nodes per km²)
2. Mean circuity (Boeing)
3. Street orientation entropy (Boeing 2019)
4. Link–node ratio
5. Pedshed ratio (Tal & Handy 2012): reachable area within 400 m and 800 m divided by the circle area, median over sampled nodes. Polygon from `c2g.create_isochrone(..., method="buffer", buffer_distance=25)`.
6. Pedestrian Route Directness: median d* / euclid over the same pairs
7. **Sevtsuk Redundancy Index and Redundancy Count**, the main baseline. Bounded:
   - 100 pairs, short and medium bins only
   - at most 10k paths per pair
   - 5 s time limit per pair
   - log how many pairs hit either limit

### 3.3 Validation notebook (`notebooks/validation.ipynb`)
1. **Stability:**
   - 5 seeds (CV < 3%)
   - ε ∈ {0.05, 0.1, 0.2}
   - `SNAP_M` ∈ {1, 5, 10}
   - `DROP_SIDEWALKS` on/off
2. **Face validity:** a 1.5 km radius for each of:
   - Eixample
   - Midtown Manhattan
   - Venice
   - Le Marais
   - Amsterdam canal ring
   - Shimokitazawa
   - Brasília Asa Sul
   - a cul-de-sac suburb
   - a campus
   - a business park

   Produce a ranked bar chart.
3. **Novelty check:** Spearman and partial correlations of RFI with each baseline. The claim holds only if RFI isn't explained by intersection density alone. Report it either way.

### 3.4 Full-solution milestones
4. **M4: Map.**
   - Implement `viz/config.py`, `viz/export.py` and `rfi map`.
   - *Done when* `map.html` opens from disk with all four layers and the `pair_id` filter works.
5. **M5: Validation notebook.**
   - Implement `baselines.py` and the notebook.
   - *Done when* the notebook runs top to bottom and produces the ranked chart and correlation table.

---

## Part 4. README
The single `README.md` at the repo root (provided) contains:
- install
- quick-start commands
- defining an area (user guide for point, place, bbox, polygon)
- one-paragraph method summary
- limitations
- parser contract (for developers)
- references

Keep it updated as milestones land.
