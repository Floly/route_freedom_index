"""Result -> four GeoDataFrames (EPSG:4326) -> standalone kepler.gl HTML. See SPEC 3.1."""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from scipy.sparse.csgraph import dijkstra
from shapely.geometry import Point, Polygon

from rfi import constants as C
from rfi.metric import pair_score
from rfi.parser.area import _project
from rfi.result import Result
from rfi.sampling import Network
from rfi.utils.log import RFILogger

log = RFILogger(__name__)

ELLIPSE_PTS = 64
SIMPLIFY_M = 1.0
COORD_PRECISION = 1e-6  # degrees, ~0.1 m


def _wgs84(gdf: gpd.GeoDataFrame, simplify_m=0.0) -> gpd.GeoDataFrame:
    """UTM GeoDataFrame -> WGS84, optionally simplified (in meters) and rounded to 6 decimals."""
    geom = gdf.geometry.simplify(simplify_m) if simplify_m else gdf.geometry
    out = gdf.assign(geometry=geom).to_crs(4326)
    out["geometry"] = shapely.set_precision(out.geometry.values, COORD_PRECISION)
    return out


def _ellipse(p, q, budget: float) -> Polygon:
    """Ellipse with foci p, q (UTM, meters) and focal sum `budget`."""
    c = np.linalg.norm(q - p) / 2
    a = budget / 2
    b = np.sqrt(max(a * a - c * c, 0.0))
    ang = np.arctan2(*(q - p)[::-1])
    th = np.linspace(0, 2 * np.pi, ELLIPSE_PTS, endpoint=False)
    x, y = a * np.cos(th), b * np.sin(th)
    mid = (p + q) / 2
    return Polygon(np.column_stack([mid[0] + x * np.cos(ang) - y * np.sin(ang),
                                    mid[1] + x * np.sin(ang) + y * np.cos(ang)]))


def _area_layer(res: Result) -> gpd.GeoDataFrame:
    n_pairs = int(sum(b["n_pairs"] for b in res.bins))
    geoms = {"interior": shapely.from_wkt(res.meta["interior_wkt"]),
             "download": shapely.from_wkt(res.meta["download_wkt"])}
    return gpd.GeoDataFrame(
        {"role": list(geoms), "rfi": res.score, "eps": res.meta["eps"], "n_pairs": n_pairs},
        geometry=list(geoms.values()), crs=4326)


def _origins_layer(res: Result, net: Network, crs) -> gpd.GeoDataFrame:
    """Sampled nodes that took part in a scored pair, with the median r of their pairs."""
    p = res.pairs
    long = pd.concat([p[["s", "r"]].rename(columns={"s": "node"}),
                      p[["t", "r"]].rename(columns={"t": "node"})])
    g = long.groupby("node")["r"].agg(median_r="median", n_pairs="size").reset_index()
    pts = [Point(*net.xy[n]) for n in g["node"]]
    return gpd.GeoDataFrame(g.drop(columns="node"), geometry=pts, crs=crs)


def _pick_demo_pairs(pairs: pd.DataFrame, seed: int) -> pd.DataFrame:
    """N_DEMO_PAIRS pairs spread evenly over the bins (seeded)."""
    rng = np.random.default_rng(seed)
    per_bin = [g.sample(frac=1, random_state=int(rng.integers(2**31))) for _, g in
               pairs.groupby("bin")]
    picked, k = [], 0
    while len(picked) < min(C.N_DEMO_PAIRS, len(pairs)):
        for g in per_bin:
            if k < len(g) and len(picked) < C.N_DEMO_PAIRS:
                picked.append(g.iloc[k])
        k += 1
    return pd.DataFrame(picked).reset_index(drop=True)


def _demo_layer(res: Result, net: Network, crs) -> gpd.GeoDataFrame:
    """Per demo pair: its ellipse and the in-ellipse edges tagged shortest/usable/unusable."""
    eps = res.meta["eps"]
    edge_geom = res.edges.geometry.to_numpy()
    rows, geoms = [], []
    for pair_id, pr in _pick_demo_pairs(res.pairs, res.meta["seed"]).iterrows():
        s, t = int(pr["s"]), int(pr["t"])
        ds, dt = dijkstra(net.csr, directed=False, indices=[s, t],
                          limit=C.dijkstra_cutoff_m(eps))
        _, _, usable, in_ell = pair_score(net, ds, dt, s, t, eps)
        detour = np.minimum(ds[net.u] + net.length + dt[net.v], ds[net.v] + net.length + dt[net.u])
        shortest = usable & (detour <= ds[t] + C.ELLIPSE_TOL_M / 2)
        status = np.where(shortest, "shortest", np.where(usable, "usable", "unusable_in_ellipse"))
        for i in np.flatnonzero(in_ell):
            rows.append({"pair_id": pair_id, "kind": "edge", "status": status[i], "r": pr["r"]})
            geoms.append(edge_geom[i])
        rows.append({"pair_id": pair_id, "kind": "ellipse", "status": "ellipse", "r": pr["r"]})
        geoms.append(_ellipse(net.xy[s], net.xy[t], (1 + eps) * ds[t]))
    out = gpd.GeoDataFrame(rows, geometry=geoms, crs=crs)
    # strings: kepler reads a numeric pair_id as a timestamp and builds a time slider
    out["pair_id"] = out["pair_id"].map(lambda i: f"{i:02d}")
    return out


def build_layers(res: Result) -> dict[str, gpd.GeoDataFrame]:
    """The four datasets of the map, all EPSG:4326."""
    crs = res.edges.crs
    net = res.network()
    streets = res.edges[["length_m", "edge_freedom", "geometry"]].copy()
    if "class" in res.edges:
        streets["class"] = res.edges["class"].astype(str)
    return {
        "area": _area_layer(res),
        "streets": _wgs84(streets, SIMPLIFY_M),
        "origins": _wgs84(_origins_layer(res, net, crs)),
        "demo": _wgs84(_demo_layer(res, net, crs), SIMPLIFY_M),
    }


def export_map(res: Result, path) -> Path:
    """Write a standalone kepler.gl HTML (needs internet for kepler's JS and the basemap)."""
    from keplergl import KeplerGl  # optional dependency: rfi[map]

    from rfi.viz.config import kepler_config

    path = Path(path)
    with log.step("map", out=str(path)):
        layers = build_layers(res)
        m = KeplerGl(height=800)
        for name, gdf in layers.items():
            m.add_data(data=gdf.copy(), name=name)
        c = shapely.from_wkt(res.meta["interior_wkt"]).centroid
        extent_m = _project(shapely.from_wkt(res.meta["interior_wkt"]), 4326,
                            res.edges.crs).area ** 0.5
        config = kepler_config(c.y, c.x, extent_m)
        path.parent.mkdir(parents=True, exist_ok=True)
        m.save_to_html(config=config, file_name=str(path), read_only=False)
    return path
