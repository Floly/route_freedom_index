"""Baseline network measures that RFI is compared against (SPEC 3.2).

All functions take the `Network` of a scored area (`Result.network()`) and, where needed, the
`Result` itself. `all_baselines(res)` returns one row of numbers per area.
"""

import time

import numpy as np
import shapely
from scipy.sparse.csgraph import dijkstra

from rfi import constants as C
from rfi.sampling import Network, sample_origins
from rfi.utils.log import RFILogger

log = RFILogger(__name__)

PEDSHED_RADII_M = (400, 800)
PEDSHED_BUFFER_M = 25
ORIENTATION_BINS = 36
REDUNDANCY_PAIRS = 100
REDUNDANCY_MAX_PATHS = 10_000
REDUNDANCY_MAX_SECONDS = 5.0


def _degree(net: Network) -> np.ndarray:
    return np.diff(net.csr.indptr)  # unique neighbours; parallel edges are collapsed in csr


def intersection_density(net: Network, area_km2: float) -> float:
    """Interior nodes with degree >= 3 per km2."""
    return float(((_degree(net) >= 3) & net.in_interior).sum() / area_km2)


def mean_circuity(net: Network) -> float:
    """Boeing's average circuity: total edge length / total straight-line length of its edges.

    Interior edges only (both endpoints inside).
    """
    keep = net.in_interior[net.u] & net.in_interior[net.v] & (net.u != net.v)
    straight = np.linalg.norm(net.xy[net.u[keep]] - net.xy[net.v[keep]], axis=1)
    ok = straight > 0
    return float(net.length[keep][ok].sum() / straight[ok].sum())


def street_bearings(net: Network) -> np.ndarray:
    """Bearings (degrees clockwise from north) of interior edges, each edge in both directions."""
    keep = net.in_interior[net.u] & net.in_interior[net.v] & (net.u != net.v)
    d = net.xy[net.v[keep]] - net.xy[net.u[keep]]
    bearing = np.degrees(np.arctan2(d[:, 0], d[:, 1])) % 360
    return np.concatenate([bearing, (bearing + 180) % 360])


def orientation_entropy(net: Network) -> float:
    """Shannon entropy (nats) of street bearings in 36 bins of 10 degrees (Boeing 2019).

    Every edge counts in both directions; bins are shifted so N/E/S/W sit at bin centres; edges
    are unweighted; interior edges only. Max is ln(36) = 3.58 (all directions equally common).
    """
    bearing = street_bearings(net)
    width = 360 / ORIENTATION_BINS
    counts = np.bincount(((bearing + width / 2) % 360 // width).astype(int),
                         minlength=ORIENTATION_BINS)
    p = counts[counts > 0] / counts.sum()
    return float(-(p * np.log(p)).sum())


def link_node_ratio(net: Network) -> float:
    """Unique links / nodes, counting only interior nodes and links between interior nodes."""
    inside = net.in_interior
    links = (inside[net.u] & inside[net.v] & (net.u != net.v))
    a, b = np.minimum(net.u[links], net.v[links]), np.maximum(net.u[links], net.v[links])
    n_links = len(np.unique(a.astype(np.int64) * len(net.xy) + b))
    return float(n_links / inside.sum())


def route_directness(pairs) -> float:
    """Pedestrian Route Directness: median d* / euclid over the scored pairs."""
    return float(np.median(pairs["dstar"] / pairs["euclid"]))


def _pedshed_area(geoms, usable_edges):
    """Area of the 25 m buffer around the given (fully reachable) edges."""
    if not usable_edges.any():
        return 0.0
    return shapely.union_all(shapely.buffer(geoms[usable_edges], PEDSHED_BUFFER_M)).area


def pedshed_ratio(net: Network, edge_geoms, seed: int = 0) -> dict:
    """Tal & Handy (2012): reachable area within 400 / 800 m over the circle area.

    Median over the sampled interior origins. Reachable = edges with both ends within the
    radius, buffered by 25 m; this reproduces `c2g.create_isochrone(method="buffer",
    buffer_distance=25)` to within 0.1% (checked on Eixample), but runs one Dijkstra row per
    origin instead of building a graph per origin (c2g: ~2 s per origin).
    `edge_geoms` must be aligned with net.u / net.v.
    """
    origins = sample_origins(net, seed)
    dist = dijkstra(net.csr, directed=False, indices=origins, limit=max(PEDSHED_RADII_M))
    out = {}
    with log.step("pedshed", origins=len(origins)):
        for r in PEDSHED_RADII_M:
            ratios = [_pedshed_area(edge_geoms, (d[net.u] <= r) & (d[net.v] <= r))
                      / (np.pi * r * r) for d in dist]
            out[f"pedshed_{r}"] = float(np.median(ratios))
    return out


def _simple_paths(net_csr, s, t, d_t, budget, max_paths, max_seconds):
    """DFS over simple s -> t paths of length <= budget.

    An edge (u, v, w) extends the path only if len + w + d_t[v] <= budget, so every branch can
    still reach t within budget. Returns (n_paths, sum_of_path_lengths, union_edge_length,
    capped, timed_out).
    """
    indptr, indices, data = net_csr.indptr, net_csr.indices, net_csr.data
    t_end = time.perf_counter() + max_seconds
    n_paths, len_sum, capped, timed_out = 0, 0.0, False, False
    used = {}  # undirected edge (a, b) -> length, for the union
    on_path = {s}
    stack = [(s, indptr[s], 0.0)]
    path, path_w = [s], []
    steps = 0
    while stack:
        steps += 1
        if steps % 2048 == 0 and time.perf_counter() > t_end:
            timed_out = True
            break
        u, k, length = stack.pop()
        if k >= indptr[u + 1]:
            on_path.discard(u)
            path.pop()
            if path_w:
                path_w.pop()
            continue
        stack.append((u, k + 1, length))
        v, w = indices[k], data[k]
        if v in on_path or length + w + d_t[v] > budget + C.ELLIPSE_TOL_M / 2:
            continue
        if v == t:
            n_paths += 1
            len_sum += length + w
            nodes = path + [v]
            for a, b, w_ab in zip(nodes[:-1], nodes[1:], path_w + [w]):
                used[(min(a, b), max(a, b))] = w_ab
            if n_paths >= max_paths:
                capped = True
                break
            continue
        on_path.add(v)
        path.append(v)
        path_w.append(w)
        stack.append((v, indptr[v], length + w))
    return n_paths, len_sum, sum(used.values()), capped, timed_out


def redundancy(net: Network, pairs, eps: float, seed: int = 0, n_pairs: int = REDUNDANCY_PAIRS):
    """Operational Sevtsuk-style redundancy of sampled pairs (short and medium bins only).

    The exact UNA formula is not published; this is our reconstruction from the tool description
    (cityform.mit.edu, UNA Rhino) and Alkhaja et al. (2023). For a pair (s, t) with shortest
    distance d* and budget B = (1 + eps) d*, let P be the simple s -> t paths of length <= B:
      RC       = |P|                              Redundancy Count
      RI_sum   = sum of len(p) over P / d*        literal reading of Alkhaja et al.
      RI_union = len(union of edges of P) / d*    unique segments; probable UNA reading
    Limits: n_pairs pairs, <= REDUNDANCY_MAX_PATHS paths and <= REDUNDANCY_MAX_SECONDS per pair.
    Pairs that hit a limit get capped=True (counts are then lower bounds).
    Returns a DataFrame with one row per pair.
    """
    import pandas as pd

    pool = pairs[pairs["bin"].isin([0, 1])]
    if pool.empty:
        raise ValueError("no pairs in the short and medium bins")
    rng = np.random.default_rng(seed)
    pick = pool.iloc[np.sort(rng.choice(len(pool), size=min(n_pairs, len(pool)), replace=False))]
    rows = []
    with log.step("redundancy", pairs=len(pick)):
        for _, p in pick.iterrows():
            s, t = int(p["s"]), int(p["t"])
            d_t = dijkstra(net.csr, directed=False, indices=t, limit=C.dijkstra_cutoff_m(eps))
            dstar = float(d_t[s])
            n, len_sum, union, capped, timed_out = _simple_paths(
                net.csr, s, t, d_t, (1 + eps) * dstar, REDUNDANCY_MAX_PATHS,
                REDUNDANCY_MAX_SECONDS)
            rows.append({"s": s, "t": t, "bin": int(p["bin"]), "dstar": dstar, "U": p["U"],
                         "A": p["A"], "r": p["r"], "rc": n, "ri_sum": len_sum / dstar,
                         "ri_union": union / dstar, "union_len": union,
                         "capped": capped or timed_out, "timed_out": timed_out})
    out = pd.DataFrame(rows)
    share = float(out["capped"].mean()) if len(out) else 0.0
    log.info("redundancy done", pairs=len(out), capped_share=round(share, 3))
    return out


def all_baselines(res, seed: int = 0):
    """(row of baselines for a scored area, per-pair redundancy DataFrame).

    Redundancy is summarised in the row by medians, with and without capped pairs.
    """
    net = res.network()
    geoms = res.edges.geometry.to_numpy()
    red = redundancy(net, res.pairs, res.meta["eps"], seed)
    row = {
        "intersection_density": intersection_density(net, res.meta["area_km2"]),
        "mean_circuity": mean_circuity(net),
        "orientation_entropy": orientation_entropy(net),
        "link_node_ratio": link_node_ratio(net),
        **pedshed_ratio(net, geoms, seed),
        "route_directness": route_directness(res.pairs),
        "rc_median": float(red["rc"].median()),
        "ri_sum_median": float(red["ri_sum"].median()),
        "ri_union_median": float(red["ri_union"].median()),
        "rc_median_uncapped": float(red.loc[~red["capped"], "rc"].median()),
        "ri_union_median_uncapped": float(red.loc[~red["capped"], "ri_union"].median()),
        "redundancy_capped_share": float(red["capped"].mean()),
    }
    return row, red
