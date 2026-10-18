"""Network container, origin sampling, Dijkstra matrix and pair selection."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra

from rfi import constants as C
from rfi.utils.log import RFILogger

log = RFILogger(__name__)


@dataclass
class Network:
    """Undirected walk network in a projected CRS (meters).

    xy: (N,2) node coordinates. u, v, length: (E,) ALL edges, parallel ones included (used for
    U and A). csr: symmetric routing matrix, parallel edges collapsed to the min length.
    """

    xy: np.ndarray
    u: np.ndarray
    v: np.ndarray
    length: np.ndarray
    in_interior: np.ndarray
    csr: csr_matrix


def from_edges(xy, u, v, length, in_interior=None) -> Network:
    xy = np.asarray(xy, dtype=np.float64)
    u, v = np.asarray(u, dtype=np.int32), np.asarray(v, dtype=np.int32)
    length = np.asarray(length, dtype=np.float64)
    n = len(xy)
    if in_interior is None:
        in_interior = np.ones(n, dtype=bool)

    keep = u != v
    a, b = np.minimum(u[keep], v[keep]).astype(np.int64), np.maximum(u[keep], v[keep])
    key = a * n + b
    uniq, first = np.unique(key, return_index=True)
    order = np.argsort(key, kind="stable")
    starts = np.searchsorted(key[order], uniq)
    w = np.minimum.reduceat(length[keep][order], starts) if len(uniq) else np.empty(0)
    ua, ub = a[first], b[first]
    csr = csr_matrix((np.concatenate([w, w]), (np.concatenate([ua, ub]), np.concatenate([ub, ua]))),
                     shape=(n, n))
    return Network(xy, u, v, length, np.asarray(in_interior, dtype=bool), csr)


def sample_origins(net: Network, seed: int) -> np.ndarray:
    """Seeded uniform sample of interior nodes, at most C.N_ORIGINS, sorted."""
    interior = np.flatnonzero(net.in_interior)
    k = min(C.N_ORIGINS, len(interior))
    return np.sort(np.random.default_rng(seed).choice(interior, size=k, replace=False))


def distance_matrix(net: Network, origins: np.ndarray, eps: float) -> np.ndarray:
    """(len(origins), N) float32 shortest-path distances, cut off at the derived cutoff."""
    with log.step("dijkstra", origins=len(origins), nodes=len(net.xy)):
        d = dijkstra(net.csr, directed=False, indices=origins, limit=C.dijkstra_cutoff_m(eps))
    return d.astype(np.float32)


def make_pairs(net, origins, dist, eps, seed, ellipse_ok=None):
    """Pairs i<j among origins, binned by Euclidean distance, with exclusions counted.

    ellipse_ok(s, t, budget) -> bool array tells whether the ellipse lies inside the download
    area; None means no boundary check (synthetic graphs). Returns (pairs DataFrame, diagnostics).
    """
    k = len(origins)
    i, j = np.triu_indices(k, 1)
    s, t = origins[i], origins[j]
    euclid = np.linalg.norm(net.xy[s] - net.xy[t], axis=1)
    dstar = dist[i, t].astype(np.float64)

    bin_id = np.full(len(i), -1)
    for b, (lo, hi) in enumerate(C.DIST_BINS):
        last = b == len(C.DIST_BINS) - 1
        bin_id[(euclid >= lo) & ((euclid <= hi) if last else (euclid < hi))] = b
    diag = {"out_of_range": int((bin_id < 0).sum()), "unreachable": 0, "too_circuitous": 0,
            "boundary": 0, "bins": []}

    _, comp = connected_components(net.csr, directed=False)
    unreachable = comp[s] != comp[t]
    circuitous = ~unreachable & ~(dstar <= C.CIRCUITY_CAP * euclid)  # inf (past cutoff) included
    ok = (bin_id >= 0) & ~unreachable & ~circuitous
    boundary = np.zeros(len(i), dtype=bool)
    if ellipse_ok is not None and ok.any():
        boundary[ok] = ~ellipse_ok(s[ok], t[ok], (1 + eps) * dstar[ok])
    ok &= ~boundary

    rng = np.random.default_rng(seed)
    keep_idx = []
    for b in range(len(C.DIST_BINS)):
        in_bin = bin_id == b
        row = {"bin": b, "candidates": int(in_bin.sum()),
               "unreachable": int((in_bin & unreachable).sum()),
               "too_circuitous": int((in_bin & circuitous).sum()),
               "boundary": int((in_bin & boundary).sum())}
        idx = np.flatnonzero(in_bin & ok)
        if len(idx) > C.MAX_PAIRS_PER_BIN:
            idx = np.sort(rng.choice(idx, size=C.MAX_PAIRS_PER_BIN, replace=False))
        row["kept"] = len(idx)
        keep_idx.append(idx)
        diag["bins"].append(row)
        reached = row["candidates"] - row["unreachable"] - row["too_circuitous"]
        if reached and row["boundary"] / reached > C.BOUNDARY_DROP_WARN:
            log.warning("boundary check dropped many pairs; RFI biased upward", bin=b,
                        dropped=row["boundary"], of=reached)
    for name in ("unreachable", "too_circuitous", "boundary"):
        diag[name] = sum(r[name] for r in diag["bins"])

    keep = np.concatenate(keep_idx)
    pairs = pd.DataFrame({"s": s[keep], "t": t[keep], "row_s": i[keep], "row_t": j[keep],
                          "bin": bin_id[keep], "euclid": euclid[keep], "dstar": dstar[keep]})
    log.info("pairs", kept=len(pairs), **{k_: diag[k_] for k_ in
             ("out_of_range", "unreachable", "too_circuitous", "boundary")})
    return pairs, diag
