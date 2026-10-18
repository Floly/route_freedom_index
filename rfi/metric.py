"""Route Freedom Index: per-pair usable/reachable ratio and the area score."""

import numpy as np

from rfi import constants as C
from rfi.sampling import Network, distance_matrix, make_pairs, sample_origins
from rfi.utils.log import RFILogger

log = RFILogger(__name__)


def pair_score(net: Network, ds, dt, s: int, t: int, eps: float):
    """U, A and the two edge masks for one pair. ds, dt: distance rows of s and t."""
    budget = (1 + eps) * float(ds[t])
    w = net.length
    ds, dt = ds.astype(np.float64), dt.astype(np.float64)
    detour = np.minimum(ds[net.u] + w + dt[net.v], ds[net.v] + w + dt[net.u])
    usable = detour <= budget + C.ELLIPSE_TOL_M / 2

    inside = (np.linalg.norm(net.xy - net.xy[s], axis=1)
              + np.linalg.norm(net.xy - net.xy[t], axis=1)) <= budget + C.ELLIPSE_TOL_M
    edge_in = inside[net.u] & inside[net.v]
    return w[usable].sum(), w[edge_in].sum(), usable, edge_in


def score_pairs(net: Network, dist, pairs, eps: float):
    """Adds U, A, r to pairs; returns (pairs, usable_count, ellipse_count, diagnostics)."""
    n = len(pairs)
    U, A = np.zeros(n), np.zeros(n)
    usable_cnt = np.zeros(len(net.u), dtype=np.int32)
    ellipse_cnt = np.zeros(len(net.u), dtype=np.int32)
    with log.step("metric", pairs=n):
        for k, (s, t, rs, rt) in enumerate(zip(pairs["s"], pairs["t"], pairs["row_s"],
                                               pairs["row_t"])):
            U[k], A[k], usable, edge_in = pair_score(net, dist[rs], dist[rt], s, t, eps)
            usable_cnt += usable
            ellipse_cnt += edge_in
    pairs = pairs.assign(U=U, A=A, r=np.where(A > 0, U / np.where(A > 0, A, 1), np.nan),
                         u_over_dstar=U / pairs["dstar"].to_numpy())
    n_bad = int((pairs["r"] > 1 + C.R_TOL).sum())
    if n_bad:
        log.warning("r > 1: bug or snapping artefact", pairs=n_bad, max_r=float(pairs["r"].max()))
    return pairs, usable_cnt, ellipse_cnt, {"r_gt_1": n_bad}


def _quartiles(x):
    return np.percentile(x, [25, 50, 75]) if len(x) else (np.nan,) * 3


def area_score(pairs) -> dict:
    """RFI = mean of per-bin median r over bins with >= MIN_PAIRS_PER_BIN pairs.

    `redundancy` is the same aggregation of U / d* (usable street length per metre of shortest
    path). It grows with trip length, so read it per bin; the mean over bins mixes lengths.
    """
    bins = []
    for b, (lo, hi) in enumerate(C.DIST_BINS):
        p = pairs[pairs["bin"] == b]
        q25, q50, q75 = _quartiles(p["r"])
        u25, u50, u75 = _quartiles(p["u_over_dstar"])
        bins.append({"bin": b, "range_m": (lo, hi), "n_pairs": len(p), "median_r": q50,
                     "iqr_r": q75 - q25, "median_u_over_dstar": u50, "iqr_u_over_dstar": u75 - u25,
                     "median_circuity": float(np.median(p["dstar"] / p["euclid"]))
                     if len(p) else np.nan, "used": len(p) >= C.MIN_PAIRS_PER_BIN})
    used = [x for x in bins if x["used"]]
    score = float(np.mean([x["median_r"] for x in used])) if used else float("nan")
    redundancy = float(np.mean([x["median_u_over_dstar"] for x in used])) if used else float("nan")
    circuity = float(np.median(pairs["dstar"] / pairs["euclid"])) if len(pairs) else float("nan")
    return {"score": score, "redundancy": redundancy, "bins": bins, "circuity": circuity}


def score_network(net: Network, eps: float = 0.1, seed: int = 0, ellipse_ok=None) -> dict:
    """Full pipeline on a Network: sample, route, score pairs, aggregate."""
    if not 0 < eps <= 0.5:
        raise ValueError(f"eps must satisfy 0 < eps <= 0.5, got {eps}")
    origins = sample_origins(net, seed)
    dist = distance_matrix(net, origins, eps)
    pairs, diag = make_pairs(net, origins, dist, eps, seed, ellipse_ok)
    pairs, usable_cnt, ellipse_cnt, mdiag = score_pairs(net, dist, pairs, eps)
    res = area_score(pairs)
    for x in res["bins"]:
        if not x["used"]:
            log.warning("bin dropped: too few pairs", bin=x["bin"], n_pairs=x["n_pairs"])
    with np.errstate(invalid="ignore", divide="ignore"):
        edge_freedom = np.where(ellipse_cnt > 0, usable_cnt / ellipse_cnt, np.nan)
    return {**res, "pairs": pairs, "edge_freedom": edge_freedom, "diagnostics": {**diag, **mdiag}}
