import numpy as np
import pytest
from scipy.sparse.csgraph import dijkstra

from rfi import constants as C
from rfi.metric import pair_score, score_network
from rfi.sampling import distance_matrix, make_pairs, sample_origins
from tests import graphs


def _rfi(net, eps=0.1, seed=0):
    return score_network(net, eps=eps, seed=seed)


def _pair_dists(net, s, t):
    d = dijkstra(net.csr, indices=[s, t])
    return d[0], d[1]


def test_ordering_grid_dropped_tree():
    g = graphs.grid()
    r_full = _rfi(g)["score"]
    r_drop = _rfi(graphs.drop_edges(g, 0.3))["score"]
    r_tree = _rfi(graphs.tree())["score"]
    assert r_full > r_drop > r_tree


def test_tree_U_equals_shortest_path():
    net = graphs.tree(n=12, seed=1)
    eps, w = 0.05, 100.0
    rng = np.random.default_rng(0)
    checked = 0
    for s, t in rng.integers(0, len(net.xy), (400, 2)):
        ds, dt = _pair_dists(net, s, t)
        dstar = ds[t]
        if s == t or not 2 * w > eps * dstar:  # precondition: every off-path edge unusable
            continue
        U, _, _, _ = pair_score(net, ds, dt, s, t, eps)
        assert U == pytest.approx(dstar)
        checked += 1
    assert checked >= 20


@pytest.mark.parametrize("spur, in_U", [(40.0, True), (60.0, False)])
def test_spur_usable_iff_detour_within_budget(spur, in_U):
    net = graphs.path_with_spur(spur)
    eps = 0.1
    s, t = 0, 10  # d* = 1000, slack = 100: spur detour is 2*spur
    ds, dt = _pair_dists(net, s, t)
    U, _, usable, _ = pair_score(net, ds, dt, s, t, eps)
    assert usable[-1] == in_U
    assert U == pytest.approx(1000.0 + (spur if in_U else 0.0))


def test_culdesacs_lower_r():
    scores = [_rfi(graphs.grid_with_culdesacs(share))["score"] for share in (0.0, 0.2, 0.5)]
    assert scores[0] > scores[1] > scores[2]


def test_duplicate_edges_changes_score_less_than_5pct():
    g = graphs.drop_edges(graphs.grid(), 0.2)
    a, b = _rfi(g)["score"], _rfi(graphs.duplicate_edges(g))["score"]
    assert abs(a - b) / a < 0.05


@pytest.mark.parametrize(
    "net",
    [
        graphs.grid(),
        graphs.drop_edges(graphs.grid(), 0.3),
        graphs.tree(),
        graphs.grid_with_culdesacs(0.3),
        graphs.jittered_grid(seed=1),
        graphs.duplicate_edges(graphs.grid(n=15)),
    ],
)
def test_r_not_above_one(net):
    res = _rfi(net)
    assert res["diagnostics"]["r_gt_1"] == 0
    assert res["pairs"]["r"].max() <= 1 + C.R_TOL


def test_usable_subset_of_ellipse_random_graph():
    """U ⊆ A edge-wise on irregular graphs, for many pairs and eps values."""
    rng = np.random.default_rng(3)
    for seed in range(3):
        net = graphs.jittered_grid(seed=seed)
        for eps in (0.05, 0.1, 0.5):
            for s, t in rng.integers(0, len(net.xy), (40, 2)):
                ds, dt = _pair_dists(net, s, t)
                if s == t or not np.isfinite(ds[t]):
                    continue
                U, A, usable, edge_in = pair_score(net, ds, dt, s, t, eps)
                assert not (usable & ~edge_in).any()
                assert U <= A * (1 + C.R_TOL)


def test_pair_exclusions_counted(monkeypatch):
    net = graphs.drop_edges(graphs.grid(n=20), 0.6, seed=2)  # many components
    monkeypatch.setattr(C, "N_ORIGINS", 100)
    origins = sample_origins(net, seed=0)
    dist = distance_matrix(net, origins, 0.1)
    assert dist.dtype == np.float32
    pairs, diag = make_pairs(net, origins, dist, 0.1, seed=0)
    assert diag["unreachable"] > 0
    assert len(pairs) <= len(C.DIST_BINS) * C.MAX_PAIRS_PER_BIN
    assert (pairs["euclid"] >= C.DIST_BINS[0][0]).all()


def test_small_bins_dropped(monkeypatch):
    monkeypatch.setattr(C, "MIN_PAIRS_PER_BIN", 10**9)
    res = _rfi(graphs.grid(n=15))
    assert np.isnan(res["score"])
    assert res["bins"] and all(not b["used"] for b in res["bins"])


def test_redundancy_is_u_over_dstar_aggregated_like_the_score():
    res = _rfi(graphs.grid())
    p = res["pairs"]
    np.testing.assert_allclose(p["u_over_dstar"], p["U"] / p["dstar"])
    used = [b for b in res["bins"] if b["used"]]
    assert used
    for b in used:
        assert b["median_u_over_dstar"] == pytest.approx(
            p.loc[p["bin"] == b["bin"], "u_over_dstar"].median())
    assert res["redundancy"] == pytest.approx(np.mean([b["median_u_over_dstar"] for b in used]))


def test_redundancy_orders_grid_above_tree_and_grows_with_trip_length():
    full, tree = _rfi(graphs.grid()), _rfi(graphs.tree())
    assert full["redundancy"] > tree["redundancy"]
    med = [b["median_u_over_dstar"] for b in full["bins"] if b["used"]]
    assert med == sorted(med) and med[-1] > med[0]  # U / d* grows with d*, so read it per bin
