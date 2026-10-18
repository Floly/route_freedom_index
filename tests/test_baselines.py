import math

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString

from rfi import baselines as bl
from tests.graphs import grid, tree


def _pairs(s, t, bin_=0):
    return pd.DataFrame({"s": [s], "t": [t], "bin": [bin_], "U": [0.0], "A": [1.0], "r": [0.0]})


def _geoms(net):
    return np.array([LineString(net.xy[[a, b]]) for a, b in zip(net.u, net.v)], dtype=object)


def test_grid_structure_baselines():
    net = grid(20, 100.0)
    assert bl.intersection_density(net, 1.0) == pytest.approx(396)  # 18**2 degree 4 + 4*18 degree 3
    assert bl.mean_circuity(net) == pytest.approx(1.0)
    assert bl.orientation_entropy(net) == pytest.approx(math.log(4))  # N, E, S, W only
    assert bl.link_node_ratio(net) == pytest.approx(2 * 20 * 19 / 400)


def test_tree_has_lower_link_node_ratio_than_grid():
    assert bl.link_node_ratio(tree(20)) < bl.link_node_ratio(grid(20)) - 0.5


def test_route_directness_is_median_ratio():
    pairs = pd.DataFrame({"dstar": [100.0, 150.0, 400.0], "euclid": [100.0, 100.0, 200.0]})
    assert bl.route_directness(pairs) == pytest.approx(1.5)


def test_pedshed_grid_is_a_share_of_the_circle():
    net = grid(25, 100.0)
    out = bl.pedshed_ratio(net, _geoms(net), seed=0)
    assert set(out) == {"pedshed_400", "pedshed_800"}
    assert all(0.1 < v < 1.5 for v in out.values())


def test_redundancy_on_3x3_grid_counts_the_six_monotone_paths():
    net = grid(3, 100.0)
    s, t = 0, 8  # opposite corners: C(4,2) = 6 shortest paths, all of length 400
    out = bl.redundancy(net, _pairs(s, t), eps=0.05, n_pairs=1).iloc[0]
    assert out["rc"] == 6 and not out["capped"]
    assert out["ri_sum"] == pytest.approx(6.0)  # six paths of length d*
    assert out["ri_union"] == pytest.approx(3.0)  # all 12 edges of 100 m over d* = 400 m


def test_redundancy_counts_a_longer_detour_only_inside_budget():
    net = grid(3, 100.0)
    s, t = 0, 2  # same side: d* = 200, detour around the square = 400
    assert bl.redundancy(net, _pairs(s, t), eps=0.5, n_pairs=1).iloc[0]["rc"] == 1
    assert bl.redundancy(net, _pairs(s, t), eps=1.0 + 1e-6, n_pairs=1).iloc[0]["rc"] > 1


def test_redundancy_cap_is_flagged(monkeypatch):
    monkeypatch.setattr(bl, "REDUNDANCY_MAX_PATHS", 5)
    net = grid(6, 100.0)
    out = bl.redundancy(net, _pairs(0, 35), eps=0.1, n_pairs=1).iloc[0]
    assert out["capped"] and out["rc"] == 5


def test_redundancy_ri_union_tracks_rfi_numerator_on_a_grid():
    """RI_union * d* is the usable length U when every in-budget walk is a simple path."""
    net = grid(5, 100.0)
    out = bl.redundancy(net, _pairs(0, 24), eps=0.05, n_pairs=1).iloc[0]
    assert out["union_len"] == pytest.approx(out["ri_union"] * out["dstar"])
    assert out["union_len"] == pytest.approx(net.length.sum())  # all 40 edges lie on monotone paths
