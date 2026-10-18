"""Synthetic graphs for tests, built as rfi.sampling.Network (same type the real parser makes)."""

import numpy as np

from rfi.sampling import from_edges


def _grid_arrays(n, spacing):
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    xy = np.column_stack([ii.ravel(), jj.ravel()]).astype(float) * spacing
    idx = np.arange(n * n).reshape(n, n)
    u = np.concatenate([idx[:, :-1].ravel(), idx[:-1, :].ravel()])
    v = np.concatenate([idx[:, 1:].ravel(), idx[1:, :].ravel()])
    return xy, u, v


def grid(n=30, spacing=100.0):
    xy, u, v = _grid_arrays(n, spacing)
    return from_edges(xy, u, v, np.full(len(u), spacing))


def jittered_grid(n=20, spacing=100.0, jitter=25.0, keep=0.85, seed=0):
    """Irregular graph; edge length >= euclid, as in real networks."""
    rng = np.random.default_rng(seed)
    xy, u, v = _grid_arrays(n, spacing)
    xy = xy + rng.uniform(-jitter, jitter, xy.shape)
    sel = rng.random(len(u)) < keep
    u, v = u[sel], v[sel]
    length = np.linalg.norm(xy[u] - xy[v], axis=1) * rng.uniform(1.0, 1.3, len(u))
    return from_edges(xy, u, v, length)


def drop_edges(net, frac, seed=0):
    rng = np.random.default_rng(seed)
    keep = rng.random(len(net.u)) >= frac
    return from_edges(net.xy, net.u[keep], net.v[keep], net.length[keep], net.in_interior)


def tree(n=30, spacing=100.0, seed=0):
    """Random spanning tree of the grid (randomised Kruskal)."""
    rng = np.random.default_rng(seed)
    xy, u, v = _grid_arrays(n, spacing)
    parent = np.arange(len(xy))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    keep = []
    for e in rng.permutation(len(u)):
        ru, rv = find(u[e]), find(v[e])
        if ru != rv:
            parent[ru] = rv
            keep.append(e)
    keep = np.array(keep)
    return from_edges(xy, u[keep], v[keep], np.full(len(keep), spacing))


def grid_with_culdesacs(share, n=30, spacing=100.0, stub=200.0, seed=0):
    """Grid plus one dead-end stub of length `stub` on a `share` fraction of nodes."""
    rng = np.random.default_rng(seed)
    xy, u, v = _grid_arrays(n, spacing)
    parents = np.flatnonzero(rng.random(len(xy)) < share)
    tips = len(xy) + np.arange(len(parents))
    tip_xy = xy[parents] + stub * np.array([0.6, 0.8])  # off-lattice direction
    return from_edges(
        np.vstack([xy, tip_xy]),
        np.concatenate([u, parents]),
        np.concatenate([v, tips]),
        np.concatenate([np.full(len(u), spacing), np.full(len(parents), stub)]),
    )


def duplicate_edges(net):
    return from_edges(
        net.xy,
        np.concatenate([net.u, net.u]),
        np.concatenate([net.v, net.v]),
        np.concatenate([net.length, net.length]),
        net.in_interior,
    )


def path_with_spur(spur, n=11, spacing=100.0):
    """Straight path of n nodes plus one perpendicular spur of length `spur` at the middle."""
    xy = np.column_stack([np.arange(n) * spacing, np.zeros(n)])
    xy = np.vstack([xy, [(n // 2) * spacing, spur]])
    u = np.concatenate([np.arange(n - 1), [n // 2]])
    v = np.concatenate([np.arange(1, n), [n]])
    length = np.concatenate([np.full(n - 1, spacing), [spur]])
    return from_edges(xy, u, v, length)
