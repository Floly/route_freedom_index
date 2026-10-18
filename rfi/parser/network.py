"""Overture Maps -> walk network (city2graph) -> rfi.sampling.Network. See SPEC 2.2."""

import hashlib
import json
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import city2graph as c2g
import geopandas as gpd
import numpy as np
import shapely

from rfi import constants as C
from rfi.parser.area import Interior, _project
from rfi.sampling import Network, from_edges
from rfi.utils.log import RFILogger

log = RFILogger(__name__)

_DROP_CLASSES = ("motorway", "trunk")  # their link variants carry the same class
_DROP_SIDEWALK_SUBCLASSES = ("sidewalk", "crosswalk")


@dataclass
class BuiltNetwork:
    net: Network
    edges: gpd.GeoDataFrame  # UTM, same order as net.u / net.v: geometry, length_m, class
    cleaning_stats: list
    meta: dict
    ellipse_ok: object  # callable(s, t, budget) -> bool array, for sampling.make_pairs


def _latest_release():
    try:
        from overturemaps.core import get_latest_release

        return get_latest_release()
    except Exception as err:  # offline or API change: fall back to the download date
        log.warning("could not read the latest Overture release", error=str(err))
        return None


def _download(polygon):
    """(segments, connectors, meta) for a WGS84 polygon, cached in C.CACHE_DIR/<sha1>/."""
    key = hashlib.sha1(shapely.to_wkt(polygon, rounding_precision=6).encode()).hexdigest()
    folder = Path(C.CACHE_DIR) / key
    files = [folder / n for n in ("segments.parquet", "connectors.parquet", "meta.json")]
    if all(f.exists() for f in files):
        log.info("cache hit", folder=str(folder))
        return (gpd.read_parquet(files[0]), gpd.read_parquet(files[1]),
                json.loads(files[2].read_text()))
    release = _latest_release()
    data = c2g.load_overture_data(
        area=polygon, types=["segment", "connector"], save_to_file=False, return_data=True,
        keep_outer_neighbors=True, release=release,
    )
    seg, con = data["segment"], data["connector"]
    if len(seg) == 0:
        raise ValueError("Overture returned no segments for this area")
    meta = {"overture_release": release, "downloaded": date.today().isoformat()}
    folder.mkdir(parents=True, exist_ok=True)
    seg.to_parquet(files[0])
    con.to_parquet(files[1])
    files[2].write_text(json.dumps(meta))
    return seg, con, meta


def _walkable(segs, drop_sidewalks):
    """Keep roads, drop motorway/trunk (and their links) and, optionally, sidewalks/crosswalks."""
    keep = segs["subtype"] == "road"
    keep &= ~segs["class"].isin(_DROP_CLASSES)
    if drop_sidewalks:
        keep &= ~segs["subclass"].isin(_DROP_SIDEWALK_SUBCLASSES)
    stats = {"before": segs["class"].value_counts().to_dict(),
             "after": segs.loc[keep, "class"].value_counts().to_dict()}
    log.info("walkable filter", before=len(segs), after=int(keep.sum()), classes=stats["after"])
    return segs[keep], stats


def make_ellipse_ok(download_utm, xy):
    """Predicate for sampling.make_pairs: is the pair's ellipse inside the download polygon?

    Tests C.ELLIPSE_BOUNDARY_PTS points on the ellipse boundary against the prepared polygon.
    """
    shapely.prepare(download_utm)
    theta = np.linspace(0, 2 * np.pi, C.ELLIPSE_BOUNDARY_PTS, endpoint=False)
    cos, sin = np.cos(theta), np.sin(theta)

    def ok(s, t, budget):
        ps, pt = xy[s], xy[t]
        mid, d = (ps + pt) / 2, pt - ps
        dist = np.linalg.norm(d, axis=1)
        axis = np.where(dist[:, None] > 0, d / np.maximum(dist, 1e-12)[:, None], [1.0, 0.0])
        normal = np.column_stack([-axis[:, 1], axis[:, 0]])
        a = np.asarray(budget, dtype=float) / 2
        b = np.sqrt(np.maximum(a**2 - (dist / 2) ** 2, 0.0))
        pts = (mid[:, None, :] + (a[:, None] * cos)[:, :, None] * axis[:, None, :]
               + (b[:, None] * sin)[:, :, None] * normal[:, None, :])
        inside = shapely.contains_xy(download_utm, pts[..., 0].ravel(), pts[..., 1].ravel())
        return inside.reshape(len(ps), -1).all(axis=1)

    return ok


def _stat(step, nodes=None, edges=None, length_m=None, **extra):
    return {"step": step, "nodes": nodes, "edges": edges,
            "length_km": None if length_m is None else round(float(length_m) / 1000, 2), **extra}


def build_network(interior: Interior, eps: float) -> BuiltNetwork:
    stats = []
    crs = interior.utm_crs
    download = interior.download_geometry(eps)

    with log.step("download", cache=C.CACHE_DIR):
        seg, con, meta = _download(download)
        stats.append(_stat("download", edges=len(seg), connectors=len(con)))

    with log.step("filter"):
        seg, wstats = _walkable(seg, C.DROP_SIDEWALKS)
        stats.append(_stat("walkable", edges=len(seg), classes=wstats["after"]))

    with log.step("project", crs=str(crs)):
        seg, con = seg.to_crs(crs), con.to_crs(crs)
        stats.append(_stat("project", edges=len(seg), length_m=seg.length.sum()))

    with log.step("split", snap_m=C.SNAP_M):
        seg = c2g.process_overture_segments(
            seg, get_barriers=False, connectors_gdf=con, threshold=C.SNAP_M)
        stats.append(_stat("split", edges=len(seg), length_m=seg.length.sum()))

    with log.step("graph"):
        nodes, edges = c2g.segments_to_graph(seg, multigraph=True, directed=False)
        stats.append(_stat("graph", len(nodes), len(edges), edges.length.sum()))
        nodes, edges = c2g.remove_isolated_components((nodes, edges))
        edges["length_m"] = edges.geometry.length
        stats.append(_stat("isolated_removed", len(nodes), len(edges), edges["length_m"].sum()))

    with log.step("routing"):
        order = nodes.index.sort_values()  # stable contiguous 0..N-1 by node_id
        nodes = nodes.loc[order]
        xy = np.column_stack([nodes.geometry.x, nodes.geometry.y])
        u = order.get_indexer(edges.index.get_level_values(0))
        v = order.get_indexer(edges.index.get_level_values(1))
        if (u < 0).any() or (v < 0).any():
            raise ValueError("edges reference nodes missing from the node table")
        interior_utm = interior.polygon_utm()
        shapely.prepare(interior_utm)
        in_interior = shapely.contains_xy(interior_utm, xy[:, 0], xy[:, 1])
        net = from_edges(xy, u, v, edges["length_m"].to_numpy(), in_interior)
        stats.append(_stat("routing", len(xy), net.csr.nnz // 2, edges["length_m"].sum(),
                           interior_nodes=int(in_interior.sum())))

    cols = [c for c in ("class", "length_m", "geometry") if c in edges.columns]
    edges_out = edges[cols].reset_index(drop=True)
    for s in stats:
        log.info("cleaning", **{k: v for k, v in s.items() if k != "classes"})
    meta = {**meta, "n_nodes": len(xy), "n_edges": len(net.u)}
    ellipse_ok = make_ellipse_ok(_project(download, 4326, crs), xy)
    return BuiltNetwork(net, edges_out, stats, meta, ellipse_ok)
