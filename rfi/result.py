"""Result container with save/load (result.json + pairs.parquet + edges.parquet)."""

import json
import math
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from rfi.parser.area import _project
from rfi.sampling import Network, from_edges
from rfi.utils.log import RFILogger

log = RFILogger(__name__)

_JSON_FIELDS = ("score", "redundancy", "bins", "circuity", "cleaning_stats", "diagnostics",
                "meta")


def _plain(x):
    """JSON-safe copy: numpy -> python, tuples -> lists, NaN/inf -> None."""
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


@dataclass
class Result:
    score: float
    bins: list  # per-bin stats: range_m, n_pairs, median_r, iqr_r, median_u_over_dstar, ..., used
    circuity: float  # median d* / euclid over all kept pairs
    pairs: pd.DataFrame
    edges: gpd.GeoDataFrame  # projected CRS, with edge_freedom
    cleaning_stats: list
    diagnostics: dict
    meta: dict
    redundancy: float = math.nan  # U / d*, aggregated like the score (per-bin median, then mean)

    def save(self, out_dir) -> Path:
        out = Path(out_dir)
        with log.step("save", out=str(out)):
            out.mkdir(parents=True, exist_ok=True)
            payload = _plain({k: getattr(self, k) for k in _JSON_FIELDS})
            (out / "result.json").write_text(json.dumps(payload, indent=2))
            self.pairs.to_parquet(out / "pairs.parquet")
            self.edges.to_parquet(out / "edges.parquet")
        return out

    def network(self) -> Network:
        """Routing network rebuilt from the node ids and coordinates stored in `edges`."""
        e = self.edges
        u, v = e["u"].to_numpy(), e["v"].to_numpy()
        xy = np.zeros((int(max(u.max(), v.max())) + 1, 2))
        xy[u] = e[["ux", "uy"]].to_numpy()
        xy[v] = e[["vx", "vy"]].to_numpy()
        interior = _project(shapely.from_wkt(self.meta["interior_wkt"]), 4326, e.crs)
        in_interior = shapely.contains_xy(interior, xy[:, 0], xy[:, 1])
        return from_edges(xy, u, v, e["length_m"].to_numpy(), in_interior)

    @classmethod
    def load(cls, out_dir) -> "Result":
        out = Path(out_dir)
        data = json.loads((out / "result.json").read_text())
        for b in data["bins"]:  # JSON null -> NaN, list -> tuple, as produced by the pipeline
            b["range_m"] = tuple(b["range_m"])
            for k in ("median_r", "iqr_r", "median_circuity", "median_u_over_dstar",
                      "iqr_u_over_dstar"):
                if b.get(k) is None:  # None: NaN saved as null, or a result saved before U/d*
                    b[k] = math.nan
        for k in ("score", "circuity", "redundancy"):
            if data.get(k) is None:
                data[k] = math.nan
        return cls(pairs=pd.read_parquet(out / "pairs.parquet"),
                   edges=gpd.read_parquet(out / "edges.parquet"), **data)
