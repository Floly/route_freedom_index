"""Result container with save/load (result.json + pairs.parquet + edges.parquet)."""

import json
import math
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from rfi.utils.log import RFILogger

log = RFILogger(__name__)

_JSON_FIELDS = ("score", "bins", "circuity", "cleaning_stats", "diagnostics", "meta")


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
    bins: list  # per-bin stats: range_m, n_pairs, median_r, iqr_r, median_circuity, used
    circuity: float  # median d* / euclid over all kept pairs
    pairs: pd.DataFrame
    edges: gpd.GeoDataFrame  # projected CRS, with edge_freedom
    cleaning_stats: list
    diagnostics: dict
    meta: dict

    def save(self, out_dir) -> Path:
        out = Path(out_dir)
        with log.step("save", out=str(out)):
            out.mkdir(parents=True, exist_ok=True)
            payload = _plain({k: getattr(self, k) for k in _JSON_FIELDS})
            (out / "result.json").write_text(json.dumps(payload, indent=2))
            self.pairs.to_parquet(out / "pairs.parquet")
            self.edges.to_parquet(out / "edges.parquet")
        return out

    @classmethod
    def load(cls, out_dir) -> "Result":
        out = Path(out_dir)
        data = json.loads((out / "result.json").read_text())
        for b in data["bins"]:  # JSON null -> NaN, list -> tuple, as produced by the pipeline
            b["range_m"] = tuple(b["range_m"])
            for k in ("median_r", "iqr_r", "median_circuity"):
                if b[k] is None:
                    b[k] = math.nan
        for k in ("score", "circuity"):
            if data[k] is None:
                data[k] = math.nan
        return cls(pairs=pd.read_parquet(out / "pairs.parquet"),
                   edges=gpd.read_parquet(out / "edges.parquet"), **data)
