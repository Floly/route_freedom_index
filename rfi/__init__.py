import json
import time

from rfi import constants as C
from rfi.metric import score_network
from rfi.parser.area import parse_area, validate_eps
from rfi.parser.network import build_network
from rfi.result import Result
from rfi.utils.log import RFILogger, configure_logging

configure_logging()
log = RFILogger(__name__)

__all__ = ["Result", "rfi"]


def _constants_snapshot() -> dict:
    return {k: getattr(C, k) for k in dir(C) if k.isupper()}


def _area_key(area):
    """JSON-safe description of the area argument, to tell whether a saved result matches it."""
    if isinstance(area, (str, dict, list, tuple)):
        return json.loads(json.dumps(area, default=str))
    return None


def rfi(area, eps: float = 0.1, seed: int = 0) -> Result:
    """Route Freedom Index of an area. See "Defining an area" in README.md for `area`."""
    t0 = time.perf_counter()
    eps = validate_eps(eps)
    with log.step("parse", eps=eps, seed=seed):
        interior = parse_area(area)
    built = build_network(interior, eps)
    res = score_network(built.net, eps, seed, built.ellipse_ok)

    net = built.net
    # node ids and coordinates travel with the edges, so a saved result can rebuild the network
    edges = built.edges.assign(
        edge_freedom=res["edge_freedom"], u=net.u, v=net.v, ux=net.xy[net.u, 0],
        uy=net.xy[net.u, 1], vx=net.xy[net.v, 0], vy=net.xy[net.v, 1])
    meta = {**built.meta, "eps": eps, "seed": seed, "constants": _constants_snapshot(),
            "area_name": interior.name, "area_km2": interior.area_km2,
            "area_input": _area_key(area), "interior_wkt": interior.polygon.wkt,
            "download_wkt": interior.download_geometry(eps).wkt,
            "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("done", score=res["score"], runtime_s=meta["runtime_s"])
    return Result(score=res["score"], bins=res["bins"], circuity=res["circuity"],
                  pairs=res["pairs"], edges=edges, cleaning_stats=built.cleaning_stats,
                  diagnostics=res["diagnostics"], meta=meta)
