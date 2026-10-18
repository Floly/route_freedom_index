"""All tunable settings. The API never exposes them; read as `C.X` inside function bodies."""

import math

D_MAX_M = 3000
DIST_BINS = ((300, 750), (750, 1500), (1500, 3000))
N_ORIGINS = 150
MAX_PAIRS_PER_BIN = 300
MIN_PAIRS_PER_BIN = 50
CIRCUITY_CAP = 1.5
BUFFER_MARGIN_M = 200
BOUNDARY_DROP_WARN = 0.10
SNAP_M = 1.0
DROP_SIDEWALKS = True
N_DEMO_PAIRS = 20
MAX_AREA_KM2_WARN = 30
ELLIPSE_TOL_M = 1e-3
R_TOL = 1e-6
ELLIPSE_BOUNDARY_PTS = 32
MIN_DIAMETER_WARN_M = 3000
CACHE_DIR = "cache"


def download_buffer_m(eps: float) -> float:
    """Half-width of the ellipse for the longest straight trip, plus margin."""
    return 0.5 * D_MAX_M * math.sqrt((1 + eps) ** 2 - 1) + BUFFER_MARGIN_M


def dijkstra_cutoff_m(eps: float) -> float:
    """Exact distances for every pair not excluded as too circuitous."""
    return (1 + eps) * D_MAX_M * CIRCUITY_CAP
