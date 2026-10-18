import logging
import os
import resource
import sys
import time
from contextlib import contextmanager

_FORMAT = "%(asctime)s %(levelname)s %(name)s | %(message)s"


def _peak_rss_mb() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**20 if sys.platform == "darwin" else peak / 2**10  # bytes vs KiB


def _fmt(fields: dict) -> str:
    return "".join(f" {k}={v}" for k, v in fields.items())


def configure_logging(level: str | None = None) -> None:
    """Configure the root handler once. Level: argument, else RFI_LOG_LEVEL, else INFO."""
    root = logging.getLogger()
    if not root.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(_FORMAT, datefmt="%H:%M:%S"))
        root.addHandler(handler)
    root.setLevel((level or os.environ.get("RFI_LOG_LEVEL", "INFO")).upper())


class RFILogger:
    """Thin wrapper over std logging. One instance per module: log = RFILogger(__name__)."""

    def __init__(self, name: str):
        self._log = logging.getLogger(name)

    def info(self, msg, **fields):
        self._log.info("%s%s", msg, _fmt(fields))

    def warning(self, msg, **fields):
        self._log.warning("%s%s", msg, _fmt(fields))

    def debug(self, msg, **fields):
        self._log.debug("%s%s", msg, _fmt(fields))

    @contextmanager
    def step(self, name: str, **fields):
        """Logs start, then elapsed seconds and peak RSS (MB) on exit."""
        self.info(f"{name} start", **fields)
        t0 = time.perf_counter()
        try:
            yield
        finally:
            self.info(
                f"{name} done", elapsed_s=round(time.perf_counter() - t0, 2),
                peak_rss_mb=round(_peak_rss_mb()),
            )
