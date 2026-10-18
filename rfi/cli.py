"""Command line: `rfi score <area args> [--eps 0.1] [--seed 0] [-v] --out DIR`."""

from pathlib import Path
from typing import Annotated

import typer

from rfi import rfi
from rfi.utils.log import configure_logging

app = typer.Typer(add_completion=False, help="Route Freedom Index.")


def _floats(text: str, n: int, what: str) -> list[float]:
    try:
        vals = [float(x) for x in text.split(",")]
    except ValueError:
        vals = []
    if len(vals) != n:
        raise typer.BadParameter(f"{what} expects {n} comma-separated numbers, got {text!r}")
    return vals


def area_from_args(place, bbox, point, radius, geojson):
    """Exactly one area form -> the object `parse_area` accepts."""
    given = [name for name, v in (("--place", place), ("--bbox", bbox), ("--point", point),
                                  ("--geojson", geojson)) if v is not None]
    if len(given) != 1:
        raise typer.BadParameter(
            f"give exactly one of --place, --bbox, --point, --geojson (got {given or 'none'})")
    if radius is not None and point is None:
        raise typer.BadParameter("--radius only goes with --point")
    if place is not None:
        return place
    if geojson is not None:
        return geojson
    if bbox is not None:
        return _floats(bbox, 4, "--bbox W,S,E,N")
    if radius is None:
        raise typer.BadParameter("--point needs --radius (meters)")
    lat, lon = _floats(point, 2, "--point LAT,LON")
    return {"point": (lat, lon), "radius_m": radius}


@app.callback()
def main():
    """Route Freedom Index."""


@app.command()
def score(
    out: Annotated[Path, typer.Option("--out", help="Output folder, e.g. results/eixample/")],
    place: Annotated[str | None, typer.Option(help='Place name, e.g. "Le Marais, Paris"')] = None,
    bbox: Annotated[str | None, typer.Option(help="W,S,E,N in WGS84")] = None,
    point: Annotated[str | None, typer.Option(help="LAT,LON (needs --radius)")] = None,
    radius: Annotated[float | None, typer.Option(help="Radius in meters, with --point")] = None,
    geojson: Annotated[Path | None, typer.Option(help="GeoJSON file, one (Multi)Polygon")] = None,
    eps: Annotated[float, typer.Option(help="Allowed detour, 0 < eps <= 0.5")] = 0.1,
    seed: Annotated[int, typer.Option(help="Seed for all random sampling")] = 0,
    verbose: Annotated[bool, typer.Option("-v", "--verbose", help="DEBUG logging")] = False,
):
    """Score an area and write result.json, pairs.parquet, edges.parquet to --out."""
    if verbose:
        configure_logging("DEBUG")
    area = area_from_args(place, bbox, point, radius, geojson)
    try:
        res = rfi(area, eps=eps, seed=seed)
    except (ValueError, TypeError) as err:
        raise typer.BadParameter(str(err)) from err
    res.save(out)
    typer.echo(f"RFI = {res.score:.4f}  (eps={eps}, seed={seed}, {res.meta['runtime_s']} s)")
    for b in res.bins:
        mark = "" if b["used"] else "  [dropped]"
        typer.echo(f"  bin {b['range_m']}: n={b['n_pairs']} median_r={b['median_r']:.4f}{mark}")
    typer.echo(f"saved to {out}")
