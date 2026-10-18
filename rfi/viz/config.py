"""kepler.gl config built in code (SPEC 3.1): layers, colors, pair_id filter, initial view."""

import math

# Sequential low -> high freedom (light yellow to dark blue)
FREEDOM_RANGE = ["#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"]
# Ordinal domain is sorted alphabetically: ellipse, shortest, unusable_in_ellipse, usable
DEMO_STATUS_COLORS = ["#ffffff", "#e31a1c", "#9e9e9e", "#1f78b4"]
BASEMAP_STYLE = "https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json"  # no token


def _color_range(name, colors, category="Custom"):
    return {"name": name, "type": "sequential", "category": category, "colors": colors}


def _vis(**over):
    base = {"opacity": 0.8, "strokeOpacity": 0.9, "thickness": 1, "strokeColor": None,
            "colorRange": _color_range("Freedom", FREEDOM_RANGE),
            "strokeColorRange": _color_range("Freedom", FREEDOM_RANGE),
            "radius": 4, "sizeRange": [0, 10], "radiusRange": [0, 50], "heightRange": [0, 500],
            "elevationScale": 5, "stroked": True, "filled": False, "enable3d": False,
            "wireframe": False, "fixedHeight": False}
    return {**base, **over}


def _layer(layer_id, data_id, label, visible, vis, channels, color=(180, 180, 180)):
    return {"id": layer_id, "type": "geojson",
            "config": {"dataId": data_id, "label": label, "color": list(color),
                       "columns": {"geojson": "_geojson"}, "isVisible": visible,
                       "visConfig": vis, "hidden": False,
                       "textLabel": [{"field": None, "color": [255, 255, 255], "size": 18,
                                      "offset": [0, 0], "anchor": "start",
                                      "alignment": "center"}]},
            "visualChannels": channels}


def _channels(**kw):
    none = {"colorField": None, "colorScale": "quantile", "strokeColorField": None,
            "strokeColorScale": "quantile", "sizeField": None, "sizeScale": "linear",
            "heightField": None, "heightScale": "linear", "radiusField": None,
            "radiusScale": "linear"}
    return {**none, **kw}


def _zoom(lat: float, extent_m: float) -> float:
    """Web-mercator zoom that fits ~1.6 * extent_m across a ~1000 px wide map."""
    m_per_px = 1.6 * max(extent_m, 1.0) / 1000
    return round(max(8.0, min(17.0, math.log2(156_543 * math.cos(math.radians(lat)) / m_per_px))), 1)


def kepler_config(lat: float, lon: float, extent_m: float) -> dict:
    """Config for datasets named area, streets, origins, demo (see rfi.viz.export)."""
    layers = [
        _layer("area", "area", "Area", True,
               _vis(strokeColor=[255, 255, 255], thickness=1.5, opacity=0.05, strokeOpacity=0.6),
               _channels()),
        _layer("streets", "streets", "Streets (edge freedom)", True,
               _vis(thickness=1.2, strokeOpacity=0.85),
               _channels(strokeColorField={"name": "edge_freedom", "type": "real"},
                         strokeColorScale="quantize")),
        _layer("origins", "origins", "Origins (median r)", False,
               _vis(stroked=False, filled=True, radius=5, opacity=0.9),
               _channels(colorField={"name": "median_r", "type": "real"},
                         colorScale="quantize")),
        _layer("demo", "demo", "Demo pair", True,
               _vis(thickness=2.5, strokeOpacity=1.0, strokeColorRange={
                   "name": "Status", "type": "qualitative", "category": "Custom",
                   "colors": DEMO_STATUS_COLORS}),
               _channels(strokeColorField={"name": "status", "type": "string"},
                         strokeColorScale="ordinal")),
    ]
    pair_filter = {"dataId": ["demo"], "id": "pair_id_filter", "name": ["pair_id"],
                   "type": "multiSelect", "value": ["00"], "enlarged": False, "plotType": None,
                   "animationWindow": "free", "yAxis": None, "speed": 1}
    return {"version": "v1", "config": {
        "visState": {"filters": [pair_filter], "layers": layers,
                     "interactionConfig": {
                         "tooltip": {"fieldsToShow": {
                             "streets": [{"name": "edge_freedom", "format": None},
                                         {"name": "class", "format": None},
                                         {"name": "length_m", "format": None}],
                             "origins": [{"name": "median_r", "format": None},
                                         {"name": "n_pairs", "format": None}],
                             "demo": [{"name": "pair_id", "format": None},
                                      {"name": "status", "format": None},
                                      {"name": "r", "format": None}],
                             "area": [{"name": "role", "format": None},
                                      {"name": "rfi", "format": None},
                                      {"name": "eps", "format": None},
                                      {"name": "n_pairs", "format": None}]},
                             "compareMode": False, "compareType": "absolute", "enabled": True},
                         "brush": {"size": 0.5, "enabled": False},
                         "geocoder": {"enabled": False}, "coordinate": {"enabled": False}},
                     "layerBlending": "normal", "splitMaps": [], "animationConfig": {
                         "currentTime": None, "speed": 1}},
        "mapState": {"bearing": 0, "dragRotate": False, "latitude": lat, "longitude": lon,
                     "pitch": 0, "zoom": _zoom(lat, extent_m), "isSplit": False},
        "mapStyle": {"styleType": "rfi_basemap", "topLayerGroups": {}, "visibleLayerGroups": {},
                     "threeDBuildingColor": [9, 7, 15], "mapStyles": {"rfi_basemap": {
                         "id": "rfi_basemap", "label": "Carto Dark Matter",
                         "url": BASEMAP_STYLE, "icon": "", "layerGroups": []}}}}}
