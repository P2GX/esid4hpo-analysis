#!/usr/bin/env python3
"""
Figure 1, panel 2 - map of ESID4HPO participating countries.

Europe-primary map with a North America inset (lower-left) so the US and
Canada are visible too. Involved countries in BIH navy/teal on light grey.

One-time data download (Natural Earth 50m admin-0 countries):
  curl -sL -o ne_50m_admin_0_countries.geojson \
    https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson

Deps:  pip install geopandas shapely pyproj
Run:   python fig1_map.py
Outputs: _work/figures/map_eu_inset.png / .pdf / .svg  (SVG editable in PowerPoint)
"""

import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
from pyproj import Transformer

HERE = pathlib.Path(__file__).resolve().parent
GEOJSON = HERE / "ne_50m_admin_0_countries.geojson"
OUTDIR = HERE / "_work" / "figures"

# ISO ADM0_A3 codes. Add "GRC" to INVOLVED_EU once Matentzoglu is confirmed.
INVOLVED_EU = ["DEU", "GBR", "NLD", "AUT", "ITA", "CZE"]
INVOLVED_NA = ["USA", "CAN"]

C_INV    = "#4e7e96"   # BIH navy/teal (matches 'new terms')
C_LAND   = "#e7eaec"   # light grey land
C_BORDER = "#ffffff"   # white borders
MUTED    = "#5f7078"

EU_PROJ = "+proj=laea +lat_0=52 +lon_0=12 +datum=WGS84"
NA_PROJ = "+proj=laea +lat_0=50 +lon_0=-96 +datum=WGS84"

# lon/lat windows framing each map
EU_WIN = (-12, 26, 36, 60)     # UK / Iberia .. Poland ; Italy .. North Sea
NA_WIN = (-128, -56, 25, 60)   # continental US + southern/central Canada


def _extent(proj, lon0, lon1, lat0, lat1):
    t = Transformer.from_crs("EPSG:4326", proj, always_xy=True)
    xs, ys = [], []
    for lon in (lon0, lon1):
        for lat in (lat0, lat1):
            x, y = t.transform(lon, lat)
            xs.append(x); ys.append(y)
    return min(xs), max(xs), min(ys), max(ys)


def main():
    if not GEOJSON.exists():
        raise SystemExit(f"{GEOJSON} not found - download it (see module docstring).")
    gdf = gpd.read_file(GEOJSON)
    gdf = gdf[gdf["ADM0_A3"] != "ATA"]
    gdf["hit"] = gdf["ADM0_A3"].isin(INVOLVED_EU + INVOLVED_NA)

    fig = plt.figure(figsize=(7.6, 7.0))

    # main: Europe
    axE = fig.add_axes([0, 0, 1, 1])
    gE = gdf.to_crs(EU_PROJ)
    gE[~gE.hit].plot(ax=axE, color=C_LAND, edgecolor=C_BORDER, linewidth=0.5)
    gE[gE.hit].plot(ax=axE, color=C_INV, edgecolor=C_BORDER, linewidth=0.6)
    xmin, xmax, ymin, ymax = _extent(EU_PROJ, *EU_WIN)
    axE.set_xlim(xmin, xmax); axE.set_ylim(ymin, ymax)
    axE.set_axis_off(); axE.set_aspect("equal")

    # inset: North America (lower-left)
    axN = fig.add_axes([0.015, 0.02, 0.34, 0.30])
    gN = gdf.to_crs(NA_PROJ)
    gN[~gN.hit].plot(ax=axN, color=C_LAND, edgecolor=C_BORDER, linewidth=0.4)
    gN[gN.hit].plot(ax=axN, color=C_INV, edgecolor=C_BORDER, linewidth=0.5)
    nx0, nx1, ny0, ny1 = _extent(NA_PROJ, *NA_WIN)
    axN.set_xlim(nx0, nx1); axN.set_ylim(ny0, ny1)
    axN.set_aspect("equal")
    axN.set_xticks([]); axN.set_yticks([])
    for s in axN.spines.values():
        s.set_edgecolor("#c4c9ce"); s.set_linewidth(1.0)
    axN.set_facecolor("white")
    axN.text(0.5, 1.04, "North America", transform=axN.transAxes,
             ha="center", va="bottom", fontsize=9, color=MUTED)

    fig.patch.set_facecolor("white")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUTDIR / f"map_eu_inset.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    print(f"saved -> {OUTDIR / 'map_eu_inset.(png|pdf|svg)'}")


if __name__ == "__main__":
    main()
