#!/usr/bin/env python3
"""Compare the grid orientation of the WRF/ERF and ROMS/REMORA grids over the ocean.

Both models give the angle of the grid's x axis counterclockwise from east, so that
u_east = u*cos(a) - v*sin(a): WRF as COSALPHA/SINALPHA in wrfinput, ROMS as `angle`
(radians) in the grid file. A vector the driver passes without rotation is turned by
the difference of the two angles at the same point. This script reports that
difference at every ocean rho point, matched to the nearest WRF mass point.

usage: compare_grid_angles.py wrfinput_d01 oneyear_grd.nc
"""
import sys

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree


def unit_vectors(lon_deg, lat_deg):
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.column_stack([np.cos(lat) * np.cos(lon),
                            np.cos(lat) * np.sin(lon),
                            np.sin(lat)])


def summary(name, deg):
    p = np.percentile(deg, [0, 1, 50, 99, 100])
    print(f"{name:34s} min {p[0]:8.3f}  p1 {p[1]:8.3f}  median {p[2]:8.3f}"
          f"  p99 {p[3]:8.3f}  max {p[4]:8.3f}  (deg)")


def main(wrf_path, roms_path):
    with Dataset(wrf_path) as w:
        attrs = {k: getattr(w, k) for k in ("MAP_PROJ", "TRUELAT1", "TRUELAT2",
                                            "STAND_LON", "CEN_LAT", "CEN_LON", "DX")
                 if hasattr(w, k)}
        wlon = np.asarray(w["XLONG"][0]).ravel()
        wlat = np.asarray(w["XLAT"][0]).ravel()
        walpha = np.arctan2(np.asarray(w["SINALPHA"][0]),
                            np.asarray(w["COSALPHA"][0])).ravel()
    with Dataset(roms_path) as r:
        mask = np.asarray(r["mask_rho"][:]) > 0.5
        rlon = np.asarray(r["lon_rho"][:])[mask]
        rlat = np.asarray(r["lat_rho"][:])[mask]
        rangle = np.asarray(r["angle"][:])[mask]

    print("WRF projection:", attrs)
    print(f"ROMS ocean rho points: {rangle.size}")

    dist, idx = cKDTree(unit_vectors(wlon, wlat)).query(unit_vectors(rlon, rlat))
    # Chord on the unit sphere -> km; points far beyond a WRF cell lie outside WRF.
    dist_km = dist * 6371.0
    dx_km = float(attrs.get("DX", 2000.0)) / 1000.0
    inside = dist_km < dx_km
    print(f"inside the WRF domain (nearest WRF point within DX): "
          f"{inside.sum()} of {inside.size}")

    a_wrf = walpha[idx][inside]
    a_roms = rangle[inside]
    delta = np.angle(np.exp(1j * (a_roms - a_wrf)))  # wrapped to (-pi, pi]

    summary("WRF alpha at ocean points", np.degrees(a_wrf))
    summary("ROMS angle", np.degrees(a_roms))
    summary("ROMS angle - WRF alpha", np.degrees(delta))
    d = np.abs(delta)
    print(f"cross-axis leak of an unrotated vector, sin|delta|: "
          f"median {np.median(np.sin(d)):.3f}, max {np.sin(d.max()):.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
