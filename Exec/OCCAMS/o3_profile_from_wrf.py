"""Build an nlay-length erf.o3vmr profile matching WRF o3input=2 (CAM ozone climatology).

Layout of ozone.formatted, verified from the data: (month, lat, plev), plev fastest,
top (0.28 hPa) to bottom. Values are volume mixing ratio.
Layer pressures come from the reference wrfinput (P+PB, domain mean per level),
because the ERF deck uses WRF levels. Output is bottom-up, as ERF indexes ilay = k.
"""
import sys
import numpy as np
import netCDF4 as nc

ref = "/pscratch/sd/k/kehinson/oracle/OCCAMS_for_Seahorce/base_drag_202108/"
lat_c = float(sys.argv[1]) if len(sys.argv) > 1 else 38.5    # CEN_LAT of the reference run
day_of_year = int(sys.argv[2]) if len(sys.argv) > 2 else 213  # 2021-08-01

lat = np.loadtxt(ref + "ozone_lat.formatted")
plev = np.loadtxt(ref + "ozone_plev.formatted") * 100.0       # hPa -> Pa
oz = np.loadtxt(ref + "ozone.formatted").reshape(12, lat.size, plev.size)

# Linear in time between mid-month values, as WRF does, wrapping Dec <-> Jan.
mid = np.array([15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5])
m1 = (np.searchsorted(mid, day_of_year) - 1) % 12
m2 = (m1 + 1) % 12
w2 = ((day_of_year - mid[m1]) % 365) / ((mid[m2] - mid[m1]) % 365)
prof_t = (1 - w2) * oz[m1] + w2 * oz[m2]

if not lat[0] <= lat_c <= lat[-1]:
    sys.exit(f"latitude {lat_c} outside the climatology's {lat[0]}..{lat[-1]}")
j = min(max(np.searchsorted(lat, lat_c), 1), lat.size - 1)
wl = (lat_c - lat[j - 1]) / (lat[j] - lat[j - 1])
prof = (1 - wl) * prof_t[j - 1] + wl * prof_t[j]

w = nc.Dataset(ref + "wrfinput_d01")
p_lay = (w.variables["P"][0] + w.variables["PB"][0]).mean(axis=(1, 2))  # bottom-up

o3 = np.exp(np.interp(np.log(p_lay), np.log(plev), np.log(prof)))

print(f"# months {m1+1}/{m2+1} weight {w2:.2f}; lat {lat_c}; nlay {o3.size}; "
      f"p {p_lay[0]/100:.0f}..{p_lay[-1]/100:.0f} hPa; o3 {o3.min():.3g}..{o3.max():.3g}")
print("erf.o3vmr = " + " ".join(f"{v:.3e}" for v in o3))
