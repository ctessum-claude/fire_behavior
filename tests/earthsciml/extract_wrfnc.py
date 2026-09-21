#!/usr/bin/env python3
"""Extract one frame of test7's coupled WRF+CFBM output (wrf.nc) into the flat text format read by
feedback_driver.F90 (record: `<name> <kind> <rank> <extents>` then the values, Fortran order).
Fields: fire-mesh fluxes FGRNHFX/FGRNQFX, atmosphere RHO, z_at_w = (PH+PHB)/g, dz8w, MUT, C1H, C2H,
unstaggered u_phy/v_phy, ZNT, FZ0, UF, VF, FXLAT/FXLONG, XLAT/XLONG, and the WRF-computed GRNHFX/GRNQFX,
RTHFRTEN/RQVFRTEN for comparison.  Usage: extract_wrfnc.py wrf.nc frame out.flat"""
import sys, numpy as np
from netCDF4 import Dataset
G = 9.81
fn, t, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
cfbm = sys.argv[4] if len(sys.argv) > 4 else None   # optional esm_dump_fire_flux_<n>.json: use CFBM's own fire-mesh fluxes
d = Dataset(fn)
def v(name): return np.asarray(d.variables[name][t], dtype=np.float64)
def tr3(a):  # netCDF (k, j, i) -> Fortran (i, k, j)
    return np.transpose(a, (2, 0, 1))
def tr2(a):  # (j, i) -> (i, j)
    return np.transpose(a, (1, 0))
ph = tr3(v('PH') + v('PHB')) / G            # (i, kstag, j) 45 levels
zw = ph
dz8w = zw[:, 1:, :] - zw[:, :-1, :]           # (i, k, j) 44 layers
U = v('U'); V = v('V')                        # (k, j, i_stag), (k, j_stag, i)
u_phy = tr3(0.5 * (U[:, :, :-1] + U[:, :, 1:]))
v_phy = tr3(0.5 * (V[:, :-1, :] + V[:, 1:, :]))
fgh, fgq = tr2(v('FGRNHFX')), tr2(v('FGRNQFX'))
if cfbm:
    import json
    j = json.load(open(cfbm))['vars']
    def arr(name):
        a = np.asarray(j[name]['data'], dtype=np.float64).reshape(j[name]['shape'], order='F')
        return a[5:85, 5:85]   # fire_flux dumps carry the memory halo ifms=-4..85; keep 1..80
    fgh, fgq = arr('grnhft'), arr('grnqft')
recs = {
 'fgrnhfx': fgh, 'fgrnqfx': fgq,
 'rho': tr3(v('RHO')), 'z_at_w': zw, 'dz8w': dz8w, 'mu': tr2(v('MUT')),
 'c1h': v('C1H'), 'c2h': v('C2H'),
 'u_phy': u_phy, 'v_phy': v_phy, 'znt': tr2(v('ZNT')), 'fz0': tr2(v('FZ0')),
 'uf': tr2(v('UF')), 'vf': tr2(v('VF')), 'fxlat': tr2(v('FXLAT')), 'fxlong': tr2(v('FXLONG')),
 'xlat': tr2(v('XLAT')), 'xlong': tr2(v('XLONG')),
 'grnhfx_wrf': tr2(v('GRNHFX')), 'grnqfx_wrf': tr2(v('GRNQFX')),
 'rthfrten_wrf': tr3(v('RTHFRTEN')), 'rqvfrten_wrf': tr3(v('RQVFRTEN')),
 'itimestep': np.asarray(d.variables['ITIMESTEP'][t]), 'xtime': np.asarray(d.variables['XTIME'][t]),
}
for k in ['DX', 'DY', 'CEN_LAT', 'CEN_LON', 'TRUELAT1', 'TRUELAT2', 'STAND_LON', 'DT']:
    recs[k.lower()] = np.asarray(getattr(d, k), dtype=np.float64)
with open(out, 'w') as f:
    for name, a in recs.items():
        a = np.asarray(a)
        dims = ' '.join(str(n) for n in a.shape)
        flat = a.reshape(-1, order='F')
        if np.issubdtype(a.dtype, np.integer):
            f.write(f"{name} int {a.ndim} {dims}\n" + ' '.join(str(int(x)) for x in flat) + '\n')
        else:
            f.write(f"{name} real {a.ndim} {dims}\n" + ' '.join(f"{float(x):.9g}" for x in flat) + '\n')
print('wrote', out, {k: np.asarray(a).shape for k, a in recs.items() if np.asarray(a).ndim})
