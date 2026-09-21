#!/usr/bin/env python3
"""Extract one atmosphere TILE of the coupled Last Chance Gulch WRF+CFBM run (domain 5, ifire = 1,
sr = 4) into the flat text format read by lastchance_feedback_driver.F90.

The atmospheric density is not in wrfout, so it is reconstructed the way phy_prep does,
rho = (1 + qv) / alt with alt the WRF dry EOS inverse density
alt = r_d theta (1 + rv/rd qv) / p0 * (p0/p)^(cv/cp); the reconstruction is checked against the
mass-coordinate identity -(c1h mut + c2h) dnw = rho_d g dz8w and the check is printed.

Usage: extract_tile.py <time-string> <i0> <j0> <half> <out.flat>
       i0, j0 are 0-based WRF (west_east, south_north) indices of the target column;
       the tile is (2*half+1)^2 atmosphere cells centred on it.
"""
import sys, numpy as np
from netCDF4 import Dataset

RUN = '/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/lastchance/runs/full/'
G, RD, P0, RV = 9.81, 287.0, 1.0e5, 461.6
CP = 7.0 * RD / 2.0
CV = CP - RD
SR = 4

t, i0, j0, half, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
w = Dataset(RUN + 'wrfout_d05_' + t)
V = lambda n: np.asarray(w.variables[n][0], dtype=np.float64)

isl = slice(i0 - half, i0 + half + 1)
jsl = slice(j0 - half, j0 + half + 1)
n = 2 * half + 1
nz = len(w.dimensions['bottom_top'])

zw = (V('PH') + V('PHB'))[:, jsl, isl] / G                   # (nz+1, n, n)
dz8w = zw[1:] - zw[:-1]
p = (V('P') + V('PB'))[:, jsl, isl]
th = V('T')[:, jsl, isl] + 300.0
qv = V('QVAPOR')[:, jsl, isl]
alt = RD * th * (1.0 + (RV / RD) * qv) / P0 * (P0 / p) ** (CV / CP)
rho = (1.0 + qv) / alt
mut = (V('MU') + V('MUB'))[jsl, isl]
c1h, c2h, znw = V('C1H'), V('C2H'), V('ZNW')
dnw = znw[1:] - znw[:-1]

# consistency check of the reconstructed density against the mass coordinate
rho_d_id = -(c1h[:, None, None] * mut[None] + c2h[:, None, None]) * dnw[:, None, None] / (G * dz8w)
dev = np.abs((1.0 / alt) / rho_d_id - 1.0)[:nz - 1]
print('rho reconstruction vs mass-coordinate identity: max rel dev %.3e (centre column %.3e)'
      % (dev.max(), dev[:, half, half].max()))

f = Dataset(RUN + 'fire_output_' + t + '.nc')
fisl = slice(SR * (i0 - half), SR * (i0 + half + 1))
fjsl = slice(SR * (j0 - half), SR * (j0 + half + 1))
fgrnhfx = np.asarray(f.variables['fgrnhfx'][fjsl, fisl], dtype=np.float64)
fgrnqfx = np.asarray(f.variables['fgrnqfx'][fjsl, fisl], dtype=np.float64)

tr3 = lambda a: np.transpose(a, (2, 0, 1))   # (k, j, i) -> (i, k, j)
tr2 = lambda a: np.transpose(a, (1, 0))      # (j, i)    -> (i, j)
recs = [
    ('fgrnhfx', tr2(fgrnhfx)), ('fgrnqfx', tr2(fgrnqfx)),
    ('rho', tr3(rho)), ('dz8w', tr3(dz8w)), ('z_at_w', tr3(zw)),
    ('mu', tr2(mut)), ('c1h', c1h), ('c2h', c2h),
    ('grnhfx_wrf', tr2(V('GRNHFX')[jsl, isl])), ('grnqfx_wrf', tr2(V('GRNQFX')[jsl, isl])),
    ('rthfrten_wrf', tr3(V('RTHFRTEN')[:nz, jsl, isl])), ('rqvfrten_wrf', tr3(V('RQVFRTEN')[:nz, jsl, isl])),
]
with open(out, 'w') as fh:
    fh.write('dims int 1 4\n%d %d %d %d\n' % (n, n, nz, SR))
    fh.write('frame int 1 3\n%d %d %d\n' % (int(t.replace('-', '').replace(':', '').replace('_', '')[-6:]), i0, j0))
    for name, a in recs:
        a = np.asarray(a)
        fh.write('%s real %d %s\n' % (name, a.ndim, ' '.join(str(x) for x in a.shape)))
        fh.write(' '.join('%.9g' % x for x in a.reshape(-1, order='F')) + '\n')
print('wrote', out, 'tile %dx%dx%d centred on WRF column (i=%d, j=%d) of %s' % (n, n, nz, i0, j0, t))
