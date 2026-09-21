#!/usr/bin/env python3
"""Print the fire-cell facts of the coupled cell under WRF column (i0, j0) at a Last Chance frame."""
import sys, numpy as np
from netCDF4 import Dataset
RUN='/projects/illinois/eng/cee/ctessum/ctessum/data/eqwefic/lastchance/runs/full/'
t,i0,j0 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
f=Dataset(RUN+'fire_output_'+t+'.nc')
FGI=[0.166,0.896,0.674,3.591,0.784,1.344,1.091,1.120,0.780,2.692,2.582,7.749,13.024,1.e-7]
WEIGHT=None
sub=lambda n: np.asarray(f.variables[n][4*j0:4*j0+4, 4*i0:4*i0+4])
gh, gq = sub('fgrnhfx'), sub('fgrnqfx')
cat, ffb, fmc = sub('nfuel_cat').astype(int), sub('fuel_frac_burnt_dt'), sub('fmc_g')
ff, fa, lfn = sub('fuel_frac'), sub('fire_area'), sub('lfn')
uf, vf = sub('uf'), sub('vf')
# block order: Fortran (i fastest) within the 4x4 block -> s = ci + 4*cj + 1 with ci the x offset
jc, ic = np.unravel_index(gh.argmax(), gh.shape)   # arrays are [j, i]
s = ic + 4*jc + 1
print('block fgrnhfx (Fortran order s=1..16):', [float(gh[cj,ci]) for cj in range(4) for ci in range(4)])
print('block fgrnqfx (Fortran order s=1..16):', [float(gq[cj,ci]) for cj in range(4) for ci in range(4)])
print('coupled cell: s=%d  ci=%d cj=%d  fire (i,j) 0-based = (%d,%d)'%(s, ic, jc, 4*i0+ic, 4*j0+jc))
print('  nfuel_cat=%d  fgi=%.6g  fuel_frac_burnt_dt=%.9g  fmc_g=%.9g'%(cat[jc,ic], FGI[cat[jc,ic]-1], ffb[jc,ic], fmc[jc,ic]))
print('  fgrnhfx=%.9g  fgrnqfx=%.9g  fuel_frac=%.9g  fire_area=%.9g  lfn=%.9g'%(gh[jc,ic],gq[jc,ic],ff[jc,ic],fa[jc,ic],lfn[jc,ic]))
print('  uf=%.9g vf=%.9g'%(uf[jc,ic],vf[jc,ic]))
print('  all cats in block:', sorted(set(cat.ravel().tolist())))
