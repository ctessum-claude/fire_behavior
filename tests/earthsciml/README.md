# EarthSciML reference drivers

Standalone drivers that call this library's routines directly so that their inputs and
outputs can be dumped (`-DESM_DUMP=ON`, `ESM_DUMP=<scheme>`) and used as reference values
for the `.esm` transcriptions in the EqWeFiC project. They are *harnesses*, not part of the
model: nothing in `driver/`, `physics/`, `io/`, `share/` or `state/` calls them.

| file | what it replays | schemes dumped |
|---|---|---|
| `feedback_driver.F90` | `wrf_mod::Provide_atm_feedback` -> `Fire_tendency` and `interp_mod::Interp_profile` on one frame of the coupled test7 `wrf.nc` (fixed 19 x 19 x 44) | `wrfnc_wind` |
| `lastchance_feedback_driver.F90` | the same feedback chain on an arbitrary atmosphere TILE, geometry read from the flat file, so one binary serves any tile of any frame of a real Lambert run | `lastchance_feedback` |
| `fmc_rain_driver.F90` | `fmc_wrffire_t%Advance_moisture_classes` on an 8-cell strip whose accumulated rain spans the 0.05 mm/h threshold, at two step lengths | `fmc`, `fmc_avg` |
| `extract_wrfnc.py` | extracts one frame of test7's `wrf.nc` into the flat text format `feedback_driver` reads | - |
| `extract_tile.py` | extracts an atmosphere tile of a coupled WRF + CFBM run (`wrfout` + `fire_output`) for `lastchance_feedback_driver`; reconstructs `rho` the way `phy_prep` does and checks it against the mass-coordinate identity | - |
| `cell_facts.py` | prints the fire-cell facts (fuel category, burn rate, moisture, fluxes) of the 4 x 4 block under one WRF column | - |

Build: configure the library as usual with `-DESM_DUMP=ON`, then link a driver against the
static library and its module directory, e.g.

    f95 -O2 -ffree-line-length-none -I<build>/mod -o fmc_rain_driver \
        fmc_rain_driver.F90 <build>/libfirelib.a

`feedback_driver` and `lastchance_feedback_driver` expect a CFBM `namelist.fire` in the
working directory (only `fire_atm_feedback` is read).  The reference data these produced,
and the flat files they consume, are large and live outside git.
