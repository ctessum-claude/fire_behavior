! EarthSciML kernel driver for the WETTING (rain) branch of fire_behavior's fuel-moisture model.
! The reference runs this repo has -- test7 and the coupled Last Chance Gulch case -- are both
! rain-free (fire_rain is identically 0 for the whole 12 h fire), so `r = rain_int - 0.05` never
! goes positive there and physics/fmc_wrffire_mod.F90:338-345 is never executed.  This driver calls
! fmc_wrffire_t%Advance_moisture_classes directly on an 8-cell strip whose accumulated rain is set
! to span the branch: below, at and above the 0.05 mm/h threshold, then 0.5, 2, 8 (= saturation_rain),
! 40 and 200 mm/h.  Two calls per run: the first at the fire step (dt_moisture = 0.5 s, so
! change = dt rlag << TOL and the update takes the second-order Taylor branch, as in every existing
! dump) and the second at dt_moisture = 600 s, which pushes the 1-h and 10-h classes past TOL = 1e-2
! and exercises the EXPONENTIAL branch -- also untested until now.
!
! Usage: ESM_DUMP=fmc ESM_DUMP_CALLS=1,2 ESM_DUMP_DIR=<dir> ./fmc_rain_driver
program fmc_rain_driver
  use fuel_anderson_mod, only : fuel_anderson_t
  use fmc_wrffire_mod, only : fmc_wrffire_t
  implicit none
  integer, parameter :: NI = 8, NJ = 1, NCLS = 5
  type (fuel_anderson_t) :: fuels
  type (fmc_wrffire_t) :: fmc
  real, dimension (NI, NJ) :: rain, rain_old, t2, t2_old, q2, q2_old, psfc, psfc_old, rh_fire, nfuel_cat, fmc_g
  real, dimension (NI) :: rain_rate                 ! mm/h
  real :: dt_fire, dt_long, fuelmc_g, fuelmc_g_live
  integer :: i

  call fuels%Initialization (1.0)                   ! fuelmc_c: canopy moisture, unused here
  fuelmc_g = 0.08; fuelmc_g_live = 0.3
  call fmc%Init (fuels, fuelmc_g, fuelmc_g_live, 1, NI, 1, NJ, 0, 0.5)

  ! forcing identical to cell (18,64) of test7 step 20, so the drying branch of cell 1 can be
  ! compared straight against the existing moisture tests
  t2 = 318.25555419921875; t2_old = t2
  q2 = 5.505851469933987e-3; q2_old = q2
  psfc = 84605.171875; psfc_old = psfc
  nfuel_cat = 2.0                                   ! timber with grass understory
  rh_fire = 0.0; fmc_g = 0.0

  rain_rate = [0.0, 0.04, 0.05, 0.5, 2.0, 8.0, 40.0, 200.0]
  dt_fire = 0.5; dt_long = 600.0

  ! --- call 1: the fire step.  Taylor branch, wetting active in cells 4..8.
  rain_old = 0.0
  do i = 1, NI
    rain(i, 1) = rain_rate(i) * dt_fire / 3600.0    ! accumulated rain [mm] over this step
  end do
  call set_state (0.05, 0.05, 0.05, 0.05, 0.30)
  fmc%dt_moisture = dt_fire
  call fmc%Advance_moisture_classes (.false., 1, NI, 1, NJ, 1, NI, 1, NJ, &
      rain, t2, q2, psfc, rain_old, t2_old, q2_old, psfc_old, rh_fire, fuelmc_g)
  call fmc%Average_moisture_classes (1, NI, 1, NJ, 1, NI, 1, NJ, nfuel_cat, fmc_g)
  print '(a)', 'call 1 (dt_moisture = 0.5 s):'
  do i = 1, NI
    print '(a, i2, a, f9.4, a, 5es14.6, a, es14.6)', '  cell ', i, ' rain_int =', rain_rate(i), &
        '  fmc_gc =', fmc%fmc_gc(i, :, 1), '  fmc_g =', fmc_g(i, 1)
  end do

  ! --- call 2: a 600 s moisture step, so change = dt rlag > TOL for the fast classes.
  rain_old = 0.0
  do i = 1, NI
    rain(i, 1) = rain_rate(i) * dt_long / 3600.0
  end do
  call set_state (0.05, 0.05, 0.05, 0.05, 0.30)
  fmc%dt_moisture = dt_long
  call fmc%Advance_moisture_classes (.false., 1, NI, 1, NJ, 1, NI, 1, NJ, &
      rain, t2, q2, psfc, rain_old, t2_old, q2_old, psfc_old, rh_fire, fuelmc_g)
  call fmc%Average_moisture_classes (1, NI, 1, NJ, 1, NI, 1, NJ, nfuel_cat, fmc_g)
  print '(a)', 'call 2 (dt_moisture = 600 s):'
  do i = 1, NI
    print '(a, i2, a, f9.4, a, 5es14.6, a, es14.6)', '  cell ', i, ' rain_int =', rain_rate(i), &
        '  fmc_gc =', fmc%fmc_gc(i, :, 1), '  fmc_g =', fmc_g(i, 1)
  end do

contains

  subroutine set_state (c1, c2, c3, c4, c5)
    real, intent (in) :: c1, c2, c3, c4, c5
    fmc%fmc_gc(:, 1, :) = c1; fmc%fmc_gc(:, 2, :) = c2; fmc%fmc_gc(:, 3, :) = c3
    fmc%fmc_gc(:, 4, :) = c4; fmc%fmc_gc(:, 5, :) = c5
  end subroutine set_state

end program fmc_rain_driver
