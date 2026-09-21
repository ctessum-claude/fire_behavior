! EarthSciML reference driver: replays fire_behavior's WRF-side coupling routines
! (wrf_mod::Provide_atm_feedback -> Fire_tendency, and interp_mod::Interp_profile with the
! WRF nearest-column mapping) on one frame of the coupled WRF+CFBM output test7/wrf.nc, and
! dumps inputs/outputs as JSON through the fork's ESM_DUMP hooks.  Usage:
!   feedback_driver <frame.flat>      (run in a directory holding the CFBM namelist.fire)
program feedback_driver
  use namelist_mod, only : namelist_t
  use proj_lc_mod, only : proj_lc_t
  use interp_mod, only : Interp_profile
  use wrf_mod, only : Provide_atm_feedback
  use module_esm_dump
  implicit none
  integer, parameter :: NX = 19, NY = 19, NZ = 44, NFX = 80, NFY = 80, SR = 4
  type (namelist_t) :: nml
  type (proj_lc_t) :: proj
  character (len = 512) :: flat
  real, dimension (NFX, NFY) :: fgrnhfx, fgrnqfx, fz0, uf_wrf, vf_wrf, fxlat, fxlong, emis_smoke
  real, dimension (NX, NY) :: mu, znt, xlat, xlong, grnhfx_wrf, grnqfx_wrf
  real, dimension (NX, NZ + 1, NY) :: rho, dz8w, z_at_w, u_phy, v_phy, zeros3, rthfrten, rqvfrten, rthfrten_wrf, rqvfrten_wrf
  real, dimension (NZ + 1) :: c1h, c2h
  real, dimension (NX, NY) :: grnhfx, grnqfx, canhfx, canqfx, grnsmk, aod
  real :: dx, dy, cen_lat, cen_lon, truelat1, truelat2, stand_lon, alfg, alfc, z1can
  real, dimension (10, 10) :: uf_calc, vf_calc, z0f_blk
  integer, dimension (10, 10) :: iw_blk, jw_blk
  real :: i_real, j_real, uo, vo, dmax, dref
  integer :: i, j, k, i_wrf, j_wrf, itimestep

  call get_command_argument (1, flat)
  call nml%Initialization ('namelist.fire')
  block
    character (len = 64) :: env
    call get_environment_variable ('FIRE_WIND_HEIGHT', env)
    if (len_trim (env) > 0) read (env, *) nml%fire_wind_height
    call get_environment_variable ('FIRE_LSM_ZCOUPLING', env)
    if (len_trim (env) > 0) nml%fire_lsm_zcoupling = (trim (env) == '1')
  end block
  call read_flat (trim (flat))
  alfg = 50.0; alfc = 50.0; z1can = 15.0      ! WRF Registry defaults fire_ext_grnd, fire_ext_crwn, fire_crwn_hgt
  emis_smoke = 0.0; zeros3 = 0.0

  ! --- fire -> atmosphere: aggregation to the atmosphere grid and the theta/qv tendencies
  call Provide_atm_feedback (nml, 1, NFX, 1, NFY, 1, NFX, 1, NFY, 1, NFX, 1, NFY, &
      1, NX + 1, 1, NZ + 1, 1, NY + 1, 1, NX, 1, NZ + 1, 1, NY, 1, NX, 1, NZ, 1, NY, SR, SR, &
      emis_smoke, tracer_opt = 0, p_phy = zeros3, t_phy = zeros3, qv = zeros3, aod5502d_smoke = aod, &
      fgrnhfx = fgrnhfx, fgrnqfx = fgrnqfx, grnhfx = grnhfx, grnqfx = grnqfx, canhfx = canhfx, canqfx = canqfx, &
      grnsmk = grnsmk, alfg = alfg, alfc = alfc, z1can = z1can, rho = rho, dz8w = dz8w, z_at_w = z_at_w, &
      mu = mu, c1h = c1h, c2h = c2h, rthfrten = rthfrten, rqvfrten = rqvfrten)
  dmax = maxval (abs (rthfrten(2:NX-1, 1:NZ-1, 2:NY-1) - rthfrten_wrf(2:NX-1, 1:NZ-1, 2:NY-1)))
  dref = maxval (abs (rthfrten_wrf(2:NX-1, 1:NZ-1, 2:NY-1)))
  print '(a, es12.4, a, es12.4)', 'rthfrten: max |driver - wrf.nc| = ', dmax, '  max |wrf.nc| = ', dref
  dmax = maxval (abs (grnhfx(2:NX-1, 2:NY-1) - grnhfx_wrf(2:NX-1, 2:NY-1)))
  print '(a, es12.4, a, es12.4)', 'grnhfx: max |driver - wrf.nc| = ', dmax, '  max |wrf.nc| = ', maxval (grnhfx_wrf)

  ! --- atmosphere -> fire winds for the 10x10 fire-cell block at the domain origin (WRF path)
  proj = proj_lc_t (cen_lat = cen_lat, cen_lon = cen_lon, dx = dx, dy = dy, standard_lon = stand_lon, &
      true_lat_1 = truelat1, true_lat_2 = truelat2, nx = NX, ny = NY)
  do j = 1, 10
    do i = 1, 10
      call proj%Calc_ij (fxlat(i, j), fxlong(i, j), i_real, j_real)
      i_wrf = min (max (1, nint (i_real)), NX); j_wrf = min (max (1, nint (j_real)), NY)
      iw_blk(i, j) = i_wrf; jw_blk(i, j) = j_wrf; z0f_blk(i, j) = fz0(i, j)
      call Interp_profile (nml%fire_lsm_zcoupling, nml%fire_lsm_zcoupling_ref, nml%fire_wind_height, 1, NZ + 1, &
          u_phy(i_wrf, :, j_wrf), v_phy(i_wrf, :, j_wrf), z_at_w(i_wrf, :, j_wrf), fz0(i, j), uo, vo)
      uf_calc(i, j) = uo; vf_calc(i, j) = vo
    end do
  end do
  print '(a, es12.4, a, es12.4)', 'uf: max |Interp_profile(state at frame) - UF(frame)| = ', &
      maxval (abs (uf_calc - uf_wrf(1:10, 1:10))), '  max |UF| = ', maxval (abs (uf_wrf(1:10, 1:10)))
  if (esm_dump_want ('wrfnc_wind')) then
    call esm_dump_open ('wrfnc_wind')
    call esm_dump_var ('itimestep', itimestep)
    call esm_dump_var ('fire_wind_height', nml%fire_wind_height); call esm_dump_var ('fire_lsm_zcoupling', nml%fire_lsm_zcoupling)
    call esm_dump_var ('fire_lsm_zcoupling_ref', nml%fire_lsm_zcoupling_ref)
    call esm_dump_var ('i_wrf', iw_blk); call esm_dump_var ('j_wrf', jw_blk)
    call esm_dump_var ('u_phy', u_phy(:, 1:NZ, :)); call esm_dump_var ('v_phy', v_phy(:, 1:NZ, :)); call esm_dump_var ('z_at_w', z_at_w)
    call esm_dump_var ('z0f', z0f_blk); call esm_dump_var ('uf_calc', uf_calc); call esm_dump_var ('vf_calc', vf_calc)
    call esm_dump_var ('uf_wrfnc', uf_wrf(1:10, 1:10)); call esm_dump_var ('vf_wrfnc', vf_wrf(1:10, 1:10))
    call esm_dump_var ('rthfrten_wrfnc', rthfrten_wrf(:, 1:NZ, :)); call esm_dump_var ('rqvfrten_wrfnc', rqvfrten_wrf(:, 1:NZ, :))
    call esm_dump_var ('grnhfx_wrfnc', grnhfx_wrf); call esm_dump_var ('grnqfx_wrfnc', grnqfx_wrf)
    call esm_dump_close ()
  end if

contains

  subroutine read_flat (path)
    character (len = *), intent (in) :: path
    character (len = 64) :: name, kind
    integer :: u, rank, ios, n(3), i
    real, allocatable :: buf(:)
    open (newunit = u, file = path, status = 'old', action = 'read')
    do
      read (u, *, iostat = ios) name, kind, rank, n(1:rank)
      if (ios /= 0) exit
      allocate (buf(max (1, product (n(1:rank)))))
      if (rank == 0) then
        read (u, *) buf(1)
      else
        read (u, *) buf(1:product (n(1:rank)))
      end if
      select case (trim (name))
      case ('fgrnhfx'); fgrnhfx = reshape (buf, [NFX, NFY])
      case ('fgrnqfx'); fgrnqfx = reshape (buf, [NFX, NFY])
      case ('fz0'); fz0 = reshape (buf, [NFX, NFY])
      case ('uf'); uf_wrf = reshape (buf, [NFX, NFY])
      case ('vf'); vf_wrf = reshape (buf, [NFX, NFY])
      case ('fxlat'); fxlat = reshape (buf, [NFX, NFY])
      case ('fxlong'); fxlong = reshape (buf, [NFX, NFY])
      case ('mu'); mu = reshape (buf, [NX, NY])
      case ('znt'); znt = reshape (buf, [NX, NY])
      case ('xlat'); xlat = reshape (buf, [NX, NY])
      case ('xlong'); xlong = reshape (buf, [NX, NY])
      case ('grnhfx_wrf'); grnhfx_wrf = reshape (buf, [NX, NY])
      case ('grnqfx_wrf'); grnqfx_wrf = reshape (buf, [NX, NY])
      case ('rho'); rho(:, 1:NZ, :) = reshape (buf, [NX, NZ, NY]); rho(:, NZ + 1, :) = rho(:, NZ, :)
      case ('dz8w'); dz8w(:, 1:NZ, :) = reshape (buf, [NX, NZ, NY]); dz8w(:, NZ + 1, :) = dz8w(:, NZ, :)
      case ('u_phy'); u_phy(:, 1:NZ, :) = reshape (buf, [NX, NZ, NY]); u_phy(:, NZ + 1, :) = 0.0
      case ('v_phy'); v_phy(:, 1:NZ, :) = reshape (buf, [NX, NZ, NY]); v_phy(:, NZ + 1, :) = 0.0
      case ('z_at_w'); z_at_w = reshape (buf, [NX, NZ + 1, NY])
      case ('rthfrten_wrf'); rthfrten_wrf = reshape (buf, [NX, NZ + 1, NY])
      case ('rqvfrten_wrf'); rqvfrten_wrf = reshape (buf, [NX, NZ + 1, NY])
      case ('c1h'); c1h(1:NZ) = buf(1:NZ); c1h(NZ + 1) = c1h(NZ)
      case ('c2h'); c2h(1:NZ) = buf(1:NZ); c2h(NZ + 1) = c2h(NZ)
      case ('dx'); dx = buf(1)
      case ('dy'); dy = buf(1)
      case ('cen_lat'); cen_lat = buf(1)
      case ('cen_lon'); cen_lon = buf(1)
      case ('truelat1'); truelat1 = buf(1)
      case ('truelat2'); truelat2 = buf(1)
      case ('stand_lon'); stand_lon = buf(1)
      case ('itimestep'); itimestep = nint (buf(1))
      end select
      deallocate (buf)
    end do
    close (u)
  end subroutine read_flat

end program feedback_driver
