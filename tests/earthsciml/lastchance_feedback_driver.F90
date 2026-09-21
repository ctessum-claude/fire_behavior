! EarthSciML reference driver: replays fire_behavior's WRF-side fire -> atmosphere coupling
! (wrf_mod::Provide_atm_feedback -> Fire_tendency) on a small atmosphere TILE extracted from a
! frame of the coupled Last Chance Gulch WRF+CFBM run (domain 5, ifire = 1, sr = 4), and dumps
! inputs/outputs as JSON through the fork's ESM_DUMP hooks.
!
! Unlike data/eqwefic/fire_build/feedback_driver (fixed 19x19x44 test7 dims) the tile geometry is
! read from the flat file, so the same binary serves any tile of any frame.  Usage:
!   lastchance_feedback_driver <tile.flat>     (run in a directory holding the CFBM namelist.fire)
program lastchance_feedback_driver
  use namelist_mod, only : namelist_t
  use wrf_mod, only : Provide_atm_feedback
  use module_esm_dump
  implicit none
  type (namelist_t) :: nml
  character (len = 512) :: flat
  integer :: nx, ny, nz, sr, i, j, k, itile, jtile, iframe
  real, allocatable, dimension (:, :) :: fgrnhfx, fgrnqfx, emis_smoke
  real, allocatable, dimension (:, :) :: mu, grnhfx_wrf, grnqfx_wrf
  real, allocatable, dimension (:, :) :: grnhfx, grnqfx, canhfx, canqfx, grnsmk, aod
  real, allocatable, dimension (:, :, :) :: rho, dz8w, z_at_w, zeros3, rthfrten, rqvfrten, rthfrten_wrf, rqvfrten_wrf
  real, allocatable, dimension (:) :: c1h, c2h
  real :: alfg, alfc, z1can, dmax, dref

  call get_command_argument (1, flat)
  call nml%Initialization ('namelist.fire')
  call read_flat (trim (flat))
  alfg = 50.0; alfc = 50.0; z1can = 15.0      ! WRF Registry defaults fire_ext_grnd, fire_ext_crwn, fire_crwn_hgt
  emis_smoke = 0.0; zeros3 = 0.0

  ! Tile treated as a standalone domain: ids = 1, ide = nx + 1 (WRF e_we convention), likewise jds/jde.
  ! Provide_atm_feedback aggregates i = 2 .. nx - 1, j = 2 .. ny - 1, so the tile centre is interior.
  call Provide_atm_feedback (nml, 1, sr * nx, 1, sr * ny, 1, sr * nx, 1, sr * ny, 1, sr * nx, 1, sr * ny, &
      1, nx + 1, 1, nz + 1, 1, ny + 1, 1, nx, 1, nz + 1, 1, ny, 1, nx, 1, nz + 1, 1, ny, sr, sr, &
      emis_smoke, tracer_opt = 0, p_phy = zeros3, t_phy = zeros3, qv = zeros3, aod5502d_smoke = aod, &
      fgrnhfx = fgrnhfx, fgrnqfx = fgrnqfx, grnhfx = grnhfx, grnqfx = grnqfx, canhfx = canhfx, canqfx = canqfx, &
      grnsmk = grnsmk, alfg = alfg, alfc = alfc, z1can = z1can, rho = rho, dz8w = dz8w, z_at_w = z_at_w, &
      mu = mu, c1h = c1h, c2h = c2h, rthfrten = rthfrten, rqvfrten = rqvfrten)

  i = nx / 2 + 1; j = ny / 2 + 1                 ! tile centre = the requested WRF column
  print '(a, 2i4)', 'tile centre (i, j) =', i, j
  print '(a, 2es14.6)', 'grnhfx  driver / wrfout  = ', grnhfx(i, j), grnhfx_wrf(i, j)
  print '(a, 2es14.6)', 'grnqfx  driver / wrfout  = ', grnqfx(i, j), grnqfx_wrf(i, j)
  dmax = maxval (abs (rthfrten(i, 1:nz, j) - rthfrten_wrf(i, 1:nz, j)))
  dref = maxval (abs (rthfrten_wrf(i, 1:nz, j)))
  print '(a, es12.4, a, es12.4)', 'rthfrten: max |driver - wrfout| = ', dmax, '  max |wrfout| = ', dref
  dmax = maxval (abs (rqvfrten(i, 1:nz, j) - rqvfrten_wrf(i, 1:nz, j)))
  dref = maxval (abs (rqvfrten_wrf(i, 1:nz, j)))
  print '(a, es12.4, a, es12.4)', 'rqvfrten: max |driver - wrfout| = ', dmax, '  max |wrfout| = ', dref

  if (esm_dump_want ('lastchance_feedback')) then
    call esm_dump_open ('lastchance_feedback')
    call esm_dump_var ('nx', nx); call esm_dump_var ('ny', ny); call esm_dump_var ('nz', nz); call esm_dump_var ('sr', sr)
    call esm_dump_var ('iframe', iframe); call esm_dump_var ('itile', itile); call esm_dump_var ('jtile', jtile)
    call esm_dump_var ('ic', i); call esm_dump_var ('jc', j)
    call esm_dump_var ('alfg', alfg); call esm_dump_var ('alfc', alfc); call esm_dump_var ('z1can', z1can)
    call esm_dump_var ('fire_atm_feedback', nml%fire_atm_feedback)
    call esm_dump_var ('fgrnhfx_blk', fgrnhfx(sr * (i - 1) + 1 : sr * i, sr * (j - 1) + 1 : sr * j))
    call esm_dump_var ('fgrnqfx_blk', fgrnqfx(sr * (i - 1) + 1 : sr * i, sr * (j - 1) + 1 : sr * j))
    call esm_dump_var ('mu', mu(i, j)); call esm_dump_var ('c1h', c1h(1:nz)); call esm_dump_var ('c2h', c2h(1:nz))
    call esm_dump_var ('rho', rho(i, 1:nz, j)); call esm_dump_var ('dz8w', dz8w(i, 1:nz, j))
    call esm_dump_var ('z_at_w', z_at_w(i, 1:nz + 1, j))
    call esm_dump_var ('grnhfx', grnhfx(i, j)); call esm_dump_var ('grnqfx', grnqfx(i, j))
    call esm_dump_var ('rthfrten', rthfrten(i, 1:nz, j)); call esm_dump_var ('rqvfrten', rqvfrten(i, 1:nz, j))
    call esm_dump_var ('grnhfx_wrfout', grnhfx_wrf(i, j)); call esm_dump_var ('grnqfx_wrfout', grnqfx_wrf(i, j))
    call esm_dump_var ('rthfrten_wrfout', rthfrten_wrf(i, 1:nz, j)); call esm_dump_var ('rqvfrten_wrfout', rqvfrten_wrf(i, 1:nz, j))
    call esm_dump_close ()
  end if

contains

  subroutine read_flat (path)
    character (len = *), intent (in) :: path
    character (len = 64) :: name, kind
    integer :: u, rank, ios, n(3)
    real, allocatable :: buf(:)

    open (newunit = u, file = path, status = 'old', action = 'read')
    ! first record must be: dims int 1 4 / nx ny nz sr
    read (u, *) name, kind, rank, n(1)
    read (u, *) nx, ny, nz, sr
    read (u, *) name, kind, rank, n(1)
    read (u, *) iframe, itile, jtile
    allocate (fgrnhfx(sr * nx, sr * ny), fgrnqfx(sr * nx, sr * ny), emis_smoke(sr * nx, sr * ny))
    allocate (mu(nx, ny), grnhfx_wrf(nx, ny), grnqfx_wrf(nx, ny))
    allocate (grnhfx(nx, ny), grnqfx(nx, ny), canhfx(nx, ny), canqfx(nx, ny), grnsmk(nx, ny), aod(nx, ny))
    allocate (rho(nx, nz + 1, ny), dz8w(nx, nz + 1, ny), z_at_w(nx, nz + 1, ny), zeros3(nx, nz + 1, ny))
    allocate (rthfrten(nx, nz + 1, ny), rqvfrten(nx, nz + 1, ny), rthfrten_wrf(nx, nz + 1, ny), rqvfrten_wrf(nx, nz + 1, ny))
    allocate (c1h(nz + 1), c2h(nz + 1))
    rthfrten_wrf = 0.0; rqvfrten_wrf = 0.0
    do
      read (u, *, iostat = ios) name, kind, rank, n(1:rank)
      if (ios /= 0) exit
      allocate (buf(max (1, product (n(1:rank)))))
      read (u, *) buf(1:max (1, product (n(1:rank))))
      select case (trim (name))
      case ('fgrnhfx'); fgrnhfx = reshape (buf, [sr * nx, sr * ny])
      case ('fgrnqfx'); fgrnqfx = reshape (buf, [sr * nx, sr * ny])
      case ('mu'); mu = reshape (buf, [nx, ny])
      case ('grnhfx_wrf'); grnhfx_wrf = reshape (buf, [nx, ny])
      case ('grnqfx_wrf'); grnqfx_wrf = reshape (buf, [nx, ny])
      case ('rho'); rho(:, 1:nz, :) = reshape (buf, [nx, nz, ny]); rho(:, nz + 1, :) = rho(:, nz, :)
      case ('dz8w'); dz8w(:, 1:nz, :) = reshape (buf, [nx, nz, ny]); dz8w(:, nz + 1, :) = dz8w(:, nz, :)
      case ('z_at_w'); z_at_w = reshape (buf, [nx, nz + 1, ny])
      case ('rthfrten_wrf'); rthfrten_wrf(:, 1:nz, :) = reshape (buf, [nx, nz, ny])
      case ('rqvfrten_wrf'); rqvfrten_wrf(:, 1:nz, :) = reshape (buf, [nx, nz, ny])
      case ('c1h'); c1h(1:nz) = buf(1:nz); c1h(nz + 1) = c1h(nz)
      case ('c2h'); c2h(1:nz) = buf(1:nz); c2h(nz + 1) = c2h(nz)
      end select
      deallocate (buf)
    end do
    close (u)
  end subroutine read_flat

end program lastchance_feedback_driver
