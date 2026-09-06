  module coupling_mod

    use proj_lc_mod, only : proj_lc_t
    use interp_mod, only : Interp_horizontal_nearest, Interp_horizontal_bilinear, HINTERP_NEAREST, HINTERP_BILINEAR, Interp_profile
#ifdef ESM_DUMP
    use module_esm_dump
#endif
    use stderrout_mod, only : Stop_simulation

    implicit none

    private

    public :: Interp_horizontal, Calc_fire_wind

  contains

    subroutine Interp_horizontal (data_in, proj_data_in, ims, ime, jms, jme, ifms, ifme, jfms, jfme, &
        num_tiles, i_start, i_end, j_start, j_end, hinterp_opt, lats_out, lons_out, data_out)

      implicit none

      type (proj_lc_t), intent (in) :: proj_data_in
      integer, intent (in) :: ifms, ifme, jfms, jfme, ims, ime, jms, jme, num_tiles, hinterp_opt
      integer, dimension (num_tiles), intent (in) :: i_start, i_end, j_start, j_end
      real, dimension(ims:ime, jms:jme), intent (in) :: data_in
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lats_out, lons_out
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: data_out

      integer :: ij, ifts, ifte, jfts, jfte


      Hinterp: select case (hinterp_opt)
        case (HINTERP_NEAREST)
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
          do ij = 1, num_tiles
            ifts = i_start(ij)
            ifte = i_end(ij)
            jfts = j_start(ij)
            jfte = j_end(ij)
            call Interp_horizontal_nearest (data_in, proj_data_in, ims, ime, jms, jme, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
                lats_out, lons_out, data_out)
          end do
          !$OMP END PARALLEL DO

        case (HINTERP_BILINEAR)
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
          do ij = 1, num_tiles
            ifts = i_start(ij)
            ifte = i_end(ij)
            jfts = j_start(ij)
            jfte = j_end(ij)
            call Interp_horizontal_bilinear (data_in, proj_data_in, ims, ime, jms, jme, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
                lats_out, lons_out, data_out)
          end do
          !$OMP END PARALLEL DO

        case default
          call Stop_simulation ('The horizontal interpolation option selected does not exist.')

      end select Hinterp

    end subroutine Interp_horizontal

    subroutine Calc_fire_wind (u3d, v3d, z_at_w, z0, iims, iime, jims, jime, kims, kime, fire_lsm_zcoupling, &
        fire_lsm_zcoupling_ref, fire_wind_height, ioms, iome, joms, jome, iops, iope, jops, jope, &
        u_out, v_out, cap_winds)

      implicit none

      integer, intent (in) :: iims, iime, jims, jime, kims, kime, ioms, iome, joms, jome, iops, iope, jops, jope
      real, intent (in) :: fire_wind_height, fire_lsm_zcoupling_ref
      logical, intent (in) :: fire_lsm_zcoupling
      real, dimension(iims:iime, jims:jime, kims:kime), intent (in) :: u3d, v3d, z_at_w
      real, dimension(iims:iime, jims:jime), intent (in) :: z0
      real, dimension(ioms:iome, joms:joms), intent (out) :: u_out, v_out
      logical, intent (in), optional :: cap_winds

      real :: wspd
      integer :: i, j
      logical :: cap_winds_flag


!      print *, 'shape u3d = ', shape (u3d)
!      print *, 'shape v3d = ', shape (v3d)
!      print *, 'shape z_at_w = ', shape (z_at_w)
!      print *, 'shape u_out = ', shape (u_out)
!      print *, 'shape v_out = ', shape (v_out)
!      print *, 'shape z0 = ', shape (z0)

      if (present (cap_winds)) then
        cap_winds_flag = cap_winds
      else
        cap_winds_flag = .false.
      end if

      do j = jops, jope
        do i = iops, iope
          call Interp_profile (fire_lsm_zcoupling, fire_lsm_zcoupling_ref, fire_wind_height, kims, kime, &
              u3d(i, j, :), v3d(i, j, :), z_at_w(i, j, :), z0(i, j), u_out(i, j), v_out(i, j))
        end do
      end do

#ifdef ESM_DUMP
        ! EarthSciML instrumentation: inputs and outputs of the log-profile wind interpolation
        ! (Interp_profile) for a 10x10 block of fire cells at the tile origin; the 3-D fields are the
        ! atmosphere winds/heights already mapped onto the fire mesh (standalone path). k = 1 is the
        ! lowest layer; z_at_w has kime - kims + 1 levels (WRF's z_at_w with the top level dropped).
      if (esm_dump_want ('fire_wind')) then
        call esm_dump_open ('fire_wind')
        call esm_dump_var ('iops', iops); call esm_dump_var ('jops', jops); call esm_dump_var ('kims', kims); call esm_dump_var ('kime', kime)
        call esm_dump_var ('fire_wind_height', fire_wind_height); call esm_dump_var ('fire_lsm_zcoupling_ref', fire_lsm_zcoupling_ref)
        call esm_dump_var ('fire_lsm_zcoupling', fire_lsm_zcoupling); call esm_dump_var ('cap_winds', cap_winds_flag)
        call esm_dump_var ('u3d', u3d(iops:min (iope, iops + 9), jops:min (jope, jops + 9), :))
        call esm_dump_var ('v3d', v3d(iops:min (iope, iops + 9), jops:min (jope, jops + 9), :))
        call esm_dump_var ('z_at_w', z_at_w(iops:min (iope, iops + 9), jops:min (jope, jops + 9), :))
        call esm_dump_var ('z0', z0(iops:min (iope, iops + 9), jops:min (jope, jops + 9)))
        call esm_dump_var ('uf', u_out(iops:min (iope, iops + 9), jops:min (jope, jops + 9)))
        call esm_dump_var ('vf', v_out(iops:min (iope, iops + 9), jops:min (jope, jops + 9)))
        call esm_dump_close ()
      end if
#endif

        ! To avoid arithmatic error
      if (cap_winds_flag) then
        do j = jops, jope
          do i = iops, iope
            wspd = sqrt (u_out(i, j) ** 2 + v_out(i, j) ** 2)
            if (wspd < 0.001) then
              u_out(i, j) = sign (0.001, u_out(i, j))
              v_out(i, j) = sign (0.001, v_out(i, j))
            end if
          end do
        end do
      end if

    end subroutine Calc_fire_wind

  end module coupling_mod
