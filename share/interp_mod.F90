  module interp_mod

    use proj_lc_mod, only : proj_lc_t

#ifdef ESM_DUMP
    use module_esm_dump
#endif

    implicit none

    private

    integer, parameter :: HINTERP_NEAREST = 1, HINTERP_BILINEAR = 2, VINTERP_WINDS_FROM_3D_WINDS = 0, VINTERP_WINDS_FROM_10M_WINDS = 1

    public :: Interp_profile, Interp_horizontal_nearest, Interp_horizontal_bilinear, HINTERP_NEAREST, HINTERP_BILINEAR, &
        VINTERP_WINDS_FROM_3D_WINDS, VINTERP_WINDS_FROM_10M_WINDS

  contains

    subroutine Interp_horizontal_nearest (data_in, proj_data_in, ims, ime, jms, jme, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        lats_out, lons_out, data_out)

    ! Purpose: Nearest neighbor interpolation/extrapolation

      implicit none

      type (proj_lc_t), intent (in) :: proj_data_in
      integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, ims, ime, jms, jme
      real, dimension(ims:ime, jms:jme), intent (in) :: data_in
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lats_out, lons_out
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: data_out

      integer :: i, j, i_in, j_in
      real :: i_real, j_real


      do j = jfts, jfte
        do i = ifts, ifte
          call proj_data_in%Calc_ij (lats_out(i, j), lons_out(i, j), i_real, j_real)
          i_in = min (max (ims, nint (i_real)), ime)
          j_in = min (max (jms, nint (j_real)), ime)
          data_out(i, j) = data_in(i_in, j_in)
        end do
      end do

    end subroutine Interp_horizontal_nearest

    subroutine Interp_horizontal_bilinear (data_in, proj_data_in, ims, ime, jms, jme, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        lats_out, lons_out, data_out)

    ! Purpose: bi-linear interpolation + nearest neighbor extrapolation

      implicit none

      type (proj_lc_t), intent (in) :: proj_data_in
      integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, ims, ime, jms, jme
      real, dimension(ims:ime, jms:jme), intent (in) :: data_in
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lats_out, lons_out
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: data_out

      integer :: i, j, i0, j0, i1, j1
      real :: i_real, j_real, di, dj
#ifdef ESM_DUMP
      real, dimension (ifms:ifme, jfms:jfme) :: e_ireal, e_jreal, e_di, e_dj
      real, dimension (ifms:ifme, jfms:jfme) :: e_i0, e_j0
      logical :: e_dumping
      e_dumping = esm_dump_want ('hinterp')
      e_ireal = 0.0; e_jreal = 0.0; e_di = 0.0; e_dj = 0.0; e_i0 = 0.0; e_j0 = 0.0
#endif


      do j = jfts, jfte
        do i = ifts, ifte
          call proj_data_in%Calc_ij (lats_out(i, j), lons_out(i, j), i_real, j_real)

          i0 = max (ims, min (ime - 1, int (floor (i_real))))
          j0 = max (jms, min (jme - 1, int (floor (j_real))))
          i1 = i0 + 1
          j1 = j0 + 1

          di = max (0.0, min (1.0, i_real - real (i0)))
          dj = max (0.0, min (1.0, j_real - real (j0)))

          data_out(i, j) = (1.0 - di) * (1.0 - dj) * data_in(i0, j0) + &
              di * (1.0 - dj) * data_in(i1, j0) + &
              (1.0 - di) * dj * data_in(i0, j1) + &
              di * dj * data_in(i1, j1)
#ifdef ESM_DUMP
          e_ireal(i, j) = i_real; e_jreal(i, j) = j_real; e_di(i, j) = di; e_dj(i, j) = dj
          e_i0(i, j) = real (i0); e_j0(i, j) = real (j0)
#endif
        end do
      end do

#ifdef ESM_DUMP
      if (e_dumping) then
        call esm_dump_open ('hinterp')
        call esm_dump_var ('ims', ims); call esm_dump_var ('ime', ime); call esm_dump_var ('jms', jms); call esm_dump_var ('jme', jme)
        call esm_dump_var ('ifms', ifms); call esm_dump_var ('ifme', ifme); call esm_dump_var ('jfms', jfms); call esm_dump_var ('jfme', jfme)
        call esm_dump_var ('ifts', ifts); call esm_dump_var ('ifte', ifte); call esm_dump_var ('jfts', jfts); call esm_dump_var ('jfte', jfte)
        call esm_dump_var ('data_in', data_in); call esm_dump_var ('data_out', data_out)
        call esm_dump_var ('i_real', e_ireal); call esm_dump_var ('j_real', e_jreal)
        call esm_dump_var ('di', e_di); call esm_dump_var ('dj', e_dj)
        call esm_dump_var ('i0', e_i0); call esm_dump_var ('j0', e_j0)
        call esm_dump_var ('lats_out', lats_out); call esm_dump_var ('lons_out', lons_out)
        call esm_dump_close ()
      end if
#endif

    end subroutine Interp_horizontal_bilinear

    subroutine Interp_profile (fire_lsm_zcoupling, fire_lsm_zcoupling_ref, fire_wind_height, kfds, kfde, &
        uin, vin, z_at_w, z0f, uout, vout)

      implicit none

      real, intent (in) :: fire_wind_height, fire_lsm_zcoupling_ref
      integer, intent (in) :: kfds, kfde
      real, dimension(:), intent (in) :: uin, vin, z_at_w
      real, intent (in) :: z0f
      logical, intent (in) :: fire_lsm_zcoupling
      real, intent (out) :: uout, vout


      real, parameter :: VK_KAPPA = 0.4
      real, dimension (kfds:kfde - 1) :: altw, hgt
      integer :: k, kdmax
      real :: loght, loglast, logz0, logfwh, ht, r_nan, fire_wind_height_local, z0fc, &
          ust_d, wsf, wsf1, uf_temp, vf_temp


        ! max layer to interpolate from, can be less
      kdmax = kfde - 2
      do k = kfds, kdmax + 1
          ! altitude of the bottom w-point
!        altw(k) = phl(k) / G
        altw(k) = z_at_w(k)
      end do

      do k = kfds, kdmax
          ! height of the mass point above the ground
        hgt(k) = 0.5 * (altw(k) + altw(k + 1)) - altw(kfds)
      end do

        ! extrapolate mid-flame height from fire_lsm_zcoupling_ref?
      if (fire_lsm_zcoupling) then
        logfwh = log (fire_lsm_zcoupling_ref)
        fire_wind_height_local = fire_lsm_zcoupling_ref
      else
        logfwh = log (fire_wind_height)
        fire_wind_height_local = fire_wind_height
      end if

        ! interpolate u
      if (fire_wind_height_local > z0f)then
        do k = kfds, kdmax
          ht = hgt(k)
          if (ht >= fire_wind_height_local) then
              ! found layer k this point is in
            loght = log(ht)
            if (k == kfds) then
                ! first layer, log linear interpolation from 0 at zr
              logz0 = log(z0f)
              uout = uin(k) * (logfwh - logz0) / (loght - logz0)
              vout = vin(k) * (logfwh - logz0) / (loght - logz0)
            else
                ! log linear interpolation
              loglast = log (hgt(k - 1))
              uout = uin(k - 1) + (uin(k) - uin(k - 1)) * (logfwh - loglast) / (loght - loglast)
              vout = vin(k - 1) + (vin(k) - vin(k - 1)) * (logfwh - loglast) / (loght - loglast)
            end if
            exit
          end if
          if (k == kdmax) then
              ! last layer, still not high enough
            uout = uin(k)
            vout = vin(k)
          end if
        end do
      else
          ! roughness higher than the fire wind height
        uout = 0.0
        vout = 0.0
      end if

        ! Extrapol wind to target height
      if (fire_lsm_zcoupling) then
        uf_temp = uout
        vf_temp = vout
        wsf = max (sqrt (uf_temp ** 2.0 + vf_temp ** 2.0), 0.1)
        z0fc = z0f
        ust_d = wsf * VK_KAPPA / log(fire_lsm_zcoupling_ref / z0fc)
        wsf1 = (ust_d / VK_KAPPA) * log((fire_wind_height + z0fc) / z0fc)
        uout = wsf1 * uf_temp / wsf
        vout = wsf1 * vf_temp / wsf
      end if

    end subroutine Interp_profile

  end module interp_mod
