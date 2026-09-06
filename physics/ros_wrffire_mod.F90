#ifdef ESM_DUMP
#define ESM_PURE
#else
#define ESM_PURE pure
#endif
  module ros_wrffire_mod

    use constants_mod, only : CMBCNST, CONVERT_J_PER_KG_TO_BTU_PER_POUND
#ifdef ESM_DUMP
    use module_esm_dump
#endif
    use fuel_mod, only : fuel_t
    use namelist_mod, only : namelist_t
    use ros_mod, only : ros_t
    use state_mod, only : state_fire_t

    implicit none

    private

    public :: ros_wrffire_t
#ifdef ESM_DUMP
    public :: Esm_dump_ros_diag
      ! EarthSciML instrumentation: per-cell intermediates of Calc_ros_wrffire
      ! (module-level because Calc_ros has intent(in) this); filled on every
      ! call, dumped by Esm_dump_ros_diag from the level-set tendency dump.
    real, dimension(:, :), allocatable, save :: e_speed, e_tanphi, e_cor_wind, e_cor_slope, &
        e_umid, e_phiw, e_phis, e_ros_base, e_ros_wind, e_ros_slope
#endif

    logical, parameter :: FIRE_GROWS_ONLY = .true.
    integer, parameter :: SLOPE_FACTOR = 1.0

      ! fuelheat: fuel particle low heat content [btu/lb]
    real, parameter :: FUELHEAT = CMBCNST * CONVERT_J_PER_KG_TO_BTU_PER_POUND
    integer, parameter :: FIRE_ADVECTION = 1 ! "0 = fire spread computed from normal wind speed/slope, 1 = fireline particle speed projected on normal" "0"

    type, extends(ros_t) :: ros_wrffire_t
      real, dimension(:, :), allocatable :: bbb, ischap, betafl, phiwc, r_0
    contains
      procedure, public :: Calc_ros => Calc_ros_wrffire
      procedure, public :: Init => Init_ros_wrffire
      procedure, public :: Set_params => Set_ros_parameters_wrffire
    end type ros_wrffire_t

  contains

    ESM_PURE function Calc_ros_wrffire (this, ifms, ifme, jfms, jfme, i, j, nvx, nvy, uf, vf, dzdxf, dzdyf) result (return_value)

      implicit none

      ! m/s = (ft/min) * 0.3048 / 60.0 = (ft/min) * .00508
      ! ft/min = m/s * 2.2369 * 88.0 = m/s *  196.850

      class (ros_wrffire_t), intent (in) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme, i, j
      real, intent (in) :: nvx, nvy, uf, vf, dzdxf, dzdyf
      real :: return_value

      real :: speed, tanphi ! windspeed and slope in the direction normal to the fireline
      real :: umid, phis, phiw, spdms, umidm, excess, ros_back, cor_wind, cor_slope, ros_base, ros_wind, ros_slope
      real, parameter :: ROS_MAX = 6.0


#ifdef ESM_DUMP
      umid = 0.0; phiw = 0.0; phis = 0.0
#endif
      if (FIRE_ADVECTION /= 0) then
          ! wind speed is total speed 
        speed = sqrt (uf * uf + vf * vf) + tiny (speed)
          ! slope is total slope
        tanphi = sqrt (dzdxf * dzdxf + dzdyf * dzdyf) + tiny (tanphi)
          ! cos of wind and spread, if >0
        cor_wind =  max (0.0, (uf * nvx + vf * nvy) / speed)
          ! cos of slope and spread, if >0
        cor_slope = max (0.0, (dzdxf * nvx + dzdyf * nvy) / tanphi)
      else
          ! wind speed in spread direction
        speed = uf * nvx + vf * nvy
          ! slope in spread direction
        tanphi = dzdxf * nvx + dzdyf * nvy
        cor_wind = 1.0
        cor_slope = 1.0
      end if

      if (.not. this%ischap(i, j) > 0.0) then
          ! Rothermel
        spdms = max (speed, 0.0)
        umidm = min (spdms, 30.0)
        umid = umidm * 196.850 ! m/s to ft/min
        phiw = umid ** this%bbb(i, j) * this%phiwc(i, j)
        phis = 0.0
        if (tanphi > 0.0) phis = 5.275 * (this%betafl(i, j)) ** (-0.3) * tanphi ** 2
        ros_base = this%r_0(i, j) * 0.00508 ! ft/min to m/s
        ros_wind = ros_base * phiw
        ros_slope = ros_base * phis
      else
        spdms = max (speed, 0.0)
        ros_back = 0.03333    ! chaparral backing fire spread rate 0.033 m/s   ! param!
          ! spread rate, m/s
        ros_wind = 1.2974 * spdms ** 1.41
        ros_wind = max (ros_wind, ros_back)
        ros_slope = 0.0
        ros_base = 0.0
      end if

      ros_wind = ros_wind * cor_wind
      ros_slope = ros_slope * cor_slope

      return_value = min (ros_base + ros_wind + SLOPE_FACTOR * ros_slope, ROS_MAX)
      if (FIRE_GROWS_ONLY) return_value = max (return_value, 0.0)
#ifdef ESM_DUMP
      if (allocated (e_speed)) then
        e_speed(i, j) = speed; e_tanphi(i, j) = tanphi; e_cor_wind(i, j) = cor_wind; e_cor_slope(i, j) = cor_slope
        e_umid(i, j) = umid; e_phiw(i, j) = phiw; e_phis(i, j) = phis
        e_ros_base(i, j) = ros_base; e_ros_wind(i, j) = ros_wind; e_ros_slope(i, j) = ros_slope
      end if
#endif

    end function Calc_ros_wrffire

    subroutine Init_ros_wrffire (this, ifms, ifme, jfms, jfme)

      implicit none

      class (ros_wrffire_t), intent (in out) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme


      allocate (this%iboros(ifms:ifme, jfms:jfme))
 
      allocate (this%ischap(ifms:ifme, jfms:jfme))
      allocate (this%betafl(ifms:ifme, jfms:jfme))
      allocate (this%bbb(ifms:ifme, jfms:jfme))
      allocate (this%phiwc(ifms:ifme, jfms:jfme))
      allocate (this%r_0(ifms:ifme, jfms:jfme))
#ifdef ESM_DUMP
      if (.not. allocated (e_speed)) then
        allocate (e_speed(ifms:ifme, jfms:jfme), e_tanphi(ifms:ifme, jfms:jfme), e_cor_wind(ifms:ifme, jfms:jfme), &
            e_cor_slope(ifms:ifme, jfms:jfme), e_umid(ifms:ifme, jfms:jfme), e_phiw(ifms:ifme, jfms:jfme), &
            e_phis(ifms:ifme, jfms:jfme), e_ros_base(ifms:ifme, jfms:jfme), e_ros_wind(ifms:ifme, jfms:jfme), &
            e_ros_slope(ifms:ifme, jfms:jfme))
        e_speed = 0.0; e_tanphi = 0.0; e_cor_wind = 0.0; e_cor_slope = 0.0; e_umid = 0.0; e_phiw = 0.0; e_phis = 0.0
        e_ros_base = 0.0; e_ros_wind = 0.0; e_ros_slope = 0.0
      end if
#endif

    end subroutine Init_ros_wrffire

#ifdef ESM_DUMP
    subroutine Esm_dump_ros_diag ()
      ! append the Calc_ros_wrffire intermediates to the currently open dump
      implicit none
      if (.not. esm_dump_active ()) return
      call esm_dump_var ('ros_speed', e_speed); call esm_dump_var ('ros_tanphi', e_tanphi)
      call esm_dump_var ('ros_cor_wind', e_cor_wind); call esm_dump_var ('ros_cor_slope', e_cor_slope)
      call esm_dump_var ('ros_umid', e_umid); call esm_dump_var ('ros_phiw', e_phiw); call esm_dump_var ('ros_phis', e_phis)
      call esm_dump_var ('ros_base', e_ros_base); call esm_dump_var ('ros_wind', e_ros_wind); call esm_dump_var ('ros_slope', e_ros_slope)
      call esm_dump_var ('ros_max', 6.0); call esm_dump_var ('ros_back_chap', 0.03333)
    end subroutine Esm_dump_ros_diag
#endif

    subroutine Set_ros_parameters_wrffire (this, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        fuels, nfuel_cat, fmc_g)

      implicit none

      class (ros_wrffire_t), intent (in out) :: this
      integer, intent(in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme
      class (fuel_t), intent (in) :: fuels
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat, fmc_g


      real ::  fuelload, fueldepth, rtemp1, rtemp2, qig, epsilon, rhob, wn, betaop, e, c, &
          xifr, etas, etam, a, gammax, gamma, ratio, ir, fuelloadm, bmst
      integer:: i, j, k, kk
      character (len = 128) :: msg
#ifdef ESM_DUMP
      real, dimension(ifms:ifme, jfms:jfme) :: e_bmst, e_fuelloadm, e_fuelload, e_fueldepth, e_betaop, e_qig, e_epsilon, &
          e_rhob, e_c, e_e, e_gammax, e_a, e_ratio, e_gamma, e_wn, e_rtemp1, e_etam, e_etas, e_ir, e_xifr
      e_bmst = 0.0; e_fuelloadm = 0.0; e_fuelload = 0.0; e_fueldepth = 0.0; e_betaop = 0.0; e_qig = 0.0; e_epsilon = 0.0
      e_rhob = 0.0; e_c = 0.0; e_e = 0.0; e_gammax = 0.0; e_a = 0.0; e_ratio = 0.0; e_gamma = 0.0; e_wn = 0.0
      e_rtemp1 = 0.0; e_etam = 0.0; e_etas = 0.0; e_ir = 0.0; e_xifr = 0.0
#endif


      Loop_j: do j = jfts, jfte
        Loop_i: do i = ifts, ifte
          k = int (nfuel_cat(i, j))
          if(k == fuels%no_fuel_cat) then
            this%ischap(i, j) = 0.0
              ! set to 1.0 to prevent grid%betafl(i,j)**(-0.3) to be Inf in fire_ros
            this%betafl(i, j) = 1.0
            this%bbb(i, j) = 1.0
            this%phiwc(i, j) = 0.0
            this%r_0(i, j) = 0.0
            this%iboros(i, j) = 0.0
          else
            this%ischap(i, j) = fuels%ichap(k)
              ! Settings of fire spread parameters from Rothermel
              ! No need to recalculate if FMC does not change
            bmst = fmc_g(i, j) / (1.0 + fmc_g(i, j))
              !  fuelload without moisture
            fuelloadm = (1.0 - bmst) * fuels%fgi(k)
            fuelload = fuelloadm * (0.3048) ** 2 * 2.205 ! to lb/ft^2
            fueldepth = fuels%fueldepthm(k) / 0.3048 ! to ft
              ! packing ratio
            this%betafl(i, j) = fuelload / (fueldepth * fuels%fueldens(k))
              ! optimum packing ratio
            betaop = 3.348 * fuels%savr(k) ** (-0.8189)
              ! heat of preignition, btu/lb
            qig = 250.0 + 1116.0 * fmc_g(i, j)
              ! effective heating number
            epsilon = exp (-138.0 / fuels%savr(k))
              ! ovendry bulk density, lb/ft^3
            rhob = fuelload/fueldepth

              ! const in wind coef
            c = 7.47 * exp (-0.133 * fuels%savr(k) ** 0.55)
            this%bbb(i,j) = 0.02526 * fuels%savr(k) ** 0.54
            e = 0.715 * exp (-3.59e-4 * fuels%savr(k))
            this%phiwc(i,j) = c * (this%betafl(i, j) / betaop) ** (-e)

            rtemp2 = fuels%savr(k) ** 1.5
              ! maximum rxn vel, 1/min
            gammax = rtemp2 / (495.0 + 0.0594 * rtemp2)
              ! coef for optimum rxn vel
            a = 1.0 / (4.774 * fuels%savr(k) ** 0.1 - 7.27)
            ratio = this%betafl(i,j)/betaop
              !optimum rxn vel, 1/min
            gamma = gammax * (ratio ** a) * exp(a * (1.0 - ratio))

             ! net fuel loading, lb/ft^2
            wn = fuelload/(1 + fuels%st(k))
            rtemp1 = fmc_g(i, j) / fuels%fuelmce(k)
              ! moist damp coef
            etam = 1.0 - 2.59 * rtemp1 + 5.11 * rtemp1 ** 2 - 3.52 * rtemp1 ** 3
              ! mineral damping coef
            etas = 0.174 * fuels%se(k) ** (-0.19)
              !rxn intensity,btu/ft^2 min
            ir = gamma * wn * FUELHEAT * etam * etas
            ! irm = ir * 1055./( 0.3048**2 * 60.) * 1.e-6     !for mw/m^2
            this%iboros(i, j) = ir * 1055.0 / ( 0.3048 ** 2 * 60.0) * 1.e-3 * (60.0 * 384.0 / fuels%savr(k)) ! I_R x t_r (kJ m^-2)
              ! propagating flux ratio
            xifr = exp((0.792 + 0.681 * fuels%savr(k) ** 0.5) &
                * (this%betafl(i, j) + 0.1)) / (192.0 + 0.2595 * fuels%savr(k))

              ! r_0 is the spread rate for a fire on flat ground with no wind.
              ! default spread rate in ft/min
            this%r_0(i, j) = ir * xifr / (rhob * epsilon * qig)
#ifdef ESM_DUMP
            e_bmst(i, j) = bmst; e_fuelloadm(i, j) = fuelloadm; e_fuelload(i, j) = fuelload; e_fueldepth(i, j) = fueldepth
            e_betaop(i, j) = betaop; e_qig(i, j) = qig; e_epsilon(i, j) = epsilon; e_rhob(i, j) = rhob; e_c(i, j) = c; e_e(i, j) = e
            e_gammax(i, j) = gammax; e_a(i, j) = a; e_ratio(i, j) = ratio; e_gamma(i, j) = gamma; e_wn(i, j) = wn
            e_rtemp1(i, j) = rtemp1; e_etam(i, j) = etam; e_etas(i, j) = etas; e_ir(i, j) = ir; e_xifr(i, j) = xifr
#endif
          end if
        end do Loop_i
      end do Loop_j
#ifdef ESM_DUMP
      if (esm_dump_want ('ros_params')) then
        call esm_dump_open ('ros_params')
        call esm_dump_var ('ifts', ifts); call esm_dump_var ('ifte', ifte); call esm_dump_var ('jfts', jfts); call esm_dump_var ('jfte', jfte)
        call esm_dump_var ('ifms', ifms); call esm_dump_var ('ifme', ifme); call esm_dump_var ('jfms', jfms); call esm_dump_var ('jfme', jfme)
        call esm_dump_var ('nfuel_cat', nfuel_cat); call esm_dump_var ('fmc_g', fmc_g)
        call esm_dump_var ('n_fuel_cat', fuels%n_fuel_cat); call esm_dump_var ('no_fuel_cat', fuels%no_fuel_cat)
        call esm_dump_var ('fgi', fuels%fgi); call esm_dump_var ('fueldepthm', fuels%fueldepthm); call esm_dump_var ('weight', fuels%weight)
        call esm_dump_var ('ichap', fuels%ichap); call esm_dump_var ('fueldens', fuels%fueldens); call esm_dump_var ('savr', fuels%savr)
        call esm_dump_var ('st', fuels%st); call esm_dump_var ('se', fuels%se); call esm_dump_var ('fuelmce', fuels%fuelmce)
        call esm_dump_var ('fgi_1h', fuels%fgi_1h); call esm_dump_var ('fgi_10h', fuels%fgi_10h); call esm_dump_var ('fgi_100h', fuels%fgi_100h)
        call esm_dump_var ('fgi_1000h', fuels%fgi_1000h); call esm_dump_var ('fgi_live', fuels%fgi_live); call esm_dump_var ('waf', fuels%waf)
        call esm_dump_var ('fuelheat', FUELHEAT); call esm_dump_var ('cmbcnst', CMBCNST)
        call esm_dump_var ('convert_j_per_kg_to_btu_per_pound', CONVERT_J_PER_KG_TO_BTU_PER_POUND)
        call esm_dump_var ('bmst', e_bmst); call esm_dump_var ('fuelloadm', e_fuelloadm); call esm_dump_var ('fuelload', e_fuelload)
        call esm_dump_var ('fueldepth', e_fueldepth); call esm_dump_var ('betaop', e_betaop); call esm_dump_var ('qig', e_qig)
        call esm_dump_var ('epsilon', e_epsilon); call esm_dump_var ('rhob', e_rhob); call esm_dump_var ('c', e_c); call esm_dump_var ('e', e_e)
        call esm_dump_var ('gammax', e_gammax); call esm_dump_var ('a', e_a); call esm_dump_var ('ratio', e_ratio); call esm_dump_var ('gamma', e_gamma)
        call esm_dump_var ('wn', e_wn); call esm_dump_var ('rtemp1', e_rtemp1); call esm_dump_var ('etam', e_etam); call esm_dump_var ('etas', e_etas)
        call esm_dump_var ('ir', e_ir); call esm_dump_var ('xifr', e_xifr)
        call esm_dump_var ('ischap', this%ischap); call esm_dump_var ('betafl', this%betafl); call esm_dump_var ('bbb', this%bbb)
        call esm_dump_var ('phiwc', this%phiwc); call esm_dump_var ('r_0', this%r_0); call esm_dump_var ('iboros', this%iboros)
        call esm_dump_close ()
      end if
#endif

    end subroutine Set_ros_parameters_wrffire

  end module ros_wrffire_mod
