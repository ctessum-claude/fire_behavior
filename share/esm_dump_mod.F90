!=======================================================================
! module_esm_dump: structured (JSON) dumps of physics-kernel inputs and
! outputs for the EarthSciML EqWeFiC project (branch earthsciml-instrumented).
!
! Usage in a physics wrapper, around a kernel call:
!   use module_esm_dump
!   ...
!   dumping = esm_dump_want('ysu')          ! counts calls per scheme
!   if (dumping) call esm_dump_open('ysu')
!   if (dumping) call esm_dump_var('tx', t3d_hv)   ! any rank 0-3 real/int/logical
!   call bl_ysu_run(...)
!   if (dumping) call esm_dump_var('ttnp', rthblten_hv)
!   if (dumping) call esm_dump_close()
!
! Sub-process diagnostics inside a kernel (guarded so the kernel stays clean
! when built without the dump module, e.g. #ifdef ESM_DUMP in physics_mmm):
!   if (esm_dump_inner) call esm_dump_var('praut_raw', praut)
! The wrapper sets esm_dump_inner = dumping (or dumping .and. i == its for
! per-column kernels) around the call; esm_dump_active() tells whether a
! dump file is currently open.
!
! Control (environment variables):
!   ESM_DUMP       comma-separated scheme names to dump, or "all" (unset = off)
!   ESM_DUMP_CALLS comma-separated 1-based call indices per scheme to dump
!                  (default "1"); every call increments the scheme's counter.
!   ESM_DUMP_DIR   output directory (default ".")
! One file per dumped call: <dir>/esm_dump_<scheme>_<call>.json with
!   {"scheme": "...", "call": n, "vars": {"name": {"kind": ..., "shape": [..],
!    "data": [...]}, ...}}
! Arrays are written in Fortran (column-major) order, flattened; "shape"
! gives the extents so an extractor can reshape. Reals use ES24.16E3.
!=======================================================================
module module_esm_dump
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
  public :: esm_dump_want, esm_dump_open, esm_dump_open_file, esm_dump_close, esm_dump_var
  public :: esm_dump_active, esm_dump_inner

  ! esm_dump_inner: set by a wrapper (or driver) to let dump calls placed
  ! INSIDE a kernel (sub-process diagnostics) append to the open file.  The
  ! wrapper sets it .true. only for the call/column it wants, so kernels that
  ! are called once per column do not write duplicate keys.
  logical, save :: esm_dump_inner = .false.

  integer, parameter :: max_schemes = 32
  character(len=32) :: scheme_names(max_schemes) = ''
  integer :: scheme_calls(max_schemes) = 0
  integer :: nschemes = 0
  integer :: unit = -1
  logical :: first_var = .true.
  character(len=32) :: cur_scheme = ''
  integer :: cur_call = 0

  interface esm_dump_var
    module procedure dump_r8_0, dump_r8_1, dump_r8_2, dump_r8_3
    module procedure dump_r4_0, dump_r4_1, dump_r4_2, dump_r4_3
    module procedure dump_i_0, dump_i_1, dump_i_2
    module procedure dump_l_0, dump_l_1
    module procedure dump_str
  end interface

contains

  logical function env_list_has(varname, item, default) result(has)
    character(len=*), intent(in) :: varname, item, default
    character(len=1024) :: val
    integer :: stat, p, q
    call get_environment_variable(varname, val, status=stat)
    if (stat /= 0) val = default
    has = .false.
    if (len_trim(val) == 0) return
    if (trim(val) == 'all') then
      has = .true.; return
    end if
    p = 1
    do
      q = index(val(p:), ',')
      if (q == 0) then
        if (trim(adjustl(val(p:))) == trim(item)) has = .true.
        exit
      end if
      if (trim(adjustl(val(p:p+q-2))) == trim(item)) has = .true.
      p = p + q
      if (p > len_trim(val)) exit
    end do
  end function env_list_has

  logical function esm_dump_want(scheme) result(want)
    character(len=*), intent(in) :: scheme
    integer :: is
    character(len=16) :: callstr
    want = .false.
    if (.not. env_list_has('ESM_DUMP', scheme, '')) return
    is = 0
    do is = 1, nschemes
      if (trim(scheme_names(is)) == trim(scheme)) exit
    end do
    if (is > nschemes) then
      nschemes = nschemes + 1
      is = nschemes
      scheme_names(is) = scheme
    end if
    scheme_calls(is) = scheme_calls(is) + 1
    write(callstr, '(I0)') scheme_calls(is)
    want = env_list_has('ESM_DUMP_CALLS', trim(callstr), '1')
    if (want) then
      cur_scheme = scheme
      cur_call = scheme_calls(is)
    end if
  end function esm_dump_want

  subroutine esm_dump_open(scheme)
    character(len=*), intent(in) :: scheme
    character(len=1024) :: dir, fname
    integer :: stat
    call get_environment_variable('ESM_DUMP_DIR', dir, status=stat)
    if (stat /= 0 .or. len_trim(dir) == 0) dir = '.'
    write(fname, '(A,"/esm_dump_",A,"_",I0,".json")') trim(dir), trim(scheme), cur_call
    open(newunit=unit, file=trim(fname), status='replace', action='write')
    write(unit, '(A,A,A,I0,A)') '{"scheme": "', trim(scheme), '", "call": ', cur_call, ', "vars": {'
    first_var = .true.
  end subroutine esm_dump_open

  ! Open a dump at an explicit path (used by the standalone kernel drivers).
  subroutine esm_dump_open_file(scheme, path)
    character(len=*), intent(in) :: scheme, path
    open(newunit=unit, file=trim(path), status='replace', action='write')
    write(unit, '(A,A,A,I0,A)') '{"scheme": "', trim(scheme), '", "call": ', 0, ', "vars": {'
    first_var = .true.
  end subroutine esm_dump_open_file

  subroutine esm_dump_close()
    write(unit, '(A)') '}}'
    close(unit)
    unit = -1
    esm_dump_inner = .false.
  end subroutine esm_dump_close

  logical function esm_dump_active() result(active)
    active = unit /= -1
  end function esm_dump_active

  subroutine header(name, kind, shape)
    character(len=*), intent(in) :: name, kind
    integer, intent(in) :: shape(:)
    integer :: d
    if (.not. first_var) write(unit, '(A)') ','
    first_var = .false.
    write(unit, '(A,A,A,A,A)', advance='no') '"', trim(name), '": {"kind": "', kind, '", "shape": ['
    do d = 1, size(shape)
      if (d > 1) write(unit, '(A)', advance='no') ', '
      write(unit, '(I0)', advance='no') shape(d)
    end do
    write(unit, '(A)', advance='no') '], "data": ['
  end subroutine header

  subroutine footer()
    write(unit, '(A)', advance='no') ']}'
  end subroutine footer

  subroutine write_r8(vals)
    real(real64), intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      if (mod(n, 4) == 1 .and. n > 1) write(unit, '(A)') ''
      write(unit, '(ES24.16E3)', advance='no') vals(n)
    end do
  end subroutine write_r8

  subroutine dump_r8_0(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v
    call header(name, 'real64', [integer ::])
    call write_r8([v])
    call footer()
  end subroutine dump_r8_0
  subroutine dump_r8_1(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:)
    call header(name, 'real64', shape(v))
    call write_r8(v)
    call footer()
  end subroutine dump_r8_1
  subroutine dump_r8_2(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:,:)
    call header(name, 'real64', shape(v))
    call write_r8(reshape(v, [size(v)]))
    call footer()
  end subroutine dump_r8_2
  subroutine dump_r8_3(name, v)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: v(:,:,:)
    call header(name, 'real64', shape(v))
    call write_r8(reshape(v, [size(v)]))
    call footer()
  end subroutine dump_r8_3

  subroutine dump_r4_0(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v
    call header(name, 'real32', [integer ::])
    call write_r8([real(v, real64)])
    call footer()
  end subroutine dump_r4_0
  subroutine dump_r4_1(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:)
    call header(name, 'real32', shape(v))
    call write_r8(real(v, real64))
    call footer()
  end subroutine dump_r4_1
  subroutine dump_r4_2(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:,:)
    call header(name, 'real32', shape(v))
    call write_r8(real(reshape(v, [size(v)]), real64))
    call footer()
  end subroutine dump_r4_2
  subroutine dump_r4_3(name, v)
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: v(:,:,:)
    call header(name, 'real32', shape(v))
    call write_r8(real(reshape(v, [size(v)]), real64))
    call footer()
  end subroutine dump_r4_3

  subroutine write_i(vals)
    integer, intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      write(unit, '(I0)', advance='no') vals(n)
    end do
  end subroutine write_i
  subroutine dump_i_0(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v
    call header(name, 'int', [integer ::])
    call write_i([v])
    call footer()
  end subroutine dump_i_0
  subroutine dump_i_1(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:)
    call header(name, 'int', shape(v))
    call write_i(v)
    call footer()
  end subroutine dump_i_1
  subroutine dump_i_2(name, v)
    character(len=*), intent(in) :: name
    integer, intent(in) :: v(:,:)
    call header(name, 'int', shape(v))
    call write_i(reshape(v, [size(v)]))
    call footer()
  end subroutine dump_i_2

  subroutine write_l(vals)
    logical, intent(in) :: vals(:)
    integer :: n
    do n = 1, size(vals)
      if (n > 1) write(unit, '(A)', advance='no') ', '
      if (vals(n)) then
        write(unit, '(A)', advance='no') 'true'
      else
        write(unit, '(A)', advance='no') 'false'
      end if
    end do
  end subroutine write_l
  subroutine dump_l_0(name, v)
    character(len=*), intent(in) :: name
    logical, intent(in) :: v
    call header(name, 'bool', [integer ::])
    call write_l([v])
    call footer()
  end subroutine dump_l_0
  subroutine dump_l_1(name, v)
    character(len=*), intent(in) :: name
    logical, intent(in) :: v(:)
    call header(name, 'bool', shape(v))
    call write_l(v)
    call footer()
  end subroutine dump_l_1

  subroutine dump_str(name, v)
    character(len=*), intent(in) :: name, v
    if (.not. first_var) write(unit, '(A)') ','
    first_var = .false.
    write(unit, '(A,A,A,A,A)', advance='no') '"', trim(name), '": {"kind": "str", "shape": [], "data": ["', trim(v), '"]}'
  end subroutine dump_str

end module module_esm_dump
