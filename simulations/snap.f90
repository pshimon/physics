module snapio
  use iso_c_binding
  use iso_fortran_env, only: int32, real64
  implicit none
  interface
    function c_rename(old, new) bind(C, name="rename")
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: old(*), new(*)
      integer(c_int) :: c_rename
    end function
  end interface
contains
  ! File layout: int32 nx, int32 ny, float64 t, then u(nx,ny) float64
  ! (header = 16 bytes, data in Fortran column-major order)
  subroutine write_snapshot(step, t, u)
    integer, intent(in)       :: step
    real(real64), intent(in)  :: t, u(:,:)
    character(len=64) :: fin, tmp
    integer :: unit
    integer(c_int) :: rc

    write(fin,'(a,i6.6,a)') 'snap_', step, '.bin'
    tmp = trim(fin)//'.tmp'

    open(newunit=unit, file=trim(tmp), access='stream', &
         form='unformatted', status='replace', action='write')
    write(unit) int(size(u,1),int32), int(size(u,2),int32), t
    write(unit) u
    close(unit)

    rc = c_rename(trim(tmp)//c_null_char, trim(fin)//c_null_char)
  end subroutine
end module

program demo
  use snapio
  implicit none
  integer, parameter :: nx = 200, ny = 150
  real(real64) :: u(nx,ny), t
  integer :: step, i, j

  do step = 0, 500
    t = 0.01_real64*step
    do j = 1, ny
      do i = 1, nx
        u(i,j) = sin(0.05_real64*i - t) * cos(0.07_real64*j + 0.5_real64*t)
      end do
    end do
    if (mod(step,10) == 0) call write_snapshot(step, t, u)
  end do
end program
