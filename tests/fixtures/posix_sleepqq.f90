! Test-only bridge for the existing mock fixture's millisecond sleep.
subroutine sleepqq(milliseconds)
  use iso_c_binding, only: c_int
  implicit none
  integer, intent(in) :: milliseconds
  integer(c_int) :: ignored
  interface
    function usleep(microseconds) bind(C,name='usleep') result(status)
      import c_int
      integer(c_int), value :: microseconds
      integer(c_int) :: status
    end function
  end interface
  ignored = usleep(int(milliseconds,c_int)*1000_c_int)
end subroutine
