C Shared L-stage controller; see controller_dispatch.f.
      program ctrl_l
      implicit none
      call controller_dispatch('L')
      end
      include 'controller_dispatch.f'
