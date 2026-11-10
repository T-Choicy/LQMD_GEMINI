C Shared A-stage controller; see controller_dispatch.f.
      program ctrl_a
      implicit none
      call controller_dispatch('A')
      end
      include 'controller_dispatch.f'
