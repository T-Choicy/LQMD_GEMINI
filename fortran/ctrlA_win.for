C Shared A-stage controller; see controller_dispatch.for.
      program ctrl_a
      implicit none
      call controller_dispatch('A')
      end
      include 'controller_dispatch.for'
