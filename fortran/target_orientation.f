      subroutine read_target_orient(unit,method)
      implicit none
      integer unit,method,ios,k,n
      character(512) line
      method=1
      read(unit,'(A)',iostat=ios) line
      if(ios.lt.0) return
      if(ios.ne.0) then
        print *, 'Cannot read TARGET_ORIENT_METHOD'
        stop 1
      endif
      line=adjustl(line)
      if(len_trim(line).eq.0.or.line(1:1).eq.'!') return
      k=index(line,'!')
      if(k.gt.0) line(k:)=' '
      line=trim(adjustl(line))
      n=len_trim(line)
      if(n.eq.0) return
      k=1
      if(line(1:1).eq.'+'.or.line(1:1).eq.'-') k=2
      if(k.gt.n) then
        print *, 'Invalid TARGET_ORIENT_METHOD'
        stop 1
      endif
      if(verify(line(k:n),'0123456789').ne.0) then
        print *, 'Invalid TARGET_ORIENT_METHOD'
        stop 1
      endif
      read(line,*,iostat=ios) method
      if(ios.ne.0) then
        print *, 'Invalid TARGET_ORIENT_METHOD'
        stop 1
      endif
      if(method.ne.1.and.method.ne.2) then
        print *, 'Invalid TARGET_ORIENT_METHOD'
        stop 1
      endif
      end
