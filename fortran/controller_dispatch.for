C Shared entry for both native and Windows controller names.
      subroutine controller_dispatch(stage)
      implicit none
      character*(*) stage
      character(2048) config,arg,root,base,exe,probe,osname
      character(8192) command
      character(256) message
      character(16) action
      character(24) shell
      integer narg,i,j,k,ios,exitcode,cmdstat
      logical windows,found
      config=''
      action=''
      osname=''
      call get_environment_variable('OS',osname,status=ios)
      windows=trim(osname) == 'Windows_NT'
      if(.not.windows)then
        call get_environment_variable('SystemRoot',osname,status=ios)
        windows=ios == 0.and.len_trim(osname) > 0
      endif
      shell='pwsh'
      if(windows) shell='powershell.exe'
      narg=command_argument_count()
      if(narg > 2) error stop 'Use a config path and at most one option'
      do i=1,narg
        call get_command_argument(i,arg,status=ios)
        if(ios /= 0) error stop 'Controller argument exceeds path limit'
        select case(trim(arg))
        case('--help','-h')
          print *, 'Controller ',stage,': [config] [--prepare-only',
     &      ' | --merge-only]'
          print *, 'Reads original INPUT.in/INPUT_ANA.in or JSON.'
          print *, 'Uses shared backend, isolated shards and 25 slots.'
          print *, 'Linux requires PowerShell 7 (pwsh) and compiled',
     &      ' LQMD/exfrag; Windows uses powershell.exe.'
          return
        case('--prepare-only','-PrepareOnly')
          if(len_trim(action) /= 0) error stop 'Duplicate action'
          action=' -PrepareOnly'
        case('--merge-only','-MergeOnly')
          if(len_trim(action) /= 0) error stop 'Duplicate action'
          action=' -MergeOnly'
        case default
          if(len_trim(config) /= 0) error stop 'Duplicate config'
          config=arg
        end select
      enddo
      call get_command_argument(0,exe,status=ios)
      if(ios /= 0) error stop 'Executable path exceeds limit'
      root=''
      found=.false.
C Search the working directory first, then ancestors of the executable.
      do j=1,2
        base='.'
        if(j == 2)then
          base=exe
          k=scan(trim(base),'/'//achar(92),back=.true.)
          if(k == 0) cycle
          base=base(:k)
        endif
        do i=0,12
          probe=trim(base)//'/scripts/windows_event_controller.ps1'
          inquire(file=trim(probe),exist=found)
          if(found)then
            root=base
            exit
          endif
          base=trim(base)//'/..'
        enddo
        if(found) exit
      enddo
      if(.not.found) error stop 'Run controller within project tree'
      if(len_trim(config) == 0)then
        arg='INPUT.in'
        if(stage == 'A') arg='INPUT_ANA.in'
        config='Input file/'//trim(arg)
        inquire(file=trim(config),exist=found)
        if(.not.found)
     &    config=trim(root)//'/fortran/Input file/'//trim(arg)
      endif
      inquire(file=trim(config),exist=found)
      if(.not.found) error stop 'Controller input file not found'
      probe=trim(root)//'/scripts/windows_event_controller.ps1'
      command=trim(shell)//' -NoLogo -NoProfile -ExecutionPolicy Bypass'
     &  //' -File '//trim(controller_quote(trim(probe),windows))
     &  //' -Stage '//stage//' -Config '
     &  //trim(controller_quote(trim(config),windows))//trim(action)
      call execute_command_line(trim(command),exitstat=exitcode,
     &  cmdstat=cmdstat,cmdmsg=message)
      if(cmdstat /= 0)then
        print *, 'Cannot start ',trim(shell),': ',trim(message)
        error stop 2
      endif
      if(exitcode /= 0)then
        print *, 'Controller backend exit code: ',exitcode
        error stop 1
      endif
      contains
      function controller_quote(value,is_windows) result(quoted)
      character(*) value
      logical is_windows
      character(2052) quoted
      if(index(value,achar(10)) > 0.or.index(value,achar(13)) > 0)
     &  error stop 'Newline in controller argument'
      if(is_windows)then
        if(index(value,'"') > 0.or.index(value,'%') > 0)
     &    error stop 'Unsafe quote or percent in controller path'
        quoted='"'//value//'"'
      else
        if(index(value,achar(39)) > 0)
     &    error stop 'Single quotes in controller paths unsupported'
        quoted=achar(39)//value//achar(39)
      endif
      end function
      end subroutine
