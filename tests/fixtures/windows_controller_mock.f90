! Mock LQMD/exfrag fixture for software-only Windows controller tests.
! It does not model any nuclear physics or produce physics-valid data.
program windows_controller_mock
  use iso_fortran_env, only: error_unit, output_unit
  implicit none
  character(len=1024) :: arg0, cwd
  integer :: ios, nrmin, nrmax, nth_run, nth_noh, nevents
  integer :: model, seed, shard, event
  logical :: legacy_mode

  call get_command_argument(0, arg0)
  arg0 = arg0(scan(trim(arg0), '/'//achar(92), back=.true.) + 1:)
  if (index(lower(trim(arg0)), 'lqmd') > 0) then
    call run_lqmd_mock()
  else
    call run_exfrag_mock()
  end if

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(error_unit, '(A)') 'MOCK_FIXTURE_ERROR: '//trim(message)
    stop 7
  end subroutine fail

  subroutine read_line(path, number, value)
    character(len=*), intent(in) :: path
    integer, intent(in) :: number
    character(len=*), intent(out) :: value
    integer :: local_unit, local_ios, row
    open(newunit=local_unit, file=path, status='old', action='read', iostat=local_ios)
    if (local_ios /= 0) call fail('cannot open '//trim(path))
    do row = 1, number
      read(local_unit, '(A)', iostat=local_ios) value
      if (local_ios /= 0) call fail('cannot read requested line from '//trim(path))
    end do
    close(local_unit)
  end subroutine read_line

  subroutine read_selector()
    integer :: local_unit, local_ios, row, comment_pos
    character(len=1024) :: value
    model = -1
    seed = -1
    open(newunit=local_unit, file='Input file/m_select.in', status='old', &
         action='read', iostat=local_ios)
    if (local_ios /= 0) call fail('cannot open Input file/m_select.in')
    row = 0
    do
      read(local_unit, '(A)', iostat=local_ios) value
      if (local_ios /= 0) exit
      if (len_trim(value) == 0 .or. value(1:1) == '!' .or. value(1:1) == '#') cycle
      row = row + 1
      comment_pos = index(value, '!')
      if (comment_pos > 0) value = value(1:comment_pos - 1)
      comment_pos = index(value, '#')
      if (comment_pos > 0) value = value(1:comment_pos - 1)
      if (row == 1) read(value, *, iostat=local_ios) model
      if (row == 3) read(value, *, iostat=local_ios) seed
      if (row >= 3) exit
    end do
    close(local_unit)
    if (model < 0 .or. seed < 0) call fail('invalid m_select.in')
  end subroutine read_selector

  integer function shard_number()
    integer :: pos, local_ios, n
    call getcwd(cwd)
    pos = index(cwd, 'shard-')
    if (pos <= 0) call fail('working directory has no shard-NNNN component')
    read(cwd(pos + 6:pos + 9), '(I4)', iostat=local_ios) n
    if (local_ios /= 0) call fail('invalid shard directory name')
    shard_number = n
  end function shard_number

  subroutine run_lqmd_mock()
    integer :: out_coll, out_angm, input_unit, nfile, events_per_file
    integer :: file_index, source_model
    logical :: input_exists
    character(len=1024) :: control, source_model_text
    character(len=16) :: allow_selected
    call read_line('Input file/LQMD.in', 1, control)
    read(control, *, iostat=ios) nrmin, nrmax, nth_run, nth_noh
    if (ios /= 0 .or. nrmin < 1 .or. nrmax < nrmin) call fail('invalid LQMD event range')
    call read_selector()
    inquire(file='Input file/INPUT.in', exist=input_exists)
    legacy_mode = input_exists
    if (seed <= 0) call fail('LQMD requires a positive comparison seed')
    if (legacy_mode) then
      call get_environment_variable('MOCK_ALLOW_SELECTED_LQMD', allow_selected, status=ios)
      if (ios /= 0 .or. trim(allow_selected) /= '1') &
        call fail('legacy selected-mass LQMD requires test-only MOCK_ALLOW_SELECTED_LQMD=1')
      call read_line('Input file/INPUT.in', 1, control)
      read(control, *, iostat=ios) nfile
      if (ios /= 0 .or. nfile /= 5) call fail('legacy INPUT.in NFILE must be five')
      call read_line('Input file/INPUT.in', 2, control)
      read(control, *, iostat=ios) events_per_file
      if (ios /= 0 .or. events_per_file /= 2) call fail('legacy INPUT.in NRUN_PER must be two in this software test')
      file_index = shard_number()
      if (nrmin /= (file_index - 1) * events_per_file + 1 .or. &
          nrmax /= file_index * events_per_file) call fail('legacy LQMD.in range is not globally numbered')
      call read_line('Input file/INPUT_ANA.in', 1, control)
      if (model < 1 .or. model > 5) call fail('legacy LQMD must use selected mass model 1..5')
      open(newunit=input_unit, file='Input file/m_select.in', status='old', action='read', iostat=ios)
      if (ios /= 0) call fail('legacy model-specific m_select.in is missing')
      read(input_unit, '(A)', iostat=ios) source_model_text
      if (ios /= 0) call fail('cannot read first m_select.in data row')
      read(source_model_text, *, iostat=ios) source_model
      close(input_unit)
      if (ios /= 0 .or. source_model /= model) call fail('m_select.in model does not match this model-specific source')
    else
      if (model /= 0) call fail('legacy JSON-style LQMD must keep m_select=0')
    end if

    open(newunit=out_coll, file='coll.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fail('cannot create coll.dat')
    if (nrmin == 1) then
      write(out_coll, '(A)') ''
      write(out_coll, '(A)') 'MOCK SOFTWARE TEST COLLISION HEADER'
      write(out_coll, '(A)') '0 0'
    end if
    open(newunit=out_angm, file='initangm.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) call fail('cannot create initangm.dat')
    do event = nrmin, nrmax
      write(out_coll, '(A,I0)') 'mock-event ', event
      write(out_angm, '(I0,1X,I0,1X,I0)') event, seed, model
    end do
    close(out_coll)
    close(out_angm)
    open(unit=72, file='mock-started.txt', status='replace', action='write')
    write(72, '(A)') 'SOFTWARE TEST FIXTURE ONLY'
    close(72)
    if (legacy_mode) then
      open(unit=73, file='mock-provenance.txt', status='replace', action='write')
      write(73, '(4(I0,1X))') model, seed, nrmin, nrmax
      close(73)
    end if
    call wait_for_release()
    write(output_unit, '(A)') 'The program end normally, ok!'
  end subroutine run_lqmd_mock

  subroutine wait_for_release()
    character(len=1024) :: release_path
    character(len=32) :: sleep_text
    logical :: release_exists
    integer :: local_status, tick, sleep_ms
    call get_environment_variable('MOCK_RELEASE_FILE', release_path, status=local_status)
    if (local_status /= 0 .or. len_trim(release_path) == 0) then
      call get_environment_variable('MOCK_SLEEP_MS', sleep_text, status=local_status)
      if (local_status == 0) then
        read(sleep_text, *, iostat=local_status) sleep_ms
        if (local_status == 0 .and. sleep_ms > 0) call sleepqq(sleep_ms)
      end if
      return
    end if
    do tick = 1, 1200
      inquire(file=trim(release_path), exist=release_exists)
      if (release_exists) return
      call sleepqq(50)
    end do
    call fail('timed out waiting for test-only MOCK_RELEASE_FILE barrier')
  end subroutine wait_for_release

  subroutine find_analysis_events()
    integer :: local_unit, local_ios
    character(len=1024) :: value
    nevents = -1
    open(newunit=local_unit, file='Input file/anahyp.in', status='old', &
         action='read', iostat=local_ios)
    if (local_ios /= 0) call fail('cannot open Input file/anahyp.in')
    do
      read(local_unit, '(A)', iostat=local_ios) value
      if (local_ios /= 0) exit
      if (trim(adjustl(value)) == '%%Simulation') then
        read(local_unit, '(A)', iostat=local_ios) value
        if (local_ios /= 0) call fail('missing Simulation event row')
        read(value, *, iostat=local_ios) nevents
        exit
      end if
    end do
    close(local_unit)
    if (nevents < 1) call fail('cannot find positive events value in anahyp.in')
  end subroutine find_analysis_events

  subroutine write_cnpcs(path, heading, base)
    character(len=*), intent(in) :: path, heading
    real(8), intent(in) :: base
    integer :: out, local_ios, block, row
    character(len=64), dimension(4) :: headings
    headings = [character(len=64) :: &
      ' NUCLIDE CROSS SECTION AS A FUNCTION OF Z', &
      ' NUCLIDE CROSS SECTION AS A FUNCTION OF N', &
      ' NUCLIDE CROSS SECTION AS A FUNCTION OF A', &
      ' CROSS SECTION OF ISOTOPICS']
    open(newunit=out, file=path, status='replace', action='write', iostat=local_ios)
    if (local_ios /= 0) call fail('cannot create '//trim(path))
    do block = 1, 4
      write(out, '(A)') trim(headings(block))
      if (block == 4) then
        do row = 1, 5
          write(out, '(I12,F14.7)') row, base * dble(row)
        end do
      else
        write(out, '(I12,F14.7)') 1, base
        if (shard == 2) write(out, '(I12,F14.7)') 2, base * 3.0d0
      end if
      write(out, '(A)') ' '
    end do
    close(out)
  end subroutine write_cnpcs

  subroutine run_exfrag_mock()
    integer :: slash, local_unit, local_ios, event_rows, angle_rows, directory_model, row_event, row_seed, row_model
    logical :: legacy_analysis
    character(len=1024) :: config_line
    real(8) :: part_base
    call getcwd(cwd)
    inquire(file='Input file/INPUT.in', exist=legacy_analysis)
    directory_model = -1
    if (.not. legacy_analysis) then
      slash = index(cwd, 'model-')
      if (slash <= 0) call fail('working directory has no model-N component')
      read(cwd(slash + 6:slash + 6), '(I1)', iostat=ios) directory_model
      if (ios /= 0 .or. directory_model < 1 .or. directory_model > 5) call fail('invalid model directory')
    end if
    shard = shard_number()
    call read_selector()
    if (model < 1 .or. model > 5) call fail('analysis model outside 1..5')
    if (.not. legacy_analysis .and. model /= directory_model) &
      call fail('analysis mass selection disagrees with model-specific input directory')
    call find_analysis_events()
    open(unit=71, file='coll.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('analysis input coll.dat is missing')
    read(71, '(A)', iostat=ios) config_line
    if (ios /= 0 .or. len_trim(config_line) /= 0) call fail('analysis coll.dat lacks initial blank header record')
    read(71, '(A)', iostat=ios) config_line
    if (ios /= 0 .or. index(config_line, 'MOCK SOFTWARE TEST COLLISION HEADER') <= 0) &
      call fail('analysis coll.dat has invalid or missing header')
    event_rows = 0
    do
      read(71, '(A)', iostat=local_ios) config_line
      if (local_ios /= 0) exit
      if (index(config_line, 'mock-event ') > 0) event_rows = event_rows + 1
    end do
    close(71)
    if (event_rows /= nevents) call fail('coll.dat event count does not match analysis input')
    open(newunit=local_unit, file='initangm.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('analysis input initangm.dat is missing')
    angle_rows = 0
    do
      read(local_unit, *, iostat=local_ios) row_event, row_seed, row_model
      if (local_ios /= 0) exit
      if (legacy_analysis .and. row_model /= model) &
        call fail('analysis m_select.in model does not match the LQMD event source')
      angle_rows = angle_rows + 1
    end do
    close(local_unit)
    if (angle_rows /= nevents) call fail('initangm.dat event count does not match analysis input')
    part_base = dble(model * 10 + shard * 5)
    call write_cnpcs('CNPCS_HOT.dat', 'HOT', part_base + 1000.0d0)
    call write_cnpcs('CNPCS_COLD.dat', 'COLD', part_base)
    write(output_unit, '(A)') 'OK, the program ends normally !!'
  end subroutine run_exfrag_mock

  function lower(value) result(out)
    character(len=*), intent(in) :: value
    character(len=len(value)) :: out
    integer :: j, c
    out = value
    do j = 1, len_trim(out)
      c = iachar(out(j:j))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(j:j) = achar(c + 32)
    end do
  end function lower

end program windows_controller_mock
