C     Production mass-selection guard and inspection-only table readers.
      MODULE mass_selection
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: m_select, read_mass_selection, load_selected_table
      PUBLIC :: lookup_selected_mass, overlay_selected_masses
      PUBLIC :: inspect_selected_table, inspect_selected_mass
      PUBLIC :: require_physical_compatibility
      PUBLIC :: lookup_selected_frdm12, validate_selected_source
      PUBLIC :: native_mass_count, fallback_mass_count
      PUBLIC :: frdm12_native_count
      PUBLIC :: mass_framework, comparison_seed
      PUBLIC :: report_mass_comparison_usage
      INTEGER, SAVE :: m_select = 0
      INTEGER, SAVE :: mass_framework = 0
      INTEGER, SAVE :: comparison_seed = 314159
      INTEGER, PARAMETER :: max_z = 1000, max_a = 1000
      REAL, SAVE :: values(0:max_z,0:max_a)
      LOGICAL, SAVE :: present_mass(0:max_z,0:max_a) = .FALSE.
      LOGICAL, SAVE :: native_mass(0:max_z,0:max_a) = .FALSE.
      REAL, SAVE :: frdm_mth(0:max_z,0:max_a)
      REAL, SAVE :: frdm_emic(0:max_z,0:max_a)
      REAL, SAVE :: frdm_beta2(0:max_z,0:max_a)
      LOGICAL, SAVE :: frdm_present(0:max_z,0:max_a) = .FALSE.
      LOGICAL, SAVE :: frdm_seen(0:max_z,0:max_a) = .FALSE.
      LOGICAL, SAVE :: frdm_loaded = .FALSE.
      CHARACTER(LEN=512), SAVE :: frdm_loaded_dir = ''
      INTEGER, SAVE :: frdm12_native_count = 0
      INTEGER, SAVE :: inspect_loaded_select = -1
      INTEGER, SAVE :: native_mass_count = 0, fallback_mass_count = 0
      CHARACTER(LEN=512), SAVE :: inspect_loaded_dir = ''
      REAL, SAVE :: comparison_mass(0:max_z,0:max_a)
      LOGICAL, SAVE :: comparison_present(0:max_z,0:max_a) = .FALSE.
      LOGICAL, SAVE :: comparison_native(0:max_z,0:max_a) = .FALSE.
      LOGICAL, SAVE :: comparison_fallback_reported(0:max_z,0:max_a)
      LOGICAL, SAVE :: comparison_usage_reported(0:max_z,0:max_a)
      LOGICAL, SAVE :: comparison_loaded = .FALSE.
      CHARACTER(LEN=512), SAVE :: comparison_loaded_dir = ''
      INTEGER, SAVE :: comparison_native_count = 0
      INTEGER, SAVE :: comparison_fallback_count = 0
      INTEGER, SAVE :: comparison_ignored_light_count = 0
      INTEGER, SAVE :: comparison_fallback_lookups = 0
      INTEGER, SAVE :: comparison_fallback_detail_count = 0
      INTEGER, SAVE :: comparison_loaded_model = -1
      INTEGER, SAVE :: comparison_shared_light_lookups = 0
      INTEGER, SAVE :: comparison_native_unique = 0
      INTEGER, SAVE :: comparison_fallback_unique = 0
      INTEGER, SAVE :: comparison_shared_light_unique = 0
      INTEGER, SAVE :: comparison_usage_unique = 0
      INTEGER, SAVE :: comparison_query_calls = 0
      CONTAINS

      SUBROUTINE mass_error(message,path,line_no)
      CHARACTER(LEN=*), INTENT(IN) :: message,path
      INTEGER, INTENT(IN), OPTIONAL :: line_no
      WRITE(error_unit,'(A)') 'Mass selection: '//TRIM(message)
      WRITE(error_unit,'(A)') TRIM(path)
      IF (PRESENT(line_no)) WRITE(error_unit,'(A,I0)')
     &  'Record line: ',line_no
      STOP 1
      END SUBROUTINE mass_error

      SUBROUTINE require_physical_compatibility()
      CHARACTER(LEN=16) :: model
      CHARACTER(LEN=256) :: missing
      IF (mass_framework == 1) THEN
        IF (m_select >= 0 .AND. m_select <= 5) RETURN
      ENDIF
      IF (m_select == 0) RETURN
      SELECT CASE(m_select)
      CASE(1)
        model = 'FRDM12'
        RETURN
      CASE(2)
        model = 'HFB02'
        missing = 'validated smooth/effective-microscopic decomposition'//
     &    ' compatible with Mth'
      CASE(3)
        model = 'KTUY05'
        missing = 'KTUY05 mean even-odd term Meo: values or complete'//
     &    ' final-version formula and parameters; Esh is supplied'
      CASE(4)
        model = 'WSM14'
        missing = 'gross/smooth-mass, pairing and microscopic'//
     &    ' decomposition plus RBF correction assignment'
      CASE(5)
        model = 'LDMM'
        missing = 'verified provenance and column-5 mass/energy'//
     &    ' definition; column 4 matches old Emic numerically'
      CASE DEFAULT
        model = 'INVALID'
        missing = 'valid m_select in 0..5'
      END SELECT
      WRITE(error_unit,'(A,I0,A,A,A,A,A)')
     &  'MASS_SELECT_PHYSICS_UNSUPPORTED m_select=',m_select,
     &  ' model=',TRIM(model),' missing=',TRIM(missing)
      STOP 1
      END SUBROUTINE require_physical_compatibility

      SUBROUTINE require_inspection_mode()
      CHARACTER(LEN=16) :: model
      IF (m_select == 0) RETURN
      SELECT CASE(m_select)
      CASE(1)
        model = 'FRDM12'
      CASE(2)
        model = 'HFB02'
      CASE(3)
        model = 'KTUY05'
      CASE(4)
        model = 'WSM14'
      CASE(5)
        model = 'LDMM'
      CASE DEFAULT
        model = 'INVALID'
      END SELECT
      WRITE(error_unit,'(A,I0,A,A)')
     &  'MASS_SELECT_INSPECTION_REQUIRES_DEFAULT m_select=',m_select,
     &  ' model=',TRIM(model)
      STOP 1
      END SUBROUTINE require_inspection_mode

      SUBROUTINE read_mass_selection(config_path)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: config_path
      CHARACTER(LEN=512) :: path,line
      INTEGER :: unit_no,ios,parsed,framework_value,seed_value
      LOGICAL :: exists,framework_present
      path = 'Input file/m_select.in'
      IF (PRESENT(config_path)) path = config_path
      m_select = 0
      mass_framework = 0
      comparison_seed = 314159
      INQUIRE(FILE=TRIM(path),EXIST=exists)
      IF (.NOT.exists) RETURN
      OPEN(NEWUNIT=unit_no,FILE=TRIM(path),STATUS='OLD',
     &  ACTION='READ',IOSTAT=ios)
      IF (ios /= 0) CALL mass_error('Cannot open selector',path)
      DO
        READ(unit_no,'(A)',IOSTAT=ios) line
        IF (ios /= 0) CALL mass_error('Missing m_select value',path)
        line = ADJUSTL(line)
        IF (LEN_TRIM(line) == 0) CYCLE
        IF (line(1:1) == '!' .OR. line(1:1) == '#') CYCLE
        READ(line,*,IOSTAT=ios) parsed
        IF (ios /= 0) CALL mass_error('Invalid m_select value',path)
        EXIT
      ENDDO
      framework_value = 0
      seed_value = 314159
      framework_present = .FALSE.
      DO
        READ(unit_no,'(A)',IOSTAT=ios) line
        IF (ios < 0) EXIT
        IF (ios /= 0) CALL mass_error('Cannot read framework value',path)
        line = ADJUSTL(line)
        IF (LEN_TRIM(line) == 0) CYCLE
        IF (line(1:1) == '!' .OR. line(1:1) == '#') CYCLE
        READ(line,*,IOSTAT=ios) framework_value
        IF (ios /= 0) CALL mass_error('Invalid mass_framework',path)
        framework_present = .TRUE.
        EXIT
      ENDDO
      IF (framework_value < 0 .OR. framework_value > 1)
     &  CALL mass_error('mass_framework must be 0 or 1',path)
      IF (framework_present) THEN
        DO
          READ(unit_no,'(A)',IOSTAT=ios) line
          IF (ios < 0) EXIT
          IF (ios /= 0)
     &        CALL mass_error('Cannot read comparison seed',path)
          line = ADJUSTL(line)
          IF (LEN_TRIM(line) == 0) CYCLE
          IF (line(1:1) == '!' .OR. line(1:1) == '#') CYCLE
          READ(line,*,IOSTAT=ios) seed_value
          IF (ios /= 0)
     &        CALL mass_error('Invalid comparison seed',path)
          EXIT
        ENDDO
      ENDIF
      CLOSE(unit_no)
      IF (parsed < 0 .OR. parsed > 5)
     &  CALL mass_error('m_select must be in 0..5',path)
      IF (seed_value <= 0)
     &  CALL mass_error('comparison seed must be positive',path)
      m_select = parsed
      mass_framework = framework_value
      comparison_seed = seed_value
      CALL require_physical_compatibility()
      IF (mass_framework == 1) CALL seed_comparison_random()
      END SUBROUTINE read_mass_selection

      SUBROUTINE seed_comparison_random()
      INTEGER :: seed_size
      INTEGER, ALLOCATABLE :: seed_vector(:)
      CALL RANDOM_SEED(SIZE=seed_size)
      ALLOCATE(seed_vector(seed_size))
      seed_vector = comparison_seed
      CALL RANDOM_SEED(PUT=seed_vector)
      DEALLOCATE(seed_vector)
      END SUBROUTINE seed_comparison_random

      SUBROUTINE load_selected_table(input_dir)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: input_dir
      CHARACTER(LEN=512) :: folder,path
      CALL require_physical_compatibility()
      folder = 'Input file/'
      IF (PRESENT(input_dir)) folder = input_dir
      IF (mass_framework == 1) THEN
        IF (m_select == 5) CALL require_physical_compatibility()
        IF (comparison_loaded .AND.
     &      comparison_loaded_dir == folder .AND.
     &      comparison_loaded_model == m_select) RETURN
        CALL load_comparison_table(folder)
        RETURN
      ENDIF
      IF (m_select == 0) RETURN
      IF (frdm_loaded .AND. frdm_loaded_dir == folder) RETURN
      path = TRIM(folder)//'frdm12.dat'
      CALL read_frdm12_production(path)
      frdm_loaded = .TRUE.
      frdm_loaded_dir = folder
      END SUBROUTINE load_selected_table

      SUBROUTINE load_comparison_table(folder)
      CHARACTER(LEN=*), INTENT(IN) :: folder
      CHARACTER(LEN=512) :: path
      CHARACTER(LEN=32) :: file_name
      INTEGER :: iz,ia
      comparison_present = .FALSE.
      comparison_native = .FALSE.
      comparison_fallback_reported = .FALSE.
      comparison_usage_reported = .FALSE.
      comparison_native_count = 0
      comparison_fallback_count = 0
      comparison_ignored_light_count = 0
      comparison_fallback_lookups = 0
      comparison_fallback_detail_count = 0
      comparison_shared_light_lookups = 0
      comparison_native_unique = 0
      comparison_fallback_unique = 0
      comparison_shared_light_unique = 0
      comparison_usage_unique = 0
      comparison_query_calls = 0
      path = TRIM(folder)//'massfrdm.dat'
      CALL read_comparison_file(path,0)
      SELECT CASE(m_select)
      CASE(0)
        file_name = 'massfrdm.dat'
      CASE(1)
        file_name = 'frdm12.dat'
      CASE(2)
        file_name = 'mass-hfb02.dat'
      CASE(3)
        file_name = 'KTUY05_m246.dat'
      CASE(4)
        file_name = 'wsm14.dat'
      CASE(5)
        file_name = 'ldmm.dat'
      CASE DEFAULT
        CALL require_physical_compatibility()
      END SELECT
      IF (m_select == 5) THEN
        WRITE(*,'(A)')
     &    'LDMM_COLUMN5_INFERRED_UNVERIFIED_PROVENANCE: local ' //
     &    'ldmm.dat column-5 candidate atomic mass-excess surface; ' //
     &    'not verified native model'
      ENDIF
      IF (m_select > 0) THEN
        path = TRIM(folder)//TRIM(file_name)
        CALL read_comparison_file(path,m_select)
      ELSE
        DO ia=1,max_a
          DO iz=0,max_z
            IF (comparison_present(iz,ia))
     &        comparison_native(iz,ia) = iz >= 8 .AND.
     &        ia-iz >= 8
          ENDDO
        ENDDO
      ENDIF
      DO ia=1,max_a
        DO iz=0,MIN(ia,max_z)
          IF (.NOT.comparison_present(iz,ia)) CYCLE
          IF (iz >= 8 .AND. ia-iz >= 8) THEN
            IF (comparison_native(iz,ia)) THEN
              comparison_native_count = comparison_native_count+1
            ELSE
              comparison_fallback_count =
     &          comparison_fallback_count+1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      comparison_loaded = .TRUE.
      comparison_loaded_dir = folder
      comparison_loaded_model = m_select
      WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')
     & 'MASS_COMPARE_FRAMEWORK mass_framework=1 m_select=',m_select,
     & ' map=',COUNT(comparison_present),' native=',
     & comparison_native_count,' fallback=',
     & comparison_fallback_count,' ignored_light=',
     & comparison_ignored_light_count,' domain_records=',
     & comparison_native_count+comparison_fallback_count
      WRITE(*,'(A,A)') 'MASS_COMPARE_MODEL=',TRIM(file_name)
      WRITE(*,'(A)') 'MASS_COMPARE_COMMON_DOMAIN Z>=8 N>=8; '//
     & 'outside=massfrdm shared_light'
      END SUBROUTINE load_comparison_table

      SUBROUTINE read_comparison_file(path,model)
      CHARACTER(LEN=*), INTENT(IN) :: path
      INTEGER, INTENT(IN) :: model
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=4) :: element
      INTEGER :: unit_no,ios,line_no,iz,inn,ia,flag
      LOGICAL, SAVE :: file_seen(0:max_z,0:max_a)
      REAL :: me,mexp,eps(4),beta(4),sp,micro,binding,sigma
      REAL :: micro_fl,me_fl,alpha(3),ws4
      LOGICAL :: has_value
      file_seen = .FALSE.
      OPEN(NEWUNIT=unit_no,FILE=TRIM(path),STATUS='OLD',
     & ACTION='READ',IOSTAT=ios)
      IF (ios /= 0) CALL mass_error('Cannot open comparison table',path)
      line_no = 0
      DO
        READ(unit_no,'(A)',IOSTAT=ios) line
        IF (ios < 0) EXIT
        line_no = line_no+1
        IF (ios /= 0) CALL mass_error('Comparison read failed',path,
     &      line_no)
        IF (LEN_TRIM(line) == 0) CYCLE
        IF (VERIFY(line(1:1),'#!') == 0) CYCLE
        IF (model == 3 .AND. INDEX(ADJUSTL(line),'ZZ') == 1) CYCLE
        me = 0.
        has_value = .TRUE.
        SELECT CASE(model)
        CASE(0)
          READ(line,*,IOSTAT=ios) iz,ia,element,flag,me
        CASE(1)
          IF (LEN_TRIM(line(126:135)) == 0)
     &      CALL mass_error('Blank FRDM12 Mth field',path,line_no)
          READ(line,'(3I5,4F10.2,4F10.3,5F10.2,F10.3,2F10.2)',
     &      IOSTAT=ios) iz,inn,ia,eps,beta,sp,micro,binding,me,
     &      mexp,sigma,micro_fl,me_fl
          IF (ios == 0 .AND. ia /= iz+inn)
     &      CALL mass_error('A differs from Z+N',path,line_no)
        CASE(2)
          READ(line,'(2I4,1X,A2,1X,I1,2F10.3)',IOSTAT=ios)
     &      iz,ia,element,flag,mexp,me
          has_value = LEN_TRIM(line(24:33)) > 0
          IF (.NOT.has_value .AND. LEN_TRIM(line(14:23)) == 0)
     &      CALL mass_error('Missing HFB mass fields',path,line_no)
        CASE(3)
          READ(line,*,IOSTAT=ios) iz,inn,me,micro,alpha
          IF (ios == 0) ia = iz+inn
        CASE(4)
          READ(line,*,IOSTAT=ios) ia,iz,ws4,me
        CASE(5)
          READ(line,*,IOSTAT=ios) iz,inn,ia,micro,me
          IF (ios == 0 .AND. ia /= iz+inn)
     &        CALL mass_error('A differs from Z+N',path,line_no)
          IF (ios == 0 .AND. .NOT.ieee_is_finite(micro))
     &        CALL mass_error('Non-finite LDMM column 4',path,line_no)
        CASE DEFAULT
          CALL require_physical_compatibility()
        END SELECT
        IF (ios /= 0) CALL mass_error('Malformed mass record',path,
     &      line_no)
        IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR. ia > max_a)
     &      CALL mass_error('Nuclide index outside 0:1000',path,line_no)
        IF (iz > ia) CALL mass_error('Z exceeds A',path,line_no)
        IF (model == 1 .OR. model == 2 .OR. model == 3 .OR.
     &      model == 4) THEN
          IF (.NOT.ieee_is_finite(me) .AND. has_value)
     &      CALL mass_error('Non-finite mass excess',path,line_no)
        ELSE IF (.NOT.ieee_is_finite(me)) THEN
          CALL mass_error('Non-finite mass excess',path,line_no)
        ENDIF
        IF (file_seen(iz,ia))
     &      CALL mass_error('Duplicate nuclide',path,line_no)
        file_seen(iz,ia) = .TRUE.
        IF (model > 0 .AND.
     &      (iz < 8 .OR. ia-iz < 8))
     &      comparison_ignored_light_count =
     &        comparison_ignored_light_count+1
        IF (.NOT.has_value) CYCLE
        IF (model == 0) THEN
          comparison_mass(iz,ia) = me
          comparison_present(iz,ia) = .TRUE.
        ELSE IF (iz >= 8 .AND. ia-iz >= 8) THEN
          comparison_mass(iz,ia) = me
          comparison_present(iz,ia) = .TRUE.
          comparison_native(iz,ia) = .TRUE.
        ENDIF
      ENDDO
      CLOSE(unit_no)
      END SUBROUTINE read_comparison_file

      SUBROUTINE read_frdm12_production(path)
      CHARACTER(LEN=*), INTENT(IN) :: path
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=4) :: element
      INTEGER :: unit_no,ios,line_no,iz,inn,ia,flag
      INTEGER :: record_count
      REAL :: eps(4),beta(4),sp,emic,ebind,mth,mexp,sigma
      REAL :: emic_fl,mth_fl
      frdm_present = .FALSE.
      frdm_seen = .FALSE.
      record_count = 0
      OPEN(NEWUNIT=unit_no,FILE=TRIM(path),STATUS='OLD',
     &  ACTION='READ',IOSTAT=ios)
      IF (ios /= 0) CALL mass_error('Cannot open FRDM12 table',path)
      line_no = 0
      DO
        READ(unit_no,'(A)',IOSTAT=ios) line
        IF (ios < 0) EXIT
        line_no = line_no+1
        IF (ios /= 0) CALL mass_error('FRDM12 read failed',path,
     &      line_no)
        IF (LEN_TRIM(line) == 0) CYCLE
        IF (VERIFY(line(1:1),'#!') == 0) CYCLE
        IF (LEN_TRIM(line(56:65)) == 0)
     &      CALL mass_error('Blank FRDM12 beta2 field',path,line_no)
        IF (LEN_TRIM(line(106:115)) == 0)
     &      CALL mass_error('Blank FRDM12 Emic field',path,line_no)
        IF (LEN_TRIM(line(126:135)) == 0)
     &      CALL mass_error('Blank FRDM12 Mth field',path,line_no)
        READ(line,'(3I5,4F10.2,4F10.3,5F10.2,F10.3,2F10.2)',
     &      IOSTAT=ios) iz,inn,ia,eps,beta,sp,emic,ebind,mth,
     &      mexp,sigma,emic_fl,mth_fl
        IF (ios /= 0)
     &      CALL mass_error('Malformed FRDM12 record',path,line_no)
        IF (ia /= iz+inn)
     &      CALL mass_error('FRDM12 A differs from Z+N',path,
     &      line_no)
        IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR. ia > max_a)
     &      CALL mass_error('FRDM12 nuclide outside 0:1000',path,
     &      line_no)
        IF (iz > ia)
     &      CALL mass_error('FRDM12 Z exceeds A',path,line_no)
        IF (frdm_seen(iz,ia))
     &      CALL mass_error('Duplicate FRDM12 nuclide',path,line_no)
        IF (.NOT.ieee_is_finite(mth) .OR.
     &      .NOT.ieee_is_finite(emic) .OR.
     &      .NOT.ieee_is_finite(beta(1)))
     &      CALL mass_error('Non-finite FRDM12 tuple field',path,
     &      line_no)
        frdm_seen(iz,ia) = .TRUE.
        frdm_present(iz,ia) = .TRUE.
        frdm_mth(iz,ia) = mth
        frdm_emic(iz,ia) = emic
        frdm_beta2(iz,ia) = beta(1)
        record_count = record_count+1
      ENDDO
      CLOSE(unit_no)
      IF (record_count == 0)
     &  CALL mass_error('FRDM12 table has no records',path)
      frdm12_native_count = record_count
      END SUBROUTINE read_frdm12_production

      SUBROUTINE lookup_selected_frdm12(iz,ia,mth,emic,beta2,found)
      INTEGER, INTENT(IN) :: iz,ia
      REAL, INTENT(OUT) :: mth,emic,beta2
      LOGICAL, INTENT(OUT) :: found
      mth = 0.
      emic = 0.
      beta2 = 0.
      found = .FALSE.
      CALL require_physical_compatibility()
      IF (m_select == 0) RETURN
      IF (.NOT.frdm_loaded) CALL load_selected_table()
      IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR. ia > max_a) RETURN
      found = frdm_present(iz,ia)
      IF (.NOT.found) RETURN
      mth = frdm_mth(iz,ia)
      emic = frdm_emic(iz,ia)
      beta2 = frdm_beta2(iz,ia)
      END SUBROUTINE lookup_selected_frdm12

      SUBROUTINE validate_selected_source(iz,ia)
      INTEGER, INTENT(IN) :: iz,ia
      REAL :: mth,emic,beta
      LOGICAL :: found
      IF (m_select == 0 .AND. mass_framework == 0) RETURN
      CALL require_physical_compatibility()
      IF (mass_framework == 1) THEN
        CALL lookup_selected_mass(iz,ia,mth,found)
        IF (found) RETURN
        WRITE(error_unit,'(A,I0,A,I0)')
     &    'MASS_COMPARE_MASS_UNCOVERED Z=',iz,' A=',ia
        STOP 1
      ENDIF
      CALL lookup_selected_frdm12(iz,ia,mth,emic,beta,found)
      IF (found) RETURN
      WRITE(error_unit,'(A,I0,A,I0,A,I0,A)')
     &  'MASS_SELECT_SOURCE_UNCOVERED m_select=',m_select,
     &  ' model=FRDM12 Z=',iz,' A=',ia,
     &  ' missing=FRDM12 native initial-source mass'
      STOP 1
      END SUBROUTINE validate_selected_source

      SUBROUTINE inspect_selected_table(model,input_dir)
      INTEGER, INTENT(IN) :: model
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: input_dir
      CHARACTER(LEN=512) :: folder,path
      CHARACTER(LEN=32) :: file_name
      CALL require_inspection_mode()
      folder = 'Input file/'
      IF (PRESENT(input_dir)) folder = input_dir
      WRITE(*,'(A)') 'MASS_TABLE_INSPECTION_ONLY: format parsing; '//
     &  'physical compatibility not established'
      IF (model < 1 .OR. model > 5)
     &  CALL mass_error('Inspection model must be in 1..5',folder)
      IF (inspect_loaded_select == model .AND.
     &    inspect_loaded_dir == folder) RETURN
      SELECT CASE(model)
      CASE(1)
        file_name = 'frdm12.dat'
      CASE(2)
        file_name = 'mass-hfb02.dat'
      CASE(3)
        file_name = 'KTUY05_m246.dat'
      CASE(4)
        file_name = 'wsm14.dat'
      CASE(5)
        file_name = 'ldmm.dat'
      END SELECT
      present_mass = .FALSE.
      native_mass = .FALSE.
      native_mass_count = 0
      fallback_mass_count = 0
C     Legacy values are inspection comparisons, not production fallback.
      path = TRIM(folder)//'massfrdm.dat'
      CALL read_native_file(path,0)
      path = TRIM(folder)//TRIM(file_name)
      CALL read_native_file(path,model)
      IF (native_mass_count == 0)
     &  CALL mass_error('Selected table has no mass records',path)
      fallback_mass_count = COUNT(present_mass .AND. .NOT.native_mass)
      inspect_loaded_select = model
      inspect_loaded_dir = folder
      WRITE(*,'(A,I0,2A)') 'inspection model=',model,': ',TRIM(path)
      WRITE(*,'(A,I0,A,I0)') 'Native masses: ',native_mass_count,
     &  '; legacy fallback masses: ',fallback_mass_count
      END SUBROUTINE inspect_selected_table

      SUBROUTINE read_native_file(path,model)
      CHARACTER(LEN=*), INTENT(IN) :: path
      INTEGER, INTENT(IN) :: model
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=4) :: element
      INTEGER :: unit_no,ios,line_no,iz,inn,ia,flag
      REAL :: me,mexp,eps(4),beta(4),sp,micro,binding,sigma
      REAL :: micro_fl,me_fl,alpha(3),ws4
      LOGICAL, SAVE :: seen(0:max_z,0:max_a)
      LOGICAL :: has_value
      seen = .FALSE.
      OPEN(NEWUNIT=unit_no,FILE=TRIM(path),STATUS='OLD',
     &  ACTION='READ',IOSTAT=ios)
      IF (ios /= 0) CALL mass_error('Cannot open mass table',path)
      line_no = 0
      DO
        READ(unit_no,'(A)',IOSTAT=ios) line
        IF (ios < 0) EXIT
        IF (ios /= 0) CALL mass_error('Read failed',path,line_no+1)
        line_no = line_no+1
        IF (LEN_TRIM(line) == 0) CYCLE
        IF (VERIFY(line(1:1),'#!') == 0) CYCLE
      IF (model == 3 .AND. INDEX(ADJUSTL(line),'ZZ') == 1) CYCLE
        has_value = .TRUE.
        SELECT CASE(model)
        CASE(0)
C         Read legacy values only for the offline comparison map.
          READ(line,*,IOSTAT=ios) iz,ia,element,flag,me
        CASE(1)
C         FRDM2012: preserve empty experimental fields and trailing notes.
          READ(line,'(3I5,4F10.2,4F10.3,5F10.2,F10.3,2F10.2)',
     &      IOSTAT=ios) iz,inn,ia,eps,beta,sp,micro,binding,me,
     &      mexp,sigma,micro_fl,me_fl
          IF (ios == 0 .AND. ia /= iz+inn)
     &      CALL mass_error('A differs from Z+N',path,line_no)
          IF (LEN_TRIM(line(126:135)) == 0)
     &      CALL mass_error('Missing FRDM theoretical mass',path,line_no)
        CASE(2)
C         HFB-2: Mexp may be blank; Mth occupies its own fixed field.
          READ(line,'(2I4,1X,A2,1X,I1,2F10.3)',IOSTAT=ios)
     &      iz,ia,element,flag,mexp,me
          has_value = LEN_TRIM(line(24:33)) > 0
          IF (.NOT.has_value .AND. LEN_TRIM(line(14:23)) == 0)
     &      CALL mass_error('Missing HFB mass fields',path,line_no)
        CASE(3)
          READ(line,*,IOSTAT=ios) iz,inn,me,micro,alpha
          IF (ios == 0) ia = iz+inn
        CASE(4)
          READ(line,*,IOSTAT=ios) ia,iz,ws4,me
        CASE(5)
C         LDMM numeric layout: Z,N,A,unknown field,unverified mass field.
C         Source, units, and physical column definitions remain unverified.
          READ(line,*,IOSTAT=ios) iz,inn,ia,micro,me
          IF (ios == 0) THEN
            IF (ia /= iz+inn)
     &        CALL mass_error('A differs from Z+N',path,line_no)
          ENDIF
        END SELECT
        IF (ios /= 0) CALL mass_error('Malformed mass record',path,line_no)
        IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR. ia > max_a)
     &    CALL mass_error('Nuclide index outside 0:1000',path,line_no)
        IF (iz > ia) CALL mass_error('Z exceeds A',path,line_no)
        IF (.NOT.ieee_is_finite(me))
     &    CALL mass_error('Non-finite mass excess',path,line_no)
        IF (seen(iz,ia)) CALL mass_error('Duplicate nuclide',path,line_no)
        seen(iz,ia) = .TRUE.
C       Absent HFB Mth is not zero: keep the legacy comparison value.
        IF (.NOT.has_value) CYCLE
        values(iz,ia) = me
        present_mass(iz,ia) = .TRUE.
        IF (model > 0) THEN
          native_mass(iz,ia) = .TRUE.
          native_mass_count = native_mass_count+1
        ENDIF
      ENDDO
      CLOSE(unit_no)
      END SUBROUTINE read_native_file

      SUBROUTINE lookup_selected_mass(iz,ia,me,found,is_native)
      INTEGER, INTENT(IN) :: iz,ia
      REAL, INTENT(OUT) :: me
      LOGICAL, INTENT(OUT) :: found
      LOGICAL, INTENT(OUT), OPTIONAL :: is_native
      REAL :: emic,beta
      me = 0.
      found = .FALSE.
      IF (PRESENT(is_native)) is_native = .FALSE.
      CALL require_physical_compatibility()
      IF (mass_framework == 1) THEN
        CALL load_selected_table()
        IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR.
     &      ia > max_a) RETURN
        IF (iz > ia) RETURN
        IF (.NOT.comparison_present(iz,ia)) THEN
          IF (iz == 0 .AND. ia >= 2 .AND. ia <= 4) THEN
            me = REAL(ia)*8.0713171
            found = .TRUE.
            CALL record_comparison_usage(iz,ia,.FALSE.,.TRUE.)
          ENDIF
          RETURN
        ENDIF
        me = comparison_mass(iz,ia)
        found = .TRUE.
        IF (PRESENT(is_native))
     &      is_native = comparison_native(iz,ia)
        CALL record_comparison_usage(iz,ia,
     &       comparison_native(iz,ia),
     &       iz < 8 .OR. ia-iz < 8)
        RETURN
      ENDIF
      IF (m_select == 0) RETURN
      CALL lookup_selected_frdm12(iz,ia,me,emic,beta,found)
      IF (PRESENT(is_native)) is_native = found
      END SUBROUTINE lookup_selected_mass

      SUBROUTINE report_mass_comparison_usage()
      IF (mass_framework /= 1) RETURN
      CALL write_comparison_usage()
      END SUBROUTINE report_mass_comparison_usage

      SUBROUTINE record_comparison_usage(iz,ia,is_native,is_shared)
      INTEGER, INTENT(IN) :: iz,ia
      LOGICAL, INTENT(IN) :: is_native,is_shared
      CHARACTER(LEN=24) :: source
      comparison_query_calls = comparison_query_calls+1
      IF (.NOT.is_native) THEN
        comparison_fallback_lookups =
     &      comparison_fallback_lookups+1
        IF (is_shared) THEN
          comparison_shared_light_lookups =
     &        comparison_shared_light_lookups+1
        ENDIF
      ENDIF
      IF (comparison_usage_reported(iz,ia)) RETURN
      comparison_usage_reported(iz,ia) = .TRUE.
      comparison_usage_unique = comparison_usage_unique+1
      IF (is_native) THEN
        comparison_native_unique = comparison_native_unique+1
      ELSE IF (is_shared) THEN
        comparison_shared_light_unique =
     &      comparison_shared_light_unique+1
      ELSE
        comparison_fallback_unique = comparison_fallback_unique+1
      ENDIF
      IF (.NOT.is_native .AND.
     &    .NOT.comparison_fallback_reported(iz,ia)) THEN
        comparison_fallback_reported(iz,ia) = .TRUE.
        IF (comparison_fallback_detail_count < 20) THEN
          comparison_fallback_detail_count =
     &        comparison_fallback_detail_count+1
          source = 'massfrdm'
          IF (is_shared) source = 'massfrdm(shared_light)'
          WRITE(*,'(A,I0,A,I0,A,I0,A,A)')
     &        'MASS_COMPARE_FALLBACK model=',m_select,
     &        ' Z=',iz,' A=',ia,' source=',TRIM(source)
        ENDIF
      ENDIF
      CALL write_comparison_usage()
      END SUBROUTINE record_comparison_usage

      SUBROUTINE write_comparison_usage()
      WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')
     & 'MASS_COMPARE_USAGE native_unique=',comparison_native_unique,
     & ' fallback_unique=',comparison_fallback_unique,
     & ' shared_light_unique=',comparison_shared_light_unique,
     & ' total_queries=',comparison_usage_unique,
     & ' total_calls=',comparison_query_calls,
     & ' fallback_calls=',comparison_fallback_lookups,
     & ' shared_light_calls=',comparison_shared_light_lookups
      END SUBROUTINE write_comparison_usage

      SUBROUTINE inspect_selected_mass(iz,ia,me,found,is_native)
      INTEGER, INTENT(IN) :: iz,ia
      REAL, INTENT(OUT) :: me
      LOGICAL, INTENT(OUT) :: found
      LOGICAL, INTENT(OUT), OPTIONAL :: is_native
      me = 0.
      found = .FALSE.
      IF (PRESENT(is_native)) is_native = .FALSE.
      CALL require_inspection_mode()
      IF (inspect_loaded_select < 1 .OR.
     &    inspect_loaded_select > 5)
     &  CALL mass_error('Inspection table has not been parsed',
     &  'inspection-only')
      IF (iz < 0 .OR. iz > max_z .OR. ia < 1 .OR. ia > max_a) RETURN
      found = present_mass(iz,ia)
      IF (found) me = values(iz,ia)
      IF (PRESENT(is_native)) is_native = native_mass(iz,ia)
      END SUBROUTINE inspect_selected_mass

      SUBROUTINE overlay_selected_masses(table,zmax,amax)
      INTEGER, INTENT(IN) :: zmax,amax
      REAL, INTENT(INOUT) :: table(0:zmax,0:amax)
      INTEGER :: iz,ia
      CALL require_physical_compatibility()
      IF (mass_framework == 1) THEN
        CALL load_selected_table()
        DO ia=1,MIN(amax,max_a)
          DO iz=0,MIN(zmax,max_z)
            IF (comparison_present(iz,ia))
     &        table(iz,ia) = comparison_mass(iz,ia)
          ENDDO
        ENDDO
        RETURN
      ENDIF
      IF (m_select == 0) RETURN
      IF (.NOT.frdm_loaded) CALL load_selected_table()
      DO ia=1,MIN(amax,max_a)
        DO iz=0,MIN(zmax,max_z)
          IF (frdm_present(iz,ia))
     &      table(iz,ia) = frdm_mth(iz,ia)
        ENDDO
      ENDDO
      END SUBROUTINE overlay_selected_masses
      END MODULE mass_selection
