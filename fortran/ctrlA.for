      PROGRAM CONTROL_ANA
      IMPLICIT NONE

      real, parameter :: dkin=1., mkin=1000.
      real, parameter :: delt=0.5, TMAX=150
      integer, parameter :: TEAMAX = 10
      integer, parameter :: CROSSMAX=30, ANGMAX=100, LINEAR_MODE = 0
      integer, parameter :: exctop = 1000      ! MeV
      
      integer NTH_RUN
      integer NTH_NOH
      integer NRMIN, NRMAX

      integer npar
      real X00CUT
      integer kktot, kkmin, kkmax
      real rr0, pp0
      integer READMESON

      integer COLM, GEMINI
      integer ANGM_DEC, ANGM_LIGHT, TESTANG
      integer SURF_COAL
      integer gemini_in_out
      
      integer ECHECK
      
      integer MCFR
      integer TESTEX, OUTENE

      real DECNEXC1, DECNEXC2
      integer DECIZ      
      integer PRINTCLKIN
      integer PRINTDEC, Z11, Z22

      integer TEAE
      integer EXCCROSS, TEAMUL
      real EX1DOWN, EX1UP
      real EX2DOWN, EX2UP
      real EX3DOWN
      real EXPNT
      integer PRTEX_WRTPHAS

      integer CNPCS, IPCS, TEELECS, TECOLD

      integer LNVICS, ELEKED, DDCSS_LNKED, DDCSS_LNKED_S
      real DIFF_T
      real B_DIFF

      integer FSVT

      real alpha,batat,gmmat
      real gsur,gtau11,eta11
      real gsym,symm
      integer momp


!      Fission and IMF parameters
      integer :: IMF_OPTION_ANA, Z_IMF_MIN_ANA
      real :: T_DELAY_ANA, SIG_DELAY_ANA, KRAMER_ANA, B_SCALE
      LOGICAL :: LESTONE_ANA, SHARP_DELAY_ANA

!      Evaporation plus gamma emission
      LOGICAL :: TL_IWBC_ANA, K_SUM_ANA, POLARIZATION_ANA
      integer :: EXOTIC_INDEX_ANA, MASS_OPTION_ANA
      real :: RATIO_ANA, e2_strength_ANA, e1_strength_ANA, THRESHOLD_ANA

!      Level density
      integer :: ADEN_TYPE_ANA
      real :: ADEN_0_ANA
      real :: K0, Kinf
      real :: A_SCALE_ANA

      integer I, IIP, NFILE, DEC1, DEC2, DEC3, NFOLD, NFOLD1, NRUN_PER

      character(len=20) NAME1, NAME2, NAME3, MKDIR1, MKDIR2, NOHUP1, 
     &  NOHUP2
      character(len=100) IFORT, NOHUP4, NOHUP3, DEL_CMD, MKDIR3, MONIT
      character(len=40) NAME
      character(len=3) FOLDER

      integer T_PER, J, NSTOP, EVENT, EVENT1
      integer, ALLOCATABLE, DIMENSION(:) :: STOPPED
      
      ! TEELECS
      integer TEELECSCOUNT(ceiling(TMAX/DELT),TEAMAX,2) ! (1热/2冷碎片,时间演化元素个数,时间数)
      integer TEELESPECIES, TEELECSR
      character*10 TEELECS_LABEL(TEAMAX)
      integer TEELECS_SETINF(TEAMAX,2)

      ! DDCSS_LNKED
      integer NEVENT, NNUCL, II, N_ANG, III , N_KIN_E, NSPECIES, NANGLE
      real KIN_E, LNKED(crossmax, ceiling(mkin/dkin)), LNKEDR, DDCSSR
      real DDCSS(crossmax, ANGMAX, ceiling(mkin/dkin))
      real ANG(crossmax, ANGMAX), READER2, MOMI
      character(len=20) DATANAME
      character(len=10) SUFFIX, READER1
      character(len=100) OPENNAME
      character(len=10) NUCL1(crossmax)
      integer BINRANGE(10, 2), BIN_N(10), N_BIN, NOW, BINWID,SUMMED(100)
      real BINLOC
      real E_BIN(10, 100)

      ! CNPCS
      real ELE_COUNT(crossmax)
      real CNPCS_CZ(150), CNPCS_CN(150), CNPCS_CA(300)
      real CNPCS_HZ(150), CNPCS_HN(150), CNPCS_HA(300)
      real READ_CNPCS
      character*10 CNPCS_LABEL(crossmax)
      integer CNPCS_SETINF(crossmax,2)
      integer NFILE_ANA, CSSPECIES

      ! EXCCROSS
      integer exccross1(0:ceiling(exctop/10.)), IE
      real exccross2(300), exce
      
      open (2004, FILE = 'Input file/INPUT_ANA.in', FORM = 'FORMATTED', 
     &  STATUS ='OLD')
      read (2004,*) NFILE
      read (2004,*) NRUN_PER
      read (2004,*) NAME
      read (2004,*) T_PER

      open (2003, FILE = 'Input file/anahyp0.in', FORM = 'FORMATTED',
     &  STATUS ='OLD')
      read (2003,*) ! CONT_ANA
      read (2003,*) NTH_RUN
      read (2003,*) NTH_NOH
      read (2003,*) NRMIN, NRMAX
      read (2003,*) ! Simulation
      read (2003,*) npar
      read (2003,*) X00CUT
      read (2003,*) kktot, kkmin, kkmax
      read (2003,*) rr0, pp0
      read (2003,*) READMESON
      read (2003,*) ! Models
      read (2003,*) COLM, GEMINI
      read (2003,*) ANGM_DEC, ANGM_LIGHT, TESTANG
      read (2003,*) SURF_COAL
      read (2003,*) ! Coalescence
      read (2003,*) MCFR
      read (2003,*) TESTEX, OUTENE
      read (2003,*) ! Gemini
      read (2003,*) DECNEXC1, DECNEXC2, DECIZ      
      read (2003,*) PRINTCLKIN
      read (2003,*)  PRINTDEC, Z11, Z22
      read (2003,*)  gemini_in_out
      read (2003,*) ! EnergyCheck
      read (2003,*) ECHECK
      read (2003,*) ! TEAE
      read (2003,*) TEAE
      read (2003,*) EXCCROSS, TEAMUL
      read (2003,*) EX1DOWN, EX1UP
      read (2003,*) EX2DOWN, EX2UP
      read (2003,*) EX3DOWN
      read (2003,*) EXPNT
      read (2003,*) PRTEX_WRTPHAS
      read (2003,*) ! Yields
      read (2003,*) CNPCS, IPCS, TEELECS, TECOLD
      read (2003,*) ! Kinematics
      read (2003,*) LNVICS, ELEKED, DDCSS_LNKED, DDCSS_LNKED_S, DIFF_T
      read (2003,*) B_DIFF
      read (2003,*) ! Tasks
      read (2003,*) FSVT
      read (2003,*) ! Dynamics
      read (2003,*) alpha,batat,gmmat
      read (2003,*) gsur,gtau11,eta11
      read (2003,*) gsym,symm
      read (2003,*) momp
      read (2003,*) ! Gemini settings
      read (2003,*) ! Fission and IMF parameters
      read (2003,*) IMF_OPTION_ANA
      read (2003,*) T_DELAY_ANA
      read (2003,*) SIG_DELAY_ANA
      read (2003,*) SHARP_DELAY_ANA
      read (2003,*) Z_IMF_MIN_ANA
      read (2003,*) KRAMER_ANA
      read (2003,*) LESTONE_ANA
      read (2003,*) B_SCALE
      read (2003,*) ! Evaporation plus gamma emission
      read (2003,*) TL_IWBC_ANA
      read (2003,*) K_SUM_ANA
      read (2003,*) POLARIZATION_ANA
      read (2003,*) EXOTIC_INDEX_ANA
      read (2003,*) MASS_OPTION_ANA
      read (2003,*) e2_strength_ANA
      read (2003,*) e1_strength_ANA
      read (2003,*) THRESHOLD_ANA
      read (2003,*) ! Level density
      read (2003,*) ADEN_TYPE_ANA
      read (2003,*) ADEN_0_ANA    
      read (2003,*) K0
      read (2003,*) Kinf
      read (2003,*) A_SCALE_ANA

      read (2004,*) x00cut
      read (2004,*) alpha,batat,gmmat
      read (2004,*) gsur,gtau11,eta11
      read (2004,*) gsym,symm
      read (2004,*) momp

      open (2005, FILE = 'MONITOR_ANA.dat', STATUS = 'UNKNOWN') ! 监视器
      call RANDOM_NUMBER(MOMI) 
      write (2005,*) MOMI

!     将exfrag.for编译为名为Name的可执行的二进制文件(.exe)
      NAME1 = 'ifort -o '
      NAME3 = ' exfrag.for'
      IFORT = TRIM(NAME1)//TRIM(NAME)//TRIM(NAME3)
      call SYSTEM(IFORT) 

      do I = 1, NFILE
        call Folder001(I, FOLDER)
        IFORT = 'rm -f ./'//TRIM(NAME)//'_ana'//FOLDER
        call SYSTEM(IFORT) ! 删除原来的LQMD程序
      end do

      call SLEEP(3)

      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI

! 按事件数分为NFILE个文件夹
      do I = 1, NFILE
        call Folder001(I, FOLDER) 
        NAME2 = TRIM(NAME)//'_ana'//TRIM(FOLDER)
        IFORT = TRIM('cp '//TRIM(NAME)//' '//TRIM(NAME2)) 
        call SYSTEM(IFORT) ! copy .exe 
        MKDIR1 = 'mkdir '
        MKDIR2 = 'ana'//TRIM(FOLDER)
        MKDIR3 = MKDIR1//' '//TRIM(MKDIR2)
        ! nohup ./exfrag >ana001/nohup0001.out 2>&1 &
        NOHUP1 = 'nohup ./'
        NOHUP2 = ' &'
        NOHUP4 = '>'//'ana'//TRIM(FOLDER)//'/'//'nohup0001.out 2>&1'
        NOHUP3 = TRIM(NOHUP1)//TRIM(NAME2)
     &    //' '//TRIM(NOHUP4)//TRIM(NOHUP2)
        call SYSTEM(MKDIR3) ! 创建目录（文件夹）

        NRMIN = (I - 1)*NRUN_PER + 1
        NRMAX = I*NRUN_PER
        NTH_RUN = I
        NTH_NOH = 1
        NPAR = NFILE*NRUN_PER
        open (2, FILE = 'Input file/anahyp.in', FORM = 'FORMATTED', 
     &    STATUS='UNKNOWN')
        write (2,*) '%%CONT_ANA'
        write (2,*) NTH_RUN
        write (2,*) NTH_NOH
        write (2,*) NRMIN, NRMAX
        write (2,*) '%%Simulation'
        write (2,*) npar
        write (2,*) X00CUT
        write (2,*) kktot, kkmin, kkmax
        write (2,*) rr0, pp0
        write (2,*) READMESON
        write (2,*) '%%Models'
        write (2,*) COLM, GEMINI
        write (2,*) ANGM_DEC, ANGM_LIGHT, TESTANG
        write (2,*) SURF_COAL
        write (2,*) '%%Coalescence'
        write (2,*) MCFR
        write (2,*) TESTEX, OUTENE
        write (2,*) '%%Gemini'
        write (2,*) DECNEXC1, DECNEXC2, DECIZ      
        write (2,*) PRINTCLKIN
        write (2,*) PRINTDEC, Z11, Z22
        write (2,*) gemini_in_out
        write (2,*) 'EnergyCheck'
        write (2,*) ECHECK
        write (2,*) '%%TEAE'
        write (2,*) TEAE
        write (2,*) EXCCROSS, TEAMUL
        write (2,*) EX1DOWN, EX1UP
        write (2,*) EX2DOWN, EX2UP
        write (2,*) EX3DOWN
        write (2,*) EXPNT
        write (2,*) PRTEX_WRTPHAS
        write (2,*) '%%Yields'
        write (2,*) CNPCS, IPCS, TEELECS, TECOLD
        write (2,*) '%%Kinematics'
        write (2,*) LNVICS, ELEKED, DDCSS_LNKED, DDCSS_LNKED_S, DIFF_T
        write (2,*) B_DIFF
        write (2,*) '%%Tasks'
        write (2,*) FSVT
        write (2,*) '%%Dynamics'
        write (2,*) alpha,batat,gmmat
        write (2,*) gsur,gtau11,eta11
        write (2,*) gsym,symm
        write (2,*) momp
        write (2,*) '%%Gemini settings'
        write (2,*)  '***** fission and IMF parameters *****'
        write (2,*) IMF_OPTION_ANA
        write (2,*) T_DELAY_ANA
        write (2,*) SIG_DELAY_ANA
        write (2,*) SHARP_DELAY_ANA
        write (2,*) Z_IMF_MIN_ANA
        write (2,*) KRAMER_ANA
        write (2,*) LESTONE_ANA
        write (2,*) B_SCALE
        write (2,*)  '***** evaporation + gamma emission *****'
        write (2,*) TL_IWBC_ANA
        write (2,*) K_SUM_ANA
        write (2,*) POLARIZATION_ANA
        write (2,*) EXOTIC_INDEX_ANA
        write (2,*) MASS_OPTION_ANA
        write (2,*) e2_strength_ANA
        write (2,*) e1_strength_ANA
        write (2,*) THRESHOLD_ANA
        write (2,*)  '***** level density(能级密度) *****'
        write (2,*) ADEN_TYPE_ANA
        write (2,*) ADEN_0_ANA    
        write (2,*) K0
        write (2,*) Kinf
        write (2,*) A_SCALE_ANA
        close(2)
        call SYSTEM(NOHUP3)
        call SLEEP(1)       !防止输入文件修改过快
      end do
      
      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI
      call SLEEP(T_PER)
      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI

      NSTOP = 0
      allocate(STOPPED(NFILE)) ! 动态调整数组大小
      STOPPED(1:NFILE) = 0
      J = 0
      do
        EVENT = 0
        EVENT1 = 0
        if (J == NFILE) then
          J = 0
          NSTOP = 0
          EVENT = 0
          EVENT1 = 0
          STOPPED(1:NFILE) = 0
          print*
          print*
        end if
        J = J + 1
        call Folder001(J, FOLDER)
        MONIT = 'ana'//FOLDER//'/'//'MONITOR'//FOLDER//'_0001.dat'
        open (3, FILE = MONIT, FORM = 'FORMATTED', STATUS = 'UNKNOWN')
        do
          read (3, *, END = 33, ERR = 33) EVENT
        end do
33      close(3)

        call RANDOM_NUMBER(MOMI)
        write (2005,*) MOMI

        print*, 'NFILE =', ' ', J, ' ', 'EVENT =', ' ', EVENT
        call SLEEP(T_PER)

        call RANDOM_NUMBER(MOMI)
        write (2005,*) MOMI

        open (3, FILE = MONIT, FORM = 'FORMATTED', STATUS = 'OLD')
        do
          read (3, *, END = 34, ERR = 34) EVENT1
        end do
34      close(3)

! EVENT = 0 相当于当前计算全部舍弃
        if (EVENT == EVENT1.OR.EVENT == 0) then
          if (EVENT < NRUN_PER) STOPPED(J) = 1
          NSTOP = NSTOP + 1
        end if

        if (NSTOP == NFILE) then
          EXIT
        end if
      end do

********************************************************************
*** 读取输入文件与初始化
            
      ! TEELECS
      if (TEELECS == 1) then
        open (19, FILE = 'Input file/TEELECS_in.in', STATUS='OLD')
        TEELESPECIES = 0 ! 元素种类数
        read (19,*)
        read (19,*)
        do II = 1, TEAMAX
        read (19,*,end=228,ERR=228) TEELECS_LABEL(II),
     &    TEELECS_SETINF(II,1), TEELECS_SETINF(II,2)
          TEELESPECIES = TEELESPECIES + 1
        end do
228     continue
        do II = 1, TEELESPECIES
          do III = 1, CEILING(TMAX/DELT)
            TEELECSCOUNT(III,II,1:2) = 0
          end do
        end do     
        close(19)    
      end if
      
      ! DDCSS_LNKED
      if (DDCSS_LNKED == 1) then
        open (556, FILE = 'Input file/DDCSS_LNKED_in.in',STATUS = 'OLD')
        NSPECIES=0
        read (556,*)
        read (556,*)
        do I=1,crossmax
          read (556, *, END=2366, err=2366) READER1, READER2
          NSPECIES = NSPECIES+  1
        end do
2366    NANGLE=0
        read (556,*)
        do I = 1, ANGMAX
          read (556, *, END=2367,err=2367) READER2
          NANGLE = NANGLE + 1
        end do
2367    do I = 1, 10
          BINRANGE(I, 1:2) = 0
          BIN_N(I) = 0
        end do
        I = 0
        read (556,*)
        read (556,*)
        do
          I = I + 1
          read (556, *, END=2368, ERR=2368) BINRANGE(I, 1), 
     &      BINRANGE(I, 2), BIN_N(I)
        end do
2368    close(556)
        I = I - 1
        N_BIN = I
        do I = 1, N_BIN
          E_BIN(I, 1:100) = 0
        end do
        do I = 1, NSPECIES
          do II = 1, NANGLE
            do III = 1, ceiling(mkin/dkin)
              DDCSS(I, II, III) = 0.
            end do
          end do
        end do
        do I = 1, NSPECIES
          do II = 1, ceiling(mkin/dkin)
            LNKED(I, II) = 0.
          end do
        end do
      end if

      ! CNPCS
      if (CNPCS == 1) then
        open (19, FILE = 'Input file/CNPCS_in.in', STATUS='OLD')
        CSSPECIES = 0 ! 元素种类数
        read (19,*)
        read (19,*)
        do II = 1, TEAMAX
          read (19,*,end=289,ERR=289) CNPCS_LABEL(II),
     &      CNPCS_SETINF(II,1), CNPCS_SETINF(II,2)
            CSSPECIES = CSSPECIES + 1
        end do        
289     continue
        close (19)
        ELE_COUNT(1:crossmax) = 0.
        CNPCS_CZ(1:150) = 0.
        CNPCS_CN(1:150) = 0.
        CNPCS_CA(1:300) = 0.
        CNPCS_HZ(1:150) = 0.
        CNPCS_HN(1:150) = 0.
        CNPCS_HA(1:300) = 0.
      end if
      
      ! EXCCROSS
      EXCCROSS1(0:ceiling(EXCTOP/10.))=0
      EXCCROSS2(1:300) = 0.
********************************************************************
      
      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI
      
********************************************************************
*** 读取初步分析结果
      NEVENT = 0
      NFILE_ANA = 0
      do J = 1, NFILE
        EVENT = 0
        call Folder001(J, FOLDER)
        MONIT = 'ana'//FOLDER//'/'//'MONITOR'//FOLDER//'_0001.dat'
        open (3, FILE = MONIT, FORM = 'FORMATTED', 
     &    STATUS = 'UNKNOWN')
        do
          read (3, *, END = 6, ERR = 6) EVENT
        end do
6       close(3)
        if (STOPPED(J) == 1.OR.EVENT == 0) then
          cycle
        end if
        ! if (MOD(J, 10) == 0) cycle 意义不明
        NEVENT = NEVENT + NRUN_PER
        NFILE_ANA = NFILE_ANA + 1
        call Folder001(J, FOLDER)
                    
        if (TEELECS == 1) then
          ! read (*,*) ! 'DEBUG: === ENTERING TEELECS OUTPUT ==='
          do I = 1, 2
            DATANAME = 'TEELECS_out'
            SUFFIX = '.dat'
            OPENNAME = 'ana'//TRIM(FOLDER)//'/'//TRIM(DATANAME)//
     &        TRIM(FOLDER)//'_0001'//TRIM(SUFFIX)
            open (20, FILE = OPENNAME, STATUS = 'OLD')
            read (20,*)
            read (20,*) ! '-------------------------'
            if (I == 1) then
              read (20,*) ! 'RESIDUES' ! 热碎片
            else
              read (20,*) ! 'FRAGMENTS'! 冷碎片
            end if
            read (20,*) ! '-------------------------'
            read (20,*)
            do II = 1, TEELESPECIES
              read (20,*) TEELECS_LABEL(II)
              do III = 1, ceiling(TMAX/DELT)
                read (20,*) KIN_E, TEELECSR
                TEELECSCOUNT(III,II,I) = TEELECSR  
     &                                 + TEELECSCOUNT(III,II,I)
              end do
              read (20,*)
            end do
          end do
          close(20)
        end if
        ! DDCSS_LNKED
        if (DDCSS_LNKED == 1) then
          !读入
          DATANAME = 'DDCSS_LNKED_out'
          SUFFIX       = '.dat'
          OPENNAME = 'ana'//TRIM(FOLDER)//'/'//TRIM(DATANAME)//
     &      TRIM(FOLDER)//'_0001'//TRIM(SUFFIX)
          open (7, FILE = OPENNAME, FORM = 'FORMATTED',
     &      STATUS = 'OLD')
          read (7,*)
          read (7,*)
          read (7,*)
          do I = 1, NSPECIES
            read (7,*) NUCL1(I)
            do II = 1, NANGLE
              read (7,*) ANG(I, II)
              do III = 1, ceiling(mkin/dkin)
                read (7,*) KIN_E, DDCSSR
                DDCSS(I, II, III) = DDCSS(I, II, III)
     &                            + DDCSSR
              end do
            end do
          end do
          read (7,*)
          do I = 1, NSPECIES
            read (7,*)
            do II = 1, ceiling(mkin/dkin)
              read (7,*) KIN_E, LNKEDR
              LNKED(I, II) = LNKED(I, II) + LNKEDR
            end do
          end do
          close(7)
        end if

        ! CNPCS
        if (CNPCS == 1) then
          ! COLD FRAGMENTS
          DATANAME = 'CNPCS_COLD'
          SUFFIX = '.dat'
          OPENNAME = 'ana'//TRIM(FOLDER)//'/'//TRIM(DATANAME)//
     &      TRIM(FOLDER)//'_0001'//TRIM(SUFFIX) ! ana001/CNPCS_COLD_0001.dat
          open (7, FILE = OPENNAME, FORM = 'FORMATTED',
     &      STATUS = 'OLD')
          read (7,*)
          do I = 1, 150
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_CZ(I) = CNPCS_CZ(I) + READ_CNPCS
          end do
          read (7,*)
          read (7,*)
          do I = 1, 150
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_CN(I) = CNPCS_CN(I) + READ_CNPCS
          end do
          read (7,*)
          read (7,*)
          do I = 1, 300
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_CA(I) = CNPCS_CA(I) + READ_CNPCS
          end do
          read (7,*) 
          read (7,*)
          do I = 1, crossmax
            read (7,*) IIP, READ_CNPCS
            ELE_COUNT(I) = ELE_COUNT(I) + READ_CNPCS
          end do  
          close(7)

          ! HOT_FRAGMENTS
          DATANAME = 'CNPCS_HOT'
          SUFFIX = '.dat'
          OPENNAME = 'ana'//TRIM(FOLDER)//'/'//TRIM(DATANAME)//
     &      TRIM(FOLDER)//'_0001'//TRIM(SUFFIX)
          open (7, FILE = OPENNAME, FORM = 'FORMATTED',
     &      STATUS = 'OLD')
          read (7,*)
          do I = 1, 150
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_HZ(I) = CNPCS_HZ(I) + READ_CNPCS
          end do
          read (7,*)
          read (7,*)
          do I = 1, 150
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_HN(I) = CNPCS_HN(I) + READ_CNPCS
          end do
          read (7,*)
          read (7,*)
          do I = 1, 300
            read (7,*) KIN_E, READ_CNPCS
            CNPCS_HA(I) = CNPCS_HA(I) + READ_CNPCS
          end do
          close(7)
        end if
        
        ! EXCCROSS
        if (EXCCROSS) then
          DATANAME = 'EXCCROSS'
          SUFFIX       = '.dat'
          OPENNAME = 'ana'//TRIM(FOLDER)//'/'//TRIM(DATANAME)//
     &      TRIM(FOLDER)//'_0001'//TRIM(SUFFIX)
          open (7, FILE = OPENNAME, FORM = 'FORMATTED',
     &      STATUS = 'OLD')
          read (7,*)
          read (7,*)
          read (7,*)
          
          do I = 0, ceiling(EXCTOP/10.)
            read (7,*) IE, exce
            EXCCROSS2(I) = EXCCROSS2(I) + exce
          end do
          read (7,*)
          read (7,*)
          read (7,*)
          do I = 1, 300
            read (7,*) IE, exce
            EXCCROSS1(I) = exccross1(I) + exce
          end do
        end if
      end do
********************************************************************

      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI
      
********************************************************************
*** 输出再分析结果
      call SYSTEM('mkdir ana_results')
      
      ! TEELECS
      if (TEELECS == 1) then
        if (GEMINI == 0) then
          open (21, FILE = 'ana_results'//'/'//'TEELECS_out1.dat',
     &      STATUS = 'UNKNOWN')
        elseif (GEMINI == 1) then
          open (21, FILE = 'ana_results'//'/'//'TEELECS_out.dat',
     &      STATUS = 'UNKNOWN')
        end if
        ! write (*,*) 'DEBUG: === ENTERING TEELECS OUTPUT ==='
        do I = 1, 2
          write (21,*)
          write (21,*) '-------------------------'
          if (I == 1) then
            write (21,*) 'RESIDUES' ! 热碎片
          else
            write (21,*) 'FRAGMENTS'! 冷碎片
          end if
          write (21,*) '-------------------------'
          write (21,*)
          do II = 1, TEELESPECIES
            write (21,*) TEELECS_LABEL(II)
            do III = 1, CEILING(TMAX/DELT)
              if (TEELECSCOUNT(III,II,I) == 0) cycle
              write (21,*) III*DELT, 
     &          TEELECSCOUNT(III,II,I)/FLOAT(NEVENT)
            end do
            write (21,*)
          end do
        end do
        close(21)
      end if      
      
      ! DDCSS_LNKED
      if (DDCSS_LNKED == 1) then
        if (GEMINI == 0) then
          open (8, FILE = 'ana_results'//'/'//'DDCSS_LNKED_out1.dat',
     &      STATUS = 'UNKNOWN')
        elseif (GEMINI == 1) then
          open (8, FILE = 'ana_results'//'/'//'DDCSS_LNKED_out.dat',
     &      STATUS = 'UNKNOWN')
        end if
          
        write (8,*) ' Double Differential Cross Sections'
        do I = 1, NSPECIES
          write (8,*) NUCL1(I)
          do II = 1, NANGLE
            write (8,*) 'Theta= ', ' ' , ANG(I, II)
            if (LINEAR_MODE == 1) then
              do III = 1, ceiling(mkin/dkin)
                if (DDCSS(I, II, III) == 0.) cycle
                write (8,*) III*dkin, DDCSS(I, II, III)/NEVENT
              end do
            else
              NOW = 1
              SUMMED(1:100) = 0
              do J = 1, 10
                    E_BIN(J, 1:100) = 0.
              end do
              do III = 1, ceiling(mkin/dkin)
                if (III*DKIN.GE.(BINRANGE(NOW, 2)+1)) then
                  do J = 1, BIN_N(NOW)
                    if (E_BIN(NOW, J) == 0.) cycle
                    write (8,*) BINRANGE(NOW,1) + (J-1)*BINWID
     &                + BINWID/2.-1, E_BIN(NOW, J)/NEVENT/
     &                SUMMED(J)
                  end do
                  NOW = NOW + 1
                  SUMMED(1:100) = 0
                end if
                if (NOW.GT.N_BIN) then
                  EXIT
                end if
                BINWID = (BINRANGE(NOW,2) - 
     &            BINRANGE(NOW, 1) + 1)/BIN_N(NOW)
                BINLOC = (III - BINRANGE(NOW, 1) + 1)/FLOAT(BINWID)
                E_BIN(NOW,CEILING(BINLOC))= E_BIN(NOW, 
     &            CEILING(BINLOC)) + DDCSS(I,II,III)
                SUMMED(CEILING(BINLOC)) = SUMMED(CEILING(BINLOC))
     &            + 1
              end do
            end if
          end do
        end do
        write (8,*) ' Light Nuclide Kinetic Energy Distribution'
        do I = 1, NSPECIES
          write (8,*) NUCL1(I)
          do II = 1, ceiling(mkin/dkin)
            if (LNKED(I, II) == 0) cycle
            write (8,*) II*dkin, LNKED(I, II)/NEVENT
          end do
        end do
        close(8)
      end if

      ! CNPCS
      if (CNPCS == 1) then
        ! COLD FRAGMENTS
        if (GEMINI == 0) then
          open (8, FILE = 'ana_results'//'/'//'CNPCS_COLD1.dat', 
     &          STATUS = 'UNKNOWN')
        elseif (GEMINI == 1) then
          open (8, FILE = 'ana_results'//'/'//'CNPCS_COLD.dat', 
     &          STATUS = 'UNKNOWN')
        end if

        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF Z:'
        do I = 1, 150
          if (CNPCS_CZ(I) == 0) cycle
          write (8,*) I, CNPCS_CZ(I)/NEVENT
        end do
        write (8,*)

        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF N:'
        do I = 1, 150
          if (CNPCS_CN(I) == 0) cycle
          write (8,*) I, CNPCS_CN(I)/NEVENT
        end do
        write (8,*)

        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF A:'
        do I = 1, 300
          if (CNPCS_CA(I) == 0) cycle
          write (8,*) I, CNPCS_CA(I)/NEVENT
        end do
        write (8,*)

        write (8,*) 'CROSS SECTION OF SPECIAL ELEMETNS:'
        do I = 1, CSSPECIES
          if (ELE_COUNT(I) == 0) cycle
          write (8,*) CNPCS_LABEL(I), ELE_COUNT(I)/NEVENT
        end do
        write (8,*)    
        close(8)

        ! HOT FRAGMENTS
        if (GEMINI == 0) then
          open (8, FILE = 'ana_results'//'/'//'CNPCS_HOT1.dat', 
     &      STATUS = 'UNKNOWN')
        elseif (GEMINI == 1) then
          open (8, FILE = 'ana_results'//'/'//'CNPCS_HOT.dat', 
     &      STATUS = 'UNKNOWN')
        end if      
        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF Z:'
        do I = 1, 150
          if (CNPCS_HZ(I) == 0) cycle
          write (8,*) I, CNPCS_HZ(I)/NEVENT
        end do
        write (8,*)
        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF N:'
        do I = 1, 150
          if (CNPCS_HN(I) == 0) cycle
          write (8,*) I, CNPCS_HN(I)/NEVENT
        end do
        write (8,*)
        write (8,*), 'TOTAL CROSS SECTION AS A FUNCTION OF A:'
        do I = 1, 300
          if (CNPCS_HA(I) == 0) cycle
          write (8,*) I, CNPCS_HA(I)/NEVENT
        end do
        write (8,*)   
        close(8)
      end if

      ! EXCCROSS
      if (EXCCROSS) then
        open (8, FILE = 'ana_results'//'/'//'EXCCROSS.dat', 
     &      STATUS = 'UNKNOWN')
          write (8,*) 'EXCITATION ENERGY CROSS SECTION'
          write (8,*)
          write (8,*) 'MeV      mb/MeV'
          do I = 0, ceiling(EXCTOP/10.)
            if (EXCCROSS2(I) == 0) cycle
            write (8,*) I*10, EXCCROSS2(I)/NEVENT
          end do
          write (8,*) 'EXCITATION ENERGY AS FUNCTION OF MASS 
     &      NUMBER'
          write (8,*)
          write (8,*) 'A      mb'
          do I = 1, 300
            if (EXCCROSS1(I) == 0) cycle
            write (8,*) I, EXCCROSS1(I)/NEVENT
          end do
      end if
          
********************************************************************        
          
      call RANDOM_NUMBER(MOMI)
      write (2005,*) MOMI

      print*
      print*, 'Total number of analyzed events:'
      print*, NEVENT
            
      CONTAINS
                   
********************************************************************
*** get the folder ID 001~999     
      subroutine Folder001(KK, FOLDER)
      
      IMPLICIT NONE
      integer, intent(in) :: KK
      character(len=3), intent(out) :: FOLDER
      integer :: NFOLD, NFOLD1, DEC1, DEC2, DEC3
    
      NFOLD = KK
      NFOLD1 = NFOLD
      DEC1 = FLOOR(NFOLD1/100.0) + 48
      NFOLD1 = NFOLD1 - FLOOR(NFOLD1/100.0) * 100
      DEC2 = FLOOR(NFOLD1/10.0) + 48
      NFOLD1 = NFOLD1 - FLOOR(NFOLD1/10.0) * 10
      DEC3 = FLOOR(NFOLD1/1.0) + 48
      FOLDER = CHAR(DEC1)//CHAR(DEC2)//CHAR(DEC3)
      return
      
      end subroutine Folder001
********************************************************************
      
      END PROGRAM