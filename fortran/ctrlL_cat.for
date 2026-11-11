	PROGRAM	CONTROL
	IMPLICIT NONE
	
	INTEGER, PARAMETER :: N_CATCH = 10


	INTEGER	nrmin, nrmax, nth_run, nth_noh
	INTEGER	npr, nzpr, nta, nzta, nanti
	REAL	x00min, x00max, x00stp, x00cut,nx00
	REAL	epp, nt, dt, ITMIN, ITMAX
	REAL	alpha, beta, gamma
	REAL	gsur, gtau, eta, gsym, symm, rho0
	REAL	asym, bsym, csym, gamsy
	REAL	vpa, vpb, vpc, q01, p01, e01
	REAL	r0, ds, c1, c2
	REAL	rsur, rmin
	REAL	bta1, bta2, bta3
	INTEGER	ish, ipot, iwrite
	REAL	cdelt, cmeps, cmom, xsym, momp
	INTEGER	iflag1, iflag2, iflag3, iflag4, iflag5, iflag6,iflag7
	INTEGER	mdel, mdin, npop, nkaon
	REAL	nsym, xmu
	REAL	kfile, v01, initi, ibin, ird, crsc, ifile, gradi ,nmint
	INTEGER	IPRIN,ifdelt,iflagcro,iflagpau,ipoli
	INTEGER	OUTBARYON, OUTMESON
	REAL	GETX00CUT, GETX00CUT_IN, GETX00CUT_ST, G_ROT
	
	INTEGER	SURF_COAL1
	INTEGER	PRINT_CLUST, PRINT_POT
	REAL	NDATA
	INTEGER	K_START_SURF, K_END_SURF
	REAL	RR0, D0_SURF, COAL_PATH, CLUST_LEAST_DIST_COEF
	REAL	H0(5), MONI, MONI1
	REAL	H1_SURF, INCL4_6



	INTEGER	NFILE, NRUN_PER
	INTEGER NFOLD, NFOLD1
	INTEGER DEC1, DEC2, DEC3, DEC4
	INTEGER I, FINISH, J
	INTEGER NTH_NOH_BIN(999), LAST_EVENT(999)
	INTEGER NOHUP, NOHUP11
	INTEGER EVENT, EVENT1
	INTEGER LENG0(N_CATCH), LENG1(N_CATCH), LENG2(N_CATCH), 
     &    LENG3(N_CATCH), LENG4(N_CATCH), LENG5(N_CATCH), LENG6(N_CATCH)
	INTEGER	LEAD(N_CATCH), ENDING(N_CATCH), T_PER_ANA
	INTEGER	T_1ST, T_PER, CALL_ANA, FAIL, NFILE_ANA, RUN_PER_ANA
      INTEGER IERR

	CHARACTER(LEN = 20) :: NAME
	CHARACTER(LEN = 9) :: NAME1, NAME3
	CHARACTER(LEN = 38) :: IFORT
	CHARACTER(3) :: CP1
	CHARACTER(4) :: CP2
	CHARACTER(27) :: CP3
	CHARACTER(6) :: MKDIR1
	CHARACTER(3) :: MKDIR2
	CHARACTER(9) :: MKDIR3
	CHARACTER(3) :: FOLDER
	CHARACTER(6) :: TOUCH1
	CHARACTER(8) :: TOUCH2
	CHARACTER(17) :: TOUCH3
	CHARACTER(8) :: NOHUP1
	CHARACTER(2) :: NOHUP2
	CHARACTER(54) :: NOHUP3
	CHARACTER(24) :: NOHUP4
	CHARACTER(4) :: NNOHUP
	CHARACTER(LEN = 1000000) :: CATCH1, CATCH2, CATCH3, CATCH4 
	CHARACTER(LEN = 100) :: NAME_ANA


	INTEGER	KK,KKK
	REAL	X1, Y1, Z1, PX1, PY1, PZ1, EM1
	INTEGER	ISO1, ISY1, ID1, IDR1

      OPEN(2004, FILE = 'ERR.out', STATUS = 'UNKNOWN')
	OPEN(2003, FILE = 'Input file/LQMD0.in', STATUS = 'OLD')
	OPEN(2, FILE = 'Input file/INPUT.in', STATUS = 'OLD')

	READ(2, *)	NFILE
	READ(2, *)	NRUN_PER
	READ(2, *)	NAME
	READ(2, *)	T_1ST
	READ(2, *)	T_PER
	READ(2, *)	CALL_ANA
	READ(2,* )	FAIL

!	NAME1 = 'ifort -o '
!	NAME3 = ' LQMD.for'
!c     ifort -o 可执行文件名(NAME) LQMD.for      
!	IFORT = NAME1//NAME//NAME3
!	CALL SYSTEM(IFORT)
!	CALL SLEEP(5)

	read(2003,*) nrmin, nrmax, nth_run, nth_noh
	read(2003,*) npr, nzpr, nta, nzta, nanti
	read(2003,*) x00min, x00max, x00stp, x00cut,nx00
	read(2003,*) epp, nt, dt, ITMIN, ITMAX
	read(2003,*) alpha, beta, gamma
	read(2003,*) gsur, gtau, eta, gsym, symm, rho0
	read(2003,*) asym, bsym, csym, gamsy
	read(2003,*) vpa, vpb, vpc, q01, p01, e01
	read(2003,*) r0, ds, c1, c2
	read(2003,*) rsur, rmin
	read(2003,*) bta1, bta2, bta3
	read(2003,*) ish, ipot, iwrite
	read(2003,*) cdelt, cmeps, cmom, xsym, momp
	read(2003,*) iflag1, iflag2, iflag3, iflag4, iflag5, iflag6,iflag7
	read(2003,*) mdel, mdin, npop, nkaon
	read(2003,*) nsym, xmu
	read(2003,*) kfile, v01, initi, ibin, ird, crsc, ifile, gradi 
     &			,nmint
	read(2003,*) IPRIN,ifdelt,iflagcro,iflagpau,ipoli
	read(2003,*) OUTBARYON, OUTMESON
	read(2003,*) GETX00CUT, GETX00CUT_IN, GETX00CUT_ST, G_ROT
	READ(2003,*)
	READ(2003,*) SURF_COAL1, PRINT_CLUST, PRINT_POT, NDATA
	READ(2003,*) K_START_SURF, K_END_SURF
	READ(2003,*) RR0, D0_SURF, COAL_PATH, CLUST_LEAST_DIST_COEF
	READ(2003,*) H0(1), H0(2), H0(3), H0(4), H0(5)
	READ(2003,*) H1_SURF, INCL4_6

	CLOSE(2003)

10	format(1X,6F10.5)
11	format(1X,8I4)
12	format(1X, I3, 1X, F5.1, 1X, 5I3, 1X, F5.1, 1X, I3)
13	format(1X,4F10.5, 1X, I3)
	
	NTH_NOH_BIN(1:999) = 0


	DO	I = 1, NFILE

		NFOLD = I
		NFOLD1 = NFOLD
		DEC1 = FLOOR(NFOLD1/100.) + 48
		NFOLD1 = NFOLD1 - FLOOR(NFOLD1/100.)*100.
		DEC2 = FLOOR(NFOLD1/10.) + 48
		NFOLD1 = NFOLD1 - FLOOR(NFOLD1/10.)*10.
		DEC3 = FLOOR(NFOLD1/1.) +48
		FOLDER = CHAR(DEC1)//CHAR(DEC2)//CHAR(DEC3)

		NRMIN = (I - 1)*NRUN_PER + 1
		LAST_EVENT(I) = NRMIN
		NRMAX = I*NRUN_PER
		NTH_RUN = I
		NTH_NOH = 1
		NTH_NOH_BIN(I) = 1 

	ENDDO
      
! 对文件夹001-100进行循环，检查是否所有事件都演化成功
      I = 0
	FINISH = 0
	DO
		IF(I.EQ.NFILE)	THEN
			FINISH = 0
          DO I = 1, NFILE
				PRINT*, NTH_NOH_BIN(I)
          ENDDO
          EXIT
		ENDIF

		I = I + 1

          call Folder001(I, FOLDER)
		call Folder0001(NTH_NOH_BIN(I), NNOHUP)

		DO
			OPEN(2, FILE = FOLDER//'/'//'initangm'//FOLDER//'_'//NNOHUP
     &	        //'.dat', STATUS = 'OLD',ACTION='READ', IOSTAT=IERR)		
			IF(IERR.EQ.0)THEN
				CLOSE(2)
				! 记录提交次数，重新提交则不计
				NTH_NOH_BIN(I) = NTH_NOH_BIN(I) + 1 
            call Folder0001(NTH_NOH_BIN(I),NNOHUP)
              ELSE
                  NTH_NOH_BIN(I) = NTH_NOH_BIN(I) - 1 
				CLOSE(2)
				EXIT
			ENDIF

          ENDDO
      ENDDO
      
********************************************************************
*** 合并与删除文件
      !CALL SYSTEM('mkdir results')
	! LENG0: 110 处字符串长度
	! LENG1: 111 处字符串长度
	! LENG2: 112 处字符串长度
      ! LENG3: 113 处字符串长度
      ! LENG4: 114 处字符串长度
	! LENG5: 115 处字符串长度
	! LENG6: 116 处字符串长度

	LENG0(1) = 3    ! 'cat'
	LENG1(1) = 26   ! ' 001/initangm001_0001.dat '
	LENG2(1) = 23   ! '>> results/initangm.dat'
	LENG3(1) = 19   ! '> 001/initangm.dat '
	LENG4(1) = 18   ! ' 001/initangm.dat '
	LENG5(1) = 8    ! 'rm -f'
	LENG6(1) = 28   ! ' ./001/initangm001_0001.dat '

	LENG0(2) = 3    ! 'cat'
	LENG1(2) = 27   ! ' 001/SURF_COAL001_0001.dat '
	LENG2(2) = 24   ! '> results/SURF_COAL.dat'
	LENG3(2) = 20   ! '> 001/SURF_COAL.dat '
	LENG4(2) = 19   ! ' 001/SURF_COAL.dat '
	LENG5(2) = 8    ! 'rm -f'
	LENG6(2) = 29   ! ' ./001/SURF_COAL001_0001.dat '

	LENG0(3) = 3    ! 'cat'
	LENG1(3) = 22   ! ' 001/coll001_0001.dat '
	LENG2(3) = 19   ! '> results/coll.dat'
	LENG3(3) = 23   ! '> 001/coll.dat '  
	LENG4(3) = 14   ! ' 001/coll.dat '    
	LENG5(3) = 8    ! 'rm -f'
	LENG6(3) = 24   ! ' ./001/coll001_0001.dat '

	LENG0(4) = 3    ! 'cat'
	LENG1(4) = 26   ! ' 001/COUNTSUR001_0001.out '
	LENG2(4) = 23   ! '> results/COUNTSUR.out'
	LENG3(4) = 19   ! '> 001/COUNTSUR.out '
	LENG4(4) = 18   ! ' 001/COUNTSUR.out '
	LENG5(4) = 8    ! 'rm -f'
	LENG6(4) = 28   ! ' ./001/COUNTSUR001_0001.out '

	LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1

!	KKK = 0
!	DO
!		OPEN(5511, FILE = '001/coll001_0001.dat', Status = 'OLD')
!		read(5511,23456)  kK,x1,y1,z1,px1,py1
!     &	,pz1,em1,iso1,isy1,id1,idr1
!		if(em1.eq.0.9380)THEN
!			KKK = 0
!			PRINT*, 'Please add a filehead and save!'
!		ELSE
!			KKK = 1
!		ENDIF
!		CLOSE(5511)
!		IF(KKK.EQ.1)EXIT
!		CALL SLEEP(1)
!	ENDDO
!23456	format(1x,i4,7f9.4,4i4)

	DO I = 1, NFILE
		call Folder001(I, FOLDER)
          call Folder0001(NTH_NOH_BIN(I), NNOHUP)
          
          ! cat 003/initangm.dat 003/initangm003_0002.dat ... > 003/initangm.dat
          WRITE(2004,*) 'cat in file '//FOLDER//': '//NNOHUP
          CATCH1(1:1000000) = ' '
	    CATCH2(1:1000000) = ' '
	    CATCH3(1:1000000) = ' '
	    CATCH4(1:1000000) = ' '
110	    CATCH1(1:LENG0(1)) = 'cat'
	    CATCH2(1:LENG0(2)) = 'cat'
	    CATCH3(1:LENG0(3)) = 'cat'
	    CATCH4(1:LENG0(4)) = 'cat'
          LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1 ! 计算当前末位
		DO J = 1, NTH_NOH_BIN(I)              
		    call Folder0001(J, NNOHUP)              
111			CATCH1(LEAD(1)+(J-1)*LENG1(1) : LEAD(1)-1+J*LENG1(1)) = 
     &		' '//FOLDER//'/'//'initangm'//FOLDER//'_'//NNOHUP//'.dat '
			CATCH2(LEAD(2)+(J-1)*LENG1(2) : LEAD(2)-1+J*LENG1(2)) = 
     &		' '//FOLDER//'/'//'SURF_COAL'//FOLDER//'_'//NNOHUP//
     &		'.dat '
			CATCH3(LEAD(3)+(J-1)*LENG1(3) : LEAD(3)-1+J*LENG1(3)) = 
     &		' '//FOLDER//'/'//'coll'//FOLDER//'_'//NNOHUP//'.dat '
			CATCH4(LEAD(4)+(J-1)*LENG1(4) : LEAD(4)-1+J*LENG1(4)) = 
     &		' '//FOLDER//'/'//'COUNTSUR'//FOLDER//'_'//NNOHUP//'.out '    
			ENDING(1:N_CATCH) = LEAD(1:N_CATCH)-1+J*LENG1(1:N_CATCH)		
		ENDDO		
		LEAD(1:N_CATCH) = ENDING(1:N_CATCH) + 1 ! 计算当前末位          
113       CATCH1(LEAD(1) : LEAD(1)+LENG3(1)-1)= '> '//FOLDER//
     &        '/initangm.dat '
          CATCH2(LEAD(2) : LEAD(2)+LENG3(2)-1)='> '//FOLDER//
     &        '/SURF_COAL.dat '
          CATCH3(LEAD(3) : LEAD(3)+LENG3(3)-1)= '> '//FOLDER//
     &        '/coll.dat '
          CATCH4(LEAD(4) : LEAD(4)+LENG3(4)-1)= '> '//FOLDER//
     &        '/COUNTSUR.out '
          ENDING(1:N_CATCH) = LEAD(1:N_CATCH)-1+LENG3(1:N_CATCH)
          LEAD(1:N_CATCH) = ENDING(1:N_CATCH) + 1 ! 计算当前末位 
          
          ! rm -f ./003/initangm003_0001.dat ./003/initangm003_0002.dat ...
115       CATCH1(LEAD(1) : LEAD(1)+LENG5(1)-1)= '&& rm -f'
          CATCH2(LEAD(2) : LEAD(2)+LENG5(2)-1)= '&& rm -f'
          CATCH3(LEAD(3) : LEAD(3)+LENG5(3)-1)= '&& rm -f'
          CATCH4(LEAD(4) : LEAD(4)+LENG5(4)-1)= '&& rm -f'
          ENDING(1:N_CATCH) = LEAD(1:N_CATCH)-1+LENG5(1:N_CATCH)
          LEAD(1:N_CATCH) = ENDING(1:N_CATCH) + 1 ! 计算当前末位
          DO J = 1, NTH_NOH_BIN(I)          
		    call Folder0001(J, NNOHUP)
116			CATCH1(LEAD(1)+(J-1)*LENG6(1) : LEAD(1)-1+J*LENG6(1)) = 
     &		' ./'//FOLDER//'/'//'initangm'//FOLDER//'_'//NNOHUP//'.dat '
			CATCH2(LEAD(2)+(J-1)*LENG6(2) : LEAD(2)-1+J*LENG6(2)) = 
     &		' ./'//FOLDER//'/'//'SURF_COAL'//FOLDER//'_'//NNOHUP//
     &		'.dat '
			CATCH3(LEAD(3)+(J-1)*LENG6(3) : LEAD(3)-1+J*LENG6(3)) = 
     &		' ./'//FOLDER//'/'//'coll'//FOLDER//'_'//NNOHUP//'.dat '
			CATCH4(LEAD(4)+(J-1)*LENG6(4) : LEAD(4)-1+J*LENG6(4)) = 
     &		' ./'//FOLDER//'/'//'COUNTSUR'//FOLDER//'_'//NNOHUP//'.out '
			ENDING(1:N_CATCH) = LEAD(1:N_CATCH)-1+J*LENG6(1:N_CATCH)		
		ENDDO
          WRITE(2004,*) CATCH1(1:ENDING(1))
          !CALL SYSTEM('rm -f ./'//FOLDER//'/initangm.dat')
          !CALL SYSTEM('rm -f ./'//FOLDER//'/SURF_COAL.dat')
          !CALL SYSTEM('rm -f ./'//FOLDER//'/coll.out')
          !CALL SYSTEM('rm -f ./'//FOLDER//'/COUNTSUR.out')
		CALL SYSTEM(CATCH1)
	    !CALL SYSTEM(CATCH2)
	    !CALL SYSTEM(CATCH3)
	    !CALL SYSTEM(CATCH4)
      ENDDO

      CALL SLEEP(3)

      ! cat 001/initangm.dat 002/initangm.dat ... > results/initangm.dat
      WRITE(2004,*) 'cat nfile: ', NFILE
      CATCH1(1:1000000) = ' '
	CATCH2(1:1000000) = ' '
	CATCH3(1:1000000) = ' '
	CATCH4(1:1000000) = ' '
      CATCH1(1:LENG0(1)) = 'cat'
	CATCH2(1:LENG0(2)) = 'cat'
	CATCH3(1:LENG0(3)) = 'cat'     
	CATCH4(1:LENG0(4)) = 'cat'
      LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1
      CALL SYSTEM('rm -f ./results/initangm.dat')
      CALL SYSTEM('rm -f ./results/SURF_COAL.dat')
	CALL SYSTEM('rm -f ./results/coll.dat')
	CALL SYSTEM('rm -f ./results/COUNTSUR.out')
      DO I = 1, NFILE
          call Folder001(I, FOLDER)
          LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1
114       CATCH1(LEAD(1) : LEAD(1)-1+LENG4(1)) = 
     &		' '//FOLDER//'/'//'initangm.dat '         
          CATCH2(LEAD(2) : LEAD(2)-1+LENG4(2)) = 
     &		' '//FOLDER//'/'//'SURF_COAL.dat '
          CATCH3(LEAD(3) : LEAD(3)-1+LENG4(3)) = 
     &		' '//FOLDER//'/'//'coll.dat '
          CATCH4(LEAD(4) : LEAD(4)-1+LENG4(4)) = 
     &		' '//FOLDER//'/'//'COUNTSUR.out '      
          ENDING(1:N_CATCH) = LEAD(1:N_CATCH) + LENG4(1:N_CATCH) - 1
          LEAD(1:N_CATCH) = ENDING(1:N_CATCH) + 1
112       CATCH1(LEAD(1) : LEAD(1)+LENG2(1)-1) = 
     &        '>> results/initangm'//'.dat'
	    CATCH2(LEAD(2) : LEAD(2)+LENG2(2)-1) = 
     &        '>> results/SURF_COAL'//'.dat'
	    CATCH3(LEAD(3) : LEAD(3)+LENG2(3)-1) = 
     &        '>> results/coll'//'.dat'
	    CATCH4(LEAD(4) : LEAD(4)+LENG2(4)-1) = 
     &        '>> results/COUNTSUR'//'.out'
          ENDING(1:N_CATCH) = LEAD(1:N_CATCH) + LENG2(1:N_CATCH) - 1
          WRITE(2004,*) CATCH1(1:ENDING(1))
	    CALL SYSTEM(CATCH1)
	    CALL SYSTEM(CATCH2)
	    CALL SYSTEM(CATCH3)
	    CALL SYSTEM(CATCH4)
          CALL SYSTEM('rm -f ./'//FOLDER//'/initangm.dat')
          CALL SYSTEM('rm -f ./'//FOLDER//'/SURF_COAL.dat')
          CALL SYSTEM('rm -f ./'//FOLDER//'/coll.dat')
          CALL SYSTEM('rm -f ./'//FOLDER//'/COUNTSUR.out')
      ENDDO
      close(2004)

********************************************************************
      
	IF(CALL_ANA.EQ.1)	THEN
		OPEN(10, FILE = 'INPUT_ANA.in', STATUS = 'UNKNOWN')	
		NFILE_ANA = NFILE
		RUN_PER_ANA = NRUN_PER - FAIL
		NAME_ANA = 'ctrlA'
		T_PER_ANA = T_PER
		WRITE(10,*)	NFILE_ANA
		WRITE(10,*)	RUN_PER_ANA
		WRITE(10,*)	TRIM(NAME)
		WRITE(10,*)	T_PER_ANA
		WRITE(10,*)	x00cut
		WRITE(10,*) alpha, beta, gamma
		WRITE(10,*) gsur, gtau, eta
		WRITE(10,*)	gsym, symm
		WRITE(10,*)	NINT(momp)
		CLOSE(10)

		CALL SYSTEM('ifort -o '//TRIM(NAME_ANA)//' ctrlA.for')

		STOP ! 终止程序，即只编译不运行

		CALL SYSTEM('nohup./ '//TRIM(NAME_ANA)//' >nohup_ana 2>&1')
		DO
			OPEN(10, FILE = 'MONITOR_ANA.dat', STATUS = 'UNKNOWN')
			DO
				READ(10,*,END=339,ERR=339)	MONI
			ENDDO
339			CLOSE(10)

			CALL SLEEP(T_PER_ANA)

			OPEN(10, FILE = 'MONITOR_ANA.dat', STATUS = 'UNKNOWN')
			DO
				READ(10,*,END=340,ERR=340)	MONI1
			ENDDO
340			CLOSE(10)

			IF(MONI1.EQ.MONI)STOP
		ENDDO
      ENDIF

********************************************************************
********************************************************************
********************************************************************
      
      CONTAINS
      
********************************************************************
*** get the folder ID 001~100      
      subroutine Folder001(KK,FOLDER)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KK
      CHARACTER(LEN=3), INTENT(OUT) :: FOLDER
      INTEGER :: NFOLD, NFOLD1, DEC1, DEC2, DEC3
    
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

********************************************************************
*** get the folder ID 0001~1000      
      subroutine Folder0001(KK,NNOHUP)
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KK
      CHARACTER(LEN=4), INTENT(OUT) :: NNOHUP
      INTEGER :: NOHUP, DEC1, DEC2, DEC3, DEC4
    
	NOHUP = KK
	DEC1 = FLOOR(NOHUP/1000.) + 48
	NOHUP = NOHUP - FLOOR(NOHUP/1000.)*1000.
	DEC2 = FLOOR(NOHUP/100.) +48
	NOHUP = NOHUP - FLOOR(NOHUP/100.)*100.
	DEC3 = FLOOR(NOHUP/10.) +48
	NOHUP = NOHUP - FLOOR(NOHUP/10.)*10.
	DEC4 = FLOOR(NOHUP/1.) +48
	NNOHUP = CHAR(DEC1)//CHAR(DEC2)//CHAR(DEC3)//CHAR(DEC4)
      return
      
      end subroutine Folder0001
********************************************************************

	END