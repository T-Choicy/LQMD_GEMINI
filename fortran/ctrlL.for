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
	INTEGER NTH_NOH_BIN(999), LAST_EVENT(999), LAST_EVENT1(999)
	INTEGER NOHUP, NOHUP11
	INTEGER EVENT, EVENT1, ITIME, ITIME1
      INTEGER LAST_TIME(999),LAST_TIME1(999)
	INTEGER LENG0(N_CATCH), LENG1(N_CATCH), LENG2(N_CATCH), 
     & LENG3(N_CATCH), LENG4(N_CATCH), LENG5(N_CATCH), LENG6(N_CATCH)
	INTEGER	LEAD(N_CATCH), ENDING(N_CATCH), T_PER_ANA
	INTEGER	T_1ST, T_PER, CALL_ANA, FAIL, NFILE_ANA, RUN_PER_ANA

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

	integer :: cur_kk,last_kk,ios
      real :: cur_x1,cur_y1,cur_z1
      real :: last_x1,last_y1,last_z1
      logical :: has_data
      
	open(2, file = 'Input file/INPUT.in', status = 'OLD')
	read(2, *)	NFILE
	read(2, *)	NRUN_PER
	read(2, *)	NAME
	read(2, *)	T_1ST
	read(2, *)	T_PER
	read(2, *)	CALL_ANA
	read(2,* )	FAIL
      
	open(2003, file = 'Input file/LQMD0.in', status = 'OLD')
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
	read(2003,*) kfile, v01, initi, ibin, ird, crsc, ifile, gradi, 
     & nmint
	read(2003,*) IPRIN,ifdelt,iflagcro,iflagpau,ipoli
	read(2003,*) OUTBARYON, OUTMESON
	read(2003,*) GETX00CUT, GETX00CUT_IN, GETX00CUT_ST, G_ROT
	read(2003,*)
	read(2003,*) SURF_COAL1, PRINT_CLUST, PRINT_POT, NDATA
	read(2003,*) K_START_SURF, K_END_SURF
	read(2003,*) RR0, D0_SURF, COAL_PATH, CLUST_LEAST_DIST_COEF
	read(2003,*) H0(1), H0(2), H0(3), H0(4), H0(5)
	read(2003,*) H1_SURF, INCL4_6
	close(2003)
10	format(1X,6F10.5)
11	format(1X,8I4)
12	format(1X, I3, 1X, F5.1, 1X, 5I3, 1X, F5.1, 1X, I3)
13    format(1X,4F10.5, 1X, I3)
	
      NAME1 = 'ifort -o '
	NAME3 = ' LQMD.for'
	IFORT = NAME1//NAME//NAME3
	call SYSTEM(IFORT) ! ifort -o 可执行文件名(NAME) LQMD.for  
   
	call SLEEP(3)

********************************************************************
*** 并行多个程序    
	NTH_NOH_BIN(1:999) = 0
	do I = 1, NFILE
		call Folder001(I, FOLDER)
		MKDIR1 = 'mkdir '
		MKDIR2 = FOLDER
		MKDIR3 = MKDIR1//MKDIR2
          call SYSTEM(MKDIR3) ! mkdir创建文件夹
          
		IFORT = TRIM( 'cp '//TRIM(NAME)//' '//TRIM(NAME)//
     &	'_'//TRIM(FOLDER) ) 
          call SYSTEM(IFORT) ! cp复制可执行程序

		nrmin = (I - 1)*NRUN_PER + 1
		LAST_EVENT(I) = nrmin
          LAST_EVENT1(I) = nrmin
		nrmax = I*NRUN_PER
		NTH_RUN = I
		NTH_NOH = 1
		NTH_NOH_BIN(I) = 1
		NOHUP1 = 'nohup ./'
		NOHUP2 = ' &'
		NOHUP4 = '>'//FOLDER//'/'//'nohup0001.out 2>&1'
		NOHUP3 = NOHUP1//TRIM(NAME)//'_'//
     &    TRIM(FOLDER)//' '//NOHUP4//NOHUP2
          
		open(1, file = 'Input file/LQMD.in', status = 'UNKNOWN')
		write(1,*) nrmin, nrmax, nth_run, nth_noh
		write(1,*) npr, nzpr, nta, nzta, nanti
		write(1,*) x00min, x00max, x00stp, x00cut,NINT(nx00)
		write(1,*) epp, NINT(nt), dt, NINT(ITMIN), NINT(ITMAX)
		write(1,*) alpha, beta, gamma
		write(1,10) gsur, gtau, eta, gsym, symm, rho0
		write(1,10) asym, bsym, csym, gamsy
		write(1,10) vpa, vpb, vpc, q01, p01, e01
		write(1,10) r0, ds, c1, c2
		write(1,*) rsur, rmin
		write(1,10) bta1, bta2, bta3
		write(1,11) ish, ipot, iwrite
		write(1,13) cdelt, cmeps, cmom, xsym, NINT(momp)
		write(1,11) iflag1, iflag2, iflag3, iflag4, iflag5, iflag6
     &				,iflag7
		write(1,11) mdel, mdin, npop, nkaon
		write(1,*) NINT(nsym), xmu
		write(1,12) NINT(kfile), v01, NINT(initi), NINT(ibin),  
     &     		NINT(ird),NINT(crsc), NINT(ifile), gradi ,NINT(nmint)
		write(1,11) IPRIN,ifdelt,iflagcro,iflagpau,ipoli
		write(1,11) OUTBARYON, OUTMESON
		write(1,*) NINT(GETX00CUT), GETX00CUT_IN, GETX00CUT_ST, 
     &				NINT(G_ROT)
		write(1,*)
		write(1,*) SURF_COAL1, PRINT_CLUST, PRINT_POT, NDATA
		write(1,*) K_START_SURF, K_END_SURF
		write(1,*) RR0, D0_SURF, COAL_PATH, CLUST_LEAST_DIST_COEF
		write(1,*) H0(1), H0(2), H0(3), H0(4), H0(5)
		write(1,*) H1_SURF, INCL4_6
		close(1)
          
		call SYSTEM(NOHUP3) ! nohup执行进程
		call SLEEP(1)	! 防止修改Input file/LQMD.in过快
      enddo
********************************************************************
      
	call SLEEP(T_1ST)	!从提交到输出第一个事件的记录所需时间

********************************************************************
*** 检查是否所有事件都演化成功	
      FINISH = 0 ! success if FINISH = NFILE
	do
		if(FINISH.EQ.NFILE)exit
          print*, 'NEW TEST' ! try again if they are not equal
          FINISH = 0
          do I = 1, NFILE
              if(FINISH.EQ.NFILE)exit
		    EVENT  = 0
		    EVENT1 = 0
		    call Folder001(I, FOLDER)
		    call Folder0001(NTH_NOH_BIN(I), NNOHUP)		
		    open(2, file = FOLDER//'/'//'monitor'//FOLDER//'_'//NNOHUP//
     &	        '.dat', status = 'UNKNOWN')		
		    ! (event,itime)=(nru,k), 范围为[nrmin,nrmax]
		    ! 若所有事件正常演化，则(event,itime) = (nrmax,nt)
		    ! 若(nru,k<nt)异常,nru事件需要重跑
              ! 若(nru,k=nt)异常,nru事件无需重跑
              do
			    read(2, *, end = 20) EVENT, ITIME
		    enddo
20		    continue
		    close(2)
              if(EVENT.NE.0)then
                  LAST_EVENT(I) = EVENT
                  LAST_TIME(I) = ITIME
              endif
              
              call SLEEP(T_PER)	! 每次check之间的时间间隔

		    open(2, file = FOLDER//'/'//'monitor'//FOLDER//'_'//NNOHUP//                   
     &	        '.dat', status = 'UNKNOWN')
		    open(5, file = FOLDER//'/'//'coll'//FOLDER//'_'//NNOHUP//
     &	        '.dat', status = 'UNKNOWN')
		    close(5)
		    open(7, file = FOLDER//'/'//'SURF_COAL'//FOLDER//'_'//NNOHUP//             
     &            '.dat', status = 'UNKNOWN')
		    close(7)
		    do
			    read(2, *, end = 21) EVENT1, ITIME1
		    enddo
21		    continue
		    close(2)
              if(EVENT1.NE.0)then
                  LAST_EVENT1(I) = EVENT1
                  LAST_TIME1(I) = ITIME1
              endif
              print*, 'FILE=', I, ' ',
     &            'EVENT=',LAST_EVENT1(I), ' ',
     &            'ITIME=', LAST_TIME1(I) 
              ! 对第I个程序判断是否演化完毕:
              ! EVENT等于最后一个事件
		    !print*, 'NFILE=', ' ', I, ' ', 'EVENT=', ' ', EVENT
		    if(EVENT.EQ.I*NRUN_PER.OR.(EVENT.EQ.I*NRUN_PER - 1))then
     			    FINISH = FINISH + 1
			    cycle
              endif
              ! 未跑完的则重新在未完成事件处重新提交
              ! （1）不同时间的两次check相同（即演化太久了）
              ! （2）当前事件不是最后一个事件时，重新提交程序
		    if(EVENT1.EQ.EVENT.AND.
     &            ITIME1.EQ.ITIME.AND.
     &            (EVENT.NE.I*NRUN_PER.OR.
     &                ITIME.NE.NINT(nt)))then
  
                  call SYSTEM(TRIM('pkill -9 -f "' // TRIM(NAME) 
     &                // '_' // TRIM(FOLDER) // '"')) ! 强制关闭
            
			    nrmin = LAST_EVENT(I) ! 从异常事件接续
                  if(ITIME.EQ.NINT(nt)) nrmin = LAST_EVENT(I) + 1
			    nrmax = I*NRUN_PER
			    NTH_RUN = I
				NTH_NOH_BIN(I) = NTH_NOH_BIN(I) + 1 ! 记录提交次数，重新提交则不计
				if(EVENT.EQ.0) NTH_NOH_BIN(I) = NTH_NOH_BIN(I) - 1
				NTH_NOH = NTH_NOH_BIN(I)
				call Folder0001(NTH_NOH_BIN(I), NNOHUP)
				NOHUP1 = 'nohup ./'
				NOHUP2 = ' &'
				NOHUP4 = '>'//FOLDER//'/'//'nohup'//NNOHUP//'.out 2>&1'
				NOHUP3 = NOHUP1//TRIM(NAME)//'_'//
     &                TRIM(FOLDER)//' '//NOHUP4//NOHUP2

				open(1, file = 'Input file/LQMD.in', status = 'UNKNOWN')
				write(1,*) nrmin, nrmax, nth_run, nth_noh
				write(1,*) npr, nzpr, nta, nzta, nanti
				write(1,*) x00min, x00max, x00stp, x00cut,NINT(nx00)
				write(1,*) epp, NINT(nt), dt, NINT(ITMIN), NINT(ITMAX)
				write(1,*) alpha, beta, gamma
				write(1,10) gsur, gtau, eta, gsym, symm, rho0
				write(1,10) asym, bsym, csym, gamsy
				write(1,10) vpa, vpb, vpc, q01, p01, e01
				write(1,10) r0, ds, c1, c2
				write(1,*) rsur, rmin
				write(1,10) bta1, bta2, bta3
				write(1,11) ish, ipot, iwrite
				write(1,13) cdelt, cmeps, cmom, xsym, NINT(momp)
				write(1,11) iflag1, iflag2, iflag3, iflag4, iflag5, iflag6
     &				,iflag7
				write(1,11) mdel, mdin, npop, nkaon
				write(1,*) NINT(nsym), xmu
				write(1,12) NINT(kfile), v01, NINT(initi), NINT(ibin),  
     &     		    NINT(ird),NINT(crsc), NINT(ifile), gradi ,NINT(nmint)
				write(1,11) IPRIN,ifdelt,iflagcro,iflagpau,ipoli
				write(1,11) OUTBARYON, OUTMESON
				write(1,*) NINT(GETX00CUT), GETX00CUT_IN, GETX00CUT_ST, 
     &				NINT(G_ROT)
				write(1,*)
				write(1,*) SURF_COAL1, PRINT_CLUST, PRINT_POT, NDATA
				write(1,*) K_START_SURF, K_END_SURF
				write(1,*) RR0, D0_SURF, COAL_PATH, CLUST_LEAST_DIST_COEF
				write(1,*) H0(1), H0(2), H0(3), H0(4), H0(5)
				write(1,*) H1_SURF, INCL4_6
				close(1)
                  
			    call SYSTEM(NOHUP3) ! 再次提交程序
              endif ! 重新提交结束
          enddo       
      enddo
********************************************************************
      !open(2004, file = 'ERR.out', status = 'UNKNOWN')
      !if(NFILE*NRUN_PER*nt/float(iwrite) .GT. 500000)then
      !    close(2004)
      !    STOP ! 数据量太大，终止合并
      !endif
      
********************************************************************
*** 合并与删除文件
      call SYSTEM('mkdir results')
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
	LENG2(2) = 24   ! '>> results/SURF_COAL.dat'
	LENG3(2) = 20   ! '> 001/SURF_COAL.dat '
	LENG4(2) = 19   ! ' 001/SURF_COAL.dat '
	LENG5(2) = 8    ! 'rm -f'
	LENG6(2) = 29   ! ' ./001/SURF_COAL001_0001.dat '

	LENG0(3) = 3    ! 'cat'
	LENG1(3) = 22   ! ' 001/coll001_0001.dat '
	LENG2(3) = 19   ! '>> results/coll.dat'
	LENG3(3) = 23   ! '> 001/coll.dat '  
	LENG4(3) = 14   ! ' 001/coll.dat '    
	LENG5(3) = 8    ! 'rm -f'
	LENG6(3) = 24   ! ' ./001/coll001_0001.dat '

	LENG0(4) = 3    ! 'cat'
	LENG1(4) = 26   ! ' 001/COUNTSUR001_0001.out '
	LENG2(4) = 23   ! '>> results/COUNTSUR.out'
	LENG3(4) = 19   ! '> 001/COUNTSUR.out '
	LENG4(4) = 18   ! ' 001/COUNTSUR.out '
	LENG5(4) = 8    ! 'rm -f'
	LENG6(4) = 28   ! ' ./001/COUNTSUR001_0001.out '

	LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1

	KKK = 0
	do
		OPEN(5511, FILE = '001/coll001_0001.dat', Status = 'OLD')
		read(5511,23456)  kK,x1,y1,z1,px1,py1
     &	,pz1,em1,iso1,isy1,id1,idr1
		if(em1.eq.0.9380)then
			KKK = 0
			PRINT*, 'Please add a filehead and save!'
		else
			KKK = 1
		endif
		close(5511)
		if(KKK.EQ.1)EXIT
		call SLEEP(1)
	enddo
23456 format(1x,i4,7f9.4,4i4)
      
      open(2004, file = 'ERR.out', status = 'UNKNOWN')
	do I = 1, NFILE
		call Folder001(I, FOLDER)
          call Folder0001(NTH_NOH_BIN(I), NNOHUP)
          
          ! cat 003/initangm.dat 003/initangm003_0002.dat ... > 003/initangm.dat
          write(2004,*) 'cat in file '//FOLDER//': '//NNOHUP
          CATCH1(1:1000000) = ' '
	    CATCH2(1:1000000) = ' '
	    CATCH3(1:1000000) = ' '
	    CATCH4(1:1000000) = ' '
110	    CATCH1(1:LENG0(1)) = 'cat'
	    CATCH2(1:LENG0(2)) = 'cat'
	    CATCH3(1:LENG0(3)) = 'cat'
	    CATCH4(1:LENG0(4)) = 'cat'
          LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1 ! 计算当前末位
		do J = 1, NTH_NOH_BIN(I)              
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
		enddo		
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
          do J = 1, NTH_NOH_BIN(I)          
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
		enddo
        write(2004,*) CATCH1(1:ENDING(1))
        call SYSTEM('rm -f ./'//FOLDER//'/initangm.dat')
        call SYSTEM('rm -f ./'//FOLDER//'/SURF_COAL.dat')
        call SYSTEM('rm -f ./'//FOLDER//'/coll.dat')
        call SYSTEM('rm -f ./'//FOLDER//'/COUNTSUR.out')
		call SYSTEM(CATCH1)
        call SYSTEM(CATCH2)
        call SYSTEM(CATCH3)
        call SYSTEM(CATCH4)
      enddo

      call SLEEP(3)

      ! cat 001/initangm.dat 002/initangm.dat ... >> results/initangm.dat
      write(2004,*) 'cat nfile: ', NFILE
      CATCH1(1:1000000) = ' '
	CATCH2(1:1000000) = ' '
	CATCH3(1:1000000) = ' '
	CATCH4(1:1000000) = ' '
      CATCH1(1:LENG0(1)) = 'cat'
	CATCH2(1:LENG0(2)) = 'cat'
	CATCH3(1:LENG0(3)) = 'cat'     
	CATCH4(1:LENG0(4)) = 'cat'
      LEAD(1:N_CATCH) = LENG0(1:N_CATCH) + 1
      !call SYSTEM('rm -f ./results/initangm.dat')
      !call SYSTEM('rm -f ./results/SURF_COAL.dat')
	!call SYSTEM('rm -f ./results/coll.out')
	!call SYSTEM('rm -f ./results/COUNTSUR.out')
      do I = 1, NFILE
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
        write(2004,*) CATCH1(1:ENDING(1))
	  call SYSTEM(CATCH1)
	  call SYSTEM(CATCH3)
        call SYSTEM('rm -f ./'//FOLDER//'/initangm.dat')
        call SYSTEM('rm -f ./'//FOLDER//'/coll.dat')
        if(SURF_COAL1.eq.1)then
          call SYSTEM(CATCH2)
          call SYSTEM(CATCH4)
          call SYSTEM('rm -f ./'//FOLDER//'/SURF_COAL.dat')
          call SYSTEM('rm -f ./'//FOLDER//'/COUNTSUR.out')
        endif
      enddo
      close(2004)
      call sleep(5)
********************************************************************
*** 删除多余的角动量数据，确保角动量计算准确
      open(5512,file='results/initangm.dat',status='old',iostat=ios)
	if(ios/=0)then
		write(*,*)'Error: Unable to open input file initangm.dat'
		stop
	endif

	open(5513,file='results/initangm1.dat',status='unknown')

	has_data=.false.

	do
		read(5512,*,iostat=ios)cur_kk,cur_x1,cur_y1,cur_z1

		if(ios<0)then
			! 到达文件末尾，把最后一个编号保存的数据写入文件
			if(has_data)then
				write(5513,*)last_kk,last_x1,last_y1,last_z1
			endif
			exit
		elseif(ios>0)then
			write(*,*)'Error reading file at KK =',cur_kk
			exit
		endif

		if(.not.has_data)then
			! 读取第一行数据时的初始化
			last_kk=cur_kk
			last_x1=cur_x1
			last_y1=cur_y1
			last_z1=cur_z1
			has_data=.true.
		elseif(cur_kk/=last_kk)then
			! 检测到新编号：说明last_kk的所有行已读取完毕，写入它收尾时的最后一行
			write(5513,*)last_kk,last_x1,last_y1,last_z1
			! 更新为新编号的数据
			last_kk=cur_kk
			last_x1=cur_x1
			last_y1=cur_y1
			last_z1=cur_z1
		else
			! 同一编号的后续行：不断用最新的x1,y1,z1覆盖旧数据
			last_x1=cur_x1
			last_y1=cur_y1
			last_z1=cur_z1
		endif
      enddo
      !close(5512) 
      !close(5513)
********************************************************************
      
	if(CALL_ANA.EQ.1)	then

		open(10, file = 'Input file/INPUT_ANA.in', status = 'UNKNOWN')
		NFILE_ANA = NFILE
		RUN_PER_ANA = NRUN_PER - FAIL
		NAME_ANA = 'ctrlA'
		T_PER_ANA = T_PER
		write(10,*)	NFILE_ANA
		write(10,*)	RUN_PER_ANA
		write(10,*)	TRIM(NAME)
		write(10,*)	T_PER_ANA
		write(10,*)	x00cut
		write(10,*) alpha, beta, gamma
		write(10,*) gsur, gtau, eta
		write(10,*)	gsym, symm
		write(10,*)	NINT(momp)		
		close(10)

		call SYSTEM('ifort -o '//TRIM(NAME_ANA)//' ctrlA.for')

		STOP	!终止程序，即只编译不运行

		call SYSTEM('nohup./ '//TRIM(NAME_ANA)//' >nohup_ana 2>&1')

		do

			open(10, file = 'MONITOR_ANA.dat', status = 'UNKNOWN')
			do
				read(10,*,end=339,ERR=339)	MONI
			enddo
339			close(10)

			call SLEEP(T_PER_ANA)

			open(10, file = 'MONITOR_ANA.dat', status = 'UNKNOWN')
			do
				read(10,*,end=340,ERR=340)	MONI1
			enddo
340			close(10)

			if(MONI1.EQ.MONI)	STOP

		enddo

      endif
      
      CONTAINS
      
********************************************************************
*** get the folder ID 001~100      
      subroutine Folder001(KK, FOLDER)
      
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
      subroutine Folder0001(KK, NNOHUP)
      
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

	END PROGRAM