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


	integer :: cur_kk,last_kk,ios
      real :: cur_x1,cur_y1,cur_z1
      real :: last_x1,last_y1,last_z1
      logical :: has_data

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
				!到达文件末尾，把最后一个编号保存的数据写入文件
				if(has_data)then
					write(5513,*)last_kk,last_x1,last_y1,last_z1
				endif
				exit
			elseif(ios>0)then
				write(*,*)'Error reading file at KK =',cur_kk
				exit
			endif

			if(.not.has_data)then
				!读取第一行数据时的初始化
				last_kk=cur_kk
				last_x1=cur_x1
				last_y1=cur_y1
				last_z1=cur_z1
				has_data=.true.
			elseif(cur_kk/=last_kk)then
				!检测到新编号：说明last_kk的所有行已读取完毕，写入它收尾时的最后一行
				write(5512,*)last_kk,last_x1,last_y1,last_z1
				!更新为新编号的数据
				last_kk=cur_kk
				last_x1=cur_x1
				last_y1=cur_y1
				last_z1=cur_z1
			else
				!同一编号的后续行：不断用最新的x1,y1,z1覆盖旧数据
				last_x1=cur_x1
				last_y1=cur_y1
				last_z1=cur_z1
			endif
		enddo
      END program