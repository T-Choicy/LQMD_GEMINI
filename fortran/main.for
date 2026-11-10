
	PROGRAM	LQMD
      USE mass_selection, ONLY: read_mass_selection,
     &  overlay_selected_masses, mass_framework,
     &  report_mass_comparison_usage

!MAIN.1公共组定义与初始化----------------------------------------------------------------------------------------------
      ! 物理常数
      REAL, PARAMETER :: PI = 3.1415926           ! 圆周率 rad
      REAL, PARAMETER :: kC = 0.00140             ! 库伦势垒常数 fm GeV
      REAL, PARAMETER :: hbarc = 0.197327         ! hbar*c fm GeV
      
      INTEGER, PARAMETER :: SURF_MAX_A = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
	INTEGER, PARAMETER :: LEAD_MAX = 233
	INTEGER, PARAMETER :: N_EMIT_MAX = 40
	INTEGER, PARAMETER :: MAX_EMIT = 4
	character*50 outfile,targ, proj,iniout,rmsout, pmsout, rdeout,
     &	         aveout
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),pyms(500),
     &       pzms(500),emson(500),idms(500),icp(500),idanti(500)
      common/rp0/xi(500),yi(500),zi(500),pxi(500),pyi(500),pzi(500),
     &       emi(500),isoi(500),isyi(500),idi(500)
      common/rpi/rmst(500),pmst(500),avet(500)
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &      ncolln,ncollp,nblkn,nblkp,nanti
      common/bb/np,nzp,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/ee/mseed,nseed
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                           
      common/skyrme/al,sal,t0,t1,arang,vy,width                           
      common/gg/nzta,nzpr,nmint
      common/densi1/denstl(500),dens1(500),velcm(3)
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp
      common/deu2/px1,py1,pz1,px2,py2,pz2                                     
      common/deu3/dg(50,10,2),dtg(2),dtg1(2),dcm(50,10,2),ideu,dtgt(2),
     &      dtg1t(2)
      common/kaot/pktot,nkn 
      common/xx/alpha,beta,gamma,symm
      common/aven/tene
      common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc,numit
      common/veloc/velx(500),vely(500),velz(500)                              
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
	common/cmod/coeft,coefp,epp,ecm,ptcm
	common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
	common/data/targ, proj
	common/poli/iflagcro,iflagpau,ipoli
	common/totcol/ipuli3,iplin3,iplip3
      common/boite/ann(800),anp(800)
	common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     &	      frtwo(500)
      common/CLWF3/RTIME0
	common/CW/c1,c2
	common/wavel/sl(500)   
	common/shell/ish,ds,r0p,r0t,r0c,Eshpp,Eshpn,Eshtp,Eshtn,Eshcp,
     &             Eshcn,Ri0
	common/ths/th1,th2,th3
	common/prneck/ratio_np,nprogp,nprogn,ntargp,ntargn,rrn
	common/diss/xldis,Ediss
      common/medium/mdel,mdin,nsym,npop,nkaon
      common/SZPI/S10(600),S11(600),S01(600),S10D(600)
      common/SZWP/W10(600),W11(600),W01(600),W10D(600)
	common/matx/xmatr(4,0:400)
	common/signbd/xsNbDt(0:200), xsNbNl(0:200)
	common/dmss/XDMASS(500),XNM144(500),XNM153(500)
	common/psdc/nntag, ffnn
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),TOT_CLUST_NUM, 
     &			He4_ERROR
      common/TUNNEL/Z_RESID, A_RESID
	common/SURF1/ CLUST_LEAST_DIST, XCC_TARG, YCC_TARG, ZCC_TARG
     &,R0_SURF, D_SURF, COL_PATH, NDATA
      integer target_orient_method
      common /target_orient/ target_orient_method
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)
	common/EMIT_CANDID/EMIT_CANDID(STORE_MAX), 
     &		CLUST_BIN(2:SURF_MAX_A, STORE_MAX)
	common/SURF/TOO_MANY
	common/SURF_LOOP/NUM_RANK(SURF_MAX_A)


      dimension dtot(10),dtotcm(10)                                           
      dimension m1(100)
	dimension xncoll(2000),xncollp(2000),xncolln(2000),xnblk(2000),
     &          xnblkp(2000),xnblkn(2000)
	character date1*10,time1*10,zone1*8
	integer		values1(8),values(8),values2(8)
	integer		hour,minute,second,interval,OUTBARYON, OUTMESON
	REAL		SETIMP

	!表面并合
	INTEGER	SURF_COAL, K_START_SURF, K_END_SURF, NNTA, ZZTTA, TUNN
	REAL	XC_TARG, YC_TARG, ZC_TARG, R_TARG
	REAL	DIST_EX, RR0
	INTEGER	CLUST_CANDID(4), MMK
	REAL	RAN_E, NUCL_E, NUCL_R, E_CL, CLUST_LEAST_DIST_COEF
	REAL(KIND = 8)	E_TOT
	REAL	AGM(3), MOM(3), P_REFL, P_REFLX, P_REFLY, P_REFLZ
	REAL	CT_X, CT_Y, CT_Z
	INTEGER	LEADING_NUCLEON(LEAD_MAX), EMD_NUCL(500)
	INTEGER	CLUST_NUM(2:SURF_MAX_A), TEMP_CL, CL_NU(SURF_MAX_A)
	INTEGER	PRINT_CLUST, PRINT_POT
	LOGICAL	EMITTED
	INTEGER	N_NEU_CL, N_PRO_CL, NUM_CANDID, Z_RESID_TEMP, A_RESID_TEMP
	REAL	M_PRIM(3), ST_M(MAX_EMIT,3), M_LAST(3)
	INTEGER	Z_LAST, A_LAST, NP_STR, KKK, ST_P(MAX_EMIT)
	INTEGER	PROPER(0:STORE_MAX)
	INTEGER	TEST_1, TEST_2, TEST_3, TEST_4, TEST_5, TEST_6, TEST_7
	INTEGER	COUNT_COM(8)
	INTEGER	TEST_6HE, TEST_6LI, TEST_7LI, TEST_7BE, TEST_8HE, TEST_8LI
	INTEGER	TEST_8B
	INTEGER	TEST1_6HE, TEST1_6LI, TEST1_7LI, TEST1_7BE, TEST1_8HE
	INTEGER	 TEST1_8LI, TEST1_8B
	INTEGER	BINHE6(0:1000), BINLI6(0:1000)
	INTEGER	BINHE61(0:1000), BINLI61(0:1000)


	INTEGER	HE6, LI6, LI7, BE7, HE8, LI8, B8
	CHARACTER(LEN = 4) COLLCHAR
	CHARACTER(LEN = 9) SURFCHAR
	CHARACTER(LEN = 8) INITCHAR
	CHARACTER(LEN = 8) COUNTSUR
      CHARACTER(LEN = 7) MONICHAR
	CHARACTER(LEN = 20) COLLCHAR1
	CHARACTER(LEN = 25) SURFCHAR1
	CHARACTER(LEN = 24) INITCHAR1
	CHARACTER(LEN = 24) COUNTSUR1
      CHARACTER(LEN = 23) MONICHAR1
	CHARACTER(LEN = 4) dat_CHAR, out_CHAR
	INTEGER NTHR(3), NTHN(4), NTH_RUN, NTH_NOH
	CHARACTER(LEN = 3) NTHRCHAR
	CHARACTER(LEN = 4) NTHNCHAR
      INTEGER ERROR ! 缓存错误信息
      
	INTEGER	EMCL_N, EMCL_AA(N_EMIT_MAX), EMCL_ZZ(N_EMIT_MAX), TEMP_ZZ, 
     &		TEMP_NN, CL_I, J00, J3, EMCL_IF(N_EMIT_MAX), EMCL_NN, ZZ_1
	REAL	EMCL_E(N_EMIT_MAX), EMCL_M(N_EMIT_MAX, 3), 
     &		EMCL_A(N_EMIT_MAX, 3)
      INTEGER COUNT_2H, COUNT_3H, COUNT_3He, COUNT_4He
      
	CHARACTER*4 READ_ELS

      ERROR = 0
	HE6 = 0
	LI6 = 0
	LI7 = 0
	BE7 = 0
	HE8 = 0
	LI8 = 0
	B8 = 0

	TEST_1 = 0
	TEST_2 = 0
	TEST_3 = 0
	TEST_4 = 0
	TEST_5 = 0
	TEST_6 = 0
	TEST_7 = 0
	COUNT_COM(1:8) = 0

	E_HE6 = 0.
	E_LI6 = 0.
	BINHE6(0:1000) = 0
	BINLI6(0:1000) = 0
	BINHE61(0:1000) = 0
	BINLI61(0:1000) = 0


!MAIN.1公共组定义与初始化结束------------------------------------------------------------------------------------------

      CALL read_mass_selection
	open(1412, FILE='Input file/massfrdm.dat', status='old')

	do I=1,10000
  		read(1412,*,end=128,err=128) READ_IZ, READ_IA, READ_ELS, 
     &	READ_I1, READ_Be
		READ_XBIN(READ_IZ, READ_IA) = READ_Be ! 质量过剩(能)
	enddo
128	continue
	
	close(1412)
      IF (mass_framework == 0) THEN
        CALL overlay_selected_masses(READ_XBIN,1000,1000)
      ELSE
        WRITE(*,'(A)') 'MASS_COMPARE_SOURCE_FROZEN: READ_XBIN uses M0'
      ENDIF

!MAIN.2打开文件与日志输出----------------------------------------------------------------------------------------------
	!open(9,file='inf.out',status='unknown')
      open(2003,file='Input file/LQMD.in',status='old')
c**** arndt.in pi-n --> delta cross-section
      open(202,file='Input file/arndt.in',status='old')
	open(77, file='dmass.dat',status='unknown') !好像没有用 
	open(55,file='ImpPar.dat',status='unknown')
      if(SURF_COAL) open(20,file='error_log.txt',status='unknown') ! 存储错误信息
      !open(30,file='proj.dat',status='unknown')
      
	call date_and_time(date1,time1,zone1,values1)
	print 123, date1
	!write(9,123) date1
123	format(/t2,'the program begin at :',a10,1x,a8/)
!MAIN.2打开文件与日志输出结束------------------------------------------------------------------------------------------


!MAIN.3固定参数的设置--------------------------------------------------------------------------------------------------
	width=2.16
	vy1=-0.0024
	arang=0.8   
	xmn= 0.938
!MAIN.3固定参数的设置结束----------------------------------------------------------------------------------------------


!MAIN.4读入设置文件----------------------------------------------------------------------------------------------------
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
	read(2003,*) betap, betat, betac
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
      
	read(2003,*)	!表面并合输入
	read(2003,*)	SURF_COAL, PRINT_CLUST, PRINT_POT, NDATA
	read(2003,*)	K_START_SURF, K_END_SURF
	read(2003,*)	RR0, D_SURF, COL_PATH, CLUST_LEAST_DIST_COEF
	read(2003,*)	H0_SURF(1), H0_SURF(2), H0_SURF(3), H0_SURF(4),
     & H0_SURF(5)
	read(2003,*)	H1_SURF, INCL4_6
      call read_target_orient(2003,target_orient_method)
      if(target_orient_method.eq.2.and.initi.ne.1) then
        print *, "TARGET_ORIENT_METHOD 2 requires initi=1"
        stop 1
      endif
!MAIN.4读入设置文件结束------------------------------------------------------------------------------------------------

	COLLCHAR	= 'coll'        ! collision 碰撞数据
	SURFCHAR	= 'SURF_COAL'
	INITCHAR	= 'initangm'    ! initialization of angular momentum 角动量初始化
	COUNTSUR	= 'COUNTSUR'    ! count of surface coalescence 表面并合计数
      MONICHAR    = 'monitor'
	dat_CHAR		= '.dat'
	out_CHAR		= '.out'

	NTH_RUN1 = NTH_RUN
	NTHR(1) = FLOOR(nth_run/100.) + 48
	NTH_RUN = NTH_RUN - FLOOR(nth_run/100.)*100
	NTHR(2) = FLOOR(nth_run/10.) + 48
	NTH_RUN = NTH_RUN - FLOOR(nth_run/10.)*10
	NTHR(3) = FLOOR(nth_run/1.) + 48
	NTHRCHAR = CHAR(NTHR(1))//CHAR(NTHR(2))//CHAR(NTHR(3))

	NTHN(1) = FLOOR(nth_NOH/1000.) + 48
	NTH_NOH = NTH_NOH - FLOOR(nth_NOH/1000.)*1000
	NTHN(2) = FLOOR(nth_NOH/100.) + 48
	NTH_NOH = NTH_NOH - FLOOR(nth_NOH/100.)*100
	NTHN(3) = FLOOR(nth_NOH/10.) + 48
	NTH_NOH = NTH_NOH - FLOOR(nth_NOH/10.)*10
	NTHN(4) = FLOOR(nth_NOH/1.) + 48
	NTHNCHAR = CHAR(NTHN(1))//CHAR(NTHN(2))//CHAR(NTHN(3))//
     &			CHAR(NTHN(4))

	COLLCHAR1 = NTHRCHAR//'/'//COLLCHAR//NTHRCHAR//'_'//NTHNCHAR//
     &			dat_CHAR
	SURFCHAR1 = NTHRCHAR//'/'//SURFCHAR//NTHRCHAR//'_'//NTHNCHAR//
     &			dat_CHAR
	INITCHAR1 = NTHRCHAR//'/'//INITCHAR//NTHRCHAR//'_'//NTHNCHAR//
     &			dat_CHAR
	COUNTSUR1 = NTHRCHAR//'/'//COUNTSUR//NTHRCHAR//'_'//NTHNCHAR//
     &			dat_CHAR
      MONICHAR1 = NTHRCHAR//'/'//MONICHAR//NTHRCHAR//'_'//NTHNCHAR//
     &			dat_CHAR

	if(nth_run1.EQ.0)	then
		COLLCHAR1 = 'coll.dat'
		SURFCHAR1 = 'SURF_COAL.dat'
		INITCHAR1 = 'initangm.dat'
		COUNTSUR1 = 'COUNTSUR.dat'
        MONICHAR1 = 'MONITOR.dat'
      endif
      
	print*, COLLCHAR1
	print*, SURFCHAR1
	print*, INITCHAR1
	print*, COUNTSUR1
      print*, MONICHAR1

	if(SURF_COAL.EQ.1) open(78,FILE=SURFCHAR1,status='unknown')
	open(56,file=INITCHAR1,status='unknown')
      if(nth_run1.NE.0) open(89,file=MONICHAR1,status='unknown')
	if(SURF_COAL.EQ.1) open(57,file=COUNTSUR1,status='unknown')


!MAIN.5设置质量劈裂----------------------------------------------------------------------------------------------------
	if(momp.eq.1) then
		cmnn=cdelt*(1.+xsym)
		cmnp=cdelt*(1.-xsym)
	else
		cmnn=cmom*(1.+xsym)
		cmnp=cmom*(1.-xsym)
	endif
!MAIN.5设置质量劈裂----------------------------------------------------------------------------------------------------


!MAIN.5运行开头输出----------------------------------------------------------------------------------------------------
	print*, 'Projectile: A, Z',npr, nzpr, ' Target: A, Z',nta, nzta
	if(nrmin.gt.nrmax) stop 'Nrmin should be smaller than Nrmax'
!MAIN.5运行开头输出结束------------------------------------------------------------------------------------------------


!MAIN.6壳修正部分预备计算----------------------------------------------------------------------------------------------
	if(ish.eq.0)then                 !壳修正开关shell,ish=0 采用droplet model
		Eshpp=0.
		Eshpn=0.
		Eshtp=0.
		Eshtn=0.
		Eshcp=0.
		Eshcn=0.
	else							 !壳修正采用DTCSM，需要手动输入参数
        print*,'please input: Eshpp, Eshpn; Eshtp, Eshtn; Eshcp, Eshcn'
		read(*,*) Eshpp, Eshpn, Eshtp, Eshtn, Eshcp, Eshcn
      endif
      
c     弹核(projectile): 半径r0p, 形变参数betap, xr_proj, th1
c     靶核(target): 半径r0t, 形变参数betat, xr_targ, th2
c     复合核(compound nucleus): 半径r0c, 形变参数betac, xr_comp, th3
      
	r0p=r0*float(npr)**0.3333         
	r0t=r0*float(nta)**0.3333       
	r0c=r0*float(npr+nta)**0.3333     
	xr_proj=3.*betap**2/(8.*3.1415926)    !变形
	xr_targ=3.*betat**2/(8.*3.1415926)
	xr_comp=3.*betac**2/(8.*3.1415926)
	dis=0.32*r0
	th1=xr_proj/dis**2
	th2=xr_targ/dis**2
	th3=xr_comp/dis**2

!MAIN.6壳修正部分预备计算结束------------------------------------------------------------------------------------------


!MAIN.7关于核子核子碰撞产生共振与pi介子的设置的读入--------------------------------------------------------------------
	read(202,518) (S10(L),S11(L),S01(L),S10D(L),W10(L),W11(L),   
     &               W01(L),W10D(L),L=1,500)
      close(202)             !从arndt.in读入，delta cross-section
518	format(1X,8F10.4)

	call matrixs
	call NbarN_to_NbarDelta_dm
!MAIN.7关于核子核子碰撞产生共振与pi介子的设置的读入结束----------------------------------------------------------------


!MAIN.8输出文件设置----------------------------------------------------------------------------------------------------
	if(initi.eq.1)then
		proj='pini.dat'
		targ='tini.dat'
	endif

	outfile= COLLCHAR1  ! coll.out
!MAIN.8输出文件设置结束------------------------------------------------------------------------------------------------


!MAIN.9体系宏观特征的计算----------------------------------------------------------------------------------------------
 	! 1 GeV/u下，经典计算的速度明显超光速，因此要考虑相对论效应   CHEN SM
      ! if(npr.eq.1) epp=sqrt(xmn**2+epp**2)-xmn ! 将动量转换为动能   CHENG HG
	
      ! 质心系动能计算
      ecm = epp*float(npr*nta)/float(npr+nta)             ! 经典的质心系动能 CHENG HG
      print*, 'classical k.e. in c.m.: ', ecm
      
      ! 靶核初始四动量 p_target = (Eta1,0,0,0)
      RMASS = 0.938
      ETA1 = FLOAT(nta) * RMASS ! CHEN SM
      PZTA1 = 0.0
      ! 入射核初始四动量 p_projectile = (Epr,0,0,pzpr1)
      if(npr.NE.0)then
          EPR = FLOAT(npr)*(RMASS + epp)                  ! 总能量 Etot = Erest + T
          PZPR1 = SQRT(EPR**2 - (RMASS*npr)**2)            ! 纵向动量/动量大小 pz/p = (Etot^2 - Erest^2)^(1/2)
      else
          EPR = 0.0
          PZPR1 = 0.0
      endif
      S = (EPR + ETA1)**2 - PZPR1**2
      SQRT_S = SQRT(S)
      ecm1 = SQRT_S - ETA1 - FLOAT(npr)*RMASS               ! 相对论的质心系动能 CHEN SM
      print*, 'relativistic k.e. in c.m.: ', ecm1
	
      ! ======================================================================
	! 参数初始化与核反应模型相关物理量计算
	! ======================================================================

	fdelta = 0.5                     ! 表面弥散宽度调整因子（无量纲）

	al = 4.0*width                   ! 高斯型密度分布宽度参数（与核物质分布有关）
	sal = (al*PI)**1.5*rho0          ! 归一化因子，与密度分布积分相关（单位：质量/体积）

	Cex = -3./4.*kC*(3./PI)**0.3333  ! 库仑交换项系数（Skyrme 力参数）

	alnew = 4.*(width + fdelta**2./4.)  ! 考虑弥散后的新宽度参数
	salnew = 0.16*(alnew*PI)**1.5      ! 对应的新归一化因子

	if (ifdelt .lt. 1) then          ! 若开关 ifdelt < 1，则忽略弥散修正
		alnew = al
		salnew = sal
	endif

	t0 = -0.124                      ! Skyrme 力参数 t0（两体相互作用强度）
	t1 = 0.0705                      ! Skyrme 力参数 t1（动量依赖项）
	vy = vy1                         ! 粘性系数（或速度相关参数，vy1 由外部给定）

	! 电荷分布半径（基于质量数 A 的经验公式，1.142 为典型参数）
	rrp = 1.142*npr**0.3333          ! 入射核（弹核）电荷分布半径（fm）
	rrt = 1.142*nta**0.3333          ! 靶核电荷分布半径（fm）
	ptcm = rrp + rrt + rsur          ! 相互作用半径 = 两核半径之和 + 表面弥散修正（fm）

	np = npr + nta                   ! 体系总核子数（中子数+质子数）
	nzp = nzpr + nzta                ! 体系总质子数

	alp = 2./al                      ! 辅助参数，与 al 成反比
	bet = .346                       ! 经验常数（可能与核子-核子散射截面有关）

	! 与核内级联或蒸发相关的参数
	fdr = alp*bet/(alp + 2*bet)      ! 径向扩散系数（与动量输运相关）
	fdp = 1./4./(2*alp + bet)/0.197/0.197  ! 动量扩散系数（含 hbar*c 换算）

	fd1 = 8*(alp*bet)**1.5/(alp + 2*bet)**1.5/(2*alp + bet)**1.5

	xfdr = sqrt(3./2./fdr)           ! 径向均方位移特征值
	xfdp = sqrt(3./2./fdp)           ! 动量空间均方位移特征值

	! ======================================================================
	! 力学系数（用于核势能或流体动力学描述）
	! ======================================================================
	pota = 2.*alpha                  ! 势能展开第一项（alpha 为外部参数）
	potb = beta/(-alpha)*gamma/(gamma + 1.)   ! 势能第二项
	potc = gamma - 1.                ! 势能指数项
	potd = pota/sal/al               ! 体粘滞或密度依赖项系数（原密度）
	pote = potb/(sal)**potc          ! 指数型势项系数（原密度）
	potf = 0.00157/sal*2./al         ! 表面能相关系数
	potg = gradi/sal*rho0/(al/4)**2  ! 梯度项系数（gradi 为外部梯度参数）
	pots = 2.*symm/sal/al            ! 对称能系数（symm 为对称能强度）

	! 使用修正后的密度参数（若 ifdelt ≥ 1 则采用 alnew/salnew）
	potdnew = pota/salnew/alnew
	potenew = potb/salnew**potc
!MAIN.9体系宏观特征的计算结束------------------------------------------------------------------------------------------


!MAIN.10泡利初始化-----------------------------------------------------------------------------------------------------
      call initp 

      do 664 j = 1, 99999
		dr2 = float(j)/40. 
		dr = sqrt(dr2)
		ir = int(dr*100.+0.5)
		ih(j+1) = ir 
664   continue 

      ih(1) = 1
!MAIN.10泡利初始化结束-------------------------------------------------------------------------------------------------


!MAIN.11打开输出文件并写入体系特征-------------------------------------------------------------------------------------
	open(ifile,file=outfile,form='formatted',status='unknown')

	if(nth_run.EQ.0)	then
		if(nrmin.eq.1) write(ifile,1120) nt,x00min,x00max,epp,dt,al
		if(nrmin.eq.1) write(ifile,1130) nta,npr,nzta,nzpr
	elseif(nth_run.EQ.1.AND.nth_noh.EQ.1)	then
		if(nrmin.EQ.1)	then
			write(ifile,1120) nt,x00min,x00max,epp,dt,al
			write(ifile,1130) nta,npr,nzta,nzpr
		endif
	endif

1120  format(/1x,i4,5(f8.3,1x))
1130  format(1x,4(i3,1x))
!MAIN.11打开输出文件并写入体系特征结束---------------------------------------------------------------------------------

	print*, 'simulation of collision process!'


!MAIN.12碰撞参数循环---------------------------------------------------------------------------------------------------
      do 557 aimp=x00min,x00max,x00stp
          ! nx00=0就是固定碰撞参数，并且x00=x00min
          ! nx00=1则覆盖所给的碰撞参数范围
          if(nx00.eq.0) x00=aimp
          if(nx00.eq.0) print*, 'Ep=', epp, 'GeV  ', 'b=', x00, ' fm'
          if(nx00.gt.0) print*, 'Ep=', epp, 'GeV  ',  x00min, ' - ',
     &    x00max, ' fm'

          do ii=1,nt
              xncoll(ii)= 0.
              xncollp(ii)=0.
              xncolln(ii)=0.
              xnblk(ii)=  0.
              xnblkp(ii)= 0.
              xnblkn(ii)= 0.
          enddo

          call clallnew

          call coulgn

          !对非碰撞体系（自由演化），不计算运动学特征
          if(npr.eq.0.or.nta.eq.0.or.epp.eq.0.)then
              gam3=1.
              goto 707
          endif

		!计算动力学特征
		e1=xmn+epp                       ! 总能量                       
		pz1=sqrt(e1**2-xmn**2)           ! c.m.下的入射核子动量                                   
		bcm=pz1/(e1+float(nta)/float(npr)*xmn)                                
		bzz=bcm                                                                 
		gamcm=1./sqrt(1-bcm**2)
		gam3=gamcm                                                              
		pb=bcm*pz1                                                              
		gg=(gamcm-1)*pb/bcm/bcm-gamcm*e1                                        
		ppo=gg*bcm+pz1                                                          
		tpo=float(npr)/float(nta)*ppo 
707		continue

		!对碎裂反应（弹核很轻），采用靶核的几何截面
		if(npr.LE.nmint) then    !without initialization
			x00min=0.
			x00max=rrt
		endif

		do 4723 ijk = 1,nt 
4723			prtime(ijk) = 0

		do I=1, 500
			XDMASS(I)=0.
			XNM144(I)=0.
			XNM153(I)=0.
		enddo

		!碰撞参数随机数
		mseedd=10
		MSEEDDD=0.
		if(NRMIN.GT.1)	then
			do	I=1,NRMIN-1
				MSEEDDD=RAN(mseedd)
			enddo
		endif
		
		!GETX00CUT
		IFCOLL=1
		COLL1=0
		COLL2=0
		GETX00CUT_ROT=0

		!MAIN.12.1事件数循环---------------------------------------------------------------------
		do 556 nru=nrmin,nrmax



			TEST_6HE = 0
			TEST_6LI = 0
			TEST_7LI = 0
			TEST_7BE = 0
			TEST_8HE = 0
			TEST_8LI = 0
			TEST_8B = 0

			TEST1_6HE = 0
			TEST1_6LI = 0
			TEST1_7LI = 0
			TEST1_7BE = 0
			TEST1_8HE = 0
			TEST1_8LI = 0
			TEST1_8B = 0

			!write(9,*) 'Running NO: ', nru

			!事件数循环上的随机数种子
			mseed=28+nru*5

			!对于碰撞体系，且碰撞参数为非固定时（即nx00不为0），碰撞参数取遍范围
			if(nta.gt.0.and.nzta.gt.0.and.npr.gt.0.and.nzpr.gt.0.and.
     &epp.gt.0. .and. nx00.gt.0) then

				if(X00CUT.NE.0.)	then
                      ! 范围是[x00min, x00min+x00cut]，对于散裂反应，范围是[0, x00cut]
					x00=x00min+X00CUT*sqrt(RANFL(mseedd)) 
				else
					x00=x00min+(x00max-x00min)*sqrt(RANFL(mseedd))
				endif

			endif
	
			!GETX00CUT
			if(GETX00CUT.EQ.1.AND.IFCOLL.EQ.1)	then
				if(GETX00CUT_ROT.NE.0)	then
					GETX00CUT_IN = GETX00CUT_IN + GETX00CUT_ST
					print*
					print*,'ABANDON AND GET X00CUT AGAIN!'
					print*
				endif
				X00=GETX00CUT_IN
				GETX00CUT_IN=GETX00CUT_IN+GETX00CUT_ST
				
				COLL1=0
				COLL2=0
				IFCOLL=0
				GETX00CUT_ROT=0
			elseif(GETX00CUT.EQ.1.AND.GETX00CUT_ROT.LT.G_ROT)	then
				if(GETX00CUT_ROT.EQ.0) then
					print*
					print*, 'STEP ONE DONE'
					GETX00CUT_IN=GETX00CUT_IN-GETX00CUT_ST
				endif
				X00=GETX00CUT_IN
				GETX00CUT_ROT=GETX00CUT_ROT+1
				print*
				print*, 'ROT:', NINT(GETX00CUT_ROT)
			elseif(GETX00CUT.EQ.1.AND.GETX00CUT_ROT.EQ.G_ROT)	then
				print*
				print*, 'GETX00CUT DONE!'
				print*, 'X00CUT:',GETX00CUT_IN
				PAUSE
				GOTO	1245
			endif

			!显示并记录最终确定的碰撞参数
			print*, 'Running NO: ', nru, '  Impact param.: ', x00
			write(55,*)	nru,x00 ! fm
			
			!计算弹靶核相空间位置的库伦修正
			if(NPR.GT.0)	call coulomb_modify

			ideut=1                                                                 
			if(iflag6.eq.1) ideut=2 ! 循环两次的意义不明

			!MAIN.12.1.1 ideut循环-------------------------------------------
			do 5561 ideu=1,ideut

				do 8001 klk=1,2
					dtg(klk)=0                                                              
					dtg1(klk)=0                                                             
8001				continue 
                                                               
				if(ideu.eq.2) goto 5562
				do 8002 kkl=1,np                                                        
					iev(kkl)=1                                                              
8002				continue
				goto 5563                                                               
5562				call minspa         ! 计算iev                                                    
5563              continue
                                                          

				!选择构造初始构型的方式（给一个随机入射角度）
				if(initi.eq.0) call init 
				if(initi.eq.1) call init2 ! 散裂反应采用这个
				if(initi.eq.1) call width0


				ncoll=0
				ncolln=0
				ncollp=0
				nblock=0
				nblkn=0
				nblkp=0
				ipuli3=0
				iplin3=0
				iplip3=0
				numit =0
				rtime=0.
				call collnum
				call prot
				rtime0=10000.


				VANISH(1:500) = 0                   ! 结团核子发射状态参数
				LEADING_NUCLEON(1:LEAD_MAX) = 0
				EMD_NUCL(1:500) = 0                 ! 非结团核子发射状态参数

				EMCL_N = 0
				do JJJ = 1, 3
					EMCL_M(1:N_EMIT_MAX, JJJ) = 0.
					EMCL_A(1:N_EMIT_MAX, JJJ) = 0.
				enddo
				EMCL_E(1:N_EMIT_MAX) = 0.
				EMCL_AA(1:N_EMIT_MAX) = 0
				EMCL_ZZ(1:N_EMIT_MAX) = 0
				EMCL_IF(1:N_EMIT_MAX) = 0

                  COUNT_2H = 0
                  COUNT_3H = 0
                  COUNT_3He = 0
                  COUNT_4He = 0

				!MAIN.12.1.1.1 时域循环------------------
				do 26 k=1,nt

					!时域循环上随机数种子
					nseed=48+nru*8+3*k

					!计算当下时间点rtime，输出时的判断语句用到
					perc=float((nru-1)*nt+k)/float((nrmax-nrmin)*nt)*100.
					if(dt.gt.0.) rtime=k*dt
					if(dt.lt.0.) rtime=-k*dt

					!考虑表面并合
					if(SURF_COAL.EQ.1.AND.k.GE.K_START_SURF.AND.k.LE.K_END_SURF)then

						!求靶核的质心
						XC_TARG = 0.
						YC_TARG = 0.
						ZC_TARG = 0.
						do I_T = 1, NP
							if(idr(I_T).GT.0)		cycle
							if(VANISH(I_T).EQ.1)	cycle
							if(EMD_NUCL(I_T).EQ.1)	cycle
							XC_TARG = XC_TARG + X(I_T)
							YC_TARG = YC_TARG + Y(I_T)
							ZC_TARG = ZC_TARG + Z(I_T)
						enddo
						XC_TARG = XC_TARG/NTA
						YC_TARG = YC_TARG/NTA
						ZC_TARG = ZC_TARG/NTA

						!求靶核的均方根半径
						R_TARG = 0.	
						do I_T = 1, NP
							if(idr(I_T).GT.0)		cycle
							if(VANISH(I_T).EQ.1)	cycle
							if(EMD_NUCL(I_T).EQ.1)	cycle
						    R_TARG = R_TARG + (X(I_T) - XC_TARG)**2
							R_TARG = R_TARG + (Y(I_T) - YC_TARG)**2
							R_TARG = R_TARG + (Z(I_T) - ZC_TARG)**2
						enddo
						R_TARG = SQRT(R_TARG/NTA)

						!排除已经出射的核子后再求质心
						XCC_TARG = 0.
						YCC_TARG = 0.
						ZCC_TARG = 0.
						PXCC_TARG = 0.
						PYCC_TARG = 0.
						PZCC_TARG = 0.
						NNTA = 0
						ZZTTA = 0
						do I_T = 1, NP
							if(idr(I_T).GT.0)		cycle
							if(VANISH(I_T).EQ.1)	cycle
							if(EMD_NUCL(I_T).EQ.1)	cycle
							DIST_EX = 0.
							DIST_EX = DIST_EX + (X(I_T) - XC_TARG)**2
							DIST_EX = DIST_EX + (Y(I_T) - YC_TARG)**2
							DIST_EX = DIST_EX + (Z(I_T) - ZC_TARG)**2
							DIST_EX = SQRT(DIST_EX)
							DIST_EX = ABS(DIST_EX - R_TARG)
						
							if(DIST_EX.GE.rr0)	cycle
	
							PXCC_TARG = PXCC_TARG + PX(I_T)
							PYCC_TARG = PYCC_TARG + PY(I_T)
							PZCC_TARG = PZCC_TARG + PZ(I_T)
														
							XCC_TARG = XCC_TARG + X(I_T)
							YCC_TARG = YCC_TARG + Y(I_T)
							ZCC_TARG = ZCC_TARG + Z(I_T)
							NNTA = NNTA + 1

							if(ISY(I_T).EQ.-1)	ZZTTA = ZZTTA + 1
                          enddo
						
                          ! 用于存储预平衡发射后的剩余核信息
                          XCC_TARG = XCC_TARG/NNTA
						YCC_TARG = YCC_TARG/NNTA
						ZCC_TARG = ZCC_TARG/NNTA
						A_RESID = NNTA
						Z_RESID = ZZTTA

						! 实时求出R0_SURF
						R0_SURF = 1.4* NTA**0.333333
						CLUST_LEAST_DIST = 
     &						CLUST_LEAST_DIST_COEF*NTA**0.333333

						! 输出靶核势垒
						if(PRINT_POT.EQ.1)then

							open(21, FILE = 'NUCL_POT.dat', 
     &							STATUS = 'NEW')
							! 质子
							X(3) = XCC_TARG
							Y(3) = YCC_TARG
							Z(3) = ZCC_TARG
							PX(3) = 0.
							PY(3) = 0.
							PZ(3) = 0.
							ISO(3) = 1
							ISY(3) = -1
							WRITE(21,*) 'PROTON POTENTIAL'
							do STEP = 1, 800
								WRITE(21,*) X(3) - XCC_TARG, 
     &									GET_TOTAL_ENERGY(3)
								X(3) = X(3) + 0.05
							enddo
							WRITE(21,*)
							! 中子
							X(3) = XCC_TARG
							Y(3) = YCC_TARG
							Z(3) = ZCC_TARG
							PX(3) = 0.
							PY(3) = 0.
							PZ(3) = 0.
							ISO(3) = 0
							ISY(3) = 1
							WRITE(21,*), 'NEUTRON POTENTIAL'
							do STEP = 1, 400
								WRITE(21,*), X(3) - XCC_TARG, 
     &									GET_TOTAL_ENERGY(3)
								X(3) = X(3) + 0.05
							enddo

							PAUSE
							close(21)
							STOP
						endif

						! 检查是否有Leading nucleon
						do I_L = 1, NP

							if(VANISH(I_L).EQ.1)cycle
							if(idr(I_L).GT.0)cycle

							DIST_EX = 0. ! 临时变量

							! 只有径位在R0_SURF + D_SURF到R0_SURF + D_SURF + COL_PATH
							! 内的核子才可以作为LEADING核子
							DIST_EX = DIST_EX + (X(I_L) - XCC_TARG)**2
							DIST_EX = DIST_EX + (Y(I_L) - YCC_TARG)**2
							DIST_EX = DIST_EX + (Z(I_L) - ZCC_TARG)**2
                              DIST_EX = SQRT(DIST_EX) ! 基于剩余核的质心求得的核半径
							NUCL_R  = DIST_EX
							DIST_EX = DIST_EX - R0_SURF - D_SURF
							if(DIST_EX.LT.0.)	cycle ! R > R0_SURF + D_SURF
							DIST_EX = DIST_EX - COL_PATH
							NUCL_E = DIST_EX
							if(DIST_EX.GT.0.)	cycle ! R < R0_SURF + D_SURF + COL_PATH
					
							! 总能量大于0的才可作为leading nucleon
							DIST_EX = GET_TOTAL_ENERGY(I_L)-0.1
							E_TOT   = DIST_EX
							if(DIST_EX.LT.0.)	cycle
													
							! 远离中的核子才可作为leading nucleon
							CT_X = (X(I_L) - XCC_TARG)
							CT_Y = (Y(I_L) - YCC_TARG)
							CT_Z = (Z(I_L) - ZCC_TARG)
							DIST_EX = CT_X*PX(I_L) + CT_Y*PY(I_L) +
     &									CT_Z*PZ(I_L) ! 可能存在问题，r·p（p没有用质心系下的动量）
							if(DIST_EX.LT.0.)	cycle

							! 径位已经超级接近R0_SURF + D_SURF + COL_PATH的核子则作发射测试
							DIST_EX = DIST_EX/NUCL_R/SQRT(em(I_L)**2 + 
     &						PX(I_L)**2 + PY(I_L)**2 + PZ(I_L)**2)
							DIST_EX = NUCL_E + DIST_EX*DT
							if(DIST_EX.GE.0)	then								
								if(E_TOT.LT.50)	then
								call TUNNEL1(NINT((1-ISY(I_L))/2.), 1 
     &								, E_TOT, TUNN)
								elseif(E_TOT.GE.50)	then
									TUNN = 1
								endif

								if(TUNN.EQ.1)then	    ! 发射
									EMD_NUCL(I_L) = 1
									LEADING_NUCLEON(I_L) = 1
									cycle
								else					! 反射
									P_REFL = CT_X*PX(I_L) + 
     &	     							CT_Y*PY(I_L) + CT_Z*PZ(I_L)
									P_REFL = P_REFL/SQRT(CT_X**2 + 
     &	     							CT_Y**2 + CT_Z**2)
									P_REFLX = CT_X*P_REFL/NUCL_R
									P_REFLY = CT_Y*P_REFL/NUCL_R
									P_REFLZ = CT_Z*P_REFL/NUCL_R
									PX(I_L) = PX(I_L) - 2.*P_REFLX
									PY(I_L) = PY(I_L) - 2.*P_REFLY
									PZ(I_L) = PZ(I_L) - 2.*P_REFLZ
									LEADING_NUCLEON(I_L) = 0
									cycle
								endif
							endif
							
							! 让这个曾经的LEADING核子继续飞行，可以参与形成CLUST，但不再作为
							! LEADING核子，直到飞行到极限边缘，在那里进行发射检测，成功则发射，
							! 不成功则反射，反射后可再次作为LEADING核子
							if(LEADING_NUCLEON(I_L).EQ.1)	cycle
							LEADING_NUCLEON(I_L) = 1
							
                              ! 递归构造 recursive construction
							RECURS_CONSTR = 0
1211                      if(RECURS_CONSTR.EQ.1) then	!是一个超长的if，if在1211，else在1212，endif在1213	

							!设leading核子为0号
							TREE(0,2:SURF_MAX_A) = 0
							CLUST_STOR(0,1) = I_L
							CLUST_STOR(0,2:SURF_MAX_A) = 0

							!递归构造结团
							He4_ERROR = 0.
							CALL CONSTR_AND_TUNNEL(0)

							!按照核子数来归类和儲存结团
							do J0 = 2, SURF_MAX_A
								CLUST_BIN(J0, 1:STORE_MAX) = 0
								CLUST_NUM(J0) = 0
							enddo
							if(TOT_CLUST_NUM.NE.0) then

								do	J1 = 1, TOT_CLUST_NUM

									if(TREE(J1, 2).EQ.0)cycle

									do	J0 = 2, SURF_MAX_A

										J00 = SURF_MAX_A - J0 + 2
										if(TREE(J1, J00).NE.0)then
											CLUST_NUM(J00) = 
     &										CLUST_NUM(J00) + 1
											TEMP_CL = CLUST_NUM(J00)
											CLUST_BIN(J00, TEMP_CL) = 
     &										FLOAT(J1)
											EXIT
										endif

									enddo

								enddo

							endif

							!查看結團的構造情況
							if(PRINT_CLUST.EQ.1)then
							
								do J1 = 2, SURF_MAX_A				! 對結團核子數循環

									do J0 = 1, CLUST_NUM(J1)		! 對某一核子數的各結團循環

										CL_I = INT(CLUST_BIN(J1, J0))	! 獲取結團ID
										
										N_PRO_CL = 0
										N_NEU_CL = 0
										do J2 = 1, J1				! 獲取結團中的核子的ID
											CL_NU(J2) = 
     &										CLUST_STOR(CL_I, J2)
											if(ISY(CL_NU(J2)).EQ.-1)
     &										N_PRO_CL = N_PRO_CL + 1
											if(ISY(CL_NU(J2)).EQ.1)
     &										N_NEU_CL = N_PRO_CL + 1
										enddo

										CLUST_X = 0.				! 求結團幾何中心
										CLUST_Y = 0.
										CLUST_Z = 0.
										do J2 = 1, J1
											CLUST_X = CLUST_X +
     &										X(CL_NU(J2))
											CLUST_Y = CLUST_Y +
     &										Y(CL_NU(J2))
											CLUST_Z = CLUST_Z +
     &										Z(CL_NU(J2))
										enddo
										CLUST_X = CLUST_X/FLOAT(J1)
										CLUST_Y = CLUST_Y/FLOAT(J1)
										CLUST_Z = CLUST_Z/FLOAT(J1)

										CLUST_D=(CLUST_X-XCC_TARG)**2	!求結團離residue距離
     &										+ (CLUST_Y - YCC_TARG)**2
     &										+ (CLUST_Z - ZCC_TARG)**2
										CLUST_D = SQRT(CLUST_D)

										CLUST_R = 0.					!求結團尺寸
										do J2 = 1, J1
											CLUST_R = 
     &										(X(CL_NU(J2))-CLUST_X)**2+
     &										(Y(CL_NU(J2))-CLUST_Y)**2+
     &										(Z(CL_NU(J2))-CLUST_Z)**2
										enddo
										CLUST_R = SQRT(CLUST_R)
     &									/FLOAT(J1)

										CLUST_E = 0.					!求結團單核子能量之和
										do	J2 = 1, J1
											CLUST_E = CLUST_E +
     &									   GET_TOTAL_ENERGY(CL_NU(J2))
										enddo

										WRITE(*, 12563) CL_I, N_PRO_CL
     &									, J1, CLUST_D ,CLUST_E,CLUST_R

									enddo

								enddo

							endif

							!发射测试
							E_CL = 0.

                              ! 对结团核子数循环
							do J0 = 2, SURF_MAX_A

								J00 = SURF_MAX_A - J0 + 2 
								if(CLUST_NUM(J00).EQ.0)cycle
								call SORT_CANDID(J00,CLUST_NUM(J00)
     &                                , NUM_CANDID)

								COUNT_COM(J00) = COUNT_COM(J00) +
     &								 NUM_CANDID
                                  
                                  ! 对某一核子数的各结团循环(由大到小)
								do	II = 1, NUM_CANDID
									IIK = EMIT_CANDID(II)
									ZZ_1 = 0
                                      
                                      ! 计算结团质子数
									do	III = 1, J00
										IIJ = CLUST_STOR(IIK,III)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
                                          if(IIJ.LE.0)then ! 有待检查                                             
                                              write(20,*)'(Line:983)',
     &                                            IIJ  
                                              cycle  
                                          endif
										if(ISY(IIJ).EQ.-1)then
											ZZ_1 = ZZ_1 + 1
										endif
                                      enddo
                                      
                                      ! 记录构造的结团数量，真实发射的结团数量<=构造的结团数量
									if(J00.EQ.6)then
										if(ZZ_1.EQ.2)	 	
     &									    TEST_6HE = TEST_6HE + 1
										if(ZZ_1.EQ.3)	 	
     &									    TEST_6LI = TEST_6LI + 1
									endif
									if(J00.EQ.7)then
										if(ZZ_1.EQ.3)	 	
     &									    TEST_7LI = TEST_7LI + 1
										if(ZZ_1.EQ.4)		
     &                                        TEST_7BE = TEST_7BE + 1
									endif
									if(J00.EQ.8)then
										if(ZZ_1.EQ.2)	 	
     &									    TEST_8HE = TEST_8HE + 1
										if(ZZ_1.EQ.3)	 	
     &									    TEST_8LI = TEST_8LI + 1
										if(ZZ_1.EQ.5)	 	
     &									    TEST_8B = TEST_8B + 1
									endif
								enddo

								if(NUM_CANDID.EQ.0.)cycle
								call RANDOM_NUMBER(RAN_E)
								RAN_E = NUM_CANDID*(RAN_E + 0.00001)
                                  ! RAN_E的范围是[1,NUM_CANDID]
                                  ! 这意味着不对某核子数的各结团进行循环，而是随机取一种A=J00的结团                                  
								IIK	= EMIT_CANDID(CEILING(RAN_E))
								
                                  ! 计算当前结团质子数、当前剩余核质量数与质子数
                                  ZZ_1 = 0
                                  TEMP_ZZ = 0
                                  A_RESID_TEMP = A_RESID
								Z_RESID_TEMP = Z_RESID 
								do	III = 1, J00
                                      A_RESID = A_RESID - 1
                                      IIJ = CLUST_STOR(IIK, III)
                                      ! 检查是否溢出，溢出则跳过    CHEN SM
                                      if(IIJ.LE.0)then ! 有待检查                                             
                                          write(20,*) 
     &                                    'clust_id = 0(Line:1018)'         
                                          A_RESID = A_RESID_TEMP
                                          Z_RESID = Z_RESID_TEMP
                                          exit  
                                      endif
									if(ISY(IIJ).EQ.-1)then
										ZZ_1 = ZZ_1 + 1
                                          Z_RESID = Z_RESID - 1
										TEMP_ZZ = TEMP_ZZ + 1
									endif
                                  enddo
                                  
                                  ! 确保A_RESID>0, Z_RESID>=0
                                  if(A_RESID.LE.0)then
                                      write(20,*) 
     &                                    'A_R <= 0(Line:1045)'
                                      A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  elseif(Z_RESID.LT.0)then
                                      write(20,*) 
     &                                    'Z_R < 0(Line:1045)'
                                      A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  endif
                                  
                                  ! 若结团异常，则不进行发射
                                  if(A_RESID.EQ.A_RESID_TEMP)cycle
                                  
								if(J00.EQ.6)then
									if(ZZ_1.EQ.2)	 	
     &								    TEST1_6HE = TEST1_6HE + 1
									if(ZZ_1.EQ.3)	 	
     &								    TEST1_6LI = TEST1_6LI + 1
								endif
								if(J00.EQ.7)then
									if(ZZ_1.EQ.3)		
     &								    TEST1_7LI = TEST1_7LI + 1
									if(ZZ_1.EQ.4)		
     &								    TEST1_7BE = TEST1_7BE + 1
								endif
								if(J00.EQ.8)then
									if(ZZ_1.EQ.2)	 	
     &								    TEST1_8HE = TEST1_8HE + 1
									if(ZZ_1.EQ.3)		
     &								    TEST1_8LI = TEST1_8LI + 1
									if(ZZ_1.EQ.5)		
     &								    TEST1_8B = TEST1_8B + 1
                                  endif
                                  
								if(EMISSION_TEST(ZZ_1, J00, 
     &                                IIK, E_CL).EQ.1)then
									call GET_CLUST_KIN(IIK, E_CL, 
     &								MOM, AGM)
									do	JJ = 1, J00
                                          IIJ = CLUST_STOR(IIK, JJ)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
										if(IIJ.LE.0)then ! 有待检查                                             
											write(20,*)'(Line:1087)',
     &                                            IIJ      
                                              exit  
                                          endif
										VANISH(IIJ)= 1
									enddo
									EMCL_N = EMCL_N + 1
									EMCL_ZZ(EMCL_N) = TEMP_ZZ
									EMCL_AA(EMCL_N) = J00
									EMCL_E(EMCL_N) = E_CL
									EMCL_IF(EMCL_N) = 1
									EMCL_M(EMCL_N, 1:3) = MOM(1:3)
									EMCL_A(EMCL_N, 1:3) = AGM(1:3)
                                      if(ZZ_1.EQ.1)then
                                          if(J00.EQ.2)then
                                              COUNT_2H = COUNT_2H + 1
                                          elseif(J00.EQ.3)then
                                              COUNT_3H = COUNT_3H + 1
                                          endif
                                      elseif(ZZ_1.EQ.2)then
                                          if(J00.EQ.3)then
                                              COUNT_3He = COUNT_3He + 1
                                          elseif(J00.EQ.4)then
                                              COUNT_4He = COUNT_4He + 1
                                          endif
                                      endif
                                      
									EXIT
								else
									A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
								endif
							enddo

							!测试都失败以后，就让这个曾经的LEADING核子继续飞行，可以参与形成CLUST
							!直到飞行到极限边缘，在那里进行发射检测，成功则发射，不成功则反射


1212                      else	!是一个超长的if，if在1211，else在1212，endif在1213
							
							CLUST_STOR(0,1) = I_L
							CLUST_STOR(0,2:SURF_MAX_A) = 0
							do J1 = 1, STORE_MAX
								CLUST_STOR(J1, 1:SURF_MAX_A) = 0
							enddo
							do J1 = 1, SURF_MAX_A
								TREE(0:STORE_MAX ,J1) = 0
							enddo
							TOT_CLUST_NUM = 0.
							NUM_RANK(1) = 1				! A=1 结团数=1
							NUM_RANK(2:SURF_MAX_A) = 0
                              
                              ! 构造结团
							call CONSTR_LOOP                       

							!查看结团的构造情况
							if(PRINT_CLUST.EQ.1)then
								if(TOT_CLUST_NUM.NE.0.)print*, K
						
								do J1 = 2, SURF_MAX_A                 ! 对结团核子数循环

									do J0 = 1, NUM_RANK(J1)           ! 对某一核子数的各结团循环
	
										CL_I = TREE(J0, J1)	          ! 获取结团ID
									
										N_PRO_CL = 0
										N_NEU_CL = 0
										do J2 = 1, J1				  ! 获取结团中核子ID并计算结团总能量
											CL_NU(J2) =             
     &										CLUST_STOR(CL_I, J2)    
C 											print*, 	CL_NU(J2)
C 											PAUSE
		! ! 检查是否溢出，溢出则跳过    CHEN SM    CHEN SM
		! if (CL_NU(J2).LE.0) then ! 有待检查
		!  write(error_msg,'(A,I0,A,I0,A)')"Error: CLUST_STOR = ",
		!&    CL_NU(J2),"(Line:",1148,")"
		!  write(20,*)char(error_msg)  ! 写入错误日志文件
		!  print *, char(error_msg) ! 同时在屏幕输出
		!  cycle  
		! endif 
                                              if(ISY(CL_NU(J2)).EQ.-1)
     &										N_PRO_CL = N_PRO_CL + 1
											if(ISY(CL_NU(J2)).EQ.1)
     &										N_NEU_CL = N_PRO_CL + 1
										enddo

										CLUST_X = 0.                 ! 求结团几何中心
										CLUST_Y = 0.
										CLUST_Z = 0.
										do J2 = 1, J1
											CLUST_X = CLUST_X +
     &										X(CL_NU(J2))
											CLUST_Y = CLUST_Y +
     &										Y(CL_NU(J2))
											CLUST_Z = CLUST_Z +
     &										Z(CL_NU(J2))
										enddo
										CLUST_X = CLUST_X/FLOAT(J1)
										CLUST_Y = CLUST_Y/FLOAT(J1)
										CLUST_Z = CLUST_Z/FLOAT(J1)

										CLUST_D=(CLUST_X-XCC_TARG)**2	!求結團離residue距離
     &										+ (CLUST_Y - YCC_TARG)**2
     &										+ (CLUST_Z - ZCC_TARG)**2
										CLUST_D = SQRT(CLUST_D)

										CLUST_R = 0.					!求結團尺寸
										do J2 = 1, J1
											CLUST_R = 
     &										(X(CL_NU(J2))-CLUST_X)**2+
     &										(Y(CL_NU(J2))-CLUST_Y)**2+
     &										(Z(CL_NU(J2))-CLUST_Z)**2
										enddo
										CLUST_R = SQRT(CLUST_R)
     &									/FLOAT(J1)

										CLUST_E = 0.					! 求結團單核子能量之和
										do J2 = 1, J1
											CLUST_E = CLUST_E +
     &									   GET_TOTAL_ENERGY(CL_NU(J2))
										enddo

										WRITE(*, 12563) CL_I, N_PRO_CL
     &									, J1, CLUST_D ,CLUST_E,CLUST_R ! ID, proton number, mass number, distance bt. cluster and residue, energy, radius 

									enddo

								enddo

								if(TOT_CLUST_NUM.NE.0.)then
									print*
								endif

                              endif

							!发射测试
                              ! 结团总能量会在发射测试(EMISSION_TEST)中计算，然后返回到E_CL
							E_CL = 0. 
							VIRT = 0.
                              
                              ! 对结团核子数循环
							do J0 = 2, SURF_MAX_A                   

								if(VIRT.EQ.1)cycle
								J00 = SURF_MAX_A - J0 + 2             ! 结团核子数(由大到小)
								if(NUM_RANK(J00).EQ.0)cycle
								COUNT_COM(J00) = COUNT_COM(J00) +
     &								 NUM_RANK(J00)
                                  
                                  ! 对某一核子数的各结团循环
								do J1 = 1, NUM_RANK(J00)           

									IIK = TREE(J1, J00)               ! 结团ID
									ZZ_1 = 0                          ! 结团质子数
									CLUST_E = 0.                      ! 结团总能量
                                      
                                      ! 计算结团质子数与总能量
									do	III = 1, J00                
                                          
                                          IIJ = CLUST_STOR(IIK, III)    ! 结团中核子ID(第IIK结团的第III核子的ID)
                                          
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
										if(IIJ.LE.0)then ! 有待检查                                             
											write(20,*)'(Line:1241)',
     &                                            IIJ       
                                              exit  
                                          endif
                                          
										if(ISY(IIJ).EQ.-1)
     &									    ZZ_1 = ZZ_1 + 1
										CLUST_E = CLUST_E +
     &									    GET_TOTAL_ENERGY(IIJ)

                                      enddo

									!if(J00.EQ.6)	then
									!	if(ZZ_1.EQ.2)	print*, '6HE '
									!	if(ZZ_1.EQ.3)	print*, '6LI '
									!endif
									!if(J00.EQ.7)	then
									!	if(ZZ_1.EQ.3)	print*, '7LI'
									!	if(ZZ_1.EQ.4)	print*, '7BE'
									!endif
									!if(J00.EQ.8)	then
									!	if(ZZ_1.EQ.2)	print*, '8HE'
									!	if(ZZ_1.EQ.3)	print*, '8LI'
									!	if(ZZ_1.EQ.5)	print*, '8B'
									!endif

                                      ! 记录构造的结团数量，真实发射的结团数量<=构造的结团数量
    !                                 if(J00.EQ.6)then
   	!								if(ZZ_1.EQ.2)	
    !&                                        TEST_6HE = TEST_6HE + 1       
    !     								if(ZZ_1.EQ.3)
    !&                                        TEST_6LI = TEST_6LI + 1
    !									    if(ZZ_1.EQ.2.AND.CLUST_E.GT.0)
    !&			                            BINHE61(NINT(CLUST_E)) = 
    !&                                        BINHE61(NINT(CLUST_E)) + 1
    ! 								    if(ZZ_1.EQ.3.AND.CLUST_E.GT.0)
    !&			                            BINLI61(NINT(CLUST_E)) = 
    !&                                        BINLI61(NINT(CLUST_E)) + 1
    !									endif
    !									if(J00.EQ.7)then
    !										if(ZZ_1.EQ.3)
    !&                                        TEST_7LI = TEST_7LI + 1       
    ! 									if(ZZ_1.EQ.4)
    !&                                        TEST_7BE = TEST_7BE + 1
    !	 								endif
    ! 								if(J00.EQ.8)then
    !										if(ZZ_1.EQ.2)
    !&                                        TEST_8HE = TEST_8HE + 1       
    ! 									if(ZZ_1.EQ.3)
    !&                                        TEST_8LI = TEST_8LI + 1
    ! 									if(ZZ_1.EQ.5)
    !&                                        TEST_8B = TEST_8B + 1
    !                                 endif
                                      
								enddo

								call RANDOM_NUMBER(RAN_E)
								RAN_E = NUM_RANK(J00)*(RAN_E+0.00001) 
                                  RAN_E = CEILING(RAN_E)
                                  ! RAN_E的范围是[1,NUM_RANK(J00)]
                                  ! 这意味着不对某核子数的各结团进行循环，而是随机取一种A=J00的结团
								IIK	= TREE(RAN_E, J00)
								
                                  ! 计算当前结团质子数、当前剩余核质量数与质子数
                                  ZZ_1 = 0
                                  TEMP_ZZ = 0
                                  A_RESID_TEMP = A_RESID
								Z_RESID_TEMP = Z_RESID                                  
								do III = 1, J00
                                      A_RESID = A_RESID - 1
                                      IIJ = CLUST_STOR(IIK, III)
                                      ! 检查是否溢出，溢出则跳过    CHEN SM
									if(IIJ.LE.0)then ! 有待检查   
                                          write(20,*)'(Line:1320)',
     &                                        IIJ  
                                          A_RESID = A_RESID_TEMP
                                          Z_RESID = Z_RESID_TEMP
                                          exit  
                                      endif 
									if(ISY(IIJ).EQ.-1)then
                                          ZZ_1 = ZZ_1 + 1               ! 结团质子数
                                          TEMP_ZZ = TEMP_ZZ + 1
                                          Z_RESID = Z_RESID - 1         ! 剩余核质子数
                                      endif
                                  enddo
                                  
                                  ! 确保A_RESID>0, Z_RESID>=0
                                  ! 理论上不可能出现，一定是bug
                                  if(A_RESID.LE.0)then
                                      write(20,*) 
     &                                    'A_R <= 0(Line:1325)'
                                       A_RESID = A_RESID_TEMP
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  elseif(Z_RESID.LT.0)then
                                      write(20,*) 
     &                                    'Z_R < 0(Line:1325)'
                                      A_RESID = A_RESID_TEMP
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  endif
                                  
                                  ! 若结团异常，则不进行发射
                                  if(A_RESID.EQ.A_RESID_TEMP)cycle
                                  
								if(J00.EQ.6)then
                                      CLUST_E = 0.
                                      do III = 1, J00
                                          IIJ = CLUST_STOR(IIK, III)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
										if(IIJ.LE.0)then ! 有待检查                                             
											write(20,*)'(Line:1361)',
     &                                            IIJ           
                                              exit  
                                          endif 
                                          CLUST_E = CLUST_E +
     &                                        GET_TOTAL_ENERGY(IIJ)                                  
                                      enddo
                                      if(CLUST_E.LE.0) CLUST_E = 0.
									if(ZZ_1.EQ.2) TEST1_6HE = 	
     &								    TEST1_6HE + 1
									if(ZZ_1.EQ.3) TEST1_6LI = 	
     &								    TEST1_6LI + 1
									if(ZZ_1.EQ.2) print*, 'TEST 6HE'
									if(ZZ_1.EQ.3) print*, 'TEST 6LI'
									if(ZZ_1.EQ.2.AND.CLUST_E.GT.0)
     &				                    BINHE6(NINT(CLUST_E)) = 
     &                                    BINHE6(NINT(CLUST_E)) + 1
									if(ZZ_1.EQ.3.AND.CLUST_E.GT.0)
     &				                    BINLI6(NINT(CLUST_E)) = 
     &                                    BINLI6(NINT(CLUST_E)) + 1
								endif
								if(J00.EQ.7)then
									if(ZZ_1.EQ.3) TEST1_7LI = 	
     &								    TEST1_7LI + 1
									if(ZZ_1.EQ.4) TEST1_7BE = 	
     &								    TEST1_7BE + 1
									if(ZZ_1.EQ.3) print*, 'TEST 7LI'
									if(ZZ_1.EQ.4) print*, 'TEST 7BE'
								endif
								if(J00.EQ.8)then
									if(ZZ_1.EQ.2)	TEST1_8HE = 	
     &								TEST1_8HE + 1
									if(ZZ_1.EQ.3)	TEST1_8LI = 	
     &								TEST1_8LI + 1
									if(ZZ_1.EQ.5)	TEST1_8B = 	
     &								TEST1_8B + 1
									if(ZZ_1.EQ.2)	print*, 'TEST 8HE'
									if(ZZ_1.EQ.3)	print*, 'TEST 8LI'
									if(ZZ_1.EQ.5)	print*, 'TEST 8B'
								endif

C 								if(J00.EQ.6.or.J00.EQ.7.or.
C      &								J00.EQ.8)	PAUSE

                                  ! 进行发射测试
								if(EMISSION_TEST(ZZ_1, J00, 
     &                                IIK, E_CL).EQ.1)then
									call GET_CLUST_KIN(IIK, E_CL, 
     &								    MOM, AGM)
									do III = 1, J00
                                          IIJ = CLUST_STOR(IIK, III)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
									    if(IIJ.LE.0)then ! 有待检查                                             
										    write(20,*)'(Line:1415)',
     &                                            IIJ  
                                              Z_RESID = Z_RESID_TEMP
                                              A_RESID = A_RESID_TEMP
                                              ERROR = III
                                              exit
                                          endif 
										VANISH(IIJ) = 1
                                      enddo
                                      if(ERROR.NE.0)then
                                          do III = 1, ERROR
                                              IIJ = CLUST_STOR(IIK,III)
                                              VANISH(IIJ) = 0
                                          enddo
                                      endif

									EMCL_N = EMCL_N + 1
									EMCL_ZZ(EMCL_N) = ZZ_1
									EMCL_AA(EMCL_N) = J00
									EMCL_E(EMCL_N) = E_CL
									EMCL_IF(EMCL_N) = 1
									EMCL_M(EMCL_N, 1:3) = MOM(1:3)
									EMCL_A(EMCL_N, 1:3) = AGM(1:3)
                                      
                                      if(ZZ_1.EQ.1)then
                                          if(J00.EQ.2)then
                                              COUNT_2H = COUNT_2H + 1
                                          elseif(J00.EQ.3)then
                                              COUNT_3H = COUNT_3H + 1
                                          endif
                                      elseif(ZZ_1.EQ.2)then
                                          if(J00.EQ.3)then
                                              COUNT_3He = COUNT_3He + 1
                                          elseif(J00.EQ.4)then
                                              COUNT_4He = COUNT_4He + 1
                                          endif
                                      endif
									
                                      EXIT
                                  else
                                      Z_RESID = Z_RESID_TEMP
                                      A_RESID = A_RESID_TEMP
								endif

							enddo


							VIRTU1 = 10000000.
							do J0 = 2, SURF_MAX_A

								if(VIRT.NE.1.OR.TOT_CLUST_NUM.EQ.0.)
     &								cycle

								J00 = SURF_MAX_A - J0 + 2
								if(NUM_RANK(J00).EQ.0)	cycle

								do J1 = 1, NUM_RANK(J00)

									IIK = TREE(J1, J00)
									ZZ_1 = 0					
									SPX = 0.
									SPY = 0.
									SPZ = 0.
									SEE = 0.
                                      
                                      ! 计算结团质子数、动量、静能
									do III = 1, J00
										IIJ = CLUST_STOR(IIK, III)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
									    if(IIJ.LE.0)then ! 有待检查                                             
										   write(20,*)'(Line:1462)',
     &                                            IIJ 
                                          endif 
										if(ISY(INT(IIJ)).EQ.-1)	 
     &									    ZZ_1 = ZZ_1 +1            ! 结团质子数
										    SPX = SPX + PX(IIJ)       ! 结团动量
										    SPY = SPY + PY(IIJ)
										    SPZ = SPZ + PZ(IIJ)
										SEE = SEE + SQRT(0.938**2 +   
     &							PX(IIJ)**2 + PY(IIJ)**2 + PZ(IIJ)**2) 
                                      enddo
									SSCL = SEE**2 - SPX**2 
     &										- SPY**2 - SPZ**2
									if(SSCL.LE.0.) SSCL=(J00*0.938)**2
									SSCL = SQRT(SSCL)                 ! 结团静能(GeV)

									VIRTU = 1000.*(SSCL - J00*0.938) - 
     &								    READ_XBIN(ZZ_1, J00)          ! 计算的质量过剩能与实际的质量过剩能之差
									VIRTU = VIRTU/FLOAT(J00)
									if(VIRTU.LE.VIRTU1)then           ! 不知道用来干嘛
										VIRTU1 = VIRTU
										IIK1 = IIK
										Z_SEL = ZZ_1
										A_SEL = J00
									endif

								enddo
									
							enddo

							if(VIRT.EQ.1.AND.TOT_CLUST_NUM.NE.0.
     &						.AND.IIK1.NE.0.)then
							
								IIK = IIK1
								ZZ_1 = 0
								do	III = 1, J00 ! 有问题，J00没有外循环，一定等于2   CHEN SM
									IIJ = CLUST_STOR(IIK,III)
                                      ! 检查是否溢出，溢出则跳过    CHEN SM
                                      if(IIJ.LE.0)then ! 有待检查                                             
                                          write(20,*)'(Line:1501)',
     &                                        IIJ
                                      endif 
									if(ISY(IIJ).EQ.-1)then
										ZZ_1 = ZZ_1 + 1
                                          TEMP_ZZ = TEMP_ZZ + 1
                                          Z_RESID = Z_RESID - 1
									endif
								enddo

								if(J00.EQ.6)then
                                      CLUST_E = 0.
                                      do J2 = 1, J00
                                          CLUST_E = CLUST_E +
     &                                        GET_TOTAL_ENERGY
     &                                        (INT(CLUST_STOR(IIK, J2)))
                                      enddo
                                      if(CLUST_E.LE.0)CLUST_E =0.
									if(ZZ_1.EQ.2) TEST1_6HE = 	
     &								    TEST1_6HE + 1
									if(ZZ_1.EQ.3) TEST1_6LI = 	
     &								    TEST1_6LI + 1
									if(ZZ_1.EQ.2) print*, 'TEST 6HE'
									if(ZZ_1.EQ.3) print*, 'TEST 6LI'
									if(ZZ_1.EQ.2.AND.CLUST_E.GT.0)
     &				                    BINHE6(NINT(CLUST_E)) = 
     &                                    BINHE6(NINT(CLUST_E)) + 1
									if(ZZ_1.EQ.3.AND.CLUST_E.GT.0)
     &				                    BINLI6(NINT(CLUST_E)) = 
     &                                    BINLI6(NINT(CLUST_E)) + 1
								endif
								if(J00.EQ.7)	then
									if(ZZ_1.EQ.3)	TEST1_7LI = 	
     &								TEST1_7LI + 1
									if(ZZ_1.EQ.4)	TEST1_7BE = 	
     &								TEST1_7BE + 1
									if(ZZ_1.EQ.3)	print*, 'TEST 7LI'
									if(ZZ_1.EQ.4)	print*, 'TEST 7BE'
								endif
								if(J00.EQ.8)	then
									if(ZZ_1.EQ.2)	TEST1_8HE = 	
     &								TEST1_8HE + 1
									if(ZZ_1.EQ.3)	TEST1_8LI = 	
     &								TEST1_8LI + 1
									if(ZZ_1.EQ.5)	TEST1_8B = 	
     &								TEST1_8B + 1
									if(ZZ_1.EQ.2)	print*, 'TEST 8HE'
									if(ZZ_1.EQ.3)	print*, 'TEST 8LI'
									if(ZZ_1.EQ.5)	print*, 'TEST 8B'
								endif

C 								if(J00.EQ.6.or.J00.EQ.7.or.
C      &								J00.EQ.8)	PAUSE
				
								A_RESID = NNTA - J00
								Z_RESID_TEMP = Z_RESID
								Z_RESID = Z_SEL
                                  
                                  ! 确保A_RESID>0, Z_RESID>=0
                                  if(A_RESID.LE.0)then
                                      write(20,*) 
     &                                    'A_R <= 0(Line:1539)'
                                      A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  elseif(Z_RESID.LT.0)then
                                      write(20,*) 
     &                                    'Z_R < 0(Line:1539)'
                                      A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
                                      cycle
                                  endif
                                  
								if(EMISSION_TEST(ZZ_1, J00, 
     &                                IIK, E_CL).EQ.1)then ! CHEN SM
									call GET_CLUST_KIN(IIK, E_CL, 
     &								    MOM, AGM)
                                      ERROR = 0
									do JJ = 1, J00
                                          IIJ = CLUST_STOR(IIK,JJ)
                                          ! 检查是否溢出，溢出则跳过    CHEN SM
									    if(IIJ.LE.0)then ! 有待检查                                             
                                              write(20,*)'(Line:1578)',
     &                                            IIJ
                                              Z_RESID = Z_RESID_TEMP
                                              A_RESID = A_RESID_TEMP
                                              ERROR = JJ
                                              exit
                                          endif 
										VANISH(IIJ) = 1
                                      enddo
                                      if(ERROR.NE.0)then
                                          do JJ = 1, ERROR
                                              IIJ = CLUST_STOR(IIK,JJ)
                                              VANISH(IIJ) = 0
                                          enddo
                                      endif
                                      
									EMCL_N = EMCL_N + 1
									EMCL_ZZ(EMCL_N) = ZZ_1 ! 原为Z_TEMP, CHEN SM
									EMCL_AA(EMCL_N) = J00
									EMCL_E(EMCL_N) = E_CL
									EMCL_IF(EMCL_N) = 1
									EMCL_M(EMCL_N, 1:3) = MOM(1:3)
									EMCL_A(EMCL_N, 1:3) = AGM(1:3)
                                      if(ZZ_1.EQ.1)then
                                          if(J00.EQ.2)then
                                              COUNT_2H = COUNT_2H + 1
                                          elseif(J00.EQ.3)then
                                              COUNT_3H = COUNT_3H + 1
                                          endif
                                      elseif(ZZ_1.EQ.2)then
                                          if(J00.EQ.3)then
                                              COUNT_3He = COUNT_3He + 1
                                          elseif(J00.EQ.4)then
                                              COUNT_4He = COUNT_4He + 1
                                          endif
                                      endif
									EXIT
								else
									A_RESID = A_RESID + J00
									Z_RESID = Z_RESID_TEMP
								endif
								
							endif

1213                      endif   !是一个超长的if，if在1211，else在1212，endif在1213

                          enddo

                          !WRITE(57,*) 'timep   d    t    h    a'
		                WRITE(57,888)k, COUNT_2H, COUNT_3H,
     &	                    COUNT_3He, COUNT_4He
888                       FORMAT(1I8,4I5)

                          ! Strong decay after freeze-out (instantaneous)
						if(K.EQ.K_END_SURF)then
							!结团粒子的强衰变
							KKK = EMCL_N
							do	JJJ = 1, EMCL_N
								M_PRIM(1:3) = EMCL_M(JJJ, 1:3)
								call CLUSTER_STRONG_DECAY(EMCL_ZZ(JJJ)
     &							, EMCL_AA(JJJ), M_PRIM, Z_LAST, A_LAST 
     &							,M_LAST, ST_P, ST_M, NP_STR)
								
								SHUT1=1.
								if(NP_STR.NE.0.AND.SHUT1.EQ.0)then
									EMCL_ZZ(JJJ) = Z_LAST
									EMCL_AA(JJJ) = A_LAST
									EMCL_M(JJJ, 1:3) = M_LAST(1:3)
									EMCL_E(JJJ) = GET_ACC_EKIN(Z_LAST, 
     &								A_LAST, M_LAST)
     									do	KK = 1, NP_STR

										KKK = KKK + 1
										EMCL_ZZ(KKK) =  ST_P(KK)
										EMCL_AA(KKK) =  NINT(3.*(ST_P(
     &									KK)**2 - ST_P(KK))/2. + 1.)
										EMCL_M(KKK, 1:3) = ST_M(KK,
     &									 1:3)

										M_PRIM(1:3) = EMCL_M(KKK, 1:3)
										EMCL_E(KKK) = GET_ACC_EKIN(
     &									EMCL_ZZ(KKK), EMCL_AA(KKK)
     &									,M_PRIM)

										EMCL_A(KKK, 1:3) = 0.
										EMCL_IF(KKK) = 1					
									enddo
								endif
                              enddo
                              EMCL_N = KKK
						!endif  ! chen sm

							do JJJ = 1, EMCL_N
								WRITE(78, 23457) EMCL_ZZ(JJJ), 
     &							EMCL_AA(JJJ), 
     &							EMCL_E(JJJ), EMCL_M(JJJ, 1), 
     &							EMCL_M(JJJ,2), EMCL_M(JJJ, 3),
     &							EMCL_A(JJJ, 1), EMCL_A(JJJ, 2),
     &							EMCL_A(JJJ, 3)
                              enddo
							!if(mod(k,iwrite).eq.0) ! chen sm
                              WRITE(78, 23458) 0, 0, NRU, 0., 0., 0., 
     &											0., 0., 0.
						endif


                          endif                          
! 表面并合判断结束-----------------------------------------------------------------------------------------------------

55213					CONTINUE

5512					format(1x, 6f9.4, 1X, I2)
23458					format(1x, i1, 1X, I1, I10, 6f9.4)
23457					format(1x,i1, 1X, I1, f9.2, 6f9.4)





					!若iflag1<0，则不处理碰撞过程
					if(iflag1.lt.0) goto 2223

					if(iflag4 .eq. 1) then
						if(iflagcro.eq.0) then		!THIS
							call rcollbb
							call resdecay
							call rcollmeson
							call mesonabs
						elseif(iflagcro.eq.1) then
							call crossnew
							call rcollnew
						endif
					elseif(iflag4.eq.2) then
						call nrcoll(crsc)
					endif

C 					if(NRU.EQ.2200)	print*, 'GGG'

					xncoll(k)=xncoll(k)+float(ncoll)/float(nrmax-
     &				nrmin+1)
					xncollp(k)=xncollp(k)+float(ncollp)/float(nrmax-
     &				nrmin+1)
					xncolln(k)=xncolln(k)+float(ncolln)/float(nrmax-
     &				nrmin+1)
					xnblk(k)=xnblk(k)+float(nblock)/float(nrmax-
     &				nrmin+1)
					xnblkp(k)=xnblkp(k)+float(nblkp)/float(nrmax-
     &				nrmin+1)
					xnblkn(k)=xnblkn(k)+float(nblkn)/float(nrmax-
     &				nrmin+1)

					!相空间的平均场演化
2223					if(iflag2.eq.0)                   call	updt
					if(iflag2.eq.1 .and. iflag5.eq.0 .and. iflag3.ne.
     &				3)								  call	updtyc
					if(iflag2.eq.1 .and. iflag5.eq.1 .and. iflag3.ne.
     &				3)								  call	updtv

					if(iflag3.eq.3)					  call	updtim	! mean field of baryons
					if(iflag3.eq.3)					  call	udtmsn	! mean field of mesons
					if(iflag3.eq.3 .and. iflag7.eq.1) call  updtsh	! 在只有靶核的时候不可开壳修正，否则出现sqrt DOMAIN的错误		
C 					if(iflag3.eq.3 .and. iflag7.eq.1) call  constr
					if(iflag3.eq.3 .and. iflag7.eq.1 .and. ipot.eq.
     &				1)								  call	neck

C 					if(NRU.EQ.2200)	print*, 'HHH'
					
					!散裂反应中求初始角动量
                      if(K.EQ.2)	then
					XKTA=0.
					YKTA=0.
					ZKTA=0.
					XXTA=0.
					YYTA=0.
					ZZTA=0.
					do II=NPR+1, NP
						XKTA=XKTA+PX(II)
						YKTA=YKTA+PY(II)
						ZKTA=ZKTA+PZ(II)
						XXTA=XXTA+X(II)
						YYTA=YYTA+Y(II)
						ZZTA=ZZTA+Z(II)
					enddo
					XXTA=XXTA/NTA
					YYTA=YYTA/NTA
					ZZTA=ZZTA/NTA

					xkpr=0.
					ykpr=0.
					zkpr=0.
					xxpr=0.
					yypr=0.
					zzpr=0.
					do ii=1,npr
						xkpr=xkpr+px(ii)
						ykpr=ykpr+py(ii)
						zkpr=zkpr+pz(ii)
						xxpr=xxpr+x(ii)
						yypr=yypr+y(ii)
						zzpr=zzpr+z(ii)
					enddo
					xxpr=xxpr/npr
					yypr=yypr/npr
					zzpr=zzpr/npr
					
					totangmx=yyta*zkta-zzta*ykta+yypr*zkpr-zzpr*ykpr
					totangmy=zzta*xkta-xxta*zkta+zzpr*xkpr-xxpr*zkpr
					totangmz=xxta*ykta-yyta*xkta+xxpr*ykpr-yypr*xkpr
					totangmx=totangmx*5.067731
					totangmy=totangmy*5.067731
					totangmz=totangmz*5.067731
    					
              write(56,*) nru, totangmx, totangmy,
     &					totangmz					
                      endif
											
     !!         do I = 1, npr
     !!           E_lab = sqrt(em(i)**2 + px(i)**2 + py(i)**2 + pz(i)**2)
     !!           write(30,23456) k,x(i),y(i),z(i)
     !!&            ,px(i),py(i),pz(i),E_lab,
     !!&            iso(i),isy(i),id(i),idr(i),
     !!&            NINT(VANISH(I)), EMD_NUCL(I)                
     !!         end do
              
					if(mod(k,iwrite).eq.0)	then
						do 81	I=1,NP
							if(rtime.GE.ITMIN.AND.rtime.LE.ITMAX) then
								if(OUTBARYON.EQ.0)	GOTO	2659
								write(ifile,23456) k,x(i),y(i),z(i)
     &							,px(i),py(i),pz(i),em(i),
     &							iso(i),isy(i),id(i),idr(i),
     &							NINT(VANISH(I)), EMD_NUCL(I)
c k是时间点
c iso
c isy: -1是质子,1是中子, (1-isy)/2即核电荷数
c
                                  
2659								if(OUTMESON.EQ.0)	GOTO	81
								write(ifile,23456) k,xms(i),yms(i),
     &							zms(i),pxms(i),pyms(i),pzms(i),
     &							emson(i),idms(i),icp(i)
							endif
81						CONTINUE
					endif
23456					format(1x,i4,7f9.4,6i4)

					!GETX00CUT
					if(GETX00CUT.EQ.1)	then

						if(K.EQ.1)	then
							do	I=1, NP
								COLL1=COLL1+ID(I)
							enddo
						endif



						if(K.EQ.NT)	then

							do I=1, NP
								COLL2=COLL2+ID(I)
							enddo
							
							if(NINT(COLL1).NE.NINT(COLL2))	IFCOLL=1

						endif

                      endif
                      
                      if(nth_run1.NE.0) write(89, *) nru, k
                      
26                continue
                  !MAIN.12.1.1.1 时域循环结束--------------
                  
5561			continue
			!MAIN.12.1.1 ideut循环结束---------------------------------------
			
!		if((mod(nru,(nrmax-nrmin+1)).EQ.1).AND.NTH_RUN.EQ.0) 
!     &    WRITE(57,*) 'nru 6He 6Li 7Li 7Be 8He 8Li 8B'
!		WRITE(57,888)	nru, TEST_6HE, TEST_6LI, TEST_7LI, TEST_7BE,
!     &	    TEST_8HE, TEST_8LI, TEST_8B
!888       FORMAT(1I8,7I5)
          
556       continue
		!MAIN.12.1事件数循环结束-----------------------------------------------------------------
		
  !        print*, 'After sorting a:'
		!do	I = 2, 8
		!	print*, 'a = ', I, COUNT_COM(I)
		!enddo
		!
		!print*
		!print*, 'after sorting iso:'
		!print*, 'HE6 = ', TEST_6HE
		!print*, 'LI6 = ', TEST_6LI
		!print*, 'LI7 = ', TEST_7LI
		!print*, 'BE7 = ', TEST_7BE
		!print*, 'HE8 = ', TEST_8HE
		!print*, 'LI8 = ', TEST_8LI
		!print*, 'B8 = ', TEST_8B
  !
		!print*
		!print*, 'TESTED iso:'
		!print*, 'HE6 = ', TEST1_6HE
		!print*, 'LI6 = ', TEST1_6LI
		!print*, 'LI7 = ', TEST1_7LI
		!print*, 'BE7 = ', TEST1_7BE
		!print*, 'HE8 = ', TEST1_8HE
		!print*, 'LI8 = ', TEST1_8LI
		!print*, 'B8 = ', TEST1_8B
  !
		!do II = 1, 1000
		!	if(BINHE61(II).EQ.0)	cycle
		!	print*, II, ' FORM ', BINHE61(II)
		!enddo
		!print*
		!do II = 1, 1000
		!	if(BINLI61(II).EQ.0)	cycle
		!	print*, II, ' FORM ', BINLI61(II)
		!enddo
  !
		!print*
  !
		!do II = 1, 1000
		!	if(BINHE6(II).EQ.0)	cycle
		!	print*, II, ' TEST ', BINHE6(II)
		!enddo
		!print*
		!do II = 1, 1000
		!	if(BINLI6(II).EQ.0)	cycle
		!	print*, II, ' TEST ', BINLI6(II)
		!enddo

C 		PAUSE
		

557   continue
!MAIN.12碰撞参数循环结束-----------------------------------------------------------------------------------------------


!MAIN.13程序结束善后---------------------------------------------------------------------------------------------------
      call date_and_time(date1,time1,zone1,values2)
      hour=values2(5)-values1(5)
	if(hour.le.0) hour=24+hour

	minute=values2(6)-values1(6)
	second=values2(7)-values1(7)
	interval=60*(hour*60+minute)+second
	hour=interval/3600
	minute=(interval-hour*3600)/60

	write(*,124) date1, time1
124	format(/t2,'the program end at ',a10,1x,a8/)
12563	FORMAT('CL_No.=', I4, ' P=', I1, ' A=', I1, 
     &' D=', F4.1, ' E=', F7.1, ' R=', F3.1)

1245	print*, 'The program end normally, ok!'
!MAIN.13程序结束善后完成-----------------------------------------------------------------------------------------------
      CALL report_mass_comparison_usage()
	end

      include 'target_orientation.for'
