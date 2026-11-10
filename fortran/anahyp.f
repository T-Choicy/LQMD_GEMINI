! CSM: Restored garbled comments in this source file; program logic unchanged.
************************************************************************
ccccc reconstruction of fragments with Coalescence Model and GEMINI model
ccccc in  C.M. Frame
************************************************************************
C*      R0+P0+Iso.(isospin)+RAM(Restructured Aggregation Model)       *C
C----------------------------------------------------------------------C
C*      Iso= 1, considering isospin, 0 no considering isospin.        *C
C*      Iram=1, RAM is included,     0 no.                            *C
ccccc     This program is the routine of Coalescence Model        cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc    MZ       : is the number of protons in the system ;      cccccc
ccccc    MN       : is the number of neutrons in the system;      cccccc
ccccc    mass     : is the number of nucleons in the system;      cccccc
ccccc    NEVENT   : is the total number of events;                cccccc
ccccc    rr0, pp0 : are the parameters of Coalescence Model.      cccccc
ccccc    MG       : the max. num. of cluster in the system.       cccccc
ccccc    LB(i)    : is the charge, i.e., =1(proton),=0(neutron).  cccccc
ccccc    XJR(I,j) : is the number of proton(j=1) or neutron(j=2)  cccccc
ccccc                    or hyperon(j=3) for I-th fragments.      cccccc
ccccc    ICL(i)   : is the mass number of cluster                 cccccc              
c///////////////////////////////////////////////////////////////////////
ccccc               Z. Q. Feng in Aug 2008 - Oct 2014             cccccc
ccccc    Fragment and hyperfragment production in                 cccccc
ccccc   HICs and proton (antiproton and meson)-induced reactions  cccccc
c///////////////////////////////////////////////////////////////////////

! Restored: -------------------------------------------------------------- GEMINI ----------------------------------------------------------------------- CSM
      USE gem_kind
      USE gemsub
      USE functs
      USE angle
      USE evaporation
      USE evap
      USE light_stuff
      USE tl_stuff
      USE yrast_mod
      USE mass_stuff
      USE mass_selection, ONLY: read_mass_selection,
     &  overlay_selected_masses
      USE stat_mod
      USE spectra
! Restored: ------------------------------------------------------------ GEMINI --------------------------------------------------------------------- CSM



******************************************************************************
CCC*** Restored: **************** 0 *******************CCC CSM
******************************************************************************

******************************************************************************
****** Restored: **************** 0.1 ********************** CSM
      INTEGER, PARAMETER :: NFILES = 10, OUTFILE = 1000, 
     &      INFILE = 100000, ENFILE = 10000

! DDCSS_LNKED
      REAL, PARAMETER :: RMASS = 0.938918, PMASS = 0.138039
      REAL, PARAMETER :: dkin = 1., dthe = 0.1, mkin = 1000.
      INTEGER, PARAMETER :: crossmax = 30, ANGMAX = 100, MG = 1000
c     RMASS, PMASS: the average mass of neutron and proton, the average mass of pions
c      DDCSS-parameters in calculating double differential cross sections in spallation reactions dsigma/dEdtheta (mb/sr*MeV)
c      dkin: size of kinetic energy bin(MeV), dthe: size of polar angle bin, mkin: upper bound of analyzed kinetic energy
c      crossmax: the number of analyzed isotopes
c      ANGMAX: the number of analyzed polar angles
c     The finer the bin, the greater demand on the quality of statistics

! LNVICS
      INTEGER, PARAMETER :: ENDCLK = 500

! TEAE
      REAL, PARAMETER :: delt = 0.5, TMAX = 150
      INTEGER, PARAMETER :: exctop = 1000      ! MeV
      INTEGER, PARAMETER :: TEAMAX = 10
c      delt the time interval Restored comment. CSM
c      TMAX the maximum number of time points Restored comment. CSM
c      TEAMAX the maximum number of elements Restored comment. CSM

! IPCS(Isotope Prduction Cross Sections)
      INTEGER, parameter :: IPCSMAX=100
c      IPCSMAX the maximum number of elements Restored comment. CSM

******************************************************************************
****** Restored: **************** 0.2 ********************** CSM
      EXTERNAL ISOTAB
      COMMON /XYZ/X(MG), Y(MG), Z(MG), xpion(MG), ypion(MG), zpion(MG)
      COMMON /MOM/XK(MG), YK(MG), ZK(MG), LB(MG), isy(MG), xkpion(MG),
     &        ykpion(MG), zkpion(MG)
      COMMON XFR(MG), YFR(MG), ZFR(MG), XKFR(MG), YKFR(MG), ZKFR(MG)
      COMMON  fr(MG), tkin(MG), utwo(MG), uthr(MG), ucou(MG), umod(MG),
     &    usym(MG), ussr(MG), ugrd(MG), frp(MG), couex(MG),utau(MG),
     &    aene(MG),excen(MG),bdgd(MG)
      COMMON /RAMISO/Iso,Iram
      COMMON /NMASS/npr,nta
      COMMON /BINDEN/XBIN(0:MG,0:MG)
      CHARACTER*2 ELS
      INTEGER IC(MG)
      DIMENSION XJR(MG,3)
      DIMENSION IFLAG(MG),frxk(MG),fryk(MG),frzk(MG),fthe(MG)
      DIMENSION AVEMn(0:MG),zMULT(0:MG),aEkinZ(0:MG),atEkinZ(0:MG),
     &    aETrZ(0:MG),atETrZ(0:MG),znMULT(0:MG,0:MG),aEkinA(MG),
     &    aETrA(MG)
      DIMENSION XKinSP(50000),ICKinSP(50000),ih(100000)
      DIMENSION aveMz(0:MG),azMULT(0:MG)
      DIMENSION MCL(MG),MLLL(MG),MLL1(MG),MCLZ(MG),MLLZ(MG),MLZ1(MG)  
      CHARACTER*30 OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,rinf(50)
      REAL*8 animf,animf2,dens(MG)
      DIMENSION  Acc(0:MG), Acimf(0:MG), Alcc(0:MG), Alcimf(0:MG),
     &    Aneu(0:MG),Azbound(0:MG),Azbimf(0:MG),Aneuimf(0:MG),
     &    Aimf(0:MG)
      DIMENSION ndel(MG),ndel1(MG),ndel2(MG),ndel3(MG),ndel4(MG),
     &    nstr(MG),nstr1(MG),nstr2(MG),npion(MG),npion1(MG),npion2(MG),
     &    npion3(MG),rtime(MG),idr(MG),idms(MG),empion(MG),icp(MG),
     &    nsqp(MG),nsqp1(MG),nsqp2(MG),nsqp3(MG),anmp(MG),anmn(MG),
     &    nkzr(MG),nkpl(MG),mszr(MG),mspl(MG),mslamd(MG),msigms(MG),
     &    msigpl(MG),nazr(MG),namn(MG),msakzr(MG),msakmn(MG),em(MG),
     &    idat(MG),nlam(2,MG),nsig(2,MG),nsms(2,MG),nszr(2,MG),
     &    nspl(2,MG),nxi(2,MG),nxne(2,MG),nxch(2,MG),nsqet(MG),neta(MG)
      REAL RAPIPI(-1:1,-80:80), PTRSPI(-1:1,0:80), PTRETA(0:80),
     &    PTPION(-1:1,0:80),PTKAON(-1:1,0:80),drflow(-1:1,-80:80),
     &    epflow(-1:1,-80:80),drfltr(-1:1,0:80),epfltr(-1:1,0:80),
     &    direrr(-1:1,-80:80),ellerr(-1:1,-80:80),dispol(-1:1,-80:80),
     &    PTSQPN(-1:1,0:80), npnrap(-1:1,-80:80), NTRPION(-1:1,0:80),
     &    npntrn(-1:1,0:80), xsecpi(-1:1,0:80), PTSIG(-1:1,0:80), 
     &    RAPETA(-80:80)
      REAL nkarap(-1:1,-80:80), nkatrn(-1:1,0:80), RAPIKA(-1:1,-80:80),
     &    PTRSKA(-1:1,0:80), drkafl(-1:1,-80:80), epkafl(-1:1,-80:80),
     &    drkaer(-1:1,-80:80), elkaer(-1:1,-80:80), drkatr(-1:1,0:80),
     &    epkatr(-1:1,0:80), PTSQKA(-1:1,0:80), NTRKAON(-1:1,0:80),
     &    xseckn(-1:1,0:80), RAPATK(-1:1,-80:80), PTRATK(-1:1,0:80),
     &    xsatkn(-1:1,0:80), NTRSIG(-1:1,0:80), PTSQSG(-1:1,0:80)
      REAL RAPILM(-1:1,-80:80,0:1), PTRSLM(-1:1,0:80,0:1), 
     &    xseclm(-1:1,0:80,0:1),
     &    nlmrap(-1:1,-80:80),drlmfl(-1:1,-80:80),drlmer(-1:1,-80:80),
     &    RAPISG(-1:1,-80:80,0:1), PTRSSG(-1:1,0:80,0:1),
     &    xsecsg(-1:1,0:80,0:1),
     &    nsgrap(-1:1,-80:80),drsgfl(-1:1,-80:80),drsger(-1:1,-80:80)
      REAL RPPRLG(0:20,-80:80,0:10), PTRSCL(0:20,0:80,0:10),
     &    xseccl(0:20,0:80,0:10), RAPFRA(-80:80,0:10),
     &    xinclcs(0:80,0:10), XKINSPC(0:80,0:10), RPPRTR(0:20,-80:80), 
     &    drflcl(0:20,-80:80), epflcl(0:20,-80:80), drtrcl(0:20,0:80),
     &    eptrcl(0:20,0:80), drcler(0:20,-80:80), epcler(0:20,-80:80),
     &    PTSQNL(0:20,0:80), PTNUCL(0:20,0:80), nakrap(-1:1,-80:80),
     &    drpter(0:20,0:80), eppter(0:20,0:80), nclrap(0:20,-80:80),
     &    ncltrn(0:20,0:80), NTRNUCL(0:20,0:80),dratk(-1:1,-80:80),
     &    drerak(-1:1,-80:80),ENERGY(0:100,0:800),
     &    velocity(0:400,0:400),Aisot(0:400,0:400),
     &    velpar(0:400,-400:400)
      DIMENSION CLX(MG),CLY(MG),CLZ(MG),CLKX(MG),CLKY(MG),CLKZ(MG),
     &    CLKIN(MG),ICL(MG),apr(MG),anu(MG),rnp(MG),aprsq(MG), 
     &    anusq(MG),ICHY(MG),ISUC(MG),NNHY(0:1,0:20,MG),NCLA(0:MG),
     &    NCLZ(-1:MG),NLMA(0:MG), RATZ(MG), RTER(MG),
     &    NSGA(0:MG),NSGZ(-1:MG),NXIA(0:MG),NXIZ(-1:MG),NLMLMA(0:MG),
     &    NLMLMZ(-1:MG),NLMSGA(0:MG),NLMSGZ(-1:MG),NSGSGA(0:MG),
     &    NSGSGZ(-2:MG),NATLA(0:MG),NATLZ(-1:MG),NLMZ(-1:MG)
      REAL (KIND=r4)  :: weight
      TYPE (identity) :: cn
      INTEGER COUNTING
      INTEGER ledger_expected_a,ledger_expected_z
      INTEGER ledger_hot_a,ledger_hot_z,ledger_cold_a,ledger_cold_z
      INTEGER ledger_source_a,ledger_source_z

      INTEGER DDCSS_LNKED, COLM, GEMINI,FSVT, CNPCS, TEAE,DECIZ
      INTEGER LNVICS,TESTEX,PRINTDEC, Z11, Z22
      INTEGER gemini_in_out
      LOGICAL FIRST_WRITE_GEMINI
      INTEGER ECHECK
      LOGICAL FIRST_ECHECK
      REAL PPR1, EPR1, VPR1
      REAL SS, S1, S2, SQRT_S
      INTEGER A0, Z0, NF, MAX_FRAG_ID, MAX_FRAG_A
      REAL (KIND=r4) E0, J0

! DDCSS_LNKED
      INTEGER mnkin, mnthe, nkin, nthe, thetad, NSPECIES, NANGLE
      REAL thet, costhet, ekin, theta, diffcrsec, DECNEXC1, DECNEXC2  ! CHENG HG
      REAL cross_count(crossmax,ceiling(mkin/dkin),ceiling(pi/dthe))
      INTEGER cross_count_s(crossmax,0:3,ceiling(10/dkin),
     &    ceiling(pi/dthe))
      INTEGER      DDCSS_LNKED_SETINF(crossmax,2)
      CHARACTER*10      DDCSS_LNKED_LABEL(crossmax)
      REAL DDCSS_LNKED_ANGLE(ANGMAX),DDCSS_LNKED_SPACFAC(ANGMAX)
      REAL KINORD, DIFFUSION(MG)      
      
      REAL DIFF_T
      REAL A_DIFF, B_DIFF, C_DIFF
      ! Output the Energy and total Momentum of Each Cold Fragment
      INTEGER LNKEB(1:CROSSMAX,ceiling(mkin/dkin))
      INTEGER NEKIN   
      REAL lncross

! FSVT
      INTEGER RDMKECM, RDPT, RDYDRM

      INTEGER NKEPION(-1:1,1:80), NYPPION(-1:1,1:80)
      REAL DNKEPION(-1:1,1:80)

      INTEGER NNKEPION(-1:1,1:80), NNTRPION(-1:1,1:80),
     &    NNYPPION(-1:1,1:80)
      INTEGER SNNKEPION(-1:1,1:80), SNNTRPION(-1:1,1:80), 
     &    SNNYPPION(-1:1,1:80)
      REAL SQNOM, SQDEN, MKECM, YDRM, BPION, RPION,X00CUT

! CNPCS(complete nuclide production cross sections)
      INTEGER ele_count(CROSSMAX), CSSPECIES
      INTEGER NZISOTC(0:300), NNISOTC(0:300), NAISOTC(0:300), 
     &    NZISOTH(0:300), NNISOTH(0:300), NAISOTH(0:300)
      INTEGER CNPCS_SETINF(crossmax,2)
      CHARACTER*10 CNPCS_LABEL(crossmax)
      
! TEAE(Time Evolution of Average Excitation Energy of the Heaviest Remnant)
      INTEGER MAXXACL, IIHR, NPROHR, NNEUHR
      REAL SUMEEHR(3,ceiling(TMAX/delT),4), BDEXHR, NEXCHR,EXPNT
      REAL COMPEXC(ceiling(TMAX/delT),4), BDGDHR,MAXEXC
      INTEGER      NRANGE(3), exccross1(0:ceiling(exctop/10.))
      real exccross2(300)
      REAL XT(ceiling(TMAX/delT),MG),YT(ceiling(TMAX/delT),MG),
     &    ZT(ceiling(TMAX/delT),MG)

! IPCS(Isotope Prduction Cross Sections)
      INTEGER IPCS, IPCS_COUNT(IPCSMAX, 0:200), IPCS_Z(IPCSMAX)
      CHARACTER*10      IPCS_LABEL(IPCSMAX)

! LNVICS

! ELEKED
      INTEGER      ELEKED,ELESPECIES,ELEKEDCOUNT(1:CROSSMAX,0:100)
      CHARACTER*10 ELEKED_LABEL(1:CROSSMAX)
      INTEGER      ELEKED_SETINF(1:crossmax,2)

! TEELECS
      INTEGER      TEELECS, TEELECSCOUNT(CEILING(TMAX/DELT),TEAMAX,2) ! Restored: ( , ,1 /2 ) CSM
      INTEGER TEELESPECIES
      CHARACTER*10      TEELECS_LABEL(TEAMAX)
      INTEGER                  TEELECS_SETINF(TEAMAX,2),TECOLD
      
! ANGM_DEC
      INTEGER      ANGM_DEC, TESTANG
      INTEGER      NUMCL, ANGM_LIGHT
      REAL      TOTANGMX, TOTANGMY, TOTANGMZ
      REAL      TOTANGMX1, TOTANGMY1, TOTANGMZ1
      REAL      TOTANGMX2, TOTANGMY2, TOTANGMZ2
      REAL      ANGMCOEF
      REAL      ANGM_CLX(MG),ANGM_CLY(MG),ANGM_CLZ(MG)
      REAL      TOTANGM3

! SURF_COAL
      INTEGER   SURF_COAL, SURF_Z(20), SURF_A(20), CLUST_TYP, VANISH(MG)
      REAL      SURFE(20), SURFMOM(20, 3), SURFAGM(20, 3)
      INTEGER      EMD_NUCL(MG), EMD_CL(MG), IF_CASCADE(MG), IFSURF(MG)

! CONT_ANA
      INTEGER      NTH_RUN, NTH_NOH, NRMIN, NRMAX, NNEVENT
      INTEGER      SWITCH(NFILES), SWITCH1(NFILES), READMESON
      INTEGER      NTH_RUN1, NTHR(3), NTHN(4), NSWITCH
      CHARACTER(LEN = 3)      NTHRCHAR
      CHARACTER(LEN = 4)      NTHNCHAR
      CHARACTER(LEN = 400) :: WRITECHAR(NFILES), READCHAR(NFILES)
      CHARACTER(LEN = 10)      ::      SUFFIX(NFILES)
      CHARACTER(LEN = 100) ::      FILENAME(NFILES)
      REAL      RAN_GEN, RAN_GEN1

! IFSTAT
      INTEGER      IFSTAT(MG)
      
      integer :: max_val, i, new_id
      integer, allocatable :: mmap(:)

! 提取出的最大结团物理量
      integer :: max_val_index(1)             ! 接收 maxloc 返回值的临时数组
      integer :: max_cluster_id               ! 质量最大的碎片编号
      integer :: max_cl_mass_num              ! 最大碎片的 质量数
      integer :: max_cl_prot_num              ! 最大碎片的 质子数
      real*8  :: max_cl_spin                  ! 最大碎片的 内禀自旋
      real*8  :: max_cl_ex_energy             ! 最大碎片的 激发能
C018  integer :: task018_event,task018_pa_a,task018_pa_z
C018  integer :: task018_extra,task018_extra_id,task018_old_id
C018  integer :: task018_sum_a,task018_sum_z,task018_slot
C018  real :: task018_pa_ex
      
      real*8 :: xbcor
      
******************************************************************************
******************************************************************************
******************************************************************************

      

******************************************************************************
CCC*** Restored: **************** 1 *******************CCC CSM
******************************************************************************


******************************************************************************
****** Restored: **************** 1.1 ********************** CSM
      open(2003, file='Input file/anahyp.in', STATUS='old')      ! Restored comment. CSM
      open(1412, file='Input file/massfrdm.dat', STATUS='old')   ! Restored comment. CSM
      read(2003,*)      ! CONT_ANA
      read(2003,*)      NTH_RUN
      read(2003,*)      NTH_NOH
      read(2003,*)      NRMIN, NRMAX
      read(2003,*)      ! Simulation
      read(2003,*)      npar
      read(2003,*)      X00CUT
      read(2003,*)      kktot, kkmin, kkmax
      read(2003,*)      rr0, pp0
      read(2003,*)      READMESON
      read(2003,*)      ! Models
      read(2003,*)      COLM, GEMINI
      read(2003,*)      ANGM_DEC, ANGM_LIGHT, TESTANG
      read(2003,*)      SURF_COAL
      read(2003,*)      ! Coalescence
      read(2003,*)      MCFR
      read(2003,*)      TESTEX, OUTENE
      read(2003,*)      ! Gemini
      read(2003,*)      DECNEXC1, DECNEXC2, DECIZ      
      read(2003,*)      PRINTCLKIN
      read(2003,*)      PRINTDEC, Z11, Z22
      read(2003,*)      gemini_in_out
      read(2003,*)      ! EnergyCheck
      read(2003,*)      ECHECK
      read(2003,*)      ! TEAE
      read(2003,*)      TEAE
      read(2003,*)      EXCCROSS, TEAMUL
      read(2003,*)      EX1DOWN, EX1UP
      read(2003,*)      EX2DOWN, EX2UP
      read(2003,*)      EX3DOWN
      read(2003,*)      EXPNT
      read(2003,*)      PRTEX_WRTPHAS
      read(2003,*)      ! Yields
      read(2003,*)      CNPCS, IPCS, TEELECS, TECOLD
      read(2003,*)      ! Kinematics
      read(2003,*)  LNVICS, ELEKED, DDCSS_LNKED, DDCSS_LNKED_S, DIFF_T
      read(2003,*)      B_DIFF
      read(2003,*)      ! Tasks
      read(2003,*)      FSVT
      read(2003,*)      ! Dynamics
      read(2003,*)      alpha,batat,gmmat
      read(2003,*)      gsur,gtau11,eta11
      read(2003,*)      gsym,symm
      read(2003,*)      momp
      read(2003,*)      ! Gemini settings
      call INPUT_ANA()
      CALL read_mass_selection
      print*, 'NRMIN = ', NRMIN, 'NRMAX = ', NRMAX
      
      theta0 = 0.
      dtheta = 0.
      nmint = 4

      XBIN(0:MG,0:MG)=0.
      do I=1,10000
          read(1412,*,end=128,err=128) IZ, IA, ELS, I1, Be
          XBIN(IZ,IA) = Be
      enddo
 128  continue
      CALL overlay_selected_masses(XBIN,MG,MG)

******************************************************************************
****** Restored: **************** 1.2 ********************** CSM
! Restored comment. CSM
      SWITCH1(1:NFILES) = 0
      
      SWITCH1(1)  = ANGM_DEC
      SUFFIX(1)      = '.dat'
      FILENAME(1) = 'initangm'
      READCHAR(1) = 'results'//'/'//TRIM(FILENAME(1))//SUFFIX(1)
      
      SWITCH1(2)  = 1
      SUFFIX(2)      = '.dat'
      FILENAME(2) = 'coll'
      READCHAR(2) = 'results'//'/'//TRIM(FILENAME(2))//SUFFIX(2)
      
      SWITCH1(3)  = SURF_COAL
      SUFFIX(3)      = '.dat'
      FILENAME(3) = 'SURF_COAL'
      READCHAR(3) = 'results'//'/'//TRIM(FILENAME(3))//SUFFIX(3)

      SWITCH1(4)  = COUNTSUR
      SUFFIX(4)      = '.dat'
      FILENAME(4) = 'COUNTSUR'
      READCHAR(4) = 'results'//'/'//TRIM(FILENAME(4))//SUFFIX(4)

      if(NTH_RUN == 0)then
            do I = 1, NFILES
                  READCHAR(I) = TRIM(FILENAME(I))//SUFFIX(I)
        enddo
        READCHAR(1) = TRIM(READCHAR(1))
        READCHAR(2) = TRIM(READCHAR(2))
        READCHAR(3) = TRIM(READCHAR(3))
        READCHAR(4) = TRIM(READCHAR(4))
      endif

      do I = 1, NFILES
            if(SWITCH1(I) == 0)cycle
            open(INFILE+I, FILE=READCHAR(I), FORM='FORMATTED',
     &      STATUS='unknown')
      enddo

      if (DDCSS_LNKED == 1) open(556,FILE='Input file/DDCSS_LNKED_in.in'
     &  ,STATUS='unknown')
      if (ELEKED == 1) open(18,FILE='Input file/ELEKED_in.in',
     &  STATUS='old')
      if (TEELECS == 1) open(19, FILE='Input file/TEELECS_in.in',
     &  STATUS='old')
      if (IPCS == 1) open(15, FILE='Input file/IPCS_in.in',STATUS='old')
      
! Restored comment. CSM
      SWITCH(1:NFILES) = 0

      NTH_RUN1 = NTH_RUN
      NTHR(1) = FLOOR(NTH_RUN1/100.) + 48
      NTH_RUN1 = NTH_RUN1 - FLOOR(NTH_RUN1/100.)*100
      NTHR(2) = FLOOR(NTH_RUN1/10.) + 48
      NTH_RUN1 = NTH_RUN1 - FLOOR(NTH_RUN1/10.)*10
      NTHR(3) = FLOOR(NTH_RUN1/1.) + 48
      NTHRCHAR = CHAR(NTHR(1))//CHAR(NTHR(2))//CHAR(NTHR(3))

      NTHN(1) = FLOOR(nth_NOH/1000.) + 48
      NTH_NOH = NTH_NOH - FLOOR(nth_NOH/1000.)*1000
      NTHN(2) = FLOOR(nth_NOH/100.) + 48
      NTH_NOH = NTH_NOH - FLOOR(nth_NOH/100.)*100
      NTHN(3) = FLOOR(nth_NOH/10.) + 48
      NTH_NOH = NTH_NOH - FLOOR(nth_NOH/10.)*10
      NTHN(4) = FLOOR(nth_NOH/1.) + 48
      NTHNCHAR = CHAR(NTHN(1))//CHAR(NTHN(2))//CHAR(NTHN(3))//
     &                  CHAR(NTHN(4))
      
      NSWITCH            = 1
      SWITCH(NSWITCH)    = DDCSS_LNKED
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'DDCSS_LNKED_out'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)
      
      NSWITCH            = 2
      SWITCH(NSWITCH)    = ELEKED
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'ELEKED_out'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)

      NSWITCH            = 3
      SWITCH(NSWITCH)    = CNPCS
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'CNPCS_HOT'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)

      NSWITCH            = 4
      SWITCH(NSWITCH)    = CNPCS
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'CNPCS_COLD'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)

      NSWITCH            = 5
      SWITCH(NSWITCH)    = NTH_RUN
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'MONITOR'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)
      
      NSWITCH            = 6
      SWITCH(NSWITCH)    = TEELECS
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'TEELECS_out'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)

      NSWITCH            = 7
      SWITCH(NSWITCH)    = EXCCROSS
      SUFFIX(NSWITCH)    = '.dat'
      FILENAME(NSWITCH)  = 'EXCCROSS'
      WRITECHAR(NSWITCH) = 'ana'//NTHRCHAR//'/'//TRIM(FILENAME(NSWITCH))
     &                        //NTHRCHAR//'_'//NTHNCHAR//SUFFIX(NSWITCH)
      
      if(NTH_RUN == 0)then
            do      I = 1, NFILES
                  WRITECHAR(I) = TRIM(FILENAME(I))//TRIM(SUFFIX(I))
            enddo
      endif

      do I = 1, NFILES
            if(SWITCH(I) == 0)cycle
            open(OUTFILE+I, FILE =trim(WRITECHAR(I)), FORM='FORMATTED',
     &      STATUS='unknown')
      enddo

      if(COLM == 1)then
        open(189,FILE='MST_ACCOUNTING.dat',STATUS='REPLACE')
        write(189,*) '# event IT expected_A expected_Z hot_A hot_Z',
     &               ' cold_A cold_Z'
      endif
      if(PRTEX_WRTPHAS == 1)open(13,FILE='WRTPHAS.dat',STATUS='unknown')
      if(FSVT. EQ. 1)open(3, FILE='FSVT.dat', STATUS='unknown')
      if(TEAMUL == 1)open(6, FILE='TEAMUL.dat',STATUS='unknown')
      if(gemini_in_out == 1)then
C018 TASK-018: uncomment C018 lines with blank or & in column 6,
C018 and comment the ordinary call input below. gemini_in_out=1.
C018 TASK018_PA_ENERGY.dat has one row per event, in coll.dat order:
C018 event_id A Z recalculated_E_MeV, using the same event phase space.
C018 The selected-source logs are in the working directory.
C018  open(185,FILE='TASK018_SELECTED.dat',STATUS='REPLACE')
C018  open(186,FILE='TASK018_FRAGMENTS.dat',STATUS='REPLACE')
C018  open(187,FILE='TASK018_PA_ENERGY.dat',STATUS='OLD')
        open(182,FILE='GEMINI_INPUT.dat',STATUS='unknown')
        open(184,FILE='FRAG_CORR.dat',STATUS='unknown')
        write(184,*) '# Z0 A0 J0 E0 Zf Af Ef Jf Id WGT'
        write(184,*) '# J0 is shared time-slice angular momentum.'
        write(184,*) '# E0 is clamped nexc, not raw excen.'
        write(184,*) '# Source is MST maxloc(ICL) for each time slice.'
        write(184,*) '# Final leaf has largest valid A; ties keep first 
     &    index.'
        write(184,*) '# With kktot=1, data-row order follows event L1.'
        FIRST_WRITE_GEMINI = .TRUE.
      endif
      if(ECHECK == 1)then
        open(183, FILE='E_SH.out', STATUS='unknown')
        FIRST_ECHECK = .TRUE.
      endif
      if(LNVICS == 1)open(10, FILE='LNVICS.dat',STATUS='unknown')
      if(IPCS == 1)open(9,FILE='IPCS_out.dat',STATUS='unknown')
      
C       call RANDOM_SEED()
C       do I = 1, NINT(NRMAX*RAN_GEN)
C             call random_number (RAN_GEN1) 
C       enddo
      
******************************************************************************
****** Restored: **************** 1.3 ********************** CSM

****** Restored: **************** 1.3.1 ********************** CSM
      do I = 1, MG
            Acc(i)=0.
            Acimf(i)=0.
            Alcc(i)=0.
            Alcimf(i)=0.
            Aneu(i)=0.
            Aneuimf(i)=0.
            Azbound(i)=0.
            Azbimf(i)=0.
            Aimf(i)=0.
            NCLA(I)=0
            NLMA(I)=0       ! Hypernucleus
            RATZ(I)=0.
            RTER(I)=0.
            NSGA(I)=0
            NXIA(I)=0
            NATLA(I)=0
            NLMLMA(I)=0
            NLMSGA(I)=0
            NSGSGA(I)=0
      enddo
      
      NCLZ(-1:MG)=0
      NLMZ(-1:MG)=0
      NSGZ(-1:MG)=0
      NXIZ(-1:MG)=0
      NATLZ(-1:MG)=0
      NLMLMZ(-1:MG)=0
      NLMSGZ(-1:MG)=0
      NSGSGZ(-2:MG)=0
      
      do J=-80,80
          RAPETA(J)=0.
          do K = -1, 1
            RAPIPI(K,J)=0.
            drflow(K,J)=0.
            epflow(K,J)=0.
            direrr(K,J)=0.
            ellerr(K,J)=0.
            dispol(K,J)=0.
            npnrap(K,J)=0.
            nkarap(K,J)=0.
            RAPIKA(K,J)=0.
            RAPATK(K,J)=0.
            drkafl(K,J)=0.
            epkafl(K,J)=0.
            drkaer(K,J)=0.
            elkaer(K,J)=0.
            RAPILM(K,J,0:1)=0.
            RAPISG(K,J,0:1)=0.
            nlmrap(K,J)=0.
            nsgrap(K,J)=0.
            drlmfl(K,J)=0.
            drlmer(K,J)=0.
            drsgfl(K,J)=0.
            drsger(K,J)=0.
            dratk(K,J) =0. 
          drerak(K,J)=0.
            nakrap(K,J)=0
          if(J >= 0)then
                  PTRSPI(K,J)=0.
                  PTPION(K,J)=0.
                  PTRETA(J)=0.
                  PTKAON(K,J)=0.
                  drfltr(K,J)=0.
                  epfltr(K,J)=0.
                  PTSQPN(K,J)=0.
                  PTSQKA(K,J)=0.
                  npntrn(K,J)=0.
                  nkatrn(K,J)=0.
                  PTRSKA(K,J)=0.
                  PTRATK(K,J)=0.
                  drkatr(K,J)=0.
                  epkatr(K,J)=0.
                  PTRSLM(K,J,0:1)=0.
                  PTRSSG(K,J,0:1)=0.
                  xsecpi(K,J)=0.
                  xseckn(K,J)=0.
                  xseclm(K,J,0:1)=0.
                  xsecsg(K,J,0:1)=0.
                  xsatkn(K,J)=0.
          endif
          enddo
            do K = 0, 20
                  RPPRLG(K,J,0:10)=0.
                  RPPRTR(K,J)=0.
                  drflcl(K,J)=0.
                  epflcl(K,J)=0.
                  drcler(K,J)=0.
                  epcler(K,J)=0.
                  nclrap(K,J)=0.
                  if(J >= 0)then
                        drtrcl(K,J)=0.
                        eptrcl(K,J)=0.
                        PTSQNL(K,J)=0.
                        PTNUCL(K,J)=0.
                        drpter(K,J)=0.
                        eppter(K,J)=0.
                        ncltrn(K,J)=0.
                        PTRSCL(K,J,0:10)=0.
                        xseccl(K,J,0:10)=0.
                  endif
            enddo
          RAPFRA(J,0:10)=0.
          if(J >= 0) xinclcs(J,0:10)=0.
          if(J >= 0) XKINSPC(J,0:10)=0.
      enddo

      do J = 0, 800
            do k = 0, 100
                ENERGY(k,J)=0.
            enddo
      enddo

      do J = 0, 400
          do k = 0, 400
              velocity(k,J)=0.
              Aisot(k,J)=0.
            enddo
      enddo

      do J = -400, 400
            do k = 0, 400
              velpar(k,J)=0.            
            enddo
      enddo

****** Restored: **************** 1.3.2 ********************** CSM
      do I = 1, MG
            fthe(I)=0.
            !ANimf(I)=0.
            !ANimf2(I)=0.
            aEkinA(I)=0. 
            aETrA(I)=0.
            rtime(I)=0.
            ndel(I)=0
            ndel1(I)=0
            ndel2(I)=0
            ndel3(I)=0
            ndel4(I)=0
            nstr(I)=0
            nstr1(I)=0
            nstr2(I)=0
            npion(I)=0
            npion1(I)=0
            npion2(I)=0
            npion3(I)=0
            neta(I)=0
            nsqp(I)=0
            nsqp1(I)=0
            nsqp2(I)=0
            nsqp3(I)=0
            nsqet(I)=0
            apr(I)=0
            anu(I)=0
            rnp(I)=0
            aprsq(I)=0.
            anusq(I)=0.
            dens(I)=0.D0
            anmp(I)=0.
            anmn(I)=0.
            nkzr(I)=0
            nkpl(I)=0
            nazr(I)=0
            namn(I)=0
            mszr(I)=0
            mspl(I)=0
            msakzr(I)=0
            msakmn(I)=0
            nlam(1:2,I)=0
            nsig(1:2,I)=0
            nsms(1:2,I)=0
            nszr(1:2,I)=0
            nspl(1:2,I)=0
            nxi(1:2,I) =0
            nxne(1:2,I)=0
            nxch(1:2,I)=0
            msigms(I)=0
            msigpl(I)=0
            mslamd(I)=0
      enddo

      do I = 1, 50000
            XKinSP(I)=0.
            ICKinSP(I)=0
      enddo

      do I = 0, MG
            aveMn(I)=0.
            zmult(I)=0.
            aEkinZ(I)=0. 
            atEkinZ(I)=0. 
            aETrZ(I)=0. 
            atETrZ(I)=0. 
          do J = 0, MG
                znMULT(I,J)=0.
          enddo 
      enddo

! Restored comment. CSM
      ANIMF = 0.0d0
      ANIMF2 = 0.0d0
      DM = 0.02      ! trans. mom./kinetic energy   GeV
      DY = 0.1       ! rap. dis.
      DP = 0.5       ! polar dis.
      
      xc = 0.
      yc = 0.
      zc = 0.

! CNPCS(CNP Cross Section)
      if(CNPCS == 1)then
        CSSPECIES = 0
        open(181,FILE='Input file/CNPCS_in.in'
     &  ,STATUS='old')
        read(181,*)
        read(181,*)
        do J = 1, CROSSMAX
          read(181,*,END=278,err=278) CNPCS_LABEL(J),
     &    CNPCS_SETINF(J,1), CNPCS_SETINF(J,2)
          CSSPECIES = CSSPECIES + 1
        enddo        
278     continue
        close(181)
        ele_count = 0
        do IIP = 0, 300
          NZISOTC(IIP) = 0
          NNISOTC(IIP) = 0
          NAISOTC(IIP) = 0

          NZISOTH(IIP) = 0
          NNISOTH(IIP) = 0
          NAISOTH(IIP) = 0
        enddo
      endif

! DDCSS_LNKED(Double Diffrencial Cross Sections && Light Nuclide Kinetic Energy Distribution)
      if(DDCSS_LNKED == 1)then
            NSPECIES = 0
            read(556,*)
            read(556,*)
            do I = 1, crossmax
                  read(556,*,END=2366,err=2366) DDCSS_LNKED_LABEL(I),
     &           DDCSS_LNKED_SETINF(I,1),DDCSS_LNKED_SETINF(I,2)
                  NSPECIES=NSPECIES+1
            enddo
            
2366        read(556,*)
            NANGLE=0
            do I = 1, ANGMAX
                  read(556,*,END=2367,err=2367) DDCSS_LNKED_ANGLE(I),
     &            DDCSS_LNKED_SPACFAC(I)
                  NANGLE = NANGLE + 1
            enddo
2367            continue
            close(556)
                        
            do k = 1, crossmax                  
                  !DDCSS
                  do i = 1, ceiling(mkin/dkin)           ! CHENG HG
                        do j = 1, ceiling(pi/dthe)
                            cross_count(k,i,j) = 0
                        enddo
              enddo
              
                  do i = 1, ceiling(10/dkin)            ! CHENG HG
                        do j = 1, ceiling(pi/dthe)
                            CROSS_COUNT_S(K,0:3,I,J) = 0
                        enddo
              enddo      
              
                  !LNKED
                  do NEKINI = 1, ceiling(mkin/dkin)
                        LNKEB(K,NEKINI) = 0
                  enddo
            enddo
      endif

! ELEKED(Element Kinetic Energy Distribution)
      if(ELEKED == 1)      then

            ELESPECIES=0
            read(18,*)
            read(18,*)
            do      J=1,CROSSMAX
                  read(18,*,END=277,ERR=277)      ELEKED_LABEL(J),
     &            ELEKED_SETINF(J,1), ELEKED_SETINF(J,2)

                  ELESPECIES=ELESPECIES+1
            enddo

277            continue

            do      I = 1,ELESPECIES
                  do      J=0,100
                        ELEKEDCOUNT(I,J)=0
                  enddo
            enddo

      endif      

! TEELECS(Time Evolution of Element Cross Section)
      if(TEELECS == 1)then
            TEELESPECIES = 0 ! Restored comment. CSM
            read(19,*)
            read(19,*)
            do J = 1, TEAMAX
                  read(19,*,END=228,ERR=228)      TEELECS_LABEL(J),
     &            TEELECS_SETINF(J,1), TEELECS_SETINF(J,2)
                  TEELESPECIES = TEELESPECIES + 1
            enddo
228            continue

            do I = 1, TEELESPECIES
                  do J = 1, CEILING(TMAX/DELT)
                        TEELECSCOUNT(J,I,1:2) = 0
                  enddo
          enddo
      endif
      
! Restored: FSVT( ) CSM
      if(FSVT == 1)then
            do M = -1, 1
                  do KK = 1, 50
                        NNKEPION(M,KK)=0
                        NNTRPION(M,KK)=0
                        NNYPPION(M,KK)=0

                        SNNKEPION(M,KK)=0
                        SNNTRPION(M,KK)=0
                        SNNYPPION(M,KK)=0

                        DNKEPION(M,KK)=0
                  enddo
            enddo
      endif

! Restored: TEAE( ) CSM
      BDEXHR=0.
      NEXCHR=0.
      BDGDHR=0.
      NPROHR=0
      NNEUHR=0
      if(TEAE == 1)      then
            do T=1,ceiling(TMAX/delT)
                  SUMEEHR(1:3,T,1)=0.
                  SUMEEHR(1:3,T,2)=0.
                  SUMEEHR(1:3,T,3)=0.
                  SUMEEHR(1:3,T,4)=0.
            enddo
      endif

      COMPEXC(1:ceiling(TMAX/delT),1)=0.
      COMPEXC(1:ceiling(TMAX/delT),2)=0.
      COMPEXC(1:ceiling(TMAX/delT),3)=0.
      COMPEXC(1:ceiling(TMAX/delT),4)=0.

      NRANGE(1:3)=0
      EXCCROSS1(0:ceiling(EXCTOP/10.))=0
      EXCCROSS2(1:300) = 0.
      EXCTOT2 = 0
! Restored: IPCS( ) CSM
      if(IPCS == 1)      then
            do      I = 1, IPCSMAX
                  IPCS_COUNT(I,0:200) = 0
                  IPCS_COUNT(I,0:200) = 0
            enddo

            NIPCS=0      
            read(15,*)
            read(15,*)
            do      I = 1, IPCSMAX
                  read(15,*,END=2236,ERR=2236) IPCS_LABEL(I), IPCS_Z(I)
                  NIPCS = NIPCS + 1
            enddo
2236            continue
            close(15)
      endif

! Restored: LNVICS( ) CSM
      if(LNVICS == 1)then
            
      endif

      L11=1
      IFEV=1
      NNEVENT = 0

******************************************************************************
****** Restored: **************** 1.4 ********************** CSM
      mseed = 200
      AMU = 0.938
      AMU1 = 0.938271998
      AMU2 = 0.939565330
      ALM = 1.115
      ASG = 1.19
      AXI = 1.32

      ukin = (0.197327)**2/2./0.9315
      xmn = 8.071
      xmp = 7.289
      ccf = 3.*pi**2./2.   
      cex = -3./4.*0.00144*(3/pi)**0.3333333
      
      cdelt = 0.00176                                     ! 00176
      if (MOMP == 0) CDELT = 0.
      cmeps = 500.
      xsym = -0.65
      cmnn = cdelt*(1. + xsym)
      cmnp = cdelt*(1. - xsym)
      rho0 = 0.16                                         ! Restored comment. CSM
      Ekfac = 0.01245                                     ! 15   ! GeV

      ! use ih array to store the distance (fm)
      do j = 1, 99999                                       
        dr2 = float(j)/40.        ! fm^2, range (0.025, 2499.975)                              
        dr = sqrt(dr2)            ! fm, range (0.158, 49.999)                                  
        ir = nint(dr*100. + 0.5)  ! manual rounding to get the dr*100 (0.01 fm)                           
        ih(j+1) = ir              ! store in the ih array
      end do
      ih(1) = 1                   ! ex. ih(2) = 16, ih(3) = 22, ih(4) = 27 (0.01 fm)

******************************************************************************
****** Restored: **************** 1.5 QMD ********************** CSM
      read(INFILE+2,965) nt, x00min, x00max, epp, dt, al
      read(INFILE+2,966) nta, npr, nzta, nzpr
 965  format(/1x,I4,5(F8.3,1x))
 966  format(1x,4(I3,1x))

      PPR1 = SQRT((NPR*EPP)**2 + 2.*(NPR*0.938)*EPP)
      EPR1 = NPR*EPP + NPR*0.938
      VPR1 = PPR1/(EPR1 + NTA*0.938)

      xb00 = x00min
      xmom = sqrt(epp*(epp+2.*AMU))                       ! Restored comment. CSM
      print*, 'Apr =', npr, 'Zpr =', nzpr, 'Ata =', nta, 'Zta =', nzta
      print*, 'Kin. Ener. in Lab (GeV/u)', epp, 'Mom. Lab.(GeV/c)', xmom

************************************************************************
************************************************************************
************************************************************************


      
************************************************************************
CCC**************** 2 Global Kinematics Initialization **************CCC
************************************************************************

                                                                                
************************************************************************
******************* 2.1 Reference System Initialization ****************

      yprcm = 1.
      pzcmpr = 1.
      if (npr > 0) then
        ! Lab system
        
        ! Target initial momentum p_target = (Eta,0,0,0)
        eta = float(nta) * rmass ! CHEN SM
        pzta = 0.
        betata = 0.
        gammta = 1.
        yta = 0.
        ! Projectile initial momentum p_projectile = (Epr,0,0,pzpr)
        if (npr /= 0) then  ! baryon
          epr = float(npr) * (rmass + epp)                  ! Total energy Etot = Erest + T
          pzpr = sqrt(epr**2 - (rmass * npr)**2)            ! Momentum magnitude pz/p = (Etot^2 - Erest^2)^(1/2)
          betapr = pzpr / epr                               ! Projectile velocity beta/v = pz/Etot
          gammpr = 1. / sqrt(1. - betapr**2)                ! gamma parameter gama = (1 - beta^2)^(-1/2)
          ypr = 0.5 * log((epr + pzpr) / (epr - pzpr))      ! Rapidity y = (1/2)*ln[(E+pz)/(E-pz)] = (1/2)*ln[(1+v)/(1-v)]
        else                ! lapton
          epr = 0.
          pzpr = 0.
          betapr = 0.
          gammpr = 1.
        end if

        betacm = pzpr / (eta + epr)                         ! Center of mass beta parameter betacm = pz/Etot
        gamacm = 1. / sqrt(1. - betacm**2)                  ! Center of mass gamma parameter

        ! Center of mass system
        ! CHENG HG
        !ECM = epp*float(npr*nta)/float(npr+nta)             ! Sum of kinetic energies E' (C.M. kinetic energy)
        !EPRCM = ECM + ETA                                   ! proj. energy Etot = Erest + T
        !PZCM = SQRT(EPRCM**2 - ETA**2)                      ! C.M. momentum pz = (Etot^2 - Erest^2)^(1/2)
        !ETACM = SQRT(PZCM**2 + (FLOAT(nta) * RMASS)**2)     ! targ. energy
        !YPRCM = 0.5*ALOG((EPRCM + PZCM)/(EPRCM - PZCM))     ! proj. rapidity y1 = (1/2)*ln[(E+pz)/(E-pz)] = (1/2)*ln[(1+v)/(1-v)]
        !YTACM = 0.5*ALOG((ETACM - PZCM)/(ETACM + PZCM))     ! targ. rapidity y2 = (1/2)*ln[(E+pz)/(E-pz)] = (1/2)*ln[(1+v)/(1-v)]
        !PZCMPR = PZCM/FLOAT(npr)                            ! Momentum per nucleon in C.M. system
      
        ! CHEN SM
        ss = (epr + eta)**2 - pzpr**2                       ! Mandelstam s
        sqrt_s = sqrt(ss)                                   ! C.M. energy sqrt_s = Etot
        ecm = sqrt_s - eta - float(npr) * rmass             ! Sum of total kinetic energies E' (C.M. system)
        s1 = ss + rmass**2 * (float(npr)**2 - float(nta)**2)
        eprcm = s1 / (2 * sqrt_s)                           ! proj. energy E'pr = (s + m1^2 - m2^2)/(2*Etot)
        s2 = ss + rmass**2 * (float(nta)**2 - float(npr)**2)  
        etacm = s2 / (2 * sqrt_s)                           ! targ. energy E'ta = (s - m1^2 + m2^2)/(2*Etot)
        pzcm = sqrt(etacm**2 - eta**2)                      ! C.M. momentum pz = E'ta^2 - m2^2 = E'pr^2 - m1^2
        yprcm = 0.5 * log((eprcm + pzcm) / (eprcm - pzcm))  ! proj. rapidity y1 = (1/2)*ln[(E+pz)/(E-pz)] = (1/2)*ln[(1+v)/(1-v)]
        ytacm = 0.5 * log((etacm - pzcm) / (etacm + pzcm))  ! targ. rapidity y1 = (1/2)*ln[(E+pz)/(E-pz)] = (1/2)*ln[(1+v)/(1-v)]
        pzcmpr = pzcm / npr                                 ! Momentum per nucleon in C.M. system
      
        print *, 'Momentum per Nucleon in C.M.(GeV/c)', pzcmpr
        print *, 'Rapidities of proj. and targ. in C.M.', yprcm, ytacm

        ! Nucleon-Nucleon reference system calculations
        if (nta /= 0) then
          if (npr == 0) then  ! lapton
            el1 = pmass + epp                               ! Total energy in lab system
            pl1 = sqrt(el1*el1 - pmass**2)                  ! Momentum
            ss = (el1 + rmass)**2 - pl1*pl1
            ss1 = ss + pmass**2 - rmass**2
            etcm = ss1 / (2. * sqrt(ss))                    ! Energy per nucleon in C.M. CHEN SM
            !ETCM = 0.5*SQRT(SS)                            ! Energy per nucleon in C.M. (Two-body equal mass approximation)
            plcm = sqrt(etcm**2 - pmass**2)                 ! Momentum magnitude in C.M.
          else                ! baryon
            el1 = rmass + epp
            pl1 = sqrt(el1*el1 - rmass**2)
            ss = (el1 + rmass)**2 - pl1*pl1
            ss1 = ss
            etcm = ss1 / (2. * sqrt(ss))                    
            !ETCM = 0.5*SQRT(SS)                        
            plcm = sqrt(etcm**2 - rmass**2)
          end if
          cta = float(nta) * etcm                           ! targ. energy in C.M. system
          czta = -float(nta) * plcm                         ! targ. momentum z in C.M. system
          cetata = czta / cta                               ! targ. beta parameter
          cyta = 0.5 * log((cta + czta) / (cta - czta))     ! targ. rapidity
        else
          czta = 0.
          cta = 0.
          cetapr = 0.
        end if

        cammta = 1.0 / sqrt(1.0 - cetata**2)                ! targ. gamma parameter
        if (npr /= 0) then
          czpr = float(npr) * plcm                          ! proj. momentum z
          cpr = float(npr) * etcm                           ! proj. energy
          cetapr = czpr / cpr                               ! proj. beta parameter
          cypr = 0.5 * log((cpr + czpr) / (cpr - czpr))     ! proj. rapidity
        else
          czpr = 0.
          cpr = 0.
          cetapr = 0.
        end if
        cammpr = 1.0 / sqrt(1.0 - cetapr**2)                ! proj. gamma parameter
        print *, 'Rapidities of NN in C.M.', cypr, cyta
      end if


      r0rec = 1.16

******************************************************************************
********************** 2.2 Radius Scaling Parameter **************************
2101  r0rec = 1.16
      idatpr = 0 
      if(idatpr > 0) r0rec = 1.38                        ! PRL77(1996)1230
      rrt = r0rec*nta**0.3333333
      x00 = x00min
      dx0 = 2.                                            ! fm
      csf = 10.*(2.*pi*x00*dx0)                           ! prod. cross section in mb
      if(npr <= nmint) csf = 10.*pi*rrt**2
      !csf = 1.                                            ! for HICs

      mz = nzta + nzpr
      mn = nta + npr - mz
      mass = npr + nta

      !rr0 = 3.                                            ! fm
      !pp0 = 0.15                                          ! GeV
      rhy = 5.
      RCR = (rr0)**2.                               
      PCR = (pp0)**2. 
      RSH = (rhy)**2
      
      print*

******************************************************************************
******************************************************************************
******************************************************************************



******************************************************************************
CCC*** Restored: **************** 3 *******************CCC CSM
******************************************************************************


      NEVENT = 0
      L11 = 1
      IFEV=1
      NNEVENT = 0
******************************************************************************
****** Restored: **************** 3.1 L1 ********************** CSM
      do 9666 L1 = 1, npar
      
      if(MOD(L1,10) == 0.and.NTH_RUN == 0)write(*,*) 'Running No.', L1
      NEVENT = NEVENT + 1            ! Restored comment. CSM
      MAXEXC = 0.
      FLAGEXC = 0
      IFFRAG = 0
      
      DIFFUSION(1:MG) = -1.
      IFSTAT(1:MG) = 0

! Restored: SURF_COAL: CSM
      if(SURF_COAL == 1)then
            IFSURF(1:MG) = 0 ! Restored: , 1, 0 CSM
               VANISH(1:MG) = 0 ! Restored: , 1, 0 CSM
        ! Restored comment. CSM
        ! Restored: 1 , 20 CSM
        ! Restored: , SURF_COAL CSM
        ! Restored: gemini, CSM
            SURF_Z(1:20) = 0 
        SURF_A(1:20) = 0
            SURFE(1:20) = 0.
        do      MM = 1, 20
                  SURFMOM(MM, 1:3) = 0.
                  SURFAGM(MM, 1:3) = 0.
            enddo
            do      MM = 1, 20
3123                  continue
                  read(INFILE+3,*,ERR=3123) SURF_Z(MM), SURF_A(MM),
     &        SURFE(MM), SURFMOM(MM, 1), SURFMOM(MM, 2), SURFMOM(MM, 3),
     &            SURFAGM(MM, 1),SURFAGM(MM, 2), SURFAGM(MM, 3)
                  if(SURF_Z(MM) == 0.and.SURF_A(MM) == 0)EXIT
                  if(SURFE(MM) <= 0.)then
                        SURFE(MM) = 543.65
                        SURFMOM(MM, 1) = 0.1395
                        SURFMOM(MM, 2) = 1.4378
                        SURFMOM(MM, 3) = 0.9863
                        SURFAGM(MM, 1) = 1.8045
                        SURFAGM(MM, 2) = -7.8661
                        SURFAGM(MM, 3) = 11.2116
                  endif
            enddo
23458            format(1x,i1, f9.2, 6f9.4)
23457            format(1x,i1,7f9.4)
      endif

            
******************************************************************************
****** Restored: **************** 3.2 it ********************** CSM
      do 9668 it = 1, kktot   
      ! Restored: ktotal- coll , rtime(it) CSM
      ! Restored: it=1 k,kktot k,kktot k CSM
      ! Restored: kktot = nt/iwrite CSM
      ! Restored: , coll.out , , GEMINI t>t_freeze CSM
      ! Restored: , CSM
      ! Restored: csm CSM
      PNTMO=1
      ! Restored comment. CSM
      cxx=0.    ! Restored comment. CSM
      cyy=0.
      czz=0.
      pxx=0.    ! Restored comment. CSM
      pyy=0.
      pzz=0.
      
      nucl=0
      mdel=0    ! Restored: delta delta CSM
      mdel1=0
      mdel2=0
      mdel3=0
      mdel4=0

      mstr=0    ! Restored: star star CSM
      mstr1=0
      mstr2=0

      mpion=0   ! Restored: pi pi CSM
      mpion1=0
      mpion2=0
      mpion3=0

      mknzr=0
      mknpl=0
      makzr=0
      makmn=0
      mlamb=0
      msigm=0
      msmin=0
      mszer=0
      mspls=0

      mcsxi=0
      mxine=0
      mxich=0
      matlamb=0
      matsigm=0
      matsmin=0
      matszer=0
      matspls=0
      matcsxi=0
      matxine=0
      matxich=0

      rhoc=0.

      do N1 = 0, 20
          do N2 = 0, 80
              NTRNUCL(N1,N2) = 0.
          enddo
      enddo
      
! FSVT
      if(FSVT == 1) then
          do N1 = -1, 1
              do N2 = 1, 50
                  NTRPION(N1,N2)=0.
                  NKEPION(N1,N2)=0.
                  NYPPION(N1,N2)=0.
              enddo
          enddo
      endif

! ANGM_DEC
      if(ANGM_DEC == 1)then
            if(IT == 1)then      
                  read(INFILE+1,*) L12, TOTANGMX1,
     &                  TOTANGMY1,TOTANGMZ1
            endif
      endif
      
******************************************************************************
****** Restored: **************** 3.3 mass=nta+npr ********************** CSM

****** Restored: **************** 3.3.1 mass=nta+npr - ***** CSM
      do 2892 i = 1, mass
       
      ! Restored: LB=iso: , 1, 0 CSM
      ! Restored: isy: , -1, 1 CSM
      ! Restored: idr: >0 CSM
      ! Restored: idat: , 0, >0 CSM
      ! Restored: vanish: , 1, 0 CSM
      ! idms: 
      if(idatpr  >  0)then
 668    read(INFILE+2,9669,err=668) k,x(i),y(i),z(i),xk(i),yk(i),
     &    zk(i),em(i),LB(i),isy(i),id,idr(i),idat(i)
 670    read(INFILE+2,9669,err=670) k,xpion(i),ypion(i),zpion(i),
     &    xkpion(i),ykpion(i),zkpion(i),empion(i),idms(i),icp(i),
     &    idat(i)
      else
        if(SURF_COAL == 1)then
      ! Restored: , CSM
 299      read(INFILE+2,9670,err=299) k,x(i),y(i),z(i),xk(i),yk(i),
     &      zk(i),em(i),LB(i),isy(i),id,idr(i), VANISH(I),     
     &      EMD_NUCL(I)                                                          
 291      if(READMESON == 1)read(INFILE+2,9669,err=291)
     &      k,xpion(i),ypion(i),zpion(i),xkpion(i),ykpion(i),
     &      zkpion(i),empion(i),idms(i),icp(i) 
        else
      ! Restored: , CSM
 289      read(INFILE+2,9669,err=289) k,x(i),y(i),z(i),xk(i),yk(i),
     &      zk(i),em(i),LB(i),isy(i),id,idr(i)                  
 290      if(READMESON == 1)read(INFILE+2,9669,err=290)
     &      k,xpion(i),ypion(i),zpion(i),xkpion(i),ykpion(i),
     &      zkpion(i),empion(i),idms(i),icp(i)
        endif
      endif          
9669      format(1x,I4,7f9.4,5I4)
9670  format(1x,I4,7f9.4,6I4)
      
! Restored comment. CSM
      !if(idr(i) <= 3.and.VANISH(I) /= 1)then
      cxx=cxx+x(i)
      cyy=cyy+y(i)
      czz=czz+z(i)
      pxx=pxx+xk(i)
      pyy=pyy+yk(i)
      pzz=pzz+zk(i)
      nucl=nucl+1
      !endif

      if(idat(i) == 0.and.VANISH(I) /= 1)then             ! Restored comment. CSM
      ! Restored comment. CSM
          if(idr(i) == 1)then                             ! Restored: delta CSM
                  mdel=mdel+1
                  if(LB(i) == 2) mdel1=mdel1+1                ! ++
                  if(LB(i) == 1) mdel2=mdel2+1                ! +
                  if(LB(i) == 0) mdel3=mdel3+1                ! 0
                  if(LB(i) == -1)mdel4=mdel4+1                ! -
          elseif(idr(i) == 2)then                         ! star  
              mstr=mstr+1                                 ! Restored: mstr 1440 CSM
              if(LB(i) == 1) mstr1=mstr1+1                ! p(1440)
              if(LB(i) == 0) mstr2=mstr2+1                ! n(1440)
          elseif(idr(i) == 4)then
              mlamb=mlamb+1
          elseif(idr(i) == 5)then
              msigm=msigm+1
              if(LB(i) == -1)msmin=msmin+1
              if(LB(i) == 0) mszer=mszer+1
              if(LB(i) == 1) mspls=mspls+1
          elseif(idr(i) == 6)then
              mcsxi=mcsxi+1
              if(LB(i) == -1)mxich=mxich+1
              if(LB(i) == 0) mxine=mxine+1
          endif
      else
          if(idr(i) == 4)then
              matlamb=matlamb+1
          elseif(idr(i) == 5)then
                  matsigm=matsigm+1
                  if(LB(i) == -1)matsmin=matsmin+1
                  if(LB(i) == 0) matszer=matszer+1
                  if(LB(i) == 1) matspls=matspls+1
          elseif(idr(i) == 6)then
                  matcsxi=matcsxi+1
                  if(LB(i) == 0) matxine=matxine+1
                  if(LB(i) == 1) matxich=matxich+1
          endif
      endif

      ! Restored comment. CSM
      if(idms(i) == 1)then                                ! pion
            mpion=mpion+1
            if(icp(i) == 1) mpion1=mpion1+1                 ! pi+
            if(icp(i) == 0) mpion2=mpion2+1                 ! pi0
            if(icp(i) == -1)mpion3=mpion3+1                 ! pi-
      elseif(idms(i) == 3)then                            ! Restored: kaon CSM
            if(icp(i) == 0) mknzr=mknzr+1                   ! K(0)
            if(icp(i) == 1) mknpl=mknpl+1                   ! K(+)
      elseif(idms(i) == 4)then                            ! Restored: kaon CSM
            if(icp(i) == 0) makzr=makzr+1                   ! k(0*)
            if(icp(i) == -1)makmn=makmn+1                   ! K(-)
      endif

      widxx = 2.*sl(i)                                    ! evolution of central density
      dr2 = (xc-x(i))**2. + (yc-y(i))**2. + (zc-z(i))**2. 
      rhoc = rhoc + XEXP(-dr2/widxx)/(PI*widxx)**1.5      ! Restored: rc CSM
      ! Restored: sl(i) FUNCTION, i, L=sigma^2 , CSM
      ! Restored: i<=npr ,i>npr CSM
2892      continue
      
      if((L1 < NRMIN.OR.L1 > NRMAX).and.(NTH_RUN /= 0))goto 9667 ! Restored: , CSM


      NNEVENT = NNEVENT + 1                               ! Restored comment. CSM
      if(L1 == 1) rtime(it) = k*dt                        ! Restored comment. CSM
      if(kktot == 1) kkmax = kkmin ! CHEN SM

ccccc  9999 ccccc Restored comment. CSM
      ! Restored: [kkmin, kkmax] k CSM
      ! Restored: k , CSM
      ! Restored: , NEVENT, CSM
      if(it >= kkmin.and.it <= kkmax)then
      ndel(it) = ndel(it) + mdel
      ndel1(it) = ndel1(it) + mdel1
      ndel2(it) = ndel2(it) + mdel2
      ndel3(it) = ndel3(it) + mdel3
      ndel4(it) = ndel4(it) + mdel4
      nstr(it) = nstr(it) + mstr
      nstr1(it) = nstr1(it) + mstr1
      nstr2(it) = nstr2(it) + mstr2
      npion(it) = npion(it) + mpion
      npion1(it) = npion1(it) + mpion1
      npion2(it) = npion2(it) + mpion2
      npion3(it) = npion3(it) + mpion3
      neta(it) = neta(it) + meta
      nkzr(it) = nkzr(it) + mknzr
      nkpl(it) = nkpl(it) + mknpl
      nazr(it) = nazr(it) + makzr
      namn(it) = namn(it) + makmn

      nlam(1,it) = nlam(1,it) + mlamb
      nsig(1,it) = nsig(1,it) + msigm
      nsms(1,it) = nsms(1,it) + msmin
      nszr(1,it) = nszr(1,it) + mszer
      nspl(1,it) = nspl(1,it) + mspls
      nxi(1,it) = nxi(1,it) + mcsxi
      nxne(1,it) = nxne(1,it) + mxine
      nxch(1,it) = nxch(1,it) + mxich

      nlam(2,it) = nlam(2,it) + matlamb
      nsig(2,it) = nsig(2,it) + matsigm
      nsms(2,it) = nsms(2,it) + matsmin
      nszr(2,it) = nszr(2,it) + matszer
      nspl(2,it) = nspl(2,it) + matspls
      nxi(2,it) = nxi(2,it) + matcsxi
      nxne(2,it) = nxne(2,it) + matxine
      nxch(2,it) = nxch(2,it) + matxich

      dens(it) = dens(it) + rhoc

      ! Restored: [kkmin, kkmax] meson CSM
      if(it == kkmax) nsqp(it) = nsqp(it) + mpion**2
      if(it == kkmax) nsqp1(it) = nsqp1(it) + mpion1**2
      if(it == kkmax) nsqp2(it) = nsqp2(it) + mpion2**2
      if(it == kkmax) nsqp3(it) = nsqp3(it) + mpion3**2
      if(it == kkmax) nsqet(it) = nsqet(it) + meta**2
      if(it == kkmax) msakzr(it) = msakzr(it) + makzr**2
      if(it == kkmax) msakmn(it) = msakmn(it) + makmn**2
      if(it == kkmax) msigms(it) = msigms(it) + msmin**2
      if(it == kkmax) msigpl(it) = msigpl(it) + mspls**2
      if(it == kkmax) mslamd(it) = mslamd(it) + mlamb**2

      ! Restored: it=1 rtime(1) CSM
      if(L1 == 1.and.it == 1) KK00 = k*dt

      ! Restored comment. CSM
      cxx = cxx/float(nucl)
      cyy = cyy/float(nucl)
      czz = czz/float(nucl)
      pxx = pxx/float(nucl)
      pyy = pyy/float(nucl)
      pzz = pzz/float(nucl)

      ! Restored comment. CSM
      nump = 0
      numn = 0

****** Restored: **************** 3.3.2 ********** CSM

! Restored: !!!! CSM
      do 2893 i = 1, mass

      if(SURF_COAL == 1.and.VANISH(i) == 1)cycle
      rhoby = 0.                                          ! Restored: i CSM

      ! Restored: i CSM
      do 2894 j = 1, mass
          if(SURF_COAL == 1.and.VANISH(j) == 1)cycle
          widxx = sl(i) + sl(j)                           ! Restored: 2.*sl(i) CSM
          dr2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 
          rhoby = rhoby + xexp(-dr2/widxx)/(pi*widxx)**1.5
2894  continue

      ! Restored: i CSM
      rds = rhoby/rho0
      if(rds > 1.)then
          if(LB(i) == 1.and.idr(i) == 0) nump = nump + 1  ! Restored comment. CSM
          if(LB(i) == 0.and.idr(i) == 0) numn = numn + 1  ! Restored comment. CSM
      endif

      ! Restored comment. CSM
      x(i) = x(i) - cxx
      y(i) = y(i) - cyy
      z(i) = z(i) - czz
      xk(i) = xk(i) - pxx
      yk(i) = yk(i) - pyy
      zk(i) = zk(i) - pzz

      ! Restored: , CSM
      pbar = sqrt(xk(i)**2 + yk(i)**2 + zk(i)**2)
      embr = sqrt(em(i)**2 + xk(i)**2 + yk(i)**2 + zk(i)**2)
      trbr = sqrt(xk(i)**2 + yk(i)**2)

      ! Restored: , CSM
      ELM = SQRT(EM(i)**2 + zkpion(i)**2 + ykpion(i)**2 + xkpion(i)**2)
      BETALM = zkpion(I)/ELM
      YMR = 0.5*ALOG((1. + BETALM)/(1. - BETALM))

      ! Restored: , , CSM
      if(idms(i)  <=  0)goto 2890
      xpion(i) = xpion(i) - cxx
      ypion(i) = ypion(i) - cyy
      zpion(i) = zpion(i) - czz
      xkpion(i) = xkpion(i) - pxx
      ykpion(i) = ykpion(i) - pyy
      zkpion(i) = zkpion(i) - pzz

      ppin = sqrt(xkpion(i)**2 + ykpion(i)**2 + zkpion(i)**2)
      emp2 = sqrt(empion(i)**2 + xkpion(i)**2 + ykpion(i)**2 + 
     &zkpion(i)**2)
      trpn = sqrt(xkpion(i)**2 + ykpion(i)**2)
      if(trpn == 0.)then
            fldr = 0.
            flep = 0.
      else
          fldr = xkpion(i)/trpn                           ! Restored: pion =cos(fai) CSM
          flep = (xkpion(i)**2 - ykpion(i)**2)/trpn**2    ! Restored: pion =cos(2fai) CSM
      endif

      zpol = zkpion(i)/ppin
      IPL = zpol/DP

      ! longitudinal rapidity
      BETA = zkpion(i)/emp2
      YRAP = 0.5*ALOG((1.+BETA)/(1.-BETA))
      RYTP = YRAP                                         ! /YPRCM
      IYY = NINT(RYTP/DY)

      ! transverse momentum
      PT = SQRT(xkpion(i)**2 + ykpion(i)**2)
      IMM = NINT(PT/DM)                                   ! /PZCMPR

      ! kinetic energy
      XSEX = csf*emp2/ppin**2
      IKK =  NINT((emp2-empion(i))/DM)                    ! NINT(ppin/DM)

      if(IKK > 80)goto 2890
      if(abs(IYY) > 60.or.IMM > 60)goto 2890
      if(abs(icp(i)) >= 2)goto 2890

      
      ! Restored: FSVT( ) CSM
      if(FSVT == 1)then
          if(IDMS(I) == 1)then                            ! Restored: pi CSM
              ! Restored: pi-/pi+ CSM
              MKECM=EMP2-EMPION(I)
              RMKECM=0.005                                ! Restored: GeV),bin 1,2,...,50 CSM
              if(MKECM <= 0. .OR. MKECM > 0.25)goto 945  ! Restored: (0,0.25]GeV CSM
                RDMKECM=CEILING(MKECM/RMKECM)
              NKEPION(icp(i),RDMKECM)=NKEPION(icp(i),RDMKECM)+1
945              continue

              ! Restored: pi-/pi+ CSM
                RRPT=0.02                                   ! Restored: (GeV/c),bin 1,2,...,50 CSM
                if(PT <= 0. .OR. PT. GT. 1.)goto 950        ! Restored: (0,1]GeV/c CSM
              RDPT=CEILING(PT/RRPT)
              NTRPION(icp(i),RDPT)=NTRPION(icp(i),RDPT)+1
950              continue

              ! Restored: p-/p+ y/yproj CSM
                YDRM = abs(YRAP/YPR)
                RYDRM = 0.02                                ! Restored: y/yproj ,bin 1,2,...,50 CSM
                if(YDRM. LE.0. .OR. YDRM. GT. 1.)goto 9661  ! Restored: (0,1] CSM
              RDYDRM=CEILING(YDRM/RYDRM)
                NYPPION(ICP(I),RDYDRM)=NYPPION(ICP(I),RDYDRM)+1
9661              continue
          endif
      endif
2890  cycle
      
2893      continue

      ! Restored comment. CSM
      ! Restored comment. CSM
      if(FSVT. EQ. 1)then
            do M = -1, 1
                  do KK = 1, 50
                  ! Restored: pi-/pi+ CSM
                    NNKEPION(M,KK)=NNKEPION(M,KK)+NKEPION(M,KK)
                    NNTRPION(M,KK)=NNTRPION(M,KK)+NTRPION(M,KK)
                    NNYPPION(M,KK)=NNYPPION(M,KK)+NYPPION(M,KK) 
                    SNNKEPION(M,KK)=SNNKEPION(M,KK)+NKEPION(M,KK)**2
                    SNNTRPION(M,KK)=SNNTRPION(M,KK)+NTRPION(M,KK)**2
                    SNNYPPION(M,KK)=SNNYPPION(M,KK)+NTRPION(M,KK)**2

                  ! Restored: dM/dKEc.m. KEc.m. CSM
              DNKEPION(M,KK)=DNKEPION(M,KK)+NKEPION(M,KK)/(RMKECM*1000)
                  enddo
            enddo
      endif

      ! Restored comment. CSM
      do I = 1, MASS
            x(i)=x(i)+cxx
            y(i)=y(i)+cyy
            z(i)=z(i)+czz
            xk(i)=xk(i)+pxx
            yk(i)=yk(i)+pyy
            zk(i)=zk(i)+pzz

            if(idms(i)  <=  0)cycle
            xpion(i)=xpion(i)+cxx
            ypion(i)=ypion(i)+cyy
            zpion(i)=zpion(i)+czz
            xkpion(i)=xkpion(i)+pxx
            ykpion(i)=ykpion(i)+pyy
            zkpion(i)=zkpion(i)+pzz
      enddo
      
      ! Restored: ANGM_DEC( ) CSM
      if(ANGM_DEC == 1)then
            ANGMCOEF = 5.067731
            TOTANGMX = 0.
            TOTANGMY = 0.
            TOTANGMZ = 0.
            do i = 1, mg
                  TOTANGMX = TOTANGMX + Y(I)*ZK(I)-Z(I)*YK(I)
                  TOTANGMY = TOTANGMY + Z(I)*XK(I)-X(I)*ZK(I)
                  TOTANGMZ = TOTANGMZ + X(I)*YK(I)-Y(I)*XK(I)
        enddo
        ! Restored: hbar CSM
            TOTANGMX = TOTANGMX*ANGMCOEF
            TOTANGMY = TOTANGMY*ANGMCOEF
            TOTANGMZ = TOTANGMZ*ANGMCOEF
            if(TESTANG == 1)then
                  PRINT*, L1, K*DT
                  PRINT*, TOTANGMX, TOTANGMY, TOTANGMZ
                  PAUSE
            endif
      endif

      if(ECHECK == 1)then
        E_BEFORE = 0.
        do i = 1, mass
          if(idms(i) /= 0 .or. idr(i) > 0) cycle
          E_BEFORE = E_BEFORE +
     &      sqrt(em(i)**2 + xk(i)**2 + yk(i)**2 + zk(i)**2) - em(i)
        enddo
      endif

ccccc  2500 ccccc Restored comment. CSM
      if(COLM == 1)then
      
      ! Restored comment. CSM
      do I = 1, MG
        XJR(I,1)=0.
        XJR(I,2)=0.
        XJR(I,3)=0.
        IC(I)=0
        IFLAG(I)=0
        mlll(I)=0
        mll1(I)=0
        mllz(I)=0
        mlz1(I)=0
        CLX(I)=0.
        CLY(I)=0.
        CLZ(I)=0.
        CLKX(I)=0.
        CLKY(I)=0.
        CLKZ(I)=0.
        CLKIN(I)=0.
        ICL(I)= 0
        ICHY(I)=0
        ISUC(I)=0
        NNHY(0:1,0:20,I)=0
      enddo

****** Restored: ************* 3.3.3 (MST) ******************** CSM
      ! Restored: : IC(1:mass), ID CSM
      III = 0   
      MST = 1
      
ccccc MST 最小树法 ccccc     
      if (MST) then

      do 120 I1 = 1, mass
        if (SURF_COAL == 1 .and. VANISH(I1) == 1) cycle
        if (idr(I1) > 0) cycle
        if (idat(I1) > 0) cycle
        ! Restored: , IC CSM
        ! Restored: cascade ID, I1 I2 , IC(I1)=IC(I2) CSM
        if (IC(I1) == 0) then 
          III = III + 1 
          IC(I1) = III                             
        endif
        do 130 I2 = I1+1, mass
          if(SURF_COAL == 1.and.VANISH(I2) == 1)cycle
          if(idr(I2)  >  0)cycle
          if(idat(I2)  >  0)cycle

          ! Restored: i1 i2 CSM
          DRX = X(I1) - X(I2)                                             
          DRY = Y(I1) - Y(I2)                                              
          DRZ = Z(I1) - Z(I2)                                             
          DRR = DRX*DRX + DRY*DRY + DRZ*DRZ                                                               
          DPX = XK(I1) - XK(I2)                                        
          DPY = YK(I1) - YK(I2)                                            
          DPZ = ZK(I1) - ZK(I2)                                            
          DPP = DPX*DPX + DPY*DPY + DPZ*DPZ
          if(DRR > RCR .or. DPP > PCR)cycle
          

          if(IC(I2) == 0)then                         ! Restored: 1234 CSM
            ! Restored: i1 CSM
            i1Z = 0
            i1N = 0
            i1A = 0  
            xI1 = 0.
            yI1 = 0.
            zI1 = 0.
            do k1 = 1, mass
              if(SURF_COAL == 1.and.VANISH(K1) == 1)cycle
              !if(Iram == 1)then
              if(IC(K1) == IC(I1))then          ! I1 cluster
                xI1 = xI1 + X(K1)
                yI1 = yI1 + Y(K1)
                zI1 = zI1 + Z(K1)
                i1A = i1A + 1                   ! Restored: i1 cluster CSM
              endif              
              ! calculate Z and N of residue
              if(IC(K1) == IC(I1).and.LB(K1) == 1)then     ! feng
                i1Z = i1Z + 1                   ! Restored: i1 cluster CSM
              elseif(IC(K1) == IC(I1).and.LB(K1) == 0)then
                i1N = i1N + 1                   ! Restored: i1 cluster CSM
              endif
            enddo
            !if(Iram == 1)then
            ! Restored: i1 cluster CSM
            xI1 = xI1/float(i1A)
            yI1 = yI1/float(i1A)
            zI1 = zI1/float(i1A)
      
            ! calculate the r square in center of mass frame
            xxI1 = 0.
            yyI1 = 0.
            zzI1 = 0.
            do k1 = 1, mass
              if(SURF_COAL == 1.and.VANISH(K1) == 1)cycle
              if(IC(K1) == IC(I1))then
                xxI1 = xxI1 + (xI1-X(K1))*(xI1-X(K1))
                yyI1 = yyI1 + (yI1-Y(K1))*(yI1-Y(K1))
                zzI1 = zzI1 + (zI1-Z(K1))*(zI1-Z(K1))
              endif
            enddo
            if(i1A < 1)cycle
            !endif
                  
            !if(Iram == 1)then
            if(i1A >= 2)then
              crms1 = sqrt((xxI1+yyI1+zzI1)/float(i1A))
              DDx = xI1 - X(I2)      
              DDy = yI1 - Y(I2)      
              DDz = zI1 - Z(I2)
              DD = sqrt(DDx*DDx + DDy*DDy + DDz*DDz)
              if(DD > (crms1+rr0)) goto 130
            endif
            !endif

            !if(LB(I2) == 1) then
            !   if(ISOTAB(i1Z+1,i1N) == 0)goto 130
            !else
            !   if(ISOTAB(i1Z,i1N+1) == 0)goto 130
            !endif
            IC(I2) = IC(I1) ! Restored: I2 I1 CSM
          else                                    ! Restored: 1234 CSM
            i2Z = 0
            i2N = 0
            i2A = 0
            xI2 = 0.
            yI2 = 0.
            zI2 = 0.
            do K1 = 1, mass
            if(SURF_COAL == 1.and.VANISH(K1) == 1)cycle
              !if(Iram == 1) then
              if(IC(K1) == IC(I2).OR.IC(K1) == IC(I1))then ! Restored: I1 I2 CSM
                xI2 = xI2 + X(K1)
                yI2 = yI2 + Y(K1)
                zI2 = zI2 + Z(K1)
                i2A = i2A + 1
              endif
              !endif
              if(IC(K1) == IC(I2) .and. LB(K1) == 1)then
                I2Z = I2Z + 1
              elseif(IC(K1) == IC(I2) .and. LB(K1) == 0)then
                I2N = I2N + 1
              endif
            enddo

            !if(Iram == 1) then
            xI2 = xI2/float(i2A)
            yI2 = yI2/float(i2A)
            zI2 = zI2/float(i2A)
            xxI2 = 0.
            yyI2 = 0.
            zzI2 = 0.
            do K1 = 1, mass
              if(SURF_COAL == 1.and.VANISH(K1) == 1)cycle
              if(IC(K1) == IC(I2).OR.IC(K1) == IC(I1))then
                xxI2=xxI2+(xI2-X(K1))*(xI2-X(K1))
                yyI2=yyI2+(yI2-Y(K1))*(yI2-Y(K1))
                zzI2=zzI2+(zI2-Z(K1))*(zI2-Z(K1))
              endif
            enddo
            if(i2A  <  1)cycle
            ! Restored comment. CSM
            crms2=sqrt((xxI2+yyI2+zzI2)/i2A)
            Rc=1.14*(float(i2A))**(1./3.)
            if(crms2 > Rc)goto 130
            !endif
            !if(ISOTAB(i1Z+I2Z,i1N+I2N) == 0)goto 130
            ICK = MAX0(IC(I1),IC(I2))                                        
            ICS = MIN0(IC(I1),IC(I2))                                        
            do J = 1, mass
              if(SURF_COAL == 1.and.VANISH(J) == 1)cycle
              if(IC(J) == ICK) IC(J) = ICS
            enddo
1234      endif                                   ! 1234       
          
130     continue                                  ! 130 I2     
                
120   continue                                    ! 120 I1    
      
      else  ! Restored: MST 2501 CSM
      
      do 140 I1 = 1, mass
            if(SURF_COAL == 1.and.VANISH(I1) == 1)cycle
        if(idr(I1) > 0)goto 140
        if(idat(I1) > 0)goto 140
        ! Restored: , IC CSM
        ! Restored: cascade ID, I1 I2 , IC(I1)=IC(I2) CSM
        if(IC(I1) == 0)then 
          IC(I1) = 1  ! Restored: : , , CSM
        endif
140   continue
      
2501  endif      
ccccc  MST 2501 ccccc Restored comment. CSM

      
****** Restored: **************** 3.3.4 ******************** CSM
C028 BEGIN COMPACT CLUSTERS
C Compact membership before any per-cluster attributes are accumulated.
      if(any(IC(1:mass) < 0) .or. any(IC(1:mass) > MG))
     &  error stop 'Invalid MST cluster label'
      max_val = maxval(IC(1:mass))
      if(allocated(mmap)) deallocate(mmap)
      allocate(mmap(max_val))
      mmap = 0
      do I=1,mass
        if(IC(I) > 0) mmap(IC(I)) = -1
      enddo
      ICMAX = 0
      do I=1,max_val
        if(mmap(I) /= -1) cycle
        ICMAX = ICMAX + 1
        mmap(I) = ICMAX
      enddo
      do I=1,mass
        if(IC(I) > 0) IC(I) = mmap(IC(I))
      enddo
C028 END COMPACT CLUSTERS
      ! Restored: ( ) CSM
      do I = 1, mass
        if(SURF_COAL == 1.and.VANISH(I) == 1)cycle
        IS = IC(I)                                          
        if(IS == 0)cycle                            ! Restored: cluster ( ) CSM
        if(LB(I) == 1)then                                           
          XJR(IS,1) = XJR(IS,1) + 1.
        else                                                   
          XJR(IS,2) = XJR(IS,2) + 1.
        endif
        ! Restored comment. CSM
        CLX(IS) = CLX(IS) + X(I)                      
        CLY(IS) = CLY(IS) + Y(I)
        CLZ(IS) = CLZ(IS) + Z(I)
        CLKX(IS) = CLKX(IS) + XK(I)
        CLKY(IS) = CLKY(IS) + YK(I)
        CLKZ(IS) = CLKZ(IS) + ZK(I)
          
        ICL(IS) = ICL(IS) + 1
        ICHY(IS) = XJR(IS,1)
        NNHY(0,0,IS) = NNHY(0,0,IS) + 1
      enddo
      
      print*, 'ICMAX = ', ICMAX
      
      ! add ID of pre-cluster (emitted)
      if(SURF_COAL == 1)then
        IICMAX = ICMAX 
        do MM = 1, 20
          ICMAX = ICMAX + 1
          if(SURF_Z(MM) == 0 .and. SURF_A(MM) == 0)then
            ICMAX = ICMAX - 1
            EXIT
          else
            if(ICMAX > MG) error stop 'Surface cluster capacity'
            XJR(ICMAX, 1) = SURF_Z(MM)
            XJR(ICMAX, 2) = SURF_A(MM) - SURF_Z(MM)
          endif
          CLKX(ICMAX) = SURFMOM(MM, 1)
          CLKY(ICMAX) = SURFMOM(MM, 2)
          CLKZ(ICMAX) = SURFMOM(MM, 3)
          CLKIN(ICMAX) = SURFE(MM)*0.001
          IFSURF(ICMAX) = 1      
        enddo
      endif      
      
      ! Mass distribution statistics
      ! Count cluster masses
      !do I = 1, nucl
      !  XMTOT = XJR(I,1) + XJR(I,2)                                   
      !  if(XMTOT > 1.0)then                                          
      !    MLL = INT( XMTOT + 0.5 )                             
      !    MLLL(MLL) = MLLL(MLL) + 1               ! Number of fragments with mass A (>=2)                                                                                          
      !  endif
      !enddo                                                       
      !! Expand into fragment list MLL1 / MCL
      !NCCP1 = 0                                                  
      !do I = 1, nucl
      !  if(MLLL(I) > 0)then                                       
      !    do IKK = NCCP1+1, NCCP1+MLLL(I)      
      !      MLL1(IKK)=I                        
      !    enddo                              
      !    NCCP1=NCCP1+MLLL(I)                
      !  endif                                 
      !enddo    
      !! For A=1 "gas" (free nucleons) also counted as fragments
      !XGAS = 0.0      
      !do I = 1, nucl  
      !  XMTOT=XJR(I,1)+XJR(I,2)                             
      !  if(XMTOT > 0.and.XMTOT <= 1.0)then
      !    XGAS=XGAS+XMTOT                    
      !  endif 
      !enddo
      !! Build final fragment mass list MCL
      !NCCP= NCCP1 + INT(XGAS+0.5)
      !do I=1,NCCP1
      !  MCL(I) = MLL1(I)                                         
      !enddo                                                       
      !do I=NCCP1+1, NCCP                                         
      !  MCL(I) = 1                                                
      !enddo                                                      
      !
      !! Charge (Z) distribution statistics (similar to mass)
      !do I=1, nucl
      !  XZTOT=XJR(I,1)                     ! XJR(I,1) is the charge number Z
      !  if(XZTOT > 1.0)then 
      !    MLZ   = INT( XZTOT + 0.5 )
      !    MLLZ(MLZ) = MLLZ(MLZ) + 1 
      !  endif                                                           
      !enddo
      !NCCPZ1=0      
      !do I=1, nucl
      !  if(MLLZ(I) > 0)then 
      !    do IKK=NCCPZ1+1, NCCPZ1+MLLZ(I)
      !      MLZ1(IKK)=I 
      !    enddo       
      !    NCCPZ1=NCCPZ1+MLLZ(I)
      !  endif                   
      !enddo                
      !XGASZ = 0.0           
      !do I=1, nucl
      !  XZTOT=XJR(I,1)  
      !  if(XZTOT > 0.and.XZTOT <= 1.0)then 
      !    XGASZ=XGASZ+XZTOT                   
      !  endif                               
      !enddo
      !NCCPZ = NCCPZ1 + INT(XGASZ + 0.5)                      
      !do I=1,NCCPZ1                                               
      !  MCLZ(I) = MLZ1(I)                                     
      !enddo                                                      
      !do I = NCCPZ1+1, NCCPZ                                      
      !  MCLZ(I) = 1                                             
      !enddo                                                  

      if(ICMAX > mg)stop 'MG is too small !'

      ! Restored: cascade CSM
      EMD_CL(1:MG) = 1
      IF_CASCADE(1:MG) = 0
      if(SURF_COAL == 1)then
        do I = 1, MASS
          if(VANISH(I) == 1 .or. IC(I) == 0)cycle       
          ! EMD_CL distinguish
          if((XJR(IC(I),1) + XJR(IC(I),2)) <= 4)then    ! CSM
            IF_CASCADE(I) = 1                           ! MST cluster
            EMD_CL(IC(I)) = EMD_CL(IC(I))*EMD_NUCL(I)   ! EMD_NUCL is input, and denotes the emitted nucleon
          endif
        enddo
      endif

****** Restored: **************** 3.3.5 ******************** CSM

C028 Account only accepted nucleons and explicitly supplied surface clusters.
      ledger_expected_a = 0
      ledger_expected_z = 0
      do I=1,mass
        if(SURF_COAL == 1.and.VANISH(I) == 1) cycle
        if(IC(I) <= 0) cycle
        ledger_expected_a = ledger_expected_a + 1
        ledger_expected_z = ledger_expected_z + LB(I)
      enddo
      if(SURF_COAL == 1)then
        do MM=1,20
          if(SURF_A(MM) == 0.and.SURF_Z(MM) == 0) exit
          ledger_expected_a = ledger_expected_a + SURF_A(MM)
          ledger_expected_z = ledger_expected_z + SURF_Z(MM)
        enddo
      endif
      ledger_hot_z = NINT(SUM(XJR(1:ICMAX,1)))
      ledger_hot_a = ledger_hot_z + NINT(SUM(XJR(1:ICMAX,2)))
      if(ledger_hot_a /= ledger_expected_a.or.
     &   ledger_hot_z /= ledger_expected_z)then
        write(*,*) 'MST_HOT_MISMATCH',L1,IT,ledger_expected_a,
     &    ledger_expected_z,ledger_hot_a,ledger_hot_z
        error stop 1
      endif

! Restored: CNPCS A N Z CSM
      ! HOT FRAGMENTS (ISOTH:isotope hot)
      if(CNPCS == 1)      then
            do IIC=1,ICMAX
                  if((XJR(IIC,1) /= 0).OR.(XJR(IIC,2) /= 0))then
                      NZISOTH(XJR(IIC,1)) = NZISOTH(XJR(IIC,1)) + 1                                                
                      NNISOTH(XJR(IIC,2)) = NNISOTH(XJR(IIC,2)) + 1
              NAISOTH(XJR(IIC,1)+XJR(IIC,2)) = 
     &          NAISOTH(XJR(IIC,1)+XJR(IIC,2))+1
                  endif
            enddo
      endif

! Restored: TEELECS CSM
      if(TEELECS == 1)then
            do IIC = 1, ICMAX
                  do I = 1, TEELESPECIES
                        if(TEELECS_SETINF(I,1) == XJR(IIC,1).and.
     &                TEELECS_SETINF(I,2) == XJR(IIC,2))then
                      TEELECSCOUNT(k,I,1) = TEELECSCOUNT(k,I,1) + 1
                        endif
                  enddo
            enddo
      endif
      
      ! Restored comment. CSM
      if(SURF_COAL == 1)then
            do MM = 1, 20
                  if(SURF_Z(MM) == 0.and.SURF_A(MM) == 0)then
                        EXIT
                  else
                        XJR(ICMAX, 1) = 0
                        XJR(ICMAX, 2) = 0
                  endif
                  CLKX(ICMAX) = 0
                  CLKY(ICMAX) = 0
                  CLKZ(ICMAX) = 0
                  CLKIN(ICMAX) = 0
 
                  IFSURF(ICMAX) = 0      
          ICMAX = ICMAX - 1
            enddo
      endif      

! ANGM_DEC: calculate the remnant spin (only one big remnant)
      if (ANGM_DEC == 1) then                  
        NUMCL = 0
        do ICC = 1, ICMAX
          ! sometimes ang_light = 10, if A were greater than it, the nuclear is big remnant  
          if (ICL(ICC) > ANGM_LIGHT) NUMCL = NUMCL + 1
        end do
        if (NUMCL > 1) goto 5584 ! here is a bug: if the remnant number 1 is not the bigest one?
        ANGM_CLX(1:MG) = 0.
        ANGM_CLY(1:MG) = 0.
        ANGM_CLZ(1:MG) = 0.
        do ICC = 1, ICMAX
          if (ICL(ICC) <= 0) cycle
          ANGM_CLX(ICC) = CLY(ICC)*CLKZ(ICC) - CLZ(ICC)*CLKY(ICC)
          ANGM_CLY(ICC) = CLZ(ICC)*CLKX(ICC) - CLX(ICC)*CLKZ(ICC)
          ANGM_CLZ(ICC) = CLX(ICC)*CLKY(ICC) - CLY(ICC)*CLKX(ICC)
          ANGM_CLX(ICC) = ANGM_CLX(ICC)*ANGMCOEF/ICL(ICC)
          ANGM_CLY(ICC) = ANGM_CLY(ICC)*ANGMCOEF/ICL(ICC)
          ANGM_CLZ(ICC) = ANGM_CLZ(ICC)*ANGMCOEF/ICL(ICC)
          TOTANGMX1 = TOTANGMX1 - ANGM_CLX(ICC)
          TOTANGMY1 = TOTANGMY1 - ANGM_CLY(ICC)
          TOTANGMZ1 = TOTANGMZ1 - ANGM_CLZ(ICC)                  
        enddo
        ! expect the pion contribution of angular momentum
        do IP = 1, MASS
          ANGM_PIONX = 0.
          ANGM_PIONY = 0.
          ANGM_PIONZ = 0.
          ANGM_PIONX = YPION(IP)*zkpion(IP) - ZPION(IP)*ykpion(IP)
          ANGM_PIONY = ZPION(IP)*xkpion(IP) - XPION(IP)*zkpion(IP)
          ANGM_PIONZ = XPION(IP)*ykpion(IP) - YPION(IP)*xkpion(IP)
          ANGM_PIONX = ANGM_PIONX*ANGMCOEF
          ANGM_PIONY = ANGM_PIONY*ANGMCOEF
          ANGM_PIONZ = ANGM_PIONZ*ANGMCOEF
          TOTANGMX1 = TOTANGMX1 - ANGM_PIONX
          TOTANGMY1 = TOTANGMY1 - ANGM_PIONY
          TOTANGMZ1 = TOTANGMZ1 - ANGM_PIONZ
        enddo
        ! expect the emitted nucleon contribution
        if (SURF_COAL == 1) then
          do MM = 1, 20
            if (SURF_Z(MM) == 0 .and. SURF_A(MM) == 0) EXIT
            TOTANGMX1 = TOTANGMX1 - SURFAGM(MM, 1)
            TOTANGMY1 = TOTANGMY1 - SURFAGM(MM, 2)
            TOTANGMZ1 = TOTANGMZ1 - SURFAGM(MM, 3)
          enddo
        endif
        ! quantum correction
        TOTANGM3 = TOTANGMX1**2 + TOTANGMY1**2 + TOTANGMZ1**2
        TOTANGM3 = SQRT(1. + 4*TOTANGM3) - 1.
        TOTANGM3 = ABS(TOTANGM3/2.)
      endif
5584  continue
      
******************************************************************************
******************************************************************************
******************************************************************************
      
      
      
******************************************************************************
CCC******************* 4 GEMINI                         *******************CCC
******************************************************************************
      
      
******************************************************************************
****** Restored: **************** 4.1 (TESTEX) ********************** CSM

      ! Restored: , CSM
!      if(TESTEX == 1)      then
!            MASS1=0
!            do I=1,MASS
!                  IC(I)=0
!                  X(I)=0.
!                  Y(I)=0.
!                  Z(I)=0.
!                  XK(I)=0.
!                  YK(I)=0.
!                  ZK(I)=0.
!                  EM(I)=0.
!            enddo
!
!            open(331,FILE='tini.dat', STATUS='unknown')
!            do      I = 1,MAss
!                  read(331,77777,end=77771)      x(i),y(i),z(i),xK(i),yK(i)
!     &           ,zK(i),em(i),LB(I),ISY(I),IDII
!                  IC(I) = 1
!                  MASS1 = MASS1 + 1
!            enddo
!77777            format(1x,7f9.4,2x,2i3,2x,i3)
!77771            continue
!            close(331)
!      endif

      ! Restored comment. CSM
      do I = 1, MASS
        if(SURF_COAL == 1.and.VANISH(I) == 1)cycle
        IS = IC(I)
        if(IS == 0)cycle
        if(ICL(IS) <= 0.)cycle

        if(TESTEX == 0)then
          ACX = CLX(IS)/ICL(IS)
          ACY = CLY(IS)/ICL(IS)
          ACZ = CLZ(IS)/ICL(IS)
          ACKX= CLKX(IS)/ICL(IS)
          ACKY= CLKY(IS)/ICL(IS)
          ACKZ= CLKZ(IS)/ICL(IS)

          XFR(I)= X(I) - ACX
          YFR(I)= Y(I) - ACY
          ZFR(I)= Z(I) - ACZ
          XKFR(I)= XK(I)-ACKX
          YKFR(I)= YK(I)-ACKY
          ZKFR(I)= ZK(I)-ACKZ
        elseif(TESTEX == 1)then
          XFR(I)= X(I)
          YFR(I)= Y(I)
          ZFR(I)= Z(I)
          XKFR(I)= XK(I)
          YKFR(I)= YK(I)
          ZKFR(I)= ZK(I)
        endif

        fr(I)=0.
        frp(I)=0.
        tkin(I)=0.
        utwo(I)=0.
        uthr(I)=0.
        ucou(I)=0.
        umod(I)=0.
        usym(I)=0.
        ussr(I)=0.
        ugrd(I)=0.
        couex(I)=0.
        utau(I)=0.
        aene(I)=0.
        excen(I)=0.
        bdgd(I)=0.
      enddo

      do I1 = 1, mass
        if (SURF_COAL == 1 .and. VANISH(I1) == 1) cycle
        IS1 = IC(I1)
        if (IS1 == 0) cycle
        ! TESTEX
        if (XJR(IS1,1) <= MCFR .and. TESTEX /= 1) cycle
                
        ! two-body interaction
        do I2 = I1+1, mass
          if (SURF_COAL == 1 .and. VANISH(I2) == 1) cycle
          IS2 = IC(I2)
          if (IS2 == 0) cycle
          ! TESTEX
          if (XJR(IS2,1) <= MCFR .and. TESTEX /= 1) cycle

          al = 2.*(sl(I1) + sl(I2))
          ! TESTEX
          if (TESTEX == 1) then
            al = 2.*(2*(0.08*float(MASS1)**0.3333333 + 0.92)**2.)
          end if
          sal = (pi*al)**1.5*rho0
          dr2 = (XFR(i1)-XFR(i2))**2 + (YFR(i1)-YFR(i2))**2 
     &        + (ZFR(i1)-ZFR(i2))**2
          p2 = (XKFR(i1)-XKFR(i2))**2 + (YKFR(i1)-YKFR(i2))**2          
     &        + (ZKFR(i1)-ZKFR(i2))**2
          pbr = p2*cmeps + 1.    ! momentum dependent item correction
          t2 = xexp(-dr2/al)/sal ! Overlapping density of two body wave packets
          idr22 = dr2*40. + 1.5 
          idr2 = min0(idr22, 100000) 
          nccc = ih(idr2)
          ! C_{/tau,/tau'}
          if (LB(i1) == LB(i2)) then
            cfm = cmnn
          else
            cfm = cmnp
          end if
          if (IS2 == IS1) then
            ucou(IS1) = ucou(IS1) + LB(i1)*LB(i2)*scou(nccc,al)/2.     ! Coulumb energy
            ugrd(IS1) = ugrd(IS1) + gsur*(6./al - 4.*dr2/al**2.)*t2/2. ! Surface gradient potential energy
            ussr(IS1) = ussr(IS1) + 2*gsym*(6./al - 4.*dr2/al**2.)     ! Surface symmetry energy correction
     &                * isy(i1)*isy(i2)*t2/2.
            usym(IS1) = usym(IS1) + 2*symm*isy(i1)*isy(i2)*t2/2.       ! Volume symmetry energy
            umod(IS1) = umod(IS1) + cfm*alog(pbr)*alog(pbr)*t2/2.      ! Momemtum dependent item
            
            ucou(IS2) = ucou(IS2) + LB(i1)*LB(i2)*scou(nccc,al)/2.     ! Coulumb energy
            ugrd(IS2) = ugrd(IS2) + gsur*(6./al - 4.*dr2/al**2.)*t2/2. ! Surface gradient potential energy
            ussr(IS2) = ussr(IS2) + 2*gsym*(6./al - 4.*dr2/al**2.)     ! Surface symmetry energy correction
     &                * isy(i1)*isy(i2)*t2/2.
            usym(IS2) = usym(IS2) + 2*symm*isy(i1)*isy(i2)*t2/2.       ! Volume symmetry energy
            umod(IS2) = umod(IS2) + cfm*alog(pbr)*alog(pbr)*t2/2.      ! Momemtum dependent item
            
            ! Sum the local overlap density factor fr and the charge density factor frp
            ! at the current position of the nucleon
            fr(i1) = fr(i1) + t2
            fr(i2) = fr(i2) + t2
            frp(i1) = frp(i1) + LB(i1)*LB(i2)*t2*rho0
            frp(i2) = frp(i2) + LB(i1)*LB(i2)*t2*rho0
          end if 
        end do

        ! single-body
        tkin(IS1) = tkin(IS1) + (XKFR(I1)*XKFR(I1) + YKFR(I1)*YKFR(I1) ! Non-relativistic Fermi kinetic energy
     &           + ZKFR(I1)*ZKFR(I1))/2./em(I1)
        utau(IS1) = utau(IS1) + gtau11*fr(I1)**ETA11                   ! Skyrme kinetic energy correction term
        utwo(IS1) = utwo(IS1) + 0.5*alpha*fr(I1)                       ! Skyrme 2-body attractive potential
        uthr(IS1) = uthr(IS1) + 1./(1. + gmmat)*batat*fr(I1)**gmmat    ! Skyrme 3-body repulsive potential
        couex(IS1) = couex(IS1) + cex*frp(I1)**0.3333333               ! Coulomb exchange potential correction        
      end do
      
      ! Restored comment. CSM
      Mul_MCFR = 0
      do I = 1, ICMAX
        ! ground state
        bdgd(I) = xmp*XJR(I,1) + xmn*XJR(I,2)  
     &          - XBIN(XJR(I,1), XJR(I,2)+XJR(I,1))      ! MeV

        if (XJR(I,1) > MCFR) then
          Mul_MCFR = Mul_MCFR + 1
          aene(i) = tkin(i) + ucou(i) + ugrd(i) + ussr(i) + usym(i)
     &            + umod(i) + utwo(i) + uthr(i) + couex(i) + utau(i)
        else      ! Restored: Y CSM
          aene(i) = -bdgd(I)/1000.
        end if
        excen(i) = aene(i)*1000 + bdgd(I)      ! MeV
        if (excen(i) < 0.) then
          excen(i) = 0.
        end if
      end do

      ! Restored comment. CSM
      do I = 1, ICMAX
        ! Ekin=Etot-Erest
        ! Etot=sqrt(Erest+p^2)
        ! Erest=m0*c^2
        CLKIN(I) = SQRT(DBLE(XJR(I,1)*AMU1+XJR(I,2)*AMU2+AENE(I))**2
     &           + DBLE(CLKX(I))**2+DBLE(CLKY(I))**2+DBLE(CLKZ(I))**2)
     &           - XJR(I,1)*AMU1 - XJR(I,2)*AMU2-AENE(I) ! GeV

        IAA = XJR(I,1) + XJR(I,2)
        ! Restored: MeV CSM
        DIFFUSION(I) = 1000.*CLKIN(I)/IAA
            
        ! PRINTCLKIN.1
c        print*, 'ICMAX = ', ICMAX
        if (PRINTCLKIN == 1) then
          if (L1 == 1.and.I == 1) then
            PRINT*
            PRINT*, 'Enter the leading event from which to analyze'
            read (*,*) LEADING1
          end if

          if (L1 >= LEADING1) then
            if (I == 1) then
              PRINT*
              PRINT*,'EVENT:',L1
              PRINT*
              PRINT*,'Hot fragment kinetic energies'
              PRINT*
              PRINT*,'         ','Z','           ','N','    ',
     &          'KIN/MeV','       ','AVKIN/MeV','       ',
     &          'THETA/Deg.'
            end if

            print*, 'I = ', I, 'A = ', IAA
            if (IAA /= 0) then
              PRINT*,INT(XJR(I,1)), INT(XJR(I,2)), CLKIN(I)*1000.,
     &          1000.*CLKIN(I)/IAA,ACOSD(CLKZ(I)/
     &          SQRT(CLKX(I)**2 + CLKY(I)**2+CLKZ(I)**2))
            end if
          end if
        end if
      end do

      if (gemini_in_out == 1) then
        max_val_index = maxloc(ICL(1:ICMAX))
        max_cluster_id = max_val_index(1)
        if (max_cluster_id > 0 .and. ICL(max_cluster_id) > 0) then      
          max_cl_mass_num = ICL(max_cluster_id)      
          max_cl_prot_num = int(XJR(max_cluster_id, 1))          
          max_cl_spin = totangm3
          max_cl_ex_energy = excen(max_cluster_id)

          if (L1 == 1) write(182,*) "   ID      A      Z   ",
     &      "    J       E*/MeV"
          write(182,188) max_cluster_id, max_cl_mass_num, 
     &      max_cl_prot_num, max_cl_spin, max_cl_ex_energy      
        else
          print *, "警告：本事件未发现任何有效碎片！"
        end if      
188     format(3I7, 2F9.2)      
C018    read(187,*) task018_event,task018_pa_a,task018_pa_z,
C018 &       task018_pa_ex
C018    if(task018_event /= L1 .or. task018_pa_a /= 233 .or.
C018 &       task018_pa_z /= 91) stop 'TASK018 energy mismatch'
      end if

      ! Restored comment. CSM

      ! Restored: tinii CSM
      if (TESTEX == 1) then
        PRINT*
        PRINT*,'kinetic energy=      tkin=  ',
     &        tkin(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'coulomb potential=   ucou=  ',
     &        ucou(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'coulomb exchange=   couex=  ',
     &        couex(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'yukawa potential=    uyuk= ',
     &        uyuk11*1000./float(mass1), 'MeV/u'
        PRINT*,'isospin potential=   usym= ',
     &        usym(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'surface symmetry=    ussr= ',
     &        ussr(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'pauli potential=    upaul= ',
     &        upaul11*1000./float(mass1), 'MeV/u'
        PRINT*,'mom. dep. potential=umomd= ',
     &        umod(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'skyrme potential=    usky= ',
     &        (utwo(1)+uthr(1))*1000./float(mass1), 'MeV/u'
        PRINT*,'skyrme  (2-body)=    utwo= ',
     &        utwo(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'skyrme  (3-body)=    uthr= ',
     &        uthr(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'gradient potential=ugradi= ',
     &        ugrd(1)*1000./float(MASS1),'MeV/u'
        PRINT*,'effective mass=      utau= ',
     &        utau(1)*1000./float(mass1), 'MeV/u'
        PRINT*,'shell correction=     Esh= ',
     &        Esh11*1000./float(mass1), 'MeV/u'
        PRINT*,'average energy=      aene= ',
     &        (tkin(1)+ucou(1)+ugrd(1)+ussr(1)+usym(1)+umod(1)+
     &        utwo(1)+uthr(1)+couex(1)+utau(1))*1000.
     &        /float(mass1), 'MeV/u'
        PRINT*
        goto 52013      
      elseif (OUTENE == 1) then
        do I = 1, ICMAX
          if(XJR(I,1)  <=  MCFR)cycle
              
          IAA=XJR(I,1)+XJR(I,2)
          PRINT*,'EVENT: ',L1                  
          PRINT*
          PRINT*,'REMNANT: ',I, ' ', 'Z= ', int(XJR(I,1)),
     &            ' ', 'A= ', ICL(IC(I))
          PRINT*
          PRINT*,'kinetic energy=     tkine=  ',
     &            tkin(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'coulomb potential=   ucou=  ',
     &            ucou(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'coulomb exchange=   couex=  ',
     &            couex(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'yukawa potential=    uyuk= ',
     &             uyuk11*1000./float(IAA), 'MeV/u'
          PRINT*,'isospin potential=   usym= ',
     &            usym(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'surface symmetry=    ussr= ',
     &            ussr(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'pauli potential=    upaul= ',
     &            upaul11*1000./float(IAA), 'MeV/u'
          PRINT*,'mom. dep. potential=umomd= ',
     &            umod(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'skyrme potential=    usky= ',
     &            (utwo(I)+uthr(I))*1000./float(IAA), 'MeV/u'
          PRINT*,'skyrme  (2-body)=    utwo= ',
     &            utwo(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'skyrme  (3-body)=    uthr= ',
     &            uthr(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'gradient potential=ugradi= ',
     &            ugrd(I)*1000./float(IAA),'MeV/u'
          PRINT*,'effective mass=      utau= ',
     &            utau(I)*1000./float(IAA), 'MeV/u'
          PRINT*,'shell correction=     Esh= ',
     &            Esh11*1000./float(IAA), 'MeV/u'
          PRINT*
          PRINT*,'groung state energy= bdgd= ',
     &            bdgd(i)/float(IAA), 'MeV/u'
          PRINT*,'average energy=      aene= ',
     &            aene(i)*1000./float(IAA), 'MeV/u'
          PRINT*,'excitation energy=    exc= ',
     &            excen(i)/float(IAA), 'MeV/u'
          print*
          PRINT*,'kinetic energy=     clkin= ',
     &            CLKIN(I)*1000, 'MeV'
          PRINT*
          PAUSE
        enddo
      endif

! check the energy conver
      if (ECHECK == 1) then
        E_AFTER = 0.
        do I = 1, ICMAX
          IZCL = XJR(I,1)
          IACL = XJR(I,2)
          if(IZCL+IACL  <=  0.) cycle
          EREST = IZCL*AMU1 + IACL*AMU2
          ETOT = sqrt((EREST + aene(I))**2
     &         + CLKX(I)**2 + CLKY(I)**2 + CLKZ(I)**2)
          E_AFTER = E_AFTER + (ETOT - EREST)
        enddo
        DELE = E_BEFORE - E_AFTER
        RELE = 0.
        if(E_BEFORE  >  0.) RELE = DELE / E_BEFORE
        if(FIRST_ECHECK)then
          write(183,'(A)') '# EVENT  IT   E_kin_bef(MeV)'
     &    // '  E_kin_aft(MeV)   DeltaE(MeV)     DeltaE/E'
          FIRST_ECHECK = .FALSE.
        endif
        write(183,9919) L1, IT, E_BEFORE*1000., E_AFTER*1000.,
     &        DELE*1000., RELE
 9919   format(I8, I4, 4ES16.8)
      end if

! Restored: TEAE l CSM
      if (TEAE == 1) then
        MAXXACL = 0
        IIHR = 0   
        do I = 1, ICMAX
          if (XJR(I,1)+XJR(I,2) >= MAXXACL) then
            MAXXACL = XJR(I,1) + XJR(I,2)
            IIHR = I
          end if
        end do
        excper = excen(IIHR)/float(MAXXACL)
9665    format(3f9.4)

        ! multi-fragmentation do not count
        if (Mul_MCFR /= 1 .and. IFFRAG == 0) IFFRAG = 1

        ! find the bigest excitation over time evolution
        if (EXCPER > MAXEXC .and. IFFRAG == 0 .and.(K*DT) >= EXPNT) then
          MAXEXC = EXCPER
          AMAXEXC = MAXXACL
        end if

        COMPEXC(IT,1) = K           ! time
        COMPEXC(IT,2) = EXCPER      ! excitation
        COMPEXC(IT,3) = MAXXACL     ! mass number
        COMPEXC(IT,4) = Mul_MCFR    ! multi-fragmentation flag

        ! Restored comment. CSM
        if (IT == KKMAX) then
          ! Restored comment. CSM
          EXCCROSS2(MAXXACL) = EXCCROSS2(MAXXACL) + EXCPER*MAXXACL
          EXCTOT = ceiling(EXCPER*MAXXACL/10.)  ! excitated energy distribution
          if (EXCTOT <= ceiling(EXCTOP/10.)) 
     &      EXCCROSS1(EXCTOT) = EXCCROSS1(EXCTOT) + 1
          if (MAXEXC >= EX1DOWN .and. MAXEXC < EX1UP) FLAGEXC = 1
          if (MAXEXC >= EX2DOWN .and. MAXEXC < EX2UP) FLAGEXC = 2
          if (MAXEXC >= EX3DOWN)                      FLAGEXC = 3
          if (FLAGEXC /= 0) then
            do T = KKMIN, KKMAX
              SUMEEHR(FLAGEXC,T,1) = SUMEEHR(FLAGEXC,T,1) + COMPEXC(T,1)
              SUMEEHR(FLAGEXC,T,2) = SUMEEHR(FLAGEXC,T,2) + COMPEXC(T,2)
              SUMEEHR(FLAGEXC,T,3) = SUMEEHR(FLAGEXC,T,3) + COMPEXC(T,3)
              SUMEEHR(FLAGEXC,T,4) = SUMEEHR(FLAGEXC,T,4) + COMPEXC(T,4)
            end do
            NRANGE(FLAGEXC) = NRANGE(FLAGEXC) + 1
          end if
        end if

        !PRTEX_WRTPHAS
        if(PRTEX_WRTPHAS == 1)then
                  do      I=1,MASS
                        XT(K,I)=X(I)
                        YT(K,I)=Y(I)
                        ZT(K,I)=Z(I)
                  enddo

                  if(L1 == 1.and.IT == KKMIN)then
                        PRINT*,'INPUT RANGE(1, 2 OR 3)'
                        read(*,*)      INPUTRANGE
                  endif

                  if(IT == KKMIN)      PRINT*,'EVENT:',L1

                  if(IT == KKMAX.and.INPUTRANGE == FLAGEXC)then
                        PRINT*, 'EVENT: ',L1
                        PRINT*,      'EXC:', MAXEXC
                        PRINT*, '  TIME/fm       EXC/MeV       MAXXACL  
     &                MUL_FRAG'

                        do T = KKMIN,KKMAX
                              PRINT*, COMPEXC(T,1)*DT, COMPEXC(T,2), 
     &                        COMPEXC(T,3),COMPEXC(T,4)
                        enddo

                        if(IT == KKMAX)then
                              PRINT*
                    PRINT*,'INPUT INSTANTS TO write OUT(ENDDED WITH 0)'

2248                              read(*,*) PNTMO
                              if(PNTMO == 0.)goto      2247

                              do I = 1, MASS
                                    if(I == 1)then
                                          write(13,*)
                                          write(13,*)      'EVENT:', L1
                                      write(13,*)      'TIME:', PNTMO
                                    endif

                            write(13,9665)      XT(nint(PNTMO/dt),I),
     &                                          YT(nint(PNTMO/dt),I),
     &                                          ZT(nint(PNTMO/dt),I)
                              enddo
                              goto 2248
                        endif
                  endif
2247                  continue
            endif
      endif

******************************************************************************
****** Restored: **************** 4.2 GEMINI ********************** CSM
    
ccccc  GEMINI 1201 ccccc Restored comment. CSM
      if(GEMINI == 1)then      ! 9992
C018  if (ICL(max_cluster_id) == 232 .and.
C018 &    XJR(max_cluster_id,1) == 90 .and.
C018 &    EXCEN(max_cluster_id) > 0.0) then
C018    task018_extra = 0
C018    do N=1,mass
C018      if(IC(N) == max_cluster_id) cycle
C018      if(IC(N) <= 0 .or. LB(N) /= 1)
C018 &        stop 'TASK018 spectator is not a proton'
C018      task018_extra=task018_extra+1
C018      task018_extra_id=IC(N)
C018      if(ICL(task018_extra_id) /= 1)
C018 &        stop 'TASK018 spectator is not isolated'
C018      CLX(max_cluster_id)=CLX(max_cluster_id)+X(N)
C018      CLY(max_cluster_id)=CLY(max_cluster_id)+Y(N)
C018      CLZ(max_cluster_id)=CLZ(max_cluster_id)+Z(N)
C018      CLKX(max_cluster_id)=CLKX(max_cluster_id)+XK(N)
C018      CLKY(max_cluster_id)=CLKY(max_cluster_id)+YK(N)
C018      CLKZ(max_cluster_id)=CLKZ(max_cluster_id)+ZK(N)
C018      IC(N)=max_cluster_id
C018      ICL(max_cluster_id)=ICL(max_cluster_id)+1
C018      XJR(task018_extra_id,1)=0.
C018      XJR(task018_extra_id,2)=0.
C018      ICL(task018_extra_id)=0
C018    enddo
C018    if(task018_extra /= 1) stop 'TASK018 proton count'
C018  endif

      n_frag = 0
      ICMAX1 = ICMAX
      IMAXMASS = 0
      AMAXFRG = 0
      ZMAXFRG = 0
      NMAXFRG = 0
      EXMAXFRG = 0.
      NHVY = 0

      do 890 I = 1, ICMAX

! Restored: 8.4.2.1----------------------- ------------ CSM

        XACL=XJR(I,1)+XJR(I,2)
        if(XACL <= 0.)goto 890

        IZ= NINT(XJR(I,1)+0.01)
        IA= NINT(XACL + 0.01)
        IN= IA-IZ


        if(IZ <= MCFR)goto 890

        CL2=xmom*xmom*ICL(I)**2

        frx = CLX(I)/ICL(I)                                       !trace
        fry = CLY(I)/ICL(I)
        frz = CLZ(I)/ICL(I)

        em1 = AMU*ICL(I)

        frkx= CLKX(I)
        frky= CLKY(I)
        frkz= CLKZ(I)

        fr2 = frkx*frkx+frky*frky+frkz*frkz

        if(fr2 > CL2)goto 890

        frkk= sqrt(frkx*frkx+frky*frky+frkz*frkz)

        e1 =  sqrt(em1*em1+frkx*frkx+frky*frky+frkz*frkz)

        bx = frkx/e1
        by = frky/e1
        bz = frkz/e1

        ! Restored comment. CSM
        NEXC=EXCEN(I)
        ! Restored comment. CSM
        !if(nexc <= DECNEXC1)goto 890
        !if(nexc > DECNEXC2)nexc= DECNEXC2
        !if(IZ > DECIZ)goto 890            ! Restored: , CSM
                                                                    ! Restored: PRINTDEC CSM

        if(XACL > AMAXFRG)then
          IMAXMASS = I
          AMAXFRG = XACL
          ZMAXFRG = IZ
          NMAXFRG = IN
          EXMAXFRG = NEXC
        endif
        if(ICL(I) > ANGM_LIGHT) NHVY = NHVY + 1

! Restored: 4.2.1-------------------- ----------- CSM




! Restored: 4.2.2----------------------- GEMINI-------------- CSM
      
        do IEN=10,60
          if(IZ == IEN)then
            IEXCI = NINT(0.02*nexc)
            energy(IEXCI,IEN) = energy(IEXCI,IEN) + 1
          endif
        enddo
C             PRINT*, 'I, IZ, IN, NEXC, TOTANGM3'
C             PRINT*, I, IZ, IN, NEXC, TOTANGM3
C             PAUSE
C018  if (I == max_cluster_id .and. IA == 232 .and. IZ == 90
C018 &    .and. EXCEN(I) > 0.0) then
C018    call input(91, IN, INT(task018_pa_ex), -1.)
C018  else
C018    call input(IZ, IN, nexc, -1.)
C018  endif
        call input(IZ, IN, nexc, -1.)
        call prepare_mass
        
        if(para%tl_iwbc)then
          call prepare_Hauser_Feshbach ()
          call read_tl ()
        else
          call prepare_evap ()
        endif

        if(para%I_angle .OR. para%k_sum)call factorial ()
!            if (para%ratio /= 1.0 .and. para%ratio > 0. .and.&
!            .NOT. para%k_SUM) call prepare_deformed_emission (para%ratio,1.)

        if(para%ratio /= 1.0 .and. .NOT. para%k_SUM)
     &    call prepare_deformed_emission (para%ratio,1.)

        if(para%Z_imf_min > 2)call prepare_5li_5he ()

        if(para%imf_option >= 2)call read_sad ()

        if(para%I_angle) call quantum_ang_generator ()

        call alloc_storage

        if(ANGM_DEC == 1)then
          if(ICL(I) > ANGM_LIGHT)then
            L = TOTANGM3
            LL_SHOW = L
          else
            L = 0
          endif
        endif

        if (N_evt <= 0)  N_evt = 0
        call fusion(l,cn)
C CSM: Capture primary state before GEMINI can change cn%J or cn%K.
        A0 = cn%A
        Z0 = cn%Z
        E0 = cn%Ex
        J0 = cn%J
        call event(cn,weight)
C018  if (I == max_cluster_id .and. IA == 232 .and. IZ == 90
C018 &    .and. EXCEN(I) > 0.0) then
C018    write(185,'(4I7,2ES20.10,I7,2I7,3ES20.10)')
C018 &       L1, I, IA, IZ, EXCEN(I),task018_pa_ex,
C018 &       INT(task018_pa_ex),
C018 &       A0, Z0, E0, J0, weight
C018    do nf=1,n_frag
C018      if(frag(nf)%Z < 0 .or. frag(nf)%A <= 0) cycle
C018      write(186,*) L1, frag(nf)%Z, frag(nf)%A, weight
C018    enddo
C018  endif
        
C CSM: Select one largest valid leaf for the MST maximum source only.
        if(gemini_in_out == 1)then
          if(I == max_cluster_id)then
            MAX_FRAG_ID = 0
            MAX_FRAG_A = 0
            do nf=1,n_frag
              if(frag(nf)%Z < 0)cycle
              if(frag(nf)%A <= 0)cycle
              if(MAX_FRAG_ID == 0)then
                MAX_FRAG_ID = nf
                MAX_FRAG_A = frag(nf)%A
              else if(frag(nf)%A > MAX_FRAG_A)then
C CSM: Strict comparison keeps the first equal-mass leaf.
                MAX_FRAG_ID = nf
                MAX_FRAG_A = frag(nf)%A
              endif
            enddo
            if(MAX_FRAG_ID > 0)then
              write(184,1840) Z0,A0,J0,E0,
     &        frag(MAX_FRAG_ID)%Z,frag(MAX_FRAG_ID)%A,
     &        frag(MAX_FRAG_ID)%J,frag(MAX_FRAG_ID)%Ex,
     &        frag(MAX_FRAG_ID)%Id,weight
            endif
          endif
        endif
 1840   format(2I8,2F10.3,2I8,2F10.3,I8,F12.5)

        call de_c

! Restored: 4.2.2----------------------- GEMINI -------------- CSM


! Restored: 4.2.3------------ GEMINI QMD -------------- CSM

c--------------------------------------------------
c----- frag(nf)%Z - proton number of the i-th fragment
c----- frag(nf)%A - nucleon number of the i-th fragment
c----- frag(nf)%Ex - excitation energy of the i-th fragment (MeV)
c----- frag(nf)%j  - spin of the i-th fragment (hbar)
c----- frag(nf)%spin% theta,phi - spin of the i-th fragment (polar angles)
c----- frag(nf)%vel% v,theta,phi - velocity vector of the i-th fragment (magnitude and polar angles) (cm/ns)
c----- frag(nf)%id - branch of decay tree of the i-th fragment
c----- frag(nf)%time, frag(i)%M
c----- n_frag, n_binary
c--------------------------------------------------
            PRINTENE=1
            ledger_source_a = 0
            ledger_source_z = 0
            do nf=1,n_frag
              if(frag(nf)%A <= 0.or.frag(nf)%Z < 0) cycle
              if(frag(nf)%Z > frag(nf)%A)
     &          error stop 'Invalid GEMINI leaf charge'
              ledger_source_a = ledger_source_a + frag(nf)%A
              ledger_source_z = ledger_source_z + frag(nf)%Z
            enddo
            if(ledger_source_a /= A0.or.ledger_source_z /= Z0)then
              write(*,*) 'GEMINI_AZ_MISMATCH',L1,IT,I,A0,Z0,
     &          ledger_source_a,ledger_source_z
              error stop 1
            endif

            do  892  nf=1, n_frag

              NFRA = frag(nf)%A
              NFRZ = frag(nf)%Z               ! =-1, gamma-ray decay
              NFRN = NFRA-NFRZ

              if(NFRA <= 0)goto 892
              if(NFRZ < 0)goto 892

C                   PAUSE

C028 BEGIN APPEND LEAF
              if(ICMAX1 >= MG) error stop 'Cold fragment capacity'
              IS = ICMAX1 + 1
              ICMAX1 = IS
              XJR(IS,3)=0.
C028 END APPEND LEAF
              XJR(IS,1)= NFRZ
              XJR(IS,2)= NFRN
              ICHY(IS) = NFRZ
              NNHY(0,0,IS)=NFRA

              M=0
              do  N=1, mass
                if(SURF_COAL == 1.and.VANISH(N) == 1)cycle
                if(IC(N)  ==  I)  then
                  M = M+1
                  if(M  >  NFRA)   goto  222
                  IC(N) = IS
                endif
              enddo

              ! Restored: , nf v theta/phi pir eaa CSM
 222          VVFR = frag(nf)%vel%v                           ! velocity (cm/ns)
              VTHF = frag(nf)%vel%theta                           ! Restored comment. CSM
              VPHF = frag(nf)%vel%phi

              if(VVFR  <  0.)   VVFR=0.

              pxr = NFRA*VVFR*SIND(VTHF)*COSD(VPHF)*AMU/30.        ! nucleon velocity (cm/ns) --> momentum (GeV/c)
              pyr = NFRA*VVFR*SIND(VTHF)*SIND(VPHF)*AMU/30.
              pzr = NFRA*VVFR*COSD(VTHF)*AMU/30.
              BINFRAG=(XMP*NFRZ+XMN*NFRN-XBIN(NFRZ,NFRA))/1000
              eaa = sqrt(DBLE(AMU1*NFRZ+AMU2*NFRN-BINFRAG)**2 + 
     &            DBLE(pxr)**2+DBLE(pyr)**2 + DBLE(pzr)**2)

              call LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)! pax--momentum of centroid of evaporated particles

 333          XX    = RAN(mseed)
              YY    = RAN(mseed)
              ZZ    = RAN(mseed)

              RR    = SQRT( XX**2 + YY**2 + ZZ**2 )

              ! Restored: IS , RCL CSM
              if((RR  <  0.001) .OR. (RR  >  1.))goto 333

              ICL(IS) = NFRA
                  
              !
              RCL= 1.28*float(NFRA)**0.3333333-0.76+0.8*
     &        float(NFRA)**(-0.3333333)

              CLX(IS)=  frx + XX*RCL                               !trace
              CLY(IS)=  fry + YY*RCL
              CLZ(IS)=  frz + ZZ*RCL

              CLKX(IS)= pax
              CLKY(IS)= pay
              CLKZ(IS)= paz
              clkin(is)=pa0-nfrn*amu2-nfrz*amu1+binfrag      !GeV

              IFSTAT(IS) = 1

              ! Restored: coll.out CSM
              if(PRINTDEC == 1.and.      
     &                XJR(I,1) >= Z11.and.XJR(I,1) <= Z22)then
                if(L11 /= L1)then
                  PRINT*, 'New event!'
                  PRINT*
                  PAUSE
                  L11=L1
                endif
                  
                IAA=XJR(I,1)+XJR(I,2)
                  
                if(PRINTENE == 1)then
                  PRINT*
                  PRINT*,'Monitor the analyzed excited remnants.'
                  PRINT*
                  PRINT*,'EVENT: ',L1
                  PRINT*
                  PRINT*,'REMNANT: ',I,' ','Z = ',int(XJR(I,1)),' ',
     &                    'N = ',int(XJR(I,2))
                      PRINT*
                      PRINT*,'kinetic energy=     tkine=  ',
     &                    tkin(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'coulomb potential=   ucou=  ',
     &                    ucou(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'coulomb exchange=   couex=  ',
     &                    couex(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'yukawa potential=    uyuk= ',
     &                    uyuk11*1000./float(IAA), 'MeV/u'
                      PRINT*,'isospin potential=   usym= ',
     &                    usym(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'surface symmetry=    ussr= ',
     &                    ussr(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'pauli potential=    upaul= ',
     &                    upaul11*1000./float(IAA), 'MeV/u'
                      PRINT*,'mom. dep. potential=umomd= ',
     &                    umod(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'skyrme potential=    usky= ',
     &                    (utwo(I)+uthr(I))*1000./float(IAA), 'MeV/u'
                      PRINT*,'skyrme  (2-body)=    utwo= ',
     &                    utwo(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'skyrme  (3-body)=    uthr= ',
     &                    uthr(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'gradient potential=ugradi= ',
     &                    ugrd(I)*1000./float(IAA),'MeV/u'
                      PRINT*,'effective mass=      utau= ',
     &                    utau(I)*1000./float(IAA), 'MeV/u'
                      PRINT*,'shell correction=     Esh= ',
     &                    Esh11*1000./float(IAA), 'MeV/u'
                      PRINT*
                      PRINT*,'groung state energy= bdgd= ',
     &                    bdgd(i)/float(IAA), 'MeV/u'
                      PRINT*,'average energy=      aene= ',
     &                    aene(i)*1000./float(IAA), 'MeV/u'
                      PRINT*,'excitation energy=    exc= ',
     &                    excen(i)/float(IAA), 'MeV/u'
                      PRINT*
                      PRINT*,'kinetic energy=     clkin= ',
     &                    CLKIN(I)*1000, 'MeV'
                      PRINT*
                      PRINT*,'Decay products:'
                      print*
                      PRINT*,'         ','NF', '          ','NFRZ',
     &                    '        ','NFRN','      ','NFRKIN'
                      PRINTENE=0
            endif
                  
                  PRINT*,NF, NFRZ,NFRN,'      ',1000.*CLKIN(IS), 'MeV'
                  
                        if(NF == N_FRAG)then
                              PRINT*
                              PAUSE
                        endif                                          
                  endif

 892            continue

C028 ICMAX1 is advanced only for accepted massive leaves.

! Restored: 4.2.3------------ GEMINI QMD ---------- CSM

! Restored: 4.2.5----------------- ---------------------- CSM
 
            if(frag(n_frag)%Z > REAL(para%zcn+4)/2.and.para%tl_iwbc) 
     &      then
                  mult0 = mult0 + mult1*weight
                  num_residue = num_residue + weight
            endif

C            enddo cascades

C            enddo loop_ell

            call DE_ARRAY_tl
            call DE_ARRAY_angle

            XJR(I,1)= 0.
            XJR(I,2)= 0.

            CLX(I) = 0.
            CLY(I) = 0.
            CLZ(I) = 0.

            CLKX(I)= 0.
            CLKY(I)= 0.
            CLKZ(I)= 0.
            CLKIN(I)=0.

            ICL(I) = 0.
C018  if (I == max_cluster_id .and. IA == 232 .and. IZ == 90
C018 &    .and. EXCEN(I) > 0.0) then
C018    task018_sum_a=0
C018    task018_sum_z=0
C018    do task018_slot=1,ICMAX1
C018      task018_sum_z=task018_sum_z+NINT(XJR(task018_slot,1))
C018      task018_sum_a=task018_sum_a+NINT(XJR(task018_slot,1)
C018 &                                  +XJR(task018_slot,2))
C018    enddo
C018    if(task018_sum_a /= 233 .or. task018_sum_z /= 91)
C018 &     stop 'TASK018 final event conservation mismatch'
C018  endif

! Restored: 4.2.5----------------- ------------------ CSM

 890  continue



      ICMAX = ICMAX1


! Restored: 4.2-------------------- GEMINI ------------------------- CSM

! Restored: 4----------------------------------------- GEMINI ---------------------------------------- CSM
9992  endif      ! Restored: GEMINI CSM

      ! Restored: , ICMAX , CSM
      if(SURF_COAL == 1)then
            IICMAX = ICMAX
            do MM = 1, 20
                  ICMAX = ICMAX + 1
                  if(SURF_Z(MM) == 0.and.SURF_A(MM) == 0)then
                        ICMAX = ICMAX - 1
                        EXIT
                  else
                        if(ICMAX > MG)
     &                    error stop 'Surface cluster capacity'
                        XJR(ICMAX, 1) = SURF_Z(MM)
                        XJR(ICMAX, 2) = SURF_A(MM) - SURF_Z(MM)
                  endif
                  CLKX(ICMAX) = SURFMOM(MM, 1)
                  CLKY(ICMAX) = SURFMOM(MM, 2)
                  CLKZ(ICMAX) = SURFMOM(MM, 3)
                  CLKIN(ICMAX) = SURFE(MM)*0.001

                  IFSURF(ICMAX) = 1      
            enddo
      endif     

      ledger_cold_z = NINT(SUM(XJR(1:ICMAX,1)))
      ledger_cold_a = ledger_cold_z + NINT(SUM(XJR(1:ICMAX,2)))
      write(189,'(8I12)') L1,IT,ledger_expected_a,ledger_expected_z,
     &  ledger_hot_a,ledger_hot_z,ledger_cold_a,ledger_cold_z
      if(ledger_cold_a /= ledger_hot_a.or.
     &   ledger_cold_z /= ledger_hot_z)then
        write(*,*) 'MST_COLD_MISMATCH',L1,IT,ledger_hot_a,
     &    ledger_hot_z,ledger_cold_a,ledger_cold_z
        error stop 1
      endif

      !PRINTCLKIN.2
      if(PRINTCLKIN == 1.and.L1 >= LEADING1)then
            do      I = 1, ICMAX
                  if(I == 1)  then
                        PRINT*
                        PRINT*,'Cold fragment kinetic energies'
                        PRINT*
                        PRINT*,'         ','Z','           ','N','    ',
     &          'KIN/MeV','       ','AVKIN/MeV','       ','THETA/Deg.'
                  endif

                  IAA=XJR(I,1)+XJR(I,2)
                  if(IAA /= 0.)then
                  PRINT*, INT(XJR(I,1)),INT(XJR(I,2)),CLKIN(I)*1000.,
     &            1000.*CLKIN(I)/IAA,ACOSD(CLKZ(I)/SQRT(CLKX(I)**2+
     &                  CLKY(I)**2+CLKZ(I)**2))
                  endif
            enddo
            PAUSE
      endif

! Restored: 9992 endif ! GEMINI CSM
      
      
C             if (l1 == 116697) pause  !116697

      do IS = 1, ICMAX
          if(ICL(IS)  <=  0.)cycle
          CLX(IS) = CLX(IS)/ICL(IS)                  !trace
          CLY(IS) = CLY(IS)/ICL(IS)
          CLZ(IS) = CLZ(IS)/ICL(IS)
          E_THIS = CLKX(IS)**2 + CLKY(IS)**2 + CLKZ(IS)**2
          E_THIS = E_THIS + (ICL(IS)*0.938)**2
          E_THIS = SQRT(E_THIS)
          call LORENTZ(0., 0., -VPR1, E_THIS, CLKX(IS), CLKY(IS), 
     &    CLKZ(IS), E_THIS, CLKX(IS), CLKY(IS), CLKZ(IS))
          CLKX(IS)=CLKX(IS)/ICL(IS)                  !There is in fact an over-dividing after which CLKX(i)
          CLKY(IS)=CLKY(IS)/ICL(IS)                  !will be the per-nucleon momentum
          CLKZ(IS)=CLKZ(IS)/ICL(IS)
      enddo

      do      I = 1, ICMAX
            if(ICL(I) <= 0)      cycle
      enddo

******************************************************************************
******************************************************************************
******************************************************************************
      


! Restored: CNPCS.4 A N Z -------------------------------------------------------------------------------- CSM
      if(CNPCS. EQ. 1) then
COLD FRAGMENTS (ISOTC:isotope cold)
        do IIC=1,ICMAX
          if ((XJR(IIC,1) /= 0).OR.(XJR(IIC,2) /= 0)) then
            NZISOTC(XJR(IIC,1)) = NZISOTC(XJR(IIC,1)) + 1   ! Restored: XJR(IIC,1) CSM
            NNISOTC(XJR(IIC,2)) = NNISOTC(XJR(IIC,2)) + 1   ! Restored: XJR(IIC,2) CSM
            NAISOTC(XJR(IIC,1)+XJR(IIC,2)) =                     ! Restored: XJR(IIC,1)+XJR(IIC,2) CSM
     &      NAISOTC(XJR(IIC,1)+XJR(IIC,2)) + 1                
          endif
          ! Cross sections of isotops
          do I = 1, CSSPECIES
            if ((XJR(IIC,1) == CNPCS_SETINF(I,1)).and.
     &          (XJR(IIC,2) == CNPCS_SETINF(I,2)))
     &        ele_count(I) = ele_count(I) + 1
          enddo      
        enddo
      endif
! Restored: CNPCS.4 A N Z ---------------------------------------------------------------------------------------------- CSM



! Restored: IPCS.3 -------------------------------------------------------------------------------------------------- CSM
      if (IPCS == 1) then
        do I = 1, ICMAX
          do J = 1, NIPCS
            if (XJR(I,1) == IPCS_Z(J).and.XJR(I,2) <= 200) 
     &        IPCS_COUNT(J,XJR(I,2))=IPCS_COUNT(J,XJR(I,2))+1
          enddo
        enddo
      endif
! Restored: IPCS.3 -------------------------------------------------------------------------------------------------------------- CSM

! Restored: TEELECS.4 ----------------------------------------------------------------------------------------------------------------- CSM
      if(TEELECS == 1.and.K >= CEILING(TECOLD/DELT))then
            do      IIC=1, ICMAX
                  do I=1, TEELESPECIES
                        if(TEELECS_SETINF(I,1) == XJR(IIC,1).and.
     &                TEELECS_SETINF(I,2) == XJR(IIC,2))then
                              TEELECSCOUNT(k,I,2)=TEELECSCOUNT(k,I,2)+1
                        endif
                  enddo
            enddo
      endif
! Restored: TEELECS.4 ------------------------------------------------------------------------------------------------------------- CSM

      ! Restored: beta CSM
      BETA_T = 2*0.938*NPR*EPP*NPR + (EPP*NPR)**2
      BETA_T = -SQRT(BETA_T)/(0.938*(NPR+NTA) + NPR*EPP)
      GAMM_T = 1./SQRT(1 - BETA_T**2)

! Restored: 8.5------------------------------ ---------------------------------------- CSM
      do 894 I=1, ICMAX
      
              XACL=XJR(I,1)+XJR(I,2)
            if(XACL  <=  0.)  goto 894

            NBH= NINT(XJR(I,3)+0.1)
            NIA= NINT(XACL+0.1)
            NFR= NINT(XJR(I,1)+0.1)
            NCH= ICHY(I)
            if(NCH  <  0)  goto 894

            do IE=2,80
                  if(XJR(I,1) == IE)then
                        do IES=2,200
                              if(XACL == IES)then
                                    AISOT(IE,IES)=AISOT(IE,IES)+1
                              endif
                          enddo
                    endif
            enddo


! Restored: LNVICS.4 ------------------------------------------------------------------------------------------------------------------------- CSM
            if(LNVICS == 1. AND. ABS(CLKZ(I)*ICL(I)) <= 1.)then
            !4He
                  if(XJR(I,1) == 1. AND. XJR(I,2) == 2)then
                  
                  endif
            endif
! Restored: LNVICS.4 --------------------------------------------------------------------------------------------------------------------- CSM

            ppin=sqrt(CLKX(I)**2 + CLKY(I)**2 + CLKZ(I)**2)
            emcl=sqrt(AMU**2+CLKX(I)**2+CLKY(I)**2+CLKZ(I)**2)*ICL(I)
            trcl=sqrt(CLKX(I)**2+CLKY(I)**2)
            velo=29.98*ppin/AMU
            vpar=29.98*CLKZ(I)/AMU
      !      write(229,*)I,vpar
            if(trcl == 0.) then
            fldr=0.
            flep=0.
            else
            fldr=CLKX(I)/trcl
            flep=(CLKX(I)**2-CLKY(I)**2)/trcl**2
            endif

! Restored: DDCSS_LNKED.3 --------------------------------------------------------------------------------------------------- CSM
            if(DDCSS_LNKED. EQ. 1.and.
     &            (IFSTAT(I) == 1.OR.IFSURF(I) == 1))then
                  
              mnkin=ceiling(mkin/dkin)                  ! Restored: bin CSM
                  mnthe=ceiling(pi/dthe)                        ! Restored: bin CSM
                  if(DIFFUSION(I) >= 0.and.DIFFUSION(I) <= DIFF_T)then
                        call RANDOM_NUMBER(A_DIFF)
                        C_DIFF = 1.-EXP(B_DIFF/(DIFFUSION(I)-DIFF_T
     &                -0.000001))
                        if(A_DIFF > C_DIFF)goto 5201
                  endif

                  do      SPECIES=1,NSPECIES
                  Ekin=0.
                  if(XJR(I,1) ==  DDCSS_LNKED_SETINF(SPECIES,1).and.             
     &                XJR(I,2) ==  DDCSS_LNKED_SETINF(SPECIES,2))then           
                      if(SURF_COAL == 1.and.IF_CASCADE(I) == 1)then
                                  if(EMD_CL(I) == 0)cycle ! Restored: cascade , CSM
                            endif             
                              ! Restored: DDCSS CSM
                              Ekin = CLkin(i)
                          costhet = clkz(i)/sqrt(CLKX(I)**2+CLKY(I)**2
     &                            +CLKZ(I)**2)
                              thet = acos(costhet)+0.000000000001
                              nkin = ceiling(ekin*1000./dkin)                  ! Restored: Ekin [ceiling(ekin/dkin)-1, ceiling(ekin/dkin)) CSM
                              nthe = ceiling(thet/dthe)                            ! Restored comment. CSM
                      if(nkin > 0.and.nkin <= mnkin.and.
     &                    nthe > 0.and.nthe < mnthe)then
                          cross_count(SPECIES,nkin,nthe)=
     &                    cross_count(SPECIES,nkin,nthe)+1
                      endif
                              ! Restored: LNKED CSM
                              nekin=ceiling(1000.*clkin(i)/dkin)
                          if(nekin  >  0 .and. nekin  <=  mnkin)then
                          LNKEB(SPECIES,NEKIN)=
     &                                LNKEB(SPECIES,NEKIN)+1
                      endif
                        endif      
                  enddo
5201                  continue            
            endif
! Restored: DDCSS_LNKED.3 ----------------------------------------------------------------------------------------------- CSM

! Restored: ELEKED.3 -------------------------------------------------------------------------------------------------------------------- CSM
            if(ELEKED == 1)then
              if(DIFFUSION(I) >= 0. .and. DIFFUSION(I) <= DIFF_T)then
                        call RANDOM_NUMBER(A_DIFF)
                        C_DIFF = 1. - EXP(B_DIFF/(DIFFUSION(I)-DIFF_T-
     &                0.000001))
                        if(A_DIFF > C_DIFF)goto 52011
                  endif
                  do      J=1,ELESPECIES
                  ! Restored comment. CSM
        if(XJR(I,1) == ELEKED_SETINF(J,1).and.
     &   XJR(I,2) == ELEKED_SETINF(J,2).and.
     &  CLKIN(I)*1000 >= 0..and.CLKIN(I)*1000 <= 100.)      then
                      ELEKEDCOUNT(J,NINT(CLKIN(I)*1000))=ELEKEDCOUNT ! Restored: NINT CSM
     &                        (J,NINT(CLKIN(I)*1000))+1
                        endif
                  enddo
52011                  continue
            endif
! Restored: ELEKED.3 ---------------------------------------------------------------------------------------------------------------- CSM



 894      continue

! Restored: 8.5------------------------------ ------------------------------------- CSM

2500      endif      ! Restored: endif CSM


9999      endif      ! Restored: endif CSM
ccccc  9999 ccccc Restored comment. CSM




9668      continue     ! Restored comment. CSM

      if(NTH_RUN /= 0) write(OUTFILE+5,*)  NNEVENT ! Restored: NEVENT*kktot CSM

9667      if(L1 >= NRMAX.and.(NTH_RUN /= 0))      EXIT

9666       continue     ! Restored comment. CSM

! Restored: 8---------------------------------------------------------------- ------------------------------------------------------------------- CSM











c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c
c
c      Restored comment. CSM
c
c
c
c////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! Restored: 9---------------------------------------------------------------- ----------------------------------------------------------------------- CSM
      if(NTH_RUN == 0) then
            print*, 'Total analyzed events:', NEVENT
      else
            print*, 'Total analyzed events:', NNevent
      endif

      factor=float(NEVENT)*float(kkmax-kkmin+1)
      th1=(theta0-dtheta)/180.*PI
      th2=(theta0+dtheta)/180.*PI
      Domiga=2.*PI*(cos(th1)-cos(th2))


! Restored: DDCSS_LNKED --------------------------------------------------------------------------------------------------------------------------------- CSM
      if(DDCSS_LNKED == 1.and.DDCSS_LNKED_S /= 1)then
            !DDCSS
9815            format(1x, '\Ekin(MeV)', 4x, 'diffcrsec(mb/srMeV)')

            write(OUTFILE+1,*)
            write(OUTFILE+1,*) 'Double Differential Cross Sections'
            write(OUTFILE+1,*)


            do      SPECIES=1,NSPECIES

                  do      II=1,NANGLE

                        if(II == 1 .and. NTH_RUN == 0)      then
                      write(OUTFILE+1,*) '-----------------------------
     &-------'
                      write(OUTFILE+1,*)  DDCSS_LNKED_LABEL(SPECIES)
                      write(OUTFILE+1,*) '-----------------------------
     &-------'
                              write(OUTFILE+1,*)
                        endif

                        if(II == 1 .and. NTH_RUN /= 0)      
     &                  write(OUTFILE+1,*)  DDCSS_LNKED_LABEL(SPECIES)

                        thetad=DDCSS_LNKED_ANGLE(II)    ! Restored comment. CSM
                        if(NTH_RUN == 0)      then
                              write(OUTFILE+1,*), '\theta= ',thetad 
                        else
                              write(OUTFILE+1,*)      thetad
                        endif
                  if(NTH_RUN == 0) write(OUTFILE+1, 9815)
                        
                  theta=thetad*pi/180     ! sr
                        
                  do i=1,ceiling(mkin/dkin)
                      if(cross_count(species,i,NINT(theta/dthe)) == 0.
     &                        AND.NTH_RUN == 0) cycle
                              
                      diffcrsec=cross_count(species,i,NINT(theta/dthe))
     &                        *pi*X00CUT**2/(2*pi*sin(theta)*dkin*dthe)
                              diffcrsec=diffcrsec/NEVENT
                    if(NTH_RUN /= 0)      diffcrsec = diffcrsec*NEVENT
                              write(OUTFILE+1,*) i*dkin, 
     &                        DDCSS_LNKED_SPACFAC(II)*diffcrsec*10
                        enddo
                        if(NTH_RUN == 0) write(OUTFILE+1,*)

                  enddo

                  if(NTH_RUN == 0) write(OUTFILE+1,*)

            enddo



            !LNKED
            write(OUTFILE+1,*) 'Light Nuclide Kinetic Energy 
     &        Distribution'
            if(NTH_RUN == 0)then
              write(OUTFILE+1,*)      'Units:',' ','MeV','   ','mb/MeV'
                  write(OUTFILE+1,*)
            endif

            do      SPECIES=1,NSPECIES

                  if(NTH_RUN == 0)then
                        write(OUTFILE+1,*)
                  write(OUTFILE+1,*)'-----------------------------------
     &---'
                        write(OUTFILE+1,*) DDCSS_LNKED_LABEL(SPECIES)
                write(OUTFILE+1,*)'-----------------------------------
     &---'
                        write(OUTFILE+1,*)
                  endif

                  if(NTH_RUN /= 0)then
                        write(OUTFILE+1,*) DDCSS_LNKED_LABEL(SPECIES)
                  endif

                  do NEKIN1 = 1, ceiling(mkin/dkin)
              if(LNKEB(SPECIES,NEKIN1) == 0.and.NTH_RUN == 0)cycle
                  lncross = 10.*LNKEB(SPECIES,NEKIN1)
     &                      *pi*X00CUT**2/(NEVENT*DKIN)                  
                        if(NTH_RUN /= 0)lncross = lncross*NEVENT
                  write(OUTFILE+1,*)      NEKIN1*DKIN, lncross
                  enddo
                  if(NTH_RUN == 0) write(OUTFILE+1,*)

            enddo
      endif

! Restored: DDCSS_LNKED ----------------------------------------------------------------------------------------------------------------- CSM
      if(DDCSS_LNKED == 1.and.DDCSS_LNKED_S == 1)then

        KINORD1=0
        ORDMAG=0
        KINORD=0.

        do K=1,CROSSMAX
          do I=1,ceiling(pi/dthe)
            do J=1,ceiling(mkin/dkin)
              ORDMAG=0
              KINORD=J*DKIN
              if(J*DKIN > 10.and.J*DKIN <= 100)then
                ORDMAG=1
                KINORD=J*DKIN/10.
              endif
              if(J*DKIN > 100.and.J*DKIN <= 1000)then
                ORDMAG=2
                KINORD=J*DKIN/100.
              endif
              if(J*DKIN > 1000.and.J*DKIN <= 10000)then
                ORDMAG=3
                KINORD=J*DKIN/1000.
              endif
              KINORD1=CEILING(KINORD/DKIN)
              cross_count_s(K,ORDMAG,KINORD1,I)=
     &          cross_count_s(K,ORDMAG,KINORD1,I)+
     &          CROSS_COUNT(K,J,I)      
            enddo
          enddo
        enddo

        !DDCSS
9816    format(1x, '\Ekin(MeV)', 4x, 'diffcrsec(mb/srMeV)')
        write(OUTFILE+1,*)
        write(OUTFILE+1,*) 'Double Differential Cross Sections'
        write(OUTFILE+1,*)


        DDCROSSFAC=0.
        do SPECIES=1,NSPECIES
          do II=1,NANGLE
            if(II == 1)then
              write(OUTFILE+1,*) '---------------------------
     &          ---------'
              write(OUTFILE+1,*)  DDCSS_LNKED_LABEL(SPECIES)
              write(OUTFILE+1,*) '---------------------------
     &          ---------'
              write(OUTFILE+1,*)
            endif

            thetad=DDCSS_LNKED_ANGLE(II)
            write(OUTFILE+1,*), '\theta= ',thetad 
            theta=thetad*pi/180
            write(OUTFILE+1, 9816)

            do ORM=0,3
              do i=1,ceiling(10/dkin)
                if(cross_count_s(species,orm,i, 
     &            NINT(theta/dthe)) == 0) cycle
                DDCROSSFAC=2*pi*sin(theta)*dthe*(dkin*10**ORM)
                diffcrsec=pi*X00CUT**2*cross_count_s(species,orm,i,
     &            NINT(theta/dthe))/DDCROSSFAC
                diffcrsec=diffcrsec/NEVENT
                write(OUTFILE+1,*) i*dkin*10**ORM, 
     &            DDCSS_LNKED_SPACFAC(II)*diffcrsec*10
              enddo
            enddo
            write(OUTFILE+1,*)
          enddo
          write(OUTFILE+1,*)
        enddo
      endif

! Restored: ELEKED -------------------------------------------------------------------------------------------------------------------------- CSM
      if(ELEKED == 1.and.NTH_RUN == 1)      then
      ! Restored: OUTPUTFILE+2 21 CSM
            write(21,*) 'Element Kinetic Energy Distribution'
            write(21,*)

            do      I=1,ELESPECIES

                  write(21,*) ELEKED_LABEL(I)

                  do      J=0,100
                        if(ELEKEDCOUNT(I,J) == 0)      cycle
                        write(21,*)      J, 10.*ELEKEDCOUNT(I,J)*
     &                  PI*X00CUT**2/NEVENT
c                  write(21,*) J, ELEKEDCOUNT(I,J)/NEVENT
              enddo
                  write(21,*)

            enddo

      endif 
          
! Restored: TEELECS --------------------------------------------------------------------------------------------------------------- CSM
      if(TEELECS == 1)then
          !write(*,*) 'DEBUG: === ENTERING TEELECS OUTPUT ==='
            do      I = 1, 2
                  write(OUTFILE+6,*)
                  write(OUTFILE+6,*) '-------------------------'
                  if(I == 1)then
                  write(OUTFILE+6,*) 'RESIDUES' ! Restored comment. CSM
              else
                  write(OUTFILE+6,*) 'FRAGMENTS'! Restored comment. CSM
                  endif
                  write(OUTFILE+6,*) '-------------------------'
                  write(OUTFILE+6,*)

                  do      J = 1, TEELESPECIES
                        write(OUTFILE+6,*) TEELECS_LABEL(J)
                        do      K = 1, CEILING(TMAX/DELT)
                              COUNTS = 0
                              COUNTS = TEELECSCOUNT(K,J,I)
                      if(NTH_RUN == 0.and.COUNTS == 0)cycle
                      if(NTH_RUN /= 0) COUNTS = COUNTS*NEVENT
                      ! Restored: (mb) = * / *10 cm^2 mb CSM
!                              if(COUNTS /= 0) write(OUTFILE+6,*) k*DELT, COUNTS*10*PI*
!     &                            X00CUT**2/NEVENT 
                      ! Restored comment. CSM
                      write(OUTFILE+6,*) k*DELT, COUNTS/FLOAT(NEVENT)
                        enddo
                        write(OUTFILE+6,*)
                  enddo
            enddo
      endif

! Restored: FSVT.5 pi-/pi+ --------------------------------------------------------------------------------------------- CSM
      if(FSVT. EQ. 1) then
! Restored: FSVT.5.1 pi-/pi+ ----------------------------------------------------------- CSM
2891  FORMAT(2X, '\FSVT.5.1: Center-of-mass kinetic energy distribution
     &of pi-/pi+')
      write(3,2891)

      do M=1,50
            if(NNKEPION(1,M) <= 0) cycle
            RPION=REAL(NNKEPION(-1,M))/NNKEPION(1,M)
            
          SQNOM=REAL(SNNKEPION(1,M))/NEVENT-
     &          (REAL(NNKEPION(1,M))/NEVENT)**2
            if(SQNOM <= 0.) cycle
            APPION=1.96*SQRT(SQNOM)/SQRT(REAL(NEVENT))

            SQDEN=REAL(SNNKEPION(-1,M))/NEVENT-
     &          (REAL(NNKEPION(-1,M))/NEVENT)**2
            if(SQDEN <= 0.) cycle
            ANPION=1.96*SQRT(SQDEN)/SQRT(REAL(NEVENT))

            BPION=SQRT((APPION*NEVENT/NNKEPION(1,M))**2+(ANPION*NEVENT/
     &          NNKEPION(-1,M))**2)
            
            write(3,*), M*1000*RMKECM, RPION, RPION*BPION
      enddo

      write(3,*) ' '
      write(3,*) ' '
      write(3,*) ' '

! Restored: FSVT.5.2 pi-/pi+ ----------------------------------------------------------------- CSM
2895  FORMAT(2X, '\FSVT.5.2: Transverse momentum distribution of pi-/pi+
     &')
      write(3,2895)

      do M=1,50
            if(NNTRPION(1,M) <= 0) cycle
            RPION=REAL(NNTRPION(-1,M))/NNTRPION(1,M)
            
          SQNOM=REAL(SNNTRPION(1,M))/NEVENT-
     &          (REAL(NNTRPION(1,M))/NEVENT)**2
            if(SQNOM <= 0.) cycle
            APPION=1.96*SQRT(SQNOM)/SQRT(REAL(NEVENT))

            SQDEN=REAL(SNNTRPION(-1,M))/NEVENT-
     &          (REAL(NNTRPION(-1,M))/NEVENT)**2
            if(SQDEN <= 0.) cycle
            ANPION=1.96*SQRT(SQDEN)/SQRT(REAL(NEVENT))

            BPION=SQRT((APPION*NEVENT/NNTRPION(1,M))**2+(ANPION*NEVENT/
     &          NNTRPION(-1,M))**2)
            
            write(3,*) M*RRPT/EPP, RPION, RPION*BPION
      enddo

      write(3,*) ' '
      write(3,*) ' '
      write(3,*) ' '

! Restored: FSVT.5.3 p-/p+ y/yproj ------------------------------------------------------------------ CSM
2896      FORMAT(2X, '\FSVT.5.3: y/yproj distribution of pi-/pi+')
      write(3,2896)

      do M=1,50
            if(NNYPPION(1,M) <= 0) cycle
            RPION=REAL(NNYPPION(-1,M))/NNYPPION(1,M)
            
          SQNOM=REAL(SNNYPPION(1,M))/NEVENT-
     &          (REAL(NNYPPION(1,M))/NEVENT)**2
            if(SQNOM <= 0.) cycle
            APPION=1.96*SQRT(SQNOM)/SQRT(REAL(NEVENT))

            SQDEN=REAL(SNNYPPION(-1,M))/NEVENT-
     &          (REAL(NNYPPION(-1,M))/NEVENT)**2
            if(SQDEN <= 0.) cycle
            ANPION=1.96*SQRT(SQDEN)/SQRT(REAL(NEVENT))

            BPION=SQRT((APPION*NEVENT/NNYPPION(1,M))**2+(ANPION*NEVENT/
     &          NNYPPION(-1,M))**2)
            
            write(3,*), RYDRM*M, RPION, RPION*BPION
      enddo      

      write(3,*) ' '
      write(3,*) ' '
      write(3,*) ' '
      
! Restored: FSVT.5.4 dM/dKEc.m. KEc.m. ------------------------------------------------------------------ CSM
2897  FORMAT(2X,'\FSVT.5.4: KEc.m. distribution of dM/dKEc.m. for pi-+')
      write(3,2897)
      
      write(3,*) 'pi-'
      do M=1,50
            write(3,*), M*1000*RMKECM, DNKEPION(-1,M)/NEVENT
      enddo

      write(3,*) 'pi+'
      do M=1,50
            write(3,*), M*1000*RMKECM, DNKEPION(1,M)/NEVENT
      enddo

      write(3,*) ' '
      write(3,*) ' '
      write(3,*) ' '


      endif
! Restored: FSVT.5 pi-/pi+ ----------------------------------------------------------------------------------------- CSM



! Restored: CNPCS.5 Z N A ---------------------------------------------------------------------------------------------------------------- CSM
      if(CNPCS == 1)then

        !HOT FRAGMENTS
        write(OUTFILE+3,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF Z'
        do IIZ=1,150
          if(NTH_RUN == 0.and.NZISOTH(IIZ) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+3,*) IIZ, PI*X00CUT**2*NZISOTH(IIZ)*10/NEVENT
          else
            write(OUTFILE+3,*) IIZ, PI*X00CUT**2*NZISOTH(IIZ)*10
          endif
        enddo

        write(OUTFILE+3,*) ' '
        write(OUTFILE+3,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF N'
        do IIN=1,150
          if(NTH_RUN == 0.and.NNISOTH(IIN) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+3,*) IIN, PI*X00CUT**2*NNISOTH(IIN)*10/NEVENT
          else
            write(OUTFILE+3,*) IIN, PI*X00CUT**2*NNISOTH(IIN)*10
          endif
        enddo

        write(OUTFILE+3,*) ' '
        write(OUTFILE+3,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF A'
        do IIA=1,300
          if(NTH_RUN == 0.and.NAISOTH(IIA) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+3,*) IIA, PI*X00CUT**2*NAISOTH(IIA)*10/NEVENT
          else
            write(OUTFILE+3,*) IIA, PI*X00CUT**2*NAISOTH(IIA)*10
          endif
        enddo

        !COLD FRAGMENTS
        write(OUTFILE+4,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF Z'
        do IIZ=1,150
          if(NTH_RUN == 0.and.NZISOTC(IIZ) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+4,*) IIZ, PI*X00CUT**2*NZISOTC(IIZ)*10/NEVENT
          else
            write(OUTFILE+4,*) IIZ, PI*X00CUT**2*NZISOTC(IIZ)*10
          endif
        enddo

        write(OUTFILE+4,*) ' '
        write(OUTFILE+4,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF N'
        do IIN=1,150
          if(NTH_RUN == 0.and.NNISOTC(IIN) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+4,*) IIN, PI*X00CUT**2*NNISOTC(IIN)*10/NEVENT
          else
            write(OUTFILE+4,*) IIN, PI*X00CUT**2*NNISOTC(IIN)*10
          endif
        enddo

        write(OUTFILE+4,*) ' '
        write(OUTFILE+4,*) 'NUCLIDE CROSS SECTION AS A FUNCTION OF A'
        do IIA=1,300
          if(NTH_RUN == 0.and.NAISOTC(IIA) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+4,*) IIA, PI*X00CUT**2*NAISOTC(IIA)*10/NEVENT
          else
            write(OUTFILE+4,*) IIA, PI*X00CUT**2*NAISOTC(IIA)*10
          endif
        enddo
      
        write(OUTFILE+4,*) ' '
        write(OUTFILE+4,*) 'CROSS SECTION OF ISOTOPICS'
        do IIP=1,crossmax
          if(NTH_RUN == 0.and.ele_count(IIP) == 0)cycle
          if(NTH_RUN == 0)then
            write(OUTFILE+4,*) IIP,PI*X00CUT**2*ele_count(IIP)*10/NEVENT
          else
            write(OUTFILE+4,*) IIP, PI*X00CUT**2*ele_count(IIP)*10
          endif
        enddo
      endif
! Restored: CNPCS.5 Z N A ------------------------------------------------------------------------------------------------------------ CSM


! Restored: TEAE.5 ------------------------------------------------------------------------------------------------------- CSM
      if (TEAE == 1) then
        if (TEAMUL == 1) then
          do I = 1, 3
            if (NRANGE(I) == 0) cycle
            write (6,*) 'RANGE:',I
            write (6,*) 'EVENTS:', NRANGE(I)
            write (6,*) '  TIME/fm       EXC/MeV       MAXXACL      
     &        MUL_FRAG'
            do T = KKMIN, KKMAX
              write (6,*) SUMEEHR(I,T,1)*DT/NRANGE(I),
     &                    SUMEEHR(I,T,2)/NRANGE(I),
     &                    SUMEEHR(I,T,3)/NRANGE(I),
     &                    SUMEEHR(I,T,4)/NRANGE(I)
            end do
            write (6,*)
          end do
        endif
        if (EXCCROSS == 1) then
          write (OUTFILE+7,*) 'EXCITATION ENERGY CROSS SECTION'
          write (OUTFILE+7,*)
          write (OUTFILE+7,*) 'MeV      mb/MeV'
          do I = 0, ceiling(EXCTOP/10.)
            if (NTH_RUN == 0) then
              if (EXCCROSS1(I) == 0) cycle
              write (OUTFILE+7,*) I*10, 10.*EXCCROSS1(I)*PI*X00CUT**2
     &          /(NEVENT*10.)
            else
              write (OUTFILE+7,*) I*10, 10.*EXCCROSS1(I)*PI*X00CUT**2
     &          /(10.)              
            end if
            !write (OUTFILE+7,*) I*10, EXCCROSS(I)/(NEVENT*10.)
          end do
          write (OUTFILE+7,*) 'EXCITATION ENERGY AS FUNCTION OF MASS 
     &      NUMBER'
          write (OUTFILE+7,*)
          write (OUTFILE+7,*) 'A      mb'
          do I = 1, 300
            if (NTH_RUN == 0) then
              if (EXCCROSS2(I) == 0) cycle
              write (OUTFILE+7,*) I, EXCCROSS2(I)/float(NEVENT)
            else
              write (OUTFILE+7,*) I, EXCCROSS2(I)
            end if
          end do
        end if
      end if
! Restored: TEAE.5 --------------------------------------------------------------------------------------------------- CSM





! Restored: IPCS.4 ------------------------------------------------------------------------------------------------------------------------------ CSM
      if(IPCS == 1)      then

            write(9,*)  'Isotope Production Cross Sections'
            write(9,*)      'Unit: mb'

            do      I=1,NIPCS

                  write(9,*)
              write(9,*)      '-------------------------------------'

                  write(9,*)      IPCS_LABEL(I)

                  do      IIN=0,200
                  if(IPCS_COUNT(I,IIN) /= 0) write(9,*),      IIN, 
     &                  IPCS_COUNT(I,IIN)*10*PI*X00CUT**2/NEVENT
                  enddo

            enddo

      endif
! Restored: IPCS.4 -------------------------------------------------------------------------------------------------------------------------- CSM


! Restored: LNVICS.5 --------------------------------------------------------------------------------------------------------------------------- CSM
      if(LNVICS == 1)      then
            
      endif
! Restored: LNVICS.5 ----------------------------------------------------------------------------------------------------------------------- CSM


      if(gemini_in_out == 1)then
        close(182)
        close(184)
      endif
      if(ECHECK == 1) close(183)

52013      continue
      PRINT*, 'OK, the program ends normally !!'
c      PAUSE

      END



      

*******************************************************************************
            FUNCTION ISOTAB(IZ,IN)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Table of Isotopes, it is extracted from NUCLEAR WALLET CARDS 
cccc (Fifth edition), JULY 1995, by JAGDISH K. TULI, NATIONAL NUCLEAR DATA 
cccc CENTER for The U. S. Nuclear Data Network, Supported by The Division of 
cccc Nuclear Physics, Office of High Energy and Nuclear Physics, US Department 
cccc of Energy.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc      if isotope (IZ,IN) exists, then ISOTAB=1, or ISOTAB=0      ccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       common /RAMISO/ Iso,Iram      

      if(iso == 0) then
      isotab=1
      return
      endif

      if(IZ == 0.and.IN == 1) then
      ISOTAB=1
      else if(IZ == 1.and.IN >= 0.and.IN <= 5) then
      ISOTAB=1
      else if(IZ == 2.and.IN >= 1.and.IN <= 8) then
      ISOTAB=1
      else if(IZ == 3.and.IN >= 1.and.IN <= 8) then
      ISOTAB=1
      else if(IZ == 4.and.IN >= 2.and.IN <= 10) then
      ISOTAB=1
      else if(IZ == 5.and.IN >= 2.and.IN <= 14) then
      ISOTAB=1
      else if(IZ == 6.and.IN >= 2.and.IN <= 16) then
      ISOTAB=1
      else if(IZ == 7.and.IN >= 3.and.IN <= 17) then
      ISOTAB=1
      else if(IZ == 8.and.IN >= 4.and.IN <= 18) then
      ISOTAB=1
      else if(IZ == 9.and.IN >= 5.and.IN <= 20) then
      ISOTAB=1
      else if(IZ == 10.and.IN >= 5.and.IN <= 22) then
      ISOTAB=1
      else if(IZ == 11.and.IN >= 6.and.IN <= 24) then
      ISOTAB=1
      else if(IZ == 12.and.IN >= 7.and.IN <= 24) then
      ISOTAB=1
      else if(IZ == 13.and.IN >= 8.and.IN <= 26) then
      ISOTAB=1
      else if(IZ == 14.and.IN >= 8.and.IN <= 28) then
      ISOTAB=1
      else if(IZ == 15.and.IN >= 9.and.IN <= 31) then
      ISOTAB=1
      else if(IZ == 16.and.IN >= 10.and.IN <= 32) then
      ISOTAB=1
      else if(IZ == 17.and.IN >= 11.and.IN <= 34) then
      ISOTAB=1
      else if(IZ == 18.and.IN >= 12.and.IN <= 33) then
      ISOTAB=1
      else if(IZ == 19.and.IN >= 13.and.IN <= 35) then
      ISOTAB=1
      else if(IZ == 20.and.IN >= 14.and.IN <= 34) then
      ISOTAB=1
      else if(IZ == 21.and.IN >= 15.and.IN <= 34) then
      ISOTAB=1
      else if(IZ == 22.and.IN >= 16.and.IN <= 36) then
      ISOTAB=1
      else if(IZ == 23.and.IN >= 17.and.IN <= 38) then
      ISOTAB=1
      else if(IZ == 24.and.IN >= 18.and.IN <= 40) then
      ISOTAB=1
      else if(IZ == 25.and.IN >= 19.and.IN <= 41) then
      ISOTAB=1
      else if(IZ == 26.and.IN >= 19.and.IN <= 43) then
      ISOTAB=1
      else if(IZ == 27.and.IN >= 21.and.IN <= 45) then
      ISOTAB=1
      else if(IZ == 28.and.IN >= 22.and.IN <= 50) then
      ISOTAB=1
      else if(IZ == 29.and.IN >= 23.and.IN <= 50) then
      ISOTAB=1
      else if(IZ == 30.and.IN >= 24.and.IN <= 52) then
      ISOTAB=1
      else if(IZ == 31.and.IN >= 25.and.IN <= 53) then
      ISOTAB=1
      else if(IZ == 32.and.IN >= 26.and.IN <= 54) then
      ISOTAB=1
      else if(IZ == 33.and.IN >= 27.and.IN <= 56) then
      ISOTAB=1
      else if(IZ == 34.and.IN >= 31.and.IN <= 57) then
      ISOTAB=1
      else if(IZ == 35.and.IN >= 33.and.IN <= 59) then
      ISOTAB=1
      else if(IZ == 36.and.IN >= 33.and.IN <= 61) then
      ISOTAB=1
      else if(IZ == 37.and.IN >= 35.and.IN <= 65) then
      ISOTAB=1
      else if(IZ == 38.and.IN >= 36.and.IN <= 64) then
      ISOTAB=1
      else if(IZ == 39.and.IN >= 38.and.IN <= 66) then
      ISOTAB=1
      else if(IZ == 40.and.IN >= 40.and.IN <= 67) then
      ISOTAB=1
      else if(IZ == 41.and.IN >= 41.and.IN <= 69) then
      ISOTAB=1
      else if(IZ == 42.and.IN >= 42.and.IN <= 71) then
      ISOTAB=1
      else if(IZ == 43.and.IN >= 43.and.IN <= 72) then
      ISOTAB=1
      else if(IZ == 44.and.IN >= 43.and.IN <= 74) then
      ISOTAB=1
      else if(IZ == 45.and.IN >= 44.and.IN <= 76) then
      ISOTAB=1
      else if(IZ == 46.and.IN >= 45.and.IN <= 77) then
      ISOTAB=1
      else if(IZ == 47.and.IN >= 47.and.IN <= 80) then
      ISOTAB=1
      else if(IZ == 48.and.IN >= 49.and.IN <= 82) then
      ISOTAB=1
      else if(IZ == 49.and.IN >= 49.and.IN <= 84) then
      ISOTAB=1
      else if(IZ == 50.and.IN >= 50.and.IN <= 87) then
      ISOTAB=1
      else if(IZ == 51.and.IN >= 52.and.IN <= 88) then
      ISOTAB=1
      else if(IZ == 52.and.IN >= 54.and.IN <= 90) then
      ISOTAB=1
      else if(IZ == 53.and.IN >= 55.and.IN <= 91) then
      ISOTAB=1
      else if(IZ == 54.and.IN >= 56.and.IN <= 93) then
      ISOTAB=1
      else if(IZ == 55.and.IN >= 57.and.IN <= 96) then
      ISOTAB=1
      else if(IZ == 56.and.IN >= 58.and.IN <= 95) then
      ISOTAB=1
      else if(IZ == 57.and.IN >= 61.and.IN <= 96) then
      ISOTAB=1
      else if(IZ == 58.and.IN >= 63.and.IN <= 97) then
      ISOTAB=1
      else if(IZ == 59.and.IN >= 62.and.IN <= 96) then
      ISOTAB=1
      else if(IZ == 60.and.IN >= 67.and.IN <= 97) then
      ISOTAB=1
      else if(IZ == 61.and.IN >= 69.and.IN <= 98) then
      ISOTAB=1
      else if(IZ == 62.and.IN >= 69.and.IN <= 99) then
      ISOTAB=1
      else if(IZ == 63.and.IN >= 71.and.IN <= 100) then
      ISOTAB=1
      else if(IZ == 64.and.IN >= 73.and.IN <= 101) then
      ISOTAB=1
      else if(IZ == 65.and.IN >= 74.and.IN <= 102) then
      ISOTAB=1
      else if(IZ == 66.and.IN >= 75.and.IN <= 103) then
      ISOTAB=1
      else if(IZ == 67.and.IN >= 76.and.IN <= 105) then
      ISOTAB=1
      else if(IZ == 68.and.IN >= 77.and.IN <= 106) then
      ISOTAB=1
      else if(IZ == 69.and.IN >= 77.and.IN <= 108) then
      ISOTAB=1
      else if(IZ == 70.and.IN >= 78.and.IN <= 110) then
      ISOTAB=1
      else if(IZ == 71.and.IN >= 79.and.IN <= 113) then
      ISOTAB=1
      else if(IZ == 72.and.IN >= 82.and.IN <= 113) then
      ISOTAB=1
      else if(IZ == 73.and.IN >= 83.and.IN <= 114) then
      ISOTAB=1
      else if(IZ == 74.and.IN >= 84.and.IN <= 116) then
      ISOTAB=1
      else if(IZ == 75.and.IN >= 85.and.IN <= 118) then
      ISOTAB=1
      else if(IZ == 76.and.IN >= 86.and.IN <= 120) then
      ISOTAB=1
      else if(IZ == 77.and.IN >= 89.and.IN <= 121) then
      ISOTAB=1
      else if(IZ == 78.and.IN >= 90.and.IN <= 124) then
      ISOTAB=1
      else if(IZ == 79.and.IN >= 92.and.IN <= 126) then
      ISOTAB=1
      else if(IZ == 80.and.IN >= 94.and.IN <= 128) then
      ISOTAB=1
      else if(IZ == 81.and.IN >= 95.and.IN <= 129) then
      ISOTAB=1
      else if(IZ == 82.and.IN >= 96.and.IN <= 132) then
      ISOTAB=1
      else if(IZ == 83.and.IN >= 102.and.IN <= 133) then
      ISOTAB=1
      else if(IZ == 84.and.IN >= 106.and.IN <= 134) then
      ISOTAB=1
      else if(IZ == 85.and.IN >= 109.and.IN <= 138) then
      ISOTAB=1
      else if(IZ == 86.and.IN >= 112.and.IN <= 142) then
      ISOTAB=1
      else if(IZ == 87.and.IN >= 114.and.IN <= 145) then
      ISOTAB=1
      else if(IZ == 88.and.IN >= 116.and.IN <= 146) then
      ISOTAB=1
      else if(IZ == 89.and.IN >= 118.and.IN <= 145) then
      ISOTAB=1
      else if(IZ == 90.and.IN >= 120.and.IN <= 147) then
      ISOTAB=1
      else if(IZ == 91.and.IN >= 122.and.IN <= 147) then
      ISOTAB=1
      else if(IZ == 92.and.IN >= 126.and.IN <= 150) then
      ISOTAB=1
      else if(IZ == 93.and.IN >= 132.and.IN <= 151) then
      ISOTAB=1
      else if(IZ == 94.and.IN >= 134.and.IN <= 153) then
      ISOTAB=1
      else if(IZ == 95.and.IN >= 137.and.IN <= 153) then
      ISOTAB=1
      else if(IZ == 96.and.IN >= 136.and.IN <= 156) then
      ISOTAB=1
      else if(IZ == 97.and.IN >= 140.and.IN <= 156) then
      ISOTAB=1
      else if(IZ == 98.and.IN >= 140.and.IN <= 158) then
      ISOTAB=1
      else if(IZ == 99.and.IN >= 142.and.IN <= 158) then
      ISOTAB=1
      else if(IZ == 100.and.IN >= 142.and.IN <= 160) then
      ISOTAB=1
      else if(IZ == 101.and.IN >= 146.and.IN <= 160) then
      ISOTAB=1
      else if(IZ == 102.and.IN >= 148.and.IN <= 161) then
      ISOTAB=1
      else if(IZ == 103.and.IN >= 149.and.IN <= 162) then
      ISOTAB=1
      else if(IZ == 104.and.IN >= 149.and.IN <= 162) then
      ISOTAB=1
      else if(IZ == 105.and.IN >= 150.and.IN <= 161) then
      ISOTAB=1
      else if(IZ == 106.and.IN >= 153.and.IN <= 160) then
      ISOTAB=1
      else if(IZ == 107.and.IN >= 153.and.IN <= 159) then
      ISOTAB=1
      else if(IZ == 108.and.IN >= 155.and.IN <= 159) then
      ISOTAB=1
      else if(IZ == 109.and.IN >= 157.and.IN <= 159) then
      ISOTAB=1
      else if(IZ >= 110) then
      ISOTAB=1
      else
      ISOTAB=0
      endif

      RETURN
      END
C---------------------- ALSO IN LQMD ----------------
      function sl(i)                                  ! Restored: evolt.f,subroutine width0 CSM
      
      COMMON /NMASS/npr, nta
      c1 = 0.08
      c2 = 0.92
      if (i <= npr) then 
        sl = (c1*float(npr)**0.3333333 + c2)**2.    ! Restored comment. CSM
      else
        sl = (c1*float(nta)**0.3333333 + c2)**2.    ! Restored comment. CSM
      endif
      return
      
      end
      
      ! Exponential Functions with Strict Numerical Protection
      function xexp(x)
      
      if (x <= -87.) then
        xexp = exp(-87.)
      else if (x >= 87.) then
        xexp = exp(87.)
      else
        xexp = exp(x)
      end if
      return
      
      end

      function scou(ncc, al)

      if (ncc == 1 .or. ncc == 5000) then
        scou = 0.
        return
      end if
      c4 = sqrt(al)
      c5 = 0.00144
      dr = (float(ncc) - 1.)/100. ! fm
      t1 = dr/c4 
      call erf(t1,t3)                                                        
      scou = c5/dr*t3
      return
      
      end

C*********************************************************************
C     Error function erf
C     IN: (x) --> OUT: (er)      
      subroutine erf(x,er) 

      p = .3275911                                                              
      a1 = .254829592                                                           
      a2 = -.284496736                                                          
      a3 = 1.421413741                                                          
      a4 = -1.453152027                                                         
      a5 = 1.061405429                                                          
      t = 1./(1 + p*x)                                               
      ! Rational Approximation
      er = 1 - (a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5)*xexp(-x**2)                  
      return
      
      end

C*********************************************************************
C     Lorentz Transformation with velocity (BX,BY,BZ) IN LQMD
C     IN: (A0,AX,AY,AZ) --> OUT: (AP0,APX,APY,APZ)
      subroutine LORENTZ(BX, BY, BZ, A0, AX, AY, AZ, AP0, APX, APY, APZ)
      
      real*8 BS, G, BFACTOR, BS2
      BS = DBLE(BX)**2 + DBLE(BY)**2 + DBLE(BZ)**2
      BS2 = 1D0 - BS
      if (BS2 <= 1D-10) BS2 = 1D-10
      G = 1D0/DSQRT(BS2)
      if (BS < 1D-8) then
        BFACTOR = 0D0
      else
        BFACTOR = (G - 1D0)/BS*(BX*AX + BY*AY + BZ*AZ)
      endif
      APX = AX + BFACTOR*BX - G*A0*BX
      APY = AY + BFACTOR*BY - G*A0*BY
      APZ = AZ + BFACTOR*BZ - G*A0*BZ
      AP0 = G*(A0 - (BX*AX + BY*AY + BZ*AZ))
      return
      
      end
