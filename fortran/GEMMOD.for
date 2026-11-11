!***** GEMINI
!***** Nucl. Phys. A 483 (1988) 371
!------------------------------begin:   MODULE gem_kind-----------------------------
      MODULE gem_kind
      IMPLICIT NONE

      INTEGER, PARAMETER, PUBLIC :: r4 = 4  !single precision
      INTEGER, PARAMETER, PUBLIC :: r8 = 8  !double
	 
	INTEGER, PUBLIC :: LL_SHOW

      TYPE, PUBLIC :: info
      REAL    ::  Excn,dexcn,R_num,t_delay,sig_delay,Kramer,ratio
      REAL    ::  E2_strength, E1_strength, aden_0
      REAL    ::  aden_k0, aden_kinf
      REAL    ::  a_scale, b_scale, threshold
      LOGICAL :: diagnostics,time_flag,I_angle,quantum
      LOGICAL :: sharp_delay,tl_iwbc,polarization 
      LOGICAL :: time_wise,lestone
      LOGICAL :: k_sum,j_flag,ssum,fluct,expand,ex_flag
      INTEGER :: Zcn,Acn,lmin,lmax,N_weight,imf_option,Z_IMF_min
      INTEGER :: exotic_index,mass_option,aden_type
      CHARACTER (LEN=40) :: reaction
      CHARACTER (LEN=12) :: day
      END TYPE info

      TYPE (info), PUBLIC :: para 

      TYPE, PUBLIC :: vell
      REAL (KIND=r4) :: v,theta,phi
      END TYPE vell
      TYPE, PUBLIC :: spinn
      REAL (KIND=r4) :: theta,phi
      END TYPE spinn

      TYPE, PUBLIC :: identity
      INTEGER :: Z,A,Id,level
      REAL (KIND=r4) :: Ex,time,J,K,M,exx
      TYPE (spinn) :: spin
      TYPE (vell) :: vel 
      END TYPE identity

      INTEGER, PARAMETER, PUBLIC :: l_max_quantum = 50
      INTEGER, PUBLIC :: A_cut

      REAL (KIND=r4), PARAMETER, PUBLIC :: pi = 3.1415926
	REAL (KIND=r4), PARAMETER, PUBLIC :: r0 = 1.16
	REAL (KIND=r4), PARAMETER, PUBLIC :: sep = 2.0
	REAL (KIND=r4), PARAMETER, PUBLIC :: r0ld = 1.2249

	END MODULE gem_kind
!------------------------------end:   MODULE gem_kind-----------------------------------

!------------------------------begin:   MODULE trig-----------------------------------
	MODULE trig
! the following are simple trigometric functions with
! arguments in degrees - they exit are intrinsic functions
! on Compact (Didital) compliers and so are not needed.
! However they are not part of standard fortran95.
	IMPLICIT NONE

	PUBLIC :: SIND, COSD, TAND, ASIND, ACOSD, ATAND, ATAN2D
	PRIVATE :: SIND4, COSD4, TAND4, ASIND4, ACOSD4, ATAND4, ATAN2d4
	PRIVATE :: SIND8, COSD8, TAND8, ASIND8, ACOSD8, ATAND8, ATAN2d8

	REAL (KIND=4), PARAMETER, PRIVATE :: pi4 = 3.141592653589793238_4
	REAL (KIND=8), PARAMETER, PRIVATE :: pi8 = 3.141592653589793238_8

	INTERFACE SIND
      MODULE PROCEDURE SIND4,SIND8
	END INTERFACE SIND

	INTERFACE COSD
      MODULE PROCEDURE COSD4,COSD8
      END INTERFACE COSD

	INTERFACE TAND
      MODULE PROCEDURE TAND4,TAND8
	END INTERFACE TAND

	INTERFACE ASIND
      MODULE PROCEDURE ASIND4,ASIND8
	END INTERFACE ASIND

	INTERFACE ACOSD
      MODULE PROCEDURE ACOSD4,ACOSD8
	END INTERFACE ACOSD

	INTERFACE ATAND
      MODULE PROCEDURE ATAND4,ATAND8
	END INTERFACE ATAND

	INTERFACE ATAN2D
      MODULE PROCEDURE ATAN2D4,ATAN2D8
	END INTERFACE ATAN2D

      CONTAINS
      ELEMENTAL FUNCTION SIND4(theta) RESULT(funct_out)
	  REAL (kind=4),INTENT(IN) :: theta 
	  REAL (KIND=4) :: funct_out
      funct_out = SIN(theta*pi4/180.0)
      RETURN
      END FUNCTION SIND4

      ELEMENTAL FUNCTION SIND8(theta) RESULT(funct_out)
      REAL (kind=8),INTENT(IN) :: theta 
      REAL (KIND=8) :: funct_out

      funct_out = SIN(theta*pi8/180.0)
      RETURN
      END FUNCTION SIND8

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ELEMENTAL FUNCTION COSD4(theta) RESULT(funct_out)
        REAL (kind=4),INTENT(IN) :: theta
        REAL (KIND=4) :: funct_out

        funct_out = COS(theta*pi4/180.0)
        RETURN
        END FUNCTION COSD4

        ELEMENTAL FUNCTION COSD8(theta) RESULT(funct_out)
        REAL (kind=8),INTENT(IN) :: theta 
        REAL (KIND=8) :: funct_out

        funct_out = COS(theta*pi8/180.0)
        RETURN
        END FUNCTION COSD8

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
       ELEMENTAL  FUNCTION TAND4(theta) RESULT(funct_out)
       REAL (kind=4), INTENT(IN) :: theta 
       REAL (KIND=4) :: funct_out
     
       funct_out = TAN(theta*pi4/180.0)
     
       RETURN
       END FUNCTION TAND4

       ELEMENTAL  FUNCTION TAND8(theta) RESULT(funct_out)
       REAL (kind=8), INTENT(IN) :: theta 
       REAL (KIND=8) :: funct_out
     
       funct_out = TAN(theta*pi8/180.0)
     
       RETURN
       END FUNCTION TAND8
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ELEMENTAL FUNCTION ASIND4(x) RESULT(funct_out)
       REAL (kind=4), INTENT(IN) :: x 
       REAL (KIND=4) :: funct_out
     
       funct_out = ASIN(x)*180.0/pi4
     
       RETURN
       END FUNCTION ASIND4

       ELEMENTAL FUNCTION ASIND8(x) RESULT(funct_out)
       REAL (kind=8), INTENT(IN) :: x 
       REAL (KIND=8) :: funct_out
     
       funct_out = ASIN(x)*180.0/pi8
     
       RETURN
       END FUNCTION ASIND8

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ELEMENTAL FUNCTION ACOSD4(x) RESULT(funct_out)
       REAL (kind=4), INTENT(IN) :: x 
       REAL (KIND=4) :: funct_out
     
       funct_out = ACOS(x)*180.0/pi4
     
       RETURN
       END FUNCTION ACOSD4

       ELEMENTAL FUNCTION ACOSD8(x) RESULT(funct_out)
       REAL (kind=8), INTENT(IN) :: x 
       REAL (KIND=8) :: funct_out
     
       funct_out = ACOS(x)*180.0/pi8
     
       RETURN
       END FUNCTION ACOSD8
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ELEMENTAL FUNCTION ATAND4(x) RESULT(funct_out)
       REAL (kind=4), INTENT(IN) :: x 
       REAL (KIND=4) :: funct_out
     
       funct_out = ATAN(x)*180.0/pi4
     
       RETURN
       END FUNCTION ATAND4

       ELEMENTAL FUNCTION ATAND8(x) RESULT(funct_out)
       REAL (kind=8), INTENT(IN) :: x 
       REAL (KIND=8) :: funct_out
     
       funct_out = ATAN(x)*180.0/pi8
     
       RETURN
       END FUNCTION ATAND8
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       ELEMENTAL FUNCTION ATAN2D4(x,y) RESULT(funct_out)
       REAL (kind=4), INTENT(IN) :: x,y 
       REAL (KIND=4) :: funct_out
  
       IF (x == 0.0 .AND. y == 0.0) THEN
        funct_out = 0.0
        RETURN
       END IF      
       funct_out = ATAN2(x,y)*180.0/pi4
       RETURN
       END FUNCTION ATAN2D4
  
       ELEMENTAL FUNCTION ATAN2D8(x,y) RESULT(funct_out)
       REAL (kind=8), INTENT(IN) :: x,y 
       REAL (KIND=8) :: funct_out
 
       IF (x == 0.0 .AND. y == 0.0) THEN
        funct_out = 0.0
        RETURN
       END IF      
        funct_out = ATAN2(x,y)*180.0/pi8
       RETURN
       END FUNCTION ATAN2D8
       END MODULE trig
!------------------------------end:   MODULE trig-----------------------------------


!------------------------------begin:   MODULE sierk_stuff--------------------------
	 MODULE sierk_stuff
!    original subroutine of Arnie Sierk
!srk This subroutine returns the barrier height bfis, the ground-state
!srk energy E_rot, in MeV, and the angular momentum at which the fission
!srk barrier disappears,  Jmax,  in units of h-bar,
!srk when called with integer arguments z, the atomic number,
!srk a, the atomic mass number, and J, the angular momentum in units
!srk of h-bar, (Planck's constant divided by 2*pi).
!
!srk      The fission barrier for  J = 0  is calculated from a 7th order
!srk fit in two variables to 638 calculated fission barriers for z value
!srk from 20 to 110.  These  638 barriers are fit with an rms deviation
!srk 0.10 MeV by this 49-parameter function.
!srk If  barfit  is called with (z,a) values outside the range of the
!srk the barrier height is set to 999 MeV, and a message is
!srk printed on the default output file.
!
!srk      For J values not equal to zero, the values of
!srk J at which the barrier is  80%  and  20%  of the J=0 value are
!srk respectively fit to 20-parameter functions of  Z  and  A, over a more
!srk restricted range of  A  values, than is the case for  J = 0.
!srk The value of J where the barrier disappears, Jmax, for 61 nuclei
!srk is fit to a 35-parameter function of Z and A,  with the same range
!srk Z  and  A  values as  J-80  and J-20.
!         Once again, if an  (iz,ia) pair is outside of the range of
!    validity of the fit, the barrier value is set to  0.0  and a message
!    is printed.  These three values  (Bfis(J=0),J-80, and J-20) and the
!    constraints  of  Bfis = 0 and  d(Bfis)/dJ = 0 at J = Jmax and J = 0
!    lead to a fifth-order fit to Bfis(J) for J> L-20.  The first three
!    constraints lead to a third-order fit for the region J < J-20.
!
!
!         The ground-state energies are calculated from a 175-parameter
!    fit in Z, A, and L to 329 ground-state energies for 36 different
!    Z  and  A  values.
!srk (The range of Z and A is the same as for J-80, J-20, and J-max)
!
!srk      The calculated barriers from which the fits were
!srk made were calculated in 1983-1984 by A. J. Sierk of Los Alamos
!srk National Laboratory   Group T-9, using  Yukawa-plus-exponential double
!srk folded nuclear energy, exact Coulomb diffuseness corrections,
!srk and diffuse-matter moments of inertia. The parameters of the model
!srk are those derived by Moller and Nix in 1979:
!srk r-0 = 1.16 fm, as = 21.13 MeV, kappa-s = 2.3  a = 0.68 fm.
!srk The diffuseness of the matter and charge distributions used
!srk corresponds to a surface diffuseness parameter (defined by Myers)
!srk of 0.99 fm.  The calculated barriers for J = 0 are
!srk accurate to a little less than 0.1 MeV;  the output from this
!srk subroutine is a little less accurate.  Worst errors may be as large
!srk as 0.5 MeV; characteristic uncertainty is in the range of 0.1-0.2
!srk MeV.  
!srk      Jmax may be in error by up to 1.0 for light nuclei,
!srk but is probably good to about 0.2-0.3 for heavy nuclei which do
!srk not have prolate ground states  (Z >75 or so).
!    The values of egs are generally approximated to within
!    about 0.1-0.2 MeV;  the largest deviation is about 0.5 MeV,
!    near L-I for light nuclei.c
!         The rms deviation of Lmax from the 61 input values is 0.31
!    h-bar.  The approximate value is nearly always within
!    0.5 h-bar of the calculated one.
!
!srk written by A. J. Sierk,  LANL  T-9
!srk Version 1.0   February, 1984
!    Version 1.1   January, 1985  Improved coefficients in egs and Lmax
!    Version 1.2   September, 1985  Improved Lmax, egs coefficients
!    Version 1.21   June, 1986   minor changes made

	USE gem_kind
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: jmax_sierk, bfis_sierk, yrast_sierk, m_inertia_sierk, 
     &          bs_sierk
	PUBLIC :: fisrot
	PRIVATE :: lpoly

	REAL (KIND=r4), PRIVATE :: amin,amax
	REAL (KIND=r8), PRIVATE :: aa,zz,J20,J80,bfis,JJ,JMAX
	REAL (KIND=r8), PRIVATE :: zzz,x,y,q,qa,qb,aj,ak,a1,a2,E_rot
	REAL (KIND=r8), DIMENSION(7), PRIVATE :: pa, pz
	REAL (KIND=r8), DIMENSION(9), PRIVATE :: pl
	REAL (KIND=r8), DIMENSION(5), PRIVATE :: paa
	REAL (KIND=r8), DIMENSION(8), PRIVATE :: pzz

	REAL (KIND=r8), PRIVATE :: aizro,ai70,aimax,ai95,aimin,bizro,
     &     bi70,bimax,bi95
	REAL (KIND=r8), PRIVATE :: aimax2,ai952,ff1,ff2,fg1,fg2,aimaxh,
     &     aimidh
	REAL (KIND=r8), PRIVATE :: saizro,sai70,sai95,sai952,simax2,
     &     sbimax,sbi70,sbi95
	REAL (KIND=r8), PRIVATE :: sbizro,q1,q2,q3,q4,q5,q6,gam,saimax,
     &     bb,gam2,aa2
	REAL (KIND=r8), PRIVATE :: bb2,aa3,saimid
	REAL (KIND=r8), PRIVATE :: bb3,gam3,alpha,beta,silt,sjlt,silt2,
     &     saimin, saimx
	REAL (KIND=r8), PRIVATE :: sigt,sjgt,sigt2,f1,f3,f1m,f2,f4,f2m,
     &     aa4,bb4
	REAL (KIND=r4), PRIVATE :: Z, A


	REAL (KIND=r8), DIMENSION(5,4), PARAMETER, PRIVATE :: 
     &     emncof = RESHAPE((/
     & -9.011e2_r8,-1.40818e3_r8, 2.77e3_r8,-7.06695e2_r8, 8.89867e2_r8,
     &  1.35355e4_r8,-2.03847e4_r8, 1.09384e4_r8,-4.86297e3_r8,
     & -6.18603e2_r8, -3.26367e3_r8, 1.62447e3_r8, 1.36856e3_r8, 
     &  1.31731e3_r8, 1.53372e2_r8, 7.48863e3_r8,-1.21581e4_r8, 
     &  5.50281e3_r8,-1.33630e3_r8, 5.05367e-2_r8/), (/5,4/))

	REAL (KIND=r8), DIMENSION(5,4), PARAMETER, PRIVATE :: 
     &     elmcof = RESHAPE((/
     &  1.84542e3_r8,-5.64002e3_r8, 5.66730e3_r8,-3.15150e3_r8, 
     &  9.54160e2_r8,-2.24577e3_r8, 8.56133e3_r8,-9.67348e3_r8,
     &  5.81744e3_r8,-1.86997e3_r8, 2.79772e3_r8,-8.73073e3_r8,
     &  9.19706e3_r8,-4.91900e3_r8, 1.37283e3_r8,-3.01866e1_r8, 
     &  1.41161e3_r8,-2.85919e3_r8, 2.13016e3_r8,-6.49072e2_r8/),
     &  (/5,4/))

	REAL (KIND=r8), DIMENSION(7,5), PARAMETER, PRIVATE :: 
     &     emxcof = RESHAPE((/
     & -4.10652732e6_r8, 1.00064947e7_r8,-1.09533751e7_r8,
     &  7.84797252e6_r8,-3.78574926e6_r8, 1.12237945e6_r8,
     & -1.77561170e5_r8, 1.08763330e7_r8,-2.63758245e7_r8,
     &  2.85472400e7_r8,-2.01107467e7_r8, 9.48373641e6_r8,      
     & -2.73438528e6_r8, 4.13247256e5_r8,-8.76530903e6_r8, 
     &  2.14250513e7_r8,-2.35799595e7_r8, 1.70161347e7_r8,
     & -8.23738190e6_r8, 2.42447957e6_r8,-3.65427239e5_r8, 
     &  6.30258954e6_r8,-1.52999004e7_r8, 1.65640200e7_r8,      
     & -1.16695776e7_r8, 5.47369153e6_r8,-1.54986342e6_r8, 
     &  2.15409246e5_r8,-1.45539891e6_r8, 3.64961835e6_r8,
     & -4.21267423e6_r8, 3.24312555e6_r8,-1.67927904e6_r8, 
     &  5.23795062e5_r8,-7.66576599e4_r8/),(/7,5/))

	REAL (KIND=r8), DIMENSION(7,7), PARAMETER, PRIVATE :: 
     &     elzcof = RESHAPE((/ 
     &  5.11819909e5_r8,-1.30303186e6_r8, 1.90119870e6_r8,
     & -1.20628242e6_r8, 5.68208488e5_r8, 5.48346483e4_r8,
     & -2.45883052e4_r8,-1.13269453e6_r8, 2.97764590e6_r8,
     & -4.54326326e6_r8, 3.00464870e6_r8,-1.44989274e6_r8,  
     & -1.02026610e5_r8, 6.27959815e4_r8, 1.37543304e6_r8,
     & -3.65808988e6_r8, 5.47798999e6_r8,-3.78109283e6_r8, 
     &  1.84131765e6_r8, 1.53669695e4_r8,-6.96817834e4_r8,
     & -8.56559835e5_r8, 2.48872266e6_r8,-4.07349128e6_r8,  
     &  3.12835899e6_r8,-1.62394090e6_r8, 1.19797378e5_r8, 
     &  4.25737058e4_r8, 3.28723311e5_r8,-1.09892175e6_r8, 
     &  2.03997269e6_r8,-1.77185718e6_r8, 9.96051545e5_r8,
     & -1.53305699e5_r8,-1.12982954e4_r8, 4.15850238e4_r8,  
     &  7.29653408e4_r8,-4.93776346e5_r8, 6.01254680e5_r8,
     & -4.01308292e5_r8, 9.65968391e4_r8,-3.49596027e3_r8,
     & -1.82751044e5_r8, 3.91386300e5_r8,-3.03639248e5_r8,
     &  1.15782417e5_r8,-4.24399280e3_r8,-6.11477247e3_r8,  
     &  3.66982647e2_r8/),(/7,7/))

	REAL (KIND=r8), DIMENSION(88), PARAMETER, PRIVATE :: 
     &     egscof_a = (/
     &  -1.781665232e6_r8,-2.849020290e6_r8, 9.546305856e5_r8,
     &   2.453904278e5_r8, 3.656148926e5_r8, 4.358113622e6_r8, 
     &   6.960182192e6_r8,-2.381941132e6_r8,-6.262569370e5_r8,
     &  -9.026606463e5_r8,-4.804291019e6_r8,-7.666333374e6_r8,  
     &   2.699742775e6_r8, 7.415602390e5_r8, 1.006008724e6_r8, 
     &   3.505397297e6_r8, 5.586825123e6_r8,-2.024820713e6_r8,
     &  -5.818008462e5_r8,-7.353683218e5_r8,-1.740990985e6_r8,
     &  -2.759325148e6_r8, 1.036253535e6_r8, 3.035749715e5_r8,  
     &   3.606919356e5_r8, 5.492532874e5_r8, 8.598827288e5_r8,
     &  -3.399809581e5_r8,-9.852362945e4_r8,-1.108872347e5_r8,
     &  -9.229576432e4_r8,-1.431344258e5_r8, 5.896521547e4_r8, 
     &   1.772385043e4_r8, 1.845424227e4_r8, 4.679351387e6_r8,  
     &   7.707630513e6_r8,-2.718115276e6_r8,-9.845252314e5_r8,
     &  -1.107173456e6_r8,-1.137635233e7_r8,-1.870617878e7_r8, 
     &   6.669154225e6_r8, 2.413451470e6_r8, 2.691480439e6_r8, 
     &   1.237627138e7_r8, 2.030222826e7_r8,-7.334289876e6_r8,  
     &  -2.656357635e6_r8,-2.912593917e6_r8,-8.854155353e6_r8,
     &  -1.446966194e7_r8, 5.295832834e6_r8, 1.909275233e6_r8, 
     &   2.048899787e6_r8, 4.290642787e6_r8, 6.951223648e6_r8,
     &  -2.601557110e6_r8,-9.129731614e5_r8,-9.627344865e5_r8,  
     &  -1.314924218e6_r8,-2.095971932e6_r8, 8.193066795e5_r8, 
     &   2.716279969e5_r8, 2.823297853e5_r8, 2.131536582e5_r8, 
     &   3.342907992e5_r8,-1.365390745e5_r8,-4.417841315e4_r8,
     &  -4.427025540e4_r8,-3.600471364e6_r8,-5.805932202e6_r8,  
     &   1.773029253e6_r8, 4.064280430e5_r8, 7.419581557e5_r8, 
     &   8.829126250e6_r8, 1.422377198e7_r8,-4.473342834e6_r8,
     &  -1.073350611e6_r8,-1.845960521e6_r8,-9.781712604e6_r8,
     &  -1.575666314e7_r8, 5.161226883e6_r8, 1.341287330e6_r8,  
     &   2.083994843e6_r8, 7.182555931e6_r8, 1.156915972e7_r8,
     &  -3.941330542e6_r8/)

	REAL (KIND=r8), DIMENSION(87), PARAMETER, PRIVATE :: egscof_b = (/ 
     & -1.108259560e6_r8,-1.543982755e6_r8,-3.579820035e6_r8,
     & -5.740079339e6_r8, 2.041827680e6_r8, 5.981648181e5_r8, 
     &  7.629263278e5_r8, 1.122573403e6_r8, 1.777161418e6_r8,
     & -6.714631146e5_r8,-1.952833263e5_r8,-2.328129775e5_r8,  
     & -1.839672155e5_r8,-2.871137706e5_r8, 1.153532734e5_r8, 
     &  3.423868607e4_r8, 3.738902942e4_r8, 2.421750735e6_r8, 
     &  4.107929841e6_r8,-1.302310290e6_r8,-5.267906237e5_r8,
     & -6.197966854e5_r8,-5.883394376e6_r8,-9.964568970e6_r8,  
     &  3.198405768e6_r8, 1.293156541e6_r8, 1.506909314e6_r8,
     &  6.387411818e6_r8, 1.079547152e7_r8,-3.517981421e6_r8,
     & -1.424705631e6_r8,-1.629099740e6_r8,-4.550695232e6_r8,
     & -7.665548805e6_r8, 2.530844204e6_r8, 1.021187317e6_r8,  
     &  1.141553709e6_r8, 2.182540324e6_r8, 3.646532772e6_r8,
     & -1.228378318e6_r8,-4.813626449e5_r8,-5.299974544e5_r8,
     & -6.518758807e5_r8,-1.070414288e6_r8, 3.772592079e5_r8, 
     &  1.372024952e5_r8, 1.505359294e5_r8, 9.952777968e4_r8,  
     &  1.594230613e5_r8,-6.029082719e4_r8,-2.023689807e4_r8,
     & -2.176008230e4_r8,-4.902668827e5_r8,-8.089034293e5_r8, 
     &  1.282510910e5_r8,-1.704435174e4_r8, 8.876109934e4_r8, 
     &  1.231673941e6_r8, 2.035989814e6_r8,-3.727491110e5_r8,  
     &  4.071377327e3_r8,-2.375344759e5_r8,-1.429330809e6_r8,
     & -2.376692769e6_r8, 5.216954243e5_r8, 7.268703575e4_r8, 
     &  3.008350125e5_r8, 1.114306796e6_r8, 1.868800148e6_r8,
     & -4.718718351e5_r8,-1.215904582e5_r8,-2.510379590e5_r8,  
     & -5.873353309e5_r8,-9.903614817e5_r8, 2.742543392e5_r8,
     &  9.055579135e4_r8, 1.364869036e5_r8, 1.895325584e5_r8, 
     &  3.184776808e5_r8,-9.500485442e4_r8, -3.406036086e4_r8,
     & -4.380685984e4_r8,-2.969272274e4_r8,-4.916872669e4_r8,
     &  1.596305804e4_r8, 5.741228836e3_r8, 6.669912421e3_r8 /)

	REAL (KIND=r8), DIMENSION(5,7,5), PARAMETER, PRIVATE :: egscof = 
     &	 RESHAPE((/egscof_a,egscof_b/),(/5,7,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &     aizroc = RESHAPE((/
     &   2.34441624e4_r8,-5.88023986e4_r8, 6.37939552e4_r8,
     &  -4.79085272e4_r8, 2.27517867e4_r8,-5.35372280e3_r8,
     &  -4.19782127e4_r8, 1.09187735e5_r8,-1.24597673e5_r8, 
     &   9.93997182e4_r8,-4.95141312e4_r8, 1.19847414e4_r8,    
     &   4.18237803e4_r8,-1.05557152e5_r8, 1.16142947e5_r8,
     &  -9.00443421e4_r8, 4.48976290e4_r8,-1.10161792e4_r8,
     &  -8.27172333e3_r8, 2.49194412e4_r8,-3.39090117e4_r8, 
     &   3.33727886e4_r8,-1.98040399e4_r8, 5.37766241e3_r8,    
     &   5.79695749e2_r8,-1.61762346e3_r8, 2.14044262e3_r8,
     &  -3.55379785e3_r8, 3.25502799e3_r8,-1.15583400e3_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &     ai70c = RESHAPE((/
     &  3.11420101e4_r8,-7.54335155e4_r8, 7.74456473e4_r8,
     & -4.79993065e4_r8, 2.23439118e4_r8,-4.81961155e3_r8,
     & -7.24025043e4_r8, 1.72276697e5_r8,-1.72027101e5_r8, 
     &  1.03891065e5_r8,-4.83180786e4_r8, 1.08040504e4_r8,    
     &  7.14932917e4_r8,-1.72792523e5_r8, 1.75814382e5_r8,
     & -1.07245918e5_r8, 4.86163223e4_r8,-1.10623761e4_r8,
     & -2.87206866e4_r8, 6.76667976e4_r8,-6.50167483e4_r8, 
     &  3.67161268e4_r8,-1.74755753e4_r8, 4.67495427e3_r8,    
     &  1.67914908e4_r8,-3.97304542e4_r8, 3.81446552e4_r8,
     & -2.04628156e4_r8, 7.20091899e3_r8,-1.49978283e3_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &     ai95c = RESHAPE((/
     & -6.17201449e5_r8, 1.45561724e6_r8,-1.47514522e6_r8, 
     &  9.37798508e5_r8,-3.74435017e5_r8, 7.81254880e4_r8, 
     &  1.24304280e6_r8,-2.94179116e6_r8, 3.00170753e6_r8,
     & -1.92737183e6_r8, 7.79238772e5_r8,-1.64803784e5_r8,    
     & -1.49648799e6_r8, 3.52658199e6_r8,-3.56784327e6_r8, 
     &  2.26413602e6_r8,-9.02243251e5_r8, 1.88619658e5_r8, 
     &  7.27293223e5_r8,-1.72140677e6_r8, 1.75634889e6_r8,
     & -1.12885888e6_r8, 4.57150814e5_r8,-9.74833991e4_r8,    
     & -3.75965723e5_r8, 8.83032946e5_r8,-8.87134867e5_r8, 
     &  5.58350462e5_r8, -2.20433857e5_r8, 4.62178756e4_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &     aimaxc = RESHAPE((/
     & -1.07989556e6_r8, 2.54617598e6_r8,-2.56762409e6_r8, 
     &  1.62814115e6_r8,-6.39575059e5_r8, 1.34017942e5_r8, 
     &  2.17095357e6_r8,-5.13081589e6_r8, 5.19610055e6_r8,
     & -3.31651644e6_r8, 1.31229476e6_r8,-2.77511450e5_r8,    
     & -2.66020302e6_r8, 6.26593165e6_r8,-6.31060776e6_r8, 
     &  3.99082969e6_r8,-1.56447660e6_r8, 3.25613262e5_r8, 
     &  1.29464191e6_r8,-3.05746938e6_r8, 3.09487138e6_r8,
     & -1.97160118e6_r8, 7.79696064e5_r8,-1.63704652e5_r8,    
     & -7.13073644e5_r8, 1.67482279e6_r8,-1.67984330e6_r8, 
     &  1.05446783e6_r8,-4.10928559e5_r8, 8.43774143e4_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &      ai952c = RESHAPE((/
     & -7.37473153e5_r8, 1.73682827e6_r8,-1.75850175e6_r8, 
     &  1.11320647e6_r8,-4.41842735e5_r8, 9.02463457e4_r8, 
     &  1.49541980e6_r8,-3.53222507e6_r8, 3.59762757e6_r8,
     & -2.29652257e6_r8, 9.21077757e5_r8,-1.90079527e5_r8,    
     & -1.80243593e6_r8, 4.24319661e6_r8,-4.29072662e6_r8, 
     &  2.71416936e6_r8,-1.07624953e6_r8, 2.20863711e5_r8, 
     &  8.86920591e5_r8,-2.09589683e6_r8, 2.13507675e6_r8,
     & -1.36546686e6_r8, 5.48868536e5_r8,-1.14532906e5_r8,    
     & -4.62131503e5_r8, 1.08555722e6_r8,-1.09187524e6_r8, 
     &  6.87308217e5_r8, -2.70986162e5_r8, 5.61637883e4_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(6,5), PARAMETER, PRIVATE :: 
     &     aimax2c = RESHAPE((/
     & -1.16343311e6_r8, 2.74470544e6_r8,-2.77664273e6_r8, 
     &  1.76933559e6_r8,-7.02900226e5_r8, 1.49345081e5_r8, 
     &  2.36929777e6_r8,-5.60655122e6_r8, 5.70413177e6_r8,
     & -3.66528765e6_r8, 1.47006527e6_r8,-3.15794626e5_r8,    
     & -2.82646077e6_r8, 6.66086824e6_r8,-6.72677653e6_r8, 
     &  4.27484625e6_r8,-1.69427298e6_r8, 3.58429081e5_r8, 
     &  1.39112772e6_r8,-3.29007553e6_r8, 3.34544584e6_r8,
     & -2.14723142e6_r8, 8.61118401e5_r8,-1.84500129e5_r8,    
     & -7.21329917e5_r8, 1.69371794e6_r8,-1.69979786e6_r8, 
     &  1.07037781e6_r8, -4.20662028e5_r8, 8.80728361e4_r8/),(/6,5/))

	REAL (KIND=r8), DIMENSION(4,4), PARAMETER, PRIVATE :: 
     &     aimax3c = RESHAPE((/
     & -2.88270282e3_r8, 5.30111305e3_r8,-3.07626751e3_r8, 
     &  6.56709396e2_r8, 5.84303930e3_r8,-1.07450449e4_r8, 
     &  6.24110631e3_r8,-1.33480875e3_r8,-4.20629939e3_r8, 
     &  7.74058373e3_r8,-4.50256063e3_r8, 9.65788439e2_r8,    
     &  1.23820134e3_r8,-2.28228958e3_r8, 1.33181316e3_r8,
     & -2.87363568e2_r8/), (/4,4/))

	REAL (KIND=r8), DIMENSION(4,4), PARAMETER, PRIVATE :: 
     &     aimax4c = RESHAPE((/
     & -3.34060345e3_r8, 6.26384099e3_r8,-3.77635848e3_r8, 
     &  8.57180868e2_r8, 6.76377873e3_r8,-1.26776571e4_r8, 
     &  7.64206952e3_r8,-1.73406840e3_r8,-4.74821371e3_r8, 
     &  8.89857519e3_r8,-5.36266252e3_r8, 1.21614216e3_r8,    
     &  1.46369384e3_r8,-2.74251101e3_r8, 1.65205435e3_r8,
     & -3.74262365e2_r8/), (/4,4/))

	REAL (KIND=r8), DIMENSION(6,4), PARAMETER, PRIVATE :: 
     &     bizroc = RESHAPE((/
     &  5.88982505e2_r8,-1.35630904e3_r8, 1.32932125e3_r8,
     & -7.78518395e2_r8, 2.73122883e2_r8,-3.49600841e1_r8,
     & -9.67701343e2_r8, 2.24594418e3_r8,-2.24303790e3_r8, 
     &  1.35440047e3_r8,-4.96538939e2_r8, 6.66791793e1_r8,    
     &  1.17090267e3_r8,-2.71181535e3_r8, 2.67008958e3_r8,
     & -1.58801770e3_r8, 5.66896359e2_r8,-8.21530057e1_r8,
     & -3.83031864e2_r8, 9.05191483e2_r8,-9.30560410e2_r8, 
     &  5.96618532e2_r8,-2.34403480e2_r8, 3.97909172e1_r8/), (/6,4/))

	REAL (KIND=r8), DIMENSION(6,4), PARAMETER, PRIVATE :: 
     &     bi70c = RESHAPE((/
     &  2.32414810e3_r8,-5.42381778e3_r8, 5.40202710e3_r8,
     & -3.26923144e3_r8, 1.18318943e3_r8,-1.93186467e2_r8,
     & -4.38084778e3_r8, 1.03523570e4_r8,-1.05573803e4_r8, 
     &  6.59901160e3_r8,-2.47601209e3_r8, 4.19497260e2_r8,    
     &  4.35377377e3_r8,-1.01728647e4_r8, 1.01311246e4_r8,
     & -6.14038462e3_r8, 2.21957562e3_r8,-3.62854365e2_r8,
     & -1.84533539e3_r8, 4.41613298e3_r8,-4.59403284e3_r8, 
     &  2.95951225e3_r8,-1.14630148e3_r8, 2.02702459e2_r8/), (/6,4/))

	REAL (KIND=r8), DIMENSION(6,4), PARAMETER, PRIVATE :: 
     &     bi95c = RESHAPE((/
     &  1.55359266e3_r8,-3.58209715e3_r8, 3.50693744e3_r8,
     & -2.03992913e3_r8, 7.05498010e2_r8,-1.49075519e2_r8,
     & -2.86876240e3_r8, 6.77107086e3_r8,-6.90300614e3_r8, 
     &  4.20246063e3_r8,-1.50290693e3_r8, 3.13662258e2_r8,    
     &  2.60138185e3_r8,-5.95414919e3_r8, 5.70261588e3_r8,
     & -3.17188958e3_r8, 9.89207911e2_r8,-1.76320647e2_r8,
     & -1.75198402e3_r8, 4.16635208e3_r8,-4.25212424e3_r8, 
     &  2.59953301e3_r8,-9.09813362e2_r8, 1.51070448e2_r8/), (/6,4/))

	REAL (KIND=r8), DIMENSION(6,4), PARAMETER, PRIVATE :: 
     &     bimaxc = RESHAPE((/
     &  4.17708254e3_r8,-8.59358778e3_r8, 6.46392215e3_r8,
     & -8.84972189e2_r8,-1.59735594e3_r8, 1.39662071e3_r8,
     & -1.56318394e4_r8, 3.54574417e4_r8,-3.35945173e4_r8, 
     &  1.65495998e4_r8,-3.32021998e3_r8,-1.46150905e3_r8,    
     &  1.41292811e4_r8,-3.11818487e4_r8, 2.77454429e4_r8,
     & -1.19628827e4_r8, 1.28008968e3_r8, 1.66111636e3_r8,
     & -1.92878152e4_r8, 4.56505796e4_r8,-4.66413277e4_r8, 
     &  2.89229633e4_r8,-1.07284346e4_r8, 1.50513815e3_r8/), (/6,4/))

	REAL (KIND=r8), DIMENSION(100), PARAMETER, PRIVATE :: b_a = (/
     & -17605.486,    7149.518,   23883.053,   21651.272,    -348.402,    
     &  67617.841,  -43739.724,  -72255.975,  -79130.794,  -44677.571,   
     & -38284.438,    1629.405,   49928.158,   41804.076,   -4744.661,    
     &  32124.997,  -30258.577,  -45416.684,  -40347.682,  -20718.220,    
     & -7009.141,  -12274.997,    5703.797,    4375.711,   -2519.378,    
     &  44633.464,  -19664.988,  -59871.761,  -55426.914,    -834.542,    
     & -168222.448,  111291.480,  179939.092,  198092.569,  111347.893,   
     &  97361.913,   -9404.674, -125757.922, -108003.071,    7218.634,    
     & -79483.915,   77325.943,  113433.935,  100793.473,   51134.715,    
     &  17974.287,   28084.202,  -15187.024,  -12244.758,    4662.638,    
     & -52775.944,   26372.246,   69061.418,   66412.802,    4456.150,    
     &  191743.968, -131672.398, -204896.966, -227867.388, -126609.907,    
     & -115704.145,   22053.324,  146315.778,  131364.763,    1176.787,   
     &  89552.374,  -92216.043, -129901.381, -115369.145,  -56853.433,   
     & -21565.916,  -26595.626,   19490.177,   16758.541,   -2102.079,   
     &  43089.492,  -24205.277,  -54145.651,  -54916.885,   -6982.880,   
     & -148893.417,  105958.013,  157835.761,  178691.634,   97246.541,   
     &  94988.919,  -28100.912, -115999.016, -110379.031,  -10165.655,   
     & -68277.650,   75226.302,  100957.904,   89635.423,   41887.718,   
     &  17720.511,   14837.634,  -17476.245,  -15762.519,   -1297.194/)

	 REAL (KIND=r8), DIMENSION(100), PARAMETER, PRIVATE :: b_b = (/
     & -25673.002,   15584.030,   30291.575,   33008.500,    6161.347,
     &  83671.721,  -60540.161,  -86493.251, -101211.404,  -53525.729,   
     & -56872.290,   21837.845,   65632.326,   67209.922,   11367.245,   
     &  37417.011,  -43832.496,  -56023.110,  -49995.410,  -21521.767,   
     & -10406.752,   -4483.654,   11314.537,   10435.311,    2206.541,   
     &  11115.697,   -6672.659,  -11696.751,  -14183.580,   -3586.280,   
     & -33914.223,   23453.203,   32630.486,   40847.174,   21223.389,   
     &  24704.794,  -10226.568,  -25497.313,  -29082.833,   -7190.865,   
     & -14734.982,   17457.735,   21457.946,   19706.266,    7745.219,   
     &   4283.242,     214.103,   -5058.390,   -4715.987,   -1303.557,   
     &  -3196.271,    1579.967,    2718.777,    3997.275,    1361.541,   
     &   9152.527,   -5182.345,   -7412.465,  -10822.583,   -5804.143,   
     &  -7177.709,    2317.015,    5874.736,    8208.479,    2858.079,   
     &   3919.329,   -4068.562,   -4982.296,   -5065.564,   -1948.636,
     &  -1147.488,     265.063,    1372.992,    1333.320,     416.369,
     &    442.331,    -219.051,    -337.663,    -575.255,    -232.634,
     &  -1217.351,     590.169,     854.438,    1468.791,     811.481,
     &   1031.157,    -227.306,    -651.183,   -1174.558,    -528.960,
     &   -536.534,     413.429,     533.977,     663.064,     278.203,
     &    154.548,     -54.623,    -169.817,    -183.006,     -63.999/)

	 REAL (KIND=r8), DIMENSION(5,5,8), PARAMETER, PRIVATE :: b =
     &	 RESHAPE((/b_a,b_b/),(/5,5,8/))

	 CONTAINS
	 SUBROUTINE jmax_sierk (z_in,a_in,Jmax_s)

! DUMMY ARGUMENTS:
!
!    input: Z (REAL (KIND=4) ::) proton number of nucleus
!           A (REAL (KIND=4) ::) nucleon number of nucleus
!           J (REAL (KIND=4) ::) spin of nucleus (hbar)
!    output:
!         bfis (REAL (KIND=4) ::) - symmetric fission barrier of nucleus
!         Jmax (REAL (KIND=4) ::) - angular momentum at which bfis=0 MeV
!
!
!          C H A N G E   L O G
!
!      Date     | Name  | Description
! ----------------+-------+-----------------------------------------------------
! 21-11-89        |  RJC  | I have beautified this routine of Sierks and made all
!                           the double precision variables single precision.
! ----------------+-------+-----------------------------------------------------
! 30-8-90         |  RJC  | corrected two errors of no probable consequence
!                         ! error 1/ there were two variables Z in routine
!                         ! 2/ there was a el instead of J in if statement
! ----------------+-------+-----------------------------------------------------

	REAL (KIND=r4), INTENT(IN) :: Z_in, A_in
	REAL (KIND=r4), INTENT(OUT), OPTIONAL :: Jmax_s

	INTEGER :: i,k

	z = z_in
	A = a_in

	IF (z < 19.0 .OR. z > 111.0) THEN
       WRITE (UNIT=*, FMT=*) "z=",Z, "out of limits"
       IF (PRESENT(Jmax_s)) jmax_s = -1.0
       RETURN
      END IF

	 amin = 1.2*z + 0.01*z**2.0
	 amax = 5.8*z - 0.024*z**2.0

	 IF (a < amin .OR. a > amax) THEN
       WRITE (UNIT=*, FMT=*) "A out of limits"
       IF (PRESENT(Jmax_s)) jmax_s = -1.0
       RETURN
	 END IF

	aa = 2.5e-3_r8  *REAL(a,KIND=r8)
	zz = 1.0e-2_r8  *REAL(z,KIND=r8)

	jmax= 0.0
	CALL lpoly (zz,7,pz)
	CALL lpoly (aa,7,pa)

	DO i=1,5
	  DO k=1,7
        Jmax = Jmax + emxcof(k,i)*pz(k)*pa(i)
	  END DO
	END DO

	IF (PRESENT(jmax_s)) jmax_s = jmax
	RETURN
	END SUBROUTINE jmax_sierk

	SUBROUTINE bfis_sierk(j,bfis_s)
	REAL (KIND=r4), INTENT(IN) :: J
	REAL (KIND=r4), INTENT(OUT) :: bfis_s
	INTEGER :: i,k

	 bfis_s = 999.0
	 IF (z > 102.0 .AND. J > 0.0) THEN
       WRITE (UNIT=*, FMT=*) "z and J out of limits"
       RETURN
	 END IF

!bfis = 0
!DO i=1,7
!   DO k=1,7
!       bfis = bfis + elzcof(k,i)*pz(k)*pa(i)
!   END DO
!END DO

	 bfis = DOT_PRODUCT(pz,MATMUL(elzcof,pa))

 	 IF (J < 1) THEN
       bfis_s = bfis
	 RETURN
	 END IF

	J80 = 0.0
	J20 = 0.0

	 DO i=1,4
        DO k=1,5
       J80 = J80 + elmcof(k,i)*pz(k)*pa(i)
       J20 = J20 + emncof(k,i)*pz(k)*pa(i)
	  END DO
	 END DO

	 IF (J <= J20) THEN
	  q = 0.2/(J20**2*J80**2*(J20-J80))
	  qa = q*(4.0*J80**3-J20**3)
	  qb = -q*(4.0*J80**2-J20**2)
	  bfis = bfis*(1.0+qa*REAL(J,KIND=r8)**2.0+qb*REAL(J,KIND=r8)**3.0)
	 ELSE
	  x = J20/Jmax
	  y = J80/Jmax
	  aj = (-20.0*x**5+25.0*x**4-4.0)*(y-1.0)**2*y*y
	  ak = (-20.0*y**5+25.0*y**4-1.0)*(x-1.0)**2*x*x
	  q = 0.2/((y-x)*((1.0-x)*(1.0-y)*x*y)**2)
	  qa = q*(aj*y-ak*x)
	  qb = -q*(aj*(2.0*y+1.0)-ak*(2.0*x+1.0))
	  zzz = REAL(J,KIND=r8)/Jmax
	  a1 = 4.0*zzz**5-5.0*zzz**4 + 1.0
	  a2 = qa*(2.0*zzz+1.0)
	  bfis = bfis*(a1+(zzz-1.0)*(a2+qb*zzz)*zzz*zzz*(zzz-1.0))
	 END IF

	IF (bfis <= 0.0)  bfis = 0.0
	IF (J > Jmax)  bfis = 0.0
	bfis_s = bfis

	RETURN
	END SUBROUTINE bfis_sierk

	SUBROUTINE yrast_sierk(J,E_rot_s)
	REAL (KIND=r4), INTENT(IN) :: J
	REAL (KIND=r4), INTENT(OUT) :: E_rot_s
	INTEGER :: k,l,m

  	   E_rot_s = 999.0
	   IF (z > 102.0 .AND. J > 0.0) THEN
	   WRITE (UNIT=*, FMT=*) "z and J out of limits"
	   RETURN
	   END IF

	   IF (A < amin .OR. A > amax) THEN
	   WRITE (UNIT=*, FMT=*) "A is out of range"
	   RETURN
	   END IF

	  E_rot = 0.0
	  JJ = REAL(J,KIND=r8)/Jmax
	  CALL lpoly (jj,9,pl)

!srk Now calculate rotating ground-state energy
!      IF (J > Jmax) RETURN
	 DO k=1,5
	   DO l=1,7
         DO  m=1,5
           E_rot = E_rot + egscof(m,l,k)*pz(l)*pa(k)*pl(2*m-1)
         END DO
	   END DO
	 END DO

	IF (E_rot < 0.0) E_rot = 0.0
	E_rot_S = E_rot
	RETURN
	END SUBROUTINE yrast_sierk

	SUBROUTINE M_inertia_sierk(J,M_inertia_1,M_inertia_2,M_inertia_3)
	REAL (KIND=r4), INTENT(IN) :: J
	REAL (KIND=r4), INTENT(OUT), OPTIONAL :: m_inertia_1,m_inertia_2,
     &                m_inertia_3

	INTEGER :: l,k

	JJ = REAL(J,KIND=r8)/Jmax

	aizro = 0.0_r8
	ai70 = 0.0_r8
	aimax = 0.0_r8
	ai95 = 0.0_r8
	aimin = 0.0_r8
	bizro = 0.0_r8
	bi70 = 0.0_r8
	bimax = 0.0_r8
	bi95 = 0.0_r8
	aimax2 = 0.0_r8
	ai952 = 0.0_r8

	DO l=1,6
	  DO k = 1,5
      aizro  = aizro  + aizroc(l,k) *pz(l)*pa(k)
      ai70   = ai70   + ai70c(l,k)  *pz(l)*pa(k)
      ai95   = ai95   + ai95c(l,k)  *pz(l)*pa(k)
      aimax  = aimax  + aimaxc(l,k) *pz(l)*pa(k)
      ai952  = ai952  + ai952c(l,k) *pz(l)*pa(k)
      aimax2 = aimax2 + aimax2c(l,k)*pz(l)*pa(k)
	  END DO

	  DO  k = 1,4
      bizro = bizro + bizroc(l,k)*pz(l)*pa(k)
      bi70  = bi70  + bi70c(l,k) *pz(l)*pa(k)
      bi95  = bi95  + bi95c(l,k) *pz(l)*pa(k)
      bimax = bimax + bimaxc(l,k)*pz(l)*pa(k)
	  END DO
	END DO

	ff1 = 1.0
	ff2 = 0.0
	fg1 = 1.0
	fg2 = 0.0
	IF (Z > 70) THEN
	    aimaxh = 0.0
	    aimidh = 0.0
	 DO l=1,4
        DO k=1,4
            aimaxh = aimaxh + aimax3c(l,k)*pz(l)*pa(k)
            aimidh = aimidh + aimax4c(l,k)*pz(l)*pa(k)
        END DO
	 END DO
	    IF (z > 80) ff1 = 0.0
	    IF (z >= 80) fg1 = 0.0
	    IF (bimax > 0.95) fg1 = 0.0
	    IF (aimaxh > aimax) ff1 = 0.0
	    ff2 = 1.0 - ff1
	    fg2 = 1.0 - fg1
	    aimax = aimax*ff1 + ff2*aimaxh
	    aimax2 = aimax2*ff1 + ff2*aimidh
	END IF
	bimax = bimax*fg1 + aimidh*fg2
	saizro = MAX(aizro,0.0_r8)
	sai70 = MAX(ai70,0.0_r8)
	sai95 = MAX(ai95,0.0_r8)
	saimax = MAX(aimax,0.0_r8)
	sai952 = MAX(ai952,0.0_r8)
	simax2 = MAX(aimax2,0.0_r8)
	sbimax = MAX(bimax,0.0_r8)
	sbi70 = MAX(bi70,0.0_r8)
	sbi95 = MAX(bi95,0.0_r8)
	sbizro = MAX(bizro,0.0_r8)

	q1 = -3.148849569
	q2 =  4.465058752
	q3 = -1.316209183
	q4 =  2.26129233
	q5 = -4.94743352
	q6 =  2.68614119
	gam = - 20.0*LOG(ABS(saizro-sai95)/ABS(saizro-saimax))
	aa = q1*saizro + q2*sai70 + q3*sai95
	bb = q4*saizro + q5*sai70 + q6*sai95
	gam2 = - 20.0*LOG(ABS(saizro-sai952)/ABS(saizro-simax2))
	aa2 = q1*saizro + q2*sai70 + q3*sai952
	bb2 = q4*saizro + q5*sai70 + q6*sai952
	aa3 = q1*sbizro + q2*sbi70 + q3*sbi95
	bb3 = q4*sbizro + q5*sbi70 + q6*sbi95
	gam3 = 60.0
	alpha = pi*(jj - 0.7)
	beta = 5.0*pi*(jj - 0.9)

	silt = saizro + aa*jj**2.0 + bb*jj**4.0
	sjlt = sbizro + aa3*jj**2.0 + bb3*jj**4.0
	silt2 = saizro + aa2*jj**2.0 + bb2*jj**4.0

	IF (jj <= 0.70) THEN
	   saimin = sjlt
	   saimx = silt
	   saimid = silt2
	ELSE
	   sigt = saizro + (saimax-saizro)*EXP(gam*(jj-1.0))
	   sjgt = sbi95 + (sbimax-sbi95)*EXP(gam3*(jj-1.0))
	   sigt2 = saizro + (simax2-saizro)*EXP(gam2*(jj-1.0))
	   IF (jj <= 0.95) THEN
       f1 = silt*(COS(alpha))**2 + sigt*(SIN(alpha))**2
       saimx = f1
       f3 = sjlt*(COS(alpha))**2 + sjgt*(SIN(alpha))**2
       saimin = f3
       f1m = silt2*(COS(alpha))**2 + sigt2*(SIN(alpha))**2
       saimid = f1m
	   ELSE
       f2 = silt*(COS(beta))**2 + sigt*(SIN(beta))**2
       saimx = f2
       f4 = sjlt*(COS(beta))**2 + sjgt*(SIN(beta))**2
       saimin = f4
       f2m = silt2*(COS(beta))**2 + sigt2*(SIN(beta))**2
       saimid = f2m
	   END IF
	END IF
	IF (ff2 > 0.01 .AND. fg2 > 0.01) THEN
	   q1 = 4.001600640
	   q2 = 0.960784314
	   q3 = 2.040816327
	   aa3 =  q1*sai70 - q2*saimax - (1.0+q3)*saizro
	   bb3 = -q1*sai70 + (1.0+q2)*saimax + q3*saizro
	   aa4 =  q1*sai70 - q2*simax2 - (1.0+q3)*saizro
	   bb4 = -q1*sai70 + (1.0+q2)*simax2 + q3*saizro
	   saimx = saizro + aa3*jj**2.0 + bb3*jj**4.0
	   saimid = saizro + aa4*jj**2.0 + bb4*jj**4.0
	END IF
	IF (saimid > saimx) saimid = saimx
	saimin = MAX(saimin,0.0_r8)
	IF (PRESENT(M_inertia_1)) m_inertia_1 = saimx
	IF (PRESENT(M_inertia_2)) M_inertia_2 = saimid
	IF (PRESENT(M_inertia_3)) M_inertia_3 = saimin
	RETURN
	END SUBROUTINE M_inertia_sierk

	SUBROUTINE bs_sierk (j,bs_s)
	REAL (kind=r4), INTENT (IN) :: j
	REAL (kind=r4), INTENT (OUT) :: bs_s
	INTEGER :: iz,ia,ill
	REAL (KIND=r8) :: xa,bs

	xa = a/320.0_r8
	JJ = REAL(J,KIND=r8)/Jmax
	CALL lpoly (jj,9,pl)
	CALL lpoly (xa,5,paa)
	CALL lpoly (zz,8,pzz)

	bs = 0.0_r8
	DO iz = 1,8
	 DO ia = 1,5
       DO ill = 1,5
          bs = bs + b(ill,ia,iz)*pzz(iz)*paa(ia)*pl(2*ill-1)
       END DO
	 END DO
	END DO
	bs_s = REAL(bs)

	END SUBROUTINE bs_sierk
cc 
	SUBROUTINE lpoly (x,n,pl)
! FUNCTIONAL DESCRIPTION:
!srk this subroutine calculates the ordinary Legendre Polynomials of
!srk order 0 to n-1 of argument  x  and stores them in the vector
!srk pl.  They are calculated by recursion relation from the first two
!srk polynomials.
!
!srk written by A. J. Sierk   LANL  T-9  February,1984
!
!srk NOTE:  pl and x must be REAL (KIND=8) :: on 32-bit computers!
!
! DUMMY ARGUMENTS:
!
!    input: x (REAL (KIND=4) ::) - argument of Legendre Polynomial
!           n (INTEGER (kind=4) ::) - n-1 is highest order of polynomial
!    output: pl(1:20) (REAL (KIND=4) ::) - value of Legendre Polynomials
!

	INTEGER, INTENT(IN) :: n
	REAL (KIND=r8), INTENT(IN)  :: x
	REAL (KIND=r8), DIMENSION(:), INTENT(OUT) :: pl
	INTEGER :: i

	pl(1) = 1.0
	pl(2) = x
	DO i=3,n
	  pl(i) = (REAL(2*i-3,KIND=r8)*x*pl(i-1) 
     &   - REAL(i-2,KIND=r8)*pl(i-2))/REAL(i-1,KIND=r8)
	END DO

	RETURN
	END SUBROUTINE lpoly

      SUBROUTINE fisrot (Z,A,al,bf)
!     Cohen-Plasil-Swiatecki rotating liquid drop fission barrier
!     Provided by courtesy of Dr. F. Plasil, Oak Ridge National Lab.
!
      REAL (KIND=r4), INTENT(IN) :: a,z,al
      REAL (KIND=r4), INTENT(OUT) :: bf

      REAL (KIND=r4), PARAMETER, DIMENSION(6,11) :: x1b = RESHAPE (
     & (/0.28,.243,.221,.208,.195,.18,.211,.186,.17,.1506,.136,.12, 
     & 0.152,.131,.1155,.096,.0795,.0625,.09725,.0795,.065,.0506,.0375,
     & 0.0253,.05771,.0455,.03414,.0235,.014,.0065,.03325,.0235,.0153,
     & .0081, 0.001,.0,.01625,.009,.0032,.0,.0,.0,.0071,.0,.0,.0,.0,.0,
     & .0,.0, 0.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,0.,0.,0.,0.,0.,0./),
     & (/6,11/))

      REAL (KIND=r4), PARAMETER, DIMENSION(6,11) :: x2b = RESHAPE ( 
     & (/.18,.1695,.1515,.133,.1155,.0949,.1495,.1363,.1165,.099, 
     &.0815,.0594,.12,.1032,.0864,.0678,.0469,.028,.09,.0725,.0556,.037,
     & .019,.0057,.0625,.045,.0304,.016,.005,0.,.0406,.0264,.0151,.0052,
     & 0.,0.,.0253,.0144,.0027,0.,0.,0.,.0141,.006,0.,0.,0.,0.,.0065,
     & .0008, 0.,0.,0.,0.,.002,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./),
     & (/6,11/))


      REAL (KIND=r4), PARAMETER, DIMENSION(10,20) :: x3b = RESHAPE ( 
     & (/.0949,.0755,.0564,.0382,.0223,.0121,.00588,.00242,.00069, 
     & .0001,.0873,.0684,.049,.0306,.0162,.0074,.00267,.00055,0.,0.,
     & .0801, .061,.0418,.0235,.0108,.00373,.00071,0.,0.,0.,.073,.054,
     & .035, .0178,.0062,.00125,0.,0.,0.,0.,.0661,.047,.0284,.012,.0025,
     & 0.,0.,0., 
     & 0.,0.,.0594,.0404,.022,.0065,0.,0.,0.,0.,0.,0.,.0528,.034,.0159, 
     &.002,0.,0.,0.,0.,0.,0.,.0465,.0277,.01,0.,0.,0.,0.,0.,0.,0.,.0401, 
     & .0217,.0044,0.,0.,0.,0.,0.,0.,0.,.0339,.0158,.00024,0.,0.,0.,0., 
     & 0.,0.,0.,.028,.0106,0.,0.,0.,0.,0.,0.,0.,0.,.0219,.0064,0.,0.,0., 
     &0.,0.,0.,0.,0.,.0164,.0025,0.,0.,0.,0.,0.,0.,0.,0.,.0122,0.,0.,0., 
     & 0.,0.,0.,0.,0.,0.,.0085,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0057,0.,0., 
     &0.,0.,0.,0.,0.,0.,0.,.0035,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0016,0.,0., 
     &0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 
     &      0.,0.,0.,0.,0./),(/10,20/))

      REAL (KIND=r4) :: an,paren,eso,cx,bx,dx,by,cy,dy,b2,b1
      INTEGER :: ix,iy

      an=a-z
      paren=1.-1.7826*((an-z)/a)**2
      eso=17.9439*paren*a**.66667
      x=0.019655*z*(z/a)/paren
      y=1.9254*al**2/(paren*a**2.3333)
      ix=20.*x+.999
      cx=ix
      bx=20.*x+.999
      dx=bx-cx
      IF (x <= 0.25) THEN
         by=10.*y+.999
         IF (by.GT.9.) by=9.
         IF (by.LT.1.) by=1.
         iy=by
         cy=iy
         dy=by-cy
         b2=(x1b(ix+1,iy+1)-x1b(ix,iy+1))*dx+x1b(ix,iy+1)
         b1=(x1b(ix+1,iy)-x1b(ix,iy))*dx+x1b(ix,iy)
         bf=(b2-b1)*dy+b1
      ELSE IF (x <= .5) THEN
         by=20.*y+.999
         IF (by.GT.11.) by=11.
         IF (by.LT.1.) by=1.
         ix=ix-5
         iy=by
         cy=iy
         dy=by-cy
         b1=(x2b(ix+1,iy)-x2b(ix,iy))*dx+x2b(ix,iy)
         b2=(x2b(ix+1,iy+1)-x2b(ix,iy+1))*dx+x2b(ix,iy+1)
         bf=(b2-b1)*dy+b1
         x=MIN(x,0.95)
         ix=20.*x+.999
         ix=ix-10
         by=100.*y+.999
         IF (by.GT.19.) by=19.
         IF (by.LT.1.) by=1.
         iy=by
         cy=iy
         dy=by-cy
         b1=(x3b(ix+1,iy)-x3b(ix,iy))*dx+x3b(ix,iy)
         b2=(x3b(ix+1,iy+1)-x3b(ix,iy+1))*dx+x3b(ix,iy+1)
         bf=(b2-b1)*dy+b1
       END IF  
       bf=bf*eso
      RETURN
      END SUBROUTINE fisrot
	END MODULE sierk_stuff
!------------------------------end:   MODULE sierk_stuff-----------------------------------

!------------------------------begin:   MODULE sortsub-----------------------------------
	MODULE sortsub
	USE gem_kind
	IMPLICIT NONE
	PRIVATE

	 PUBLIC :: write_file_info, get_event, read_file_header, 
     &          open_event_file
	 PUBLIC :: start_sort, name

	INTEGER, PARAMETER, PUBLIC :: max_frag=100
	TYPE (identity), DIMENSION(max_frag), PUBLIC :: frag
	INTEGER, PUBLIC :: n_frag, n_binary
	REAL (KIND=r4), PUBLIC :: weight

	CONTAINS
	SUBROUTINE write_file_info(unit_out)
	INTEGER, INTENT(IN) :: unit_out

	WRITE (UNIT=unit_out, FMT="(TR2,A40)") para%reaction
	WRITE (UNIT=unit_out, FMT="(TR2,a19,a9)")
     &	"date of gemini run=",para%day

	WRITE (UNIT=unit_out, FMT="(TR2,A28,f10.3,a8)")
     &	"number of events per l-wave=",para%R_num,"*(2*L+1)"

	WRITE (UNIT=unit_out, FMT="(TR2,a22,i3,TR2,a4,i3,/TR2,a19,f7.2,
     &  TR1,a3)") "Compund nucleus:- Zcn=",para%Zcn, "Acn=",para%Acn,
     &  "Excitation energy =", para%Excn, "MeV"
	 IF (para%dexcn > 0) THEN
	  WRITE (UNIT=unit_out, FMT="(TR2,A31,F7.3)")
     &    "spread in excitation energy= +-",para%dexcn
	 END IF

	WRITE (UNIT=unit_out, FMT="(TR2,a47,i3,TR2,i3)")
     &	"minimum and maximum L waves simulated by Gemini="
     &    ,para%lmin, para%lmax

	 If(para%N_weight > 0)then
	 WRITE (UNIT=unit_out, FMT="(TR2,a24,I3,A)")
     & "Weighting was turned on for ",para%N_weight,"binary decays"
	 ELSE
	 WRITE (UNIT=unit_out, FMT="(TR2,A)")
     & "NO Weighting"
	 end if

	 IF (para%I_angle) THEN
	 WRITE (UNIT=unit_out, FMT="(TR2,a)")
     &	 "Angles and Velocities of all fragments calculated"
	 ELSE
	 WRITE (UNIT=unit_out, FMT="(TR2,a)")
     &    "No angles and velocities calculated"
	 END IF

	 IF (para%IMF_option == 1 .AND. para%kramer /= 0) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A,F8.3)")
     &	    "Kramers factor=",para%kramer
	 END IF

	 IF (para%t_delay > 1000.0 .and. para%ratio /= 1.0
     & .and. para%imf_option == 0) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A38)")
     & "evaporation from deformed system only"
	 ELSE IF (para%t_delay > 0.0) THEN
	  WRITE (UNIT=unit_out, FMT="(TR2,A17,TR1,F9.3,TR1,A24)")
     & "Fission delay of",para%t_delay,"x10-21 seconds included"
	  WRITE (UNIT=unit_out, FMT="(TR2,A30,TR1,F9.3,TR1,A24)")
     & "second moment of eta dependence=",para%sig_delay
     & ," included"
	  IF (para%sharp_delay) THEN
       WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &  "Fission and Imf emission totally suppressed during delay time"
	  ELSE
        WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &      "Decay width ramps up during delay time"
	  END IF
	 ELSE
	  WRITE (UNIT=unit_out, FMT="(TR2,A16)") "No fission delay"
	 END IF

	IF (para%imf_option == 0) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "No IMF emission, only evaporation"
	ELSE IF (para%imf_option == 1) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &  "evaporation and symmetric fission only"
	ELSE IF (para%IMF_option == 2) THEN
       WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "All asymmetric divisions possible"
	  ELSE IF (para%IMF_option == 3) THEN
       WRITE (UNIT=unit_out, FMT="(TR2,A)")  "All asymmetric divisions 
     & possible, but total fission yield determined from 
     & symmetry division"
	END IF

	 IF (para%imf_option == 1) THEN
	    if(para%b_scale <= 0.0) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &        "Sierks fission barriers used"
        ELSE
	  WRITE (UNIT=unit_out, FMT="(TR2,A,TR2,F5.2)")
     &  "Sierks fission barrier scaled by ",para%b_scale
	    END IF
	 ELSE IF (para%imf_option >= 2) THEN
	    if(para%b_scale <= 0.0) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &        "conditional saddle barriers used"
        ELSE
	  WRITE (UNIT=unit_out, FMT="(TR2,A,TR2,F5.2)")
     &  "conditional saddle barriers scaled by ",para%b_scale
	    END IF
	 END IF


	SELECT CASE (para%aden_type)
	CASE (0)
	  WRITE (UNIT=unit_out, FMT="(TR2,A23,F5.2)")
     &     "level density constant=",para%aden_0
	CASE (-1)
        WRITE (UNIT=unit_out, FMT="(TR2,A40)")
     &    "Toke + Swiatecki level density parameter"
	CASE (-2)
        WRITE (UNIT=unit_out, FMT="(TR2,A32)")
     &    "Ignatyuk level density parameter"
	CASE (-3)
	  WRITE (UNIT=unit_out, FMT="(TR2,A50)")
     &    "Gottschalk and Ledergerber level density parameter"
	CASE (-4)
	  WRITE (UNIT=unit_out, FMT="(TR2,A45)")
     &    "Temperature dependent level density parameter"
	CASE (-6)
	   WRITE (UNIT=unit_out, FMT="(TR2,A45)")
     &    "Lestones Temperature dependent level density parameter"
	CASE (-9)
         WRITE (UNIT=unit_out, FMT="(TR2,A,F5.2)")
     &   "Fineman's Temp dependent level-density parameter, kappa=",
     &   para%aden_0
      CASE (-22)
        WRITE (UNIT=unit_out, FMT="(TR2,A,F5.2,TR2,A,F5.2)")
     &    "U-dependent k0=",para%aden_k0,"kinf=",para%aden_kinf
	CASE DEFAULT
	  WRITE (UNIT=unit_out, FMT="(TR2,A,I5)")
     &  "aden_type=", para%aden_type
	END SELECT

	 IF (para%imf_option == 1)
     &WRITE (UNIT=unit_out, FMT="(TR2,A,F6.3)") "af/an=",para%a_scale


	 IF (para%tl_IWBC) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A29)")
     &  "Transmission coeff. from IWBC"
	ELSE
	    WRITE (UNIT=unit_out, FMT="(TR2,A35)")
     &    "Sharp cut-off Transmission coeff."
	END IF


	IF (para%ratio == 1.0 .OR. para%ratio == 0.0) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A,F5.2)")
     &   "spherical transmission coef ,ratio=",para%ratio
	 ELSE
	    WRITE (UNIT=unit_out, FMT="(TR2,A,F5.2)")
     &    "deformed transmission coef ,ratio=",para%ratio
	END IF

!IF (para%tl_iwbc) THEN
!    WRITE (UNIT=unit_out, FMT="(TR2,A,E8.2)")
!    "threshold probabilty for evaporation=",para%threshold
!END IF 


!IF (para%collective_enhancement) THEN
!    WRITE (UNIT=unit_out, FMT="(TR2,A)")
!    "rotationally-symmetric-deformed level densities"
!ELSE
!    WRITE (UNIT=unit_out, FMT="(TR2,A)")
!    "spherical level densities"
!END IF

	IF (para%lestone) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "Lestone modification to the transition state formulism"
	ELSE
	 WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &  "Standard transitional state formulism"
	END IF

	IF (para%quantum) THEN
        WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &  "Quantum treatment of spin projection"
	ELSE
	  WRITE (UNIT=unit_out, FMT="(TR2,A)")
     & "Semi-classical treatment of spin projection and emission angles"
	END IF

	IF (para%z_imf_min == 3) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)") "evaporation Z=0,1,2"
	ELSE IF (para%Z_imf_min == 4) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)") "evaporation Z=0,1,2,3"
	ELSE IF (para%Z_imf_min == 5) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)") "evaporation Z=0,1,2,3,4"
	END IF

	IF (para%mass_option == 0) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"liquid drop masses used"
	ELSE IF (para%Mass_Option == 1) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &   "experimental masses used"
	ELSE IF (para%Mass_Option == 2) THEN
  	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"experimental masses with shell fade-out used"
	END IF

	IF (para%polarization) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "Polarization of protons considered"
	END IF

	IF (para%exotic_index == 0) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "Only boring fragments evaporated"
	ELSE IF (para%exotic_index == 1) THEN
       	WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"Some mildly exciting fragments evaporated"
	ELSE IF (para%exotic_index == 2) THEN
		WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &    "Some real exotic fragments evaporated"
	ELSE 
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"Some real exotic fragments evaporated"
	END IF

	if (para%IMF_option >= 2) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A,TR2,F6.3)")
     &    "af/an=",para%a_scale
        	WRITE (UNIT=unit_out, FMT=*)
     &	"barrier scaled by",para%b_scale
	END IF 

	IF (para%time_flag) THEN
          WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"emission times are written in the event file"
	END IF

	IF (para%time_flag) THEN
	    WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &	"fragment spins written in the event file"
	END IF

	IF (para%k_sum) THEN
	 WRITE (UNIT=unit_out, FMT="(TR2,A)")
     &  "K and shape distributions considered"
	END IF 

	WRITE (UNIT=unit_out, FMT="(TR2,A,F7.2,A)")"E2 strength=",
     &	para%E2_strength, "W.u."

	WRITE (UNIT=unit_out, FMT="(TR2,A,F7.2,A)")"E1 strength=",
     &	para%E1_strength, "W.u."


	RETURN 
	END SUBROUTINE write_file_info

	SUBROUTINE get_event (l,I_read_error)

![subroutine_header_comments]
	INTEGER, INTENT(OUT) :: l
	INTEGER, INTENT(OUT) :: I_read_error
	INTEGER :: i
	INTEGER (KIND=2), DIMENSION(600) :: array
	INTEGER (KIND=2) :: n_record
	REAL (KIND=r4) :: partr4
	INTEGER :: ipt

	READ (UNIT=2, IOSTAT=I_read_error) n_record, (array(i),
     &	 i=1,n_record-1)
	IF (I_read_error /= 0) RETURN

	ipt = 1
	n_frag = array(ipt)
	ipt = ipt + 1
	n_binary = array(ipt)
	ipt = ipt + 1
	l = array(ipt)
	weight = TRANSFER(array(ipt+1:ipt+2),partr4)

	 ipt = ipt + 2

	 DO i=1,n_frag
	    frag(i)%A = array(i+ipt)
	    frag(i)%Z = IBITS(frag(i)%A,8,8)
	    IF (frag(i)%z == 255) frag(i)%z = -1
	    frag(i)%A = IBITS(frag(i)%A,0,8)
	 END DO
	 ipt = ipt + n_frag

	FORALL (i=1:n_frag) frag(i)%id = array(ipt+i)
	ipt = ipt + n_frag

	if (para%i_angle) THEN
		FORALL (i=1:n_frag)
          frag(i)%vel%v = TRANSFER(array(ipt+2*i-1:ipt+2*i),partr4)
	    END FORALL
	    ipt = ipt + 2*n_frag
		
	    FORALL (i=1:n_frag)
	    frag(i)%vel%theta = TRANSFER(array(ipt+2*i-1:ipt+2*i),partr4)
	    END FORALL
	    ipt = ipt + 2*n_frag

	    FORALL (i=1:n_frag)
          frag(i)%vel%phi = TRANSFER(array(ipt+2*i-1:ipt+2*i),partr4)
	    END FORALL
	    ipt = ipt + 2*n_frag
	END IF

	if (para%polarization) THEN
	   forall (i=1:n_frag)  frag(i)%M = REAL(array(i+ipt) )/2.0
	   ipt = ipt + n_frag
	END IF

	if (para%time_flag) THEN
	    FORALL (i=1:n_frag)
       frag(i)%time = TRANSFER(array(ipt+2*i-1:ipt+2*i),partr4)
	    END FORALL
	    ipt = ipt + 2*n_frag
	END IF 
	
	if (para%j_flag .OR. para%polarization) THEN
	   forall (i=1:n_frag)  frag(i)%J = REAL(array(i+ipt) )/2.0
	   ipt = ipt + n_frag
	END IF

	if (para%ex_flag) THEN
	   forall (i=1:n_frag)  frag(i)%exx = REAL(array(i+ipt) )/10.0
	   ipt = ipt + n_frag
	END IF

	RETURN  
	END SUBROUTINE get_event

	SUBROUTINE read_file_header(unit_in,I_read_error)

![subroutine_header_comments]
	INTEGER, INTENT(IN) :: unit_in
	INTEGER, INTENT(OUT) :: i_read_error
	REAL (KIND=r4), save :: r_num=0.0
	READ (UNIT=unit_in, IOSTAT=I_read_error) para
	para%r_num = para%r_num + r_num
	r_num = para%r_num
	
	RETURN 
	END SUBROUTINE read_file_header

	SUBROUTINE open_event_file(unit_in,I_read_error,unit_out)
![subroutine_header_comments]
	INTEGER, INTENT(IN) :: unit_in,unit_out
	INTEGER, INTENT(OUT) :: I_read_error
	INTEGER :: i
	CHARACTER(LEN=4), DIMENSION(3), PARAMETER :: 
     &	   EXT=(/".inp",".evt",".out"/)
	CHARACTER(LEN=30) :: FILE_name
	CHARACTER (LEN=34) :: evtfile,infile,outfile
	LOGICAL, SAVE :: be=.FALSE.
	
	WRITE (UNIT=*, FMT=*)"input file"

	READ (UNIT=5, FMT="(a)", IOSTAT=I_read_error) file_name
	IF (I_read_error /= 0) THEN
	   WRITE (UNIT=*, FMT=*)"WHAT!!!"
	   RETURN
	END IF

	file_name = ADJUSTL(file_name)    !remove leading blanks
	outfile = TRIM(file_name)//ext(3)
	evtfile = TRIM(file_name)//ext(2)

	WRITE (UNIT=*, FMT=*)"file=",evtfile

	OPEN (	UNIT=unit_in, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",  
     & STATUS = "OLD", position = "REWIND", ACTION = "READ", SHARED, 
     & FILE = evtfile, IOSTAT=I_read_error)   

	IF (I_read_error /= 0) THEN
	  WRITE (UNIT=*, FMT=*)"Could not open event file",evtfile
	  RETURN
	END IF

	IF (unit_out >= 1) THEN 
	   IF (.NOT. be)THEN       !if an outfile is already open, doen"t open a new one
	  OPEN (UNIT=unit_out, STATUS="NEW", ACTION="WRITE", FILE=outfile)
	     be = .TRUE.
	   END IF
	   WRITE (UNIT=unit_out,FMT=*)"evtfile=",evtfile
	END IF

	RETURN 
	END SUBROUTINE open_event_file

	SUBROUTINE start_sort(unit_in,unit_out,I_read_error)

![subroutine_header_comments]
	INTEGER, INTENT(IN) :: unit_in,unit_out
	INTEGER, INTENT(OUT) :: I_read_error

	CALL open_event_file(unit_in,I_read_error,unit_out)
	if (i_read_error /= 0)RETURN
	CALL read_file_header(unit_in,I_read_error)
	IF (i_read_error == 0) THEN
	    call write_file_info(6)
	    IF (unit_out >= 1) THEN
	    call write_file_info(unit_out)
	    END IF
	END IF
	RETURN
	END SUBROUTINE start_sort

	FUNCTION name(Z) RESULT(funct_out)
	INTEGER, INTENT(IN) :: z
	CHARACTER(LEN=2) :: funct_out
	CHARACTER(LEN=2), DIMENSION(103), PARAMETER :: element = (/
     & " H","He","Li","Be"," B"," C"," N"," O"," F","Ne",  
     & "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",  
     & "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn",  
     & "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr",  
     & "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",  
     & "Sb","Te"," I","Xe","Cs","Ba","La","Ce","Pr","Nd",  
     & "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",  
     & "Lu","Hf","Ta"," W","Re","Os","Ir","Pt","Au","Hg",  
     & "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",  
     & "Pa"," U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",  
     & "Md","No","Lw"/)
	 IF (z > 0 .AND. Z < 104) THEN
	    funct_out = element(Z)
	 ELSE IF (z == 0) THEN
	    funct_out = " n"
	 ELSE
	    funct_out = "  "
	  END IF
	RETURN
	END FUNCTION name
   
	END MODULE sortsub
!------------------------------end:   MODULE sortsub-----------------------------------

!------------------------------begin:   MODULE spectra-----------------------------------
	MODULE spectra
	USE gem_kind
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: hist_ek,hist_normalize,hist_zero,hist_zero_acc,
     &          hist_write,hist_read,hist_print

	INTEGER, PARAMETER, PRIVATE :: Nmax = 50       ! energy bins in histrogram
	INTEGER, PARAMETER, PRIVATE :: Nmodes = 3      ! number of spectra
	INTEGER, PARAMETER, PRIVATE :: jmax = 80       ! number of Jbins
	REAL (KIND=r4) , PARAMETER, PRIVATE :: deltaE = 1. ! histrogram bin size
	REAL (KIND=r4), PRIVATE, DIMENSION(Nmodes,0:Nmax) :: ek_spectra
	REAL (KIND=r4), PRIVATE, DIMENSION(Nmodes,0:jmax,0:Nmax) :: 
     &		ek_spectra_acc

	CONTAINS
!subroutine to increment spectra for each particle type
	   SUBROUTINE hist_ek(ek,de,mode,yield)
	   REAL (KIND=r4), INTENT(IN) :: ek,de,yield
	   INTEGER, INTENT(IN) :: mode
	   INTEGER :: imax,imin,i,Imode
	   REAL (KIND=r4) :: r_max,r_min

	   if (ek < 0.) THEN 
	     WRITE (*,*)"negative kinetic energy in hist_ek, ek= ",ek 
	     RETURN
	   END IF

	   Imode = mode
	   if (mode == 7) Imode = 3
	   if (Imode > Nmodes) return
	   if (Imode < 1) return

	   r_max = (ek+de/2.)/deltaE
	   imax = FLOOR(r_max)
	   r_min = (ek-de/2.)/deltaE
	   imin = FLOOR(r_min)
	   imin = MAX(0,imin)

	   IF (imin == imax) THEN
	      ek_spectra(Imode,imin) = ek_spectra(Imode,imin) + yield
	      RETURN
	   END IF 

	   DO i=imin,imax
	     if (i > Nmax) EXIT
	     IF (i == imin) THEN
	       ek_spectra(Imode,i) = ek_spectra(Imode,i) +
     &        (real(i)+deltaE-r_min)/de*yield
	     ELSE IF ( i == imax) THEN
       ek_spectra(Imode,i) = ek_spectra(Imode,i) + 
     &        (r_max - real(i))/de*yield
	     ELSE
       ek_spectra(Imode,i) = ek_spectra(Imode,i) + deltaE/de*yield
	     END IF 
	   END DO
   
	   RETURN
	   END SUBROUTINE hist_ek

!normalize and accumulate spectra
	   SUBROUTINE hist_normalize(gamma_total,j)
	   REAL (KIND=r4), INTENT(IN) :: gamma_total
	   INTEGER, INTENT(IN) :: j
	   INTEGER :: i

	   if (j < 0 .OR. j > jmax) RETURN

	   DO i=1,Nmodes
	      ek_spectra_acc(i,j,:) = ek_spectra_acc(i,j,:) +
     &       ek_spectra(i,:)/gamma_total
	   END DO
	   RETURN
	   END SUBROUTINE hist_normalize

!zero spectra
	   SUBROUTINE hist_zero
	   ek_spectra = 0.
	   RETURN
	   END SUBROUTINE hist_zero

!zero accumulated spectra
	   SUBROUTINE hist_zero_acc
	   ek_spectra_acc = 0.
	   RETURN
	   END SUBROUTINE hist_zero_acc

	   SUBROUTINE hist_write(filename)
	   CHARACTER (LEN=44), INTENT(IN) :: filename

       OPEN (UNIT=2, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",  
     &       ACTION = "WRITE", STATUS = "NEW", FILE = filename)

	   WRITE (UNIT=2) ek_spectra_acc
	   CLOSE (UNIT=2)

	   RETURN
	   END SUBROUTINE hist_write


	   SUBROUTINE hist_read(filename)
	   CHARACTER (LEN=44), INTENT(IN) :: filename

       OPEN (UNIT=2, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",  
     &	  ACTION = "READ", STATUS = "OLD", FILE = filename)

	   READ (UNIT=2) ek_spectra_acc
	   CLOSE (UNIT=2)

	   RETURN
	   END SUBROUTINE hist_read

	   SUBROUTINE hist_print(plb,lmin,lmax,dif)
	   real(kind=r4), INTENT(IN) :: plb      ! pi-lambdabar sqaured in mb
	   REAL (KIND=r4), INTENT(IN):: lmax     ! 
	   REAL (KIND=r4), INTENT(IN) :: dif     ! diffuseness
	   INTEGER, INTENT(IN) :: lmin 
	   INTEGER :: i,l
	   REAL(Kind=r4) :: total_n,total_p,total_a,fermi

	 DO i=0,Nmax
      total_n = 0.
      total_p = 0.
      total_a = 0.
      DO l=lmin,jmax
        IF (dif == 0) THEN
           if (l <= lmax) THEN
              fermi = 1.
           ELSE
              fermi = 0.
           END IF 
         ELSE 
            fermi = plb/(1.+exp((real(l)-lmax)/dif))
         END IF
        total_n = total_n + ek_spectra_acc(1,l,i)*fermi
        total_p = total_p + ek_spectra_acc(2,l,i)*fermi
        total_a = total_a + ek_spectra_acc(3,l,i)*fermi
      END DO
      ! make per steradian
      total_n = total_n/4./3.14159
      total_p = total_p/4./3.14159
      total_a = total_a/4./3.14159
      WRITE (UNIT=80, FMT=*) (REAL(i)+.5)*deltaE,total_n
      WRITE (UNIT=81, FMT=*) (REAL(i)+.5)*deltaE,total_p
      WRITE (UNIT=82, FMT=*) (REAL(i)+.5)*deltaE,total_a
       END DO 
	  RETURN
	  END SUBROUTINE hist_print
	END MODULE spectra
!------------------------------end:   MODULE spectra-----------------------------------

!------------------------------begin:   MODULE light_stuff-----------------------------------
	MODULE light_stuff
	USE gem_kind
	IMPLICIT NONE
	PRIVATE
	PUBLIC :: write_mult

	INTEGER, PARAMETER, PUBLIC :: N_mode = 51
	REAL (KIND=r4), DIMENSION(-2:N_mode), PUBLIC :: mult0, mult1
	REAL (KIND=r4), DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: Z1 = (/ 
     &  0.0,                        !4n
     &  0.0,                        !3n
     &  0.0,                        !2n
     &  0.0,                        !n
     &  1.0,                        !1H
     &  1.0,                        !2H
     &  1.0,                        !3H
     &  1.0,                        !4H
     &  2.0,                        !3He
     &  2.0,2.0,2.0,2.0,2.0,        !4He
     &  2.0,                        !5He
     &  2.0,2.0,                    !6He
     &  2.0,                        !7He
     &  2.0,2.0,                    !8He
     &  3.0,                        !5Li
     &  3.0,3.0,3.0,3.0,3.0,3.0,    !6li
     &  3.0,3.0,3.0,                !7Li
     &  3.0,3.0,3.0,                !8Li
     &  3.0,                        !9Li
     &  4.0,4.0,                    !6Be
     &  4.0,4.0,4.0,                !7Be
     &  4.0,4.0,                    !8Be
     &  4.0,4.0,4.0,4.0,4.0,4.0,4.0,!9Be
     &  4.0,4.0,4.0,4.0,            !10Be
     &  5.0,                        !10B
     &  6.0/)                        !12C

	REAL (KIND=r4), DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: A1 =(/
     &  4.0,                         !4n
     &  3.0,                         !3n
     &  2.0,                         !2n
     &  1.0,                         !n
     &  1.0,                         !1H
     &  2.0,                         !2H
     &  3.0,                         !3H
     &  4.0,                         !4H
     &  3.0,                         !3He
     &  4.0,4.0,4.0,4.0,4.0,         !4He
     &  5.0,                         !5He
     &  6.0,6.0,                     !6He
     &  7.0,                         !7He
     &  8.0,8.0,                     !8He
     &  5.0,                         !5Li
     &  6.0,6.0,6.0,6.0,6.0,6.0,     !6Li
     &  7.0,7.0,7.0,                 !7Li
     &  8.0,8.0,8.0,                 !8Li
     &  9.0,                         !9li
     &  6.0,6.0,                     !6Be
     &  7.0,7.0,7.0,                 !7Be
     &  8.0,8.0,                     !8Be
     &  9.0,9.0,9.0,9.0,9.0,9.0,9.0, !9Be
     &  10.0,10.0,10.0,10.0,         !10Be
     &  10.0,                        !10B
     &  12.0/)                       !12C

	REAL (KIND=r4), DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: S =(/
     &  0.0,                          !4n
     &  1.5,                          !3n
     &  0.0,                          !2n
     &  0.5,                          !n
     &  0.5,                          !1H
     &  1.0,                          !2H
     &  0.5,                          !3H
     &  2.0,                          !4H
     &  0.5,                          !3He
     &  0.0,0.0,0.0,2.0,2.0,          !4He
     &  1.5,                          !5He
     &  0.0,2.0,                      !6He
     &  1.5,                          !7He
     &  0.0,2.0,                      !8He
     &  1.5,                          !5Li
     &  1.0,3.0,0.0,2.0,2.0,1.0,      !6Li
     &  1.5,0.5,3.5,                  !7Li
     &  2.0,1.0,3.0,                  !8Li
     &  1.5,                          !9Li
     &  0.0,2.0,                      !6Be
     &  1.5,0.5,3.5,                  !7Be
     &  0.0,2.0,                      !8Be
     &  1.5,0.5,2.5,0.5,2.5,1.5,3.5,  !9Be
     &  0.0,2.0,2.0,1.0,              !10Be
     &  3.0,                          !10B
     &  0.0/)                         !12C

	REAL (KIND=r4), DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: def =(/
     &  32.28,                   !4n
     &  24.21,                   !3n
     &  16.14,                   !2n
     &  8.07,                    !n
     &  7.28,                    !1H
     &  13.13,                   !2H
     &  14.949,                  !3H
     &  25.86,                   !4H
     &  14.93,                   !3He
     &  2.424,22.63,23.434,24.264,25.75,   !4He
     &  11.3860,                 !5He
     &  17.5922,19.39,           !6He
     &  26.109,                  !7He
     &  31.598,34.198,           !8He
     &  11.6780,                 !5Li
     &  14.0854,16.368,17.648,18.40,19.46,19.73, !6Li
     &  14.9067,15.384,19.537,   !7Li
     &  20.9453,21.9261,23.200,  !8Li
     &  24.9538,                 !9Li
     &  18.3735,20.04,           !6Be
     &  15.7686,16.1977,20.338,  !7Be
     &  4.9416,7.9816,           !8Be
     &  11.3475,13.04,13.777,14.12,14.40,16.047,18.10, !9Be
     &  12.606,15.974,18.5643,18.565, !10Be
     &  12.051,                  !10B
     &  0.0/)                    !12C
    
	INTEGER, DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: level_no = (/
     &  1,                       !4n
     &  1,                       !3n
     &  1,                       !2n
     &  0,                       !n
     &  0,                       !1H
     &  0,                       !2H
     &  0,                       !3H
     &  1,                       !4H   ground state unstable
     &  0,                       !3He
     &  0,1,2,3,4,               !4He
     &  1,                       !5He  ground state unstable
     &  0,1,                     !6He
     &  1,                       !7He  ground state unstable
     &  0,1,                     !8He
     &  1,                       !5Li  ground state unstable
     &  0,1,2,3,4,5,             !6Li
     &  0,1,2,                   !7Li
     &  0,1,2,                   !8Li
     &  0,                       !9Li
     &  1,2,                     !6Be  ground state unstable
     &  0,1,2,                   !7Be
     &  1,2,                     !8Be  ground state unstable
     &  0,1,2,3,4,5,6,           !9be
     &  0,1,2,3,                 !10Be
     &  0,                       !10B
     &  0/)                      !12C

	INTEGER, DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: mode_gs = (/
     &  -2,                      !4n
     &  -1,                      !3n
     &  0,                       !2n
     &  1,                       !n
     &  2,                       !1H
     &  3,                       !2H
     &  4,                       !3H
     &  5,                       !4H   ground state unstable
     &  6,                       !3He
     &  7,7,7,7,7,               !4He
     &  8,                       !5He  ground state unstable
     &  9,9,                     !6He
     &  10,                      !7He  ground state unstable
     &  11,11,                   !8He
     &  12,                      !5Li  ground state unstable
     &  13,13,13,13,13,13,       !6Li
     &  14,14,14,                !7Li
     &  15,15,15,                !8Li
     &  16,                      !9Li
     &  17,17,                   !6Be  ground state unstable
     &  18,18,18,                !7Be
     &  19,19,                   !8Be  ground state unstable
     &  20,20,20,20,20,20,20,    !9be
     &  21,21,21,21,             !10Be
     &  22,                      !10B
     &  23/)                     !12C

	INTEGER, DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: exotic = (/
     &  3,                       !4n
     &  3,                       !3n
     &  3,                       !2n
     &  0,                       !n
     &  0,                       !1H
     &  1,                       !2H
     &  1,                       !3H
     &  3,                       !4H   ground state unstable
     &  1,                       !3He
     &  0,3,3,3,3,               !4He
     &  1,                       !5He  ground state unstable
     &  2,2,                     !6He
     &  2,                       !7He  ground state unstable
     &  2,2,                     !8He
     &  2,                       !5Li  ground state unstable
     &  0,0,0,0,0,0,             !6Li
     &  0,0,0,                   !7Li
     &  1,1,1,                   !8Li
     &  1,                       !9Li
     &  0,0,                     !6Be  ground state unstable
     &  0,0,0,                   !7Be
     &  0,0,                     !8Be  ground state unstable
     &  0,0,0,0,0,0,0,           !9be
     &  0,0,0,0,                 !10Be
     &  0,                       !10B
     &  0/)                      !12C

!INTEGER, DIMENSION(-2:N_mode), PARAMETER, PUBLIC :: exotic = (/
!3,                       &!4n
!3,                       &!3n
!3,                       &!2n
!0,                       &!n
!0,                       &!1H
!0,                       &!2H
!0,                       &!3H
!3,                       &!4H   ground state unstable
!0,                       &!3He
!0,3,3,3,3,               &!4He
!0,                       &!5He  ground state unstable
!0,3,                     &!6He
!0,                       &!7He  ground state unstable
!3,3,                     &!8He
!0,                       &!5Li  ground state unstable
!0,3,0,3,3,3,             &!6Li
!0,0,3,                   &!7Li
!3,3,3,                   &!8Li
!3,                       &!9Li
!3,3,                     &!6Be  ground state unstable
!0,0,3,                   &!7Be
!0,0,                     &!8Be  ground state unstable
!3,3,3,0,0,0,0,           &!9be
!3,3,3,3,                 &!10Be
!3,                       &!10B
!3/)                       !12C

!states with a nonzero level are unstable      
       REAL(KIND=r4), DIMENSION(-2:N_mode),PARAMETER,PUBLIC :: width=(/
     &  0.0,                                   !4n
     &  0.0,                                   !3n
     &  0.0,                                   !2n
     &  0.0,                                   !n
     &  0.0,                                   !1H
     &  0.0,                                   !2H
     &  0.0,                                   !3H
     &  0.0,                                   !4H   
     &  0.0,                                   !3He
     &  0.0,0.5,0.84,2.01,5.01,                !4he
     &  -0.60,                                 !5He  ground state unstable
     &  0.0,0.113,                             !6He
     &  0.16,                                  !7He
     &  0.0,0.8,                               !8He
     &  -1.5,                                  !5Li  ground state unstable
     &  0.0,0.024,0.0,1.7,0.54,1.5,            !6Li
     &  0.0,0.0,0.093,                         !7Li
     &  0.0,0.0,0.033,                         !8Li
     &  0.0,                                   !9Li
     &  0.096,1.16,                            !6Be  ground state unstable
     &  0.0,0.0,0.175,                         !7Be
     &  6.8e-6,-1.5,                           !8Be  ground state unstable
     &  0.0,0.15,7.7e-4,1.08,0.282,0.743,1.54, !9Be
     &  0.0,0.0,0.0,0.0,                       !10Be
     &  0.0,                                   !10B
     &  0.0/)                                  !12C

!for state with width /= 0, when choosing the
!energy of the state, there is a lower limit given by the
!q_value for the decay.
	REAL(KIND=r4), DIMENSION(-2:N_mode), PARAMETER,PUBLIC :: limit=(/
     &  0.0,                                     !4n
     &  0.0,                                     !3n
     &  0.0,                                     !2n
     &  0.0,                                     !n
     &  0.0,                                     !1H
     &  0.0,                                     !2H
     &  0.0,                                     !3H
     &  -2.9,                                    !4H
     &  0.0,                                     !3He
     &  0.0,-0.39,-1.19,-2.025,-3.5,             !4He
     &  -0.89,                                   !5He
     &  0.0,-0.975,                              !6He
     &  -0.444,                                  !7He
     &  0.0,-1.0,                                !8He
     &  -1.96,                                   !5Li
     &  0.0,-0.71,0.0,-2.83,-1.67,-4.17,         !6li
     &  0.0,0.0,-2.16,                           !7Li
     &  0.0,0.0,-0.22,                           !8Li
     &  0.0,                                     !9Li
     &  -1.371,-3.04,                            !6Be
     &  0.0,0.0,-2.9824,                         !7Be
     &  -0.0917,-3.13,                           !8Be
     &  0.0,0.0,-0.75,-1.12,-1.35,-3.034,-5.094, !9Be
     &  0.0,0.0,0.0,0.0,                         !10Be
     &  0.0,                                     !10B
     &  0.0/)                                    !12C

	REAL(KIND=r4), DIMENSION(100), SAVE, PUBLIC, TARGET :: he5 =
     &  (/7.9083718E-02, 0.4554695,     0.6319224,     0.8052081,   
     &   0.9870652,      1.184213,      1.403126,      1.649059,   
     &   1.928399,      2.246720,      2.609853,      3.021339,   
     &   3.480712,      3.987466,      4.518053,      5.046234,   
     &   5.525939,      5.903042,      6.116350,      6.132256,   
     &   5.951068,      5.605902,      5.151142,      4.645151,   
     &   4.132445,      3.645012,      3.200817,      2.806499,   
     &   2.461934,      2.163713,      1.907133,      1.686299,   
     &   1.496093,      1.332298,      1.190793,      1.068257,   
     &   0.9616849,     0.8686457,     0.7872887,     0.7156538,   
     &   0.6523586,     0.5963300,     0.5465656,     0.5021940,   
     &   0.4625382,     0.4269564,     0.3949688,     0.3661532,   
     &   0.3400630,     0.3163927,     0.2948750,     0.2753005,   
     &   0.2574341,     0.2411040,     0.2261390,     0.2124056,   
     &   0.1997746,     0.1881840,     0.1774403,     0.1675267,   
     &   0.1583473,     0.1498465,     0.1419592,     0.1346347,   
     &   0.1278281,     0.1214848,     0.1155677,     0.1100507,   
     &   0.1041386,     0.1002591,     9.5705785E-02, 9.1447316E-02,
     &   8.7451912E-02, 8.3710834E-02, 8.0197766E-02, 7.6900125E-02,
     &   7.3795922E-02, 7.0873909E-02, 6.8130739E-02, 6.5537177E-02,
     &   6.3092850E-02, 6.0781863E-02, 5.8603425E-02, 5.6547668E-02,
     &   5.4600596E-02, 5.2753977E-02, 5.1012747E-02, 4.9361911E-02,
     &   4.7792416E-02, 4.6307988E-02, 4.4893898E-02, 4.3564241E-02,
     &   4.2350385E-02, 4.1146707E-02, 3.9985135E-02, 3.8882740E-02,
     &   3.7839983E-02, 3.6844742E-02, 3.5895180E-02, 3.4998152E-02 /)

	REAL (KIND=r4), DIMENSION(100), SAVE, PUBLIC, TARGET :: li5 =     !li5 line shape
     &  (/1.4504250E-03, 3.8533162E-02, 9.6412919E-02, 0.1632309,   
     &  0.2356376,     0.3133758,     0.3971918,     0.4883010,    
     &  0.5881912,     0.6985758,     0.8214914,     0.9586236,    
     &   1.112319,      1.284578,      1.476477,      1.688190,    
     &   1.918175,      2.161553,      2.408900,      2.645893,    
     &   2.853843,      3.010771,      3.098557,      3.104661,    
     &   3.029880,      2.885232,      2.691240,      2.468492,    
     &   2.234132,      2.005620,      1.790282,      1.593327,    
     &   1.416697,      1.260246,      1.122674,      1.002306,    
     &  0.8972153,     0.8054150,     0.7252306,     0.6550402,    
     &  0.5934471,     0.5392347,     0.4916269,     0.4493487,    
     &  0.4118552,     0.3785081,     0.3487624,     0.3221603,    
     &  0.2982841,     0.2768121,     0.2574504,     0.2399428,    
     &  0.2241266,     0.2097653,     0.1966244,     0.1846383,    
     &  0.1736885,     0.1636661,     0.1544695,     0.1460203,    
     &  0.1382399,     0.1312022,     0.1244449,     0.1183150,    
     &  0.1126390,     0.1073720,     0.1024795,     9.7927205E-02,
     &  9.3707457E-02, 8.9770079E-02, 8.6072780E-02, 8.2613520E-02,
     &  7.9382330E-02, 7.6355636E-02, 7.3518194E-02, 7.0855774E-02,
     &  6.8356417E-02, 6.6001333E-02, 6.3788287E-02, 6.1706606E-02,
     &  5.9736867E-02, 5.7883032E-02, 5.6134474E-02, 5.4475829E-02,
     &  5.2908845E-02, 5.1428229E-02, 5.0023757E-02, 4.8627682E-02,
     &  4.7432180E-02, 4.6227124E-02, 4.5088600E-02, 4.4007234E-02,
     &  4.2974968E-02, 4.1992918E-02, 4.1058376E-02, 4.0168889E-02,
     &  3.9323967E-02, 3.8509648E-02, 3.7737995E-02, 3.7002672E-02 /)

	REAL (KIND=r4), DIMENSION(100), SAVE, PUBLIC, TARGET :: be8 = 
     & (/1.0039957E-09, 3.2721155E-06, 6.3356660E-05, 3.6410498E-04,
     &  1.1814651E-03, 2.8016411E-03, 5.4829111E-03, 9.4491532E-03,
     &  1.4896543E-02, 2.2000713E-02, 3.0931385E-02, 4.1849297E-02,
     &  5.4934327E-02, 7.0381209E-02, 8.8400558E-02, 0.1092571,   
     &  0.1332446,     0.1607178,     0.1921125,     0.2279434,   
     &  0.2688520,     0.3155786,     0.3690677,     0.4304372,   
     &  0.5010681,     0.5826480,     0.6772066,     0.7871007,   
     &  0.9152740,      1.065109,      1.240304,      1.444788,   
     &   1.681980,      1.953733,      2.258750,      2.589431,   
     &   2.928975,      3.248340,      3.509043,      3.668774,   
     &   3.698671,      3.593336,      3.374747,      3.081249,   
     &   2.754354,      2.426771,      2.119425,      1.842805,   
     &   1.600313,      1.391263,      1.212180,      1.059609,   
     &  0.9296771,     0.8189389,     0.7242035,     0.6429818,   
     &  0.5730698,     0.5129786,     0.4602561,     0.4146703,   
     &  0.3747103,     0.3396650,     0.3088302,     0.2815771,   
     &  0.2574029,     0.2358631,     0.2166860,     0.1995172,   
     &  0.1841138,     0.1702102,     0.1577836,     0.1461642,   
     &  0.1363510,     0.1271027,     0.1186560,     0.1109695,   
     &  0.1039595,     9.7560890E-02, 9.1707118E-02, 8.6336568E-02,
     &  8.1412978E-02, 7.6889671E-02, 7.2727934E-02, 6.8895988E-02,
     &  6.5399133E-02, 6.2097840E-02, 5.9084486E-02, 5.6295421E-02,
     &  5.3709153E-02, 5.1367849E-02, 4.9127694E-02, 4.7070879E-02,
     &  4.5098834E-02, 4.3350108E-02, 4.1727263E-02, 4.0129568E-02,
     &  3.8671881E-02, 3.7332911E-02, 3.6058463E-02, 3.4895122E-02 /)


	CONTAINS
	    SUBROUTINE write_mult (num_residue)
	    REAL (KIND=r4), INTENT(IN) :: num_residue
	    REAL (KIND=r4) :: mult
	    INTEGER :: l

	    if (num_residue <= 0.) RETURN
	    WRITE (UNIT=*, FMT=*) "Primary multiplicities"
	    if (para%tl_iwbc) THEN
         WRITE (UNIT=*, FMT=*) "num_residue=",num_residue
         WRITE (UNIT=*, FMT=*) "  evaportion mode,  mult(ER)"
         DO l=-2,n_mode
         if (mult0(l) > 0) WRITE (UNIT=*, FMT=*) l, mult0(l)/num_residue
         END DO
	    END IF 

	    WRITE (UNIT=*, FMT=*) "Secondary multiplicities"


    !neutron mult
	    mult = mult0(-2)*4. + mult0(-1)*3. + mult0(0)*2.
	    mult = mult + mult0(1) + mult0(5) + mult0(12) + mult0(14)*2. 
     &      + mult0(15) + mult0(17)*2. + mult0(23) + mult0(30) 
     &      + mult0(41) + mult0(42) + mult0(43) + mult0(44) + mult0(45)

	    WRITE (UNIT=*, FMT=*) "neutron=",mult/num_residue

    !proton
	    mult = mult0(2) + mult0(8) + mult0(9) + mult0(10) + mult0(11)
     &     + mult0(18) + mult0(23) + mult0(32)*2. + mult0(33)*2.
	    WRITE (UNIT=*, FMT=*) "proton=",mult/num_residue

    !deuteron
	    mult = mult0(3) + mult0(20) + mult0(22) + mult0(24)
	    WRITE (UNIT=*, FMT=*) "deuteron=",mult/num_residue          


    !triton
	    mult = mult0(4) + mult0(5) + mult0(27)
	    WRITE (UNIT=*, FMT=*) "triton=",mult/num_residue

    !3he
	    mult = mult0(6) + mult0(36)
	    WRITE (UNIT=*, FMT=*) "3he=",mult/num_residue

    !alpha
	  mult = mult0(7) + mult0(12) + mult0(14) + Mult0(18) + mult0(20)
     &     + mult0(22) + mult0(23) + mult0(24) + mult0(27) + mult0(32)
     &     + mult0(33) + mult0(36) + mult0(37)*2. + mult0(38)*2.
     &     + mult0(41)*2. + mult0(42)*2. + mult0(43)*2. + mult0(44)*2.
     &     + mult0(45)*2.
	    WRITE (UNIT=*, FMT=*) "alpha=",mult/num_residue

    !5he
	    WRITE (UNIT=*, FMT=*) "5He=",mult0(12)/num_residue

    !6he
	    mult = mult0(13) + mult0(15) + mult0(17)     
	    WRITE (UNIT=*, FMT=*) "6He=",mult/num_residue

    !7HE
	    mult = mult0(15) + mult0(17)
	    WRITE (UNIT=*, FMT=*) "7He=",mult/num_residue

    !5Li
	    WRITE (UNIT=*, FMT=*) "5Li=",mult0(18)/num_residue

    !6li
	    mult = mult0(19) + mult0(21)
	    WRITE (UNIT=*, FMT=*) "6Li=",mult/num_residue	
	    WRITE (UNIT=*, FMT=*) "6LI*=",mult0(20)/num_residue

    !7Li
	    mult = mult0(25) + mult0(26)
	    WRITE (UNIT=*, FMT=*) "7Li=",mult/num_residue
	    WRITE (UNIT=*, FMT=*) "7LI*=",mult0(27)/num_residue

    !8LI
	    mult = mult0(29) + mult0(29)
	    WRITE (UNIT=*, FMT=*) "8Li=",mult/num_residue

    !7Be
	    mult = mult0(34) + mult0(35)
	    WRITE (UNIT=*, FMT=*) "7Be=",mult/num_residue
	    WRITE (UNIT=*, FMT=*) "7Be*=",mult0(36)/num_residue

    !8Be
	  mult = mult0(37) + mult0(42) + mult0(43) + mult0(44) + mult0(45)
	    WRITE (UNIT=*, FMT=*) "8Be=",mult/num_residue
	    WRITE (UNIT=*, FMT=*) "8Be*=",mult0(38)/num_residue

    !9Be
	    mult = mult0(39) + mult0(40)
	    WRITE (UNIT=*, FMT=*) "9Be=",mult/num_residue
	    WRITE (UNIT=*, FMT=*) "9Be*=",mult0(41)/num_residue

    !10 Be
	    mult = mult0(46) + mult0(47) + mult0(48) + mult0(49)
	    WRITE (UNIT=*, FMT=*) "10Be=",mult/num_residue
	    END SUBROUTINE write_mult
	END MODULE light_stuff
!------------------------------end:   MODULE light_stuff-----------------------------------

!------------------------------begin:   MODULE mass_stuff-----------------------------------
	MODULE mass_stuff
	USE gem_kind
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: getmass, prepare_mass, prepare_mass_tf, massfr, 
     &          droplet_mass
	PUBLIC :: get_pairing

	INTEGER, PARAMETER, PRIVATE :: n_mass = 4200
	INTEGER, DIMENSION(0:100), PRIVATE :: chz
	REAL (KIND=r4), DIMENSION(n_mass), PRIVATE :: frmass,exmass,beta2

	INTEGER, DIMENSION(0:100), PARAMETER, PUBLIC :: cha1 = (/ 
     & 1,  1,  3,  5,  6,  8,  9, 11, 13, 14, 16, 18, 19, 22, 22,  
     & 24, 25, 27, 29, 30, 32, 34, 35, 37, 38, 40, 42, 44, 46, 48, 
     & 51, 53, 55, 57, 59, 61, 63, 66, 68, 70, 72, 74, 77, 79, 81, 
     & 83, 86, 88, 90, 92, 94, 97, 99, 101,103,106,108,110,113,115,
     & 118,120,123,125,128,130,133,136,138,141,143,146,149,151,154,
     & 156,159,162,165,167,170,173,175,179,183,184,190,193,196,199,
     & 203,206,206,207,212,222,226,233,230,231,232/)

	INTEGER, DIMENSION(0:100), PARAMETER, PUBLIC :: cha2 = (/  
     & 1,  3,  7, 10, 14, 17, 19, 21, 24, 27, 31, 34, 37, 39, 41,  
     & 44, 48, 51, 54, 57, 60, 64, 67, 70, 73, 76, 79, 82, 85, 88, 
     & 91, 94, 97,100,103,106,109,111,114,117,120,123,125,128,131, 
     & 133,136,139,141,144,146,149,151,154,156,159,161,164,166,168,
     & 171,173,175,177,180,182,184,186,188,191,193,195,197,199,201,
     & 203,205,207,209,211,213,210,214,215,218,219,226,229,230,232,
     & 236,238,240,241,246,247,251,251,255,256,258/)

	CONTAINS

	SUBROUTINE getmass(iZ,iA,mass_exp,mass_ld,beta_two)
! 
! FUNCTIONAL DESCRIPTION:       
! 
! This subroutine looks up the mass excess (MeV) of a nucleus (Z,A).
! 
! DUMMY ARGUMENTS:
! 
! input  Iz: INTEGER :: - proton number of nucleus
!        Ia: INTEGER :: - nucleon number of nucleus
!        type: logical =1 mass defect from Yakawa+exponential model
!                           =0 experimental value
!
! output mass: REAL (KIND=r4) ::   mass_exp = mass defect with pairing and d
!                                  shell correction(MeV)
!                                  mass_ld = mass defect from liquid drop
!                                  term only
!      only experimental mass for Z < 11 are included in table at
!      present
! 

	INTEGER, INTENT(IN) :: iZ,iA
	REAL (KIND=r4), INTENT(OUT), OPTIONAL :: mass_ld,mass_exp,beta_two
	INTEGER :: ind
	REAL (KIND=r4) :: m_ld, m_exp, beta22

	 beta22 = 0.

	 IF (iZ > 100) THEN
	   CALL massfr(REAL(iZ),REAL(iA),m_ld)
	   m_exp = m_ld
	 ELSE IF (iA >= cha1(iZ) .AND. iA <= cha2(iZ)) THEN
	   ind = chz(iz) - cha1(iz) + iA
	   if (ind > n_mass)THEN
       WRITE (UNIT=*, FMT=*) "increase n_mass in mass.f90"
       STOP
	 END IF 

	   m_ld = frmass(ind)
	   beta22 = beta2(ind)

	   IF (iZ <= 100) THEN
         m_Exp = exmass(ind)
	   ELSE
        WRITE (UNIT=*, FMT=*) "No exp mass for Z,A =",iZ,iA
         m_Exp = frmass(ind)
	   END IF
	  ELSE
	   CALL massfr(REAL(iZ),REAL(iA),m_ld)
	   m_exp = m_ld
	  END IF

	IF (PRESENT(mass_ld)) mass_ld = m_ld
	IF (PRESENT(mass_exp)) mass_exp = m_exp
	IF (PRESENT(beta_two)) beta_two = beta22

	RETURN
	END SUBROUTINE getmass

	SUBROUTINE prepare_mass
	INTEGER :: izz,I_READ_ERROR,ind
	INTEGER :: Inn,iaa
	REAL (KIND=r4) :: x1,x2,x3,x4,x5,x6,x7,x8,x9,m1,m2,m3,m4,m5,m6
	REAL (KIND=r4) :: pairing

! the moller table starts at Z=8, to get
! lighter elements, use the thomas-fermi table
	CALL prepare_mass_tf

	 chz(0) = 1
	 DO izz=1,100
	 chz(izz) = chz(izz-1) + cha2(izz-1) - cha1(Izz-1) + 1
	 END DO

	OPEN ( UNIT=56, FILE = "tbl/mass.tbl", STATUS = "OLD", ACTION= "READ")
C	 WRITE (*,*) "opened mass.tbl"

	 DO
        READ (UNIT=56, FMT="(3I5,9F10.3,6F10.2)", IOSTAT=I_READ_ERROR)
     &   IZZ,Inn,IAA,X1,X2,X3,X4,X5,X6,X7,X8,X9,M1,M2,M3,M4,M5,M6

	    IF (I_READ_ERROR /= 0) EXIT
	    IF (izz > 100) EXIT
	    IF (iAA >= cha1(iZZ) .AND. iAA <= cha2(iZZ)) THEN
          ind = chz(izz) - cha1(izz) + iAA

          IF (ind > n_mass)THEN
             WRITE (UNIT=*, FMT=*) "increase n_mass in mass.f90"
             STOP
          END IF 

        CALL get_pairing(izz,iaa,pairing)
             

        frmass(ind) = m2 - m1 - pairing
           
        exmass(ind) = m2

        beta2(ind) = x6

!       IF (m3 /= 0) THEN
!            exmass(ind) = m3
!       ELSE
!            exmass(ind) = m2
!       END IF

	    END IF 
	 END DO
	 CLOSE (UNIT=56)
	 END SUBROUTINE prepare_mass

	 SUBROUTINE prepare_mass_tf
	 INTEGER :: i, i_READ_Error,ind
	 INTEGER :: Inn,izz,iaa
	 REAL (KIND=r4) :: x1,x2,x3,x4,x5,m1,m2
	 REAL (KIND=r4) :: shell,pairing
	 LOGICAL :: be

	 INQUIRE (FILE = "gdat/mass_tf.gdat",  EXIST = be)

!open file and read plm2 numbers
	 if (be) THEN
	   OPEN (UNIT=2, FILE = "gdat/mass_tf.gdat", STATUS = "OLD",    
     &         ACTION = "READ", FORM = "UNFORMATTED")

	   READ (UNIT=2) chz,frmass,exmass
	   CLOSE (UNIT=2)
C	   WRITE (*,*) "opened mass_tf.gdat"
	   RETURN
	 END IF 

	chz(0) = 1
	 DO izz=1,100
	    chz(izz) = chz(izz-1) + cha2(izz-1) - cha1(Izz-1) + 1
	 END DO

	 OPEN (UNIT=56, FILE="tbl/mass_tf.tbl", STATUS="OLD", ACTION="READ")
C	 WRITE (*,*) "opened mass_tf.tbl"
! skip introductory remarks in file
	 DO i=1,66
	    READ (UNIT=56, FMT=*)
	 END DO

	 DO
	    READ (UNIT=56,
     &    FMT="(I3,I4,I4,F5.1,F5.1,F8.2,F6.2,F5.2,F6.2,F8.2,F7.2,F7.2)" 
     &    ,IOSTAT=I_READ_ERROR)                                         
     &    IZZ,Inn,IAA,x1,x2,x3,shell,pairing,x4,x5,M1,m2              

	    IF (I_READ_ERROR /= 0) EXIT
	    IF (izz > 100) EXIT
	    IF (iAA >= cha1(iZZ) .AND. iAA <= cha2(iZZ)) THEN
          ind = chz(izz) - cha1(izz) + iAA

       IF (ind > n_mass)THEN
           WRITE (UNIT=*, FMT=*) "increase n_mass in mass.f90"
           STOP
       END IF 

       CALL get_pairing(izz,iaa,pairing)
 

       frmass(ind) = m1 - shell - pairing

!      if you want to backshift level densities by the
!      congruence then uncomment the following line
!      and the liune in the GET_PAIRING subroutine      
!      frmass(ind) = m1 - shell - pairing - x4

        exmass(ind) = m1

!       IF (m2 /= 0) THEN
!           exmass(ind) = m2
!       ELSE
!           exmass(ind) = m1
!       END IF
	    END IF 
	 END DO
	 CLOSE (UNIT=56)

!having got all the masses, put them into a file
	 OPEN (UNIT=2, FILE = "gdat/mass_tf.gdat", STATUS = "NEW", 
     &      ACTION = "WRITE", FORM = "UNFORMATTED")
 
	WRITE (UNIT=2) chz,frmass,exmass
	CLOSE (UNIT=2)
	WRITE (*,*) "created mass_tf.gdat"
	END SUBROUTINE prepare_mass_tf

	 PURE SUBROUTINE massfr(Z,AA,mass)       
! 
! FUNCTIONAL DESCRIPTION:       
! 
! Calculates macroscopic finite range model masses of spherical
! nucleus using formula of Krappe, Nix, and Sierk 
! Reference- (Phys Rev C20, 992 (1979))
! modified to use the parameters of Moller + Nix Nucl. Phys. A361(1981)
! 117. Pairing correction term for odd-odd nuclei
! is included, as this is the most appropriate ground state for hot nuclei
! where shell and pairing effects have washed out.
! 
! DUMMY ARGUMENTS:
! 
! input       Z : REAL (KIND=r4) :: number of protons
!             N : REAL (KIND=r4) :: number of neutrons
!     
! output  mass : REAL (KIND=r4) :: mass defect in MeV
! 
! 
	REAL (KIND=r4), INTENT(IN) :: Z,AA
	REAL (KIND=r4), INTENT(OUT) :: mass

	REAL (KIND=r4) :: N,A13,A23,Fiss,Enz,Evol,x,fact,Es,Ec,Ewig
	REAL (KIND=r4) :: ad,Ecd,akf,Ecpf,Ea0,Eca,deltau,deltal,Ep,I

	REAL (KIND=r4), PARAMETER :: mn = 8.071431
	REAL (KIND=r4), PARAMETER :: mh = 7.289034
	REAL (KIND=r4), PARAMETER :: e2 = 1.4399764
	REAL (KIND=r4), PARAMETER :: b = 0.99
	REAL (KIND=r4), PARAMETER :: w = 36.
	REAL (KIND=r4), PARAMETER :: ael = 1.433e-5
	REAL (KIND=r4), PARAMETER :: r0 = 1.16
	REAL (KIND=r4), PARAMETER :: a = 0.68
	REAL (KIND=r4), PARAMETER :: as = 21.13
	REAL (KIND=r4), PARAMETER :: ks = 2.3
	REAL (KIND=r4), PARAMETER :: av = 15.9937
	REAL (KIND=r4), PARAMETER :: kv = 1.927
	REAL (KIND=r4), PARAMETER :: rp = 0.8
	REAL (KIND=r4), PARAMETER :: c0 = 4.4
	REAL (KIND=r4), PARAMETER :: ca = 0.212

	IF(AA .LE. 0.) THEN
	mass=0.
	RETURN
	ENDIF

	N = AA - Z
	A13 = AA**(1./3.)
	A23 = AA**(2./3.)

! relative neutron excess
	I = (N-Z)/AA
	fiss = Z**2./A13

! neutron-proton terms
	Enz = mn*N + mh*Z

! Volume energy
	Evol = -av*(1.-kv*I**2.)*AA

! Surface energy
	x=a/(r0*A13)
	fact=1.-3.*x**2. + (1.+1./x)*(2.+3.*x+3.*x**2.)*exp(-2./x)
	Es = as*(1.-ks*I**2.)*A23*fact

! Coulomb energy
	fact = fiss-0.76361*Z**1.3333/A13
	Ec = 0.6*e2/r0*fact

! Wigner term
	Ewig = w*abs(I)-ael*Z**2.39

! correction to Coulomb energy for diffuse surface
! see Davies & Nix Phys. Rev. C14 (1976) 1977
	ad=0.7071*b
	x = ad/(r0*A13)
	fact= 1 - 1.875*x+2.625*x**3.
     &      -.75*exp(-2./x)*(1.+4.5*x+7.*x**2.+3.5*x**3.)
	Ecd = -3.*Z**2.*e2*ad**2./(r0*a13)**3.*fact

! correction to coulomb energy due to proton form factor
	akf=1./r0*(7.06858*z/AA)**.33333
	x=(rp*akf)**2.
	Ecpf=-0.125*rp**2.*e2/r0**3. 
     &     *(3.0208-0.113541667*x+0.0012624*x**2.)*Z**2./AA

!A0 term
	Ea0=c0

! Charge asymmetry term
	Eca=ca*(Z**2.-N**2.)/AA

! pairing term for odd-odd nuclei
	deltau = 12./sqrt(AA)
	deltal = 20./AA
	Ep =deltau - 0.5 * deltal

! add all terms
	mass = Enz+Evol+Es+Ec+Ewig+Ecd+Ecpf+Ea0+Eca+Ep

	RETURN
	END SUBROUTINE massfr

	PURE SUBROUTINE droplet_mass(Z,A,mass)
 
! 
! gives the Droplet model masses from the work of W.D. Myers
! "Droplet model of atomic nuclei" (IFI, Plenum, New York, 1977)

	REAL (KIND=r4), INTENT(IN) :: Z,A
	REAL (KIND=r4), INTENT(OUT) :: mass

	REAL (KIND=r4) :: N,delta,epsilons,A3,A23,A43
	REAL (KIND=r4) :: Ev,Es,Ec

	REAL (KIND=r4), PARAMETER :: a1=15.96
	REAL (KIND=r4), PARAMETER :: a2=20.69
	REAL (KIND=r4), PARAMETER :: J=36.8
	REAL (KIND=r4), PARAMETER :: c1=0.73219
	REAL (KIND=r4), PARAMETER :: Q=17.
	REAL (KIND=r4), PARAMETER :: K=240.
	REAL (KIND=r4), PARAMETER :: L=100.
	REAL (KIND=r4), PARAMETER :: c2=0.00016302
	REAL (KIND=r4), PARAMETER :: c3=1.28846
	REAL (KIND=r4), PARAMETER :: c4=.55911
	REAL (KIND=r4), PARAMETER :: c5=.00049274

	N = A - Z
	A3 = A**0.33333333
	A23 = A**0.6666666
	A43 = A**1.3333333
	delta = ((N-Z)/A+3./16.*c1/Q*Z/A23)/(1.+9./4.*J/Q/A3)
	epsilons = (-2*a2/A3+L*delta**2.+c1*Z**2./A43)/K
	Ev = (-a1+j*delta**2.-0.5*K*epsilons**2.)*A
	Es = (a2+9./4.*J**2./Q*delta**2.)*A23 
	Ec =  c1*Z**2./A3 - c2*Z**2.*A3 - c3*Z**2./A 
     &    - c4/2**0.33333*Z - c5*Z**2.
	mass = Ev + Es + Ec

! Wigner term
	IF (NINT(n) == NINT(Z) .AND. MODULO(INT(N),2) == 1) THEN
	    mass = mass + 30./A
	ELSE
	    mass = mass + 30.* ABS(N-Z)/A
	END IF

! odd-even
! IF (MODULO(INT(Z),2) == 0 .AND. MODULO(INT(N),2) == 0) THEN
!   mass = mass - 12./A**0.5 +10/A
!ELSE IF (MODULO(INT(Z),2) == 1 .AND. MODULO(INT(N),2) == 1) THEN
!   mass = mass  + 12./A**0.5 - 10/A
!ELSE
!   mass = mass + 10/A
!END IF

	mass = mass + N*8.07169 + Z*7.28922

! remove atomic binding
! mass = mass - 0.00001433*Z**2.39
	RETURN
	END SUBROUTINE droplet_mass

	PURE SUBROUTINE get_pairing(iz,ia,pairing)
  
![subroutine_header_comments]
	INTEGER, INTENT(IN) :: iz,ia
	REAL (KIND=r4), INTENT(OUT) :: pairing

	INTEGER :: inn,oez,oen

	inn = ia - iz
	oez = MODULO(iz,2)
	oen = MODULO(inn,2)
	if (inn == iz .AND. oez == 1) THEN
	   pairing = 4.8/REAL(inn)**0.33333333 + 4.8/REAL(iz)**0.33333333 
     &   - 6.6/REAL(ia)**0.6666666 + 30/REAL(Ia)
	ELSE IF (oez == 1 .AND. oen == 1) THEN
	   pairing = 4.8/REAL(inn)**0.33333333 + 4.8/REAL(iz)**0.33333333         
     &   - 6.6/REAL(ia)**0.6666666
	ELSE IF (oez == 1 .AND. oen == 0) THEN
	   pairing = 4.8/REAL(iz)**0.33333333
	ELSE IF (oez == 0 .AND. oen == 1) THEN
	   pairing = 4.8/REAL(inn)**0.33333333
	ELSE
	   pairing = 0                                       
	END IF 

! zero pairing for odd-odd
! pairing = pairing - 4.8/REAL(inn)**0.33333333 &
! - 4.8/REAL(iz)**0.33333333                   & 
! + 6.6/REAL(ia)**0.6666666

! congruence energy
! if you want to backshift by the congruence energy then
! uncomment the follow line
! pairing = pairing  - 10.*EXP(-4.2*REAL(ABS(iz-inn))/REAL(ia))
  
	RETURN 
	END SUBROUTINE get_pairing

	END MODULE mass_stuff
!------------------------------end:   MODULE mass_stuff-----------------------------------

!------------------------------begin:   MODULE functs-----------------------------------
	MODULE functs
	USE gem_kind 
	USE trig
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: gauss, cl_go, factorial, exp_decay, transform, ran_j
	PUBLIC :: get_little_a, OPERATOR(+), zbeta
	PRIVATE :: addvectors

	INTERFACE OPERATOR (+)
	    MODULE PROCEDURE addvectors
	END INTERFACE

	REAL (KIND=r8), DIMENSION(0:300), PUBLIC  :: FCT
	REAL (KIND=r8), PARAMETER, PRIVATE  :: fact=50.0

	CONTAINS
	    FUNCTION gauss () RESULT (FUNCT_out)

    !+
    ! 
    ! FUNCTIONAL DESCRIPTION:
    ! 
    ! gaussian random number generator centered
    ! on zeroed with standard deviation
    ! of unity
    ! 

	    REAL (kind=r4) :: funct_out
	    REAL (kind=r4), SAVE :: x=0.0,theta=0.0
	    REAL (KIND=r4) :: ran
	    LOGICAL, SAVE  :: first = .true.
    
	    IF (first) THEN
        CALL RANDOM_NUMBER(ran)
        x = (-2.0*LOG(RAN+1.0e-37))**0.5
        CALL RANDOM_NUMBER(ran)
        theta = 360.0*RAN
        funct_out = x*COSD(theta)
        first = .FALSE.
	    ELSE
        funct_out = x*SIND(theta)
        first = .TRUE.
	    END IF
	    END FUNCTION gauss

	    PURE FUNCTION cl_go(J1,J2,J,M1,M2,M) RESULT(funct_out)

    !calculates clebsch Gordan coeeficients
    ! J's and M's must be twice the real values

	    INTEGER, INTENT(IN) :: J1,J2,J,M1,M2,M
	    REAL (kind=r4) :: funct_out

	    INTEGER :: Ja,Jb,Jc,Jd,Je,Jf,Jg,Jh,Ji,Jj
	    INTEGER :: fase,zmin,zmax,z
	    REAL (KIND=r8) :: cg,f1,f2,cfg

	    funct_out = 0.0
	    IF (M1+M2 /= M) RETURN
	    IF (ABS(M1) > J1) RETURN
	    IF (ABS(M2) > J2) RETURN
	    IF (ABS(M) > J) RETURN
	    IF (J > J1+J2) RETURN
	    IF (J < ABS(J1-J2)) RETURN

	    IF (j == 0 .AND. J1 == 0 .AND. j2 == 0)THEN
        funct_out = 1.0
        RETURN
	    END IF

	    Zmin = MAX(0, -J+J2-M1, -J+J1+M2)
	    Zmax = MIN(J1+J2-J, J2+M2, J1-M1)

	    cg = 0.0
	    DO z=zmin,zmax,2
        Ja = Z/2 
        Jb = (J1 + J2 - J - Z)/2 
        Jc = (J1 - M1 - Z)/2 
        Jd = (J2 + M2 - Z)/2 
        Je = (J - J2 + M1 + Z)/2 
        Jf = (J - J1 - M2 + Z)/2 
        fase = ((-1)**(z/2))
        f2 = REAL(fase, KIND=r8)
        cg = cg + f2/(FCT(Ja)*FCT(Jb)*FCT(Jc)*FCT(Jd)*FCT(Je)*FCT(Jf))
	    END DO

	    Ja = (J1 + J2 - J)/2 
	    Jb = (J1 - J2 + J)/2 
	    Jc = (-J1 + J2 + J)/2 
	    Jd = (J1 + M1)/2 
	    Je = (J1 - M1)/2 
	    Jf = (J2 + M2)/2 
	    Jg = (J2 - M2)/2 
	    Jh = (J + M)/2 
	    Ji = (J - M)/2 
	    Jj = (J1 + J2 + J + 2)/2 
	    F1 = REAL(J + 1, KIND=r8)
	    cfg= F1*FCT(Ja)*FCT(Jb)*FCT(Jc)*FCT(Jd)*FCT(Je)*FCT(Jf)* 
     &	     FCT(Jg)*FCT(Jh)*FCT(Ji)/FCT(Jj)
		if(cfg .LT. 0.)  cfg=0.
	    cg = SQRT(cfg)*cg
	    funct_out = REAL(cg/fact**0.5)

	    RETURN
	    END FUNCTION cl_go

	    SUBROUTINE factorial ()
	    INTEGER :: i

	    FCT(0) = 1.0_r8
	    FCT(1) = 1.0_r8/fact
	    DO i=2,300
	        FCT(i) = FCT(i-1)*REAL(i, KIND=r8)/fact
	    END DO
	    RETURN
	    END SUBROUTINE factorial

	    FUNCTION exp_decay(decay_width) RESULT(funct_out)
!   for a decay width in MEV, returns in zs,
!   a random time appropriate for exponential decay law
	    REAL (kind=r4), INTENT(IN) :: decay_width
	    REAL (kind=r4) :: funct_out
	    REAL (kind=r4) :: ran
	
	    CALL RANDOM_NUMBER(ran)
	    funct_out =  - 0.65824*LOG(RAN+1.0e-37)/decay_width
	    RETURN
	    END FUNCTION exp_decay

	 PURE SUBROUTINE transform(theta1,phi1,theta2,phi2,theta3,phi3)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine performs a rotational transformation.
! theta1,phi1 are the coordinates (spherical angles in
! degrees)of a unit vector in the original coordinate systems.
! the z  axis is made to rotate in the phi=phi2 plane by an
! angle theta2. The coordinates of the vector in the new
! reference frame are theta3,phi3
!
!
! To calculate the angle between two vectors
! theta1,phi1 and theta2,phi2
! use CALL transform(theta1,phi1,-theta2,phi2,theta3,phi3)
! i.e. theta3 in the angle between them, note the negative sign
! on theta2
!
!      
!
! DUMMY ARGUMENTS:
!
! input: theta1,phi1 (REAL (KIND=r4) ::) theta and phi angles (degrees) of initial unit
!        vector
!
!      theta2,phi2 (REAL (KIND=r4) ::) rotational angles
!
! OUTPUT: theta3,phi3 (REAL (KIND=r4) ::) theta and phi angles(degrees) of unit vector
!            after rotational transform
!
	REAL (KIND=r4), INTENT(IN) :: theta2,phi2
	REAL (KIND=r4), INTENT(IN) :: theta1,phi1
	REAL (KIND=r4), INTENT(OUT) :: theta3,phi3
	REAL (KIND=r4) :: xp,yp,zp
	REAL (KIND=r4) :: xt,yt,zt

! rotate vector in x-y plane by -phi2
! and find cartesian coordinates
	xp = SIND(theta1)*COSD(phi1-phi2)
	yp = SIND(theta1)*SIND(phi1-phi2)
	zp = COSD(theta1)

! rotate vector in z-x plane by theta2
	zt = zp*COSD(theta2) - xp*SIND(theta2)
	xt = xp*COSD(theta2) + zp*SIND(theta2)
	yt = yp

! rotate vector in x-y plane back by +phi2 and find spherical coordinates
	IF (xt == 0.0 .AND. yt == 0.0) THEN
	   phi3 = 0
	ELSE
	   phi3 = ATAN2d(yt,xt) + phi2
	   if (phi3 > 360.0) phi3 = phi3 - 360.0
	   if (phi3 < 0.0) phi3 = phi3 + 360.0
	END IF 

	if (zt > 1.0) THEN
	  theta3 = 0
	ELSE IF(zt < -1.0) THEN
	  theta3 = 180
	ELSE
	   theta3 = ACOSD(zt)
	END IF
  
	RETURN
	END SUBROUTINE transform

	PURE FUNCTION addvectors (v1,v2) RESULT(v3)
!
! FUNCTIONAL DESCRIPTION:
!
! The subroutine adds together two vectors "1"+"2" specified in spherical
! polar coordinates also giving the answer "3" in polar coordinates
!
! DUMMY ARGUMENTS:
!
! input: r1,theta1,phi1 (type (vell)) magnitude, theta and phi angle (degrees) of
!                                vector "1"
!        r2,theta2,phi2 (type (vell)) - vector "2"
!
! output: r3,theta3,phi3 (type (vell) vector "3" = "1" + "2"

	TYPE (vell), INTENT(IN) :: v1,v2
	TYPE (vell) :: v3
	REAL (KIND=r4) :: x,y,z

	z = v1%v*COSD(v1%theta) + v2%v*COSD(v2%theta)
	x = v1%v*SIND(v1%theta)*COSD(v1%phi) + 
     &    v2%v*SIND(v2%theta)*COSD(v2%phi)
	y = v1%v*SIND(v1%theta)*SIND(v1%phi) + 
     &    v2%v*SIND(v2%theta)*SIND(v2%phi)

	v3%v = (x**2.0+y**2.0+z**2.0)**0.5

	IF (v3%v == 0.0) THEN
	   v3%theta = 0.0
	   v3%phi = 0.0
	ELSE
	   v3%theta = ACOSD(z/v3%v)
	   IF (y == 0.0 .AND. x == 0.0) THEN
         v3%phi = 0.0
	   ELSE
       v3%phi = ATAN2D(y,x)
       IF (v3%phi < 0.0) v3%phi = v3%phi + 360.0
         END IF
	END IF

	RETURN
	END FUNCTION addvectors

	SUBROUTINE ran_j(x,ac,U,aden,rua,xmax)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine chooses in the Monte Carlo fashion the angular momentum
! associated with one or other of the angular momentum bearing modes; i.e.
! twisting, tilting, bending etc.
!
! DUMMY ARGUMENTS:
!
! input: ac (REAL (KIND=r4) ::) energy of modes = ac*x**2./2. MeV
!        U (REAL (KIND=r4) ::) - thermal excitation energy at saddle-point
!      aden (REAL (KIND=r4) ::) - little "a" level density parameter MeV-1
!      rua (REAL (KIND=r4) ::) - 2.*(aden*U)**.5
!        xmax (REAL (KIND=r4) ::) - maximum value of x possible, i.e. for the tilting
!      mode the maximum projection is the total angular momentum
!
! output: x angular momentum (hbar)
!
! IMPLICIT INPUTS:
!
!    none
!
! IMPLICIT OUTPUTS:
!
!    none
!
!
! SIDE EFFECTS:
!
!    none
!
!
	REAL (KIND=r4) :: y1,y2,xmax1,xmax2,ran
	REAL (KIND=r4), INTENT(OUT) :: x
	REAL (KIND=r4), INTENT(IN) :: aden,rua,xmax,U,ac

	xmax1 = (2.0*U/ac)**0.5
	xmax2 = MIN(xmax,xmax1)

	DO
	   CALL RANDOM_NUMBER(ran)
	   x = xmax2*(1.0-2.0*RAN)
	   y1 = exp(2.0*(aden*(U-0.5*ac*x**2.0))**0.5-rua)
	   CALL RANDOM_NUMBER(y2)
	   IF (y2 < y1) EXIT
	END DO

	RETURN
	END subroutine ran_j

	 SUBROUTINE get_little_a(A,Z,J,U,surf,shell_corr,E_rot,aden,
     &    entropy, pair_energy,aa_scale)

! 
! calculates the level density parameter "litte a" aden
! aden   = little a
! aden_1 = little at zero temperature, if little a is
!          temperature dependent
! 
	REAL (KIND=r4), INTENT(IN) :: A,Z,J,surf,aa_scale,pair_energy,
     &     shell_corr,E_rot
	REAL (KIND=r4), INTENT(INOUT) :: U
	REAL (KIND=r4), INTENT(OUT) :: aden,entropy
	REAL (KIND=r4) :: aden_1,aden_2,us,facte,E_sym,kappa
	REAL (KIND=r4), PARAMETER :: E_fade = 18.52

	INTEGER :: i
!backshift excitation energy

	U =  U + pair_energy
	IF (para%mass_option == 0) THEN
	    U = U + shell_corr           !liquid drop mass is used
	ELSE IF (para%mass_option == 2) THEN

	    u = u + shell_corr*(1.0-EXP(-(U+E_rot)/e_fade))   !shell corrections
                                                     !fade out

!     u = u + shell_corr/(1.0+EXP(-(U+E_rot-55.)/5.))  !shell corrections
                                                     !fade out
!    u = u + shell_corr*(1.0-EXP(-U/e_fade))  !shell corrections
                                             !fade out
	END IF

	E_sym = 0.0   

	SELECT CASE (para%aden_type)
	CASE (0)
	    aden = A/para%Aden_0
	CASE (-1)  !Toke + Swiatecki - one sphere
	    aden = A/14.61*(1+4.0*surf/A**0.3333333)
	CASE (-2)  !Ignatyuk et al
	    aden = A*0.073 + 0.095*surf*A**0.6666666        
!    aden = A*0.073 + 0.095*A**0.6666666
	CASE (-3) !Gottschalk + Ledergerber
	    aden = A*0.0999373 - 0.045934*surf*A**0.6666666
	CASE (-4) ! Ormand et al temp dependent
	    aden_1 = A/8.0
	    aden = aden_1/1.4*(1.0+0.4*EXP(0.1-1.213*U/A))
	case (-5)
	    aden_1 = A/8.0*
     &	 ((Z/A)**0.33333333+((A-z)/A)**0.333333333)/2.0**0.3333333333
	    aden = aden_1/1.4*(1.0+0.4*EXP(0.1-1.213*U/A))
	CASE (-6) ! lestone   
!   aden_1 = (1.6*A+1.8*surf*A**0.666666666)/15.5
	    aden_1 = (1.6*A+1.8*A**0.666666666)/15.5
	    aden_2 = 0.5*A/15.5

!   no temperature dependent symmerty energy
!    us = (U-E_sym)/441.0/A**0.33333333
!    facte = 1.0-EXP(-Us*8.0)
!    us = Us*A/(aden_1-aden_2*facte)
!    aden = aden_1 - aden_2*(1.0-EXP(-us))

!   with temperature dependent symmetry energy
	   DO i=1,2
      us = (U-E_sym)/441/A**0.33333333
      facte = 1.0-EXP(-Us*8.0)
      us = Us*A/(aden_1-aden_2*facte)
      aden = aden_1 - aden_2*(1.0-EXP(-us))
      E_sym = 0.82247*(1./aden-1.0/aden_1)*(A-2*Z)**2.0
	   END DO
	   U = U - E_sym

	CASE (-7)
	    aden = A/11.5 + 0.117*surf*A**0.66666*EXP(-U/A*5/3.0)
	CASE (-8)
	    aden = A/15 + 0.34*surf*A**0.66666*EXP(-U/A*5/3.0)
	CASE (-9) !Fineman
	    aden = A/(8.0+para%aden_0*U/A)
	CASE (-10)
	     aden = A/(7.0 + para%aden_0*U/A)
	CASE (-11)
	     aden = A/(6.0 + para%aden_0*U/A)
	CASE (-12)
	     aden = A/(8.0+para%aden_0*(U/A)**2.)
	CASE (-13)
	     aden = (a-z)/para%aden_0
	CASE (-14)
	     aden = z/para%aden_0
	CASE (-15) !grimes CASE C
	    aden = A/8.7873*EXP(-0.0493*(z-zbeta(A))**2.)
	CASE (-16) !grimes case B
	    aden = A/9.009*EXP(-0.000641*(A-2*z)**2.)
	CASE (-17) !grimes case B
	    aden = A/12.024*EXP(0.000461*(A-2*z)**2.)
	CASE (-18) !grimes case B
!    aden = A/8*(1.-0.1*U/A+0.01*(A-2*z)*U/A)
!    aden = A/8*(1.-0.2*U/A+0.01*(A-2*z)*U/A)
!    aden = A/8*(1.-0.3*U/A+0.02*(A-2*z)*U/A)
!    aden = A/8*(1.-0.4*U/A+0.02*(A-2*z)*U/A)

!     aden = A/8.*(1-0.0125*(A-2*z)*exp(-u/A/0.31))
!     aden = A/10.*(1-0.02*(A-2*z)*exp(-u/A/0.31))
!     aden = A/12.*(1-0.02*(A-2*z)*exp(-u/A/0.31))
!     aden = A/12.-A/12.*0.03*(A-2*z)*exp(-u/A/0.31))
!     aden = A/12.-A/12.*0.04*(A-2*z)*exp(-u/A/0.31)
!     aden = A/12.-A/12.*0.04*(A-2*z)*exp(-u/A/0.51)
!     aden = A/12.-A/12.*0.04*(A-2*z)*exp(-u/A/0.21)

!     aden = A/12. -0.5*(A-2.*Z)*exp(-u/A/0.31) + 0.2*(A-2.*Z)
!     aden = A/14. -0.5*(A-2.*Z)*exp(-u/A/0.31) + 0.2*(A-2.*Z)
!     aden = A/14. -0.5*(A-2.*Z)*exp(-u/A/0.51) + 0.2*(A-2.*Z)
!     aden = A/16. -0.5*(A-2.*Z)*exp(-u/A/0.51) + 0.2*(A-2.*Z)
!     aden = A/16. -0.4*(A-2.*Z)*exp(-u/A/0.81) + 0.2*(A-2.*Z)


! good
!      aden = 0.075*A + 0.0625*A*EXP(-U/A/0.51) &
!      -0.5*(A-2.*Z)*exp(-u/A/0.51) + 0.2*(A-2.*Z)

!      aden = 0.075*A + 0.075*A*EXP(-U/A/0.51) &
!      -0.6*(A-2.*Z)*exp(-u/A/0.51) + 0.2*(A-2.*Z)

!fig25
!      aden = 0.069*A + 0.081*A*EXP(-U/A/0.6) &
!      -0.65*(A-2.*Z)*exp(-u/A/0.6) + 0.25*(A-2.*Z)

!      aden = 0.1005*A + 0.081*A*EXP(-U/A/0.4) &
!      -0.65*(A-2.*Z)*exp(-u/A/0.4) 


!fig27
!       kappa =  0.2*(Z-zbeta(A))
!       aden = A/(8.+ kappa*U/A)

!fig28
        kappa = 0.1*(z-zbeta(A))
!       kappa = 0.1*(z-zbeta(A))-.4

!fig29
!       kappa= 0.4

!       kappa = 8.2e-5*Z**2
!        kappa = 1.6e-4*z**2
!        kappa = 0.6 + (z-70.)*0.1
       aden = A/(8.+ kappa*(U/A)**2.)
	CASE (-19)      !grimes case B
	    aden = A**(2./3.)*(Z**(1./3)+(A-Z)**(1./3.))/2.**(2./3.)
     &           /para%aden_0
	CASE (-20) !my calculations
	    aden = 0.861 + 0.08361*A + (A-2.*Z)**2*(6.229e-2/A-8.23e-5)
	CASE (-21) !my calculations
	    aden = A/(8.+para%aden_0*U/A*J/50.)
C CSM: U-dependent k evolves from k0 at ground state to kinf at high U.
      CASE (-22) ! U-dependent k: k0 (gs) -> kinf (high U)
        aden_2 = EXP(-0.00493*EXP(0.0332*A)*U/A)
        aden_1 = para%aden_kinf -
     &    (para%aden_kinf-para%aden_k0)*aden_2
        aden = A/aden_1
	END SELECT

	aden = aden*aa_scale

	IF (U > 0.0) THEN
	    entropy = 2.0*(Aden*U)**0.5
	ELSE
	    entropy = 0.0
	END IF

	RETURN  
	END SUBROUTINE get_little_a

	ELEMENTAL FUNCTION zbeta(A) RESULT(funct_out)
	REAL (KIND=r4), INTENT(IN) :: A
	REAL (KIND=r4) :: funct_out 
	funct_out = A/(1.98+0.0155*A**(2./3.))
	RETURN
	END FUNCTION zbeta

	END MODULE functs
!------------------------------end:   MODULE functs-----------------------------------



!------------------------------begin:   MODULE tl_stuff-----------------------------------
	MODULE tl_stuff
	USE gem_kind
	USE trig
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: get_tl,get_tl_half, prepare_deformed_emission, read_tl
	PUBLIC :: clear_look_up, prepare_plm2, prepare_interpolation, 
     &          get_tll
	PUBLIC :: sle,inverse,prepare_tl
	PUBLIC :: DE_ARRAY_tl   !suj

	INTEGER, ALLOCATABLE, DIMENSION(:), PRIVATE :: tl1
	LOGICAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE :: look_up
	REAL (KIND=r4), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: tl_table

	INTEGER, PARAMETER, PRIVATE :: n = 3
	INTEGER, PRIVATE :: i,I_def,ix
	INTEGER, PRIVATE :: kk_max,kk
	INTEGER, PARAMETER, PRIVATE :: N_tl = 22

	REAL (KIND=r4),PRIVATE :: eek
	REAL (KIND=r4),PRIVATE :: ee_tl_up,ee_tl_down
	REAL (KIND=r4),PRIVATE :: r_max,eee_tl

	REAL (KIND=r4), DIMENSION(4), PRIVATE :: eetl
	REAL (KIND=r4), DIMENSION(3,3), PRIVATE :: mat2, m_inv_low, 
     &                m_inv_high
	REAL (KIND=r4), DIMENSION(3), PRIVATE :: coef_low, coef_high

!REAL (KIND=r4),PARAMETER,PRIVATE :: x1 = 0.693 
!REAL (KIND=r4),PARAMETER,PRIVATE :: x3 = 1.58 
!REAL (KIND=r4),PARAMETER,PRIVATE :: x4 = 2.08 

	REAL (KIND=r4),PARAMETER,PRIVATE :: x1 = 0.8
	REAL (KIND=r4),PARAMETER,PRIVATE :: x2 = 1.00
	REAL (KIND=r4),PARAMETER,PRIVATE :: x3 = 1.3 
	REAL (KIND=r4),PARAMETER,PRIVATE :: x4 = 1.7 
	REAL (KIND=r4),PARAMETER,PRIVATE :: d_beta_k = 1.0
	REAL (KIND=r4),PARAMETER,PRIVATE :: d_beta_d = 10.0

	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:),PRIVATE,TARGET :: 
     &     lm_weight,ds
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:),PRIVATE :: plm2
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:), PRIVATE :: factor
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:,:), PRIVATE,TARGET :: 
     &     ratio_r
	REAL (KIND=r4), POINTER, DIMENSION(:),PRIVATE :: ds0,beta0
	REAL (KIND=r4), POINTER, DIMENSION(:,:),PRIVATE :: ratio_r0,
     &     lm_weight0
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:),PRIVATE :: ee_tl,xi,tl_s
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:),PRIVATE,TARGET :: beta

	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:,:,:),PRIVATE :: tl_def,
     &     tl_pol
	REAL (KIND=r4), ALLOCATABLE,DIMENSION(:,:,:), PRIVATE :: tl_coef
	REAL (KIND=r4), DIMENSION(0:l_max_quantum,3), PRIVATE :: c1,c2,c3,
     &     c4

	INTEGER, PRIVATE :: n_points
	INTEGER, PRIVATE :: i_beta
	INTEGER, PRIVATE :: i_beta_max,mode0
	INTEGER, DIMENSION(10), PRIVATE :: n_points_b
	INTEGER, DIMENSION(0:l_max_quantum), PRIVATE :: n_points_l,Lbeta 

	REAL (KIND=r4), DIMENSION(0:L_max_quantum), PRIVATE :: d_beta
	REAL (KIND=r4), PRIVATE :: al1,al2
	REAL (KIND=r4), PRIVATE :: scale,expand0
	LOGICAL, PRIVATE :: deformed, expanded, modified

	CONTAINS

	SUBROUTINE get_tl(l,ek,tl,tl_up,tl_down)
  
	INTEGER, INTENT(IN) :: l
	REAL (KIND=r4), INTENT(IN) :: ek
	REAL (KIND=r4), DIMENSION(0:l_max_quantum), INTENT(OUT) :: tl
	REAL (KIND=r4), INTENT(OUT) :: tl_up,tl_down
	REAL (KIND=r4) :: tl_min
!check to see to transmission coefficients are in look_up table

	if (l > l_max_quantum) THEN
	     WRITE (*,*) "l>l_amx_quantum",l    !rjc
	END IF

	if (modified) THEN
	    ix = INT(ek+1)
	  if (ix > tl1(l)) THEN
        tl(0:l) = 1.0
        RETURN
	    ELSE if (ix <= 50 .AND. l <= l_max_quantum) THEN
         IF (look_up(ix,l)) THEN
            IF (para%k_sum) THEN
                kk = l*(l+1)/2 + 1                           
                tl(0:l) = tl_table(kk:kk+l,ix)
            ELSE
                tl(0) = tl_table(l+1,ix)
            END IF 
        RETURN
        END IF 
	    END IF 
	END IF
	eek = ek
	IF (Mode0 <= 1) eek = ek + 1.0
	eek = eek*scale

	IF (modified) THEN

! ******** ratio=x1 ***********************
	    call get_tll(c1(l,:),eek,eetl(1))

!***** ratio=x2 ***************************
	    call get_tll(c2(l,:),eek,eetl(2))

!************ ratio=x3 ********************
	    call get_tll(c3(l,:),eek,eetl(3))

!    prepare for quadratic interpolation
!    for ratio <1.3 low
!    for ratio >1.3 high

	     coef_low = MATMUL(m_inv_low,eetl(1:3))

!************ ratio=x4 ***********************
	    if (r_max > x3) THEN
        call get_tll(c4(l,:),eek,eetl(4))
        coef_high= MATMUL(m_inv_high,eetl(2:4))
	     ELSE
        eetl(4) = -7.0
	     END IF 

	    if (deformed) THEN
 
         if (para%k_sum) THEN
            i_beta = Lbeta(l)
         ELSE
            i_beta = 1
         END IF
    
         ALLOCATE (ee_tl(n_points_b(i_beta)))
         ALLOCATE (tl_s(n_points_b(i_beta)))

         ratio_r0 => ratio_r(1:n_points_b(i_beta),:,i_beta)
!         ds0 => ds(1:n_points_b(i_beta),i_beta)

!        quadratic interpolation

         WHERE (ratio_r0(:,2) > x3)
             ee_tl = MATMUL(ratio_r0,coef_high)
         ELSEWHERE
             ee_tl = MATMUL(ratio_r0,coef_low)
         END WHERE          


         WHERE (ee_tl > 65.0)
             tl_s = 0.0
         ELSEWHERE
             tl_s = 1.0/(1.0+exp(ee_tl))
         END WHERE


!         WRITE (UNIT=17, FMT=*) eetl
!
!         DO i=1,9
!           WRITE (UNIT=17, FMT=*) ratio_r0(i,2),ee_tl(i),tl_s(i)
!         END DO 

         if (para%k_sum) THEN
            kk = l*(l+1)/2 + 1                           
            lm_weight0 => lm_weight(1:n_points_b(i_beta),kk:kk+l)
            tl(0:l) = MATMUL(tl_s,lm_weight0(:,1:l+1))
            tl_min = MINVAL(tl(0:l))
         ELSE
            lm_weight0 => lm_weight(1:n_points_b(i_beta),1:1)
             tl(0:0) = MATMUL(tl_s,lm_weight0(:,1:1))
            tl_min = tl(0)
         END IF 

        DEALLOCATE (ee_tl)
        DEALLOCATE (tl_s)

	    ELSE  if (expanded) THEN    ! no deformation, just expanded
         IF (expand0 <= x3) THEN
             eee_tl = coef_low(1)*expand0**2 + coef_low(2)*expand0 
     &       + coef_low(3)
         ELSE
             eee_tl = coef_high(1)*expand0**2 + coef_high(2)*expand0 
     &       + coef_high(3)
         END IF              

        IF (eee_tl  > 65.0) THEN
           tl(0) = 0.0
        ELSE 
            tl(0) = 1.0/(1.0+exp(eee_tl))
        END IF
        tl_min = tl(0)
        if (para%k_sum .AND. l>0 ) tl(1:l) = tl(0)
	   END IF                             

!   load transmission coefficients in lookup table
	    if (tl_min >= 0.99) THEN
         tl1(l) = ix
	    ELSE IF (ix <= 50) THEN
         if (para%k_sum) THEN
           kk = l*(l+1)/2 + 1                           
           tl_table(kk:kk+l,ix) = tl(0:l)
         ELSE
           tl_table(l+1,ix) = tl(0)
         END IF 
         look_up(ix,l) = .TRUE.
	     END IF 


	ELSE  ! not deformed or expanded
	    call get_tll(c2(l,:),eek,eee_tl)
	    IF (eee_tl  > 65.0) THEN
        tl(0) = 0.0
	    ELSE 
        tl(0) = 1.0/(1.0+exp(eee_tl))
	    END IF

	    if (para%k_sum .AND. l>0) THEN
        tl(1:l) = tl(0)
        RETURN
	    END IF 

	  IF ((mode0 == 2 .or. mode0 == 1) .AND. para%polarization) THEN
        call get_tll(c1(l,:),eek,ee_tl_up)
        call get_tll(c3(l,:),eek,ee_tl_down)

        IF (ee_tl_up  > 65) THEN
            tl_up = 0.0
        ELSE 
            tl_up = 1.0/(1.0+exp(ee_tl_up))
        END IF

        IF (ee_tl_down  > 65) THEN
            tl_down = 0.0
        ELSE 
            tl_down = 1.0/(1.0+exp(ee_tl_down))
        END IF
	    END IF
	END IF
	RETURN
	END SUBROUTINE get_tl

	SUBROUTINE get_tl_half(mode,iz2,l_min,ek_min)
	INTEGER, INTENT(IN) :: mode,iz2, l_min
	REAL   (KIND=r4), INTENT(OUT) :: ek_min
	REAL (KIND=r4) :: c1,c2,c3
!   find barrier i.e. where tl = .5

	    IF (TL_COEF(mode,IZ2,7) == 0) THEN
        WRITE (UNIT=*, FMT=*) "mode=",mode,"IZ2=",IZ2
	    END IF

	    c1 = tl_coef(mode,iz2,7)
	    c2 = tl_coef(mode,iz2,1) + REAL(l_min+1)*tl_coef(mode,iz2,2)   
     &        +REAL(l_min+1)**2*tl_coef(mode,iz2,3)

	    c3 = tl_coef(mode,iz2,4) + REAL(l_min+1)*tl_coef(mode,iz2,5)
     &        + REAL(l_min+1)**2*tl_coef(mode,iz2,6)
	    ek_min = c2**2-4.0*c1*c3 

	    IF (c1 == 0) THEN
          ek_min = -c3/c2
	    ELSE IF (ek_min < 0) THEN  ! Tl never gets to .5, This should never
                             ! happen, but it occasionally does for
                             ! light daughter nuclei and very large
                             ! orbital angular momentum. Such large
                             ! angular momentum contribute little
                             ! to the decay width, so I don"t think this
                             ! is a problem, and this IF statement is
                             ! really here to stop the code from
                             ! crashing
        ek_min = (c3/c1)**0.5
	    ELSE
        ek_min = (-c2-ek_min**0.5)/(2.0*c1)      ! Tl=.5 at this kinetic
	    END IF
	    if (mode <= 1) ek_min = ek_min - 1
	    if (ek_min <= 1) ek_min = 1.0

	RETURN
	END SUBROUTINE get_tl_half
    

	SUBROUTINE prepare_deformed_emission (ratio,expand)
!get radius at angle beta to symmetry axis

	REAL (KIND=r4), INTENT(IN) :: ratio,expand
	INTEGER :: j_beta,ell,m

	deformed = .FALSE.
	expand0 = 0.0
	expanded = .FALSE.

	IF (expand /= 0.0 .AND. expand /= 1.) THEN
	   expand0 = expand
	   expanded = .TRUE. 
	END IF 

	if (ratio == 1.0 .OR. ratio == 0.0) THEN
	  r_max = expand
	  RETURN
	END IF 

	deformed = .TRUE.

	IF (ratio < 0.0) STOP !not possible

	if (para%k_sum) THEN
	  i_beta_max = 10
	ELSE
	  i_beta_max = 1
	END IF 
 
	DO j_beta=1,i_beta_max
	    if (n_points_b(j_beta) == 0) CYCLE    
	    ALLOCATE (xi(n_points_b(j_beta)))
	    ALLOCATE (factor(n_points_b(j_beta)))
	    beta0 => beta(1:n_points_b(j_beta),j_beta)
	    ds0 => ds(1:n_points_b(j_beta),j_beta)
	    ratio_r0 => ratio_r(1:n_points_b(j_beta),:,j_beta)

    !the angle xi is defined by
    !x = ratio**(2.0/3.0)*COS(xi)    = a*COS(xi)
    !y = 1.0/ratio**(1.0/3.0)*SIN(xi) = c*SIN(xi)

    
	    xi = ATAN(ratio*TAN(beta0))

    !radius
	 factor = (SIN(xi)**2 + (ratio*COS(xi))**2)**0.5/ratio**(1.0/3.0)

    !surface area
	    ds0 = 2.0*pi*ratio**(1.0/3.0)*SIN(xi)*COS(xi)**2
     & 	     /COS(beta0)**2*(COS(xi)**2+(ratio*SIN(xi))**2)**0.5

	    if (expanded) factor = factor*expand

	    ratio_r0(:,1) = factor**2
	    ratio_r0(:,2) = factor
!   ratio_r0(:,3) = 1.0 this is now done in prepare_interpolation

	    DEALLOCATE (xi)
	    DEALLOCATE (factor)
	END DO

	r_max = MAXVAL(ratio_r(:,2,1))

	if ( .NOT. para%k_sum) THEN
	    lm_weight(:,1) = ds(:,1)/SUM(ds)
	    RETURN
	END IF 

!neeed to calculate tl are a function of projection
!on symmetry axis
!now need to get |d(l,m,m")|**2 coefficients
!to get m projection 

	kk = 0
	DO  ell=0,l_max_quantum
	    n_points = n_points_l(ell)
	    DO m=0,ell
        kk = kk + 1
        lm_weight(1:n_points,kk) = plm2(1:n_points,kk)
     &        *ds(1:n_points,Lbeta(ell))
        lm_weight(1:n_points,kk) = lm_weight(1:n_points,kk)
     &        /SUM(lm_weight(1:n_points,kk))
	    END DO
	END DO

	RETURN
	END SUBROUTINE prepare_deformed_emission

	SUBROUTINE read_tl ()
	LOGICAL :: be
	INTEGER :: mode,iii,jj
!----- normal transmission coeff ----------

	ALLOCATE (tl_coef(0:n_tl,140,7))

	 INQUIRE (FILE = "gdat/tl.gdat", EXIST = be)

	if (be) THEN
	  OPEN (UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &    ACTION = "READ", FILE = "gdat/tl.gdat", STATUS = "OLD")
C	    WRITE (*,*) "opened tl.gdat"
	    READ (UNIT=8) tl_coef

	    CLOSE (UNIT=8)
	ELSE
	    OPEN (UNIT=8, FILE= "tl.tbl", STATUS = "OLD", ACTION = "READ")
C	    WRITE (*,*) "opened tl.tbl"

	    READ (UNIT=8, FMT=*)
	    DO mode=0,n_tl
       READ (UNIT=8, FMT=*) 
       DO jj=1,100
          READ (UNIT=8, FMT="(7(x,F10.3))") tl_coef(mode,jj,:)
       END DO
       READ (UNIT=8, FMT=*)
	    END DO 
	    close (UNIT=8)

	    OPEN (UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &  ACTION = "WRITE", FILE = "gdat/tl.gdat",  STATUS = "NEW")
    
	    WRITE (UNIT=8) tl_coef

	    CLOSE (UNIT=8)
	    WRITE (*,*) "created tl.gdat"
	 END IF 

!-----  transmission coeff for deformed nuclei------

	IF (para%ratio /= 1 .OR. para%K_sum .or. para%expand) THEN
	    ALLOCATE (tl_def(n_tl,100,7,3))

	    INQUIRE (FILE = "gdat/tl_def.gdat", EXIST = be)

	    if (be) THEN
        OPEN (UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &  ACTION = "READ", FILE = "gdat/tl_def.gdat", STATUS = "OLD")
C        WRITE (*,*) "opened tl_def.gdat"
        READ (UNIT=8) tl_def

        CLOSE (UNIT=8)
	    ELSE
        OPEN(UNIT=8, FILE= "tbl/tl_def.tbl", STATUS = "OLD",
     &   ACTION = "READ")
C        WRITE (*,*) "opened tl_def.tbl"
         READ (UNIT=8, FMT=*)
         DO mode=1,n_tl
             READ (UNIT=8, FMT=*) 
             DO iii = 1,3
             READ (UNIT=8,FMT=*) 
             DO jj=1,100
                READ (UNIT=8, FMT="(7(x,F10.3))") tl_def(mode,jj,:,iii)
             END DO
             READ (UNIT=8, FMT=*)
             END DO
         END DO 
        close (UNIT=8)

        OPEN ( UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &  ACTION = "WRITE", FILE = "gdat/tl_def.gdat", STATUS = "NEW")
 
        WRITE (UNIT=8) tl_def

        CLOSE (UNIT=8)
        WRITE (*,*) "created tl_def.gdat"
	    END IF 
	 END IF 
!----- transmission coeff with l.s term for neutons and protons

	IF (para%polarization) THEN
	    ALLOCATE (tl_pol(2,100,7,2))

	    INQUIRE (FILE = "gdat/tl_pol.gdat", EXIST = be)

	    if (be) THEN
        OPEN ( UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",  
     &  ACTION = "READ", FILE = "gdat/tl_pol.gdat", STATUS = "OLD")
C        WRITE (*,*) "opened tl_pol.gdat"
        READ (UNIT=8) tl_pol

        CLOSE (UNIT=8)
	    ELSE
        OPEN (UNIT=8, FILE= "tbl/tl_pol.tbl",STATUS = "OLD",
     &   ACTION = "READ")
C        WRITE (*,*) "opened tl_pol.tbl"
         READ (UNIT=8, FMT=*)
         DO mode=1,2
             READ (UNIT=8, FMT=*)
             DO iii = 1,2
             READ (UNIT=8,FMT=*) 
             DO jj=1,100
                READ (UNIT=8, FMT=*) tl_pol(mode,jj,:,iii)
             END DO
             READ (UNIT=8, FMT=*)
             END DO
         END DO 
        close (UNIT=8)

        OPEN ( UNIT=8, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &  ACTION = "WRITE", FILE = "gdat/tl_pol.gdat", STATUS = "NEW")

        WRITE (UNIT=8) tl_pol

        CLOSE (UNIT=8)
        WRITE (*,*) "created tl_pol.gdat"
	    END IF 
	END IF

	if (para%ratio < 0) r_max = -para%ratio
	if (para%k_sum .OR. para%ratio /= 1.0 .or. para%expand) THEN
	    call prepare_interpolation ()
	    ALLOCATE (look_up(50,0:l_max_quantum))
	    ALLOCATE (tl1(0:l_max_quantum))
	    ALLOCATE (tl_table(kk_max,50))

	    if (para%k_sum) THEN
        ALLOCATE (lm_weight(i_def,kk_max))
        ALLOCATE (plm2(i_def,kk_max))
        CALL prepare_plm2 ()
	    ELSE
        ALLOCATE (lm_weight(i_def,1))
	    END IF
	END IF 
	RETURN
	END SUBROUTINE read_tl

	SUBROUTINE clear_look_up ()
	if (.not. modified) RETURN
	look_up = .FALSE.
	tl1 = 10000
	RETURN
	END SUBROUTINE clear_look_up

	SUBROUTINE prepare_plm2 ()
	LOGICAL :: be
	INTEGER :: lll,ii_def,ell,m,kk,k
	REAL (KIND=r8), DIMENSION(0:l_max_quantum) :: coef, coef_old, 
     &	coef_old_old
	REAL (KIND=r8), DIMENSION(:), ALLOCATABLE :: add
	REAL (KIND=r8) :: yy,xx

!see if file containing plm2 values exists
	INQUIRE (FILE = "gdat/plm2.gdat", EXIST = be)

!open file and read plm2 numbers
	 if (be) THEN
	   OPEN ( UNIT=2,   FILE = "gdat/plm2.gdat", STATUS = "OLD",  
     &	  ACTION = "READ", FORM = "UNFORMATTED")   
C	   WRITE (*,*) "opened plm2.gdat"
	   READ (UNIT=2) lll,ii_def
        if (lll == l_max_quantum .AND. ii_def == i_def)THEN
        READ (UNIT=2) plm2
        CLOSE (UNIT=2)
        RETURN
	   END IF  
	   CLOSE (UNIT=2)
	 END IF 

	 ALLOCATE (add(i_def))

	plm2 = 0.0
	kk = 0
	DO ell=0,l_max_quantum
	   i_beta = Lbeta(ell)
	    DO m=0,ell
        kk = kk + 1
        IF (m == 0) THEN
            SELECT CASE (ell)
            CASE (0) 
                coef(0) = 1.0
            CASE (1) 
                coef(0) = 0.0
                coef(1) = 1.0
            CASE DEFAULT
                DO k=0,ell
       coef(k) = -REAL(ell-1,KIND=r8)*coef_old_old(k)/REAL(ell,KIND=r8)
        IF (k > 0) coef(k) = coef(k) +
     &          REAL(2*ell-1,KIND=r8)*coef_old(k-1)/REAL(ell,KIND=r8)
                END DO
            END SELECT
            coef_old_old = coef_old
            coef_old = coef
        ELSE       !derivative of Pl
            DO k=0,ell-m
                coef(k) = coef(k+1)*REAL(k+1,KIND=r8)
            END DO
            coef(ell-m+1) = 0.0
        END IF

    
        loop_beta4: DO i=1,n_points_l(ell)
            xx = COS(beta(i,i_beta))
            yy = 0.0
            DO k=0,ell-m
                yy = yy + coef(k)*xx**REAL(k,KIND=r8)
            END DO
            add(i) = REAL(SIN(beta(i,i_beta)),KIND=r8)
     &            **REAL(2*m,KIND=r8)*yy**2
        END DO loop_beta4
        plm2(:,kk) =  REAL(add(:)/SUM(add),KIND=r4)
	    END DO
	END DO

	DEALLOCATE (add)


!having calculated the plm2"s, put them into a file
	OPEN (UNIT=2, FILE = "gdat/plm2.gdat", STATUS = "NEW", ACTION ="WRITE",
     &     FORM = "UNFORMATTED")
	WRITE (UNIT=2) l_max_quantum,i_def
	WRITE (UNIT=2) plm2
	CLOSE (UNIT=2)
	WRITE (*,*) "created plm2.gdat"
	RETURN
	END SUBROUTINE prepare_plm2

	SUBROUTINE prepare_interpolation ()
	INTEGER :: ell

!quadratic interpolation
	mat2(1,1) =  x1**2
	mat2(1,2) =  x1
	mat2(1,3) =  1.0
	mat2(2,1) = x2**2
	mat2(2,2) = x2
	mat2(2,3) = 1.0
	mat2(3,1) = x3**2
	mat2(3,2) = x3
	mat2(3,3) = 1.0

	CALL inverse(mat2,m_inv_low)

	mat2(1,1) =  x2**2
	mat2(1,2) =  x2
	mat2(1,3) =  1.0
	mat2(2,1) = x3**2
	mat2(2,2) = x3
	mat2(2,3) = 1.0
	mat2(3,1) = x4**2
	mat2(3,2) = x4
	mat2(3,3) = 1.0

	CALL inverse(mat2,m_inv_high)

!cubic interpolation
!mat3(1,1) =  x1**3.0
!mat3(1,2) =  x1**2
!mat3(1,3) =  x1
!mat3(1,4) = 1.0
!mat3(2,1) = x2**3.0
!mat3(2,2) = x2**2
!mat3(2,3) = x2
!mat3(2,4) = 1.0
!mat3(3,1) = x3**3.0
!mat3(3,2) = x3**2
!mat3(3,3) = x3
!mat3(3,4) = 1.0
!mat3(4,1) = x4**3.0
!mat3(4,2) = x4**2
!mat3(4,3) = x4
!mat3(4,4) = 1.0

!CALL inverse(mat3,m3_inv)

	IF (para%k_sum) THEN
	    n_points_b = 0
	    DO ell=0,l_max_quantum
       IF (ell == 0) THEN
           i_beta = 10
       ELSE
           i_beta = INT(50.0/REAL(ell))
       END IF 
       if (i_beta > 10) i_beta = 10
       if (i_beta <= 1) i_beta = 1
       d_beta(ell) = REAL(i_beta)
       Lbeta(ell) = i_beta
       n_points_l(ell) = INT((90.0-d_beta(ell)/2.0)/d_beta(ell)) + 1
       n_points_b(i_beta) =  n_points_l(ell)
	    END DO
	    I_def = n_points_l(l_max_quantum)
	    ALLOCATE (beta(i_def,10))
	    ALLOCATE (ds(i_def,10))
	    ALLOCATE (ratio_r(i_def,3,10))

	     DO i_beta=1,10
        DO i=1,n_points_b(i_beta)
           beta(i,i_beta) = ((REAL(i)-0.5)*REAL(i_beta)-0.01)*pi/180.0
        END DO 
	     END DO

	    kk_max = l_max_quantum*(l_max_quantum+1)/2+l_max_quantum+1
	ELSE
	    IF (para%ratio > 0.0) THEN
        i_def = INT((90.0-d_beta_d/2.0)/d_beta_d) + 1    
        n_points_b(1) = i_def
        ALLOCATE (beta(i_def,1))
        ALLOCATE (ds(i_def,1))
        ALLOCATE (ratio_r(i_def,3,1))

        DO i=1,i_def
           beta(i,1) = ((REAL(i)-0.5)*d_beta_d-0.01)*pi/180.0
        END DO 
	    ELSE  ! gaussian distribution of Coulomb barriers
        i_def = 11
        n_points_b(1) = i_def
        ALLOCATE (ds(i_def,1))
        ALLOCATE (ratio_r(i_def,3,1))
	    END IF 
	    kk_max = l_max_quantum + 1
    
	END IF     

	ratio_r(:,3,:) = 1.0

	RETURN
	END SUBROUTINE prepare_interpolation

	PURE SUBROUTINE get_tll(c,ek,eetl)
   
![subroutine_header_comments]
	REAL (KIND=r4), DIMENSION(3), INTENT(IN) :: c
	REAL (KIND=r4), INTENT(IN) :: ek
	REAL (KIND=r4), INTENT(OUT) :: eetl
	REAL (KIND=r4) :: emax, d2e
!the fits to the transmission coefficients can occosionally
! have problems. check to see if the fit has a maximum
! for positive ek, if so set tl=0 for ek below this value

	if (c(1) == 0.0) THEN
	    emax = 0.0
	    d2e = 0.0
	ELSE
	  emax = c(3)/c(1)
	    IF (emax > 0.0) THEN     !real stationary points
        emax = emax**0.5
        d2e = 2.0*c(3)/emax**3.0
	    ELSE
        emax = 0.0              !no real maximum
        d2e = 0.0
	    END IF
	END IF                          

	if (ek < emax .AND. d2e < 0.0) THEN !fit has a maximum at
                                         !positive ek, set tl=0 below
                                         !this energy
	  eetl = 100.0
	ELSE IF (ek > emax .and. d2e > 0.0) THEN !fit has minimum at
                                          !positive ek, set eetl at
                                          !it minimum value above this
                                          !energy
	  eetl = c(1)*emax + c(2) + c(3)/emax
	ELSE
	  eetl = c(1)*ek + c(2) + c(3)/ek
	END IF 
	
	RETURN
	END SUBROUTINE get_tll

	PURE SUBROUTINE inverse(m,m_inv)
!finds the inverve of a square matrix

	REAL (KIND=r4), INTENT(IN), DIMENSION(:,:)  :: m
	REAL (KIND=r4), INTENT(OUT), DIMENSION(:,:) :: m_inv

	INTEGER :: i,n
	REAL (KIND=r4), DIMENSION(SIZE(m,1)) :: y
	REAL (KIND=r4), DIMENSION(SIZE(m,1),SIZE(m,2)) :: mm

	n = SIZE(m,1)

	DO i=1,n
	    y = 0.0
	    y(i) = 1.0
	    mm = m
	    CALL sle(mm,y)
	    m_inv(:,i) = y
	END DO

	RETURN
	END SUBROUTINE inverse


	PURE SUBROUTINE sle(m,y)
         
!solves a system of linear equations by row reduction method
	REAL (KIND=r4), INTENT(INOUT), DIMENSION(:,:) :: m
	REAL (KIND=r4), INTENT(INOUT), DIMENSION(:)   :: y

	INTEGER :: i,j,n
	REAL (KIND=r4) :: const
	n = SIZE(y)
	
	DO j=1,n
	    const = m(j,j)
	    m(j,j:n) = m(j,j:n)/const
	    y(j) = y(j)/const
	    DO i=1,n
        IF (i /= j) THEN
            const = m(i,j)
            m(i,j:n) = m(i,j:n) - const*m(j,j:n)
            y(i) = y(i) - const*y(j)
        END IF
	    END DO
	END DO
	RETURN
	END SUBROUTINE sle

	SUBROUTINE prepare_tl(mode,iz2,ia2)
!transmission were fit to calculations for residues
!which lie on the evaporation attracter line (a_eal) for each Z value.
!However there is a small A dependence on the Coulomb barriers and
!transmission coefficients
!and this subroutine is a simple scaling of kinetic energy which
!I a find accounts for this reasonably well
	INTEGER, INTENT(IN) :: mode,iz2,ia2
	REAL (KIND=r4) :: a_eal,delta_a
	INTEGER :: l

	mode0 = mode
	scale = 1.0
!correction for A (daughter) dependence of coulomb barriers
	IF (iz2 > 20) THEN
	SELECT CASE (mode)
	  CASE(2)
	     a_eal = 2.045*REAL(iz2) + 3.57e-3*REAL(iz2)**2
	     delta_a = REAL(ia2) - a_eal
	     scale = 1.+delta_a*(-5.412976E-04-9.176213E-02/REAL(iz2)
     &	     + 3.916151E-06*REAL(iz2))
	     scale = 1./scale
	  CASE(7)
	     a_eal = 2.045*REAL(iz2) + 3.57e-3*REAL(iz2)**2
	     delta_a = REAL(ia2) - a_eal
	     scale = 1.+delta_a*(-8.185244E-04-8.234510E-02/REAL(iz2)
     &	     + 1.018610E-05*REAL(iz2))
	     scale = 1./scale
	  END SELECT
	END IF 

	modified = deformed .OR. expanded
	if (mode < 1) modified = .FALSE.
	if (iz2 < 20) modified = .FALSE.

	   DO l=0,l_max_quantum
	    al1 = REAL(l+1)
	    al2 = al1**2
   
	    c2(l,1) = tl_coef(mode,iz2,7)

	    c2(l,2) = tl_coef(mode,iz2,1) + tl_coef(mode,iz2,2)*al1    
     &  +tl_coef(mode,iz2,3)*al2

	    c2(l,3) = tl_coef(mode,iz2,4) + tl_coef(mode,iz2,5)*al1 
     &    +tl_coef(mode,iz2,6)*al2
	   END DO 

	    IF (modified) THEN
       DO l=0,l_max_quantum
        al1 = REAL(l+1)
        al2 = al1**2

        c1(l,1) = tl_def(mode,iz2,7,1)
        c1(l,2) = tl_def(mode,iz2,1,1) + tl_def(mode,iz2,2,1)*al1 
     &   +tl_def(mode,iz2,3,1)*al2

        c1(l,3) = tl_def(mode,iz2,4,1) + tl_def(mode,iz2,5,1)*al1 
     &   +tl_def(mode,iz2,6,1)*al2

        c3(l,1) = tl_def(mode,iz2,7,2)

        c3(l,2) = tl_def(mode,iz2,1,2) + tl_def(mode,iz2,2,2)*al1
     &  + tl_def(mode,iz2,3,2)*al2

        c3(l,3) = tl_def(mode,iz2,4,2) + tl_def(mode,iz2,5,2)*al1 
     &  + tl_def(mode,iz2,6,2)*al2

        c4(l,1) = tl_def(mode,iz2,7,3)

        c4(l,2) = tl_def(mode,iz2,1,3) + tl_def(mode,iz2,2,3)*al1
     &  + tl_def(mode,iz2,3,3)*al2

        c4(l,3) = tl_def(mode,iz2,4,3) + tl_def(mode,iz2,5,3)*al1
     &        + tl_def(mode,iz2,6,3)*al2
	     END DO 
	   END IF 

	   IF ((mode == 1 .OR. mode==2) .AND. para%polarization) THEN
       DO l=0,l_max_quantum
        al1 = REAL(l+1)
        al2 = al1**2

        c1(l,1) = tl_pol(mode,iz2,7,1)

        c1(l,2) = tl_pol(mode,iz2,1,1) + tl_pol(mode,iz2,2,1)*al1
     &        + tl_pol(mode,iz2,3,1)*al2

        c1(l,3) = tl_pol(mode,iz2,4,1) + tl_pol(mode,iz2,5,1)*al1
     &        + tl_pol(mode,iz2,6,1)*al2

        c3(l,1) = tl_pol(mode,iz2,7,2)

        c3(l,2) = tl_pol(mode,iz2,1,2) + tl_pol(mode,iz2,2,2)*al1
     &        + tl_pol(mode,iz2,3,2)*al2

        c3(l,3) = tl_pol(mode,iz2,4,2) + tl_pol(mode,iz2,5,2)*al1
     &        + tl_pol(mode,iz2,6,2)*al2
	   END DO
	END IF
	RETURN
	END SUBROUTINE prepare_tl

!suj->
	SUBROUTINE DE_ARRAY_tl

	IF( ALLOCATED (ee_tl)) DEALLOCATE (ee_tl)
	IF( ALLOCATED (tl_s)) DEALLOCATE (tl_s)
	IF( ALLOCATED (xi)) DEALLOCATE (xi)
	IF( ALLOCATED (factor)) DEALLOCATE (factor)
	IF( ALLOCATED (tl_coef)) DEALLOCATE (tl_coef)
	IF( ALLOCATED (tl_def)) DEALLOCATE (tl_def)
	IF( ALLOCATED (tl_pol)) DEALLOCATE (tl_pol)
	IF( ALLOCATED (look_up)) DEALLOCATE (look_up)
	IF( ALLOCATED (tl1)) DEALLOCATE (tl1)
	IF( ALLOCATED (tl_table)) DEALLOCATE (tl_table)
	IF( ALLOCATED (lm_weight)) DEALLOCATE (lm_weight)
	IF( ALLOCATED (plm2)) DEALLOCATE (plm2)
	IF( ALLOCATED (lm_weight)) DEALLOCATE (lm_weight)
	IF( ALLOCATED (beta)) DEALLOCATE (beta)
	IF( ALLOCATED (ds)) DEALLOCATE (ds)
	IF( ALLOCATED (ratio_r)) DEALLOCATE (ratio_r)
	IF( ALLOCATED (beta)) DEALLOCATE (beta)
	IF( ALLOCATED (ds)) DEALLOCATE (ds)
	IF( ALLOCATED (ratio_r)) DEALLOCATE (ratio_r)
	IF( ALLOCATED (ds)) DEALLOCATE (ds)
	IF( ALLOCATED (ratio_r)) DEALLOCATE (ratio_r)

	RETURN
	END SUBROUTINE DE_ARRAY_tl

	END MODULE tl_stuff
!------------------------------end:   MODULE tl_stuff-----------------------------------

!------------------------------begin:   MODULE stat_mod-----------------------------------
	MODULE stat_mod
	USE gem_kind
	USE mass_stuff
	IMPLICIT NONE
	PRIVATE 

	PUBLIC :: alloc_storage, get_level_density,collective_enhancement

	INTEGER, PARAMETER, PUBLIC :: max_channels = 2000
	INTEGER, PARAMETER, PUBLIC :: n_light_ch = 49
	INTEGER, PARAMETER, PUBLIC :: n_light_ksum = 250000
	INTEGER, PARAMETER, PUBLIC :: n_light_normal = 24000
	INTEGER, PARAMETER, PUBLIC :: Z_shell = 2

	TYPE, PUBLIC :: dau_info
	    INTEGER :: ll,level
	    REAL (KIND=r4) :: e_kinetic,l_p_s,excit
	    REAL (KIND=r4) :: Ex_residue,J_residue
	    REAL (KIND=r4) :: k_residue,s_particle
	    REAL (KIND=r4) :: K_particle
	END TYPE dau_info

	TYPE (dau_info), DIMENSION(n_light_ch), PUBLIC :: dau
	
	INTEGER, DIMENSION(max_channels), PUBLIC :: Za
	INTEGER, DIMENSION(max_channels), PUBLIC :: Aa
	INTEGER, PUBLIC :: I_channels
	INTEGER, DIMENSION(60), PUBLIC :: evap_mode
	INTEGER, DIMENSION(50), PUBLIC :: list_HF,list_HF0
	INTEGER, PUBLIC :: N_light

	REAL (KIND=r4), DIMENSION(max_Channels), PUBLIC :: U_therm,gamma,
     &	  time
	REAL (KIND=r4), PUBLIC :: Mass0,Bs
	REAL (KIND=r4), PUBLIC :: pe,ee,m_inertia_0

	TYPE, PUBLIC :: evap_store
	    INTEGER :: spin,kk2
	    REAL (KIND=r4) :: gam,lps,energy
	END TYPE evap_store
	TYPE (evap_store), DIMENSION(:), ALLOCATABLE, PUBLIC :: storage


	CONTAINS
	    SUBROUTINE alloc_storage ()
        IF (para%k_sum) THEN
            n_light = n_light_ksum
        ELSE
            n_light = n_light_normal
        END IF      
        ALLOCATE (storage(0:n_light))
	    END SUBROUTINE alloc_storage

	    SUBROUTINE get_level_density(U,J,aden,M_Inertia,
     &        entropy,level_density)
      
    ! 
    ! FUNCTIONAL DESCRIPTION: 
    ! 
    ! returns the daughter level density divided by the
    ! parent level density
    ! 
    ! FORMAL PARAMETERS:
    !  
    ! U - thermal excitation energy above yrast line
    ! A - mass number of nucleus
    ! entropy - entropy 
    ! M_inertia - moment of inertia
    ! aden - little "a" level density parameter
    ! J - spin of nucleus
    !  
    !-
	    REAL (KIND=r4), INTENT(IN) :: U,J
	    REAL (KIND=r4), INTENT(IN) :: aden,M_inertia,entropy
	    REAL (KIND=r4), INTENT(OUT) :: level_density

	    REAL (KIND=r4) :: sigma_2,temp
	    REAL (KIND=r4) :: ee_ek,level_density_spherical

        
        ee_ek = entropy - ee
    !   choice of level density


        temp = (U/aden)**0.5 
        sigma_2 = m_inertia*temp/40.848



    !   normal fermi gas level densities for spherical nucleus
        level_density_spherical = (2.0*J+1.0)/pe*EXP(ee_ek)
     &        /(1.0+U**1.25*sigma_2**1.5)/aden**0.25

        level_density = level_density_spherical

	
	    RETURN
	    END SUBROUTINE get_level_density


	    SUBROUTINE collective_enhancement(IZ,IA,U,aden,inertia,
     &	level_density)
	    INTEGER, INTENT(IN) :: iz,ia
	    REAL (KIND=r4), INTENT(IN) :: inertia,U,aden
	    REAL (KIND=r4), INTENT(INOUT) :: level_density
	    REAL (KIND=r4) :: beta2,inertia2,f,e_cr,d_cr,temp,sigma2_perp

	    RETURN

	    CALL getmass(iz,ia,beta_two=beta2)    

!   this follows Hansen and Jenssen. 
	    IF (beta2 == 0.0) RETURN
	    inertia2 = inertia *(1.0+beta2/3.0)
	    E_cr = 120.0*beta2**2*REAL(IA)**0.33333
	    d_cr = 1400.*beta2**2/REAL(ia)**0.6666666

	    f = (U-e_cr)/d_cr
	    IF (f > 65.0) THEN
	    f = 0.0
        RETURN
	    ELSE
        f = 1.0/(1.0+EXP(f))
	    END IF 

	    temp = (U/aden)**0.5 
	    sigma2_perp = inertia2*temp/40.848
	    level_density = ((sigma2_perp-1.0)*f + 1.0)*level_density
    
	    END SUBROUTINE collective_enhancement

	END MODULE stat_mod
!------------------------------end:   MODULE stat_mod-----------------------------------



!------------------------------begin:   MODULE yrast_mod-----------------------------------
	MODULE yrast_mod
 
	USE gem_kind
	USE sierk_stuff
	USE mass_stuff

	IMPLICIT NONE
	PRIVATE

	PUBLIC :: sadfit, read_sad, yrastld, yrast, saddle, def_rot, de_c
	PRIVATE :: cubic

	REAL (KIND=r4), ALLOCATABLE, DIMENSION(:,:,:,:,:), PRIVATE :: C
	REAL (KIND=r4), DIMENSION(150), PRIVATE :: sad_array
	REAL (KIND=r4), PRIVATE :: Bs,Bc
	REAL (KIND=r4), PRIVATE :: br_perp,br_para

	CONTAINS
	SUBROUTINE sadfit(Z0,A0,J0)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine returns the mass-asymmetric saddle
! point energy for the nucleus z0,a0 with angular
! momentum j0. the saddle point energy is obtained from
! a fits to saddle pint energies calculated for two
! spheroids with a yukawa-plus-exponential nuclear force
! of Krappe, Sierk, AND Nix Phys. Rev C20 (1979)992
! using the parameters of Moller and Nix (Nucl. Phy.
! A361 (1981)117). the saddle point energies were
! calculated with the program sadnew of sierk for nuclei
! which lie along the line a=2.08*z+0.0029*z**2.  (see
! charity et al lbl22447) the program calculates the
! fissility parameter and y parameter in the finite range
! model. and uses these in a bi-cubic interpolation to
! obtain the saddle point energies. the barriers are
! then normalized so the symmetric barriers reproduce the
! the sierk barriers calculated by barfit which were
! obtained from a more realistic shape parameterization.
! Bob Charity     5 SEPT 1987
!
! DUMMY ARGUMENTS:
!
!input: Z0 (REAL (KIND=r4) ::)- proton number of decaying nucleus
!       A0 (REAL (KIND=r4) ::)- nucleon number of decaying nucleus
!       J0 (REAL (KIND=r4) ::)- spin of decaying nucleus
!
! IMPLICIT INPUTS:
!
!    none
!
! IMPLICIT OUTPUTS:
!
! sad_array(1:N_array) - finite range model conditional saddle-point energies
!                      less the masses of the two fragments and their Coulomb
!                      energy
!
! SIDE EFFECTS:
!
!    none
!
!
	REAL (KIND=r4), INTENT(IN) :: Z0,A0,J0

	INTEGER :: ii,n_array
	INTEGER :: Iy,Ix,Ia
	INTEGER :: Ia1,Ia2,Iz1,Iz2
	LOGICAL :: first = .TRUE.
	REAL (KIND=r4) :: r1,r2,A1,A2,Z1,Z2
	REAL (KIND=r4) :: E_coul
	REAL (KIND=r4) :: ac,rz,cs,delcs,esz,emz,zsqoa,x,tz,elz,elzohb
	REAL (KIND=r4) :: y,sad0,ess,bfis,sadf
	REAL (KIND=r4) :: alpha,Dkx,Dky,Dka,Rx,Ry,Ra,k1,k2,k3,k4
	REAL (KIND=r4) :: c0,c1,c2,c3
	REAL (KIND=r4) :: mass1,mass2

	REAL (KIND=r4), PARAMETER :: srznw=1.16
	REAL (KIND=r4), PARAMETER :: asnw=21.13
	REAL (KIND=r4), PARAMETER :: aknw=2.3
	REAL (KIND=r4), PARAMETER :: hbarc=197.32858
	REAL (KIND=r4), PARAMETER :: alfinv=137.035982
	REAL (KIND=r4), PARAMETER :: bb=0.99
	REAL (KIND=r4), PARAMETER :: um=931.5016
	REAL (KIND=r4), PARAMETER :: elm=0.51104
	REAL (KIND=r4), PARAMETER :: spdlt=2.99792458e23

	REAL (KIND=r4), DIMENSION(6) :: ky = (/0.0,0.02,0.05,0.10,0.15,
     &      0.20/)
	REAL (KIND=r4), DIMENSION(11):: ka = (/0.0,.1,.2,.3,.4,.5,.6,.7,
     &	  .8,.9,1./)
	REAL (KIND=r4), DIMENSION(8) :: kx = (/0.0,0.08299,0.1632,0.2435, 
     &      0.3217,0.402, 0.51182,0.617422/)
 
!/*----find fissility of nucleus---------------------------*/
	N_array = INT(A0/2.0)
	Ia = 0
	ac = 0.6*hbarc/alfinv/srznw
	rz = srznw*A0**.333333333333333
	cs = asnw * (1.0-aknw*((A0-2.0*Z0)/A0)**2.0)
	delcs = 45.0*hbarc/(8.0*srznw*alfinv)*
     &        ((bb/(1.4142*srznw))**3.0)*(Z0/A0)**2.0
	cs = cs + delcs
	esz = cs*A0**(2.0/3.0)
	emz = um*A0 - elm*Z0
	zsqoa = Z0**2.0/A0
	x = ac*zsqoa/(2.0*cs)

	IF (x > 0.61743) THEN
	   if (first) THEN
C       WRITE (UNIT=*, FMT=*) "No barriers available for this nucleus"	!hg cheng
C       WRITE (UNIT=*, FMT=*) "Z=",Z0,"A=",A0
C       WRITE (UNIT=*, FMT=*) "using scaled 194Hg barriers"
      first = .FALSE.
	   END IF 
	   x = 0.61743
	END IF

!/*----find y fissility parameter-------------------------*/
	IF (esz == 0.0) THEN
	   WRITE (UNIT=*, FMT=*) "esz=0 in sad, Z=",Z0,"A=",A0
	   STOP
	END IF

C 	IF (emz/esz <= 0.) THEN				!Hg Cheng
C 	 WRITE (UNIT=*, FMT=*) "square root of negative for a=",A0,"z=",Z0
C 	 WRITE (UNIT=*, FMT=*) "emz=",emz,"esz=",esz
C 	 STOP
C 	END IF

	tz = rz*sqrt(abs(emz/esz))/spdlt	!Hg Cheng （本来没有abs）
	elz = esz*tz
	elzohb = elz/6.582173e-22
	y = 1.25*(J0/elzohb)**2.0

!/*----------find normalization factor-------------------*/
	IF (Z0 >= 19.0) THEN
	   sad0 = -5.5286e-2 + 190.03*x + 235.189*x**2.0
     &         - 2471.249*x**3.0 + 4266.905*x**4.0 - 2576.048*x**5.0
	   CALL jmax_sierk(Z0,A0)   !initialize sierk sub
	   CALL bfis_sierk(0.0,bfis)
	   sadf = bfis/sad0
	ELSE
	   ess = 2376*x - 6062.4*x**2.0 + 8418.3*x**3.0
	   sadf = (esz/ess)**2.12
	END IF

!/*-----select subroutines for cubic-spline interpolation--*/


!/*find nearest knots for interpolation */
	Iy = 2
	DO
	   IF (y <= ky(Iy)) EXIT
	   Iy = Iy +1
	   IF (Iy > 5) EXIT
	END DO
	Iy = Iy - 1

	ix = 2
	DO
	   IF (x <= kx(Ix)) EXIT
	   Ix = Ix + 1
	   IF (Ix > 6) EXIT
	END DO
	Ix = Ix - 1
! /* now y is between knots ky(iy-1) & ky(iy) and x is  */
! /* between knots kx(ix) & kx(ix+1)       */
	Dkx = kx(Ix+1) - kx(Ix)
	Rx = (x - kx(Ix))/Dkx
	Dky = ky(Iy+1) - ky(Iy)
	Ry = (y - ky(iy))/Dky

! /*--- fill array with saddle point energies-------*/
	DO ii = N_array,1,-1
	   alpha = (A0 - 2.0 * REAL(ii))/A0
	   IF (alpha >= ka(Ia + 1)) THEN
	       DO
           IF (alpha < ka(Ia +1)) EXIT
           Ia = Ia + 1
       END DO
       k1 = cubic( C(Iy,Ix,1,Ia,1), C(Iy,Ix+1,1,Ia,1),
     &        C(Iy,Ix,2,Ia,1), C(Iy,Ix+1,2,Ia,1),Dkx,Rx)
       k2 = cubic( C(Iy,Ix,1,Ia+1,1), C(Iy,Ix+1,1,Ia+1,1),
     &        C(Iy,Ix,2,Ia+1,1), C(Iy,Ix+1,2,Ia+1,1),Dkx,Rx)
       k3 = cubic( C(Iy,Ix,1,Ia,2), C(Iy,Ix+1,1,Ia,2),
     &        C(Iy,Ix,2,Ia,2), C(Iy,Ix+1,2,Ia,2),Dkx,Rx)
       k4 = cubic( C(Iy,Ix,1,Ia+1,2), C(Iy,Ix+1,1,Ia+1,2),
     &        C(Iy,Ix,2,Ia+1,2), C(Iy,Ix+1,2,Ia+1,2),Dkx,Rx)

       Dka = ka(Ia+1) -ka(Ia)

       C0 = k1
       C1 = Dka*k3
       C2 = 3.0*(k2-k1) - (k4+2.0*k3)*Dka
       C3 = 2.0*(k1-k2) + (k3+k4)*Dka

           IF (J0 > 0.1) THEN
           k1 = cubic( C(Iy+1,Ix,1,Ia,1), C(Iy+1,Ix+1,1,Ia,1),
     &           C(Iy+1,Ix,2,Ia,1), C(Iy+1,Ix+1,2,Ia,1),Dkx,Rx)
           k2 = cubic( C(Iy+1,Ix,1,Ia+1,1),                
     &             C(Iy+1,Ix+1,1,Ia+1,1), C(Iy+1,Ix,2,Ia+1,1),
     &             C(Iy+1,Ix+1,2,Ia+1,1),Dkx,Rx)
           k3 = cubic( C(Iy+1,Ix,1,Ia,2), C(Iy+1,Ix+1,1,Ia,2),
     &             C(Iy+1,Ix,2,Ia,2), C(Iy+1,Ix+1,2,Ia,2),Dkx,Rx)  
           k4 = cubic( C(Iy+1,Ix,1,Ia+1,2),       
     &             C(Iy+1,Ix+1,1,Ia+1,2), C(Iy+1,Ix,2,Ia+1,2),
     &             C(Iy+1,Ix+1,2,Ia+1,2),Dkx,Rx)                   
           C0 = (k1-C0)*Ry + C0
           C1 = (Dka*k3-C1)*Ry + C1
           C2 = (3.0*(k2-k1) - (k4+2.*k3)*Dka - C2)*Ry + C2
           C3 = (2.0*(k1-k2) + (k3+k4)*Dka - C3)*Ry + C3
           END IF
       END IF
	   Ra = (alpha - ka(Ia))/Dka
	   sad_array(ii) = (C0 + C1*Ra + C2*Ra**2.0 + C3*Ra**3.0)*sadf
	END DO

	DO ii = para%Z_imf_min,N_array
	   A1 = REAL(ii)
	   A2 = A0 - A1
	   Z1 = A1/A0 * Z0
	   Z2 = Z0 - Z1
	   ia1 = A1
	   ia2 = A2
	   iz1 = Z1
	   iZ2 = Z2
	   CALL getmass(iz1,ia1,Mass_LD=mass1)
	   CALL getmass(Iz1+1,ia1,mass_LD=x)
	   Mass1 = Mass1 + (x-Mass1)*(Z1-REAL(iz1))
	   CALL getmass(iz2,ia2,Mass_LD=mass2)
	   CALL getmass(iz2+1,ia2,mass_LD=x)
	   Mass2 = Mass2 + (x-Mass2)*(Z2-REAL(iz2))

	   r1 = r0*A1**0.33333
	   r2 = r0*A2**0.33333
	   E_coul = Z1*Z2*1.44/(r1 + r2 + sep)
	   sad_array(ii) = sad_array(ii) - Mass1 - Mass2 - E_coul
	END DO
	RETURN
	END SUBROUTINE sadfit


	SUBROUTINE read_sad ()
	INTEGER :: i1,i2,i3,i4,i5
	LOGICAL :: be

	ALLOCATE (c(6,8,2,11,2))
	INQUIRE (FILE = "gdat/sad.gdat", EXIST = be)

	 if (be) THEN
	 OPEN (UNIT=8, FILE = "gdat/sad.gdat", ACTION = "READ", STATUS = "OLD",
     &    ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED")
	    READ (UNIT=8) c
	    CLOSE (UNIT=8)
C	    WRITE (*,*) "opnened sad.gdat"
	ELSE
	 OPEN ( UNIT=8, FILE= "tbl/sad.tbl", STATUS = "OLD", ACTION = "READ")
C	    WRITE (*,*) "opened sad.tbl"
	    READ (UNIT=8, FMT=*)
	    DO i5=1,2
	     DO i4=1,11   
	        DO i3=1,2 
            READ (UNIT=8,FMT=*)
            do i2=1,8
                READ (UNIT=8, FMT=*) (c(i1,i2,i3,i4,i5),i1=1,6)
            END DO
        END DO
	     END DO 
	    END DO

	    CLOSE (UNIT=8)

	    OPEN ( UNIT=8, FILE = "gdat/sad.gdat", ACTION = "WRITE",
     &    STATUS = "NEW", ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED")
	    WRITE (UNIT=8) c
	    CLOSE (UNIT=8)
	    WRITE (*,*) "created sad.gdat"
	END IF
	END SUBROUTINE read_sad
 

	subroutine de_c
	IF( ALLOCATED (c)) DEALLOCATE (c)
	return
	end subroutine de_c


	SUBROUTINE yrastld (Z,A,l,E_rot)
!
! FUNCTIONAL DESCRIPTION:
!
! Finds, by interpolation, the ROTATING LIQUID DROP MODEL yrast line.
! Modified from routine of Plasil o.r.n.l
!
!
! DUMMY ARGUMENTS:
!
! inputs:
! Z REAL (KIND=r4) :: - proton number of nucleus
! A REAL (KIND=r4) :: - nucleon number of nucleus
! l REAL (KIND=r4) :: - spin of nucleus
! outputs:
! E_rot REAL (KIND=r4) :: - rotational energy of nucleus
!
! IMPLICIT INPUTS:
!
!    none
!
! IMPLICIT OUTPUTS:
!
!    none
!
!
! SIDE EFFECTS:
!
!    none
!
!
 
	INTEGER :: ix,iy
	REAL (KIND=r4), INTENT(IN) :: Z,A,l
	REAL (KIND=r4), INTENT(OUT) :: E_rot
	REAL (KIND=r4) :: N
	REAL (KIND=r4) :: paren,eso,x,ero,y,cx,bx,dx,by,cy,dy,h1,h2
	REAL (KIND=r4) :: Hf

	REAL (KIND=r4), DIMENSION(6,11), PARAMETER :: x1h = RESHAPE((/
     &  0.0,0.0,0.0,0.0,0.0,0.0,-0.0057,-0.0058,-0.006,-0.0061,-0.0062,
     & -0.0063,-0.0193,-0.0203,-0.0211,-0.022,-0.023,-0.0245,-0.0402, 
     & -0.0427,-0.0456,-0.0497,-0.054,-0.0616,-0.0755,-0.0812,-0.0899,
     & -0.0988,-0.109,-0.12,-0.1273,-0.1356,-0.147,-0.1592,-0.1745,   
     & -0.1897,-0.1755,-0.1986,-0.2128,-0.2296,-0.251,-0.26,-0.255,   
     & -0.271,-0.291,-0.301,-0.327,-0.335,-0.354,-0.36,-0.365,-0.372, 
     & -0.403,-0.42,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/),
     &  (/6,11/))

	REAL (KIND=r4), DIMENSION(6,11), PARAMETER :: x2h = RESHAPE((/ 
     & 0.0,0.0,0.0,0.0,0.0,0.0,-0.0018,-0.0019,-0.00215,-0.0024,-0.0025,
     & -0.003,-0.0063,-0.00705,-0.0076,-0.0083,-0.0091,-0.0095,-0.015,  
     & -0.0158,-0.0166,-0.0192,-0.0217,-0.025,-0.0245,-0.0254,-0.029,   
     & -0.0351,-0.0478,-0.0613,-0.0387,-0.0438,-0.0532,-0.0622,-0.0845, 
     & -0.0962,-0.0616,-0.0717,-0.0821,-0.0972,-0.1123,-0.1274,-0.0793, 
     & -0.1014,-0.1138,-0.1262,-0.1394,-0.1526,-0.12,-0.134,-0.1503,    
     & -0.1666,-0.1829,-0.1992,-0.1528,-0.171,-0.1907,-0.2104,-0.2301,  
     & -0.2498,0.0,0.0,0.0,0.0,0.0,0.0/),(/6,11/))

	REAL (KIND=r4), DIMENSION(10,20), PARAMETER :: x3h = RESHAPE((/          
     &  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.00012,-0.00014,
     & -0.00016, -0.00018,-0.0002,-0.00024,-0.00029,-0.00036,-0.00065,
     & -0.00089,-0.00047, -0.0005,-.00058,-.00065,-.00074,-.00085,
     & -.00101,-.00124,-.00138,-.00178,-0.001,-0.00105,-0.00124,
     & -0.00138,-0.00156,-0.00179,-0.00275,-0.00292,-0.003,-0.003,
     & -0.00176,-0.0019,-0.00211,-0.00235,-0.00263,-0.00298,-0.00449,
     & -0.0053,-0.0053,-0.0053,-0.003,-0.00308,-0.00318,-0.00352,
     & -0.00392,-0.00417,-0.0062,-0.0062,-0.0062,-0.0062,-0.00374,
     & -0.0041, -0.00444,-0.00488,-0.00521,-0.00545,-0.0066,-0.0066,
     & -0.0066,-0.0066, -0.0053,-0.0055,-0.00585,-0.0064,-0.00695,
     & -0.007,-0.007,-0.007,-0.007, -0.007,-0.00632,-0.007,-0.00742,
     & -0.00792,-0.00856,-0.009,-0.009,-0.009, -0.009,-0.009,-0.0079,
     & -0.0085,-0.01022,-0.0119,-0.012,-0.012,-0.012, -0.012,-0.012,
     & -0.012,-0.00944,-0.0102,-0.0142,-0.0182,-0.019,-0.019,     
     & -0.019,-0.019,-0.019,-0.019,-0.0112,-0.0133,-0.0182,-0.0238,
     & -0.024,-0.024,-0.024,-0.024,-0.024,-0.024,-0.01303,-0.0178,
     & -0.0226,-0.0274, -0.028,-0.028,-0.028,-0.028,-0.028,-0.028,
     & -0.0165,-0.0254,-0.0343, -0.0343,-0.034,-0.034,-0.034,-0.034,
     & -0.034,-0.034,-0.0203,-0.033,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,
     & -0.04,-0.04,-0.025,-0.0406,-0.046,-0.047,-0.047,-0.047,-0.047,
     & -0.047,-0.047,-0.047,-0.03036,-0.0482,-0.048,-0.048,-0.048,
     & -0.048,-0.048,-0.048,-0.048,-0.048,-0.0363,-0.0558,        
     & -0.056,-0.056,-0.056,-0.056,-0.056,-0.056,-0.056,-0.056,-0.04234,        
     & -0.0634,-0.064,-0.064,-0.064,-0.064,-0.064,-0.064,-0.064,
     & -0.064,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/),(/10,20/))

	N = A - Z
	paren = 1.0 - 1.7826*((N-Z)/A)**2 ! 对称能归一化修正因子
	if (paren <= 0.) paren = 0.1  ! （极端中子过剩），强制设为0.1避免数学错误
	eso = 17.9439*paren*A**.66667 ! 液滴模型的表面能项
	x = 0.019655*Z*(Z/A)/paren  ! 形变参数（与电荷分布相关）
	ero = 34.548*l**2/A**1.6667 ! 刚性转子转动能（与 \(I^2/A^{5/3}\) 成正比）
	y = 1.9254*l**2/(paren*A**2.3333) ! 转动-形变耦合参数

	ix = 20.*x + 1.0
	cx = ix
	bx = 20.*x + 1.0
	dx = bx - cx

	IF (x <= 0.25) THEN
	  by = 10.*y + 1.
	  IF (by > 9.0) by = 9.0
	  IF (by < 1.0) by = 1.0
	  iy = by
	  cy = iy
	  dy = by - cy
	  h1 = (x1h(ix+1,iy)-x1h(ix,iy))*dx+x1h(ix,iy)
	  h2 = (x1h(ix+1,iy+1)-x1h(ix,iy+1))*dx+x1h(ix,iy+1)
	  hf = (h2-h1)*dy+h1
	ELSE IF (x < 0.5) THEN
	  by = 20.0*y + 1.0
	  IF (by > 11.0) by = 10.0
	  IF (by < 1.0) by = 1.0
	  ix = ix - 5
	  iy = by
	  cy = iy
	  dy = by - cy
	  h1 = (x2h(ix+1,iy)-x2h(ix,iy))*dx+x2h(ix,iy)
	  h2 = (x2h(ix+1,iy+1)-x2h(ix,iy+1))*dx+x2h(ix,iy+1)
	  hf = (h2-h1)*dy + h1
	ELSE
	  IF (x > 0.9499999) x=0.9499999
	  ix = 20.0*x + 1.0
	  ix = ix - 10
        ix = min(max(ix,0),9) ! 限制ix范围在0~9，防止数组越界
	  by = 100.0*y + 1.0
	  IF (by > 19.0) by = 19.0
	  IF (by < 1.0) by = 1.0
	  iy = by
	  cy = iy
	  dy = by - cy
	  h1 = (x3h(ix+1,iy)-x3h(ix,iy))*dx+x3h(ix,iy)
	  h2 = (x3h(ix+1,iy+1)-x3h(ix,iy+1))*dx+x3h(ix,iy+1)
	  hf = (h2-h1)*dy+h1
	  END IF
	E_rot =  ero + hf*eso
	RETURN
	END SUBROUTINE yrastld

	SUBROUTINE yrast(Z,A,J,K,ratio,E_rot)

!+
! 
! FUNCTIONAL DESCRIPTION: 
! 
! There are two options which are hard wired into the code for when
! you have an angular momentum greater then the value at which
! the fission barrier vanishes. (jmax)
! if kink=.true. the moment-of-inertia at jmax is used to extrapolate
! up from the E_rot value at jmax. As Sierks momentum of inertia
! changes rapidly in the last few hbar, this option produces a kink
! in the "yrast line". Alternatively, if kink=.false. one can
! start the extrapolation at jmax-dj, in a region where the
! moment-of-inertia is more constant and so the "yrast line"
! has no kink and parabolic (nearly) for all angular momenta.
! These yrast lines are only useful for evaporation, as
! equilibrium fission is not defined as the fission barrier
! is zero. If one has a fission delay turned on, then the
! code will ignore fission, otherwise the fission yield
! predicted by gemini is not well defined.
! The kink=.true. option may be relavent if one assumed
! the object evaporating particles is deformed or di-nuclear
! formed in an initially  fast-fission or deep-inelastic
! reaction, the evaporation of light particles during a fission delay
! may remove enough angular momentum for the fission competition to
! be treated statistically after the delay.
! The kink=.FAlse. option may be more applicable for a asymmetric
! entrance channel where a more compact (spherical) object is formed.
! again if enough angular momentum is removed during the delay,
! the fission can be considered statistically.
!  When the deformed delay option is used, the yrast line is the
!  rotation plus deformation energy of a spheroid, ratio is the
!  ratio of the symmetry axis to the other axis.
!-
	REAL (KIND=r4), INTENT(IN) :: Z,A,J,K,ratio
	REAL (KIND=r4), INTENT(OUT) :: E_rot

	REAL (KIND=r4) :: mom_inertia_1
	REAL (KIND=r4) :: E_def,E_yrast,jmax
	REAL (KIND=r4), PARAMETER :: dj = 3
	LOGICAL, SAVE :: kink = .true.

	IF (ratio > 0.0)THEN
!  rotational and deformation energy for rotating spheroid
!  taken from rotating liquid drop model.
!   The Bs and Bc quantities were calculated with formulas from
!   R. Nix"s thesis.
	    IF (ratio /= 1.0) THEN
        E_def = (Bs-1.0)*17.9439* 
     &       (1.0-1.7826*((A-2.0*Z)/A)**2.0)*A**0.66666666 
     &        + (bc-1.0)*0.7053*Z**2.0/A**0.33333333

	    ELSE
        E_def = 0.0
	    END IF
	    E_rot = (Br_perp*(J*(J+1.0)-K**2.0)+Br_para*K**2.0)
     &            *34.540/A**(5.0/3.0)

	    E_yrast = E_rot + e_def
	    E_rot = E_yrast
	    RETURN
	END IF

	IF (z < 19) THEN       !rotating liquid drop yrast line
	    CALL yrastld(Z,A,J,E_rot)
	    RETURN
	ELSE   !rotating finite range model
	    CALL jmax_sierk(Z,A,Jmax) !sierk's yrast line
	    IF (J <= Jmax-dj) THEN
        CALL yrast_sierk(j,E_rot)
	    ELSE IF (kink .AND. J <= jmax) THEN
        CALL yrast_sierk(j,E_rot)
	    ELSE  !there is no well defined yrast line here, but we approx one.
          !by using the moment of inertia at Jmax to extrapolate up.
        CALL m_inertia_sierk(jmax,mom_inertia_1)

        IF (kink) THEN
         CALL yrast_sierk(jmax,E_rot)
         mom_inertia_1 = mom_inertia_1*0.4*r0ld**2.0*A**1.6666666666
         E_rot = E_rot + (J**2.0-jmax**2.0)/(2.0*Mom_inertia_1)*40.848
        ELSE
!           !using the moment of inertia at jmax to extrapolate up,
!           !produces a kink in the yrast line, a less kinked yrast
!           !line can be obtained using the reduced moment of inertia
!           !in the following line - of course it is not clear what the
!           correct yrast line is for such large angular momentum       
            CALL yrast_sierk(jmax-dj,E_rot)
            mom_inertia_1 = (5.0 + mom_inertia_1)/6.0
            mom_inertia_1 = mom_inertia_1*0.4*r0ld**2.0*A**1.6666666666
      E_rot = E_rot + (J**2.0-(jmax-dj)**2.0)/(2.0*Mom_inertia_1)*40.848
        END IF
	    END IF
	END IF
	RETURN
	END SUBROUTINE yrast

	    PURE FUNCTION cubic(a,b,c,d,e,f) RESULT(funct_out)
    !
    ! FUNCTIONAL DESCRIPTION:
    !
    !    Cubic equation for cubic spline interpolation used in subroutine SADFIT
    !
    ! DUMMY ARGUMENTS:
    !
    !    input: a,b,c,d,e (REAL (KIND=4) ::) coefficients
    !           f (REAL (KIND=4) ::) value at which cubic polynomial is evaluated
    !
	     REAL (KIND=r4), INTENT(IN) :: a,b,c,d,e,f
	    REAL (KIND=r4) :: funct_out

	    funct_out = a + f*(e*c+f*(3.0*(b-a)-(d+2.0*c)*e 
     &                +f*(2.0*(a-b)+(c+d)*e)))

	    RETURN
	    END FUNCTION cubic


	    SUBROUTINE saddle(Z1,A1,Z2,A2,E_bar,m_inertia_eff)
!
! FUNCTIONAL DESCRIPTION:
!
! calculates the condition saddle-point deformation + rotational energy for a
! particular binary decay configuration.
!
! DUMMY ARGUMENTS:
!
! INPUT:
!      Z1 (REAL (KIND=r4) ::) proton number of lighter emitted complex fragment
!      A1 (REAL (KIND=r4) ::) nucleon number of lighter emitted complex fragment
!
!      Z2 (REAL (KIND=r4) ::) proton number of heavier emitted complex fragment
!      A2 (REAL (KIND=r4) ::) nucleon number of heavier emitted complex fragment
!
!  output:
!      E_bar (REAL (KIND=r4) ::) - conditional saddle-point deformation plus rotational
!                     energy
!
! IMPLICIT INPUTS:
!
! sad_array - array containing saddle-point energies from finite range model
!             less the mass of the emitted fragments and their coulomb energy
!
!
	    REAL (KIND=r4), INTENT(IN) :: A1,A2,Z1,Z2
	    REAL (KIND=r4), INTENT(OUT) :: E_bar,m_inertia_eff

	    INTEGER :: Ia1,Ia2,Iz1,Iz2,i
	    REAL (KIND=r4) :: Mass1,Mass2
	    REAL (KIND=r4) :: r1,r2,A_red,E_coul
	    REAL (KIND=r4) :: M_inertia_1,M_inertia_2,M_inertia_orbit
	    REAL (KIND=r4) :: M_inertia_saddle,m_inertia_k

!    find  the saddle point energy for the configuration of same mass
!    asymmetry but where the z/a of the fragments is the same
!    as the parent nucleus

	    ia1 = A1
	    ia2 = A2
	    iz1 = Z1
	    Iz2 = Z2
	    i = MIN(ia1,ia2)

!   calculate the mass of the fragments
	    CALL getmass(iz1,ia1,MASS_LD=mass1)
	    CALL getmass(iz2,ia2,MASS_LD=mass2)

    ! calculate coulomb energy
	    R1 = r0*A1**0.3333333
	    R2 = r0*A2**0.3333333
	    M_inertia_1 = 0.4*A1*R1**2
	    M_inertia_2 = 0.4*A2*R2**2
	    A_red = A1*A2/(A1+A2)
	    M_inertia_orbit = A_red*(R1+R2+Sep)**2
	    M_inertia_saddle = M_inertia_1 + M_inertia_2 + M_inertia_orbit
	    E_Coul = 1.44*Z1*Z2/(R1+R2+sep)

C CSM: enable conditional saddle barrier scaling
	    IF (para%b_scale > 0.0 .AND. para%b_scale /= 1.0
     &          .AND. sad_array(i) > 0.0) THEN
	       E_bar = para%b_scale*sad_array(i) + Mass1 + Mass2 + E_coul
	    ELSE
	       E_bar = sad_array(i) + Mass1 + Mass2 + E_coul
	    END IF
 
	    m_inertia_k = m_inertia_1 + m_inertia_2
	    M_Inertia_eff = 1.0/(1.0/M_inertia_k - 1.0/M_inertia_saddle )

	    RETURN
	    END SUBROUTINE saddle

	    SUBROUTINE def_rot (ratio,bss)
	    REAL (KIND=r4), INTENT(IN) :: ratio
	    REAL (KIND=r4), INTENT(OUT) :: bss
	    REAL (KIND=r4) :: ex

	    IF (ratio == 1.0) THEN
        bs = 1.0
        bc = 1.0
        br_perp = 1.0
        br_para = 1.0
	    ELSE
        IF (ratio > 1.0) THEN      !Prolate 
            ex = (1.0-1.0/ratio**2.0)**0.5
            Bs = 0.5*(1.0-ex**2.0)**0.3333333*
     &           (1.0+ASIN(ex)/ex/(1.0-ex**2.0)**0.5)
            Bc = 0.5*(1.0-ex**2.0)**0.3333333*LOG((1.0+ex)/(1.0-ex))/ex
        ELSE !oblate
            ex = (1.0/ratio**2.0-1.0)**0.5
            Bs = 0.5*(1.0+ex**2.0)**0.333333333*
     &         (1.0+LOG(ex + (1.0+ex**2.0)**0.5)/ex/(1.0+ex**2.0)**0.5)
            Bc = (1.0+EX**2.0)**0.33333*ATAN(ex)/ex
        END IF 
        Br_perp = 2.0/(1.0/ratio**0.66666666+ratio**1.33333)
        Br_para = ratio**0.6666666
	    END IF    
	    bss = bs
	    RETURN
	    END SUBROUTINE def_rot 

	END MODULE yrast_mod
!------------------------------end:   MODULE yrast_mod-----------------------------------



!------------------------------begin:   MODULE angle-----------------------------------
	MODULE angle
	USE gem_kind
	USE trig
	USE stat_mod
	USE sierk_stuff
	USE functs
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: sphere_lp, sphere_imf
	PUBLIC :: quantum_lp, quantum_ang_generator,find_theta,quantum_imf
	PUBLIC :: fast_fission_angle
	PUBLIC :: DE_ARRAY_angle    !suj

	REAL (KIND=r4), ALLOCATABLE, DIMENSION(:,:), PRIVATE :: q_ang
	REAL (KIND=r8), ALLOCATABLE, DIMENSION(:), PRIVATE :: sin2,cos2
	REAL (KIND=r4), ALLOCATABLE, DIMENSION(:), PRIVATE :: sin1
	INTEGER, PARAMETER :: n_angle = 180
	CONTAINS

	SUBROUTINE sphere_imf (parent,daughter,U_thermal)
!
! FUNCTIONAL DESCRIPTION:
!
! calculates the spins and excitation energy of the two
! fragments produced in a binary decay. for methods of
! determining spins see schmitt & pacheco nUCL. pHYS.
! a379(1982)313. the thermal excitation energy is divided in
! proportion to the masses of the fragments, but allowing for
! thermal fluctuations.
! When I_angle = .TRUE., the emission angles and velocities are
! also calculated
!
!
! DUMMY ARGUMENTS:
!
! input: I_chan - exit-channel number of binary decay
!        k_chan - subscript number of decaying system.
!
! IMPLICIT INPUTS:
!
! U_therm(I_chan) - thermal excitation energy of exit-channel
! J(k_chan) - spin of decaying system
! A(N_frag),A(N_frag-1) - nucleon number of binary fragments.
! Z(N_frag),Z(N_frag-1) - proton number of the two binary fragments
! theta_j(k_chan),phi_j(K_chan) - polar angles of the spin vector of the
!                                 decaying system
! v(k_chan),theta_v(k_chan),phi_v(k_chan) - magnitude and polar angles of the
!                                 velocity vector of the decaying system
!
! IMPLICIT OUTPUTS:
!
! J(N_frag),J(N_frag-1) - spins of the two binary fragments
! EX(N_frag),EX(N_frag-1) - excitation energies of these fragments
! theta_j(N_frag),phi_j(N_frag) - polar angles of the spin vector of the
! theta_j(N_frag-1),phi_j(N_frag-1)   binary fragments
! v(N_frag),theta_v(N_frag),phi_v(N_frag) magnitude and polar angles
! v(N_frag-1),theta_v(N_frag-1),phi_v(N_frag-1) of the velocity vectors of
!                                                  the fragments
!
! SIDE EFFECTS:
!
!    none
!
!

	REAL, INTENT(IN) :: u_thermal
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter

	REAL (KIND=r4) :: A1,A2,Z1,Z2,r1,r2,A0,Ared,Aden,aden1,aden2
	REAL (KIND=r4) :: theta,c,s,ss,cc,J0
	REAL (KIND=r4) :: rua,fact10,fact1,fact3,fact4
	REAL (KIND=r4) :: a_bending,a_wriggling
	REAL (KIND=r4) :: a_tilting,a_twisting
	REAL (KIND=r4) :: de,sig,demax,demin,yy,y
	REAL (KIND=r4) :: theta_J0_cm,phi_J0_cm
	REAL (KIND=r4) :: theta_j0_r,phi_j0_r
	REAL (KIND=r4) :: theta_symmetry_j0,phi_symmetry_j0
	REAL (KIND=r4) :: dphi,Ek,ratio_l,vrel,v1_cm,v2_cm
	REAL (KIND=r4) :: theta_symmetry_cm,phi_symmetry_cm
	REAL (KIND=r4) :: theta_s1_r,phi_s1_r
	REAL (KIND=r4) :: theta_s1_j0,phi_s1_j0
	REAL (KIND=r4) :: theta_s1_cm,phi_s1_cm
	REAL (KIND=r4) :: theta_s2_j0,phi_s2_j0
	REAL (KIND=r4) :: theta_s2_cm,phi_s2_cm
	REAL (KIND=r4) :: theta_s2_r,phi_s2_r
	REAL (KIND=r4) :: R1x,R1y,R1z,R2x,R2y,R2z
	REAL (KIND=r4) :: omega_t,ratio_1,ratio_2,m_inertia_saddle
	REAL (KIND=r4) :: S1x,S1y,S1z,S2x,S2y,S2z,S1,S2,m_inertia_orbit
	REAL (KIND=r4) :: M_inertia_12,M_inertia_1,M_inertia_2
	REAL (KIND=r4) :: K  ! spin projection along disintegration axis
	REAL (KIND=r4) :: J_twisting ! angular momentum in twisting mode
	REAL (KIND=r4) :: J_bending_x !angular momentum in one bending mode
	REAL (KIND=r4) :: J_bending_z !angular momentum in other bending mode
	REAL (KIND=r4) :: J_wriggling_x !angular momentum in one wriggling mode
	REAL (KIND=r4) :: J_wriggling_z !angular momentum in other wriggling mode
	REAL (KIND=r4) :: temp ! temperature of saddle point configuration.
	REAL (KIND=r4) :: U  !thermal energy of saddle point configuration
	REAL (KIND=r4) :: U1 ! thermal excitation energy of fragment 1
	REAL (KIND=r4) :: U2 ! thermal excitation energy of fragment 2
	REAL (KIND=r4) :: Ex1 ! total excitation energy of fragment 1
	REAL (KIND=r4) :: Ex2 ! total excitation energy of fragment 2
	REAL (KIND=r4) :: J0xz ! projection of J0 perpendicular to symmetry axis
	REAL (KIND=r4) :: J0x  ! x component of J0 in rotating coordinate system
	REAL (KIND=r4) :: J0z  ! z component of J0
	REAL (KIND=r4) :: l,lz,lx ! components of orbital angular momentum
	REAL (KIND=r4) :: surf,entropy1,entropy2,ran

	u = u_thermal
	J0 = parent%J
	A1 = daughter(1)%A
	A2 = daughter(2)%A
	Z1 = daughter(1)%Z
	Z2 = daughter(2)%Z
	r1 = r0*A1**0.333333
	r2 = r0*A2**0.333333
	M_inertia_1 = 0.4*A1*r1**2.0
	M_inertia_2 = 0.4*A2*r2**2.0
	A0 = A1 + A2
	Ared = A1*A2/A0
	M_inertia_orbit = Ared*(r1+r2+sep)**2.0

	U1 = U*A1/(A1+A2)
	U2 = U - U1
	surf = 1.0

	CALL get_little_a(A1,Z1,0.,U1,surf,0.0,0.0,aden1,entropy1,0.0,1.0)
	CALL get_little_a(A2,Z2,0.,U2,surf,0.0,0.0,aden2,entropy2,0.0,1.0)

	aden = aden1 + aden2
	temp = SQRT(U/aden)
	rua  = entropy1 + entropy2

	M_inertia_12 = M_inertia_1 + M_inertia_2
	M_inertia_saddle = M_inertia_12 + M_inertia_orbit

	 fact10 = M_inertia_orbit*(M_inertia_2-M_inertia_1)
     &     /(2.0*M_inertia_1*M_inertia_2)
	fact1 = (fact10**2.0+1.0)**0.5
	theta = atand(fact10+fact1)
	fact3 = M_inertia_12/(M_inertia_1*M_inertia_2)
	c = COSD(theta)
	s = SIND(theta)
	ss = s*s
	cc = c*c
	fact4 = 2.0*s*c/M_inertia_orbit
	 a_bending = (cc/M_inertia_1+ss/M_inertia_2 
     &           +1.0/M_inertia_orbit-fact4)*41.566
	 a_wriggling = (ss/M_inertia_1+cc/M_inertia_2 
     &           +1.0/M_inertia_orbit+fact4)*41.566 
	 a_tilting = M_inertia_orbit/M_inertia_12/M_inertia_saddle*41.566
	a_twisting = fact3*41.566

	de = u + 1.0

	DO 
!   choose k_tilting
	    CALL ran_j(K,a_tilting,U,aden,rua,J0)

!   choose j_twisting
	    sig = 5.0*(temp/a_twisting)**0.5
	    CALL ran_j(J_twisting,a_twisting,U,aden,rua,sig)

!   choose j_bending_x
	    sig = 5.0*(temp/a_bending)**0.5
	    CALL ran_j(J_bending_x,a_bending,U,aden,rua,sig)

!   choose j_bending_z
	    CALL ran_j(J_bending_z,a_bending,U,aden,rua,sig)

!   choose j_wriggling_x
	    sig = 5.0*(temp/a_wriggling)**0.5
	    CALL ran_j(J_wriggling_x,a_wriggling,U,aden,rua,sig)

!   choose j_wriggling_z
	    CALL ran_j(J_wriggling_z,a_wriggling,U,aden,rua,sig)

	    de = 0.5*(a_tilting*K**2.0+a_twisting*J_twisting**2.0 
     &    + a_bending*(J_bending_x**2.0+J_bending_z**2.0)       
     &    + a_wriggling*(J_wriggling_x**2.0+J_wriggling_z**2.0))

	    IF (de <= U) EXIT

	END DO

!calculate components of r1, depolarizing mode of fragment 1
	R1x = J_bending_x*c + J_wriggling_x*s
	R1y = J_twisting
	R1z = J_bending_z*c + J_wriggling_z*s

!       calculate components of r2
	R2x = J_wriggling_x*c - J_bending_x*s
	R2y = -J_twisting
	R2z = J_wriggling_z*c - J_bending_z*s

	CALL RANDOM_NUMBER(ran)
	omega_t = ran*360.0
	J0xz = (J0**2.0-K**2.0)**0.5
	J0x = J0xz*sind(omega_t)
	J0z = J0xz*cosd(omega_t)
	ratio_1 = M_inertia_1/M_inertia_saddle

!calculate components of s1
	S1x = J0x*ratio_1 + R1x
	S1y = M_inertia_1/M_inertia_12*K + R1y
	S1z = J0z*ratio_1 + R1z
	s1 = (s1x**2.0 + S1y**2.0 + s1z**2.0)**0.5

!calculate components of s2
	ratio_2 = M_inertia_2/M_inertia_saddle
	S2x = J0x*ratio_2 + R2x
	S2y = M_inertia_2/M_inertia_12*K + R2y
	S2z = J0z*ratio_2 + R2z
	S2 = (S2x**2.0 + S2y**2.0 + s2z**2.0)**0.5

	U = U - de
	temp = (U/aden)**0.5
	U = max(0.0,U-temp)
	temp = (U/aden)**0.5
	sig = (2.0*temp**3.0*aden1*aden2/aden)**0.5
	U1 = U*A1/A0
	U2 = U - U1
	rua = 2.0*(aden*U)**0.5
	demax = MIN(5.0*sig,U1)
	demin = MIN(5.0*sig,U2)

	DO
	    CALL RANDOM_NUMBER(ran)
	    de = (demin+demax)*RAN - demax
	  y = exp(2.0*(aden1*(U1+de))**0.5 + 2.0*(aden2*(U2-de))**0.5-rua)
	    CALL RANDOM_NUMBER(yy)
	    IF (yy <= y) EXIT
	END DO

	Ex1 = U1 + de + S1**2.0*20.8993/M_inertia_1
	Ex2 = U2 - de + s2**2.0*20.8993/M_inertia_2

	daughter(1)%Ex = Ex1
	daughter(1)%J = nint(s1)
	IF (para%ssum .AND. MODULO(daughter(1)%A,2) /= 0) THEN
	    CALL RANDOM_NUMBER(ran)
	    IF (RAN > 0.5) THEN
	        daughter(1)%j = daughter(1)%j + 0.5
	    ELSE
        daughter(1)%j = ABS(daughter(1)%j - 0.5)
	    END IF
	END IF

	daughter(2)%Ex = Ex2
	daughter(2)%J = nint(s2)
	IF (para%ssum .AND. MODULO(daughter(2)%a,2) /= 0) THEN
	    CALL RANDOM_NUMBER(ran)
	    IF (RAN > 0.5) THEN
        daughter(2)%j = daughter(2)%j + 0.5
	    ELSE
        daughter(2)%j = ABS(daughter(2)%j - 0.5)
	    END IF
	END IF

	daughter(1)%k = 0.0
	daughter(2)%k = 0.0

	IF (.NOT. para%I_angle) RETURN
! this part of the subroutines calculates the final
! centre-of-mass velocity vectors of the complex fragments and
! their spin vectors.  vectors are specified as:
!                x_y_z
! x is polar angular coordinate for vector (theta or phi)
! y  identifies the vector
! options  :    V1 -  velocity of fragment 1
!             S1 -  spin of fragment 1
!             V2 -  velocity of fragment 2
!             S2 -  spin of fragment 2
!            symmetry - specifies the symmetry axis of
!            saddle-scission point configuration
!             J0 - spin vector of fragment 0 (decaying system)
!             V0 - velocity vector of fragment 0
! Z is frame in which polar coordinates are specified
!             r - rotation frame of saddle point configuration
!             see Schmitt and Pacheco Nucl Phys. A379 (1982) 313
!             J0 - frame with Vector J0 as Z axis
!             cm - centre of mass frame of reaction.


! spin vector of compound system
	theta_J0_cm = parent%spin%theta
	phi_j0_cm = parent%spin%phi


! in the frame of the di-nuclear complex the total angular
! momentum vector processes. The theta and phi angles at
! scission are. see schmitt and pacheco for coordinate system.
	IF (j0 == 0.0) THEN
	   CALL RANDOM_NUMBER(ran)
	   phi_j0_r = 360.0*RAN
	   CALL RANDOM_NUMBER(ran)
	   theta_J0_r = ACOSD(1.0-2.0*RAN)
	ELSE
	   theta_j0_r = ACOSD(J0z/J0)
	   phi_J0_r = ATAN2D(K,J0x)
	END IF

! the symmetric axis occurs at y=0,z=0 in the above coordinate
! system i.e. theta_symmetry_r=90., phi_symmetry_r =0., find
! angles of the symmetry axis relative to the j0 vector
! (.i.e. z axis parallel to J0)
	   CALL transform(90.0,0.0,-theta_j0_r,phi_j0_r 
     &	 ,theta_symmetry_J0,phi_symmetry_J0)

! randomize the phi coordinate.
	CALL RANDOM_NUMBER(ran)
	dphi = 360.0*RAN
	phi_symmetry_J0 = MODULO(phi_symmetry_J0+dphi,360.0)

! find symmetry axis in lab frame
	 CALL transform(theta_symmetry_J0,phi_symmetry_J0, 
     &	theta_J0_cm,phi_J0_cm, theta_symmetry_cm,phi_symmetry_cm)

! in the rotating frame of dinuclear system, the spin vector of
! fragment 1 is s1X,s1Y,s1Z, relative to di nuclear coordinate system
! mentioned above
	IF (s1 == 0) THEN
	   CALL RANDOM_NUMBER(ran)
	   theta_s1_r = ACOSD(1.0-2.0*RAN)
	   CALL RANDOM_NUMBER(ran)
	   phi_s1_r = 360.0*RAN
	ELSE
	   theta_s1_r = ACOSD(S1z/S1)
	   phi_s1_r = ATAN2D(S1y,S1x)
	END IF

! rotate to frame of J0 vector (.i.e. z axis parallel to J0)
	 CALL transform(theta_s1_r,phi_s1_r,  
     &	-theta_J0_r,phi_j0_r,theta_s1_J0,phi_s1_J0)
	 phi_s1_J0 = MODULO(phi_s1_j0+dphi,360.0)
! rotate to lab frame      
	 CALL transform(theta_s1_J0,phi_s1_j0, 
     &	theta_J0_cm,phi_J0_cm,theta_s1_cm,phi_s1_cm)

!  do the same for fragment 2
	 IF (s2 == 0.0) THEN
	   CALL RANDOM_NUMBER(ran)
	   theta_s2_r = ACOSD(1.0-2.0*RAN)
	   CALL RANDOM_NUMBER(ran)
	   phi_s2_r = 360.0*RAN
	 ELSE
	   theta_s2_r = ACOSD(S2z/S2)
	   phi_s2_r = ATAN2D(S2y,S2x)
	 END IF

! rotate to frame of J0 vector (.i.e. z axis parallel to J0)
	CALL TRANSFORM(theta_s2_r,phi_s2_r, 
     &	-theta_j0_r,phi_J0_r,theta_s2_J0,phi_s2_J0)
	phi_s2_J0 = MODULO(phi_s2_j0+dphi,360.0)
! rotate to lab frame      
	CALL transform(theta_s2_J0,phi_s2_j0, 
     &    	theta_J0_cm,phi_J0_cm,theta_s2_cm,phi_s2_cm)

! calculate velocities in frame of emitting system
! first find the kinetic energy of separation
! Ek = 1.44*Z1*Z2/(r1+r2+sep)
! generlization of Viola see Hinde Nucl Phys A472, 318 (1987).
	 Ek = 0.755*Z1*Z2/(A1**0.333333+A2**0.333333333) + 7.3

! find orbital angular momentum
	ratio_l = M_inertia_orbit/M_inertia_saddle
	lx = J0x*ratio_l - R1x - R2x
	lz = J0z*ratio_l - R1z - R2z
	l = (lx**2.0 + lz**2.0)**0.5
!ek = ek + 20.8993*l**2.0/M_inertia_orbit

	Vrel = (2.0*Ek/Ared)**0.5*0.9794
	V1_cm = A2/A0*Vrel
	V2_cm = Vrel-V1_cm

! in the c-o-m frame the velocity vectors are
	daughter(1)%vel = vell(V1_cm,theta_symmetry_cm,phi_symmetry_cm) + 
     &	parent%vel 

! fragment 2 get emitted 180 degrees away.
	theta_symmetry_cm = 180.0 - theta_symmetry_cm
	phi_symmetry_cm = MODULO(phi_symmetry_cm+180.0,360.0)

	daughter(2)%vel = vell(V2_cm,theta_symmetry_cm,phi_symmetry_cm) + 
     &	 parent%vel

	daughter(1)%spin%theta = theta_s1_cm
	daughter(2)%spin%theta = theta_s2_cm
	daughter(1)%spin%phi = phi_s1_cm
	daughter(2)%spin%phi = phi_s2_cm

	RETURN
	END SUBROUTINE sphere_imf


	SUBROUTINE sphere_lp(parent,daughter,ll,l_p_s,e_kinetic)
!
! FUNCTIONAL DESCRIPTION:
!
! calculates the spins and excitation energy of the residual
! system after a light particle evaporation & the emission
! angles of the fragments if i_angle = .true.
!
! DUMMY ARGUMENTS:
!
! input: I_chan - exit-channel number of binary decay
!        k_chan - subscript number of decaying system.
!
! IMPLICIT INPUTS:
!
! J(k_chan) - spin of decaying system
! j(N_frag) - spin of residue
! ll(I_chan) - orbital angular momentum
! e_kinetic(I_chan) - energy of emitted particle
! A(N_frag),A(N_frag-1) - nucleon number of binary fragments.
! Z(N_frag),Z(N_frag-1) - proton number of the two binary fragments
! theta_j(k_chan),phi_j(K_chan) - polar angles of the spin vector of the
!                                 decaying system
! v(k_chan),theta_v(k_chan),phi_v(k_chan) - magnitude and polar angles of the
!                                 velocity vector of the decaying system
!
! IMPLICIT OUTPUTS:
!
! J(N_frag),J(N_frag-1) - spins of the two binary fragments
! EX(N_frag),EX(N_frag-1) - excitation energies of these fragments
! theta_j(N_frag),phi_j(N_frag) - polar angles of the spin vector of the
! theta_j(N_frag-1),phi_j(N_frag-1)   binary fragments
!
! v(N_frag),theta_v(N_frag),phi_v(N_frag) magnitude and polar angles
! v(N_frag-1),theta_v(N_frag-1),phi_v(N_frag-1) of the velocity vectors of
!                                                  the fragments
!
!
! SIDE EFFECTS:
!
!    none
!
!

	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
	INTEGER, INTENT(IN) :: ll
	REAL (KIND=r4), INTENT(IN) :: l_P_s,e_kinetic

	REAL (KIND=r4) :: Z1,Z2,A1,A2,J0,l
	REAL (KIND=r4) :: theta,phi,v1,v2,vrel,s2,l_Plus_s,s1
	REAL (KIND=r4) :: theta_e,theta_ee,theta_eee,phi_e,phi_ee,phi_eee
	REAL (KIND=r4) :: theta_l,phi_l,theta_s1,phi_s1,theta2,phi2,
     &				  theta3,phi3
	REAL (KIND=r4) :: ran

	A1 = daughter(1)%A
	A2 = daughter(2)%A
	Z1 = daughter(1)%Z
	z2 = daughter(2)%Z
	J0 = parent%J
	s2 = daughter(2)%J
	l = REAL(ll)
	IF (para%tl_iwbc) THEN
	   l_plus_s = l_P_s
	   s1 = daughter(1)%J
	ELSE
	   l_plus_s = l
	END IF

! find orientation of residue spin vector with respect to s0
! vector (i.e. z axis parallel to s0)
	IF (s2 > 0.0 .AND. j0 > 0.0) THEN
	    theta = ACOSD((j0**2.0 + s2**2.0 - l_plus_s**2.0)/(2.0*s2*j0))
	ELSE
	    CALL RANDOM_NUMBER(ran)
	    theta = ACOSD(1.0-2.0*RAN)
	    END IF
	CALL RANDOM_NUMBER(ran)
	phi = 360.0*RAN
! rotate to lab frame      
	 CALL transform (theta,phi,parent%spin%theta,  
     &	parent%spin%phi, daughter(2)%spin%theta,daughter(2)%spin%phi)

! find direction of emitted fragment
! first the orientation of the vector l_plus_s with respect to s0
! (i.e. z axis parallel to s0)
	IF (l_plus_s > 0.0 .AND. J0 > 0.0) THEN
	 theta=ACOSD((j0**2.0 + l_plus_s**2.0-s2**2.0)/(2.0*l_plus_s*j0))
	ELSE
	    CALL RANDOM_NUMBER(ran)
	    theta = ACOSD(1.0-2.0*RAN)
	  END IF
	phi = MODULO(180.0+phi,360.0)

	IF (para%tl_iwbc .AND. s1 > 0.0) THEN
!   find orientation of l vector with respect to l_plus_s vector
!   (i.e. z axis parallel to l_plus_s)
	    IF (l_plus_s > 0.0 .AND. l > 0.0) THEN
        theta_l = ACOSD((l_plus_s**2.0 + l**2.0 - s1**2.0) 
     &       /(2.0*l_plus_s*l))        
	    ELSE
        CALL RANDOM_NUMBER(ran)
        theta_l = ACOSD(1.0-2.0*RAN)
        END IF
	    CALL RANDOM_NUMBER(ran)
	    phi_l = 360.0*RAN
!           rotate to frame of s0 (i.e. z axis parallel to s0)
	    CALL transform (theta_l,phi_l,theta,phi,theta2,phi2)

!           orientation of s1 vector  with respect to l_plus_s
!           (i.e. z axis parallel to l_plus_s)
	    IF (l_plus_s > 0.0 .AND. s1 > 0.0) THEN
        theta_s1 = ACOSD((l_plus_s**2.0 + s1**2.0 - l**2.0) 
     &       /(2.0*l_plus_s*s1))
	    ELSE
        CALL RANDOM_NUMBER(ran)
        theta_s1 = ACOSD(1.0-2.0*RAN)
        END IF
	    phi_s1 = MODULO(180.0+phi_l,360.0)
!   rotate to frame of s0 vector (i.e. z axis parallel to s0)
	   CALL transform (theta_s1,phi_s1,theta,phi,theta3,phi3)
!   rotate to lab frame      
	   CALL transform (theta3,phi3, parent%spin%theta,parent%spin%phi,
     &   daughter(1)%spin%theta,daughter(1)%spin%phi)
	    theta = theta2
	    phi = phi2
	END IF

	if (l <= l_max_quantum) THEN
!   a better approximation can be obtained if theta
!   is chosen from the distribution for m=l
	   CALL find_theta(NINT(l),NINT(l),theta_e)
	ELSE
!  the particle is emitted perpendicular to the l vector
	   theta_e = 90.0
	END IF 

	 CALL RANDOM_NUMBER(ran)
	 phi_e = 360.0*RAN
! rotate into the frame of the s0 vector (i.e. z axis parallel to s0)
	CALL transform (theta_e,phi_e,theta,phi,theta_ee,phi_ee)
! rotate into laboratory frame
	CALL transform (theta_ee,phi_ee, parent%spin%theta,parent%spin%phi   
     &                ,theta_eee,phi_eee)

	IF (e_kinetic < 0.0) THEN
	  WRITE (UNIT=*, FMT=*) "A1,A2,Z1,Z2,j0,s2,l",A1,A2,Z1,Z2,j0,s2,l
!    CALL crash ()
	END IF

	vrel = (2.0*e_kinetic*(A1+A2)/A1/A2)**0.5*0.9794
	v1 = A2/(A2+A1)*vrel

	daughter(1)%vel = vell(v1,theta_eee,phi_eee) + parent%vel

	v2 = vrel - v1
	theta_eee = 180.0-theta_eee
	phi_eee = MODULO(180.0+phi_eee,360.0)

	daughter(2)%vel = vell(v2,theta_eee,phi_eee) + parent%vel

	RETURN
	END SUBROUTINE sphere_lp
 

	SUBROUTINE quantum_lp(parent,daughter,ll,l_p_s,e_kinetic)

![subroutine_header_comments]
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
	INTEGER, INTENT(IN) :: ll
	REAL (KIND=r4) :: l_p_s,e_kinetic
	REAL (KIND=r4) :: s0,m0,s2,m2,l,l_plus_s,s1,m_lps
	REAL (KIND=r4) :: m_lps_min,m_lps_max,m2_max,m2_min,m_s1_min,
     &	              m_s1_max
	REAL (KIND=r4) :: prob,prob1,ml_max,ml_min,ml,m_s1,ran
	REAL (KIND=r4) :: clebsch_gordan,theta,phi,A1,A2,vrel,v1,v2

	s0 = parent%J
	m0 = parent%M
	s2 = daughter(2)%J
	l = REAL(ll)
	 IF (para%tl_iwbc) THEN
	   l_plus_s = l_P_s
	   s1 = daughter(1)%J
	 ELSE 
	   l_plus_s = l
	 END IF

! find m2 and m(l_plus_s)
! first find the limits for these
	m_lps_min = - l_plus_s
	m_lps_max = l_plus_s
	m2_max = m0 - m_lps_min
	IF (m2_max > s2) m2_max = s2
	m2_min = m0 - m_lps_max
	IF (m2_min < -s2) m2_min = -s2

	if (m2_max == m2_min) THEN
	    m2 = m2_max
	    m_lps = m0 - m2
	ELSE
	    CALL RANDOM_NUMBER(prob1)
	    prob = 0.0
	    m2 = m2_max + 1
	    DO 
        m2 = m2 - 1
        m_lps = m0 - m2

        clebsch_gordan=Cl_go(NINT(s2*2),NINT(l_plus_s*2.0),NINT(S0*2.0),
     &        NINT(m2*2.0),NINT(m_lps*2.0),NINT(m0*2.0))

        if (clebsch_gordan > 1.0) THEN
        WRITE (UNIT=*, FMT=*) "clebsch > 1", clebsch_gordan
        WRITE (UNIT=*, FMT=*) "s0,m0=",s0,m0,NINT(S0*2.0),NINT(m0*2.0)
        WRITE (UNIT=*, FMT=*) "s2,m2=",s2,m2,NINT(s2*2.0),NINT(m2*2.0)
        WRITE (UNIT=*, FMT=*) "lps,m_lps",l_plus_s,m_lps,
     &	    NINT(l_plus_s*2.0),NINT(m_lps*2.0)
        END IF 

        prob = prob + clebsch_gordan**2.0
        IF (prob1 <= prob .OR. m2 <= m2_min) EXIT
	    END DO
	END IF 
	daughter(2)%M = m2

	 IF (daughter(1)%z < 0) THEN      !gamma decay
	    daughter(2)%vel = parent%vel
	    daughter(1)%vel = parent%vel
	    daughter(1)%m = 0.0
	    RETURN
	 END IF

	m_lps = m0 - m2
! first find the limits for these
	m_s1_min = - s1
	m_s1_max = s1
	ml_max = m_lps - m_s1_min 
	IF (ml_max > l) ml_max = l
	ml_min = m_lps - m_s1_max
	IF (ml_min < -l) ml_min = -l

	CALL RANDOM_NUMBER(prob1)
	prob = 0.0
	ml = ml_max + 1.0
	DO 
	    ml = ml - 1.0
	    m_s1 = m_lps - ml
	 clebsch_gordan = Cl_go(NINT(l*2),NINT(s1*2.0),NINT(l_plus_s*2.0), 
     &    NINT(ml*2.0),NINT(m_s1*2.0),NINT(m_lps*2.0))
	    prob = prob + clebsch_gordan**2.0
	    IF (prob1 <= prob .OR. ml <= ml_min) EXIT
	END DO
	daughter(1)%m = m_s1

	CALL find_theta(ll,NINT(ml),theta)

	CALL RANDOM_NUMBER(ran)
	phi = 360.0*RAN

	A2 = daughter(2)%A
	A1 = daughter(1)%A

	if (e_kinetic < 0.0) THEN
	    WRITE (UNIT=*, FMT=*) "ekinetic < 0 "
!    call crash ()
	END IF 

	vrel = (2.0*e_kinetic*(A1+A2)/A1/A2)
	if (vrel < 0) THEN
	    WRITE (UNIT=*, FMT=*) "ekinetic=", e_kinetic
	    WRITE (UNIT=*, FMT=*) "A1=",A1
	    WRITE (UNIT=*, FMT=*) "A2=",A2
!    call crash ()
	END IF 
	vrel = vrel**0.5*0.9794
	v1 = A2/(A2+A1)*vrel
	daughter(1)%vel= vell(v1,theta,phi) + parent%vel

	v2 = vrel - v1
	theta = 180.0-theta
	phi = MODULO(180.0+phi,360.0)
	daughter(2)%vel = vell(v2,theta,phi) + parent%vel

	RETURN 
	END SUBROUTINE quantum_lp
   

	SUBROUTINE quantum_ang_generator ()
![subroutine_header_comments]
	INTEGER :: l,m,kk,k,i
	LOGICAL :: be

	REAL (KIND=r8) :: coef(0:l_max_quantum)
	REAL (KIND=r8) :: coef_old(0:l_max_quantum)
	REAL (KIND=r8) :: coef_old_old(0:l_max_quantum)
	REAL (KIND=r8) :: theta,x,y,sumation
	REAL (KIND=r8) :: add(90)

!REAL (KIND=r4) :: q_ang(l_max_quantum*(l_max_quantum+1)/2+l_max_quantum+1,90)

! no need to allocate this array unless we are going to calculate
! angle
! fission arrays
	IF (para%quantum) THEN
	    ALLOCATE (sin2(n_angle))
	    ALLOCATE (cos2(n_angle))
	    ALLOCATE (sin1(n_angle))

	    DO i = 1,n_angle
        theta = (REAL(i) - 0.5)*180.0/REAL(n_angle)
        sin2(i) = SIND(REAL(theta/2.,KIND=r8))
        cos2(i) = COSD(REAL(theta/2.,KIND=r8))
        sin1(i) = SIND(theta)
	    END DO
	END IF

!light particle
	ALLOCATE (q_ang(l_max_quantum*(l_max_quantum+1)/2+l_max_quantum+1,
     &	      90))

!see if file containing plm2 values exists
	INQUIRE (FILE = "gdat/q_ang.gdat",  EXIST = be)

!open file and read plm2 numbers
	 if (be) THEN
	   OPEN (UNIT=2, FILE = "gdat/q_ang.gdat", STATUS = "OLD", 
     &	ACTION = "READ", FORM = "UNFORMATTED")   
C	   WRITE (*,*) "opened q_ang.gdat"
	   READ (UNIT=2) l
        if (l == l_max_quantum)THEN
        READ (UNIT=2) q_ang
        CLOSE (UNIT=2)
        RETURN
	  END IF  
	   CLOSE (UNIT=2)
	END IF 

	kk = 0
	DO l=0,l_max_quantum
	    DO m=0,l
        kk = kk + 1
        sumation = 0.0
        IF (m == 0) THEN
            SELECT CASE (l)
            CASE (0) 
                coef(0) = 1.0
            CASE (1)
                coef(0) = 0.0
                coef(1) = 1.0
            CASE DEFAULT
                DO k=0,l
                    coef(k) = -REAL(l-1,KIND=r8)*coef_old_old(k)
     &       /REAL(l,KIND=r8)
                    IF (k > 0) coef(k) = coef(k) +              
     &              REAL(2*l-1,KIND=r8)*coef_old(k-1)/REAL(l,KIND=r8)
                END DO
            END SELECT
            coef_old_old = coef_old
            coef_old = coef
        ELSE                               !derivative of Pl
            DO k=0,l-m
                coef(k) = coef(k+1)*REAL(k+1,KIND=r8)
            END DO
            coef(l-m+1) = 0.0
        END IF
        DO i=1,90
            theta = 90.5 - REAL(i,KIND=r8)
            x = COSD(theta)
            y = 0.0
            DO k=0,l-m
                y = y + coef(k)*x**REAL(k,KIND=r8)
            END DO
            add(i) = SIND(theta)**REAL(2*m+1,KIND=r8)*y**2.0
            sumation = sumation + add(i) 
        END DO
        DO i=1,90
            q_ang(kk,i) = REAL(add(i)/sumation,KIND=r4)
            IF (i > 1) q_ang(kk,i) = q_ang(kk,i) + q_ang(kk,i-1)
        END DO
	    END DO
	END DO

!having calculated the q_ang's, put them into a file
	OPEN (UNIT=2, FILE = "q_ang.gdat", STATUS="NEW", ACTION="WRITE",
     &     FORM = "UNFORMATTED")   
	WRITE (UNIT=2) l_max_quantum
	WRITE (UNIT=2) q_ang
	CLOSE (UNIT=2)
	WRITE (*,*) "created q_ang.gdat"

	RETURN
	END  SUBROUTINE quantum_ang_generator

	SUBROUTINE find_theta(l,m,angle)
	INTEGER, INTENT(IN) ::l,m
	REAL (KIND=r4), INTENT(OUT) :: angle
	INTEGER :: i, j, mm
	REAL (KIND=r4) :: prob, signs
	IF (l > l_max_quantum) THEN
	    WRITE (UNIT=*, FMT=*) "l=",l_max_quantum," too large for qm 
     &	 angular dist"
!    CALL crash ()
	    STOP
	END IF

	mm = ABS(m)
	i = l*(l+1)/2 + l + 1 - (l-mm)
	CALL RANDOM_NUMBER(prob)
	j = 1
	DO
	   IF (prob <= q_ang(i,j)) EXIT
	   j = j + 1
	   IF (j >= 90) EXIT 
	END DO
	signs = 1.0
	CALL RANDOM_NUMBER(prob)
	IF (prob > 0.5) signs = -1.0
	CALL RANDOM_NUMBER(prob)
	angle = 90.0 + signs*(REAL(J)-0.5) + prob - 0.5
	END SUBROUTINE find_theta


	 SUBROUTINE quantum_imf (parent,daughter,U_saddle,e_kinetic,
     &                        U_scission)

	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
	REAL (KIND=r4) , INTENT(IN) :: U_saddle,e_kinetic
	REAL (KIND=r4) , INTENT(IN) :: U_scission

	INTEGER :: i,jj,ix,if1,if2,if3,if4,if5
	REAL (KIND=r4) :: j0,m,k,theta,phi,prob,fract,entropy,temp_saddle
	REAL (KIND=r4) :: k_start,Z0,a0,U,E_rot_k,surf,aden,ef,sigma_2
	REAL (KIND=r4) :: mom_inertia_1,mom_inertia_3,m_inertia_k,
     &	              m_inertia_eff
	REAL (KIND=r4) :: gam_k,x,ared,vrel,v1_cm,v2_cm
	REAL (KIND=r4) :: m_inertia_saddle,delta_k
	REAL (KIND=r4) :: z1,a1,jmax,U1,U2,aden1,aden2,ran,demax,demin,
     &	              de,y,yy
	REAL (KIND=r8) :: signs,total,temp,sig,rua,m_inertia
	REAL (KIND=r8), DIMENSION(n_angle) :: dist
	j0 = parent%J
	m = ABS(parent%M)   !note as d(J,M,K)**2 = d(J,-M,-K)**2. and we must consider
                  !both +-K, we can deal with only positive M values
                  !d(j,K,M0,theta) = (-1)**(K-M)*d(J,-K,-M,theta)
                  !d(J,K,m,theta) = (-1)**(k-M)*d(J,M,K,theta)
                  !d(J,K,M,theta) = d(J,M,K,-theta)
	A0 = parent%A
	Z0 = parent%Z

!must choose a K value to calculate angular distribution
	IF (j0 == 0.0 ) THEN
	     K = 0.0
	ELSE IF (J0 == 0.5) THEN
	     K = 0.5
	ELSE
	 CALL m_inertia_sierk(j0,mom_inertia_1,m_inertia_3=mom_inertia_3)
	 m_inertia_saddle = mom_inertia_1*0.4*r0ld**2.0*A0**1.6666666666
       m_inertia_k = mom_inertia_3*0.4*r0ld**2.0*A0**1.6666666666
       M_Inertia_eff = 1.0/(1.0/M_inertia_k - 1.0/M_inertia_saddle)
       CALL bs_sierk(j0,surf)

       IF(MODULO(NINT(A0),2) /= 0 .AND. para%ssum) THEN
         k = 0.5
       ELSE
         k = 0.0
       END IF
       k_start = k

       jj = 0

       DO
         E_rot_k = 20.424*K**2.0/M_inertia_eff
         U = U_saddle - E_rot_k
         IF (U <= 0.0) EXIT
         CALL  get_little_a(A0,Z0,0.,U,surf,0.0,0.0,aden,entropy,
     &         0.0,para%a_scale)
         ef = entropy
         IF (ef-ee > 69.0) EXIT
         temp_saddle = (U/aden)**0.5
         sigma_2  = m_inertia_0*temp_saddle/40.848
         gam_k = 1./(1.0+U**1.25*sigma_2**1.5)
     &         /aden**0.25*(2.0*j0+1.0)/pe*EXP(ef-ee)
!        K can be positive or negative
         if (k > 0.0)gam_k = 2.0*gam_k
         jj = jj + 1
         dist(jj) =  gam_k
         IF (jj > 1) dist(jj) = dist(jj-1)+ dist(jj)
         IF (gam_k <= dist(1)/100.0) EXIT
         k = k + 1.0
         if (k > J0) EXIT
       END DO

       if (jj == 0) THEN
         WRITE (UNIT=*,FMT=*) "jj = 0 in quantum imf"
         WRITE (UNIT=*,FMT=*) "j0=",j0,"m=",m,"A0=",A0,"Z0=",Z0
         WRITE (UNIT=*,FMT=*) "U_saddle",U_saddle
       END IF 

       i = 1
       CALL RANDOM_NUMBER(prob)
      
       DO
         fract = dist(i)/dist(jj)
         IF (fract > prob) EXIT
         i = i + 1
       END DO
       k = REAL(i-1) + K_start
	END IF

!post-scission reorientation effect
! B. Back et al , Phys Rev C 32 , 195 (1985).
! k gets decreased by this amount
	IF (ABS(k) > 0.5) THEN
	    delta_K = k - 
     &    k/(1.0+122.0*(J0*(J0+1.0)-k**2.0)/(a0**1.333333*z0**2.0))

	    k = k - AINT(delta_k)
	    delta_k = delta_k - AINT(delta_k)

	    IF (delta_k < 1.0) THEN
       CALL RANDOM_NUMBER(prob)
       if (prob < delta_k) k = k - 1.0
	    ELSE
       WRITE (UNIT=*, FMT=*) "delta_k =",delta_k
       stop
	    END IF 
	END IF 

!one can set m to 0 or 0.5 to see difference in angular distributions
!m = j0 - AINT(j0)


!construct angular distribution
	dist = 0.0
	ix = MAX(0,NINT(m-k))
	ix = ix - 1
	DO
       ix = ix + 1
       x = REAL(ix)
       if1 = NINT(j0-k-x)
       IF (if1 < 0) EXIT

       if2 = NINT(j0+m-x)
       IF (if2 < 0) EXIT

       if3 = NINT(x+k-m)
       if4 = NINT(k-m+2.0*x)
       if5 = NINT(2.0*j0-k+m-2.0*x)

        if (MODULO(ix,2) == 1) THEN
        signs = -1.0_r8
       ELSE
        signs = 1.0_r8
       END IF 

       dist = dist + signs*sin2**REAL(if4,KIND=r8)
     &    *cos2**REAL(if5,KIND=r8)/(fct(if1)*fct(if2)*fct(if3)*fct(ix))
	 END DO

	dist = dist**2.0*sin1
	total = SUM(dist)
	dist = dist/total

	DO i = 2,n_angle
	    dist(i) = dist(i) + dist(i-1)
	END DO

	CALL RANDOM_NUMBER(prob)

	i = 1
	DO
	   IF (dist(i) >= prob) EXIT
	   i = i + 1
	   IF (i > n_angle) EXIT
	END DO

	CALL RANDOM_NUMBER(prob)
	theta = (REAL(i) - 1.0 + prob)*REAL(n_angle)/180.0

!we have choosen only positive K.
!However there is a 50% chance that K is negative
	CALL RANDOM_NUMBER(prob)
	IF (prob > 0.5) theta = 180.0 - theta

	CALL RANDOM_NUMBER(prob)
	phi = 360.0*prob 


!The following is an approximation of the angular distribution
!from Wheeler, 
!DO i = 1,180
!    theta = REAL(i) - 0.5
!    dist(i) = SIND(theta)**2.0/(((j0+0.5)*SIND(theta))**2.0&
!    - m**2.0 - k**2.0 + 2.0*m*k*COSD(theta))**0.5
!    total = total + dist(i)
!END DO
!
!dist(1) = dist(1)/total
!
!DO i = 2,180
!    dist(i) = dist(i)/total + dist(i-1)
!END DO
!
!DO i = 2,180
!    dist(i) = dist(i) + dist(i-1)
!END DO

	  ared = REAL(daughter(2)%A*daughter(1)%A)
     &         /REAL(daughter(2)%A+daughter(1)%A)
	 Vrel = (2.0*E_kinetic/Ared)**0.5*0.9794
	 V1_cm = REAL(daughter(1)%A)/A0*Vrel
	 V2_cm = Vrel-V1_cm

! in the c-o-m frame the velocity vectors are
	 daughter(1)%vel = vell(V1_cm,theta,phi) + parent%vel

! fragment 2 get emitted 180 degrees away.
	 theta = 180.0 - theta
	 phi = MODULO(phi+180.0,360.0)
	 daughter(2)%vel = vell(V2_cm,theta,phi) + parent%vel

	daughter(1)%Ex = 0.0
	daughter(2)%ex = 0.0

	IF (para%imf_option /= 1) RETURN
!for symmetric fission - we allow post fission evaporation
! devide the excitation energy at scission betweem the fragements
! the fragments have equal temperatures on average, but I allow
! thermal fluctuations

!fragment 1
	surf = 1.0
	z1 = daughter(1)%z
	A1 = daughter(1)%A
	IF (para%ssum .AND. MODULO(daughter(1)%a,2) /= 0)THEN
	   daughter(1)%j = ANINT(2./9.*j0+0.5) - 0.5
	   daughter(1)%m = 0.5
	ELSE
	   daughter(1)%j = ANINT(2./9.*j0)
	   daughter(1)%m = 0.
	END IF
	CALL jmax_sierk(z1,a1,jmax)

	IF (jmax == -1.) THEN
	    M_inertia = 0.4*A1*(1.2*A1**0.333333)**2
	    daughter(1)%ex = 20.424*daughter(1)%j**2/m_inertia
	ELSE
	    CALL yrast_sierk(daughter(1)%j,daughter(1)%ex)
	END IF 

	U1 = U_scission*A1/A0
	CALL get_little_a(A1,Z1,0.,U1,surf,0.0,0.0,aden1,entropy,0.0,1.0)

!fragment2
	z1 = daughter(2)%z
	A1 = daughter(2)%A
	IF (para%ssum .AND. MODULO(daughter(2)%a,2) /= 0)THEN
	   daughter(2)%j = ANINT(2./9.*j0+0.5) - 0.5
	   daughter(2)%m = 0.5
	ELSE
	   daughter(2)%j = ANINT(2./9.*j0)
	   daughter(2)%m = 0.
	END IF
	CALL jmax_sierk(z1,a1,jmax)

	IF (jmax == -1.) THEN
	    M_inertia = 0.4*A1*(1.2*A1**0.333333)**2
	    daughter(2)%ex = 20.424*daughter(2)%j**2/m_inertia
	ELSE
	    CALL yrast_sierk(daughter(2)%j,daughter(2)%ex)
	END IF 

	U2 = U_scission-U1
	CALL get_little_a(A1,Z1,0.,U2,surf,0.0,0.0,aden2,entropy,0.0,1.0)

!now for fluctuations
	aden = aden1 + aden2
	temp = (U_scission/aden)**0.5
	sig = (2.0*temp**3.0*aden1*aden2/aden)**0.5
	rua = 2.0*(aden*U_scission)**0.5
	demax = MIN(REAL(5.0*sig,KIND=r4),U1)
	demin = MIN(REAL(5.0*sig,KIND=r4),U2)
	
	DO
	    CALL RANDOM_NUMBER(ran)
	    de = (demin+demax)*RAN - demax
	    y = exp(2.0*(aden1*(U1+de))**0.5 
     &	     +2.0*(aden2*(U2-de))**0.5-rua)
	    CALL RANDOM_NUMBER(yy)
	    IF (yy <= y) EXIT
	END DO

	daughter(1)%ex = daughter(1)%ex + U1 + de
	daughter(2)%ex = daughter(2)%ex + U2 - de
	daughter(1)%k = 0.0
	daughter(2)%k = 0.0

	RETURN
	END SUBROUTINE quantum_imf

	SUBROUTINE fast_fission_angle(parent,daughter,ek,U_scission)
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter
	REAL (KIND=r4), INTENT(IN) :: ek,u_scission

	REAL (KIND=r4) :: vv,v1,v2,theta,phi,ran,rua,temp,sig,demax,demin
	REAL (KIND=r4) :: A1,A2,z1,z2,j0,a0,jmax,entropy,y,yy
	REAL (KIND=r4) :: de,aden1,aden2,surf,aden,u1,u2

	z1 = daughter(1)%z
	A1 = daughter(1)%A
	z2 = daughter(2)%z
	A2 = daughter(2)%A
	j0 = parent%j
	A0 = parent%a

	CALL RANDOM_NUMBER(ran)
	theta = 180.0*ran
	CALL RANDOM_NUMBER(ran)
	phi = 360.0*ran
	vv = (2.0*ek*A0/(A1*A2))**0.5*0.9784
	v1 = A2/A0*vv
	v2 = vv - v1

	daughter(1)%vel = vell(V1,theta,phi) + parent%vel
	theta = 180.0 - theta
	phi = MODULO(phi+180.0,360.0)
	daughter(2)%vel = vell(V2,theta,phi) + parent%vel

	surf = 1.0
	IF (para%ssum .AND. MODULO(daughter(1)%a,2) /= 0)THEN
	   daughter(1)%j = ANINT(2./9.*j0+0.5) - 0.5
	   daughter(1)%m = 0.5
	ELSE
	   daughter(1)%j = ANINT(2./9.*j0)
	   daughter(1)%m = 0.
	END IF
	CALL jmax_sierk(z1,a1,jmax)
	CALL yrast_sierk(daughter(1)%j,daughter(1)%ex)
	U1 = U_scission*A1/A0
	CALL get_little_a(A1,Z1,0.,U1,surf,0.0,0.0,aden1,entropy,0.0,1.0)

!fragment2
	IF (para%ssum .AND. MODULO(daughter(2)%a,2) /= 0)THEN
	   daughter(2)%j = ANINT(2./9.*j0+0.5) - 0.5
	   daughter(2)%m = 0.5
	ELSE
	   daughter(2)%j = ANINT(2./9.*j0)
	   daughter(2)%m = 0.
	END IF
	CALL jmax_sierk(z2,a2,jmax)
	CALL yrast_sierk(daughter(2)%j,daughter(2)%ex)
	U2 = U_scission-U1
	CALL get_little_a(A2,Z2,0.,U2,surf,0.0,0.0,aden2,entropy,0.0,1.0)

!now for fluctuations
	aden = aden1 + aden2
	temp = (U_scission/aden)**0.5
	sig = (2.0*temp**3.0*aden1*aden2/aden)**0.5
		IF (aden < 0. .OR. u_scission <0 ) 
     &	WRITE (*,*)'aden=',aden,'U=',U_scission
	rua = 2.0*(aden*U_scission)**0.5
	demax = MIN(5.0*sig,U1)
	demin = MIN(5.0*sig,U2)

	DO
	    CALL RANDOM_NUMBER(ran)
	  de = (demin+demax)*RAN - demax
	  y = exp(2.0*(aden1*(U1+de))**0.5 +2.0*(aden2*(U2-de))**0.5-rua)
	    CALL RANDOM_NUMBER(yy)
	    IF (yy <= y) EXIT
	END DO

	daughter(1)%ex = daughter(1)%ex + U1 + de
	daughter(2)%ex = daughter(2)%ex + U2 - de
	daughter(1)%k = 0.0
	daughter(2)%k = 0.0
	END SUBROUTINE fast_fission_angle

!suj
	SUBROUTINE DE_ARRAY_angle

	IF( ALLOCATED (storage) ) DEALLOCATE (storage)
	IF( ALLOCATED (sin2 ) ) DEALLOCATE (sin2 )
	IF( ALLOCATED (cos2 ) ) DEALLOCATE (cos2 )
	IF( ALLOCATED (sin1 ) ) DEALLOCATE (sin1 )
	IF( ALLOCATED (q_ang) ) DEALLOCATE (q_ang)

	RETURN
	END SUBROUTINE DE_ARRAY_angle
!suj

	END MODULE angle
!------------------------------end:   MODULE angle-----------------------------------



!------------------------------begin:   MODULE evap-----------------------------------
	MODULE evap
!
! FUNCTIONAL DESCRIPTION:
!
! Calculates the decay-widths in MeV for light particle emission Z<3 from the
! Hauser-Feshbach formalism. The summation over both spin and orbital angular
! momentum is performed.
!
! DUMMY ARGUMENTS:
!
! input: mode (INTEGER ::)
!        = 1 FOR NEUTRON EMISSION
!        = 2 FOR PROTON
!        = 3     DEUTERON
!        = 4     TRITON
!        = 5     HELIUM-3
!        = 6     ALPHA PARTICLE
!        = 7     Lithium-6
!        = 8     Lithium-7
!
!      Z0 (REAL (KIND=4) ::) proton number of decaying nucleus
!      A0 (REAL (KIND=4) ::) nucleon number of decaying nucleus
!      Ex0 (REAL (KIND=4) ::) excitation energy (MeV) of decaying nucleus
!      J0 (REAL (KIND=4) ::) spin of decaying nucleus
!
! IMPLICIT INPUTS:
!
!    none
!
! IMPLICIT OUTPUTS:
!
! I_Channels - number of open light particle channels
! gamma(1:I_channels) - decay widths for each light particle channel
! Za(1:_I_channels) - proton number of light fragment
! Aa(1:I_channels)  - nucleon number of light fragment
! Ex_residue(1:I_channels) - excitation energy of daughter nucleus

	USE gem_kind
	USE mass_stuff
	USE functs
	USE yrast_mod
	USE stat_mod
	IMPLICIT NONE
	PRIVATE 

	PUBLIC :: evaporation_simple,prepare_evap

	INTEGER, PRIVATE :: mode_gs,Iz2,Ia2,S0,s2,s2_best,jj
	INTEGER, PRIVATE :: i
	INTEGER, PRIVATE :: l,lms,l_min,s2_min

	REAL (KIND=r4), PRIVATE :: Z2,A2,Z23,A23,const,Ec,pair_energy
	REAL (KIND=r4), PRIVATE :: separation_c,r2,level_density
	REAL (KIND=r4), PRIVATE :: mass2,mass2ld,shell_corr,ba,entropy
	REAL (KIND=r4), PRIVATE :: E_rot,E_sad,J2,U2,ef0,aden
	REAL (KIND=r4), PRIVATE :: ran,temp,ef,m_inertia_orbit,m_inertia_2
	REAL (KIND=r4), PRIVATE :: fract,prob,prob2,ek,ranges
	REAL (KIND=r4), PRIVATE :: m_inertia_saddle

	INTEGER, PARAMETER, PRIVATE :: N_mode=15
	REAL (KIND=r4), DIMENSION(n_mode), PARAMETER, PRIVATE :: Z1 = 
     & (/0.0,1.0,1.0,1.0,2.0,2.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0/)

	REAL (KIND=r4), DIMENSION(N_mode), PARAMETER, PRIVATE :: A1 = 
     & (/1.0,1.0,2.0,3.0,3.0,4.0,6.0,6.0,6.0,6.0,6.0,6.0,7.0,7.0,7.0/)

	REAL (KIND=r4), DIMENSION(N_mode), PARAMETER, PRIVATE :: S = 
     & (/0.5,0.5,1.0,0.5,0.5,0.0,1.0,3.0,0.0,2.0,2.0,1.0,1.5,0.5,3.5/)

	REAL (KIND=r4), DIMENSION(N_mode), PARAMETER, PRIVATE :: DEF = 
     &	(/8.07,7.28,13.13,14.949,14.93,2.424, 
     &	14.0854,16.368,17.648,18.40,19.46,19.73,     !6Li
     &    14.9067,15.384,19.537/)                      !7Li

	INTEGER, DIMENSION(N_mode), PARAMETER, PRIVATE :: level_no = 
     &    (/0,0,0,0,0,0,0,1,2,3,4,5,0,1,2/)

! separation for Coulomb barrier
	REAL (KIND=r4), DIMENSION(5,2:8), PARAMETER, PRIVATE :: 
     &	sep_c = RESHAPE(
     &	(/ 4.509, 1.466, -0.5413, -3.218,   0.0,     !proton
     &    6.723, 1.205, -0.5938, -0.2433,  0.0,        !deuteron
     &    6.319, 1.346, -0.5901, -0.8494,  0.0,        !triton
     &    5.812, 1.542, -0.8719,  0.7053,  0.0,        !3He
     &    5.935, 1.308, -0.5129, -0.4928,  0.0,        !alpha particle
     &    6.425, 1.543, -0.8504,  0.3892,  0.0,		 !6Li
     &    6.245, 1.536, -0.8019,  1.5555,  0.0/),		 !7Li
     &	(/5,7/))

!angular momentum dependence of Coulomb barrier, i.e. (l(l+1)) dependence
	REAL(KIND=r4), DIMENSION(6,8), PARAMETER,PRIVATE :: vc0l=RESHAPE(
     & (/-0.0111, -0.02063 , -0.8885, 1.508, 0.03129 , -0.02062 ,  !neutron
     &  -0.0111, -0.02063 , -0.8885, 1.508, 0.03129 , -0.02062 ,   !proton
     &  -0.0001, -469.3   , -2.555,  5.238, 0.1851  , -0.6560 ,   !deuteron
     &  -0.0011, -24.61   , -3.773,  6.793, 0.2540  , -0.5061 ,   !triton
     &  -0.0001, -376.1   , -3.886,  6.426, 0.1821  , -0.4884 ,   !3He
     &  -0.0061, -0.6657  , -3.841,  6.913, 0.2497  , -0.2416 ,   !alpha
     &   0.0021, -99.80   , -3.291,  5.537,-0.3208  ,  1.616721,  !6Li
     &   0.0021, -89.63   , -5.523,  8.092,-0.0324  ,  1.35290/), !7Li
     &	(/6,8/))

	CONTAINS
	SUBROUTINE evaporation_simple(mode,Z0,A0,Ex0,J0)

	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN) :: Z0,A0,EX0,J0

	IF (mode > n_mode) RETURN

	z2 = z0 - z1(mode)
	a2 = a0 - A1(mode)
	IF (z1(mode) > z2) RETURN
	IF (A1(mode) > A2) RETURN
	IF (z2 >= A2) RETURN
	A23 = A2**0.3333333333
	Z23 = Z2**0.3333333333
!-------------------------------------------------------------------------

	SELECT CASE(mode)
	CASE (-2:6)
	   mode_gs = mode
	CASE (7:12)
	   mode_gs = 7
	CASE (13:15)
	   mode_gs = 8
	END SELECT

	CALL get_pairing(NINT(z2),NINT(A2),pair_energy)
	if (z0 < Z_shell) pair_energy = 0 

	IF (mode == 1) THEN
	   Ec = 0.0
	ELSE
! find radius for calculation of Coulomb barrier
	   separation_c = sep_c(1,mode_gs) + sep_c(2,mode_gs)*A23
     &    + sep_c(3,mode_gs)*Z23 + sep_c(4,mode_gs)/Z2 
	   Ec = Z1(mode)*Z2*1.44/separation_c
	END IF

!find radius for calculation of l dependence term of the coulomb barrier
	  const = 1.0/(vc0l(2,mode_gs)*(1.0-EXP(-vc0l(1,mode_gs)*A2)) 
     &   + vc0l(3,mode_gs)*Z23 + vc0l(4,mode_gs)*A23 
     &   + vc0l(5,mode_gs)*z23**2.0 + vc0l(6,mode_gs)*A23**2.0)

	M_inertia_Orbit = 20.424/const
	R2 = r0*A2**0.3333333
	M_inertia_2 = 0.4*A2*R2**2
	M_inertia_saddle = M_inertia_2 + M_inertia_orbit

! calculate the mass defect of the residual system
	iz2 = z2
	ia2 = a2
	CALL getmass(iZ2,iA2,Mass2,mass2ld)
	IF (para%mass_option /= 1 .AND. z0 > Z_shell) THEN
	    shell_corr = mass2 - mass2ld - pair_energy
	ELSE
	    shell_corr = 0.0
	END IF

! calculate the binding energy
	ba = Mass2 + def(mode) - Mass0

	IF (A0 < 8) THEN
	  E_rot = J0**2.0*20.8993/M_inertia_saddle
	  E_sad = Ba + E_rot + Ec
	  J2 = J0
	ELSE
	  J2 = J0*M_inertia_2/M_inertia_saddle
	  CALL yrast(Z2,A2,J2,0.0,0.0,E_rot)
	  E_rot = E_rot + (J0-J2)**2.0*20.8993/M_inertia_orbit
	  E_sad = Ba + E_rot + Ec
	END IF
	s2_best = INT(J2) + 1

	U2 = Ex0 - E_sad
	IF (U2 < 0) RETURN
	jj = 0
	s0 = INT(J0)


	 CALL get_little_a(A2,Z2,j2,U2,bs,shell_corr,E_rot,aden,entropy,
     &      pair_energy,1.0)

	ef0 = entropy

	IF (U2 > 40) THEN
	    CALL light_particle_quick (mode,Ex0)
	ELSE
	    CALL light_particle (mode,Ex0)
	END IF

	IF (jj == 0) RETURN
	IF (storage(jj)%gam < 0) RETURN
	I_channels = I_channels + 1
	Za(I_channels) = NINT(Z1(mode))
	Aa(I_channels) = NINT(A1(mode))
	gamma(I_channels) = storage(jj)%gam

	dau(I_channels)%Ex_residue= Ex0- Ba- Ec- REAL(l*(l+1))*const - ek 
	dau(I_channels)%J_residue = s2
	dau(I_channels)%s_particle = s(mode)
	dau(I_channels)%level = level_no(Mode)
	IF (para%I_angle) THEN
	   dau(I_channels)%ll = l
	   dau(I_channels)%e_kinetic = ek + ec + REAL(l*(l+1))*const
	   dau(I_channels)%l_p_s = l
	   dau(I_channels)%s_particle = 0
	END IF
	RETURN

	END SUBROUTINE evaporation_simple

	SUBROUTINE light_particle (mode,Ex0)

	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN) :: EX0
	REAL (KIND=r4) :: eff

! loop over lms = l - s2
! l = orbital angular momentum of the emitted particle
! s2 = spin of residual
! lms is confined to limits -J0 to J0
	lms = -s0
	eff = ef0
	s2_min = s0 + 1

	DO
	   IF (lms > s0) EXIT
	   IF ((ef0-eff) > 4.0 .AND. s2_min < s2_best) EXIT
!  find minimum values of J2 and l which are consistent
!  with this value of lms
	   l_min = INT(REAL(s0+lms+1)/2.0)
	   s2_min = l_min - lms
	   s2 = s2_min
	   ef = ef0
! loop over spin of residual nucleus
	   DO
       IF (ef0-ef > 4) EXIT
       l = s2 + lms                              !orbital angular momentum
       CALL yrast(Z2,A2,REAL(s2),0.0,0.0,E_rot)  !rotational energy of residue
       U2 = Ex0 - Ba - Ec - E_rot - REAL(l*(l+1))*const
       if (u2 > 0)CALL get_little_a(A2,Z2,REAL(s2),U2,bs,shell_corr,
     &   E_rot, aden,entropy,pair_energy,1.0)

       IF (U2 > 0.0) THEN
           ef = entropy
           temp = (U2/aden)**0.5
           IF (s2 == s2_min) eff = ef
           IF (ee-ef < 65.0 .AND. u2 > 0.0) THEN

               CALL get_level_density(U2,REAL(s2),aden,
     &              M_inertia_2,entropy,level_density)

               jj = jj + 1
               IF (jj > n_light) THEN
               WRITE (UNIT=*, FMT=*) "increase N_light in stat_mod.f90"
                   STOP 
               END IF
     
               storage(jj)%gam = level_density*temp*(2.0*s(mode)+1.0)
               storage(jj)%energy = REAL(l)
               storage(jj)%spin = s2
       IF (jj > 1) storage(jj)%gam = storage(jj)%gam + storage(jj-1)%gam
           END IF
       ELSE
           ef = ef0 - 5.0
           IF (s2 == s2_min) eff = ef0 - 5.0
       END IF
       s2 = s2 + 1
       END DO
       lms = lms + 1
	END DO
	IF (jj == 0) RETURN
	IF (storage(jj)%gam < 0.0) RETURN

	i = 1

	CALL RANDOM_NUMBER(prob)

	DO 
	   fract = storage(i)%gam/storage(jj)%gam
	   IF (fract > prob) EXIT
	   i = i + 1
	END DO

	l = NINT(storage(i)%energy)
	s2 = storage(i)%spin
	CALL yrast(Z2,A2,REAL(s2),0.0,0.0,E_rot)          ! rotational energy

	U2 = Ex0 - Ba - Ec - E_rot - REAL(l*(l+1))*const

	CALL get_little_a(A2,Z2,REAL(s2),U2,bs,shell_corr,E_rot,aden,
     &	 entropy, pair_energy,1.0)
	temp = (U2/aden)**0.5

	ef = entropy
	prob2 = 1.0
	prob = 0.0
	Ranges = MIN(4.0*temp,U2)

	DO
	   IF (prob2 < prob) EXIT
	   CALL RANDOM_NUMBER(ran)
	   ek = ranges*RAN
	   prob = EXP(2.0*(aden*(U2-ek))**0.5 - ef)
     &	 *(1.0+U2**2.0)/(1.0+(U2-ek)**2.0)
	   CALL RANDOM_NUMBER(prob2)
	END DO

	RETURN
	END SUBROUTINE light_particle

	SUBROUTINE light_particle_quick (mode,Ex0)
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN) :: EX0

	REAL (KIND=r4) :: t
	REAL (KIND=r4) :: fact1,fact2,fact3,bit,x_min,integrate,ss1,ss2
	REAL (KIND=r4) :: M_inertia_star,bb(200),s2_max,ss,eg

	M_inertia_star = 1.0/(1.0/M_inertia_orbit+1.0/M_inertia_2)

! loop over lms = l - s2
! l = orbital angular momentum of the emitted particle
! s2 = spin of residual
! lms is confined to limits -J0 to J0
	lms = -s0
	ef = ef0
	s2_min = s0 + 1

	DO
	  IF (lms > s0) EXIT
	  IF ((ef0-ef) > 4.0 .AND. s2_min < s2_best) EXIT
!  find minimum vales of J2 and l which are consistent
!  with this value of lms
	   l_min = INT(REAL(s0+lms+1)/2.0)
	   s2_min = l_min - lms
	   CALL yrast(Z2,A2,REAL(s2_min),0.0,0.0,E_rot) ! rotational energy of residue
	   U2 = Ex0 - Ba - Ec - E_rot - REAL(l_min*(l_min+1))*const
	   IF (U2 > 0.01) THEN

       CALL get_little_a(A2,Z2,REAL(s2_min),U2,bs,shell_corr,e_rot,aden,
     &                  entropy, pair_energy,1.0)
      ef = entropy
      IF (ee-ef < 65.0 .AND. u2 > 0.0) THEN
          temp = (U2/aden)**0.5
          fact1 = (M_inertia_2*(REAL(l_min)+0.5)
     &         +M_inertia_orbit*(REAL(s2_min)+0.5))**2.0    
     &         /M_inertia_2/M_inertia_orbit/M_inertia_saddle/temp*20.424
          fact2 = (M_inertia_star*temp/20.424)**0.5
          bit = (REAL(lms)+0.5)*M_inertia_star  
     &         /M_inertia_orbit + M_inertia_star/2.0/M_inertia_2
          x_min = (REAL(s2_min)-0.5+bit)/fact2

!         now we are going to approximate the integral
!         from x_min to infinity of EXP(-x**2.0) using
!         the approximation for the normal probability
!         from "Handbook of Mathematical functions",
!         by Abramowitz and Stegan, Page932, accurate
!         to 7.5e-8
!
          t = 1.0/(1.0+0.2316419*1.4142135*x_min)
          fact3 = t*0.31938 - t**2.0*0.356563782
     &    + t**3.0*1.781477937 - t**4.0*1.821255978 
     &    + t**5.0*1.330274429
          fact3 = fact3*0.3989423*1.77245385

          integrate = EXP(fact1-x_min**2.0)*(fact2**2.0 
     &         + fact2*(1.0-2.0*bit)*fact3)

!         check to see if the integration along s2 gives a
!         value greater than that for s2_min
          IF (integrate/(2.0*REAL(s2_min)+1.0) < 1.0) 
     &        integrate = 2.0*REAL(s2_min)+1.0
          jj = jj + 1

          IF (jj > n_light) THEN
              WRITE (UNIT=*, FMT=*) "increase N_light"
              STOP 
          END IF

          CALL get_level_density(U2,0.0,aden,M_inertia_2,
     &         entropy,level_density)

        storage(jj)%gam = temp*(2.0*s(mode)+1.0)*integrate*level_density

          bb(jj) = bit
          storage(jj)%spin = lms
       IF (jj > 1) storage(jj)%gam = storage(jj)%gam + storage(jj-1)%gam
	      END IF
	   ELSE
       ef = 0.0
	   END IF
	   lms = lms + 1
	END DO
	IF (jj == 0) RETURN
	IF (storage(jj)%gam <= 0.0) RETURN	!hg cheng

	fract = 0
	i = 0
	
C 		pause

C 		print*, 'jj=',jj, storage(jj)%gam		!hg cheng

!find the values of (l - s2) form the stored particle widths
	CALL RANDOM_NUMBER(prob)
	DO
	  IF  (fract > prob) EXIT

		IF(I >= JJ)	THEN					!HG CHENG: a bug when joining QMD

C 			PRINT*, JJ, 'I>=JJ'
			JJ = 0
C 			PAUSE
			RETURN
		ENDIF

	  i = i + 1
C 		print*, 'i',i, storage(i)%gam				!hg cheng

		if(storage(i)%gam < 0.0) then				!HG CHENG, A BUG
			jj = 0
C 			PRINT*,'storage(i)%gam <0 '				!HG CHENG
C 			PAUSE		
			return
		endif

	  fract = storage(i)%gam/storage(jj)%gam
	END DO
	lms = storage(i)%spin


C 		pause

! now find the range over which we will choose s2 the residue's spin
	l_min = INT(REAL(s0+lms+1)/2.0)
	s2_min = l_min - lms
	CALL yrast(Z2,A2,REAL(s2_min),0.0,0.0,E_rot)               ! rotational energy
	U2 = Ex0 - Ba - Ec - E_rot - REAL(l_min*(l_min+1))*const

	 CALL get_little_a(A2,Z2,REAL(s2_min),U2,bs,shell_corr,E_rot,aden,
     &      entropy,pair_energy,1.0)

	temp = (U2/aden)**0.5
	fact2 = (M_inertia_star*temp/20.454)**0.5
	s2_max = ((8.0*fact2**2.0+(2.0*bb(i)-1)**2.0)**0.5 
     &       -(2.0*bb(i)+1.0))/4.0
	s2_max = MAX(REAL(s2_min)-.5,s2_max)
	ef = ((s2_max+bb(i))/fact2)**2.0
	eg = 2.0*s2_max+1.0
! find range over which to pick s2
	ss = s2_max
	prob = 1.0
	DO
	  IF (prob < 0.05) EXIT
	  ss = ss + fact2/2.0
	  prob = (2.0*ss+1.0)/eg*EXP(ef-((ss+bb(i))/fact2)**2.0)
	END DO
	ss2 = ss
	ss1 = REAL(s2_min) - .5

! now select the final spin of the residue and the angular momentum
! of the particle with the provision that the final thermal excitation
! energy of the residue is positive
	DO
	   DO  
       CALL RANDOM_NUMBER(ran)
       ss = ss1 + (ss2-ss1)*RAN
       prob = (2.0*ss+1.0)/eg*EXP(ef-((ss+bb(i))/fact2)**2.0)
       CALL RANDOM_NUMBER(prob2)
       IF (prob2 < prob) EXIT
	   END DO
	   s2 = NINT(ss)
	   l = s2 + lms
	   CALL yrast(Z2,A2,REAL(s2),0.0,0.0,E_rot)   ! rotational energy
	   U2 = Ex0 - Ba - Ec - E_rot - REAL(l*(l+1))*const
	   IF (U2 >= 0.0) EXIT
	END DO

!now choose the kinetic energy of the emitted particle

	CALL get_little_a(A2,Z2,REAL(s2),U2,bs,shell_corr,e_rot,aden,
     &     entropy, pair_energy,1.0)

	temp = (U2/aden)**0.5
	ef = entropy
	Ranges = MIN(4.0*temp,U2)

	DO 
	  CALL RANDOM_NUMBER(ran)
	  ek = ranges*RAN
	  prob = EXP(2.0*(aden*(U2-ek))**0.5 - ef) 
     &	 *(1.0+U2**2.0)/(1.0+(U2-ek)**2.0)
	  CALL RANDOM_NUMBER(prob2)
	  IF (prob2 < prob) EXIT
	END DO

	RETURN
	END SUBROUTINE light_particle_quick


	SUBROUTINE prepare_evap ()
	INTEGER :: i

	i = 1
	DO
	   if (i > N_mode) EXIT
	   if (z1(i) >= para%z_imf_min) EXIT 
	   list_hf0(i) = i
	   i = i + 1
	END DO
	i = i + 1
	list_hf0(i) = -1

	END SUBROUTINE prepare_evap
	END MODULE evap
!------------------------------end:   MODULE evap-----------------------------------



!------------------------------begin:   MODULE evaporation-----------------------------------
	MODULE evaporation
	USE gem_kind
	USE light_stuff
	USE trig
	USE yrast_mod
	USE tl_stuff
	USE functs
	USE mass_stuff
	USE stat_mod
	USE spectra
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: Hauser_feshbach, prepare_hauser_feshbach,prepare_5li_5he

	INTEGER, PARAMETER, PUBLIC :: n_sub = 800
!INTEGER, PARAMETER, PUBLIC :: n_sub = 2000  !used for calculation with 12C

	TYPE (evap_store), DIMENSION(0:N_sub), PUBLIC :: sub_store

	INTEGER, PRIVATE :: Iz2,Ia2
	INTEGER, PRIVATE :: l,l_min,l_max,j_sub

	LOGICAL, PRIVATE :: odd,cg_store

	INTEGER, PRIVATE :: jj,kk,kkk
	REAL (KIND=r4), PRIVATE :: sum_tl,gam_jj
	REAL (KIND=r4), PRIVATE :: U2,aden,gam_s2_max,U2_min,ek_min,ee_ek,
     &	  sum_ek
	REAL (KIND=r4), PRIVATE :: gam_old
	REAL (KIND=r4), PRIVATE :: Z2,A2,r2,r1,Ba,E_rot,de,A23
	REAL (KIND=r4), PRIVATE :: gam_s2
	REAL (KIND=r4), PRIVATE :: fract,prob,ek
	REAL (KIND=r4), PRIVATE :: gam_ek_max
	REAL (KIND=r4), PRIVATE :: s0,s2,l_plus_s_min,l_plus_s_max,
     &	  s2_start,l_plus_s
	REAL (KIND=r4), PRIVATE :: ex1,ran,s1
	REAL (KIND=r4), PRIVATE :: shell_corr,pair_energy
	REAL (KIND=r4), PRIVATE :: entropy,level_density
	REAL (KIND=r4), PRIVATE :: k2,L_plus_s_min_k2
	REAL (KIND=r4),PRIVATE :: k2_min,k2_max,k2_start,gam_k2_max,sum_k2
	REAL (KIND=r4), PRIVATE :: M_inertia_2,M_inertia_orbit

	REAL (KIND=r4), PARAMETER, PRIVATE :: lower_ek = 0.001
	REAL (KIND=r4), PARAMETER, PRIVATE :: lower_ell = 0.001
	REAL (KIND=r4), PARAMETER, PRIVATE :: lower_s2 = 0.01
	REAL (KIND=r4), PARAMETER, PRIVATE :: lower_k2 = 0.01

	CONTAINS
	SUBROUTINE Hauser_Feshbach(mode,z0,a0,Ex0,J0,k0,ratio)
! 
! FUNCTIONAL DESCRIPTION:       
! 
! Calculates the decay-widths in MeV for light particle emission Z<3 from the
! Hauser-Feshbach formalism. The summation over spin and orbital angular
! momentum is performed.
! 
! DUMMY ARGUMENTS:
! 
! input: mode (INTEGER ::) 
!       Z0 (REAL (KIND=r4) ::) proton number of decaying nucleus
!      A0 (REAL (KIND=r4) ::) nucleon number of decaying nucleus
!      Ex0 (REAL (KIND=r4) ::) excitation energy (MeV) of decaying nucleus
!      J0 (REAL (KIND=r4) ::) spin of decaying nucleus
! 
! IMPLICIT INPUTS:
! 
! none
! 
! IMPLICIT OUTPUTS:
! 
! I_Channels - number of open light particle channels
! gamma(1:I_channels) - decay widths for each light particle channel
! Za(1:_I_channels) - proton number of light fragment
! Aa(1:I_channels)  - nucleon number of light fragment
! Ex_residue(1:I_channels) - excitation energy of daughter nucleus
! 
! 
! SIDE EFFECTS:
! 
!  none
! 

	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN)  :: J0,K0,Z0,A0,Ex0,ratio

	INTEGER :: try,i
	REAL (KIND=r4), DIMENSION(:), POINTER :: line_shape
	REAL (KIND=r4) :: Mass2,mass2ld

	z2 = z0 - z1(mode)
	a2 = a0 - A1(mode)
	IF (z1(mode) > z2) RETURN
	IF (A1(mode) > A2) RETURN
	IF (z2 >= A2) RETURN
	iz2 = z2
	ia2 = a2

	if (para%ssum) THEN  !explicit consideration of evaporated particles spin
	  s1 = s(mode)
	ELSE                 !only take into account the spin by a (2s+1) factor
                     !this is much faster especially for para%k_sum
	  s1 = 0.
	END IF 

	IF (MODULO(NINT(a2),2) /= 0 .AND. para%ssum) THEN    ! odd nuclei have 1/2 integer
	   odd = .TRUE.
	ELSE
	   odd = .FALSE.
	END IF

	call get_pairing(NINT(z2),NINT(A2),pair_energy)
	if (z0 <= z_shell) pair_energy = 0.0

	IF (z1(mode) > z2) RETURN
	IF (A1(mode) > A2) RETURN
	IF (z2 >= A2) RETURN
	A23 = A2**0.33333333
	R2 = r0*A23
	r1 = r0*a1(mode)*.33333333
	M_inertia_2 = 0.4*A2*R2**2
	M_inertia_orbit = A2*A1(mode)/A0*(r1+R2+2.0)**2.0
!_______________________________________________________________ 
! calculate the mass defect of the residual system
! experimental and finite range
!--------------------------------------------------------------

	CALL getmass(iZ2,iA2,Mass2,mass2ld)  
	IF (para%mass_option /= 1 .AND. z0 > Z_shell) THEN
	    shell_corr = mass2 - mass2ld - pair_energy
	ELSE
	    shell_corr = 0.0
	END IF

! prepares the get_tl routine to find scaling for
! the A dependence of the transmission coeff

	CALL prepare_tl(mode_gs(mode),iz2,ia2)
	call clear_look_up ()

!______________________________________________________________
! if the emitted particle is unstable, it has a finite
! width so first we choose a energy for this state.
!------------------------------------------------------------

	IF (width(mode) > 0.) THEN
	   ex1 = -1000000.0
	   DO
       IF (ex1 >= limit(mode)) EXIT
       CALL RANDOM_NUMBER(RAN)
       ex1 =  width(mode)/2.0*TAND(179.998*RAN-89.999)
	   END DO
	ELSE IF (width(mode) == 0.0) THEN
	   ex1 = 0.0
	ELSE
	   SELECT CASE (NINT(z1(mode)))
	   CASE (2)
        line_shape => he5
        de = 0.04
	   CASE (3)
        line_shape => li5
        de = 0.07
	   CASE (4)
        line_shape => be8
        de = 0.07        
	   END SELECT

	   i = 1
	   CALL RANDOM_NUMBER(prob)
	   DO
       IF (prob <= line_shape(i)) EXIT 
       i = i + 1
       IF (i > 100) EXIT
	   END DO
	   CALL RANDOM_NUMBER(RAN)
	   ex1 = (REAL(i)+RAN-0.5)*de 
	   if (ex1 == 0.0) ex1 = 0.01
	   ex1 = ex1 + limit(mode)
	END IF

!---------------------------------------------------------------- 
!  calculate the binding energy
!-----------------------------------------------------------------
	ba = Mass2 + def(mode) - Mass0 + ex1

	de = 1.0
!de = 0.25
	jj = 0
	s0 = J0

! _____________________________________________________________
! summation over spin of residual system the summation starts off at the
! spin of the decaying system and decreases the spin while that is possible
! and then comes back and sums over the larger spins
!--------------------------------------------------------------
	s2_start = AINT(M_inertia_2/(M_inertia_orbit+M_inertia_2)*J0)
	IF (odd) s2_start = s2_start + 0.5
	s2 = s2_start
	gam_s2_max = 0.0
	gam_old = 0.0
	gam_s2 = 0.0

	call s2_loop (mode,z0,a0,Ex0,k0,ratio)

!construct e-j map of daughter nuclei for topdraw
!if (mode == 7) call EJ_map (jj,Ex0,Ba)      

	IF (jj == 0) RETURN
	IF ( storage(jj)%GAM == 0.0) RETURN

	I_channels = I_channels + 1
	Za(I_channels) = NINT(Z1(mode))
	Aa(I_channels) = NINT(A1(mode))
	gamma(I_channels) = storage(jj)%gam
	 if (.NOT. para%ssum) gamma(I_Channels) =
     &	gamma(I_channels)*(2.0*s(mode)+1.0)
	 evap_mode(I_channels) = mode
!WRITE (UNIT=*, FMT=*) z1(mode),a1(mode),s1,def(mode),gamma(I_channels)
	 IF (mode == 0) gamma(I_channels) = gamma(I_channels)/3.    !dineutron
                                                       !is not a real reasonance
!IF (mode == 7) CALL low_energy_alpha(jj)

!choose the spin of the residual system 
	fract = 0.0
	i = 0
	CALL RANDOM_NUMBER(prob)
	DO
	    i = i + 1
	    fract = storage(i)%gam/storage(jj)%gam
	    IF (fract > prob) EXIT
	    if (i == jj) EXIT
	END DO
	kk = jj
	kkk = i

	s2 = REAL(storage(i)%spin)/2.0
	k2 = REAL(storage(i)%kk2)/2.0
	CALL yrast(Z2,A2,s2,k2,ratio,E_rot)
	U2_min = Ex0 - Ba - E_rot
	ek = storage(i)%energy 

	dau(I_channels)%Ex_residue = Ex0 - Ba - ek
	dau(I_channels)%J_residue = s2
	dau(I_channels)%k_residue = k2
	dau(I_channels)%s_particle = s1
	dau(I_channels)%level = level_no(mode)
	dau(I_channels)%excit = ex1

	IF (.NOT. para%I_angle) RETURN

!now find an l and l+s1 values
!consider particle"s spin in determining
!the l (orbital angular momentum limits)
	l_plus_s_min = ABS(S0-S2)        !minimum value of (l+s) vector
	l_plus_s_max = S0 + S2           !maximum value of (l+s) vector
	if (para%k_sum) l_plus_s_min = MAX(l_plus_s_min,ABS(k0-k2))

	l_min = NINT(l_plus_s_min-s1)    ! minimum value of l vector
	l_max = NINT(l_plus_s_max+s1)    !maximum value of l vector

	l_max = MIN(l_max,l_max_quantum) 

	IF (l_min < 0) l_min = NINT(MAX(s1-l_plus_s_max,0.0))


	try = 0
	cg_store = .FALSE.
	loop_ell2: DO
	    try = try + 1

	    IF (try >= 4) THEN
        ek = storage(i)%energy 
        U2 = u2_min - ek
	    ELSE
        u2 = -1.0
        ek = -1.0
        DO
           CALL RANDOM_NUMBER(ran)
           ek = storage(i)%energy + (1.0-2.0*ran)*de/2.0
           U2 = u2_min - ek
           IF (U2 > 0.0 .AND. ek > 0.0) EXIT
        END DO
	    END IF 

	    IF (.NOT.(ratio == 1.0 .OR. ratio == 0.0) .OR. para%k_sum)
     &     call clear_look_up ()

	    CALL ell_loop (mode,k0,jj)
	    IF(jj > 0 .AND. sum_tl > 0) EXIT loop_ell2
	    IF(try >= 5) EXIT loop_ell2
	END DO loop_ell2

	dau(I_channels)%Ex_residue = Ex0 - Ba - ek
	dau(I_channels)%e_kinetic = ek 

	DO i=2,jj
	     sub_store(i)%gam = sub_store(i)%gam + sub_store(i-1)%gam
	END DO 

	IF (sub_store(jj)%gam == 0.0) THEN
	  WRITE (UNIT=*, FMT=*) "sub_store=0, jj=",jj
	  WRITE (UNIT=*, FMT=*) "mode=",mode
	  WRITE (UNIT=*, FMT=*) "tries=",try
	END IF 

	IF (jj == 0 .OR. sum_tl == 0.0) THEN
	   WRITE (UNIT=*, FMT=*) "z0=",Z0
	   WRITE (UNIT=*, FMT=*) "A0=",A0
	   WRITE (UNIT=*, FMT=*) "s0=",S0
	   WRITE (UNIT=*, FMT=*) "k0=",k0
	   WRITE (UNIT=*, FMT=*) "Ex0=",Ex0
	   WRITE (UNIT=*, FMT=*) "Erot=",E_rot
	   WRITE (UNIT=*, FMT=*) "U2=",u2
	   WRITE (UNIT=*, FMT=*) "ek=",ek
	   WRITE (UNIT=*, FMT=*) "s2=",s2
	   WRITE (UNIT=*, FMT=*) "jj=",kk
	   WRITE (UNIT=*, FMT=*) "kkk=",kkk
	   WRITE (UNIT=*, FMT=*) "mode=",mode
	   WRITE (UNIT=*, FMT=*) "de=",de
	   WRITE (UNIT=*, FMT=*) "jj=",jj
	   WRITE (UNIT=*, FMT=*) "sum_tl=",sum_tl
	END IF

!choose the spin of the residual system 
	i = 1      
	CALL RANDOM_NUMBER(prob)
	DO 
	    fract = sub_store(i)%gam/sub_store(jj)%gam
	    IF (fract > prob) EXIT
	    i = i + 1
	END DO
	l = sub_store(i)%spin

	dau(I_channels)%ll  = l
	IF (para%K_sum .OR.
     &	    para%polarization .AND. (mode == 2 .OR. mode == 1)) THEN
	    dau(I_channels)%l_p_s =  sub_store(i)%lps
	    dau(I_channels)%k_particle = REAL(sub_store(i)%kk2)/2.0
	ELSE
	    l_plus_s_max = MIN(l+s1,l_plus_s_max)
	    l_plus_s_min = MAX(ABS(l-s1),l_plus_s_min)
	    CALL RANDOM_NUMBER(RAN)
	    l_plus_s = INT(RAN*(l_plus_s_max + 1 - l_plus_s_min))
     &    + l_plus_s_min
	    IF (l_plus_s > l_plus_s_max) l_plus_s = l_plus_s_max
	    dau(I_channels)%l_p_s = l_Plus_s
	END IF
	RETURN

	END SUBROUTINE Hauser_Feshbach


	SUBROUTINE s2_loop (mode,z0,a0,Ex0,k0,ratio)
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0,ratio
	INTEGER :: increase_s2

! sum over the spin of the daughter
	 increase_decrease_S2: DO increase_s2 = -1,1,2      !increase or decrease S2
	    S2 = s2_start
	    loop_s2: DO                                     !summation over s2 residue spin

!   consider particle"s spin in determining
!   the l (orbital angular momentum limits)
        l_plus_s_min = ABS(S0-S2)   !minimum value of (l+s) vector
        l_plus_s_max = S0 + S2       !maximum value of (l+s) vector

        l_max = NINT(l_plus_s_max+s1)      !maximum value of l vector
        l_max = MIN(l_max,l_max_quantum) 
        l_min = NINT(l_plus_s_min-s1)      ! minimum value
                                           ! of l 
        IF (l_min < 0) THEN
            l_min = NINT(MAX(s1-l_plus_s_max,0.0))
        END IF

		IF(L_MIN > L_MAX)	L_MIN = L_MAX	!HG CHENG FOR QMD


        k2_if: if (para%k_sum) THEN
            l_plus_s_min_k2 = l_plus_s_min
            CALL k2_loop (mode,z0,a0,Ex0,k0,ratio)
            sum_ek = sum_k2
            if (u2_min <= de/2.0) EXIT loop_s2
        ELSE k2_if
            CALL yrast(Z2,A2,s2,0.0,ratio,E_rot)
            U2_min = Ex0 - Ba - E_rot
            CALL ek_loop (mode,z0,a0,Ex0,k0)
            IF (U2_min <= de/2.0) EXIT loop_s2
        END IF  k2_if

        gam_old = gam_s2 
        gam_s2 = sum_ek


        gam_s2_max = MAX(gam_s2_max,gam_s2)
        if (gam_s2_max == 0) EXIT loop_s2

        IF (gam_s2 < lower_s2*gam_s2_max) EXIT loop_s2


    !   Forces the exit of s2 loop incase
    !   there is a run away down to s2=0
        IF(.NOT. para%k_sum .AND. increase_s2 == -1 .AND. gam_s2>gam_old
     &     .AND. gam_old > 0 .AND. gam_s2 < gam_s2_max)  EXIT loop_s2

        s2 = s2 + increase_s2
        IF (S2 < 0) EXIT loop_s2
	    END DO loop_s2
	    s2_start = s2_start - increase_s2
	    if (s2_start < 0) EXIT increase_decrease_s2
	END DO increase_decrease_s2
	RETURN
	END SUBROUTINE s2_loop

	SUBROUTINE ek_loop (mode,z0,a0,Ex0,k0)
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0
	REAL (KIND=r4) :: ekH,deH,recoil
	INTEGER :: increase_ek

	sum_ek = 0.0
	gam_ek_max = 0.0

!find energy at which tl=0.5 approximately
	 ek_min = REAL((l_min+1)*l_min)*20.454/
     &	m_inertia_orbit + 1.44*z1(mode)*z2/(R2+R1+2)

!start just below the barrier
	ek_min = 0.8*ek_min
	ek_min = REAL(NINT(ek_min+de/2.0)) - de/2.0
	ek_min = MAX(ek_min,de/2.0)
	IF (ek_min >= U2_min) RETURN !exit loop2

! integrate over kinetic energy of the emitted particle
	increase_decrease_ek: DO increase_ek=1,-1,-2  !increase

! or decrease the kinetic energy
	    ek = ek_min 
	    U2 = U2_min - ek
	    loop_ek: DO !loop over ek
!       calculate the sum of the transmission coefficients
        CALL ell_loop(mode,k0,j_sub)
        IF (sum_tl <= 0.0) EXIT loop_ek

!       calculate the level density parameter
        CALL get_little_a(A2,Z2,S2,U2,bs,shell_corr,E_rot,
     &       aden,entropy,pair_energy,1.0)

        IF (u2 <= 0) EXIT loop_ek
        ee_ek = entropy - ee
        IF (ABS(ee_ek) > 65.0 ) EXIT loop_ek


        CALL get_level_density(U2,s2,aden,M_Inertia_2,
     &       entropy,level_density)

        CALL collective_enhancement(iz2,ia2,U2,aden,m_inertia_2,
     &	   level_density)

        gam_jj = level_density*de*sum_tl
        
        if (gam_jj <= 0) EXIT loop_ek
        jj = jj + 1

!        check if array is already fill
        IF (jj > N_light) THEN
            WRITE (UNIT=*, FMT=*) "z0=",Z0
            WRITE (UNIT=*, FMT=*) "A0=",A0
            WRITE (UNIT=*, FMT=*) "s0=",S0
            WRITE (UNIT=*, FMT=*) "Ex0=",Ex0
            WRITE (UNIT=*, FMT=*) "mode=",mode
        WRITE (UNIT=*, FMT=*) "increase N_light_normal or n_light_ksum"
        WRITE (UNIT=*, FMT=*) "n_light=",n_light,"in stat_mod.f90"
            STOP 
        END IF

        !histograms for n,p alpha
        !recoil = (A0-A1(mode))/A0
        !deH = de*recoil
        !ekH = ek*recoil
        !if (A0 > 0.75*para%Acn)CALL hist_ek(ekH,deH,mode,gam_jj)


        storage(jj)%gam = gam_jj
        gam_ek_max = MAX(storage(jj)%gam,gam_ek_max)
        sum_ek = sum_ek + storage(jj)%gam
        
        storage(jj)%spin = NINT(2.0*s2)
        storage(jj)%kk2 =  NINT(2.0*k2)
        storage(jj)%energy = ek
       IF (jj > 1) storage(jj)%gam = storage(jj)%gam + storage(jj-1)%gam
       IF (gam_jj < lower_ek*gam_ek_max) EXIT loop_ek

        ek = ek + increase_ek*de
        IF (ek <= 0) EXIT loop_ek
        U2 = U2_min - ek
        IF (U2 < 0.01) EXIT loop_ek
	   END DO loop_ek                ! end loop over ek
	    ek_min = ek_min - increase_ek*de
	    if (ek_min <= 0) EXIT increase_decrease_ek
	END DO  increase_decrease_ek
	RETURN
	END SUBROUTINE ek_loop

	SUBROUTINE ell_loop (mode,k0,j_sub)
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN)  :: K0
	INTEGER, INTENT(INOUT) :: j_sub

!INTEGER, PARAMETER :: n_store = 400
	INTEGER, PARAMETER :: n_store = 1000      !used for calculation with 12C
	INTEGER :: i_store,i_cg
	INTEGER :: kl,kl_min,kl_max

	REAL (KIND=r4) :: weight_s,tlw
	REAL (KIND=r4) :: max_l_plus_s,min_l_plus_s
	REAL (KIND=r4) :: k1,l_plus_s,k_l_plus_s
	REAL (KIND=r4) :: clebsch_gordan1,clebsch_gordan2
	REAL (KIND=r4) :: tlw_max,tl_up,tl_down
	REAL (KIND=r4), DIMENSION(N_store) :: clebsch_gordan_store
	REAL (KIND=r4), DIMENSION(0:l_max_quantum) :: tl


!sum the transmission coefficients over the allowed range of orbital angular momentum
	if (l_min > l_max) THEN
	    WRITE (*,*) "l_min > lmax",l_min,l_max 
	    STOP
	END IF 
	l = l_min
	sum_tl = 0.0
	tlw_max = 0.0
	j_sub = 0
	if (u2 <= 0.0)RETURN
	k_l_plus_s = k0 - k2

	i_cg = 0

	loop_ell: DO  !summation over l orbital angular momentum

!   go get the transmission coefficients
	    CALL get_tl(l,ek,tl,tl_up,tl_down)

!   find the number of a particular l from the quantum sum rules
	    max_l_plus_s = MIN(l+s1,l_plus_s_max)
	    min_l_plus_s = MAX(ABS(l-s1),l_plus_s_min)

	    if (para%k_sum) THEN

        kl_min = MAX(NINT(k0-k2-s1),-l)
        kl_max = MIN(NINT(k0-k2+s1),l)

        weight_s = 0.0
        l_plus_s = min_l_plus_s
        loop_l_plus_s: DO 

            i_cg = i_cg + 1
            if (cg_store .AND. i_cg <= i_store) THEN
                clebsch_gordan1 = clebsch_gordan_store(i_cg)
            ELSE
                clebsch_gordan1 = Cl_go(NINT(s2*2),NINT(l_plus_s*2.0),
     &    NINT(S0*2.0),NINT(k2*2.0),NINT(k_l_plus_s*2.0),NINT(k0*2.0))
                if (i_cg <= n_store)
     &          clebsch_gordan_store(i_cg) = clebsch_gordan1
            END IF 

            if (clebsch_gordan1 == 0) EXIT loop_l_plus_s

            loop_kell: DO kl = kl_min,kl_max
                k1 = k0-k2-REAL(kl)

                i_cg = i_cg + 1
                if (cg_store .AND. i_cg <= i_store) THEN
                    clebsch_gordan2 = clebsch_gordan_store(i_cg)
                ELSE
       Clebsch_gordan2 = Cl_go(l*2,NINT(s1*2.0),NINT(l_plus_s*2.0),
     &                    kl*2,NINT(k1*2.0),NINT(k_l_plus_s*2.0))
       if (i_cg <= n_store) clebsch_gordan_store(i_cg) = clebsch_gordan2
                END IF 

                weight_s = weight_s + tl(ABS(kl))*
     &                (clebsch_gordan1*clebsch_gordan2)**2.0

                j_sub = j_sub + 1
                sub_store(j_sub)%gam = tl(ABS(kl))*
     &          (clebsch_gordan1*clebsch_gordan2)**2.0
                sub_store(j_sub)%spin = l
                sub_store(j_sub)%lps = l_plus_s                   
                sub_store(j_sub)%kk2 = NINT(2.0*k1)
            END DO loop_kell
            l_plus_s = l_plus_s + 1.0
            IF (l_plus_s > max_l_plus_s) EXIT loop_l_plus_s
        END DO loop_l_plus_s
        
        tlw = weight_s

	    ELSE

        weight_s = max_l_plus_s - min_l_plus_s + 1
        tlw = tl(0)*weight_s
        j_sub = j_sub + 1

        IF ((mode == 2 .or. mode == 1)
     &        .AND. weight_s >= 1.0 .AND. para%polarization) THEN
            IF (weight_s == 2) THEN
                tlw = tl_up + tl_down
                sub_store(j_sub)%gam = tl_up
                sub_store(j_sub)%spin = l
                sub_store(j_sub)%lps = l + s1
                j_sub = j_sub + 1
                sub_store(j_sub)%gam = tl_down
                sub_store(j_sub)%spin = l
                sub_store(j_sub)%lps = l - s1
            ELSE IF (max_l_plus_s == l + s1) THEN
                tlw = tl_up
                sub_store(j_sub)%gam = tl_up
                sub_store(j_sub)%spin = l
                sub_store(j_sub)%lps = l + s1
            ELSE
                tlw = tl_down
                 sub_store(j_sub)%gam = tl_down
                 sub_store(j_sub)%spin = l
                 sub_store(j_sub)%lps = l - s1
            END IF
        ELSE
            sub_store(j_sub)%gam = tl(0)*weight_s
            sub_store(j_sub)%spin = l
        END IF
       END IF 

       sum_tl = sum_tl + tlw

       tlw_max = MAX(tlw,tlw_max)

       IF (tlw < lower_ell*tlw_max) EXIT loop_ell

       l = l + 1

       IF (l > l_max) EXIT loop_ell
	 END DO loop_ell         !finish sum over orbital angular momentum

	if (j_sub > n_sub) THEN
	    WRITE (UNIT=*, FMT=*) "increase n_sub in light_particle.f90"
	    STOP
	END IF 

	if (.NOT.  para%k_sum) RETURN
	if (cg_store .AND. i_cg > i_store) THEN
	   if (I_cg <= n_store)THEN
       i_store = i_cg
	   ELSE
       i_store = n_store
       WRITE (UNIT=*, FMT=*) 'i_cg=',i_cg
	   END IF 
	ELSE if (.NOT. cg_store) THEN
	   if (I_cg <= n_store)THEN
       i_store = i_cg
	   ELSE
       i_store = n_store
	   END IF 
	   cg_store = .TRUE.
	END IF 

	RETURN
	END SUBROUTINE ell_loop

	SUBROUTINE k2_loop (mode,z0,a0,Ex0,k0,ratio)
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN)  :: K0,Z0,A0,Ex0,ratio

	INTEGER :: increase_k2

!we have limited the maxiumum ell possible - can change this
!in gemini.inc, therefore we have further restrictions on
!l_plus_s

	l_plus_s_max  = min(l_plus_s_max,l_max_quantum+s1)
	k2_max = k0 + l_plus_s_max
	if (k2_max > s2) k2_max = s2
	k2_min = k0 - l_plus_s_max
	if (k2_min < -s2) k2_min = -s2

	IF (S0 == 0.0) then
	    K2_START = 0.0
	else
	    k2_start = k0/s0*s2
	end if     

	if (odd) THEN
	    k2_start = AINT(k2_start+0.5)-0.5
	ELSE
	    k2_start = AINT(k2_start)
	END IF
	if (k2_start > k2_max) k2_start = k2_max
	if (k2_start < k2_min) k2_start = k2_min

	gam_k2_max = 0.0
	sum_k2 = 0.0
	increase_decrease_k2: DO increase_k2 = -1,1,2  !increase or
	                                               !decrease k2     
	    k2 = k2_start
	    loop_k2: DO !summation over k2

	    CALL yrast(Z2,A2,s2,k2,ratio,E_rot)
	    U2_min = Ex0 - Ba - E_rot

    !there is also a minimum value of l_plus_s for which is
    !possible to account for the differences between k0 and k2
	    l_plus_s_min = MAX(l_plus_s_min_k2,ABS(k0-k2)) 
	    if (l_plus_s_min_k2 > l_plus_s_max) EXIT loop_k2 

	    l_max = NINT(l_plus_s_max+s1)           !maximum value of l vector
    !l_max futher restricted by max L for which transmission coeff
    !are calculated - l_max_quantum - can be chaned in gemini.inc
	    l_max = MIN(l_max,l_max_quantum) 

	    l_min = NINT(l_plus_s_min-s1)       ! minimum value of l
	    IF (l_min < 0) THEN
		l_min = NINT(MAX(s1-l_plus_s_max,0.0))
	    END IF

	    cg_store = .FALSE.
	    CALL ek_loop (mode,z0,a0,Ex0,k0)
	    IF (U2_min <= 0.1) RETURN          !EXIT loop_s2
    
	    sum_k2 = sum_k2 + sum_ek

!    if (mode == 7) call hist_2d(INT(s2),INT(k2),sum_ek,1)
	    gam_k2_max = MAX(gam_k2_max,sum_ek)
	    if (gam_k2_max == 0) EXIT loop_k2
	    IF (sum_ek < lower_k2*gam_k2_max) EXIT loop_k2

	    k2 = k2 + increase_k2
	    IF (k2 < k2_min) EXIT loop_k2
	    if (k2 > k2_max) EXIT loop_k2
	END DO loop_k2
	k2_start = k2_start - increase_k2
	if (k2_start > k2_max) EXIT increase_decrease_k2
	if (k2_start < k2_min) EXIT increase_decrease_k2

	END DO increase_decrease_k2
	END SUBROUTINE k2_loop



	SUBROUTINE prepare_Hauser_Feshbach ()
	INTEGER :: jj,index
	jj = 0
	DO index=1,n_mode
	    IF (Z1(index) >= para%Z_imf_min)EXIT
	    IF (para%exotic_index >= exotic(index))THEN
        jj = jj + 1
        list_HF0(jj) = index
	 END IF  
	END DO
	jj = jj + 1
	list_hf0(jj) = -1

	RETURN
	END SUBROUTINE prepare_Hauser_Feshbach

	SUBROUTINE prepare_5li_5he ()
!initialization
	REAL (KIND=r4) :: sum_he5,sum_li5,sum_be8
	INTEGER :: i

	sum_he5 = SUM(he5)
	sum_li5 = SUM(li5)
	sum_be8 = SUM(be8)
	DO i=1,100
	    he5(i) = he5(i)/sum_he5
	    li5(i) = li5(i)/sum_li5
	    be8(i) = be8(i)/sum_be8
	    IF (i > 1) THEN
        he5(i) = he5(i) + he5(i-1)
        li5(i) = li5(i) + li5(i-1)
        be8(i) = be8(i) + be8(i-1)
	    END IF
	END DO
	RETURN
	END SUBROUTINE prepare_5li_5he

	END MODULE evaporation
!------------------------------end:   MODULE evaporation-----------------------------------



!------------------------------begin:   MODULE nonstatistical_mod-----------------------------------
	MODULE nonstatistical_mod
	USE gem_kind
	USE angle
	USE functs
	USE trig
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: non_statistical 
	PRIVATE :: prepare_nonstatistical

	REAL (KIND=r4), DIMENSION(68), SAVE :: BE6_1 = 
     &  (/2.0270086E-09, 5.8145673E-07, 7.1854870E-06, 3.2172113E-05,  
     &    8.9211062E-05, 1.8850228E-04, 3.3525753E-04, 5.3003017E-04,  
     &    7.6970848E-04, 1.0486309E-03, 1.3595871E-03, 1.6944309E-03,  
     &    2.0449837E-03, 2.4032623E-03, 2.7617388E-03, 3.1135539E-03,  
     &    3.4528391E-03, 3.7744523E-03, 4.0741344E-03, 4.3484448E-03,  
     &    4.5947176E-03, 4.8111109E-03, 4.9964017E-03, 5.1498408E-03,  
     &    5.2712206E-03, 5.3608082E-03, 5.4191947E-03, 5.4473905E-03,  
     &    5.4465602E-03, 5.4181116E-03, 5.3636283E-03, 5.2847988E-03,  
     &    5.1833903E-03, 5.0613396E-03, 4.9207220E-03, 4.7632651E-03,  
     &    4.5908946E-03, 4.4055441E-03, 4.2090644E-03, 4.0033292E-03,  
     &    3.7900696E-03, 3.5710568E-03, 3.3479610E-03, 3.1227483E-03,  
     &    2.8964968E-03, 2.6708802E-03, 2.4473441E-03, 2.2272691E-03, 
     &    2.0119594E-03, 1.8026819E-03, 1.6005816E-03, 1.4068034E-03,  
     &    1.2223840E-03, 1.0483334E-03, 8.8556873E-04, 7.3493982E-04,  
     &    5.9724285E-04, 4.7317994E-04, 3.6334124E-04, 2.6820163E-04,  
     &    1.8810487E-04, 1.2316174E-04, 7.3201008E-05, 3.7638565E-05,  
     &    1.5269745E-05, 3.9710330E-06, 3.5455668E-07, 2.4267577E-10/)

	 REAL (KIND=r4), DIMENSION(41), SAVE :: he6_1 = 
     &  (/5.9025350E-04, 1.4816680E-03, 2.4142121E-03, 3.2942002E-03,  
     &  4.0769968E-03, 4.7419108E-03, 5.2821259E-03, 5.6987219E-03,  
     &  5.9972950E-03, 6.1864676E-03, 6.2764324E-03, 6.2776804E-03,  
     &  6.2013729E-03, 6.0580783E-03, 5.8579575E-03, 5.6113214E-03,  
     &  5.3261835E-03, 5.0112591E-03, 4.6740063E-03, 4.3213200E-03,  
     &  3.9594155E-03, 3.5939862E-03, 3.2301033E-03, 2.8722088E-03,  
     &  2.5243657E-03, 2.1901692E-03, 1.8727884E-03, 1.5749880E-03,  
     &  1.2991506E-03, 1.0474186E-03, 8.2146755E-04, 6.2264630E-04,  
     &  4.5202073E-04, 3.1013964E-04, 1.9707170E-04, 1.1218890E-04,  
     &  5.3955449E-05, 1.9486508E-05, 3.9684455E-06, 1.6708732E-07,  
     &  8.2843659E-15/)

	 REAL (KIND=r4), DIMENSION(85), SAVE :: be9_2 = 
     &  (/8.2305390E-26, 2.8723518E-18, 6.4047232E-15, 6.3987297E-13, 
     &  1.4879972E-11, 1.5212331E-10, 9.2700458E-10, 3.9760013E-09, 
     &  1.3261308E-08, 3.6661554E-08, 8.7709815E-08, 1.8716629E-07, 
     &  3.6419763E-07, 6.5686294E-07, 1.1119278E-06, 1.7841808E-06, 
     &  2.7345291E-06, 4.0288242E-06, 5.7347315E-06, 7.9200827E-06, 
     &  1.0650383E-05, 1.3986260E-05, 1.7982296E-05, 2.2684633E-05, 
     &  2.8128959E-05, 3.4341094E-05, 4.1336891E-05, 4.9120306E-05, 
     &  5.7680518E-05, 6.7000161E-05, 7.7047647E-05, 8.7784290E-05, 
     &  9.9158489E-05, 1.1111212E-04, 1.2357711E-04, 1.3648247E-04, 
     &  1.4974950E-04, 1.6329350E-04, 1.7702780E-04, 1.9086317E-04, 
     &  2.0470416E-04, 2.1845852E-04, 2.3203198E-04, 2.4532984E-04, 
     &  2.5825755E-04, 2.7072892E-04, 2.8265314E-04, 2.9394246E-04, 
     &  3.0451489E-04, 3.1429125E-04, 3.2319827E-04, 3.3116009E-04, 
     &  3.3811448E-04, 3.4400052E-04, 3.4876211E-04, 3.5234989E-04, 
     &  3.5471984E-04, 3.5583484E-04, 3.5566138E-04, 3.5417802E-04, 
     &  3.5136752E-04, 3.4721725E-04, 3.4172178E-04, 3.3488654E-04, 
     &  3.2672478E-04, 3.1725483E-04, 3.0650699E-04, 2.9451671E-04, 
     &  2.8132813E-04, 2.6699918E-04, 2.5159851E-04, 2.3520501E-04, 
     &  2.1790304E-04, 1.9980530E-04, 1.8102973E-04, 1.6171398E-04, 
     &  1.4202381E-04, 1.2212848E-04, 1.0225191E-04, 8.2647246E-05, 
     &  6.3622749E-05, 4.5560981E-05, 2.8965014E-05, 1.4560422E-05, 
     &  3.6027143E-06/)

	  TYPE (vell), PARAMETER :: blank_vell = vell(0.0,0.0,0.0)
	  TYPE (spinn), PARAMETER :: blank_spinn = spinn(0.0,0.0)
	  TYPE (identity), PARAMETER :: neutron = 
     &  identity(0,1,0,0,0.0,0.0,0.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: proton = 
     &  identity(1,1,0,0,0.0,0.0,0.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: deuteron = 
     &  identity(1,2,0,0,0.0,0.0,1.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: triton = 
     &  identity(1,3,0,0,0.0,0.0,0.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: he3 = 
     &  identity(2,3,0,0,0.0,0.0,0.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: alpha = 
     &  identity(2,4,0,0,0.0,0.0,0.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: M1 = 
     &identity(-1,11,0,0,0.0,0.0,1.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: E2 = 
     & identity(-1,2,0,0,0.0,0.0,2.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: E1 = 
     & identity(-1,1,0,0,0.0,0.0,1.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: he6 = 
     &  identity(2,6,0,0,0.0,0.0,0.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: he7 = 
     &  identity(2,7,0,0,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: li6 = 
     &  identity(3,6,0,0,0.0,0.0,1.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: li7 = 
     &  identity(3,7,0,0,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: li8 = 
     &  identity(3,8,0,0,0.0,0.0,2.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: be7 = 
     &  identity(4,7,0,0,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: be8 = 
     &  identity(4,8,0,1,0.0,0.0,0.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: be9 = 
     &  identity(4,9,0,0,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: be10 = 
     & identity(4,10,0,0,0.0,0.0,0.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: li5 = 
     &  identity(3,5,0,1,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: he5 = 
     &  identity(2,5,0,1,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: n2 = 
     &  identity(0,2,0,1,0.0,0.0,0.0,0.0,0.0,0.0,blank_spinn,blank_vell)

	  TYPE (identity), PARAMETER :: n3 = 
     &  identity(0,3,0,1,0.0,0.0,1.5,0.0,0.0,0.0,blank_spinn,blank_vell)

	CONTAINS

	SUBROUTINE prepare_nonstatistical ()
	    REAL (KIND=r4) :: tot
	    INTEGER :: i
    
	    tot = SUM(be6_1)
	    DO I=1,68
        be6_1(i) = be6_1(i)/tot
        IF (i > 2) be6_1(i) = be6_1(i) + be6_1(i-1)
	    END DO

	    tot = SUM(be9_2)
	    DO I=1,85
        be9_2(i) = be9_2(i)/tot
        IF (i > 2) be9_2(i) = be9_2(i) + be9_2(i-1)
	    END DO

	    tot = SUM(he6_1)
	    DO I=1,41
        He6_1(i) = He6_1(i)/tot
        IF (i > 2) He6_1(i) = He6_1(i) + He6_1(i-1)
	    END DO
	END SUBROUTINE prepare_nonstatistical 


	SUBROUTINE non_statistical(parent,daughter,n_extra)
        
!+
! 
! Does the decay of specific know levels of light nuclei.
!
! 3-body decays, He6(1.8 Mev), Be6, Li6(5.36) treats as
! two 2-body decays following Lane + Thomas,
! Rev Modern Physics, vol 30 1958 page=257
! R matrix parameters for He5 and Li5 intermediates
! from Stammbach and Walters, nucl phys A180 (1960) page=225.
! FUNCTIONAL DESCRIPTION: 
! 
!-
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), DIMENSION(3), INTENT(out) :: daughter
	INTEGER, INTENT(out) :: n_extra
	INTEGER :: i,j,ll
	LOGICAL, save :: first = .TRUE.
	REAL (KIND=r4) :: gamma_width,decay_time,vrel,v1,v2,ared,theta,phi
	REAL (KIND=r4) :: e,prob,ran
	REAL (KIND=r4), DIMENSION(3,3) :: v_f
	REAL (KIND=r4), DIMENSION(3) :: vv_f,p_f
	REAL (KIND=r4) :: e_kinetic,l_p_s

	IF (first) THEN
	    CALL prepare_nonstatistical ()
	    first = .FALSE.
	END IF

	IF (parent%level <= 0) RETURN


	frag_Z: SELECT CASE (parent%Z)  ! select z
	    CASE (0)           ! dineutrons etc
	    frag_A_n: SELECT CASE (parent%A)
	    CASE (2)       !dineutron dineutron dineutron
        daughter(1) = neutron 
        daughter(2) = neutron
        gamma_width = 1
        e_kinetic = 0.1
        n_extra = 2
        ll = 0
        l_p_s = 0.5
	    CASE (3)            !trineutron trineutron trineutron
        daughter(1) = neutron
        daughter(2) = n2
        gamma_width = 1
        e_kinetic = 0.1
        n_extra = 2
        ll = 1
        l_p_s = 1.5
	    CASE (4)               ! 4n 4n 4n 4n 4n 4n 4n 4n 4n 4n 4n
        daughter(1) = neutron
        daughter(2) = n3
        gamma_width = 1
        e_kinetic = 0.1
        n_extra = 2
        ll = 1
        l_p_s = 1.5
	    END SELECT frag_a_n
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
	CASE (1)
! 4H  4H  4H  4H  4H  4H  4H  4H  4H  4H  4H  4H  4H
	    daughter(1) = neutron
	    daughter(2) = triton
	    gamma_width = 1
	    e_kinetic = 2.9 + parent%ex
	    n_extra = 2
	    ll = 1
	    l_p_s = 1.5

!HEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEH
	CASE (2)  !He fragments
	    frag_a_He: SELECT CASE (parent%A)

	    CASE (4) ! 4He 4He 4He 4He 4He 4He 4He 4He 4He 4He 4He 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = proton
            daughter(2) = triton
            gamma_width = 0.5
            e_kinetic = 0.3950 + parent%ex
            n_extra = 2
            ll = 0
            l_p_s = 0.5
        CASE (2)
            n_extra = 2
            gamma_width = 0.8
            ll = 0
            l_p_s = 0.5
            e_kinetic = 0.0
            CALL RANDOM_NUMBER(ran)
            IF (RAN <= 0.2381) THEN
                e_kinetic = 0.4320 + parent%ex
            END IF
            IF (e_kinetic > 0) THEN
                daughter(1) = neutron
                daughter(2) = he3
            ELSE
                daughter(1) = proton
                daughter(2) = triton
                e_kinetic = 1.195 + parent%ex
            END IF
        CASE (3)
            n_extra = 2
            gamma_width = 2.01
            ll = 1
            l_p_s = 1.5
            e_kinetic = 0.0
            CALL RANDOM_NUMBER(ran)
            IF (RAN <= 0.3731) THEN
                e_kinetic = 1.262 + parent%ex
            END IF
            IF (e_kinetic > 0) THEN
                daughter(1) = neutron
                daughter(2) = he3
            ELSE
                daughter(1) = proton
                daughter(2) = triton
                e_kinetic = 2.025 + parent%ex
            END IF
        CASE (4)
            n_extra = 2
            gamma_width = 5.01
            ll = 1
            l_p_s = 1.5
            e_kinetic = 0.0
            CALL RANDOM_NUMBER(ran)
            IF (RAN <= 0.4741) THEN
                e_kinetic = 2.75 + parent%ex
            END IF
            IF (e_kinetic > 0) THEN
                daughter(1) = neutron
                daughter(2) = he3
            ELSE
                daughter(1) = proton
                daughter(2) = triton
                e_kinetic = 3.515 + parent%ex
            END IF

       END SELECT
	    CASE (5) ! 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 5He 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = neutron
            daughter(2) = alpha
            gamma_width = 0.6
            e_kinetic = parent%ex + 0.89
            if (e_kinetic <= 0) WRITE (UNIT=*, FMT=*) "he5",parent%ex
            n_extra = 2
            ll = 1
            l_p_s = 1.5
        END SELECT
	    CASE (6) ! 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 6He 
        SELECT CASE (parent%level)
        CASE (1)
            CALL RANDOM_NUMBER(prob)
            i = 1
            DO
                IF (prob <= he6_1(i)) EXIT
                i = i + 1
                IF (i >= 41) EXIT
            END DO
            e_kinetic = REAL(i)*0.02
            daughter(1) = neutron
            daughter(2) = he5
            daughter(2)%ex = 0.8245 - e_kinetic -0.89
            gamma_width = 0.113
            ll = 1
            l_p_s = 0.5 
            n_extra = 2
        END SELECT
	    CASE (7) ! 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 7He 
        daughter(1) = neutron
        daughter(2) = he6
        gamma_width = 0.16
        e_kinetic = 0.445 + parent%ex
        n_extra = 2
        ll = 1
        l_p_s = 1.5
	  CASE (8) ! 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 8He 
        daughter(1) = neutron
        daughter(2) = he7
        gamma_width = 0.8
        e_kinetic = 1.01 + parent%ex
        n_extra = 2
        ll = 1
        l_p_s = 1.5
	    END SELECT frag_A_He
!LILILILILILILILILILILILILILILILILILILILILILILILILILILI
	 CASE (3)  !li fragments                     
	    frag_A_Li: SELECT CASE (parent%A)
	    CASE (5) ! 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 5Li 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = proton
            daughter(2) = alpha
            gamma_width = 1.5
            e_kinetic = 1.96 + parent%ex
            n_extra = 2
            ll = 1
            l_p_s = 1.5
        END SELECT 
	    CASE (6) ! 6LI 6Li 6LI 6Li 6LI 6Li 6LI 6Li 6LI 6Li 6LI 6Li 6LI 6Li 6LI 6Li 
        SELECT CASE(parent%level)
        CASE(2)
            daughter(1) = E2
            daughter(2) = li6
            gamma_width = 8.2e-6
            e_kinetic = 0.0
            n_extra = 2
            ll = 2
            l_p_s = 2
        CASE (4)
    !       the distribution of kinetic energies between the
    !       protons and the 5He is well approximated by a gaussian  
            DO
                e_kinetic = 0.8875 + 0.23872*Gauss()
                IF (e_kinetic > 0.0 .AND. e_kinetic <= 1.66) EXIT
            END DO
            
            daughter(1) = proton
            daughter(2) = he5
            daughter(2)%ex = 1.667 - e_kinetic - 0.89
            ll = 1
            l_p_s = 0.5
            gamma_width = 0.54
            n_extra = 2
        CASE DEFAULT
            daughter(1) = deuteron
            daughter(2) = alpha
            IF (parent%level == 1) THEN
                gamma_width = 0.024
                e_kinetic = 2.186 +parent%ex - 1.4753
                ll = 2
                l_p_s = 3 
            ELSE IF (parent%level == 3) THEN
                gamma_width = 1.7
                e_kinetic = 4.31 + parent%ex - 1.4753
                ll = 2
                l_p_s = 2
            ELSE IF (parent%level == 5) THEN
                gamma_width = 1.5
                e_kinetic = 5.65 + parent%ex - 1.4753
                ll = 0
                l_p_s = 1
            END IF
            n_extra = 2
        END SELECT
	    CASE (7) ! 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 7Li 
        SELECT CASE (parent%level)
          CASE (1)
            daughter(1) = M1
            daughter(2) = li7
            gamma_width = 0.00626e-6
            e_kinetic = 0.0
            n_extra = 2
            ll = 1
            l_p_s = 1
        CASE DEFAULT
            daughter(1) = triton
            daughter(2) = alpha
            gamma_width = 0.093
            e_kinetic = 4.63 - 2.4681 + parent%ex
            n_extra = 2
            ll = 3
            l_p_s = 3.5
        END SELECT
	    CASE (8) ! 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 8Li 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = M1
            daughter(2) = li8
            gamma_width = 3.8e-8
            e_kinetic = 0.0
            n_extra = 2
            ll = 1
            l_p_s = 1
        CASE (2)
            daughter(1) = neutron
            daughter(2) = li7
            gamma_width = 0.033
            e_kinetic = 0.22
            n_extra = 2
            ll = 1
            l_p_s = 1.5
        END SELECT
	    END SELECT frag_A_Li
!BEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBEBE
	  CASE (4) ! Be fragments
	    frag_a_Be: SELECT CASE (parent%A)
	    CASE (6) ! 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 6Be 
        SELECT CASE (parent%level)
        CASE (1) 
            CALL RANDOM_NUMBER(prob)
            i = 1
            DO
                IF (prob <= be6_1(i)) EXIT
                i = i + 1
            END DO
            e_kinetic = REAL(i)*0.02
            daughter(1) = proton
            daughter(2) = li5
            daughter(2)%ex = 1.3706 - e_kinetic - 1.96
            gamma_width = 0.096
            ll = 1
            l_p_s = 1.5 
            n_extra = 2
        CASE (2)
            daughter(1) = proton
            daughter(2) = proton
            daughter(3) = alpha
            gamma_width = 1.16
            e_kinetic = 3.041
            n_extra = 3
        END SELECT
	    CASE (7) ! 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 7Be 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = M1
            daughter(2) = be7
            gamma_width = 2.4e-9
            e_kinetic = 0.0
            n_extra = 2
            ll = 1
            l_p_s = 1
        CASE (2)
            daughter(1) = he3
            daughter(2) = alpha
            gamma_width = 0.175
            e_kinetic = 2.9824 + parent%ex
            n_extra = 2
            ll = 3
            l_p_s = 3.5
         END SELECT
	    CASE (8) ! 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 8Be 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = alpha
            daughter(2) = alpha
            gamma_width = 6.8e-6
            e_kinetic = 0.09178 + parent%ex
            n_extra = 2
            ll = 0
            l_p_s = 0
        CASE (2) 
            daughter(1) = alpha
            daughter(2) = alpha
            gamma_width = 1.5
            e_kinetic = 3.1318 + parent%ex
            n_extra = 2
            ll = 2
            l_p_s = 2
        END SELECT
	  CASE (9)       ! 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 9Be 
        SELECT CASE (parent%level)
        CASE (1) 
            daughter(1) = E1
            daughter(2) = Be9
            gamma_width = 0.15
            e_kinetic = 0.0
            n_extra = 2
            ll = 1
            l_p_s = 1
        CASE (3:6) 
            daughter(1) = neutron
            daughter(2) = be8
            SELECT CASE (parent%level)
            CASE (3)
                gamma_width = 1.08
                e_kinetic = 1.1345 + parent%ex
                ll = 1
                l_p_s = 0.5
            CASE (4)
                gamma_width = 0.282
                e_kinetic = 1.3845 + parent%ex
                ll = 2
                l_p_s = 2.5
            CASE (5)
                gamma_width = 0.743
                e_kinetic = 4.7 - 1.6655 + parent%ex
                ll = 1
                l_p_s = 1.5
            CASE (6)
                gamma_width = 1.54
                e_kinetic = 6.76 - 1.6655 + parent%ex
                ll = 3
                l_p_s = 3.5
            END SELECT
            n_extra = 2
        CASE (2)
            daughter(1) = alpha
            daughter(2) = he5
            gamma_width = 7.7e-4
            CALL RANDOM_NUMBER(prob)
            i = 1
            DO
                IF (prob <= be9_2(i)) EXIT
                i = i + 1
            END DO
            e_kinetic = REAL(i)*0.01
            daughter(2)%ex = 0.8563 - e_kinetic - 0.89
            ll = 1
            l_p_s = 1
            n_extra = 2
        END SELECT
	  CASE (10) ! 10Be 10Be 10Be 10Be 10Be 10Be 10Be 10Be 10Be 10Be 10Be 10Be 
        SELECT CASE (parent%level)
        CASE (1)
            daughter(1) = E2
            daughter(2) = Be10
            gamma_width = 0.15
            e_kinetic = 0.0
            n_extra = 2
            ll = 2
            l_p_s = 2
        CASE (2)
            daughter(1) = E2
            daughter(2) = Be10
            gamma_width = 0.15
            e_kinetic = 0.0
            n_extra = 2
            ll = 2
            l_p_s = 2
        CASE (3)
            daughter(1) = E1
            daughter(2) = Be10
            gamma_width = 0.15
            e_kinetic = 0.0
            n_extra = 2
            ll = 1
            l_p_s = 1
        END SELECT
	  END SELECT frag_A_Be
!??????????????????????????????????????????????????
	CASE DEFAULT
	    WRITE (UNIT=*, FMT=*) "decay of unstable fragment not known"
	    WRITE (UNIT=*, FMT=*) "z=",parent%z
	    WRITE (UNIT=*, FMT=*) "A=",parent%A
	    WRITE (UNIT=*, FMT=*) "level=",parent%level
!    CALL crash ()
	    STOP
	END SELECT frag_z

	decay_time = exp_decay(gamma_width)
	DO j=1,n_extra
	    daughter(j)%id = parent%id
	    daughter(j)%time = parent%time + decay_time
	END DO

!I have made this very simple at the moment,
!we ignore recoil from gamma"s
!otherwise decay is isotropic for classical angular
!distributions.

	daughter(1)%id = parent%id
	daughter(2)%id = parent%id
	daughter(3)%id = parent%id
       
	IF (para%i_angle) THEN
	    IF (e_kinetic == 0) THEN   ! gamma ray emission
        daughter(1)%vel = parent%vel
	    ELSE
        SELECT CASE (n_extra)
        CASE (2)  ! 2body decay
            IF (para%quantum) THEN
                if(e_kinetic <= 0) THEN
                   WRITE (UNIT=*, FMT=*) "ek<0 for non_statistical",
     &	parent%Z, parent%A, parent%level, parent%ex, e_kinetic
!                   CALL crash ()
                END IF
            CALL quantum_lp(parent,daughter(1:2),ll,l_p_s,e_kinetic)
            ELSE
                CALL RANDOM_NUMBER(ran)
                theta = ACOSD(1.0-2.0*RAN)
                CALL RANDOM_NUMBER(ran)
                phi = 360.0*RAN
                ared = REAL(daughter(1)%A*daughter(2)%A)
     &                /REAL(daughter(1)%A+daughter(2)%A)
                vrel = (2.0*e_kinetic/ared)**0.5*0.9794
                v1 = vrel*ared/REAL(daughter(1)%A)
                v2 = vrel - v1
                daughter(1)%vel = vell(v1,theta,phi) + parent%vel 

                theta = 180.0-theta
                phi = MODULO(phi+180.0,360.0)

                daughter(2)%vel = vell(v2,theta,phi) + parent%vel 
                daughter(1)%ex = 0.
                daughter(2)%ex = 0.
            END IF                                        
        CASE (3)  ! 3-body decay

!           zero momentum
            p_f = 0.0

!           choose initial velocities , isotropic       
            DO i=1,3
                CALL RANDOM_NUMBER(ran)
                theta = ACOSD(1.0-2.0*RAN)
                CALL RANDOM_NUMBER(ran)
                phi = 360.0*RAN
                CALL RANDOM_NUMBER(ran)
                vv_f(i) = RAN**0.33333333/REAL(daughter(i)%A)
                v_f(i,1) = vv_f(i)*SIND(theta)*COSD(phi)
                v_f(i,2) = vv_f(i)*SIND(theta)*SIND(phi)
                v_f(i,3) = vv_f(i)*COSD(theta)
!               find total momentum     
                DO j=1,3
                    p_f(j) = p_f(j) + v_f(i,j)*REAL(daughter(i)%A)
                END DO
            END DO
!           find centre of mass velocity        
            p_f = p_f/REAL(parent%A)

!           remove centre of mass motion  and find total kinetic energy
            e = 0                   
            DO i=1,3
                vv_f(i) = 0.0
                DO j=1,3
                    v_f(i,j) = v_f(i,j) - p_f(j)
                    vv_f(i) = vv_f(i) + v_f(i,j)**2.0
                END DO
                e = e + 0.5*REAL(daughter(i)%A)*vv_f(i)**2.0*1.042
                vv_f(i) = vv_f(i)**0.5
            END DO
!           find angles and scale velocity to get desired total
!           kinetic energy      
            DO i=1,3
                theta = ACOSD(v_f(i,3)/vv_f(i))
                phi = ATAN2D(v_f(i,2),v_f(i,1))
                vv_f(i) = vv_f(i)*(e_kinetic/e)**0.5

!               add to velocity of parent       
                daughter(i)%vel = vell(vv_f(i),theta,phi) + parent%vel 
                daughter(i)%ex = 0.
           END DO
        END SELECT
	    END IF
	END IF

	RETURN 
	END SUBROUTINE non_statistical
	END MODULE nonstatistical_mod
!------------------------------end:   MODULE nonstatistical_mod-----------------------------------


!------------------------------begin:   MODULE statistical_mod-----------------------------------
	MODULE statistical_mod
	USE gem_kind
	USE functs
	USE light_stuff
	USE mass_stuff
	USE angle
	USE yrast_mod
	USE evaporation
	USE tl_stuff
	USE stat_mod
	USE nonstatistical_mod
	USE spectra
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: decay,decay_time_wise,decay_prob_wise,decay_imf,decay_lp
	PUBLIC :: decay_widths,gamma_ray,imf,complex_frag
	PUBLIC :: fission_width,fluctuate,viola,init_decay,xray,
     &	      normalise_fission

	INTEGER, PUBLIC :: I_light
	INTEGER, DIMENSION(60), PRIVATE :: Iz_max,Iz_min
	INTEGER, PRIVATE :: zm
	INTEGER, PRIVATE :: M_weight
	INTEGER, PUBLIC :: n_fission
	LOGICAL, PRIVATE :: first_weight,k_weight,X_decay

	REAL (KIND=r4), DIMENSION(60), PRIVATE :: gamma_z
	REAL (KIND=r4), PRIVATE :: shell_corr0, pair_energy0, collective
	REAL (KIND=r4), PRIVATE :: gamma_gamma, E_rot_yrast
	REAL (KIND=r4), PRIVATE :: fact,fact_old
	REAL (KIND=r4), PUBLIC  :: entropy0,U0,U_scission

	CONTAINS

	SUBROUTINE decay(parent,daughter,n_extra,ratio,weight,n_binary)
	TYPE (identity), INTENT(INOUT) :: parent
	TYPE (identity), INTENT(OUT), DIMENSION(3) :: daughter
	INTEGER, INTENT(OUT) :: n_extra
	INTEGER, INTENT(INOUT) :: n_binary
	REAL (KIND=r4), INTENT(INOUT) :: ratio
	REAL (KIND=r4), INTENT(INOUT) :: weight
	LOGICAL :: fast_fission
	REAL (KIND=r4) :: E_scission,Ek,u_thermal
	INTEGER :: iz1,iz2,ia1,ia2

	 IF (parent%level > 0) THEN             ! decay from a known level
	    CALL non_statistical(parent,daughter,n_extra)
	    RETURN
	 ELSE
	    CALL statistical(parent,daughter(1:2),n_extra,ratio,weight)

	    IF (ratio > 0.0) THEN
!      if deformed and decay is after delay time        
       IF (daughter(2)%time > para%t_delay .AND. n_extra == 2) THEN

!         ignore this decay as it happened after the delay time during
!         which the system was deformed, turn off deformation and
!         calculate decay from non-deformed object

          ratio = 0.0
          n_extra = 0

!         after delay if entropy at saddle point
!         is greater than entropy at ground state
!         we have a fission

          CALL fission_width(REAL(parent%z), REAL(parent%A),parent%Ex,
     &          parent%J,fast_fission)
        
          IF (fast_fission) THEN
  !            fast fission
               n_fission = n_fission + 1

             CALL viola(parent%z,parent%a,E_scission,ek,iz1,iz2,ia1,ia2)
             U_thermal = parent%ex - E_scission

               daughter(1)%Z = iz1
               daughter(2)%Z = iz2
               daughter(1)%A = ia1
               daughter(2)%A = ia2

               daughter(1)%time = para%t_delay
               daughter(2)%time = para%t_delay
               daughter(1)%level = 0
               daughter(2)%level = 0
       
               IF (para%quantum)THEN
                   CALL fast_fission_angle(parent,daughter,ek,U_thermal)
               ELSE
                   CALL sphere_imf(parent,daughter(1:2),U_thermal)
!                  daughter(1)%vel%v = 0.0
!                  daughter(2)%vel%v = 0.0
               END IF
               n_extra = n_extra + 2

           ELSE
!              normal statistical decay
            Call statistical (parent,daughter(1:2),n_extra,ratio,weight)
!              correct time of emission 
       daughter(1)%time = daughter(1)%time - parent%time + para%t_delay
       daughter(2)%time = daughter(1)%time
           END IF


!      if deformed and no decay possible 
       ELSE IF (n_extra == 1 .and. para%t_delay < 1000.0) THEN
           !another try at decaying after the
           !system relaxes its shape
            ratio = 0.0
            CALL statistical (parent,daughter(1:2),n_extra,ratio,weight)
            IF (parent%level == 1) THEN
      daughter(1)%time  = daughter(1)%time - parent%time  + para%t_delay
       daughter(2)%time  = daughter(1)%time
            END IF 
        END IF
       END IF

       if (n_extra > 1) THEN
        IF (daughter(1)%a > 4 .OR. daughter(1)%level >= 1) THEN
             daughter(1)%id = ISHFT(parent%id,1)
             daughter(2)%id = IBSET(daughter(1)%id,0)
             n_binary = n_binary + 1
        ELSE
             daughter(1)%id = parent%id
             daughter(2)%id = parent%id
        END IF
       END IF

       IF (X_decay) CALL xray(parent,daughter,n_extra)

       END IF
	 END SUBROUTINE decay


	SUBROUTINE statistical(parent,daughter,n_extra,ratio,weight)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine chooses the decay channel of a nucleus taking into
! account the decay widths and the weighting factors. it then calculates
! the weight of the event. if kk=0, then no complex fragment decay was
! possible. the decaying nucleus is parent%Z, parent%A  weighting is
! considered in choosing the exit channels if k_weight= .TRUE.
!
! DUMMY ARGUMENTS:
!
!
! IMPLICIT INPUTS:
!
! I_channels - total number of exit channels open
! I_light - number of light particle evaporative channels open
! K_weight - flag to signal weighting is included is choosing exit channel
! gamma(1:I_channel) - consecutive tally of total decay widths at each
!                      channel number
!                    gamma(I_light) - total decay width for light particle
!            emission
!            gamma(I_channel) -total decay width for complex fragment
!            emission
!  Za(1:I_channel) - proton number of lightest fragment for each exit channel
!  Aa(1:I_channel) - nucleon number of lightest fragment for each exit
!                    channel
!  fact - factor by which decay widths are increased when weighting in on
!  weight - running weight of event
!  gamma_z(3:Zm) - decay width as a function of fragment Z_value
!  Zm - number of Z-values
!
!                           channel
! weight - new running weight of event
! fact - new weighting factor
!
!
!
	TYPE (identity), INTENT(INOUT) :: parent
	TYPE (identity), INTENT(OUT), DIMENSION(2) :: daughter
	INTEGER, INTENT(OUT) :: n_extra
	REAL (KIND=r4), INTENT(IN) :: ratio
	REAL (KIND=r4), INTENT(INOUT) :: weight

	REAL (KIND=r4) :: gamma_total
	n_extra = 0
	gamma_total = 0.0
	gamma_gamma = 0.0


! calculate the partial decay widths for all possible channels
	CALL decay_widths(parent,ratio,gamma_total)

C 	PRINT*, I_LIGHT, I_CHANNELS, 'AA'	!HG CHENG

! if the gamma ray decay width represents 99% of the total decay width
! then we assume the system will gamma decay to the ground state
	IF (gamma_total > 0.0 .AND. gamma_gamma >= 0.99*gamma_total)
     &	gamma_total = 0.0

	IF (gamma_total <= 0.0)   THEN     !no decay of the compound nucleus possible
	   n_extra = 1
	   daughter(1) = parent
	   daughter(1)%exx = daughter(1)%ex
	   daughter(1)%ex = 0.0 
	   RETURN
	END IF 

	n_extra = 2

	IF (para%time_wise) THEN
	   CALL decay_time_wise (parent,daughter(1:2))
	ELSE
	   CAll decay_prob_wise (parent,daughter(1:2),weight)
	END IF

	RETURN
	END SUBROUTINE statistical

	   SUBROUTINE INIT_decay ()
!  initilization of flag and parameters for statistical
	   list_hf = list_hf0
	   first_weight = .TRUE.
	   IF (para%N_weight > 0) K_weight = .TRUE.
	   m_weight = para%n_weight
	   x_decay = .TRUE.
	   fact = 1.0
	   fact_old = 1.0
	   n_fission = 0
	   RETURN
	   END SUBROUTINE INIT_decay
 
	SUBROUTINE decay_time_wise(parent,daughter)

![subroutine_header_comments]
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(OUT),DIMENSION(2) :: daughter
	INTEGER :: i,minimum_channel
	REAL (KIND=r4) :: decay_time
	REAL (KIND=r4) :: minimum_time,eta,t_star,t_statistical
	REAL (KIND=r4) :: t_delay_eta

	minimum_time = 1.0e30
	minimum_channel = 0
! light charged particles and gamma rays - no dynamical delays      
	DO i=1,I_light

	   time(i) = exp_decay(gamma(i))

	   IF (time(i) <= minimum_time) THEN
       minimum_Channel = i
       minimum_time = time(I)
	   END IF
	END DO

! complex fragments and fission - dynamical delays
	DO i=I_light+1,I_Channels
	   eta = REAL(parent%A-2*Aa(i))/REAL(parent%A)
	   t_delay_eta = para%t_delay*EXP(-eta**2.0/
     &				 (2.0*para%sig_delay**2.0))

	   t_statistical = exp_decay(gamma(i))

	   IF (parent%time > t_delay_eta) THEN  !elapsed time already
                                            !greater than fission
                                            !delay time, so ignore it.
	       time(i) = t_statistical
	   ELSE
       IF (para%sharp_delay) THEN !fission decay width is zero up to
                             !t_delay_eta
           time(i) = t_statistical + t_delay_eta - parent%time 
       ELSE !fission decay widths ramps up until  t_delay_eta

           t_star = (exp_decay(gamma(i))*t_delay_eta*2.0 + 
     &          parent%time**2.0)**0.5 - parent%time
           IF (t_star <= t_delay_eta) THEN ! emission during
                                             ! ramping period
               time(i) = t_star   
           ELSE                              !emission after
                                             !ramping period
               time(i) = t_statistical + t_delay_eta - parent%time
           END IF
       END IF
      END IF

       IF (time(i) <= minimum_time) THEN
       minimum_Channel = i
       minimum_time = time(I)
       END IF
	END DO

	i = minimum_channel
	decay_time = parent%time + minimum_time

C 	PRINT*, I, I_LIGHT, 'DEC'	!HG CHENG

	IF (i <= I_light) THEN
!   we have a light particle evaporation!
!   put the parameters of the exit channel in the array of the
!   decay chain. The decay channel is i

       CALL decay_lp(parent,daughter,i,decay_time)
	ELSE
!   we have a complex fragment emission.
!   put the parameters of the exit channel in the array of the
!   decay chain the decay channel is i

       CALL decay_imf (parent,daughter,i,decay_time)
	END IF

	RETURN
	END SUBROUTINE decay_time_wise
        
	SUBROUTINE decay_prob_wise(parent,daughter,weight)

![subroutine_header_comments]

	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(OUT),DIMENSION(2) :: daughter

	INTEGER :: NZ,ii,ij,iz
	INTEGER :: I,I_max,I_min
	REAL (KIND=r4) :: ran
	REAL (KIND=r4) :: gamma_tot_weight,gamma_tot_real,gamma0
	REAL (KIND=r4) :: prob_l,prob1,prob2,prob3,prob_tot,x,y
	REAL (KIND=r4) :: tau,tau1,tau2,delta_time,fract
	REAL (KIND=r4) :: decay_time,weight
	REAL (KIND=r4), PARAMETER :: fact0 = 0.3

!create running sum of partial decay widths
	IF (I_light > 1) THEN
	   DO i = 2,I_light
         gamma(I) = gamma(i) + gamma(i-1)
	   END DO
	END IF

	IF (I_channels > I_light+1) THEN
	    DO i = I_light+2,I_channels
          gamma(I) = gamma(I) + gamma(I-1)
	   END DO
	END IF

!check to see if weighting is possible
	IF (I_channels == I_light .OR. I_light == 0) THEN
	   IF (para%n_weight == 1) THEN
	      K_weight = .FALSE.
	      fact_old = 1.0
	   END IF
	   fact = 1.0
	ELSE
	   fact = fact_old
	END IF

!determine the weighting factor at start of decay
!(the time delay is not considered here)
	IF (K_weight .AND. first_weight) THEN
	   fact = fact0*gamma(I_light)/(1.0-fact0)/gamma(I_channels)
	   IF (fact <= 1.0) fact = 1.0
	   fact_old = fact
	   first_weight = .FALSE.
	END IF

!calculate total decay width
	gamma_tot_weight = 0.0
	gamma_tot_real = 0.0

	IF (I_light > 0) THEN
	   gamma_tot_weight = gamma_tot_weight + gamma(I_light)
	   gamma_tot_real = gamma_tot_real + gamma(I_light)
	END IF

	IF (I_channels > I_light) THEN
	   gamma_tot_weight = gamma_tot_weight + fact*gamma(I_channels)
	   gamma_tot_real = gamma_tot_real + gamma(I_channels)
	END IF

!calculate the probability for light particle decay(including weighting and time delay)
	IF (I_light == 0) THEN
	   prob_l = 0.0
	ELSE IF (I_channels > I_light) THEN
	   IF (parent%time > para%t_delay .OR. para%t_delay == 0) THEN
       prob_l = gamma(I_light)/gamma_tot_weight
	   ELSE
       tau = (para%t_delay - parent%time)*gamma(I_light)/0.65824
       prob_l = (gamma_tot_weight - fact*gamma(I_channels)*EXP(-tau))/
     &	gamma_tot_weight
	   END IF
	ELSE
	   prob_l = 1.0
	END IF

! generate random numbers
	CALL RANDOM_NUMBER(prob1)
	CALL RANDOM_NUMBER(prob2)
	I = 0

!choose the exist channel
	IF (prob1 <= prob_l) THEN
!  we have a light particle evaporation!
!  Find which decay channel the random number generator has selected

	   DO 
       I = I + 1
       fract = gamma(i)/gamma(I_light)
       if (fract > prob2) EXIT
	   END DO

!  if weighting turned on, calculate the running weight
      IF (K_weight) THEN
       IF (parent%time > para%t_delay) THEN
           weight = weight*gamma_tot_weight/gamma_tot_real
       ELSE
           weight = weight/gamma_tot_real/prob_l*(gamma(I_light)+
     &			  gamma(I_channels)*(1.0-EXP(-tau)))
       END IF
      END IF

!  calculate the time since the formation of the compound nucleus at
!  which this light particle emission occurred
	   IF (parent%time > para%t_delay) THEN
       decay_time = parent%time + exp_decay(gamma_tot_real)
	   ELSE
       tau1 = 0.65824/gamma(I_light)
       tau2 = 0.65824/gamma_tot_real
       delta_time = para%t_delay - parent%time
       prob1 = tau1*(1.0-EXP(-delta_time/tau1))
       prob2 = tau2*EXP(-delta_time/tau2)
       prob_tot = prob1 + prob2
       CALL RANDOM_NUMBER(ran)
       IF (RAN > prob1/prob_tot) THEN
           CALL RANDOM_NUMBER(ran)
           decay_time = para%t_delay - tau2*LOG(RAN+1.0e-37)
       ELSE
           x = EXP(-delta_time/tau1)
           CALL RANDOM_NUMBER(ran)
           y = RAN*(1.0-x)+x
           decay_time = parent%time - tau1*LOG(y)
       END IF
       END IF

!  put the parameters of the exit channel in the array of the
!  decay chain. The decay channel is i
       CALL decay_lp(parent,daughter,i,decay_time)

	ELSE
!  we have a complex fragment emission.
!  choose the exit channel.
!  With weighting first choose the charge split

	   IF (K_weight) THEN
      IF (parent%time > para%t_delay) THEN
          weight = weight*gamma_tot_weight/fact/gamma_tot_real
      ELSE
          weight = weight*gamma(I_channels)*EXP(-tau)
     &             /gamma_tot_real/(1.0-prob_l)
      END IF
	   END IF

	   IF (K_weight .and. para%N_weight == 1) THEN
!      count how many charge splits are possible
       NZ = 0
       DO I = 3,Zm
          IF (gamma_z(I) > 0.0) NZ = NZ + 1
       END DO

       II = INT(prob2*REAL(NZ)) + 1
       IF (II > NZ) II = NZ

!      find which channel was selected
       IJ = 0
       IZ = para%Z_imf_min - 1

       DO
          IZ = IZ + 1
          IF (gamma_z(iz) > 0) IJ = IJ + 1
          IF (IJ >= II) EXIT
       END DO

!      Next choose the mass split.
!      Find channels which correspond to chosen charge split.
       I_max = Iz_max(iz)
       I_min = Iz_min(iz)

       IF (I_min == I_light+1) THEN
           gamma0 = 0.0
       ELSE
           gamma0 = gamma(I_min-1)
       END IF

!      channels i_min+1 to i_max are associated with our chosen
!      charge split
       CALL RANDOM_NUMBER(prob3)
       I = I_min - 1
       DO
           IF (fract > prob3) EXIT
           I = I + 1
           IF (I >= I_max) EXIT
           fract = (gamma(i)-gamma0)/gamma_z(iz)
       END DO
       weight = weight*REAL(NZ)*gamma_Z(iz)/gamma(I_channels)
       K_weight = .FALSE. ! after complex fragment emission weighting
       fact = 1.0          ! turned off
       fact_old = 1.0
	   ELSE
!      without weighting choose exist channel according to the calculated decay widths
       IF (M_weight > 0) THEN
           M_weight = M_weight - 1
           IF (M_weight == 0) THEN
              K_weight = .FALSE.
              fact = 1.0
              fact_old = 1.0
           END IF
       END IF

       I = I_light
       DO
          I = I + 1
          fract = gamma(i)/gamma(I_channels)
          IF (fract > prob2) EXIT
       END DO
	   END IF

!  calculate the time since the formation of the compound nucleus for the saddle crossing
	   IF (parent%time > para%t_delay) THEN
        decay_time = parent%time + exp_decay(gamma_tot_real)
	   ELSE
        decay_time = para%t_delay + exp_decay(gamma_tot_real)
	   END IF

!  put the parameters of the exit channel in the array of the decay chain the decay channel is i
	   CALL decay_imf (parent,daughter,i,decay_time)

	END IF
	RETURN

	END SUBROUTINE decay_prob_wise


	SUBROUTINE decay_imf (parent,daughter,i,decay_time)
! an imf emission or fission has been decided this subroutine loads the information about
! the decay products to the frag array
	INTEGER, INTENT(IN) :: i
	REAL (KIND=r4), INTENT(IN) :: decay_time
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter

	LOGICAL, SAVE :: first = .TRUE.

	n_fission = n_fission + 1
	daughter(1)%Z = Za(i)
	daughter(1)%A = Aa(i)
	daughter(2)%Z = parent%Z - daughter(1)%Z
	daughter(2)%A = parent%A - daughter(1)%A
	daughter(1)%level = 0
	daughter(2)%level = 0

C 	  if(u_therm(i) .LE. 0.1)  then      ! Zhaoqing Feng
C 	  daughter(1)   = parent
C 	  daughter(1)   = parent
C 	  daughter(2)%Z = parent%Z - daughter(1)%Z
C 	  daughter(2)%A = parent%A - daughter(1)%A
C 	  daughter(1)%level = 0
C 	  daughter(2)%level = 0
C 	  daughter(2)%time = 0.
C 	  daughter(1)%time = 0.
C  	  return
C 	  endif

	IF (para%quantum)THEN
	    IF (para%imf_option == 1) THEN
       CALL quantum_imf(parent,daughter,u_therm(i),dau(i)%l_p_s,
     &				  U_scission)
	    ELSE
        IF (first) THEN
          WRITE (UNIT=*, FMT=*) "no quantum mechanical treatment for
     &		 IMF spin available yet"
       WRITE (UNIT=*, FMT=*) "set quantum=.FALSE. for semi-classical"
          first = .FALSE.
        END IF
        daughter(2)%vel%v = 0.0
        daughter(1)%vel%v = 0.0
        daughter(2)%vel%theta = 0.0
        daughter(1)%vel%theta = 0.0
        daughter(2)%vel%phi = 0.0
        daughter(1)%vel%phi = 0.0
        daughter(1)%exx = 0.0
        daughter(2)%exx = 0.0
        daughter(2)%Ex = 0.0
        daughter(1)%Ex = 0.0
        daughter(2)%level = 0
        daughter(1)%level = 0
        daughter(2)%J= 0.0
        daughter(1)%J= 0.0
	   END IF
	ELSE

       CALL sphere_imf(parent,daughter,u_therm(i))
!    daughter(2)%vel%v = 0.0
!    daughter(1)%vel%v = 0.0
	END IF

	daughter(2)%time = decay_time
	daughter(1)%time = decay_time

	RETURN
	END SUBROUTINE decay_imf

	SUBROUTINE decay_lp (parent,daughter,i,decay_time)
!	 light_particle emission or gamma-ray emission has been decided
!	 This subroutine loads the information about the decay products to the frag array
c
	INTEGER, INTENT(IN) :: i
	REAL (KIND=r4), INTENT(IN) :: decay_time
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT),DIMENSION(2) :: daughter

	daughter(1)%Z = Za(i)
	daughter(1)%A = Aa(i)
	if (za(i) >= 0) THEN                       ! because za=-1 for gamma rays
	   daughter(2)%Z = parent%Z - Za(i)
	   daughter(2)%A = parent%A - Aa(i)
	   daughter(1)%exx = 0.0
	ELSE
	   daughter(2)%Z = parent%Z 
	   daughter(2)%A = parent%A 
	   daughter(1)%exx = dau(i)%excit
	END IF    
 
	daughter(1)%level = dau(i)%level
	daughter(1)%J = dau(i)%s_particle
	daughter(1)%k = dau(i)%k_particle
	daughter(1)%ex = dau(i)%excit
	daughter(2)%Ex = dau(i)%Ex_residue
	daughter(2)%exx = 0.0
	daughter(2)%J = dau(i)%J_residue
	daughter(2)%k = dau(i)%k_residue
	daughter(1)%time = decay_time
	daughter(2)%time = decay_time
	daughter(2)%level = 0

	IF (para%I_angle) THEN
	   IF (para%quantum) THEN         !quantum approach
          CALL quantum_lp(parent,daughter,
     &         dau(i)%ll,dau(i)%l_p_s,dau(i)%e_kinetic)
	   ELSE             !semi classical approach
       IF (daughter(1)%z < 0) THEN          ! gamma emission
           daughter(2)%vel = parent%vel     !velocity vector
           daughter(2)%spin%theta = parent%spin%theta
           daughter(2)%spin%phi = parent%spin%phi
           daughter(1)%vel%v = 0.0
           daughter(1)%vel%theta = 0.0
           daughter(1)%vel%phi = 0.0
       ELSE
           CALL sphere_lp(parent, daughter, dau(i)%ll, dau(i)%l_p_s,
     &					dau(i)%e_kinetic)
       END IF
	   END IF
	END IF
	END SUBROUTINE decay_lp

	SUBROUTINE decay_widths(parent,rratio,gamma_total)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine calculates all the decay widths for both complex
! fragment decay and light particle emission.  decay widths are
! found in array gamma.  gamma(i_light) is the total decay width
! for light particle emission and gamma(i_channels) is the total
! decay width for complex fragment emission. if i_channels=i_light there is no complex fragment decay channels.
! the variable kk is set to zero if no decay is possible the
! calculated decay width are in units of MeV.
!
! DUMMY ARGUMENTS:
!
! input: k (INTEGER ::) - channel number in Z,A,I_decay,Ex,J arrays of nucleus
!                       whose decay is being simulated.
!
! IMPLICIT INPUTS:
!
! Z(k) - proton number of decaying nucleus
! A(k) - nucleon number of decaying nucleus
! Ex(k) - excitation energy of decaying nucleus
! J(k) - spin of decaying nucleus
!
! IMPLICIT OUTPUTS:
!
! I_channels - total number of exit channels open
! I_light - number of light particle evaporative channels open
! gamma(1:I_channel) - consecutive tally of total decay widths at each
!                      channel number
!                      gamma(I_light) - total decay width for light particle
!            emission
!            gamma(I_channel) -total decay width for complex fragment
!            emission
!  Za(1:I_channel) - proton number of lightest fragment for each exit channel
!  Aa(1:I_channel) - nucleon number of lightest fragment for each exit channel
!  gamma_z(3:Zm)   - decay width as a function of fragment Z_value
!
!
! SIDE EFFECTS:
!
!    none
!
!
	USE evap
	TYPE (identity), INTENT(INOUT) :: parent
	REAL (KIND=r4), INTENT(IN) :: rratio
	REAL (KIND=r4), INTENT(out) :: gamma_total            !total decay width
	INTEGER :: mode
	INTEGER :: i,Iz0,Ia0,I_evap
	REAL (KIND=r4) :: A0,Z0,J0,k0
	REAL (KIND=r4) :: Ex0,aden,temp,mass0ld
	REAL (KIND=r4) :: prob,expand
	REAL (KIND=r4) :: level_density,gamma_tot
	REAL (KIND=r4) :: es0,ec0,x,homega,ratio

	Z0 = parent%Z
	A0 = parent%A
	Ex0 = parent%Ex
	J0 = parent%J
	K0 = parent%k
	gamma_total = 0.0

!calculate the mass defect of the compound nucleus
	iz0 = Z0
	ia0 = A0

	CALL get_pairing(NINT(Z0),NINT(A0),pair_energy0)

!get experimental mass and finite range mass
	CALL getmass(iZ0,iA0,Mass0,mass0ld) 
	if (iz0 > Z_shell) THEN
	   shell_corr0 = mass0 - mass0ld - pair_energy0
	ELSE
	   shell_corr0 = 0
	END IF 

	ratio = rratio

	IF (para%fluct) THEN
	     CALL fluctuate(Z0,A0,Ex0,J0,k0,ratio)
	     parent%k = k0
	END IF 

	IF (ratio > 0.0) CALL def_rot(ratio,bs)

	CALL yrast(Z0,A0,J0,K0,ratio,E_rot_yrast)   ! find yrast energy


	U0 = Ex0 - E_rot_yrast                      ! thermal energy above the yrast line

	IF (para%k_sum) THEN
	  expand = 1.0
	  CALL prepare_deformed_emission (ratio,expand)
	ELSE IF (para%expand) THEN
	   expand = 1.1
	   if (A0 > 170) expand=1.15
!  expand = 1.0 + 0.2*u0/A0 
	  CALL prepare_deformed_emission (ratio,expand)
	END IF 

	IF (U0 <= 0.0) THEN         ! this is not possible
	   RETURN
	END IF


	 CALL get_little_a(A0,Z0,j0,U0,bs,shell_corr0,E_rot_yrast,
     &    	aden,entropy0,pair_energy0,1.0)


	IF (U0 <= 0.0) THEN         ! this is not possible
	   RETURN
	END IF


	Es0 = 17.9439*(1.0-1.7826*((A0-2.0*Z0)/A0)**2.0)*A0**0.66666666
	Ec0 = 0.7053*Z0**2.0/A0**0.33333333
	x = ec0/2.0/es0
	if(x>1.) x=1.    !suj
	homega = 36.34*((1.0-x)/A0)**0.5
	temp = (U0/aden)**0.5
	collective = homega/temp
	
!prepare pe and ee so we can use get_level_density
	pe = 1.0
	ee = entropy0

	M_inertia_0 = 0.4*r0**2.0*A0**1.6666666

	  CALL get_level_density(U0,J0,aden,M_inertia_0,entropy0,
     &  level_density)
	  CALL collective_enhancement(iz0,ia0,U0,aden,m_inertia_0,
     &  level_density)

!final command for level density preparation
	  pe = level_density*6.2831

	I_channels = 0
	if (k0 >= 0) THEN
!calculate the partial decay widths for light particle emission
	   IF (para%tl_iwbc .AND. z0 > 10) THEN  ! transmission coeff from
                                            ! incoming wave
                                            !  boundary condition model
   !  the fit to the IWBC transmission is not very good for emssion
   !  from very light system, so for z0<10 we do not use it

	   mode = 1
   !CALL hist_zero() !histgrams for n,p and alpha 
	   DO
       IF (list_HF(mode) <= 0) EXIT
       CALL Hauser_Feshbach (list_HF(mode),Z0,A0,Ex0,J0,K0,ratio)
       mode = mode + 1
	   END DO

! *************** weight channels to get more statistics ********
! the multiplicities will be wrong, but the shape of the spectra
! should be approximately right
!         DO i=1,i_channels
!           if (evap_mode(i) == 34 .OR. evap_mode(i) == 35)THEN
!               gamma(i) = gamma(i)*10.0
!           ELSE if (evap_mode(i) == 13 .OR. evap_mode(i) == 15)THEN
!               gamma(i) = gamma(i)*4.0
!           ELSE if (evap_mode(i) > 18) THEN
!               gamma(i) = gamma(i)*2.0
!           END IF 
!         END DO 
!*************************************************************

        if (para%threshold > 0.0 .AND. n_fission == 0) THEN
            gamma_tot = SUM(gamma(1:I_channels))
            mode = 1
            DO i=1,i_channels

               if (gamma(i)/gamma_tot > para%threshold) THEN
                   list_hf(mode) = evap_mode(i)
                   mode = mode  + 1
c      if be8 gound state is kept - then keep the first excited state
                ELSE IF (evap_mode(i) == 38) THEN  
                   IF(list_hf(mode-1) == 37) THEN 
                     list_hf(mode) = evap_mode(i) 
                     mode = mode  + 1
                   END IF 
               END IF 
            END DO
            list_hf(mode) = -1
        END IF 
       
       ELSE			! sharp cut-off transmission coeff.
       j0 = INT(j0)    ! evaporation simple assumes integer spin
       parent%j = j0
       mode = 1
      DO
          IF (list_HF(mode) <= 0) EXIT
          CALL evaporation_simple(mode,Z0,A0,Ex0,J0)
          mode = mode + 1
      END DO
	    END IF
	 END IF 

!WRITE (UNIT=*, FMT=*) gamma(1)+gamma(2)+gamma(3)+gamma(4)+gamma(5)+gamma(6)
!call write_array
!CALL WRITE_spectra2
!WRITE (UNIT=*, FMT=*) gamma(1)/SUM(gamma(1:3)),gamma(2)/SUM(gamma(1:3))&
!,gamma(3)/SUM(gamma(1:3))
!WRITE (UNIT=*, FMT=*) gamma(1),gamma(2),gamma(3)
!STOP

	if (para%k_sum) ratio = 0.0

	i_evap = I_Channels

! calculate the decay widths for gamma ray emission if near the yrast line
	DO mode=1,2
	   CALL gamma_ray(mode,Z0,A0,Ex0,J0,ratio)
	END DO


	I_light = I_channels
	IF (rratio <= 0.0) THEN   ! no IMF emission from deformed compound nucleus during fission delay

!  Complex fragment decay. or fission
	   SELECT CASE(para%imf_option)
	   CASE (1)
       CALL fission_width(Z0,A0,Ex0,J0)
	   CASE (2) 
!      Calculate barriers for imf emission, then decay widths
       IF (A0 > A_cut) THEN
           CALL sadfit(Z0,A0,J0)
           CALL IMF(Z0,A0,Ex0,J0)
       END IF
	   CASE(3)
!      Calculate barriers for imf emission, then decay widths
       IF (A0 > A_cut) THEN
           CALL sadfit(Z0,A0,J0)
           CALL IMF(Z0,A0,Ex0,J0)
           CALL normalise_fission(Z0,A0,Ex0,J0)
       END IF

	  END SELECT

	END IF


	IF (I_channels > Max_channels) THEN
	   WRITE (UNIT=*, FMT=*) "dimension of channels exceeded"
	   STOP
	END IF

! if no decay modes are possible or only gamma decay is possible
! then turn the decay off
	 IF (I_channels <= 0 .OR. I_channels == 1 .AND. Aa(1) == 0) THEN
	   gamma_tot = 0.0
	   RETURN
	 END IF

	if (para%tl_iwbc) THEN
	    gamma_tot = SUM(gamma(1:I_channels))
	    DO i=1,I_evap
        prob = gamma(i)/gamma_tot
        mult1(evap_mode(i)) = mult1(evap_mode(i)) + prob
	    END DO
	END IF

	gamma_total = 0.0 
	IF (I_channels > 0) gamma_total = SUM(gamma(1:I_channels))


!histograms for n,p, alpha's
!CALL hist_normalize(gamma_total*para%r_num,int(j0))

	IF (gamma_total <= 0.0) THEN
   ! this step is not strictly correct as to get to this part of the
   ! program I_channels must be greater than zero and hence at one of
   ! the possible exit-channels is open. However because of probably
   ! some underflow error, the total decay width was calculated to be
   ! zero. I stop the decay of the fragment here, but in reality the
   ! fragment should be allowed to decay.

	   RETURN
	END IF

	RETURN
	END SUBROUTINE decay_widths
 
	SUBROUTINE gamma_ray(mode,Z0,A0,Ex0,J0,ratio)
!+
!
! FUNCTIONAL DESCRIPTION:
! calculates the decay width for gamma emission
!
!-
	INTEGER, INTENT(IN) :: mode
	REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0,ratio
	INTEGER :: jj,i
	REAL (KIND=r4) :: U2_max,U2,de,aden,sum_gam_max,e,E_rot
	REAL (KIND=r4) :: s0,s0_min,s0_max
	REAL (KIND=r4) :: sum_gam,gam_max,ef
	REAL (KIND=r4) :: const,fract,prob
	REAL (KIND=r4) :: entropy,e_GDR,gamma_GDR,Fe1
	REAL (KIND=r4) :: level_density,ran
	LOGICAL, PARAMETER :: GDR = .FALSE.

	REAL (KIND=r4), PARAMETER :: wue(2)= (/6.8e-8,4.9e-14/)
	REAL (KIND=r4), SAVE :: inhibition(2) = (/1.0e-3,50.0/)

	inhibition(1) = para%E1_strength
	inhibition(2) = para%E2_strength
	de = 1
	if (u0 <= 15.0) de = 0.2

!giant dipole strength function
!see D.R. Chakrabarty et al. Phys Rev. C36 (1987) 1886
!for details
	E_GDR = 77.0/A0**0.33333333
	gamma_GDR = 4.8+0.0026*U0**1.6
	Fe1 = 2.09e-5*Z0*(A0-Z0)/A0*gamma_GDR

!const = A0**(REAL(2*mode)/3.0)*wue(mode)
	const = A0**(REAL(2*mode)/3.0)*wue(mode)*3.14159*2.0
	const = const*inhibition(mode)
	s0 = J0
	s0_min = ABS(s0-REAL(mode))
	s0_max = s0 + REAL(mode)
	s0 = s0_min
	sum_gam_max = 0.0
	jj = 0
	loop_s0_sum: DO 
	  CALL yrast(Z0,A0,s0,0.0,ratio,E_rot)
	  U2_max = Ex0 - E_rot
	  IF (U2_max > 0.01) THEN
	      e = 0.1
	      sum_gam = 0.0
	      gam_max = 0.0
	      loop_e_sum: DO 
          U2 = U2_max - e
          IF (U2 > 0.01) THEN
              CALL get_little_a(A0,Z0,s0,U2,bs,shell_corr0,E_rot,
     &             aden,entropy,pair_energy0,1.0)
              
              ef = entropy -  ee
              IF (ef <= 65.0 .AND. U2 > 0.0) THEN
                  jj = jj + 1
                  IF (jj > n_light) THEN
             PRINT*,  "increase N_light in stat_mod from ",n_light
                        STOP 
                  END IF

             CALL get_level_density(U2,s0,aden,M_inertia_0,entropy,
     &			level_density)

            storage(jj)%gam = e**(2*mode+1)*const*de*level_density

!                 if GDR strength function is used
!                 do it here
                  if (mode == 1 .AND. GDR)THEN
           storage(jj)%gam = storage(jj)%gam/const*
     &     Fe1*e/((e**2.0-E_GDR**2.0)**2.0+gamma_GDR**2.0*e**2.0)
                  END IF


                  storage(jj)%spin = NINT(s0*2.0)
                  storage(jj)%energy = e
                  gam_max = MAX(storage(jj)%gam,gam_max)
                  sum_gam = sum_gam + storage(jj)%gam
                  IF (jj > 1) storage(jj)%gam = 
     &            storage(jj-1)%gam + storage(jj)%gam
                  IF (storage(jj)%gam < 0.01*gam_max) EXIT loop_e_sum
              ELSE
                  EXIT loop_e_sum
              END IF
          ELSE
              EXIT loop_e_sum
          END IF
          e = e + de
      END DO loop_e_sum
      sum_gam_max = MAX(sum_gam,sum_gam_max)
      IF (sum_gam < 0.01*sum_gam_max) EXIT loop_s0_sum
	 ELSE
      EXIT loop_s0_sum
	  END IF
	  s0 = s0 + 1
	  IF (s0 > s0_max) EXIT loop_s0_sum
	END DO loop_s0_sum
	IF (jj == 0) RETURN
	IF (storage(jj)%gam <= 0) RETURN
	  I_channels = I_channels + 1
	IF (I_channels > n_light_ch) THEN
	    WRITE (*,*) " increase n_light_ch in stat_mod.f90",i_channels
	    STOP
	END IF 
	Za(I_channels) = -1
	Aa(I_channels) = mode
	gamma(I_channels) = storage(jj)%gam
	gamma_gamma = gamma_gamma + storage(jj)%gam
	i = 0
	CALL RANDOM_NUMBER(prob)
	DO
	   i = i + 1
	   fract = storage(i)%gam/storage(jj)%gam
	   IF (fract > prob) EXIT
	END DO
	dau(I_channels)%J_residue = REAL(storage(i)%spin)/2.0
	dau(I_channels)%s_particle = REAL(mode)
	CALL RANDOM_NUMBER(ran)
	e = storage(i)%energy + de*(0.5 - RAN)
	IF (e < 0.0) e = 0.0
	dau(I_channels)%Ex_residue = Ex0 - e
	dau(I_channels)%level = 0
	dau(I_channels)%ll = mode
	dau(I_channels)%l_p_s = mode
	dau(I_channels)%excit = e
	RETURN
	END SUBROUTINE gamma_ray
 

	SUBROUTINE IMF(z0,a0,Ex0,J0)
!
! FUNCTIONAL DESCRIPTION:
!
! Calculated the decay-widths for complex fragment emission Z>2 from the
! transition state formalism of Moretto.
!
! DUMMY ARGUMENTS:
!
! input: Z0 (REAL (KIND=r4) ::) proton number of decaying nucleus
!        A0 (REAL (KIND=r4) ::) nucleon number of decaying nucleus
!        Ex0 (REAL (KIND=r4) ::) excitation energy of decaying nucleus
!        J0 (REAL (KIND=r4) ::) spin of decaying nucleus
!
! IMPLICIT INPUTS:
!
! arrays cha1,cha2 - maximum and minimum mass for each Z-values between which
!                       the emitted complex fragment must lie.
!
! IMPLICIT OUTPUTS:
!
! Zm - number of distinct exit channels classified by Z-value.
! gamma_z(3:Zm) - decay widths as a function of the smaller complex fragment Z-value.
! I_channels - total number of open decay channels.
! gamma(1:I_channels) - decay width for each exit channel
!
! SIDE EFFECTS:
!
!    none
!
!

	REAL (KIND=r4), INTENT(IN) :: Z0, A0, Ex0, J0

	INTEGER :: ia0,iz0,Iz1,Iz2,Ia_max,Ia_min,Ia
	INTEGER :: kl
	REAL (KIND=r4) :: Z1,Z2,A1,A2

	iz0 = NINT(z0)
	ia0 = NINT(A0)
	zm = INT(z0/2.0)
	iz1 = para%z_imf_min - 1
	DO
	  IF (iz1 >= zm) EXIT
	  iz1 = iz1 + 1
	  iz2 = iz0 - iz1
	  z1 = REAL(iz1)
	  z2 = REAL(iz2)
	  iA_max = MIN(cha2(iz1),ia0 - cha1(iz2))
	  iA_min = MAX(cha1(iz1),ia0 - cha2(iz2))


! if z0 is even there will be some double counting
! unless a1 is restricted less then or equal to "a0/2".
	  IF (iz1 == iz2) iA_max  = MIN(iA_max,INT(A0/2.0)) 
	  gamma_z(iz1) = 0.0
	  DO ia = iA_min,iA_max
	     a1 = REAL(ia)
	     a2 = REAL(ia0 - ia)
!    calculate the decay width for the binary division
	     CALL complex_frag(A0,Z1,A1,Z2,A2,Ex0,J0,kl)
!    WRITE(8,*)z1,a1,z2,a2,gamma(I_channels)
	     IF (kl == 1) THEN
         IF (gamma_z(iz1) == 0.0) Iz_min(iz1) = I_channels
         gamma_z(iz1) = gamma_z(iz1) + gamma(I_channels)
         Iz_max(iz1) = I_channels
	   END IF
	  END DO
	END DO
	RETURN
	END SUBROUTINE IMF


!this subroutine checks to see if there is a fission peak (maximum) in
!the charge distribution, if so then its scales all channels in the 
!peak so that the total fission yield is given by the one-dimensional 
!Bohr-wheeler transitional state formalism 
	 SUBROUTINE Normalise_fission(Z0,A0,Ex0,J0)
	REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0
	INTEGER :: iz1,Zmin,i,iz
	REAL (KIND=r4):: gamma_min
	REAL (KIND=r4):: gam_tot,gam_sym,factor

!find the minimum in the z distribution
	gamma_min = 9.e32
	iz1 = para%z_imf_min
	DO iz = iz1,zm-1
       if (gamma_z(iz) < gamma_min) THEN
       gamma_min = gamma_z(iz)
       zmin = iz
       END IF 
	END DO

	if (zmin == zm-1) return                  !no fission peak, so do not modify

!calculate total width from minimum to symmetry
	gam_tot = 0.
	DO iz = zmin,zm
	   gam_tot = gam_tot + gamma_z(iz)
	END DO

	i = i_channels
	CALL fission_width(Z0,A0,Ex0,J0)
	IF (i == i_channels) THEN
	  factor = 0.                 !fission wasn't possible
	ELSE
	  gam_sym = gamma(I_channels)
	  I_channels = i
	  factor = gam_sym/gam_tot     !factor to multiply fission widths by
	END IF 
	DO i=I_light+1,I_channels
	   if(Za(i) >= zmin) gamma(i) = gamma(i)*factor
	END DO
	END SUBROUTINE Normalise_fission

	SUBROUTINE complex_frag(A0,Z1,A1,Z2,A2,Ex0,J0,kl)
!
! FUNCTIONAL DESCRIPTION:
!
! this subroutine calculates the decay width for a particular binary decay
! exit channel.
! the variable kl is set to zero if this exit channel is not open
!
! DUMMY ARGUMENTS:
!
! INPUT:
!      A0 (REAL (KIND=r4) ::) nucleon number of decaying nucleus
!      Ex0 (REAL (KIND=r4) ::) excitation energy of decaying nucleus
!      J0 (REAL (KIND=r4) ::) spin of decaying nucleus
!      Z1 (REAL (KIND=r4) ::) proton number of lighter emitted complex fragment
!      A1 (REAL (KIND=r4) ::) nucleon number of lighter emitted complex fragment
!
!      Z2 (REAL (KIND=r4) ::) proton number of heavier emitted complex fragment
!      A2 (REAL (KIND=r4) ::) nucleon number of heavier emitted complex fragment
!
!  output:
!
!    kl (INTEGER ::) = 0 if exit channel not open
!                    =1 when open
!
! IMPLICIT INPUTS:
!
! I_channel - previous channel number
!
! IMPLICIT OUTPUTS:
!
! I_channels - new channel number
! gamma(I_channels) - decay width of channel
! Za(I_channels) - z-value of lightest emitted complex fragment
! Aa(I_channels - A-value of lightest emitted complex fragment
! U_therm(I_channels) - thermal excitation energy a saddle-point
!
! SIDE EFFECTS:
!
!    none
!
!
	INTEGER, INTENT(OUT) :: kl
	REAL (KIND=r4), INTENT(IN) :: A0,Z1,A1,Z2,A2,J0,Ex0
	REAL (KIND=r4) :: E_sad,U_saddle,aden,ef,temp,sigma_2,gam_k,
     &				  gamma_extra
	REAL (KIND=r4) :: M_Inertia_eff,E_rot_k,U,temp_saddle,sigma_k_2
	REAL (KIND=r4) :: surf,entropy,gamm,k

!subroutine saddle calculates the  saddle point energy of the saddle point configuration (E_sad)
	 CALL saddle(Z1,A1,Z2,A2,E_sad,m_inertia_eff)

!In principle  E_rot_saddle >= E_rot_yrast.
!or equivalently the fission barrier is always positive.
!Let us just check this in case the saddle point
!interpolation or more importantly, the extrapolation
!(if you are going to spins where the symmetric fission barrier is zero) is wrong.

	IF (E_sad < E_rot_yrast) E_sad = E_rot_yrast 
	U_saddle = Ex0 - E_sad

	Kl = 1
	IF (U_saddle <= 0.0) THEN
	  kl = 0
	  RETURN
	END IF

!we are using experimental masses for the ground state
!however, barrier is relative to liquid drop mass,
!so we remove shell and pairing terms from the excitation energy

	if (Z1+Z2 > Z_shell) U_saddle = U_saddle+shell_corr0+pair_energy0

	IF (U_saddle <= 0.0) THEN
	  kl = 0
	  RETURN
	END IF


	surf = (A1**0.66666+A2**0.66666)/A0**0.6666666      
	CALL get_little_a(A0,Z1+Z2,j0,U_saddle,surf,0.0,0.0,aden,entropy,
     &				  0.0,1.0)
	ef = entropy

!one may consider including the following line.
!As U_saddle is <= U_yrast, the only time the condition is met
!is a large U values when Aden_saddle > aden_yrast.
!IF (ef > ee) ef = ee

	IF (ee-ef > 65.0) THEN
	   kl = 0
	   RETURN
	END IF

	temp = (U_saddle/aden)**0.5
	sigma_2 = M_inertia_0*temp/40.848
	 gamm = temp/(1.0+U_saddle**1.25*sigma_2**1.5)/aden**0.25*
     &		(2.0*j0+1.0)/pe*EXP(ef-ee)


!include lestones modification to Bohr-Wheeler formulism by
!including summation over K - the projection of the angular
!momentum on the symmetry axis, the previous decay width is
!just for k=0, so now add contribution from other K values
	IF (para%lestone .AND. J0 > 0) THEN  
!    m_inertia_k = m_inertia_1 + m_inertia_2
!    M_Inertia_eff = 1.0/(1.0/M_inertia_k - 1.0/M_inertia_saddle )


        IF (MODULO(NINT(A0),2) /= 0 .AND. para%ssum) THEN
            gamma_extra = 0.0
            k = 0.5
        ELSE
            gamma_extra = gamm
            k = 1.0
        END IF

       loop_k: DO 
        E_rot_k = 20.424*K**2./M_inertia_eff
        U = U_saddle - E_rot_k
        IF (U <= 0) THEN
            EXIT loop_k
        ELSE
      CALL get_little_a(A0,Z1+Z2,j0,U,surf,0.0,0.0,aden,entropy,0.0,1.0)
            
            ef = entropy
            IF (ef-ee > 69) THEN
                EXIT loop_k
            ELSE
                temp_saddle = (U/aden)**0.5
                sigma_k_2 = temp_saddle*m_inertia_0/40.848
                gam_k = temp_saddle/(1.0+sigma_k_2**1.5*U**1.25)
     &                /aden**0.25/pe*EXP(ef-ee)*(2.0*j0+1.0)
!                K can be positive or negative
                 gam_k = 2.0*gam_k
                 gamma_extra = gamma_extra + gam_k
                 IF (gam_k <= gamma_extra/100.0) EXIT loop_k
            END IF
        END IF
        k = k + 1
        IF (INT(k) > INT(J0)) EXIT loop_k
        END DO loop_k
        gamm =  gamma_extra/(2.0*J0+1.0)
	END IF

	if (gamm <= 0) RETURN
	I_channels = I_channels + 1
	Za(I_channels) = nint(Z1)
	Aa(I_channels) = nint(A1)
	U_therm(I_channels) = Ex0 - E_sad
	gamma(I_Channels) = gamm
	
	RETURN
	END SUBROUTINE complex_frag
 
	SUBROUTINE fission_width(Z0,A0,Ex0,J0,fast_fission)
!
! FUNCTIONAL DESCRIPTION:
!
! Calculated the decay-widths for fission decay by the transition state formalism.
!
! DUMMY ARGUMENTS:
!
! input: Z0 (REAL (KIND=r4) ::) proton number of decaying nucleus
!        A0 (REAL (KIND=r4) ::) nucleon number of decaying nucleus
!        Ex0 (REAL (KIND=r4) ::) excitation energy of decaying nucleus
!        J0 (REAL (KIND=r4) ::) spin of decaying nucleus
!
! IMPLICIT INPUTS:
!
!
! IMPLICIT OUTPUTS:
!
! Zm - number of distinct exit channels classified by Z-value.
! gamma_z(3:Zm) - decay widths as a function of the smaller complex fragment
!                 Z-value.
! I_channels - total number of open decay channels.
! gamma(1:I_channels) - decay width for each exit channel
!
! SIDE EFFECTS:
!
! none
!
!
	USE sierk_stuff
	REAL (KIND=r4), INTENT(IN) :: Z0,A0,Ex0,J0
	LOGICAL, INTENT(OUT), OPTIONAL :: fast_fission

	INTEGER :: i,jj,iz1,ia1,iz2,ia2
	REAL (KIND=r4) :: barrier,E_rot,U_saddle,U_saddle0
	REAL (KIND=r4) :: aden,ef,Ek,fract,prob,k,j00,jmax
!REAL (KIND=r4) :: omega_0,beta
	REAL (KIND=r4) :: E_saddle,temp_saddle,sigma_2
	REAL (KIND=r4) :: U,M_Inertia_eff,E_rot_k,gam_k,gamma_extra,
     &				  M_inertia_k
	REAL (KIND=r4) :: surf,entropy,mom_inertia_1,mom_inertia_3
	REAL (KIND=r4) :: E_scission,entropy_sp,m_inertia_saddle
	REAL (KIND=r4), PARAMETER :: de = 0.5
	LOGICAL, SAVE :: first = .TRUE.

	IF (PRESENT(fast_fission)) fast_fission = .TRUE. !suj
	entropy_sp = 0
!Sierks routine does work if Z is less than 19
	IF (z0 < 19) RETURN       
!find the saddle-point energy from Sierk"s finite range model fits
	CALL yrast(Z0,A0,J0,0.0,0.0,E_rot)
	CALL Jmax_sierk(Z0,A0,jmax)
	if (j0 > jmax) THEN
	  barrier = 0.
	  j00 = jmax 
	ELSE
	  CALL bfis_sierk(J0,barrier)
	  if (para%b_scale > 0) barrier = barrier*para%b_scale
	  j00 = j0
	END IF 
	E_saddle = E_rot + Barrier
	U_saddle = Ex0 - E_saddle


!we are using experimental masses for the ground state
!however, barrier is relative to liquid drop mass,
!so we remove shell and pairing terms from the excitation energy

	if (z0 > Z_shell) U_saddle = U_saddle + shell_corr0 + pair_energy0
	U_saddle0 = U_saddle

	IF (U_saddle <= 0.0001) THEN
	   RETURN
	END IF

!surf = 2.0**0.3333333
	CALL bs_sierk(J00,surf)

	CALL get_little_a(A0,Z0,J0,U_saddle,surf,0.0,0.0,
     &  	 aden,entropy_SP,0.0,para%a_scale)

	IF (PRESENT(fast_fission)) fast_fission = .FALSE.    
	IF (entropy_sp >= entropy0 - 1.0) THEN
	    IF (PRESENT(fast_fission)) fast_fission = .TRUE.
	    IF (first) THEN
       WRITE (UNIT=*, FMT=*) "entropy at sp > entropy at gs - 1"
       WRITE (UNIT=*, FMT=*) "use of transmission state form not valid"
       WRITE (UNIT=*, FMT=*) "starting at JCN=",J0
       WRITE (UNIT=*, FMT=*) " "
       first = .FALSE.
       END IF
	END IF

	ef = entropy_sp

	IF (ee-ef > 69.0) THEN
		  RETURN
	END IF


	ek = 0.0001
	jj = 0
	loop_ek: DO 
	   U_saddle = U_saddle0 - ek
	   jj = jj + 1 

!   normal fermi gas level densities for spherical nucleus
	  CALL get_little_a(A0,Z0,j0,U_saddle,surf,0.0,0.0,aden,entropy,
     &    0.0,para%a_scale)

       ef = entropy
       Temp_saddle = ((U_saddle)/aden)**0.5          !here
       sigma_2  = m_inertia_0*temp_saddle/40.848

       storage(jj)%gam = de/(1.0+U_saddle**1.25*sigma_2**1.5)
     &            /aden**0.25*(2.0*j0+1.0)/pe*EXP(ef-ee)*collective


!   include lestones modification to Bohr-Wheeler formulism by
!   including summation over K - the projection of the angular
!   momentum on the symmetry axis, the previous decay width is
!   just for k=0, so now add contribution from other K values

       IF (para%lestone .AND. J0 > 0.0) THEN  
       CALL m_inertia_sierk(j00,mom_inertia_1,m_inertia_3=mom_inertia_3)
        m_inertia_saddle = mom_inertia_1*0.4*r0ld**2.0*A0**1.6666666666
        m_inertia_k = mom_inertia_3*0.4*r0ld**2.0*A0**1.6666666666
        M_Inertia_eff = 1.0/(1.0/M_inertia_k - 1.0/M_inertia_saddle)

        IF (MODULO(NINT(A0),2) /= 0 .AND. para%ssum) THEN
            gamma_extra = 0.0
            k = 0.5
        ELSE
            gamma_extra = storage(jj)%gam
            k = 1.0
        END IF

        loop_k: DO 
            E_rot_k = 20.424*K**2.0/M_inertia_eff
            U = U_saddle - E_rot_k
            IF (U <= 0.0) THEN
                EXIT loop_k
            ELSE
            CALL  get_little_a(A0,Z0,j0,U,surf,0.0,0.0,aden,entropy,
     &            0.0,para%a_scale)

                ef = entropy
                IF (ef-ee > 69.0) EXIT loop_k
                temp_saddle = (U/aden)**0.5
                sigma_2  = m_inertia_0*temp_saddle/40.848
                gam_k = de/(1.0+U**1.25*sigma_2**1.5)
     &                /aden**0.25*(2.0*j0+1.0)/pe*EXP(ef-ee)*collective
!                K can be positive or negative
                gam_k = 2.0*gam_k
                gamma_extra = gamma_extra + gam_k
                IF (gam_k <= gamma_extra/100.0) EXIT loop_k
            END IF
            k = k + 1
            IF (INT(k) > INT(J0)) EXIT loop_k
        END DO loop_k
        storage(jj)%gam =  gamma_extra/(2*j0+1)
       END IF 


       storage(jj)%energy = ek
       ek = ek + de
       IF (jj > 1) storage(jj)%gam = storage(jj)%gam + storage(jj-1)%gam 
       IF (storage(jj)%gam < storage(1)%gam/100.0) EXIT loop_ek
       IF (U_Saddle0 - ek <= 0.0) EXIT loop_ek
	END DO loop_ek       

	if (storage(jj)%gam == 0) RETURN
!at this point we have a non-zero fission decay width      
	I_channels = I_channels + 1
	Za(I_channels) = INT(Z0/2.0)
	Aa(I_channels) = INT(A0/2.0)
	gamma(I_channels) = storage(jj)%gam
	if (para%imf_option == 3) return ! only interested in the gamma not the rest

!now find the kinetic energy in the fission degree of freedom at the saddle point

	i = 1
	CALL RANDOM_NUMBER(prob)
      
	DO
	  fract = storage(i)%gam/storage(jj)%gam
	  IF (fract > prob) EXIT
	  i = i + 1
	END DO

	ek = storage(i)%energy

	IF (I_channels > n_light_ch) THEN
	   WRITE (UNIT=*, FMT=*) "N_light_ch too small in stat_mod.f90"
	   STOP 
	END IF 
	U_therm(I_channels) = U_saddle0 - ek
	dau(I_channels)%E_kinetic = ek


!IF (para%Kramer) THEN
!    omega_0 = 2.207-2.5e-4*J0**2.0
!    beta = beta2*temp_saddle**2.0
!    eta = beta/omega_0/2.0
!    Kramer_fact = (1.0+eta**2.)**0.5 - eta
!    gamma(I_channels) = gamma(I_channels)*Kramer_fact
!END IF

	 IF (para%Kramer /= 0) gamma(I_Channels) = 
     &	 gamma(I_Channels)*para%kramer

	gamma_z(Za(I_channels)) = gamma(I_channels)

	gamma_z(3:Za(I_channels)-1) = 0.0

	Iz_max(Za(I_channels)) = I_channels
	Iz_min(Za(I_channels)) = I_channels
	Zm = Za(I_channels)

! U_thermal is the energy above the the saddle point
! energy - in old versions of GEMINI this thermal energy
! is subdivided between the two fission fragments - the extra
! energy dissapated during the saddle to scission transition
! is ignored. This is expected to be good for light systems
! where the saddle and scission points are almost degerate.
! let us try and add this extra energy.

	CALL viola(INT(Z0),INT(A0),E_scission,Ek,iz1,iz2,ia1,ia2)
	U_scission = Ex0 - E_scission
	dau(I_channels)%l_p_s = ek

	 If ( .NOT. (para%i_angle .AND. para%quantum))
     &	U_therm(I_channels) = U_scission

	RETURN
	END SUBROUTINE fission_width


	SUBROUTINE fluctuate(Z0,A0,Ex0,J0,k0,ratio)

!       THis program considers the effect of thermal shape
!       fluctuations on the decay mode of the excited system.
!       The shape considers are spheriodal deformation and the
!       deformation parameter is the dimensional relative
!       quadpole moment (Q_rel) related to the absolute moment by
!       Q_abs = Q_rel*3./4./pi*Z*e*R0
!       where R0 =r0*A**0.33333
!       however I will also use the deformation quantities
!       ABS(ex**2) = 1.- 1./ratio**2. (prolate)
!                 = 1./ratio**2 - 1  (oblate)
!       where ratio is ratio of symmetry axis(c) to other axes(a).
!       where ratio = c/a
!       ratio >1 prolate
!       ratio <1 oblate
!       Q_rel = 8./15.*pi*(c/R0)**2.*ex**2.
!             = 8./15.*pi/R0**2.*(c**2.-a**2.)
!             = 8./15.*pi*(c_red**2.-a_red**2), where c_red = c/R0
	USE sierk_stuff
	INTEGER :: increase_q,jj,i
	INTEGER, PARAMETER :: n_fluct = 10000
	LOGICAL :: odd,sphere_less

	REAL (KIND=r4), INTENT(IN) :: A0,Z0,Ex0,j0
	REAL (KIND=r4), INTENT(OUT) :: k0,ratio
	REAL (KIND=r4) :: k0_start,k0_stop,fract,ran,delta_k0
	REAL (KIND=r4) :: q,q_start,Es0,Ec0,er0_sphere,e_def_sphere
	REAL (KIND=r4) :: E_def,x,entropy,aden,try,bc,bs
	REAL (KIND=r4) :: dbc,dbs,d2bc,d2bs,g,dg,d2g,e_def_old
	REAL (KIND=r4) :: fact1,fact2,fact3,fact4,fact6,ex
	REAL (KIND=r4) :: c_red,a_red,br_perp,br_para
	REAL (KIND=r4) :: dbr_perp,dbr_para,d2br_perp,d2br_para
	REAL (KIND=r4) :: uu0,uu00,mom_integral,er0_perp,er0_para
	REAL (KIND=r4) :: de,d2e,b,dqdx,e_rot,ee_ek
	REAL (KIND=r4) :: level_density,temp
	REAL (KIND=r4) :: prob_max,prob_tot
	REAL (KIND=r4) :: barrier,mom_inertia_1,mom_inertia_3
	REAL (KIND=r4) :: e_saddle,m_inertia_k
	REAL (KIND=r4) :: m_inertia_eff,m_Inertia_saddle
	REAL (KIND=r4), DIMENSION(n_fluct) :: prob,rat,kkk
	REAL (KIND=r4), PARAMETER :: delta_q = 0.2
	REAL (KIND=r4), PARAMETER :: q_max = 5.0
!REAL (KIND=r4)  :: hbar1,hbar2,db_dx,det,d2et
	COMPLEX (KIND=r8) :: stuff,stuff2,stuff3,stuff5,stuff6

	k0 = -1.0
!get fission barrier
	CALL yrast(Z0,A0,J0,0.0,0.0,E_rot)
	CALL Jmax_sierk(Z0,A0)
	CALL bfis_sierk(J0,barrier)
	E_saddle = barrier + E_rot
	CALL m_inertia_sierk(j0,mom_inertia_1,m_inertia_3=mom_inertia_3)
	m_inertia_saddle = mom_inertia_1*0.4*r0ld**2.0*A0**1.6666666666
	m_inertia_k = mom_inertia_3*0.4*r0ld**2.0*A0**1.6666666666
	M_Inertia_eff = 1.0/(1.0/M_inertia_k - 1.0/M_inertia_saddle)


	odd = .FALSE.
	if (MODULO(NINT(a0),2) /= 0 .AND. para%ssum) odd = .TRUE.

	Es0 = 17.9439*(1.0-1.7826*((A0-2.*Z0)/A0)**2.0)*A0**0.66666666
	Ec0 = 0.7053*Z0**2.0/A0**0.33333333
	er0_sphere = 34.540/A0**(5.0/3.0)

	sphere_less = .FALSE.
	prob_max = 0.0
	jj = 0
	q_start = 0
	loop_increase_decrease_q: DO increase_q = -1,1,2
	    q = q_start        
	    loop_q: DO 

!         include the following line and one only has emission from a sphere
!        if (q /= q_start) EXIT loop_increase_decrease_q

        if (q == 0) THEN              !sphere
            bc = 1.0
            bs = 1.0
            br_perp = 1.0
            br_para = 1.0
            x = 0.0

            ratio = 1.0
!           derrivatives with respect to x=ex**2.
            de = 0.0
            d2e = es0*4.0/45.0 - ec0*2.0/45.0

            dbr_perp = -1.0/6.0
            d2br_perp = -7.0/18.0 
            dbr_para = 1.0/3.0
            d2br_para = 4.0/9.0 

            e_rot = j0*(j0+1.0)*er0_sphere 
            uu0 = ex0 - e_rot
            CALL get_little_a(A0,Z0,j0,UU0,bs,0.0,E_rot,
     &            aden,entropy0,0.0,1.0)

           IF (UU0 <= 0.0) THEN                 ! this is not possible
               q = q + increase_q*delta_q
               e_def_old = e_rot
               e_def_sphere = e_rot
               sphere_less = .TRUE.
               CYCLE loop_q
           END IF

           !prepare pe and ee so we can use get_level_density
           pe = 1.0
           ee = entropy0
          
           M_inertia_0 = 0.4*r0**2.0*A0**1.6666666

           CALL get_level_density(uU0,J0,aden,M_inertia_0,
     &          entropy0,level_density)

           !final command for level density preparation
           pe = level_density*6.2831

            k0_stop = j0
            if (odd) THEN
                k0_start = 0.5
            ELSE
                k0_start = 0.0
            END IF
            delta_k0 = 1.0

        ELSE

            try = 0.25 - q**2.0*q*(5.0/8.0/pi)**3.0
            if (try > 0) THEN
                try = try**0.5
                if (q < 0)THEN
                c_red = (0.5+try)**(1.0/3.0) - ABS(0.5-try)**(1.0/3.0)
                ELSE
                c_red = (0.5+try)**(1.0/3.0) + ABS(0.5-try)**(1.0/3.0)
                END IF                        
            ELSE

                 stuff = CSQRT(CMPLX (try))

                 stuff2 = cmplx(0.5) + stuff
                 stuff3 = stuff2**CMPLX (1.0/3.0)
!                        WRITE (UNIT=*, FMT=*) stuff2,stuff3,stuff3**(3.0)
                 stuff5 = cmplx(0.5) - stuff
                 stuff6 = stuff5**CMPLX (1.0/3.0)
!                        WRITE (UNIT=*, FMT=*) stuff5,stuff6,stuff6**(3.0)

  
!                       WRITE (UNIT=*, FMT=*) stuff3+stuff6
!                       WRITE (UNIT=*, FMT=*) -0.5*(stuff3+stuff6)&
!                       +CMPLX(0.0,3.0**0.5/2.0)*(stuff3-stuff6)
!                       WRITE (UNIT=*, FMT=*) -0.5*(stuff3+stuff6)&
!                       -CMPLX(0.0,3.0**0.5/2.0)*(stuff3-stuff6)

                c_red = REAL(stuff3+stuff6)
                if (c_red <= 0) THEN
                     WRITE (UNIT=*, FMT=*) "cubic solution is neg?"
                END IF 
            END IF 



            a_red = (1.0/C_red)**0.5
            ratio = c_red/a_red
            x = 1.0-1.0/ratio**2.0
            ex = ABS(x)**0.5

            IF (q > 0) THEN !prolate

!               Coulomb, first and second derrivative
                fact1 = (1.0-ex**2.)**0.33333333333
                fact2 = LOG((1.0+ex)/(1.0-ex)) 
                Bc = 0.5*fact1*fact2/ex


                dBc = -fact1/2.0/ex**2.0*fact2 - fact2/3.0/fact1**2.0
     &                + 1.0/fact1**2.0/ex
       

                d2Bc =  fact2/3.0/ex/fact1**2.0 + fact1/ex**3.0*fact2
     &          -4.0/9.0*ex/fact1**5.0*fact2 - 2.0/ex**2.0/fact1**2.0
     &          + 2.0/3.0/fact1**5.0 
                 
!               surface, first and second derrivates
                G = ASIN(ex)/ex/fact1**1.5
                fact3 = 1.0 + G
                Bs = fact1/2.0*fact3 
                

                fact4 = (1.0-ex**2.0)
                dG =   1.0/ex/fact4 - ASIN(ex)/ex**2.0/fact4**0.5
     &                + ASIN(Ex)/fact4**1.5 
                dBs = -ex/3.0/fact1**2.0*fact3 + fact1*dG/2.0

                d2G = -2.0/ex**2.0/fact4 + 3.0/fact4**2.0 
     &                + 2.0*ASIN(Ex)/ex**3.0/fact4**0.5 
     &                - ASIN(ex)/ex/fact4**1.5
     &                + 3.0*ex*ASIN(ex)/fact4**2.5
                d2Bs = -(3.0+ex**2.0)/9.0/fact1**5.0*fact3
     &                - 2.0*ex/3.0/fact1**2.0*dG
     &                + fact1/2.0*d2G

!               rotational perp
                fact6 = 2.0-ex**2.0
                br_perp = 2.0*fact1**2.0/fact6

                dbr_perp = -8.0/3.0*ex/fact1/fact6
     &                + 4.0*ex*fact1**2./fact6**2.0

                d2br_perp = -8.0/3.0/fact1/fact6
     &                -16.0/9.0*ex**2.0/fact1**4.0/fact6
     &                -32.0/3.0*ex**2./fact1/fact6**2.0
     &                + 4.0*fact1**2.0/fact6**2.0
     &                + 16.0*ex**2.0*fact1**2.0/fact6**3.

!               rotational para
                br_para = 1.0/fact1

                dbr_para = 2.0*ex/3.0/fact1**4.0

                d2br_para = 2.0/3.0/fact1**4.0 +
     &                8.0/9.0*ex**2.0/fact1**7.0


!               derrivative of deformation energy with respect to ex
                de = dBs*Es0 + dBc*Ec0 
                d2E = d2Bs*Es0 + d2Bc*Ec0               

!               change variable to x = ex**2                            
!               valculate new derrivatives
                d2e = d2e*0.25/x - de*0.25/x**1.5
                de = de*0.5/ex

                d2br_perp = d2br_perp*0.25/x - dbr_perp*0.25/x**1.5
                dbr_perp = dbr_perp*0.5/ex

                d2br_para = d2br_para*0.25/x - dbr_para*0.25/x**1.5
                dbr_para = dbr_para*0.5/ex

                k0_stop = j0
                if (odd) THEN
                    k0_start = 0.5
                ELSE
                    k0_start = 0.0
                END IF
                delta_k0 = 1.0


                if (sphere_less) THEN
                   E_def = (Bs-1.0)*Es0 + (bc-1.0)*Ec0
                   k0 = k0_start
                   Er0_perp = (j0*(j0+1.0)-k0**2.0)*er0_sphere
                   Er0_para = k0**2.0*er0_sphere
                   E_rot = br_para*er0_para + br_perp*er0_perp
                   E_def = e_def + e_rot
                   if (e_def > e_def_old) THEN
                       EXIT loop_increase_decrease_q
                   END IF 
                   if (e_def > ex0) THEN
                      q = q + increase_q*delta_q
                      E_def_old = e_def
                      CYCLE loop_q
                   ELSE
                      prob_max = 0.0
                      sphere_less = .FALSE.
                   END IF 
                END IF 

            ELSE   !oblate
!               Coulomb, first and second derrivative
                fact1 = (1.0+ex**2.0)**0.33333333333
                fact2 = ATAN(ex)
                Bc = fact1*fact2/ex

            dBc = -bc/ex + 2.0/3.0*fact2/fact1**2.0 + 1.0/ex/fact1**2.
       

                d2Bc = bc/ex**2.0 - dbc/ex  -2.0/3.0/fact1**5.0
     &         - 8.0/9.0*ex*fact2/fact1**5.0 - 1.0/ex**2.0/fact1**2.0  


!              surface
               fact3 = (1.0+ex**2.0)**0.5
               g = LOG(ex+fact3)/ex/fact3
               Bs = fact1/2.0*(1.0+g)

               dg = (1.0+ex/fact3)/(ex+fact3)/ex/fact3 - g/ex
     &               -g*ex/(1.0+ex**2.0)

               d2g = - dg/ex +g/ex**2. - g/(1.0+ex**2.0) 
     &               -dg*ex/(1.0+ex**2.0)
     &              + 2.0*g*ex**2.0/(1.0+ex**2.0)**2.0
     &         + (1.0/fact3 - ex**2.0/fact3**3.0)/(ex+fact3)/ex/fact3
     &         - (1.0+ex/fact3)**2.0/ex/fact3/(ex+fact3)**2.0
     &         - (1.0+ex/fact3)/(ex+fact3)/ex**2.0/fact3
     &         - (1.0+ex/fact3)/(ex+fact3)/fact3**3.0
      
               dbs =  (1.0+g)/3.0/fact1**2.0*ex + fact1/2.0*dg
      
             d2bs = (1.0+g)/3.0/fact1**2.0 + 2.0*ex*dg/3.0/fact1**2.0
     &          - 4.0/9.0*ex**2.0*(1.0+g)/fact1**5.0 +fact1/2.0*d2g

!              rotational_para
               Br_para = 1.0/fact1
               dBr_para = -2.0*ex/3.0/fact1**4.0
         d2Br_para = -2.0/3.0/fact1**4.0 + 16.0/9.0*ex**2.0/fact1**7.0
      
!               rotational perp
                fact6 = 2.0+ex**2.0
                br_perp = 2.0*fact1**2.0/fact6

                dbr_perp = 8.0/3.0*ex/fact1/fact6
     &                - 4.0*ex*fact1**2.0/fact6**2.0

                d2br_perp = 8.0/3.0/fact1/fact6
     &               -16.0/9.0*ex**2.0/fact1**4.0/fact6
     &               -32.0/3.0*ex**2.0/fact1/fact6**2.0
     &               - 4.0*fact1**2.0/fact6**2.0
     &               + 16.0*ex**2.0*fact1**2.0/fact6**3.0

!               derrivative of deformation energy with
!               respect to ex
                de = dBs*Es0 + dBc*Ec0 
                d2E = d2Bs*Es0 + d2Bc*Ec0               

!               change variable to x = -ex**2                           
!               valculate new derrivatives
                d2e = -d2e*0.25/x - de*0.25/(-x)**1.5
                de = -de*0.5/ex

                d2br_perp = -d2br_perp*0.25/x - dbr_perp*0.25/(-x)**1.5
                dbr_perp = -dbr_perp*0.5/ex

                d2br_para = -d2br_para*0.25/x - dbr_para*0.25/(-x)**1.5
                dbr_para = -dbr_para*0.5/ex


                k0_start = j0
                if (odd) THEN
                    k0_stop = 0.5
                ELSE
                    k0_stop = 0.0
                END IF
                delta_k0 = -1.0

                if (sphere_less) THEN
                   E_def = (Bs-1.0)*Es0 + (bc-1.0)*Ec0
                   k0 = k0_start
                   Er0_perp = (j0*(j0+1.0)-k0**2.0)*er0_sphere
                   Er0_para = k0**2.0*er0_sphere
                   E_rot = br_para*er0_para + br_perp*er0_perp
                   E_def = e_def + e_rot
                   if (e_def > e_def_old) THEN
                        EXIT loop_q
!                       WRITE (UNIT=*, FMT=*) "min E_def+e_rot=",E_def_old,&
!                      " is larger than ex0=",ex0
!                       STOP
                   END IF 
                   if (e_def >= ex0) THEN
                      q = q + increase_q*delta_q
                      E_def_old = e_def
                      CYCLE loop_q
                   ELSE
                       uu0 = ex0 - e_def

  !prepare pe and ee so we can use get_level_density
                       pe = 1.
                       ee = entropy0
          
                       M_inertia_0 = 0.4*r0**2.0*A0**1.6666666

                       CALL get_level_density(uU0,J0,aden,M_inertia_0,
     &                      entropy0,level_density)


                       !final command for level density preparation
                        pe = level_density*6.2831
                       sphere_Less = .FALSE.
                   END IF 
                END IF 
            END IF 
        END IF                

!        deformation energy
         E_def = (Bs-1.0)*Es0 + (bc-1.0)*Ec0  
         uU0 = ex0 - e_def
         if (uU0 <= 0) EXIT loop_q

!        inertia parameter for x=ex**2.
         B = (3.0-x)/90.0/(1.0-x)**2.66666666666
!         db_dx = -1.0/90.0/(1.0-x)**2.66666666666&
!         + 4.0*(3.0-x)/135.0/(1.0-x)**3.66666666666


         dqdx = 8.0/15.0*pi/(1.0-x)**0.6666666
     &         + 2.0/3.0*8.0/15.0*pi*x/(1.0-x)**1.6666666

         k0 = k0_start
         loop_k: DO 
             Er0_perp = (j0*(j0+1.0)-k0**2.)*er0_sphere
             Er0_para = k0**2.0*er0_sphere

             E_rot = br_para*er0_para + br_perp*er0_perp

!           if fission unstable - forget
!           for prolate shapes check to see if rotation
!            plus deformation is greater than Sierk"s saddle point
!           values                                
             if (q > 0 .AND. e_def + e_rot
     &             > E_saddle +  20.424*K0**2./M_inertia_eff) THEN
                k0_start = k0 + 1.0
             ELSE

                 if (e_rot >= uU0) EXIT loop_k
                 uU00 = uu0 - e_rot
                 
!                        calculate the level density parameter
                 CALL get_little_a(A0,Z0,j0,uU00,bs,0.0,E_rot,
     &                aden,entropy,0.0,1.0)


                 IF (uu00 <= 0.0) EXIT loop_k
                 ee_ek = entropy - ee
                 IF (ABS(ee_ek) > 65.0 ) EXIT loop_k


                 CALL get_level_density(Uu00,j0,aden,M_Inertia_0,
     &                entropy,level_density)

                 temp = (uu00/aden)**0.5

                 mom_integral = (5.0/4.0*B/er0_sphere*temp)**0.5
!                 det = de + dbr_perp*er0_perp + dbr_para*er0_para
!                 d2et = d2e + d2br_perp*er0_perp&
!                 + d2br_para*er0_para

!                hbar correction
!                 hbar1 =  er0_sphere*4.0/5.0/b*det**2.0/24.0/temp**3.
!                 hbar2 = -er0_sphere/15.0/b*&
!                 (d2et - 0.5*det/b**0.5*db_dx)/temp**2.0

                 jj = jj + 1 
                 if (jj > n_fluct) THEN
                     WRITE (UNIT=*, FMT=*) "n_fluct too small"
                     STOP 
                 END IF 
!                 prob(jj) = mom_integral*(1.+hbar1+hbar2)*level_density
                 prob(jj) = mom_integral*level_density*br_para**0.5
!                this prob is calculated for x space - we want q space so
!                multiply by dxdq         
                 kkk(jj) = k0
                 rat(jj) = ratio
                 prob(jj) = prob(jj)/dqdx*delta_q
                  if (k0 /= 0) prob(jj) = 2.0*prob(jj)
                 prob_tot = prob_tot + prob(jj)
                 prob_max = MAX(prob(jj),prob_max)
                 if (prob(jj) <= 0.01*prob_max) EXIT loop_k
             END IF 
             k0 = k0 + delta_k0
             IF (k0 > k0_stop .AND. delta_k0 > 0.0) EXIT loop_k
             IF (k0 < k0_stop .AND. delta_k0 < 0.0) EXIT loop_k
         END DO loop_k
        if (k0 == k0_start .and. prob(jj) <= .01*prob_max) EXIT loop_q

        q = q + increase_q*delta_q
        if (abs(q) > q_max) EXIT loop_q
        END DO loop_q
	   q_start = q_start - increase_q*delta_q
	    e_def_old = e_def_sphere
	    if (e_def_sphere > ex0)sphere_less = .TRUE.
	END DO loop_increase_decrease_q

	if (jj == 0)  RETURN             !not possible

	if (jj > n_fluct) THEN
	    WRITE (UNIT=*, FMT=*) "increase n_fluct in fluctuate"
	END IF 

	DO i=2,jj
	   prob(i) = prob(i) + prob(i-1)
	END DO

	CALL RANDOM_NUMBER(ran)
	i = 1
	DO 
       fract = prob(i)/prob(jj)
       IF (fract > ran) EXIT
	  i = i + 1
	END DO

	k0 = kkk(i)
	ratio = rat(i) 
 
	END SUBROUTINE fluctuate


	SUBROUTINE viola(Z,A,E_scission,ek_viola,z1,z2,a1,a2)
	INTEGER, INTENT(IN) :: Z,A
	REAL (KIND=r4), INTENT(OUT) :: E_scission,ek_viola
	INTEGER, INTENT(OUT) :: z1,z2,a1,a2
	REAL (KIND=r4) :: mass1,mass2,Qfiss

	z1 = Z/2
	a1 = a/2
	CALL getmass(z1,A1,MASS_EXP=mass1)

	a2 = A - a1
	z2 = z - z1
	CALL getmass(z2,a2,MASS_EXP=mass2)

	Qfiss = mass0 - mass1 - mass2
	ek_viola = 0.1189*REAL(z)**2.0/REAL(a)**0.3333333333 + 7.3
	E_scission = Ek_viola - Qfiss

	END SUBROUTINE viola

	SUBROUTINE xray(parent,daughter,n_extra)
	TYPE (identity), INTENT(IN) :: parent
	TYPE (identity), INTENT(INOUT), DIMENSION(3) :: daughter
	INTEGER, INTENT(INOUT) :: n_extra
	REAL (KIND=r4) :: time_k_xray,gamma_k_xray

	IF (n_extra == 1) THEN
!  if K_vacancy not filled when residue at this point
!  let it decay now.
       IF (parent%Z >= para%Z_imf_min .AND. X_decay) THEN
      n_extra = n_extra + 1
      daughter(2) = daughter(1)
      daughter(1)%z = -1
      daughter(1)%a = parent%Z
      daughter(1)%ex = 0.0
      daughter(1)%exx = 0.0
      daughter(1)%level = 0
      daughter(1)%J = 0.0
      daughter(1)%m = 0.0
      daughter(1)%id = parent%id
      gamma_k_xray = 1.076e-11*REAL(parent%z)**3.93
      time_k_xray = exp_decay(gamma_k_xray)
      daughter(1)%time = time_k_xray + parent%time
      X_decay = .FALSE.
       END IF
       RETURN
	END IF

!K xray decay of vacancy produced in CN
	IF (daughter(1)%Z >= para%Z_imf_min) THEN
	  X_decay = .FALSE.
	  RETURN
	END IF 

	IF (x_decay .AND. parent%z > para%Z_imf_min) THEN
	   gamma_k_xray = 1.076e-11*REAL(parent%z)**3.93
	   time_k_xray = exp_decay(gamma_k_xray)

	   IF (daughter(1)%time-parent%time < time_k_xray) RETURN

!  insert a k xray into decay sequence
!  after the last parent nucleus
	   N_extra = n_extra + 1
	   daughter(3) = daughter(2)
	   daughter(2) = daughter(1)
	   daughter(1)%z = -1
	   daughter(1)%a = parent%Z
	   daughter(1)%ex = 0.0
	   daughter(1)%exx = 0.0
	   daughter(1)%J = 0.0
	   daughter(1)%id = parent%id
	   daughter(1)%level = 0
	   daughter(1)%time = time_k_xray + parent%time
	   x_decay = .FALSE.
	END IF
	END SUBROUTINE xray

	END MODULE statistical_mod
!------------------------------end:   MODULE statistical_mod-----------------------------------



!------------------------------begin:   MODULE gemsub-----------------------------------
	module gemsub
	USE gem_kind
	USE trig
	USE evaporation
	USE statistical_mod
	USE light_stuff
	USE angle
	IMPLICIT NONE
	PRIVATE

	PUBLIC :: event, decay_event, dmp, crash
	PUBLIC :: fusion, input, INPUT_ANA

	INTEGER, PARAMETER, PUBLIC :: max_frag = 100
	TYPE (identity), DIMENSION(max_frag), PUBLIC :: frag
	CHARACTER (LEN=44), PRIVATE :: evtfile, infile, outfile
	CHARACTER (LEN=44), PUBLIC ::  spectrafile

	INTEGER, PARAMETER, PUBLIC :: Amax = 4
	INTEGER, PUBLIC :: N_evt, Num_l, Max_l
	INTEGER, PUBLIC :: N_frag, N_binary



!	Fission and IMF parameters
	INTEGER	::	IMF_OPTION_ANA, Z_IMF_MIN_ANA
	REAL ::		T_DELAY_ANA, SIG_DELAY_ANA, KRAMER_ANA, B_SCALE
	LOGICAL ::	LESTONE_ANA, SHARP_DELAY_ANA

!	Evaporation plus gamma emission
	LOGICAL ::	TL_IWBC_ANA, K_SUM_ANA, POLARIZATION_ANA
	INTEGER ::	EXOTIC_INDEX_ANA, MASS_OPTION_ANA
	REAL ::		RATIO_ANA, e2_strength_ANA, e1_strength_ANA,
     &			THRESHOLD_ANA

!	Level density
	INTEGER ::	ADEN_TYPE_ANA
      REAL :: ADEN_K0_ANA, ADEN_KINF_ANA
	REAL ::		ADEN_0_ANA
	REAL ::		A_SCALE_ANA




	CONTAINS

	SUBROUTINE INPUT_ANA()

!	Fission and IMF parameters
	READ(2003,*)
	READ(2003,*)	IMF_OPTION_ANA
	READ(2003,*)	T_DELAY_ANA
	READ(2003,*)	SIG_DELAY_ANA
	READ(2003,*)	SHARP_DELAY_ANA
	READ(2003,*)	Z_IMF_MIN_ANA
	READ(2003,*)	KRAMER_ANA
	READ(2003,*)	LESTONE_ANA
	READ(2003,*)	B_SCALE
!	Evaporation plus gamma emission
	READ(2003,*)
	READ(2003,*)	TL_IWBC_ANA
	READ(2003,*)	K_SUM_ANA
	READ(2003,*)	POLARIZATION_ANA
	READ(2003,*)	EXOTIC_INDEX_ANA
	READ(2003,*)	MASS_OPTION_ANA
	READ(2003,*)	e2_strength_ANA
	READ(2003,*)	e1_strength_ANA
	READ(2003,*)	THRESHOLD_ANA
!	Level density
	READ(2003,*)
	READ(2003,*)	ADEN_TYPE_ANA
	READ(2003,*)	ADEN_0_ANA	
      READ(2003,*) ADEN_K0_ANA
      READ(2003,*) ADEN_KINF_ANA
	READ(2003,*)	A_SCALE_ANA

	RETURN
	END SUBROUTINE INPUT_ANA


	SUBROUTINE event(cn,weight)
!
! FUNCTIONAL DESCRIPTION:
!
! This is THE subroutine of GEMINI. This subroutine can be used in other
!  programs as an event simulator. It simulates the decay of a compound
!  nucleus and follows the decay chain down until particle emission is
!  impossible.
!
! DUMMY ARGUMENTS:
!
!    none
!
! IMPLICIT INPUTS:
!
! frag(1)%Z - proton number of compound nucleus
! frag(1)%A - nucleon number of compound nucleus
! frag(1)%Ex - excitation energy of compound nucleus (MeV)
! frag(1)%j  - spin of compound nucleus (hbar)
! frag(1)%spin% theta,phi - spin of compound nucleus 
!                            (polar angles)
! frag(1)%vel% v,theta,phi - velocity vector of compound nucleus
!                            (magnitude and polar angles) (cm/ns)
! flags:
!    I_weight - weighting on
!    I_angle - fragment velocity vectors are calculated
!    A_cut - minimum nucleon number for complex fragment decay
!
!
! IMPLICIT OUTPUTS:
! arrays:
! n_frag - total number of fragments in decay tree
! frag(1:N_frag)%Z - proton numbers of all fragments produced in decay chain
! frag(1:N_frag)%A - nucleon numbers . . .
! frag(1:n_frag)%Ex - excitation energies ....
! frag(1:n_frag)%spin - spin vectors .....
! frag(1:n_frag)%vel - velocity vectors .....
! frag(1:n_frag)%id - branch of decay tree
!
!
! final fragments are distinguished by Ex = -1.
!
! SIDE EFFECTS:
!
! uses a lot of cpu time
!

	TYPE (identity), INTENT(INOUT) :: cn
	REAL (KIND=r4), INTENT(out) :: weight
	REAL (KIND=r4) :: ratio

	N_FRAG = 0
	N_binary = 0
	weight = 1.0

	CALL init_decay ()

	mult1 = 0.0

	IF (para%k_sum) THEN
	   ratio = 0.0
	ELSE
	   ratio = para%ratio
	END IF 

!the following recursive subroutine decays the cn and all daughter fragments

	CALL decay_event(cn,ratio,weight)

	RETURN
	END SUBROUTINE event

	RECURSIVE SUBROUTINE decay_event(parent,ratio,weight)
	    TYPE (identity), INTENT(INOUT) :: parent
	    REAL (KIND=r4), INTENT(INOUT) :: ratio,weight

	    TYPE (identity), DIMENSION (3) :: daughter
	    INTEGER :: n_extra,i
		
C 		PRINT*, 'PARENT%Z, PARENT%A, PARENT%EX, PARENT%LEVEL'	!HG CHENG
C 		PRINT*, PARENT%Z, PARENT%A, PARENT%EX, PARENT%LEVEL

C 		PAUSE !HG CHENG 

	    CALL decay(parent,daughter,n_extra,ratio,weight,n_binary)

	     DO i=1,n_extra
        IF ((daughter(i)%level > 0 .AND. daughter(i)%ex .GE. 0.0).OR.		!HG CHENG
     &     (daughter(i)%A > amax .AND. daughter(i)%ex > 0.0 ))THEN
			
C 		PAUSE !HG CHENG

           CALL decay_event(daughter(i),ratio,weight)
        ELSE
           n_frag = n_frag + 1
C 		PRINT*, N_FRAG,DAUGHTER(I)%Z, DAUGHTER(I)%A, DAUGHTER(I)%EX !HG CHENG
           frag(n_frag) = daughter(i)
        END IF
	     END DO 

	END SUBROUTINE decay_event

	SUBROUTINE dmp(L,weight)
!
! FUNCTIONAL DESCRIPTION:
!
!      Compresses the arrays Z,A,etc and writes the
!      event in the event file
!
! DUMMY ARGUMENTS:
!
!    input: L (INTEGER ::) - l-wave
!
! IMPLICIT INPUTS:
!
!    N_frag - number of particle in arrays
!    Max_frag - maximum possible number in arrays
!    Z,A,Ex,I_decay,theta_v,phi_v,v - arrays containing information
!                                   about each fragment
!    I_open - flag indicating when event is open
!
! IMPLICIT OUTPUTS:
!
!    none
!
!
! SIDE EFFECTS:
!
!    increases length of event file
!
!
	INTEGER, INTENT(IN) :: L
	REAL (KIND=r4), INTENT(IN) :: weight
	INTEGER :: i
	REAL (KIND=r4) :: phi0,vx,vy,vz,vx_new,vy_new,vz_new
	REAL (KIND=r4) :: ran
	INTEGER (KIND=2) :: array(300),part2(2)=(/0,0/)
	INTEGER :: ipt,n_records
	INTEGER :: tot_a
! compress event arrays, so we only write out residual fragments.

	IF (N_frag > Max_frag) THEN
	  WRITE (UNIT=*, FMT=*) "Parameter Max_frag is too small. 
     &		(see gemsub.f90)"
	  STOP
	END IF

	CALL RANDOM_NUMBER(ran)
	phi0 = 360.0*ran

	tot_a = 0
	DO i = 1, N_frag
	   if (frag(i)%z >= 0) tot_a = tot_a + frag(i)%A
	   frag(i)%A = frag(i)%A + ISHFT(frag(i)%Z, 8)

!  rotate velocity vectors from z axis normal to reaction
!  plane to new z axis along the beam axis
	   if (para%quantum .AND.
     &	   .not. para%polarization .AND. para%i_angle) THEN
       vz = COSD(frag(i)%vel%theta)
       vx = COSD(frag(i)%vel%phi)*SIND(frag(i)%vel%theta)
       vy = SIND(frag(i)%vel%phi)*SIND(frag(i)%vel%theta)

       vz_new = vx
       vx_new = vz
       vy_new = vy

       frag(i)%vel%theta = ACOSD(vz_new)		!很明显是角度
       if (vy_new /= 0.0 .AND. vx_new /= 0.0) THEN
           frag(i)%vel%phi = ATAN2D(vy_new,vx_new)
!          randomize the reaction plane in phi angle
           frag(i)%vel%phi = MODULO(frag(i)%vel%phi + phi0,360.0)
           if (frag(i)%vel%phi < 0.0) frag(i)%vel%phi =
     &         frag(i)%vel%phi + 360.0
       ELSE
           frag(i)%vel%phi = 0.0
       END IF
	  END IF       
	END DO

	if (tot_a /= para%acn) THEN
	   stop "mass not conserved"
	END IF 


!open the event file at the end to append new event

	OPEN ( UNIT=2, ACCESS = "SEQUENTIAL", FORM = "UNFORMATTED",
     &       STATUS = "OLD", ACTION = "WRITE", POSITION = "APPEND",
     &	   FILE = evtfile)


	ipt = 2
	array(ipt) = INT(N_frag,KIND=2)
	ipt = ipt + 1
	array(ipt) = INT(N_binary,KIND=2)
	ipt = ipt + 1
	array(ipt) = INT(l,KIND=2)

	array(ipt+1:ipt+2) = TRANSFER(weight,part2)
	ipt = ipt + 2

	FORALL (i=1:n_frag) array(i+ipt) = INT(frag(i)%A,kind=2)
	ipt = ipt + n_frag

	FORALL (i=1:n_frag) array(i+ipt) = INT(frag(i)%id,kind=2)
	ipt = ipt + n_frag

	if (para%i_angle) THEN
       FORALL (i=1:n_frag) array(ipt+i*2-1:ipt+i*2) =
     &      TRANSFER(frag(i)%VEL%V,part2)
       ipt = ipt + n_frag*2
       FORALL (i=1:n_frag) array(ipt+i*2-1:ipt+i*2) =
     &      TRANSFER(frag(i)%VEL%theta,part2)
       ipt = ipt + n_frag*2
       FORALL (i=1:n_frag) array(ipt+i*2-1:ipt+i*2) =
     &      TRANSFER(frag(i)%VEL%phi,part2)
       ipt = ipt + n_frag*2
	END IF

	if (para%polarization) THEN
	   forall (i=1:n_frag) array(i+ipt) = INT(frag(i)%M*2.0,kind=2)
	   ipt = ipt + n_frag
	END IF

	if (para%time_flag) THEN
	     FORALL (i=1:n_frag) array(ipt+i*2-1:ipt+i*2) =
     &     TRANSFER(frag(i)%time,part2)
	     ipt = ipt + n_frag*2
	END IF 


	if (para%j_flag .OR. para%polarization) THEN
	     FORALL (i=1:n_frag) array(i+ipt) = INT(frag(i)%j*2.0,kind=2)
	     ipt = ipt + n_frag
	END IF 

	if (para%ex_flag) THEN
	     FORALL (i=1:n_frag) array(i+ipt) = INT(frag(i)%exx*10.)
	     ipt = ipt + n_frag
	END IF

	n_records = ipt
	array(1) = INT(n_records,kind=2)

	WRITE (UNIT=2) (array(i),i=1,n_records)
	PRINT*, 'number of fragments:', array(2)
	CLOSE (UNIT=2)

	RETURN
	END SUBROUTINE dmp

	SUBROUTINE fusion(L,CN)
!
! FUNCTIONAL DESCRIPTION:
!
! determines the z,a,excitation energy for each l wave and how many events to simulate
!
! DUMMY ARGUMENTS:
!
! input: L (INTEGER ::) l-wave number in hbar
!
! IMPLICIT INPUTS:
!
! Zcn - proton number of compound nucleus
! Acn - nucleon number of compound nucleus
! Excn - excitation energy of compound nucleus
!
! IMPLICIT OUTPUTS:
!
! Z(1) - proton number of compound nucleus
! A(1) - nucleon number of compound nucleus
! Ex(1) - excitation energy of compound nucleus
!
! J(1),theta_j(1),phi_j(1) - spin vector of compound nucleus, magnitude and
!                            polar angles
!
! V(1),theta_v(1),phi_v(1) - velocity vector of compound nucleus, magnitude
!                              and polar angles
!
! SIDE EFFECTS:
!
!    none
!
!
	INTEGER, INTENT(IN) :: L
	TYPE (identity), INTENT(out) :: cn       !compound nucleus

	LOGICAL :: odd
	REAL (KIND=r4) :: ran

	cn%Z = para%Zcn
	cn%A = para%Acn
	cn%Ex = para%Excn
	cn%id = 1
	cn%k = 0.

	CALL RANDOM_NUMBER(ran)
	if (para%dexcn > 0.0) THEN
	    cn%EX = cn%EX + para%dexcn*(1.0-2.0*ran)
	END IF
	cn%ex  = MAX(cn%ex,0.0_r4)

	IF (para%tl_iwbc) THEN  ! explicitly consider evaporated fragment spins
  ! values, must allow for non-integer spins
	   IF (MODULO(cn%A,2) /= 0 .AND. para%ssum) THEN   
!      odd A = have 1/2 integer value spins
       odd = .TRUE.                              
       CALL RANDOM_NUMBER(ran)
       IF (ran > 0.5) THEN  
                                    
          cn%J = REAL(L) + 0.5
       ELSE
           cn%J = ABS(REAL(L)-0.5)
       END IF
	   ELSE ! even A = integer values spins
       cn%J = REAL(L)
       odd = .FALSE.
	   END IF
	ELSE  !ignore evaporated fragment"s spin value, except for a degeneracy in the calculating gamma), hence we
!can use only integer value spins.
	   cn%J = REAL(L)
	END IF

	if (para%k_sum) THEN
!    cn%k = cn%j
	    if (odd) THEN
        cn%k =0.5
	    ELSE
        cn%k =0.0
	    END IF
	ELSE IF (para%ratio > 0.0) THEN
	     if (para%ratio > 1.0) THEN           !prolate
         if (odd) THEN
            cn%k =0.5
         ELSE
            cn%k =0.0
         END IF
	     ELSE                                 !oblate
         cn%k = cn%j
	     END IF 
	END IF 

	cn%time = 0.0
	cn%level = 0
	IF (para%quantum) THEN
	   cn%M = cn%J
	   cn%vel%theta = 90.0
	   cn%vel%phi = 0.0
	ELSE
	  cn%spin%theta = 90.0
	  cn%spin%phi = 0.0
	  cn%vel%theta = 0.0
	  cn%vel%phi = 0.0
	END IF

	cn%vel%v = 0.0

	RETURN
	END SUBROUTINE fusion

	SUBROUTINE input(ICZ,ICN,NEXC, ratio_ana)
!
! FUNCTIONAL DESCRIPTION:
! read input parameters for calculations
!
! DUMMY ARGUMENTS:
!
! none
!
! IMPLICIT INPUTS:
!
! none
!
! IMPLICIT OUTPUTS:
!
! Zcn - proton number of compound nucleus
! Acn - nucleon number of compound nucleus
! Excn - excitation energy of compound nuclei
! Lmin,lmax - minimum and maximum l-wave simulated in calculation
! R_num - gives the number of simulations per l-wave i.e. 2l(l+1)*R_num
! I_weight - flag for weighting
! I_angle - flag to calculate angles
! A_cut - minimum A for complex fragment emission
! reaction - title of simulation
! day - date
! Max_l - maximum l-wave already simulated if restarting
! Num_l - number of simulations already made for the above l-wave
!
! SIDE EFFECTS:
!
! none
!
!
	CHARACTER(LEN=40) :: condition
	CHARACTER(LEN=40) :: file_name
	CHARACTER(LEN=4), DIMENSION(4), PARAMETER :: 
     &              ext=(/".evt",".inp",".out",".spe"/)
	INTEGER :: L,I_read_error, ICZ, ICN, ICA, NEXC
	INTEGER (KIND=2), DIMENSION(4) :: parti2
	LOGICAL :: new=.FALSE., continue=.FALSE. ,more=.FALSE.

	REAL :: ratio_ana

	ICA=ICZ+ICN
	file_name = 'frag'

	file_name = ADJUSTL(file_name)       !remove leading blanks
	infile  = TRIM(file_name)//ext(2)
	evtfile = TRIM(file_name)//ext(1)
	spectrafile = TRIM(file_name)//ext(4)

	condition = 'new'
	condition = ADJUSTL(condition)       !remove leading blanks

	Max_L = -1

!*************** INPUT PARAMETERS *********************
	    para%Zcn= ICZ
	    para%Acn= ICA
	    para%EXcn=NEXC

	    para%dEXcn=0.
	    para%Lmin= 0.
	    para%Lmax= 0.
	    para%R_Num=100
	    para%N_weight=0.

!**************** output choices ***********************
	    para%diagnostics= .false.
	    para%time_flag=.false.
	    para%j_flag =  .false.
	    para%ex_flag=  .true.           !.false.
	    para%I_angle=  .true.           !.false.  velocity
	    para%quantum=  .false.
!**************** fission and IMF parameters ***********
	    para%imf_option =	IMF_OPTION_ANA
	    para%t_delay =		T_DELAY_ANA
	    para%sig_delay =	SIG_DELAY_ANA
	    para%sharp_delay =	SHARP_DELAY_ANA
	    para%Z_imf_min =	Z_IMF_MIN_ANA
	    para%kramer =		KRAMER_ANA
	    para%lestone =		LESTONE_ANA
	    para%b_scale =		B_SCALE
!**************** evaporation + gamma emission***********
	    para%tl_iwbc =		TL_IWBC_ANA
	    para%k_sum =		K_SUM_ANA
	    para%polarization =	POLARIZATION_ANA
	    para%exotic_index =	EXOTIC_INDEX_ANA
	    para%ratio =		RATIO_ANA
	    para%mass_option =	MASS_OPTION_ANA
	    para%e2_strength =	e2_strength_ANA
	    para%e1_strength =	e1_strength_ANA
	    para%threshold =	THRESHOLD_ANA
!**************** level density *************************
	    para%aden_type =	ADEN_TYPE_ANA
	    para%aden_0 =		ADEN_0_ANA
!   para%collective_enhancement
          para%aden_k0 = ADEN_K0_ANA
          para%aden_kinf = ADEN_KINF_ANA
	    para%a_scale =		A_SCALE_ANA



	   para%ssum = para%tl_iwbc
	   para%fluct = .FALSE.
	   if (para%k_sum) para%ssum = .FALSE.
	   if (para%k_sum) para%fluct = .TRUE.
	   if (para%ratio < 0.) THEN
       para%expand = .true.
       para%ratio = 0.
	   ELSE
       para%expand = .FALSE.
	   END IF 

	 IF (para%polarization .AND. (para%ratio /= 0. .AND. 
     &	 para%ratio /= 1.)) THEN
       WRITE (UNIT=*,FMT=*) "polarization option only available for"
       WRITE (UNIT=*,FMT=*) "spherical emisster, i.e. ratio = 0 or 1"
       STOP
	   END IF 

!the minimum z considered from IMF emission is allowed to be either 3 , 4 or 5.
	   para%z_imf_min = MAX(para%Z_imf_min,3)
	   para%Z_imf_min = MIN(para%Z_imf_min,5)

	   A_cut = 0
	   IF (para%z_imf_min == 3) THEN
	   A_cut = MAX(A_cut,8)       !minimum values of A_cut set to 8 by default
	   ELSE
	   A_cut = MAX(A_cut,12)      !minimum values of A_cut set to 12 by default
	   END IF


	   IF (para%E2_strength <= 0)para%e2_strength = 50.0
	   IF (para%E1_strength <= 0)para%e1_strength = 0.1

!decide to select decay channels by calculating times or channel probabilities.

	IF (para%imf_Option == 1) THEN
	    para%time_wise = .FALSE.
	ELSE IF (para%t_delay == 0.0) THEN
	    para%time_wise = .FALSE.
	ELSE IF (para%sig_delay == 0.0) THEN
	    para%time_wise = .FALSE.
	ELSE
	    para%time_wise = .TRUE.
	END IF

	IF (para%imf_option >= 2 .and. para%quantum) THEN
	   PRINT*, "imf_option>=2,quantum=.true. not possible"
	   STOP
	END IF

	IF (para%sig_delay > 0.0 .AND. para%N_weight > 0) THEN
	   PRINT*, "sig_delay>0 and N_weight>0 option not possible"
	   STOP
	END IF


	IF (para%sig_delay /= 0.0 .AND. para%ratio > 0.0) THEN
	   WRITE (UNIT=*, FMT=*) "for sig_delay /= 0 , ratio must <= zero"
	   STOP
	END IF

	 IF (.NOT. (para%kramer == 0.0 .OR. para%kramer == 1.0)
     &     .AND. para%IMF_option >= 2) THEN
	   WRITE (UNIT=*, FMT=*) "Kramers can only be non zero for 
     &	  IMF_option=1"
	   STOP
	END IF


	IF (para%ratio > 0.0 .AND. para%Z_imf_min >= 6) THEN
	    WRITE (UNIT=*, FMT=*) "Deformed transmission coefficients not 
     &	available"
	    WRITE (UNIT=*, FMT=*) "for  B and heavier fragments, 
     &	set z_imf_min=5"
	    STOP
	END IF

	IF (para%ratio > 0.0 .AND. para%t_delay == 0.0) THEN
	    WRITE (UNIT=*, FMT=*) "option ratio>0.0 and t_delay>0.0"
	    WRITE (UNIT=*, FMT=*) "does not compute"
	    STOP
	END IF

	RETURN
	END SUBROUTINE input

	    SUBROUTINE crash()
 ! write out event information when a error has been detected
	    INTEGER :: i

	 OPEN (UNIT=78, FILE="gem.crash", ACTION="WRITE", STATUS="NEW")       

	    DO i=1,n_frag
        WRITE (UNIT=78, FMT=*) "fragment",i
        WRITE (UNIT=78, FMT=*) "z,a=",frag(i)%Z,frag(i)%A
        WRITE (UNIT=78, FMT=*) "id,level=",frag(i)%id,frag(i)%level
        WRITE (UNIT=78, FMT=*) "Ex,time=",frag(i)%ex,frag(i)%time
        WRITE (UNIT=78, FMT=*) "spin(j,theta,phi)=",
     &       frag(i)%j,frag(i)%spin%theta,frag(i)%spin%phi
        WRITE (UNIT=78, FMT=*) "J,M=",frag(i)%j,frag(i)%M
        WRITE (UNIT=78, FMT=*) "vel(v,theta,phi)="
     &       ,frag(i)%vel%v,frag(i)%vel%theta,frag(i)%vel%phi

        WRITE (UNIT=78, FMT=*) 

	    END DO 
	    WRITE (UNIT=*, FMT=*) "see gem.crash"
	    STOP 
	    END SUBROUTINE crash

	END MODULE gemsub
!------------------------------end:   MODULE gemsub-----------------------------------
