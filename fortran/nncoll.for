cc	hadron-hadron collisions
c///////////////////////////////////////////////////////////////////////////////
     	  subroutine prodrsn(prx,pry,prz,ppx,ppy,ppz,srt,i1,i2,sels,stot,
     &                     kdp,icm,kst,xmd)
c******* production of resonances delta(1232), N(1440), N(1535)
c******* delta, idr=1; N(1440), idr=2; N(1535), idr=3
c******* sels -- total elastic cross section, sins -- total inelastic c.s in NN collisions
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00
        common/ee/mseed,nseed
        COMMON/SZPI/S10(600),S11(600),S01(600),S10D(600)
        COMMON/SZWP/W10(600),W11(600),W01(600),W10D(600)
	  common/dmss/XDMASS(500),XNM144(500),XNM153(500)

	kdp=0
	kst=0
	pi=3.141592654
	RMASS = 0.938
	EMASS = 0.547
	DLMAS = 1.232
	DMMIN = 1.077
	DMMAX = srt-RMASS
	skprd = 5.     !mb for kaon
	DSMTH = 1.1
	NS = 1000.
	DM = (DMMAX-DMMIN)/NS

	ppx=0.
	ppy=0.
	ppz=0.

	is1=iso(i1)
	is2=iso(i2)
	iy1=isy(i1)
	iy2=isy(i2)
c....ietanq gives different parametrizations of SIGQ0 in relcol
      IETANQ=1
c      SIGNQ -> N(1535)-CREATION
c      SIGNPI-> S-STATE PION PRODUCTION
c --- CROSS SECTION FOR N(1535) PRODUCTION
      CUTQ = 2.*RMASS + EMASS
      SRT2 = SRT - CUTQ
      SIGQ0= 0.0
      IF(SRT2 .GT. 0.0) THEN
        IF(IETANQ .EQ. 0) SIGQ0= 0.34 * SRT2 / (0.253 + SRT2**2)
        IF(IETANQ .EQ. 1) SIGQ0= 0.4  * SRT2 / (0.552 + SRT2**2)
        IF(IETANQ .EQ. 2) SIGQ0= 0.204* SRT2 / (0.058 + SRT2**2)
      END IF

        LL = MAX1(500.*(srt - 2.015), 0.) + 1
        LL = MIN0(LL,500)
	  SSDP=S10D(LL)

	IF(iso(i1).eq.iso(i2)) THEN
	 sins=1.33333*xmd*sigmand(srt)+sigmans(srt,iso(i1),iso(i2))+
     &      1.66667*sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2)
       SIGNQ = SIGQ0
       SIGNPI= SSDP
	ELSE
	 sins=0.66667*xmd*sigmand(srt)+sigmans(srt,iso(i1),iso(i2))+
     &      sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2)
       SIGNQ = 3.0 * SIGQ0
       SIGNPI= 0.5 * SSDP
	ENDIF
	 srem=sels+sins+SIGNQ+SIGNPI+skprd
	 rsnq=(sels+sins+SIGNQ)/stot
	 rspi=(sels+sins+SIGNQ+SIGNPI)/stot
       x1=RANFL(mseed)
       if(x1 .GT. srem/stot) return
        pr=sqrt(prx*prx+pry*pry+prz*prz)
        c2=prz/pr
        rsig=sels/stot
        if(x1.gt.rsig) goto 3

        as=(3.65*(srt-1.866))**6  
        a=6.*as/(1.+as) 
        ta=-2.*pr*pr 
      x2=RANFL(mseed) 
        t1=alog((1-x2)*XEXP(a*ta)+x2)/a
        c1=1.-t1/ta
        t1=2.*pi*RANFL(nseed)
        goto 4

  3     if(srt.lt.2.015)  return 
                  
!!        em1=0.75*(srt-2.015)+1.077
!	  BREIT-WIGNER DISTRIBUTION
c  	  if(DMMAX .LE. DMMIN)  RETURN    
c 33	  DLMTS = DMMIN + (DMMAX-DMMIN)*RANFL(nseed)
c	  XM = (DLMTS-DLMAS)/0.115
c	  PRDL =1./(1.+4.*XM**2)
c	  IF(PRDL .LT. RANFL(mseed))  GOTO 33
c	  IF(DLMTS .LT. DSMTH)       THEN             !smearing delta mass
c	  PSM=(DLMTS-DMMIN)/(DSMTH-DMMIN)	
c	  IF(PSM .LT. RANFL(mseed))    GOTO 33
c	  ENDIF

        SS = srt*srt
	  PRDLmax = 0.
	  DO  I=1,NS
	  DLMTS = DMMIN + DM*I
        pr2   = (SS+DLMTS*DLMTS-RMASS*RMASS)**2/(4.*SS)-DLMTS**2
	  if(pr2 .LE. 0.)   goto 33
	  pr = sqrt(pr2)
	  XM    = 0.46*DLMAS**2/((DLMTS**2-DLMAS**2)**2+0.115**2*DLMAS**2)
	  PRDL  = pr*XM*DLMTS
	  IF(PRDL .GT. PRDLmax)  PRDLmax=PRDL
	  ENDDO
 33	  DLMTS = DMMIN + (DMMAX-DMMIN)*RANFL(nseed)
        pr2   = (SS+DLMTS*DLMTS-RMASS*RMASS)**2/(4.*SS)-DLMTS**2
	  if(pr2 .LE. 0.)   goto 33
	  pr = sqrt(pr2)
	  XM    = 0.46*DLMAS**2/((DLMTS**2-DLMAS**2)**2+0.115**2*DLMAS**2)
	  PRDL  = pr*XM*DLMTS
	  RS    = PRDL/PRDLmax
	  IF(RS .LT. RANFL(mseed))   GOTO 33
	  IF(DLMTS .LT. DSMTH)       THEN             !smearing delta mass
	  PSM=(DLMTS-DMMIN)/(DSMTH-DMMIN)	
	  IF(PSM .LT. RANFL(mseed))  GOTO 33
	  ENDIF

	  em1 = DLMTS
	  em2 = 1.44
	  em3 = 1.535
	  IF((SRT.GT.2.18) .AND. (SRT.LT.2.68))  THEN
 44	  DLMTS = DMMIN + 0.2 + (DMMAX-DMMIN)*RANFL(nseed)
	  XM = (DLMTS-em2)/0.2
	  PRDL =1./(1.+4.*XM**2)
	  IF(PRDL .LT. RANFL(mseed))  GOTO 44
	  ENDIF
	  em2 = DLMTS    
	  IF((SRT.GT.2.32) .AND. (SRT.LT.2.68)) THEN
 55	  DLMTS = DMMIN + 0.3 + (DMMAX-DMMIN)*RANFL(nseed)
	  XM = (DLMTS-em3)/0.15
	  PRDL =1./(1.+4.*XM**2)
	  IF(PRDL .LT. RANFL(mseed))  GOTO 55
	  ENDIF
	  em3 = DLMTS

	IF(iso(i1).eq.iso(i2)) THEN
	 rsignd=(sels+1.33333*xmd*sigmand(srt))/stot
	 rsigns=(sels+1.33333*xmd*sigmand(srt) +
     &         sigmans(srt,iso(i1),iso(i2)))/stot
	 rsigdd=(sels+1.33333*xmd*sigmand(srt) +sigmans(srt,iso(i1),
     &  iso(i2))+1.66667*sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2))/stot

        if((x1.gt.rsig) .and. (x1.le.rsignd)) then
	if(iso(i1).eq.1) then
	if(RANFL(nseed).lt.0.75) then
	em(i1)=  em1
	idr(i1)= 1
	iso(i1)= 2
	isy(i1)= 0
	em(i2)=  0.938
	idr(i2)= 0
	iso(i2)= 0
	isy(i2)= 1
	else
	em(i1)=  em1
	idr(i1)= 1
	iso(i1)= 1
	isy(i1)= 0
	em(i2)=  0.938
	idr(i2)= 0
	iso(i2)= 1
	isy(i2)= -1
	endif
	elseif(iso(i1).eq.0) then
	if(RANFL(nseed).lt.0.75) then
	em(i1)=  em1
	idr(i1)= 1
	iso(i1)= -1
	isy(i1)= 0
	em(i2)=  0.938
	idr(i2)= 0
	iso(i2)= 1
	isy(i2)= -1
	else
	em(i1)=  em1
	idr(i1)= 1
	iso(i1)= 0
	isy(i1)= 0
	em(i2)=  0.938
	idr(i2)= 0
	iso(i2)= 0
	isy(i2)= 1
	endif
	endif
	  elseif((x1.gt.rsignd).and.(x1.le.rsigns)) then
	em(i1)=  em2
	idr(i1)= 2
	isy(i1)= 0
	  elseif((x1.gt.rsigns).and.(x1.le.rsigdd)) then
	em(i1)=  em1
	idr(i1)= 1
	em(i2)=  em1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
 	   if(RANFL(nseed).lt.0.6) then
	    if(iso(i1).eq.1) then
	iso(i1)= 0
	iso(i2)= 2
	    elseif(iso(i1).eq.0) then
	iso(i1)= 1
	iso(i2)= -1
	    endif
	   else
	    if(iso(i1).eq.1) then
	iso(i1)= 1
	iso(i2)= 1
	    elseif(iso(i1).eq.0) then
	iso(i1)= 0
	iso(i2)= 0
	    endif
	   endif
	  elseif((x1.gt.rsigdd).and.(x1.le.rsnq)) then
	em(i1)=  em3
	idr(i1)= 3
	isy(i1)= 0
	  elseif((x1.gt.rsnq).and.(x1.le.rspi)) then
          icm=0
 10       icm=icm+1
          IF (idms(icm).NE.0) GOTO 10
          IF (icm.GE.np) THEN
          WRITE(*,'(10X,''***  TOO MANY PIONS IN prodrsn  ***'')')
          ENDIF
	idms(icm)=1
      IF ((iso(i1)+iso(i2)).EQ.1) icp(icm)= 0
      IF ((iso(i1)+iso(i2)).EQ.0) icp(icm)= -1
      IF ((iso(i1)+iso(i2)).EQ.2) icp(icm)= 1
      IF ((icp(icm).eq.1) .or. (icp(icm).eq.-1)) emson(icm)=.13957 
	IF (icp(icm).eq.0) emson(icm)=.13496
      IF ((iso(i1)+iso(i2)).NE.1) THEN
         IF(iso(i1).EQ.0) iso(i2)= 1
         IF(iso(i1).EQ.1) iso(i2)= 0
      ENDIF
	if(iso(i2).eq.1) isy(i2)= -1
	if(iso(i2).eq.0) isy(i2)= 1
	  goto 50
	  else
	 kst=1
	 kssb= 0
	 call prkbyby(srt,i1,i2,skprd,kssb)
	 if(kssb.ne.0) return
	 call prbbkkb(srt,i1,i2,skprd)
	 return
	  endif

	ELSE

	rsignd=(sels+0.66667*xmd*sigmand(srt))/stot
	rsigns=(sels+0.66667*xmd*sigmand(srt)+sigmans(srt,iso(i1),
     &        iso(i2)))/stot
	rsigdd=(sels+0.66667*xmd*sigmand(srt)+sigmans(srt,iso(i1),
     &        iso(i2))+sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2))/stot
	xdd=sdd1/sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2)

        if((x1.gt.rsig).and.(x1.le.rsignd)) then
	    if(RANFL(nseed).lt.0.5) then
	if(iso(i1).eq.1) then
	em(i1)=  em1
	idr(i1)= 1
	isy(i1)= 0
	else
	em(i2)=  em1
	idr(i2)= 1
	isy(i2)= 0
	endif
	    else
	if(iso(i1).eq.0) then
	em(i1)=  em1
	idr(i1)= 1
	isy(i1)= 0
	else
	em(i2)=  em1
	idr(i2)= 1
	isy(i2)= 0
	endif
	    endif
	  elseif((x1.gt.rsignd).and.(x1.le.rsigns)) then
	    if(RANFL(nseed).lt.0.5) then
	if(iso(i1).eq.1) then
	em(i1)=  em2
	idr(i1)= 2
	isy(i1)= 0
	else
	em(i2)=  em2
	idr(i2)= 2
	isy(i2)= 0
	endif
	    else
	if(iso(i1).eq.0) then
	em(i1)=  em2
	idr(i1)= 2
	isy(i1)= 0
	else
	em(i2)=  em2
	idr(i2)= 2
	isy(i2)= 0
	endif
	    endif
	  elseif((x1.gt.rsigns).and.(x1.le.rsigdd)) then
	em(i1)=  em1
	idr(i1)= 1
	em(i2)=  em1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
	   if(RANFL(nseed).lt.xdd) then
	iso(i1)= 0
	iso(i2)= 1
	   else
	iso(i1)= -1
	iso(i2)= 2
	   endif
	  elseif((x1.gt.rsigdd).and.(x1.le.rsnq)) then
	em(i1)=  em3
	idr(i1)= 3
	isy(i1)= 0
  	  elseif((x1.gt.rsnq).and.(x1.le.rspi)) then
          icm=0
 20       icm=icm+1
          IF (idms(icm).NE.0) GOTO 20
          IF (icm.GE.np) THEN
          WRITE(*,'(10X,''***  TOO MANY PIONS IN prodrsn  ***'')')
          ENDIF
	idms(icm)=1
      IF ((iso(i1)+iso(i2)).EQ.1) icp(icm)= 0
      IF ((iso(i1)+iso(i2)).EQ.0) icp(icm)= -1
      IF ((iso(i1)+iso(i2)).EQ.2) icp(icm)= 1
      IF ((icp(icm).eq.1) .or. (icp(icm).eq.-1)) emson(icm)=.13957 
	IF (icp(icm).eq.0) emson(icm)=.13496
      IF ((iso(i1)+iso(i2)).NE.1) THEN
         IF(iso(i1).EQ.0) iso(i2)= 1
         IF(iso(i1).EQ.1) iso(i2)= 0
      ENDIF
	if(iso(i2).eq.1) isy(i2)= -1
	if(iso(i2).eq.0) isy(i2)= 1
	  goto 50
	  else
	 kst=1
	 kssb= 0
	 call prkbyby(srt,i1,i2,skprd,kssb)
	 if(kssb.ne.0) return
	 call prbbkkb(srt,i1,i2,skprd)
	 return
	  endif
	ENDIF

!******************************************************************
        if((x1.gt.rsig) .and. (x1.le.rsignd))  then
	  IM = INT((em1-DMMIN)/0.01) + 1
	  IF(IM.LE.200 .AND. IM.GE.1)	XDMASS(IM) = XDMASS(IM) + 1
	  elseif((x1.gt.rsignd).and.(x1.le.rsigns)) then
	  IM = INT((em2-DMMIN)/0.01) + 1
	  IF(IM.LE.300 .AND. IM.GE.1)	XNM144(IM) = XNM144(IM) + 1
	  elseif((x1.gt.rsigdd).and.(x1.le.rsnq)) then
	  IM = INT((em3-DMMIN)/0.01) + 1
	  IF(IM.LE.300 .AND. IM.GE.1)	XNM153(IM) = XNM153(IM) + 1
	  endif
!******************************************************************

	 SS=srt*srt
       pr2=(SS+em(i1)*em(i1)-em(i2)*em(i2))**2/(4.*SS)-em(i1)**2
	 if(pr2.LE.0.) then
	 iso(i1)=is1
	 iso(i2)=is2
	 isy(i1)=iy1
	 isy(i2)=iy2
	 idr(i1)=0
	 idr(i2)=0
	 em(i1)=  0.938
	 em(i2)=  0.938
	 return
	 endif

	 pr=sqrt(pr2)
       c1=1.-2.*RANFL(mseed)
       t1=2.*pi*RANFL(mseed)

 4      s2=sqrt(1.-c2*c2)
        if(prx.eq.0. .or. pry.eq.0.) goto 5
        t2=atan2(pry,prx)
 5     if(c1.gt.1)  c1=1.
	 if(c1.lt.-1) c1=-1.
	  s1=sqrt(1.-c1*c1)
        ct1=cos(t1)
        st1=sin(t1)
        ct2=cos(t2)
        st2=sin(t2)

        prz=pr*(c1*c2-s1*s2*ct1)
        ss =c2*s1*ct1+s2*c1 
        prx=pr*(ss*ct2-s1*st1*st2)
        pry=pr*(ss*st2+s1*st1*ct2)

	  return

* --- SET THE NEW MOMENTUM COORDINATES IN S-STATE PRODUCTION CASE
50      PMASS=emson(icm)
* ----- MAXIMAL PIONENERGY
        EMM=(srt**2-4.*RMASS**2+PMASS**2)/2./srt
        ENN=PMASS+(EMM-PMASS)*RANFL(nseed)
	  IF((ENN*ENN-PMASS**2).LE.0.) THEN
	  idms(icm)=0
	  icp(icm) =0
	  emson(icm)=0.
	  RETURN
	  ENDIF
        OM=SQRT(ENN*ENN-PMASS**2)
        DELTA=srt-ENN
* ----- PION ANGLE WINK IN CM  - FRAME
        ANG=2.*RANFL(nseed)-1.
        WINK=ACOS(ANG)
* ----  NUCLEON ANGLES  THET AND PHI
        THET=ACOS(2.*RANFL(nseed)-1.)
        PHI=2.*PI*RANFL(nseed)
        pxms(icm)=OM*SIN(WINK)
        pyms(icm)=0.0
        pzms(icm)=OM*COS(WINK)
        xms(icm)=0.5*(x(i1)+x(i2))
        yms(icm)=0.5*(y(i1)+y(i2))
        zms(icm)=0.5*(z(i1)+z(i2))
        SQ=DELTA**2*(DELTA**2-4.*RMASS**2-OM**2)
        SQ=SQ/(DELTA**2-(OM*COS(THET))**2)
        QQ=SQ
	  IF(SQ.LE.0.) THEN
	  idms(icm)=0
	  icp(icm) =0
	  emson(icm)=0.
	  RETURN
	  ENDIF
        SQ=SQRT(SQ)
        QA=THET
        QG=PHI
        QT=WINK
        QN1=-(COS(QT)*SIN(QA)*COS(QG)+SIN(QT)*COS(QA))
        QN2= SIN(QA)*SIN(QG)
        QN3=-(-SIN(QT)*SIN(QA)*COS(QG)+COS(QT)*COS(QA))
* ----- LENGHTH OF THE RELATIV NUCLEON VECTOR K
	prx = -pxms(icm)/2.+QN1*SQ/2.
	pry = -pyms(icm)/2.+QN2*SQ/2.
	prz = -pzms(icm)/2.+QN3*SQ/2.
      ppx = -pxms(icm)/2.-QN1*SQ/2.
      ppy = -pyms(icm)/2.-QN2*SQ/2.
      ppz = -pzms(icm)/2.-QN3*SQ/2.

  	  knn=1
  	  kdp=1
        return
        end
c///////////////////////////////////////////////////////////////////////////////
        	subroutine collres(srt,i1,i2,snnel,xmd)
c***** collisions of resonance-resonance or resonance-nucleon
c***** sels --- elastic cross section for  NR->NR
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00 
        common/ee/mseed,nseed
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/densi1/denstl(500),dens1(500),velcm(3)

	xpi=31.415926
	pi= 3.1415926
	sgnd= 6.           !mb
	skprd=5.

	RMASS = 0.938
	EMASS = 0.547
      DML = 1.115
      DMS = 1.189
      DMN = 0.938
	DMD = 1.232
	DN1 = 1.44
	DN2 = 1.535
c....ietanq gives different parametrizations of SIGQ0 in relcol
      IETANQ=1

	 is1=iso(i1)
	 is2=iso(i2)
	 id1=idr(i1)
	 id2=idr(i2)
	 if(idr(i1).GT.idr(i2)) then
	 is1=iso(i2)
	 is2=iso(i1)
	 id1=idr(i2)
	 id2=idr(i1)
	 endif

	  if(srt.lt.2.015) return
        if(id1.eq.0 .and. id2.eq.0) return

        CUTQ = 2. * RMASS + EMASS
        SRT2 = srt - CUTQ
        SIGQ0= 0.
        IF(SRT2 .GT. 0.0) THEN
          IF(IETANQ .EQ. 0) SIGQ0= 0.34 * SRT2 / (0.253 + SRT2**2)
          IF(IETANQ .EQ. 1) SIGQ0= 0.4  * SRT2 / (0.552 + SRT2**2)
          IF(IETANQ .EQ. 2) SIGQ0= 0.204* SRT2 / (0.058 + SRT2**2)
        END IF

        x1=x(i1)                                                                
        y1=y(i1)                                                                
        z1=z(i1)
        x2=x(i2)                                                                
        y2=y(i2)                                                                
        z2=z(i2)
        r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
        pz1=pz(i1)                                                              
        pz2=pz(i2)                                                              
        px1=px(i1)                                                              
        px2=px(i2)                                                              
        py1=py(i1)                                                              
        py2=py(i2)   
	  em1=em(i1)                                                              
        em2=em(i2)                                                                                              
        e1=sqrt(em(i1)*em(i1)+pz1*pz1+px1*px1+py1*py1) 
        e2=sqrt(em(i2)*em(i2)+pz2*pz2+px2*px2+py2*py2)
        bx=velcm(1)
        by=velcm(2)
        bz=velcm(3)
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
	  prf=sqrt(srt*srt/4.-.938*.938)
   	  pr2=((srt**2+em(i1)**2-em(i2)**2)/(2.*srt))**2-em(i1)**2
	  if(pr2.le.0.) return
	  cfi=prf*prf/pr2/4.
	  cfpw=em2*prf*prf/sqrt(pr2)/4.

	  cdabs = 4.*fdb(srt)               !fdb-0.065 modification on absorption of delta NPA552(1993)549
!!	  cdabs = 2.                        !xmd
!!	  cdpaw = fpw(srt)

* ------- CALCULATION OF THE CROSS SECTION FOR NORMAL DETAILED BALANCE
         IF((is1+is2) .EQ. 1) THEN
         SIGQN=3.0*SIGQ0
         ELSE
         SIGQN=SIGQ0
         END IF

	  prfdn2=((srt**2+DMN**2-DMD**2)/(2.*srt))**2-DMN**2
	  cfdn=prfdn2/pr2
	  if(prfdn2.LE.0.) cfdn=0.
	  prfnn1=((srt**2+DMN**2-DN1**2)/(2.*srt))**2-DMN**2
	  cfnn1=prfnn1/pr2
	  if(prfnn1.LE.0.) cfnn1=0.
	  prfnn2=((srt**2+DMN**2-DN2**2)/(2.*srt))**2-DMN**2
	  cfnn2=prfnn2/pr2
	  if(prfnn2.LE.0.) cfnn2=0.

      if(id1.eq.0 .and. id2.eq.1) then
	if(is1.eq.0) then                                         !ND->ND
	if(is2.eq.-1) sels= 9.*snnel/32.
	if(is2.eq.0)  sels=13.*snnel/32.
	if(is2.eq.1)  sels=17.*snnel/32.
	if(is2.eq.2)  sels=21.*snnel/32.
	else
	if(is2.eq.-1) sels=21.*snnel/32.
	if(is2.eq.0)  sels=17.*snnel/32.
	if(is2.eq.1)  sels=13.*snnel/32.
	if(is2.eq.2)  sels= 9.*snnel/32.
	endif

	  CROSdabs= cfi*sigmand(srt)/cdabs
!!	  CROSdabs= cfpw*sigmand(srt)/cdpaw   
	  if(is2.eq.2 .or. is2.eq.-1) then                         !ND->NN   sigmand(srt), sgnd
	  sdnn=CROSdabs
	  elseif(is2.eq.is1) then
	  sdnn=0.33333*CROSdabs
	  else
	  sdnn=0.66667*CROSdabs
	  endif
	  sns1=0.25*cfnn1*sigmans(srt,is1,is2)                     !ND -> NN*(1440)
	  sns2=0.25*cfnn2*SIGQN                                    !ND -> NN*(1535)
	  CROSSnds=sels+sdnn+sns1+sns2+skprd
	elseif(id1.eq.0 .and. id2.eq.2) then
	  sels= snnel                                              !NN*(1440) -> NN*(1440)
        ssnn=4.*cfi*sigmans(srt,is1,is2)                         !NN*(1440) -> NN
        ssnd=2.*cfdn*sigmans(srt,is1,is2)                        !NN*(1440) -> ND
        snn2=cfnn2*sigmans(srt,is1,is2)                          !NN*(1440) -> NN*(1535)
        CROSSnds=sels+ssnn+ssnd+snn2
	elseif(id1.eq.0 .and. id2.eq.3) then
	  sels= snnel                                              !NN*(1535) -> NN*(1535)
	  ssnn= SIGQN*prf**2/pr2                                   !NN*(1535) -> NN
        ssnd=2.*cfdn*SIGQN                                       !NN*(1535) -> ND
        snn1=cfnn1*SIGQN                                         !NN*(1535) -> NN*(1440)
        CROSSnds=sels+ssnn+ssnd+snn1
	elseif(id1.eq.0 .and. id2.eq.4) then                       !LN -> LN, SN
        TTX=px1-px2
        TTY=py1-py2
        TTZ=pz1-pz2
        PKB=SQRT(TTX**2+TTY**2+TTZ**2)
	  CROSSnds=CHYPN(srt,PKB,1)+CHYPN(srt,PKB,2)
	elseif(id1.eq.0 .and. id2.eq.5) then                       !SN -> LN, SN
        TTX=px1-px2
        TTY=py1-py2
        TTZ=pz1-pz2
        PKB=SQRT(TTX**2+TTY**2+TTZ**2)
	  CROSSnds=CHYPN(srt,PKB,3)+CHYPN(srt,PKB,4)
	elseif(id1.eq.1 .and. id2.eq.1) then
	  csdd=cfi*sigddnn(srt,is1,is2)/4.                         !DD - >NN
	  css1=0.25*sigmans(srt,is1,is2)                           !DD -> NN*(1440)
	  css2=0.25*SIGQN                                          !DD -> NN*(1535)
	  CROSSnds=csdd+css1+css2+skprd
	else
	 return
	endif

	dint=CROSSnds/xpi

        b2=.0001+bx*bx+by*by+bz*bz 
        pb=bx*px1+by*py1+bz*pz1                                                 
        gg=(gam-1.)*pb/b2-gam*e1                                                
        pzr=gg*bz+pz1                                                           
        pxr=gg*bx+px1                                                           
        pyr=gg*by+py1                                                           
        rb=bx*(x1-x2)+by*(y1-y2)+bz*(z1-z2)                                     
        gg=(gam-1.)*rb/b2                                                       
        dx=gg*bx+x1-x2                                                          
        dy=gg*by+y1-y2                                                          
        dz=gg*bz+z1-z2
        dr=sqrt(dx*dx+dy*dy+dz*dz)
        pr=sqrt(pxr*pxr+pyr*pyr+pzr*pzr)
        dzz=(pxr*dx+pyr*dy+pzr*dz)/pr                                           
        c00=dr*dr-dzz*dzz                           !!  feng
	  if(c00.le.0.) c00=0.
        bb=sqrt(c00)
        if(bb.gt.sqrt(dint)) return
        n2=n2+1                                                                 
        e1r=sqrt(em(i1)*em(i1)+pr*pr)                                        
        e2r=sqrt(em(i2)*em(i2)+pr*pr)                                           
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.                                         
        if(dd.lt.abs(dzz))   return
        c2=pzr/pr
c-------------------------------------------------------
      if(id1.eq.0 .and. id2.eq.1) then
	rat0=sels/CROSSnds
	rat1=(sels+sdnn)/CROSSnds
	rat2=(sels+sdnn+sns1)/CROSSnds
	rat3=(sels+sdnn+sns1+sns2)/CROSSnds
      xx=RANFL(nseed)

	 if(xx.GT.rat0) goto 11         !ND --> ND
	  if(is1.eq.0) then
	    if(is2.eq.0 .AND. RANFL(mseed).le.0.923) then
		iso(i1)=1
	    iso(i2)=-1
	    isy(i1)=-1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	    if(is2.eq.1 .AND. RANFL(mseed).le.0.941) then
		iso(i1)=1
	    iso(i2)=0
	    isy(i1)=-1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	    if(is2.eq.2 .AND. RANFL(mseed).le.0.571) then
		iso(i1)=1
	    iso(i2)=1
	    isy(i1)=-1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	  else
	    if(is2.eq.1 .AND. RANFL(mseed).le.0.923) then
		iso(i1)=0
	    iso(i2)=2
	    isy(i1)=1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	    if(is2.eq.0 .AND. RANFL(mseed).le.0.941) then
		iso(i1)=0
	    iso(i2)=1
	    isy(i1)=1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	    if(is2.eq.-1 .AND. RANFL(mseed).le.0.571) then
		iso(i1)=0
	    iso(i2)=0
	    isy(i1)=1
	    isy(i2)=0
		idr(i1)=0
	    idr(i2)=1
          em(i1)=DMN
          em(i2)=DMD
	    endif
	  endif
	   pr=sqrt(pr2)
 	 goto 420
 11	   if((is1+is2).GT.2 .OR. (is1+is2).LT.0) then
	   if(xx.GT.rat3) goto 44
	   return
	   endif
	   if(xx.GT.rat1) goto 22
	   if(is1+is2.eq.2) then        !n+delta++ -->p+p, p+delta+ -->p+p 
	   iso(i1)=1
	   iso(i2)=1
	   isy(i1)=-1
	   isy(i2)=-1
	   elseif(is1+is2.eq.1) then    !p+delta0 -->p+n, n+delta+ -->p+n
	   if(iso(i1).eq.0) isy(i1)=1
	   if(iso(i2).eq.0) isy(i2)=1
	   if(iso(i1).eq.1) isy(i1)=-1
	   if(iso(i2).eq.1) isy(i2)=-1
	   elseif(is1+is2.eq.0) then    !n+delta0 -->n+n, p+delta- -->n+n
	   iso(i1)=0
	   iso(i2)=0
	   isy(i1)=1
	   isy(i2)=1
	   endif
	 idr(i1)=0
	 idr(i2)=0
       em(i1)=.938
       em(i2)=.938
       pr=prf
 	 goto 420
 22	 if(xx.GT.rat2) goto 33
   	  pr2=((srt**2+DN1**2-DMN**2)/(2.*srt))**2-DN1**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	   if(is1+is2 .NE. 1) then
	   if(is1.eq.1)	 iso(i1)=1
	   if(is1.eq.1)	 iso(i2)=1
	   if(is1.eq.0)	 iso(i1)=0
	   if(is1.eq.0)	 iso(i2)=0
	   endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=2
       em(i1)= DMN
       em(i2)= DN1
 	 goto 420
 33	 if(xx.GT.rat3) goto 44
        pr2=((srt**2+DN2**2-DMN**2)/(2.*srt))**2-DN2**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	   if(is1+is2 .NE. 1) then
	   if(is1.eq.1)	 iso(i1)=1
	   if(is1.eq.1)	 iso(i2)=1
	   if(is1.eq.0)	 iso(i1)=0
	   if(is1.eq.0)	 iso(i2)=0
	   endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=3
       em(i1)= DMN
       em(i2)= DN2
 	 goto 420
 44	 kssb= 0
	 call prkbyby(srt,i1,i2,skprd,kssb)                      ! kaon production by baryon-baryon collisions
	 if(kssb.ne.0) return
	 call prbbkkb(srt,i1,i2,skprd)                           ! kaon-kaonb production by baryon-baryon collisions
	 return
	elseif(id1.eq.0 .and. id2.eq.2) then
	rat0=sels/CROSSnds
	rat1=(sels+ssnn)/CROSSnds
	rat2=(sels+ssnn+ssnd)/CROSSnds
      xx=RANFL(nseed)

	 if(xx.LT.rat0) then
	   pr=sqrt(pr2)
	 elseif(xx.GE.rat0 .AND. xx.LT.rat1) then
	   if(iso(i1).eq.0) isy(i1)=1
	   if(iso(i2).eq.0) isy(i2)=1
	   if(iso(i1).eq.1) isy(i1)=-1
	   if(iso(i2).eq.1) isy(i2)=-1
	 idr(i1)=0
	 idr(i2)=0
       em(i1)=.938
       em(i2)=.938
       pr=prf
	 elseif(xx.GE.rat1 .AND. xx.LT.rat2) then
	if(prfdn2.LE.0.) return
	  if(RANFL(mseed).LE.0.5) then
	  iso(i1)=0
	  iso(i2)=is1+is2-iso(i1)
	  else
	  iso(i1)=1
	  iso(i2)=is1+is2-iso(i1)
	  endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=1
       em(i1)= DMN
       em(i2)= DMD
       pr=sqrt(prfdn2)
	 else
	if(prfnn2.LE.0.) return
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=3
       em(i1)= DMN
       em(i2)= DN2
       pr=sqrt(prfnn2)
	 endif
 	 goto 420
	elseif(id1.eq.0 .and. id2.eq.3) then
	rat0=sels/CROSSnds
	rat1=(sels+ssnn)/CROSSnds
	rat2=(sels+ssnn+ssnd)/CROSSnds
      xx=RANFL(nseed)

	 if(xx.LT.rat0) then
	   pr=sqrt(pr2)
	 elseif(xx.GE.rat0 .AND. xx.LT.rat1) then
	   if(iso(i1).eq.0) isy(i1)=1
	   if(iso(i2).eq.0) isy(i2)=1
	   if(iso(i1).eq.1) isy(i1)=-1
	   if(iso(i2).eq.1) isy(i2)=-1
	 idr(i1)=0
	 idr(i2)=0
       em(i1)=.938
       em(i2)=.938
       pr=prf
	 elseif(xx.GE.rat1 .AND. xx.LT.rat2) then
	if(prfdn2.LE.0.) return
	  if(RANFL(mseed).LE.0.5) then
	  iso(i1)=0
	  iso(i2)=is1+is2-iso(i1)
	  else
	  iso(i1)=1
	  iso(i2)=is1+is2-iso(i1)
	  endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=1
       em(i1)= DMN
       em(i2)= DMD
       pr=sqrt(prfdn2)
	 else
	if(prfnn1.LE.0.) return
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=2
       em(i1)= DMN
       em(i2)= DN1
       pr=sqrt(prfnn1)
	 endif
 	 goto 420
	elseif(id1.eq.0 .and. id2.eq.4) then
	  xx=RANFL(nseed)
	  rat=CHYPN(srt,PKB,1)/CROSSnds
	  if(xx.lt.rat) then
	  pr=sqrt(pr2)
	  else
   	  pr2=((srt**2+DMS**2-DMN**2)/(2.*srt))**2-DMS**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	 if(iso(i1)+iso(i2).eq.1) then
	 if(RANFL(mseed).lt.0.5) then
	 iso(i1)=1
	 isy(i1)=-1
	 iso(i2)=0
	 else
	 iso(i1)=0
	 isy(i1)=1
	 iso(i2)=1
	 endif
	 else
	 if(RANFL(mseed).lt.0.5) then
	 iso(i1)=0
	 isy(i1)=1
	 iso(i2)=0
	 else
	 iso(i1)=1
	 isy(i1)=-1
	 iso(i2)=-1
	 endif
	 endif
	  idr(i1)=0
	  em(i1)=DMN
	  idr(i2)=5
	  em(i2)=DMS
	  isy(i2)=0
 	  endif
	  goto 420
	elseif(id1.eq.0 .and. id2.eq.5) then
	  xx=RANFL(nseed)
	  rat=CHYPN(srt,PKB,3)/CROSSnds
	  if(xx.lt.rat) then
   	  pr2=((srt**2+DMS**2-DMN**2)/(2.*srt))**2-DMS**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	 if(iso(i1)+iso(i2).eq.2) then
	 iso(i1)=1
	 isy(i1)=-1
	 iso(i2)=1
	 elseif(iso(i1)+iso(i2).eq.1) then
	 if(RANFL(mseed).lt.0.5) then
	 iso(i1)=1
	 isy(i1)=-1
	 iso(i2)=0
	 else
	 iso(i1)=0
	 isy(i1)=1
	 iso(i2)=1
	 endif
	 elseif(iso(i1)+iso(i2).eq.0) then
	 if(RANFL(mseed).lt.0.5) then
	 iso(i1)=1
	 isy(i1)=-1
	 iso(i2)=-1
	 else
	 iso(i1)=0
	 isy(i1)=1
	 iso(i2)=0
	 endif
	 elseif(iso(i1)+iso(i2).eq.-1) then
	 iso(i1)=0
	 isy(i1)=1
	 iso(i2)=-1
	 else
	 return
	 endif
	  idr(i1)=0
	  em(i1)=DMN
	  idr(i2)=5
	  em(i2)=DMS
	  isy(i2)=0
	  else
   	  pr2=((srt**2+DML**2-DMN**2)/(2.*srt))**2-DML**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	 if(iso(i1)+iso(i2).eq.1) then
	 iso(i1)=1
	 isy(i1)=-1
	 elseif(iso(i1)+iso(i2).eq.0) then
	 iso(i1)=0
	 isy(i1)=1
	 else
	 return
	 endif
	  idr(i1)=0
	  em(i1)=DMN
	  iso(i2)=0
	  idr(i2)=4
	  em(i2)=DML
	  isy(i2)=0
 	  endif
	  goto 420
	elseif(id1.eq.1 .and. id2.eq.1) then
	rat1=csdd/CROSSnds
	rat2=(csdd+css1)/CROSSnds
	rat3=(csdd+css1+css2)/CROSSnds
      xx=RANFL(nseed)
	 if((is1+is2).GT.2 .OR. (is1+is2).LT.0) then
	 if(xx.GT.rat3) goto 88
	 return
	 endif

       if(xx.GT.rat1) goto 66
        if(is1+is2.eq.2) then
	  iso(i1)=1
	  iso(i2)=1
	  isy(i1)=-1
	  isy(i2)=-1
	  elseif(is1+is2.eq.1) then
	  iso(i1)=1
	  iso(i2)=0
	  isy(i1)=-1
	  isy(i2)=1
	  elseif(is1+is2.eq.0) then
	  iso(i1)=0
	  iso(i2)=0
	  isy(i1)=1
	  isy(i2)=1
	  endif
	 idr(i1)=0
	 idr(i2)=0
       em(i1)=.938
       em(i2)=.938
       pr=prf
 	 goto 420
 66	 if(xx.GT.rat2) goto 77
   	  pr2=((srt**2+DN1**2-DMN**2)/(2.*srt))**2-DN1**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	   if(is1+is2 .EQ. 1) then
		if(RANFL(mseed) .LT. 0.5) then
	    iso(i1)=1
	    iso(i2)=0
	    else
	    iso(i1)=0
	    iso(i2)=1
	    endif
	   else
	   if(is1.eq.1)	 iso(i1)=1
	   if(is1.eq.1)	 iso(i2)=1
	   if(is1.eq.0)	 iso(i1)=0
	   if(is1.eq.0)	 iso(i2)=0
	   endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=2
       em(i1)= DMN
       em(i2)= DN1
 	 goto 420
 77  	 if(xx.GT.rat3) goto 88
        pr2=((srt**2+DN2**2-DMN**2)/(2.*srt))**2-DN2**2
	  if(pr2.le.0.) return
	  pr=sqrt(pr2)
	   if(is1+is2 .EQ. 1) then
		if(RANFL(mseed) .LT. 0.5) then
	    iso(i1)=1
	    iso(i2)=0
	    else
	    iso(i1)=0
	    iso(i2)=1
	    endif
	   else
	   if(is1.eq.1)	 iso(i1)=1
	   if(is1.eq.1)	 iso(i2)=1
	   if(is1.eq.0)	 iso(i1)=0
	   if(is1.eq.0)	 iso(i2)=0
	   endif
       if(iso(i1).eq.0) isy(i1)=1
       if(iso(i1).eq.1) isy(i1)=-1
	 isy(i2)=0
	 idr(i1)=0
	 idr(i2)=3
       em(i1)= DMN
       em(i2)= DN2
 	 goto 420
 88	 kssb= 0
	 call prkbyby(srt,i1,i2,skprd,kssb)                      ! kaon production by baryon-baryon collisions
	 if(kssb.ne.0) return
	 call prbbkkb(srt,i1,i2,skprd)                           ! kaon-kaonb production by baryon-baryon collisions
	 return
	else
      return
	endif
c---------------------------------------------------

 420  c1=1.-2.*RANFL(mseed)
      t1=2.*pi*RANFL(mseed)        
      s2=sqrt(1.-c2*c2)
        if(pxr.eq.0. .or. pyr.eq.0.) goto 430
        t2=atan2(pyr,pxr)
 430    s1=sqrt(1.-c1*c1)
        ct1=cos(t1)                                                             
        st1=sin(t1)                                                             
        ct2=cos(t2)                                                             
        st2=sin(t2)                                                             
        pzr=pr*(c1*c2-s1*s2*ct1)                                                 
        ss=c2*s1*ct1+s2*c1                                                      
        pxr=pr*(ss*ct2-s1*st1*st2)                                               
        pyr=pr*(ss*st2+s1*st1*ct2) 
        eaa=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        ebb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)
        
        CALL LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)
        CALL LORENTZ(-bx,-by,-bz,ebb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)
        px(i1)= pax 
        py(i1)= pay 
        pz(i1)= paz 
        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz

    	return
     	end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine prkbyby(srt,i1,i2,sref,kssb)
c****** baryon-baryon collisions for the production of kaon 
c****** NN --> NYK, NN --> deltaYK, Ndelta --> NYK, Ndelta --> deltaYK, delta delta --> deltaYK
c****** Y=lambda, sigma+, 0, -, K=0, +.
c****** i1 nucleon or delta, i2 hyperon
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00
        common/densi1/denstl(500),dens1(500),velcm(3)
        common/ee/mseed,nseed
        common/medium/mdel,mdin,nsym,npop,nkaon

	kssb=0
	if(sref.LE.0.) return

	PI=3.141592654
      DMK0=0.498               !GeV
      DMK= 0.494
      DML= 1.115
      DMS= 1.189
      DMN= 0.938
	DMD= 1.232

!	in-medium correction on the threshold energies
	  if(nkaon.EQ.0) then
	  ak= 0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3   Vk= 27.8 MeV (rho0)
	  bk= 0.333        !GeVfm^3
	  endif

	  xp= (x(i1)+x(i2))/2.
	  yp= (y(i1)+y(i2))/2.
	  zp= (z(i1)+z(i2))/2.
	  rho = densit(xp,yp,zp)
	  xkmed = sqrt(DMK**2-ak*denscal(rho)+(bk*rho)**2) + bk*rho
	  xk0med= sqrt(DMK0**2-ak*denscal(rho)+(bk*rho)**2)+ bk*rho
	  DMK0= xk0med
	  DMK=  xkmed
c***********************************************************************

	SNLK0=DMN+DML+DMK0
	SNLK= DMN+DML+DMK
	SNSK0=DMN+DMS+DMK0
	SNSK= DMN+DMS+DMK
	SDLK0=DMD+DML+DMK0
	SDLK= DMD+DML+DMK
	SDSK0=DMD+DMS+DMK0
	SDSK= DMD+DMS+DMK

	XMK=0.
	XMR=0.
	XMH=0.

	ID1=idr(i1)
	ID2=idr(i2)
	IS1=iso(i1)
	IS2=iso(i2)
	icm= 0
 88	icm=icm+1
	IF(idms(icm).NE.0) goto 88
      IF(icm.GT.np) STOP 'TOO MANY MESONS IN prkbyby!'
	xx=RANFL(nseed)
	rati=0.

	IF(ID1.GE.2 .OR. ID2.GE.2) RETURN
	IF(ID1+ID2 .GE. 1) GOTO 18

c****** NN --> NLK0
	IF(srt.LE.SNLK0) GOTO 11
	sig=signlk(srt,IS1,IS2,0,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 4
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK0
	XMR=DMN
	XMH=DML
   	ENDIF
c****** NN --> NLK+
 11  	IF(srt.LE.SNLK) GOTO 12
	sig=signlk(srt,IS1,IS2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 4
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK
	XMR=DMN
	XMH=DML
	ENDIF
c****** NN --> NSK0
 12	rati=ratf
	IF(srt.LE.SNSK0) GOTO 13
	sig=signsk(srt,IS1,IS2,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 5
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK0
	XMR=DMN
	XMH=DMS
	ENDIF
c****** NN --> NSK+
 13	rati=ratf
	IF(srt.LE.SNSK) GOTO 14
	sig=signsk(srt,IS1,IS2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 5
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK
	XMR=DMN
	XMH=DMS
	ENDIF
c****** NN --> DLK0
 14	rati=ratf
	IF(srt.LE.SDLK0) GOTO 15
	sig=sigdlk(srt,IS1,IS2,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK0
	XMR=DMD
	XMH=DML
 	ENDIF
c****** NN --> DLK+
 15	rati=ratf
	IF(srt.LE.SDLK) GOTO 16
	sig=sigdlk(srt,IS1,IS2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK
	XMR=DMD
	XMH=DML
 	ENDIF
c****** NN --> DSK0
 16	rati=ratf
	IF(srt.LE.SDSK0) GOTO 17
	sig=sigdsk(srt,IS1,IS2,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK0
	XMR=DMD
	XMH=DMS
 	ENDIF
c****** NN --> DSK+
 17	rati=ratf
	IF(srt.LE.SDSK) GOTO 121
	sig=sigdsk(srt,IS1,IS2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK
	XMR=DMD
	XMH=DMS
 	ENDIF
	GOTO 121

 18	CONTINUE
	IF(ID1+ID2 .GE. 2) GOTO 28

	kn=IS1
	kd=IS2
	IF(ID1.EQ.1) THEN
	kn=IS2
	kd=IS1
	ENDIF
c****** ND --> NLK0
	IF(srt.LE.SNLK0) GOTO 21
	sig=sndnlk(srt,kn,kd,0,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 4
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK0
	XMR=DMN
	XMH=DML
   	ENDIF
c****** ND --> NLK+
 21  	IF(srt.LE.SNLK) GOTO 22
	sig=sndnlk(srt,kn,kd,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 4
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK
	XMR=DMN
	XMH=DML
	ENDIF
c****** ND --> NSK0
 22	rati=ratf
	IF(srt.LE.SNSK0) GOTO 23
	sig=sndnsk(srt,kn,kd,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 5
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK0
	XMR=DMN
	XMH=DMS
	ENDIF
c****** ND --> NSK+
 23	rati=ratf
	IF(srt.LE.SNSK) GOTO 24
	sig=sndnsk(srt,kn,kd,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 5
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	XMK=DMK
	XMR=DMN
	XMH=DMS
	ENDIF
c****** ND --> DLK0
 24	rati=ratf
	IF(srt.LE.SDLK0) GOTO 25
	sig=snddlk(srt,kn,kd,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK0
	XMR=DMD
	XMH=DML
 	ENDIF
c****** ND --> DLK+
 25	rati=ratf
	IF(srt.LE.SDLK) GOTO 26
	sig=snddlk(srt,kn,kd,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK
	XMR=DMD
	XMH=DML
	ENDIF
c****** ND --> DSK0
 26	rati=ratf
	IF(srt.LE.SDSK0) GOTO 27
	sig=snddsk(srt,kn,kd,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK0
	XMR=DMD
	XMH=DMS
 	ENDIF
c****** ND --> DSK+
 27	rati=ratf
	IF(srt.LE.SDSK) GOTO 121
	sig=snddsk(srt,kn,kd,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK
	XMR=DMD
	XMH=DMS
 	ENDIF
	GOTO 121

 28	CONTINUE
 	IF(ID1+ID2 .GE. 3) RETURN

	kd1=IS1
	kd2=IS2
	IF(IS1.GT.IS2) THEN
 	kd1=IS2
	kd2=IS1
	ENDIF
c****** DD --> DLK0
	IF(srt.LE.SDLK0) GOTO 31
	sig=cdddlk(srt,kd1,kd2,0,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK0
	XMR=DMD
	XMH=DML
   	ENDIF
c****** DD --> DLK+
 31	IF(srt.LE.SDLK) GOTO 32
	sig=cdddlk(srt,kd1,kd2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 4
	XMK=DMK
	XMR=DMD
	XMH=DML
 	ENDIF
c****** DD --> DSK0
 32	rati=ratf
	IF(srt.LE.SDSK0) GOTO 33
	sig=cdddsk(srt,kd1,kd2,0,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK0
	XMR=DMD
	XMH=DMS
 	ENDIF
c****** DD --> DSK+
 33	rati=ratf
	IF(srt.LE.SDSK) GOTO 121
	sig=cdddsk(srt,kd1,kd2,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GT.rati .and. xx.LT.ratf) THEN
	idms(icm)=3
	icp(icm)= 1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 5
	XMK=DMK
	XMR=DMD
	XMH=DMS
 	ENDIF

 121	CONTINUE
	IF(idms(icm).EQ.0) RETURN

	PPK=(srt*srt-(XMH+XMR+XMK)**2)*(srt*srt-(XMH+XMR-XMK)**2)/
     &    (4.*srt*srt)
	IF(PPK.LE.0.) RETURN
C..generate K momentum in NN cms
      CALL EVENT(PPK,XMK,PKX,PKY,PKZ,PK0)
c..beta-factor for lorentz-boost to kaon-system:
      BETX=PKX/(srt-PK0)
      BETY=PKY/(srt-PK0)
      BETZ=PKZ/(srt-PK0)
 35    XX    = 1. - 2.*RANFL(mseed)
       YY    = 1. - 2.*RANFL(mseed)
       ZZ    = 1. - 2.*RANFL(mseed)
       RR    = SQRT( XX**2 + YY**2 + ZZ**2 )
      IF((RR .LT. 0.001) .OR. (RR .GT. 1.))             GOTO 35
      S=(srt-PK0)**2-PKX**2-PKY**2-PKZ**2
      PA2 = (S**2+XMH**4+XMR**4-2.*S*XMH**2-2*S*XMR**2-2.*XMH**2*XMR**2)
     &      /(4.*S)
      IF(PA2.LE.0) RETURN
      PA = SQRT(PA2)
      PXX= PA * XX / RR
      PYY= PA * YY / RR
      PZZ= PA * ZZ / RR
      ENUR= SQRT(XMR**2+PA2)
      EHYP= SQRT(XMH**2+PA2)
* --------- TRANSFORMATION INTO kaon-system:
c***** nucleon or delta
      CALL LORENTZ(BETX,BETY,BETZ,ENUR,PXX,PYY,PZZ,
     &     EENUR,PPX,PPY,PPZ)
c***** lambda or sigma
      CALL LORENTZ(BETX,BETY,BETZ,EHYP,-PXX,-PYY,-PZZ,
     &     EEHYP,PHYPX,PHYPY,PHYPZ)
C.........Lorentz-transformation to calculational frame
	bx= -velcm(1)
	by= -velcm(2)
	bz= -velcm(3)
      CALL LORENTZ(bx,by,bz,PK0,PKX,PKY,PKZ,P0K,PXK,PYK,PZK)   ! kaon
      CALL LORENTZ(bx,by,bz,ENUR,PPX,PPY,PPZ,U0,UX,UY,UZ)      ! nucleon or delta
      CALL LORENTZ(bx,by,bz,EHYP,PHYPX,PHYPY,PHYPZ,            ! lambda or sigma
     &     W0,WX,WY,WZ)

	 px(i1)= UX
	 py(i1)= UY
	 pz(i1)= UZ
	 em(i1)= XMR

	 px(i2)= WX
	 py(i2)= WY
	 pz(i2)= WZ
	 em(i2)= XMH

       xxx=(x(i1) + x(i2))/2.
       yyy=(y(i1) + y(i2))/2.
       zzz=(z(i1) + z(i2))/2.
	 xms(icm)=  xxx
	 yms(icm)=  yyy
	 zms(icm)=  zzz
	 pxms(icm)= PXK
	 pyms(icm)= PYK
	 pzms(icm)= PZK
	 emson(icm)=XMK

	  kssb=1
        return
        end

c///////////////////////////////////////////////////////////////////////////////
     	subroutine prkpnby(srt,i1,i2,sref,kssb)
c****** pion(eta,omga,rho)-baryon collisions for the production of kaon 
c****** PN --> YK, PD --> YK
c****** Y=lambda, sigma+, 0, -, K=0, +.
c****** i1=pion; i2=N,D
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00 
        common/densi1/denstl(500),dens1(500),velcm(3)
        common/ee/mseed,nseed
        common/medium/mdel,mdin,nsym,npop,nkaon

	kssb=0
	if(sref.LE.0.) return

      DMK0=0.498
      DMK= 0.494
      DMN= 0.938
      DML= 1.115
      DMS= 1.189
	DMD= 1.232

!	in-medium correction on the threshold energies
	  if(nkaon.EQ.0) then
	  ak= 0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3   Vk= 27.8 MeV (rho0)
	  bk= 0.333        !GeVfm^3
	  endif

	  rho = densit(x(i2),y(i2),z(i2))
	  xkmed = sqrt(max(0.,DMK**2-ak*denscal(rho)+(bk*rho)**2)) +bk*rho
	  xk0med= sqrt(max(0.,DMK0**2-ak*denscal(rho)+(bk*rho)**2))+bk*rho
        DMK0=xk0med
        DMK= xkmed
c***********************************************************************
	SLK0= DML+DMK0
	SLK=  DML+DMK
	SSK0= DMS+DMK0
	SSK=  DMS+DMK

	ID1=idms(i1)
	ID2=idr(i2)
	IS1=icp(i1)
	IS2=iso(i2)
	XMK=0.
	XMH=0.
	xx=RANFL(nseed)
	rati=0.

	IF(ID1 .EQ. 2) GOTO 30    !etaN --> YK
	IF(ID1 .EQ. 8) GOTO 40    !rhoN --> YK
	IF(ID1 .EQ. 9) GOTO 30    !omegaN --> YK
	IF(ID1 .GE. 2) RETURN
	IF(ID2 .GE. 1) GOTO 15
c****** PN --> LK0
	IF(srt.LE.SLK0) GOTO 11
	sig=spnlk(srt,IS1,IS2,0)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK0
	XMH=DML
   	ENDIF
c****** PN --> LK+
 11  	IF(srt.LE.SLK) GOTO 12
	sig=spnlk(srt,IS1,IS2,1)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK
	XMH=DML
	ENDIF
c****** PN --> SK0
 12	rati=ratf
	IF(srt.LE.SSK0) GOTO 13
	sig=spnsk(srt,IS1,IS2,0,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	ENDIF
c****** PN --> SK+
 13	rati=ratf
	IF(srt.LE.SSK) GOTO 25
	sig=spnsk(srt,IS1,IS2,1,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK
	XMH=DMS
	ENDIF
	GOTO 25

 15	CONTINUE
 	IF(ID2.GE.2) RETURN
c****** PD --> LK0
	IF(srt.LE.SLK0) GOTO 21
	sig=cpdlk(srt,IS1,IS2,0)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK0
	XMH=DML
   	ENDIF
c****** PD --> LK+
 21  	IF(srt.LE.SLK) GOTO 22
	sig=cpdlk(srt,IS1,IS2,1)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK
	XMH=DML
	ENDIF
c****** PD --> SK0
 22	rati=ratf
	IF(srt.LE.SSK0) GOTO 23
	sig=cpdsk(srt,IS1,IS2,0,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	ENDIF
c****** PD --> SK+
 23	rati=ratf
	IF(srt.LE.SSK) GOTO 25
	sig=cpdsk(srt,IS1,IS2,1,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK
	XMH=DMS
	ENDIF
	GOTO 25

  30	CONTINUE
  	IF(ID2 .GE. 1) RETURN
	XMP=0.138
	IF(ID1 .EQ. 2) XMM=0.548
	IF(ID1 .EQ. 9) XMM=0.783

	S=srt**2
      P2PN = (S**2+DMN**4+XMP**4-2.*S*DMN**2-2*S*XMP**2-
     &        2.*DMN**2*XMP**2)/(4.*S)
      P2MM = (S**2+DMN**4+XMM**4-2.*S*DMN**2-2*S*XMM**2-
     &        2.*DMN**2*XMM**2)/(4.*S)

      IF(P2PN .LE. 0) RETURN
      IF(P2MM .LE. 0) RETURN
      PPN = SQRT(P2PN)
      PMM = SQRT(P2MM)

c****** etaN --> LK, omegaN --> LK
	IF(IS2 .EQ. 0) THEN
	IF(srt.LE.SLK0)  GOTO 31
	sig=spnlk(srt,0,0,0)*PPN/PMM
	rati=sig/sref
	  IF(xx.LT.rati) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK0
	XMH=DML
	  ENDIF
	ELSE
	IF(srt.LE.SLK)  GOTO 31
	sig=spnlk(srt,0,1,1)*PPN/PMM
	rati=sig/sref
	  IF(xx.LT.rati) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK
	XMH=DML
   	  ENDIF
	ENDIF

 31	rati=ratf
	IF(IS2 .EQ. 0) THEN
	sig=spnsk(srt,0,0,0,ich)*PPN/PMM
	rat1=rati+sig/sref
	rat2=rati+3.*sig/sref

	  IF(xx.GE.rati .and. xx.LT.rat1) THEN
	IF(srt.LE.SLK0) RETURN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	  ELSEIF(xx.GE.rat1 .and. xx.LT.rat2) THEN
	IF(srt.LE.SLK) RETURN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= -1
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	  ELSE
	  RETURN
	  ENDIF
	ELSE
	sig=spnsk(srt,0,1,1,ich)*PPN/PMM
	rat1=rati+sig/sref
	rat2=rati+3.*sig/sref

	  IF(xx.GE.rati .and. xx.LT.rat1) THEN
	IF(srt.LE.SLK) RETURN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 5
	XMK=DMK
	XMH=DMS
	  ELSEIF(xx.GE.rat1 .and. xx.LT.rat2) THEN
	IF(srt.LE.SLK0) RETURN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 1
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	  ELSE
	  RETURN
	  ENDIF
	ENDIF
	GOTO 25

40	CONTINUE
  	IF(ID2 .GE. 1) RETURN
	XMP=0.138
	XMM=0.775

	S=srt**2
      P2PN = (S**2+DMN**4+XMP**4-2.*S*DMN**2-2*S*XMP**2-
     &        2.*DMN**2*XMP**2)/(4.*S)
      P2MM = (S**2+DMN**4+XMM**4-2.*S*DMN**2-2*S*XMM**2-
     &        2.*DMN**2*XMM**2)/(4.*S)

      IF(P2PN .LE. 0) RETURN
      IF(P2MM .LE. 0) RETURN
      PPN = SQRT(P2PN)
      PMM = SQRT(P2MM)

c****** rhoN --> LK0
	IF(srt.LE.SLK0) GOTO 41
	sig=spnlk(srt,IS1,IS2,0)*PPN/PMM
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK0
	XMH=DML
   	ENDIF
c****** rhoN --> LK+
 41  	IF(srt.LE.SLK) GOTO 42
	sig=spnlk(srt,IS1,IS2,1)*PPN/PMM
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 4
	XMK=DMK
	XMH=DML
	ENDIF
c****** rhoN --> SK0
 42	rati=ratf
	IF(srt.LE.SSK0) GOTO 43
	sig=spnsk(srt,IS1,IS2,0,ich)*PPN/PMM
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 0
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK0
	XMH=DMS
	ENDIF
c****** rhoN --> SK+
 43	rati=ratf
	IF(srt.LE.SSK) GOTO 25
	sig=spnsk(srt,IS1,IS2,1,ich)*PPN/PMM
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=3
	icp(i1)= 1
	iso(i2)= ich
	idr(i2)= 5
	XMK=DMK
	XMH=DMS
	ENDIF


 25	CONTINUE
	IF(idms(i1).EQ.1) RETURN

 26   XX       = 1. - 2.*RANFL(mseed)
      YY       = 1. - 2.*RANFL(mseed)
      ZZ       = 1. - 2.*RANFL(mseed)
      RR       = SQRT( XX**2 + YY**2 + ZZ**2 )
      IF((RR .LT. 0.001) .OR. (RR .GT. 1.) )             GOTO 26
	S=srt**2
      PA2 = (S**2+XMH**4+XMK**4-2.*S*XMH**2-2*S*XMK**2-2.*XMH**2*XMK**2)
     &      /(4.*S)
      IF(PA2.LE.0) RETURN
      PA = SQRT(PA2)
      PXX= PA * XX / RR
      PYY= PA * YY / RR
      PZZ= PA * ZZ / RR
      EKAN= SQRT(XMK**2+PA2)
      EHYP= SQRT(XMH**2+PA2)
	bx= -velcm(1)
	by= -velcm(2)
	bz= -velcm(3)
c***** kaon
      CALL LORENTZ(bx,by,bz,EKAN,PXX,PYY,PZZ,
     &     Ekaon,PKX,PKY,PKZ)
c***** lambda or sigma
      CALL LORENTZ(bx,by,bz,EHYP,-PXX,-PYY,-PZZ,
     &     EEHYP,PHYPX,PHYPY,PHYPZ)

	 pxms(i1)=PKX
	 pyms(i1)=PKY
	 pzms(i1)=PKZ
	 emson(i1)=XMK

	 px(i2)= PHYPX
	 py(i2)= PHYPY
	 pz(i2)= PHYPZ
	 em(i2)= XMH

	  kssb=1
        return 
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine prkknhy(srt,i1,i2,sref)
c****** kaon-hyperon collisions
c****** KL --> PN,PD; KS --> PN,PD
c****** Y=lambda, sigma+, 0, -, K=0, +.
c****** i1=pion; i2=N,D
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/ee/mseed,nseed
        common/densi1/denstl(500),dens1(500),velcm(3)
        
        emch=.13957 
	  emne=.13496
	  emav=(2.*emch+emne)/3.
	  enuc=.938

	  if(idms(i1) .NE. 3) return
	  if(idr(i2) .LT. 4) return
	  if(idr(i2) .GT. 5) return
	  ichkn=icp(i1)
	  icht= icp(i1)+iso(i2)

   	  pin2=((srt**2+emson(i1)**2-em(i2)**2)/(2.*srt))**2-emson(i1)**2
   	  pfl2=((srt**2+emav**2-enuc**2)/(2.*srt))**2-emav**2
	  if(pin2 .LE. 0.) return
	  rat=0.
 	  xx=RANFL(nseed)

	if(idr(i2) .GT. 4) goto 11
c****** LK --> PN 
	  csf= 3.*pfl2/pin2
	if(icht.eq.0) then
	  rati=csf*spnlk(srt,0,0,ichkn)/sref
	  ratf=csf*(spnlk(srt,0,0,ichkn)+spnlk(srt,-1,1,ichkn))/sref
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 0
	XMP=emne
	ELSEIF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=1
	icp(i1)= -1
	iso(i2)= 1
	idr(i2)= 0
	XMP=emch
   	ENDIF
	endif

	if(icht.eq.1) then
	  rati=csf*spnlk(srt,1,0,ichkn)/sref
	  ratf=csf*(spnlk(srt,1,0,ichkn)+spnlk(srt,0,1,ichkn))/sref ! ¸üÕý
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 0
	XMP=emch
	ELSEIF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=1
	icp(i1)= 0
	iso(i2)= 1
	idr(i2)= 0
	XMP=emne
   	ENDIF
	endif
	goto 12

c****** SK --> PN
	  csf= pfl2/pin2
 11	if(icht.eq.-1)  then
	  rati=csf*spnsk(srt,-1,0,ichkn,inc)/sref
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= -1
	iso(i2)= 0
	idr(i2)= 0
	XMP=emch
	ENDIF
	endif

	if(icht.eq. 0) then
	rati=csf*spnsk(srt,-1,1,ichkn,inc)/sref
	ratf=csf*(spnsk(srt,-1,1,ichkn,inc)+spnsk(srt,0,0,ichkn,inc))/sref
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= -1
	iso(i2)= 1
	idr(i2)= 0
	XMP=emch
	ELSEIF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=1
	icp(i1)= 0
	iso(i2)= 0
	idr(i2)= 0
	XMP=emne
	ENDIF
	endif

	if(icht.eq. 1) then
	rati=csf*spnsk(srt,0,1,ichkn,inc)/sref
	ratf=csf*(spnsk(srt,0,1,ichkn,inc)+spnsk(srt,1,0,ichkn,inc))/sref
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= 0
	iso(i2)= 1
	idr(i2)= 0
	XMP=emne
	ELSEIF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(i1)=1
	icp(i1)= 1
	iso(i2)= 0
	idr(i2)= 0
	XMP=emch
	ENDIF
	endif

	if(icht.eq. 2) then
	rati=csf*spnsk(srt,1,1,ichkn,inc)/sref
	IF(xx.LT.rati) THEN
	idms(i1)=1
	icp(i1)= 1
	iso(i2)= 1
	idr(i2)= 0
	XMP=emch
	ENDIF
	endif
	IF(idms(i1).NE.1)  RETURN

 12   XX       = 1. - 2.*RANFL(mseed)
      YY       = 1. - 2.*RANFL(mseed)
      ZZ       = 1. - 2.*RANFL(mseed)
      RR       = SQRT( XX**2 + YY**2 + ZZ**2 )
      IF((RR .LT. 0.001) .OR. (RR .GT. 1.) )             GOTO 12
	S=srt**2
      PA2 = (S**2+enuc**4+XMP**4-2.*S*enuc**2-2*S*XMP**2-
     &       2.*enuc**2*XMP**2)/(4.*S)
      IF(PA2.LE.0) RETURN
      PA = SQRT(PA2)
      PXX= PA * XX / RR
      PYY= PA * YY / RR
      PZZ= PA * ZZ / RR
      EPI= SQRT(XMP**2+PA2)
      ENU= SQRT(enuc**2+PA2)
	bx= velcm(1)
	by= velcm(2)
	bz= velcm(3)
c***** kaon
      CALL LORENTZ(bx,by,bz,EPI,PXX,PYY,PZZ,
     &     EPN,PPX,PPY,PPZ)
c***** lambda or sigma
      CALL LORENTZ(bx,by,bz,ENU,-PXX,-PYY,-PZZ,
     &     ENL,PNX,PNY,PNZ)

	 pxms(i1)=PPX
	 pyms(i1)=PPY
	 pzms(i1)=PPZ
	 emson(i1)=XMP

	 px(i2)= PNX
	 py(i2)= PNY
	 pz(i2)= PNZ
	 em(i2)= enuc
	 if(iso(i2).eq.0) isy(i2)= 1
	 if(iso(i2).eq.1) isy(i2)= -1

        return 
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine prbbkkb(srt,i1,i2,sref)
c****** baryon-baryon collisions for the production of kaon-kaonbar
c****** NN --> NNKK^, Ndelta --> NNKK^, delta delta --> NNKK^, Ndelta --> NdeltaKK^, delta delta --> delta deltaKK^
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00  
        common/densi1/denstl(500),dens1(500),velcm(3)
        common/ee/mseed,nseed
        common/medium/mdel,mdin,nsym,npop,nkaon

	if(sref.LE.0.) return
	xx=(x(i1)+x(i2))/2.
	yy=(y(i1)+y(i2))/2.
	zz=(z(i1)+z(i2))/2.
	rho=densit(xx,yy,zz)

	PI=3.141592654
      DMN= 0.938
	DMD= 1.232
      DMK0=0.498               !GeV
      DMK= 0.494

!	in-medium correction on the threshold energies
	  if(nkaon.EQ.0) then
	  ak= 0.
	  atk=0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3   Vk= 27.8 MeV (rho0)
	  atk=0.31         !Vanti-K= -100.3 MeV
	  bk= 0.333        !GeVfm^3
	  endif
	  xkmed = sqrt(DMK**2-ak*denscal(rho)+(bk*rho)**2) + bk*rho
	  xk0med= sqrt(DMK0**2-ak*denscal(rho)+(bk*rho)**2)+ bk*rho
	  xakmd = sqrt(DMK**2-atk*denscal(rho)+(bk*rho)**2) - bk*rho
	  xak0md= sqrt(DMK0**2-atk*denscal(rho)+(bk*rho)**2)- bk*rho

	DMKB0=xak0md
	DMKB= xakmd
	DMK0= xk0med
	DMK=  xkmed

	SNNKK1=2.*DMN+DMK0+DMKB0       !NNK^0KB^0
	SNNKK2=2.*DMN+DMK0+DMKB        !NNK^0KB
	SNNKK3=2.*DMN+DMK+DMKB0        !NNKKB^0
	SNNKK4=2.*DMN+DMK+DMKB         !NNKKB
	SNDKK1=DMN+DMD+DMK0+DMKB0      !NDK^0KB^0
	SNDKK2=DMN+DMD+DMK0+DMKB       !NDK^0KB
	SNDKK3=DMN+DMD+DMK+DMKB0       !NDKKB^0
	SNDKK4=DMN+DMD+DMK+DMKB        !NDKKB
	SDDKK1=DMD+DMD+DMK0+DMKB0      !DDK^0KB^0
	SDDKK2=DMD+DMD+DMK0+DMKB       !DDK^0KB
	SDDKK3=DMD+DMD+DMK+DMKB0       !DDKKB^0
	SDDKK4=DMD+DMD+DMK+DMKB        !DDKKB

	XMK1=0.
	XMK2=0.
	XMB1=0.
	XMB2=0.

	ID1=idr(i1)
	ID2=idr(i2)
	IS1=iso(i1)
	IS2=iso(i2)

	icm1= 0
 98	icm1=icm1+1
	IF(idms(icm1).NE.0) goto 98
	icm2= 0
 99	icm2=icm2+1
	IF(idms(icm2).NE.0) goto 99

      IF(icm1.GT.np .OR. icm2.GT.np) STOP 'TOO MANY MESONS IN prkbyby!'
	xx=RANFL(nseed)
	rati=0.

	IF(ID1.GE.2 .OR. ID2.GE.2) RETURN
	IF(ID1+ID2 .GE. 1) GOTO 15

c****** NN --> NNK0KB0
	IF(srt.LE.SNNKK1) GOTO 11
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,1,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
   	ENDIF
c****** NN --> NNK0KB
 11  	IF(srt.LE.SNNKK2) GOTO 12
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,2,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** NN --> NNKKB0
 12	rati=ratf
	IF(srt.LE.SNNKK3) GOTO 13
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,3,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** NN --> NNKKB
 13	rati=ratf
	IF(srt.LE.SNNKK4) GOTO 505
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,4,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
	GOTO 505

 15	CONTINUE
	IF(ID1+ID2 .GE. 2) GOTO 28

	kn=IS1
	kd=IS2
	IF(ID1.EQ.1) THEN
	kn=IS2
	kd=IS1
	ENDIF

c****** ND --> NNK0KB0
	IF(srt.LE.SNNKK1) GOTO 21
	sig=cnnkkb(srt,kn,kd,DMKB0,DMKB,1,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
   	ENDIF
c****** ND --> NNK0KB
 21  	IF(srt.LE.SNNKK2) GOTO 22
	sig=cnnkkb(srt,kn,kd,DMKB0,DMKB,2,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** ND --> NNKKB0
 22	rati=ratf
	IF(srt.LE.SNNKK3) GOTO 23
	sig=cnnkkb(srt,kn,kd,DMKB0,DMKB,3,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** ND --> NNKKB
 23	rati=ratf
	IF(srt.LE.SNNKK4) GOTO 24
	sig=cnnkkb(srt,kn,kd,DMKB0,DMKB,4,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** ND --> NDK0KB0
 24	rati=ratf
	IF(srt.LE.SNDKK1) GOTO 25
	sig=cndtkkb(srt,kn,kd,DMKB0,DMKB,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 1
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	isy(i2)=0
	XMK1=DMK0
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMD
	ENDIF
c****** ND --> NDK0KB
 25	rati=ratf
	IF(srt.LE.SNDKK2) GOTO 26
	sig=cndtkkb(srt,kn,kd,DMKB0,DMKB,2,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 1
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	isy(i2)=0
	XMK1=DMK0
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMD
	ENDIF
c****** ND --> NDKKB0
 26	rati=ratf
	IF(srt.LE.SNDKK3) GOTO 27
	sig=cndtkkb(srt,kn,kd,DMKB0,DMKB,3,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 1
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	isy(i2)=0
	XMK1=DMK
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMD
	ENDIF
c****** ND --> NDKKB
 27	rati=ratf
	IF(srt.LE.SNDKK4) GOTO 505
	sig=cndtkkb(srt,kn,kd,DMKB0,DMKB,4,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 1
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	isy(i2)=0
	XMK1=DMK
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMD
	ENDIF
	GOTO 505

 28	CONTINUE
	IF(ID1+ID2 .GE. 3) RETURN

c****** DD --> NNK0KB0
	IF(srt.LE.SNNKK1) GOTO 31
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,1,ich1,ich2)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
   	ENDIF
c****** DD --> NNK0KB
 31  	IF(srt.LE.SNNKK2) GOTO 32
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,2,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** DD --> NNKKB0
 32	rati=ratf
	IF(srt.LE.SNNKK3) GOTO 33
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,3,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB0
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** DD --> NNKKB
 33	rati=ratf
	IF(srt.LE.SNNKK4) GOTO 34
	sig=cnnkkb(srt,IS1,IS2,DMKB0,DMKB,4,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 0
	idr(i2)= 0
	if(ich1.eq.1) isy(i1)=-1
	if(ich1.eq.0) isy(i1)=1
	if(ich2.eq.1) isy(i2)=-1
	if(ich2.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB
	XMB1=DMN
	XMB2=DMN
	ENDIF
c****** DD --> DDK0KB0
 34	rati=ratf
	IF(srt.LE.SDDKK1) GOTO 35
	sig=cddnnkkb(srt,IS1,IS2,DMKB0,DMKB,1,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
	XMK1=DMK0
	XMK2=DMKB0
	XMB1=DMD
	XMB2=DMD
	ENDIF
c****** DD --> DDK0KB
 35	rati=ratf
	IF(srt.LE.SDDKK2) GOTO 36
	sig=cddnnkkb(srt,IS1,IS2,DMKB0,DMKB,2,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
	XMK1=DMK0
	XMK2=DMKB
	XMB1=DMD
	XMB2=DMD
	ENDIF
c****** DD --> DDKKB0
 36	rati=ratf
	IF(srt.LE.SDDKK3) GOTO 37
	sig=cddnnkkb(srt,IS1,IS2,DMKB0,DMKB,3,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
	XMK1=DMK
	XMK2=DMKB0
	XMB1=DMD
	XMB2=DMD
	ENDIF
c****** DD --> DDKKB
 37	rati=ratf
	IF(srt.LE.SDDKK2) GOTO 505
	sig=cddnnkkb(srt,IS1,IS2,DMKB0,DMKB,4,ich1,ich2)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i1)=ich1
	iso(i2)=ich2
	idr(i1)= 1
	idr(i2)= 1
	isy(i1)= 0
	isy(i2)= 0
	XMK1=DMK
	XMK2=DMKB
	XMB1=DMD
	XMB2=DMD
	ENDIF

 505	CONTINUE
	IF(idms(icm1).EQ.0 .OR. idms(icm2).EQ.0) RETURN

	XM1=XMK1+XMB1
	XM2=XMK2+XMB2
	SS=srt*srt
	PP2=(SS-(XM1+XM2)**2)*(SS-(XM1-XM2)**2)/(4.*SS)
	IF(PP2.LE.0.) RETURN
	SSKB1=XM1**2+PP2
	SSKB2=XM2**2+PP2
C*****generate B1K1 momentum in NN cms
	  CALL  TWOPARTDECAY(PERR,SSKB1,XMK1,XMB1,PPX,PPY,PPZ,P01,P02)
        IF(PERR.EQ.1)THEN
        PRINT*, 'SSBARPRODUCION:PERR1'
        RETURN
        ENDIF
C*****generate B2K2 momentum in NN cms
        CALL  TWOPARTDECAY(PERR,SSKB2,XMK2,XMB2,PXX,PYY,PZZ,P001,P002)
        IF(PERR.EQ.1) THEN
        PRINT*, 'SSBARPRODUCTION:PERR2'
        RETURN
        ENDIF
C.........Lorentz-transformation to calculational frame
	bx= -velcm(1)
	by= -velcm(2)
	bz= -velcm(3)

      CALL LORENTZ(bx,by,bz,P01,PPX,PPY,PPZ,PK0,PKX,PKY,PKZ)    ! kaon1
       xxx=(x(i1) + x(i2))/2.+(RANFL(nseed)-0.5)*0.4
       yyy=(y(i1) + y(i2))/2.+(RANFL(nseed)-0.5)*0.4
       zzz=(z(i1) + z(i2))/2.+(RANFL(nseed)-0.5)*0.4
	 xms(icm1)=  xxx
	 yms(icm1)=  yyy
	 zms(icm1)=  zzz
	 pxms(icm1)= PKX
	 pyms(icm1)= PKY
	 pzms(icm1)= PKZ
	 emson(icm1)=XMK1
      CALL LORENTZ(bx,by,bz,P02,-PPX,-PPY,-PPZ,W0,WX,WY,WZ)     ! baryon1
	 px(i1)= WX
	 py(i1)= WY
	 pz(i1)= WZ
	 em(i1)= XMB1
      CALL LORENTZ(bx,by,bz,P001,PXX,PYY,PZZ,PK0,PKX,PKY,PKZ)    ! kaon2
       xxx=(x(i1) + x(i2))/2.+(RANFL(nseed)-0.5)*0.4
       yyy=(y(i1) + y(i2))/2.+(RANFL(nseed)-0.5)*0.4
       zzz=(z(i1) + z(i2))/2.+(RANFL(nseed)-0.5)*0.4
	 xms(icm2)=  xxx
	 yms(icm2)=  yyy
	 zms(icm2)=  zzz
	 pxms(icm2)= PKX
	 pyms(icm2)= PKY
	 pzms(icm2)= PKZ
	 emson(icm2)=XMK2
      CALL LORENTZ(bx,by,bz,P002,-PXX,-PYY,-PZZ,W0,WX,WY,WZ)     ! baryon2
	 px(i2)= WX
	 py(i2)= WY
	 pz(i2)= WZ
	 em(i2)= XMB2

        return 
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine prpnkkb(srt,i1,i2,sref,kkb)
c****** pion-baryon collisions for the production of kaon-kaonbar
c****** PN --> NKK^B, PD --> NKK^B
c****** i1=pion; i2=N,D
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00  
        common/densi1/denstl(500),dens1(500),velcm(3)
        common/ee/mseed,nseed
        common/medium/mdel,mdin,nsym,npop,nkaon

	 kkb=0
	 if(sref.LE.0.) return

	xx=(xms(i1)+x(i2))/2.
	yy=(yms(i1)+y(i2))/2.
	zz=(zms(i1)+z(i2))/2.
	rho=densit(xx,yy,zz)

	PI=3.141592654
      DMN= 0.938
	DMD= 1.232
      DMK0=0.498               !GeV
      DMK= 0.494

!	in-medium correction on the threshold energies
	  if(nkaon.EQ.0) then
	  ak= 0.
	  atk= 0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3   Vk= 27.8 MeV (rho0)
	  atk=0.31         !Vanti-K= -100.3 MeV
	  bk= 0.333        !GeVfm^3
	  endif
	  xkmed = sqrt(DMK**2-ak*denscal(rho)+(bk*rho)**2) + bk*rho
	  xk0med= sqrt(DMK0**2-ak*denscal(rho)+(bk*rho)**2)+ bk*rho
	  xakmd = sqrt(DMK**2-atk*denscal(rho)+(bk*rho)**2) - bk*rho
	  xak0md= sqrt(DMK0**2-atk*denscal(rho)+(bk*rho)**2)- bk*rho

	DMKB0=xak0md
	DMKB= xakmd
	DMK0= xk0med
	DMK=  xkmed

	SNKK1=DMN+DMK0+DMKB0       !NK^0KB^0
	SNKK2=DMN+DMK0+DMKB        !NK^0KB
	SNKK3=DMN+DMK+DMKB0        !NKKB^0
	SNKK4=DMN+DMK+DMKB         !NKKB

	IS1=icp(i1)
	IS2=iso(i2)
	XMK1=0.
	XMK2=0.

	icm1= 0
 100	icm1=icm1+1
	IF(idms(icm1).NE.0) goto 100
	icm2= 0
 101	icm2=icm2+1
	IF(idms(icm2).NE.0) goto 101
	xx=RANFL(nseed)

	IF(idr(i2).GE.1) goto 15
c****** PN --> NK0KB0
	IF(srt.LE.SNKK1) GOTO 11
	sig=cpbkkb(srt,IS1,IS2,DMKB0,DMKB,1,ich)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB0
   	ENDIF
c****** PN --> NK0KB
 11  	IF(srt.LE.SNKK2) GOTO 12
	sig=cpbkkb(srt,IS1,IS2,DMKB0,DMKB,2,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB
	ENDIF
c****** PN --> NKKB0
 12  	IF(srt.LE.SNKK3) GOTO 13
	sig=cpbkkb(srt,IS1,IS2,DMKB0,DMKB,3,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB0
	ENDIF
c****** PN --> NKKB
 13  	IF(srt.LE.SNKK4) GOTO 506
	sig=cpbkkb(srt,IS1,IS2,DMKB0,DMKB,4,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB
	ENDIF
	GOTO 506

 15	continue
	IF(idr(i2).GE.2) RETURN

c****** PD --> NK0KB0
	IF(srt.LE.SNKK1) GOTO 21
	sig=cspdkkb(srt,IS1,IS2,DMKB0,DMKB,1,ich)
	rati=sig/sref
	IF(xx.LT.rati) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= 0
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB0
   	ENDIF
c****** PD --> NK0KB
 21  	IF(srt.LE.SNKK2) GOTO 22
	sig=cspdkkb(srt,IS1,IS2,DMKB0,DMKB,2,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 0
	icp(icm2)= -1
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK0
	XMK2=DMKB
	ENDIF
c****** PD --> NKKB0
 22  	IF(srt.LE.SNKK3) GOTO 23
	sig=cspdkkb(srt,IS1,IS2,DMKB0,DMKB,3,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= 0
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB0
	ENDIF
c****** PD --> NKKB
 23  	IF(srt.LE.SNKK4) GOTO 506
	sig=cspdkkb(srt,IS1,IS2,DMKB0,DMKB,4,ich)
	ratf=rati+sig/sref
	IF(xx.GE.rati .and. xx.LT.ratf) THEN
	idms(icm1)=3
	idms(icm2)=4
	icp(icm1)= 1
	icp(icm2)= -1
	iso(i2)= ich
	idr(i2)= 0
	if(ich.eq.1) isy(i2)=-1
	if(ich.eq.0) isy(i2)=1
	XMK1=DMK
	XMK2=DMKB
	ENDIF

 506	CONTINUE
	IF(idms(icm1).EQ.0 .OR. idms(icm2).EQ.0) RETURN

	PPN=(srt*srt-(XMK1+XMK2+DMN)**2)*(srt*srt-(XMK1+XMK2-DMN)**2)/
     &    (4.*srt*srt)
	IF(PPN.LE.0.) RETURN
C..generate N momentum in PB cms system
      CALL EVENT(PPN,DMN,PNX,PNY,PNZ,PN0)
c..beta-factor for lorentz-boost to nucleon-system:
      BETX=PNX/(srt-PN0)
      BETY=PNY/(srt-PN0)
      BETZ=PNZ/(srt-PN0)
 31    XX    = 1. - 2.*RANFL(mseed)
       YY    = 1. - 2.*RANFL(mseed)
       ZZ    = 1. - 2.*RANFL(mseed)
       RR    = SQRT( XX**2 + YY**2 + ZZ**2 )
      IF((RR .LT. 0.001) .OR. (RR .GT. 1.))             GOTO 31

      S=(srt-PN0)**2-PNX**2-PNY**2-PNZ**2
      PA2 = (S**2+XMK1**4+XMK2**4-2.*S*XMK1**2-2*S*XMK2**2-
     &       2.*XMK1**2*XMK2**2)/(4.*S)
      IF(PA2.LE.0) RETURN
      PA = SQRT(PA2)
      PXX= PA * XX / RR
      PYY= PA * YY / RR
      PZZ= PA * ZZ / RR
      EKN= SQRT(XMK1**2+PA2)
      EKB= SQRT(XMK2**2+PA2)
* --------- TRANSFORMATION INTO nucleon-system:
c***** kaon
      CALL LORENTZ(BETX,BETY,BETZ,EKN,PXX,PYY,PZZ,
     &     EKNN,PPX,PPY,PPZ)
c***** antikaon
      CALL LORENTZ(BETX,BETY,BETZ,EKB,-PXX,-PYY,-PZZ,
     &     EKBN,PKBX,PKBY,PKBZ)
C.........Lorentz-transformation to calculational frame
	bx= -velcm(1)
	by= -velcm(2)
	bz= -velcm(3)
      CALL LORENTZ(bx,by,bz,PN0,PNX,PNY,PNZ,P0N,PXN,PYN,PZN)          ! nucleon
      CALL LORENTZ(bx,by,bz,EKNN,PPX,PPY,PPZ,U0,UX,UY,UZ)             ! kaon
      CALL LORENTZ(bx,by,bz,EKBN,PKBX,PKBY,PKBZ,W0,WX,WY,WZ)          ! antikaon

	 px(i2)= PXN
	 py(i2)= PYN
	 pz(i2)= PZN
	 em(i2)= DMN

	 xms(icm1)=  xms(i1)
	 yms(icm1)=  yms(i1)
	 zms(icm1)=  zms(i1)
	 pxms(icm1)= UX
	 pyms(icm1)= UY
	 pzms(icm1)= UZ
	 emson(icm1)=XMK1

       xxx=(xms(i1) + x(i2))/2.+(RANFL(nseed)-0.5)*0.4
       yyy=(yms(i1) + y(i2))/2.+(RANFL(nseed)-0.5)*0.4
       zzz=(zms(i1) + z(i2))/2.+(RANFL(nseed)-0.5)*0.4
	 xms(icm2)=  xxx
	 yms(icm2)=  yyy
	 zms(icm2)=  zzz
	 pxms(icm2)= WX
	 pyms(icm2)= WY
	 pzms(icm2)= WZ
	 emson(icm2)=XMK2

	  kkb=1
        return
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine pnpnkkb(srt,i1,i2,sref)
c****** pion-pion collisions for the production of kaon-kaonbar
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00  
        common/densi1/denstl(500),dens1(500),velcm(3)
        common/ee/mseed,nseed
        common/medium/mdel,mdin,nsym,npop,nkaon

      if(sref.LE.0.) return

	xx=(xms(i1)+xms(i2))/2.
	yy=(yms(i1)+yms(i2))/2.
	zz=(zms(i1)+zms(i2))/2.
	rho=densit(xx,yy,zz)

      DMK0=0.498               !GeV
      DMK= 0.494
!	in-medium correction on the threshold energies
	  if(nkaon.EQ.0) then
	  ak= 0.
	  atk=0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3   Vk= 27.8 MeV (rho0)
	  atk=0.31         !Vanti-K= -100.3 MeV
	  bk= 0.333        !GeVfm^3
	  endif
	  xkmed = sqrt(DMK**2-ak*denscal(rho)+(bk*rho)**2) + bk*rho
	  xk0med= sqrt(DMK0**2-ak*denscal(rho)+(bk*rho)**2)+ bk*rho
	  xakmd = sqrt(DMK**2-atk*denscal(rho)+(bk*rho)**2) - bk*rho
	  xak0md= sqrt(DMK0**2-atk*denscal(rho)+(bk*rho)**2)- bk*rho

	DMKB0=xak0md
	DMKB= xakmd
	DMK0= xk0med
	DMK=  xkmed

	SKK1=DMK0+DMKB0       !K^0KB^0
	SKK2=DMK0+DMKB        !K^0KB
	SKK3=DMK+DMKB0        !KKB^0
	SKK4=DMK+DMKB         !KKB
	xx=RANFL(nseed)
	npch=icp(i1)+icp(i2)
	if(npch.lt.-1 .OR. npch.gt.1) return

	IF(npch.eq.-1) THEN
	sig=cpnpnkkb(srt,DMKB0,DMKB,2)
	rati=sig/sref
	IF(xx.GT.rati) RETURN
	idms(i1)=3
	idms(i2)=4
	icp(i1) =0
	icp(i2) =-1
	XMK1=DMK0
	XMK2=DMKB
	ELSEIF(npch.eq.0) THEN
	sig=cpnpnkkb(srt,DMKB0,DMKB,1)+cpnpnkkb(srt,DMKB0,DMKB,4)
	rati=sig/sref
	IF(xx.GT.rati) RETURN
	rkz=cpnpnkkb(srt,DMKB0,DMKB,1)/sig
	yy=RANFL(mseed)
	if(yy.lt.rkz) then
	idms(i1)=3
	idms(i2)=4
	icp(i1) =0
	icp(i2) =0
	XMK1=DMK0
	XMK2=DMKB0
	else
	idms(i1)=3
	idms(i2)=4
	icp(i1) =1
	icp(i2) =-1
	XMK1=DMK
	XMK2=DMKB
	endif
	ELSE
	sig=cpnpnkkb(srt,DMKB0,DMKB,3)
	rati=sig/sref
	IF(xx.GT.rati) RETURN
	idms(i1)=3
	idms(i2)=4
	icp(i1) =1
	icp(i2) =0
	XMK1=DMK
	XMK2=DMKB0
	ENDIF

	PP2=(srt*srt-(XMK1+XMK2)**2)*(srt*srt-(XMK1-XMK2)**2)/
     &    (4.*srt*srt)
	IF(PP2.LE.0.) RETURN
	PR=SQRT(PP2)
 111   XX    = 1. - 2.*RANFL(mseed)
       YY    = 1. - 2.*RANFL(mseed)
       ZZ    = 1. - 2.*RANFL(mseed)
       RR    = SQRT( XX**2 + YY**2 + ZZ**2 )
      IF((RR .LT. 0.001) .OR. (RR .GT. 1.))             GOTO 111
      PXX= PR * XX / RR
      PYY= PR * YY / RR
      PZZ= PR * ZZ / RR
      EKN= SQRT(XMK1**2+PP2)
      EKB= SQRT(XMK2**2+PP2)
	bx= -velcm(1)
	by= -velcm(2)
	bz= -velcm(3)
      CALL LORENTZ(bx,by,bz,EKN,PXX,PYY,PZZ,U0,UX,UY,UZ)             ! kaon
      CALL LORENTZ(bx,by,bz,EKB,-PXX,-PYY,-PZZ,W0,WX,WY,WZ)          ! antikaon

	 pxms(i1)= UX
	 pyms(i1)= UY
	 pzms(i1)= UZ
	 emson(i1)=XMK1
	 pxms(i2)= WX
	 pyms(i2)= WY
	 pzms(i2)= WZ
	 emson(i2)=XMK2

        return
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine stranex(prx,pry,prz,srt,i1,i2,sref)
c****** strange exchange reactions
c****** PY --> KbN
c****** i1=pion; i2=lambda, sigma
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/ee/mseed,nseed

	IF(sref.LE.0.) RETURN
      IF(srt.LE.1.436) RETURN

	 PI=3.141592654
	 DMN= 0.938
       DMK0=0.498               !GeV
       DMK= 0.494
	 DMKA=(DMK0+DMK)/2.
       DML= 1.115
       DMS= 1.189
	 DMPC=0.13957
	 DMP0=0.13496
	 DMPA=(2.*DMPC+DMP0)/3.

        pr=sqrt(prx*prx+pry*pry+prz*prz)
        c2=prz/pr

        TTX=pxms(i1)-px(i2)
        TTY=pyms(i1)-py(i2)
        TTZ=pzms(i1)-pz(i2)
        PKB=SQRT(TTX**2+TTY**2+TTZ**2)
	  xx=RANFL(nseed)
	  yy=RANFL(mseed)
        SS=srt*srt
	  prilm2=(SS-(DMPA+DML)**2)*(SS-(DMPA-DML)**2)/(4.*SS)
	  prisg2=(SS-(DMPA+DMS)**2)*(SS-(DMPA-DMS)**2)/(4.*SS)
	  prf2=  (SS-(DMKA+DMN)**2)*(SS-(DMKA-DMN)**2)/(4.*SS)
	  IF(prilm2.LE.0.) RETURN
	  IF(prisg2.LE.0.) RETURN
	  IF(prf2.LE.0.)   RETURN
	  factlam=prf2/prilm2
	  factsig=prf2/prisg2
        RPL=2.*factlam*CROKB(srt,PKB,4)/sref
        RPS=2.*factsig*CROKB(srt,PKB,5)/sref+RPL
        IF(xx.LE.RPL) THEN
	  IF(idr(i2).NE.4) RETURN
	if(yy.LT.0.5) then
	ickb= 0
	XMKB= DMK0
	else
	ickb= -1
	XMKB= DMK
	endif
	icn= icp(i1)-ickb
	if(icn.lt.0 .OR. icn.gt.1) return
	  ELSEIF(xx.GT.RPL .AND. xx.LE.RPS) THEN
	  IF(idr(i2).NE.5) RETURN
	if(yy.LT.0.5) then
	ickb= 0
	XMKB= DMK0
	else
	ickb= -1
	XMKB= DMK
	endif
	icn= icp(i1)+iso(i2)-ickb
	if(icn.lt.0 .OR. icn.gt.1) return
	  ELSE
	  RETURN
	  ENDIF

	  pr=sqrt(prf2)
	  emson(i1)=XMKB
	  idms(i1)=4
	  icp(i1)=ickb
	  em(i2)=DMN
	  idr(i2)=0.
	  iso(i2)=icn
	  if(iso(i2).eq.1) isy(i2)=-1
	  if(iso(i2).eq.0) isy(i2)=1

        c1=1.-2.*RANFL(mseed)
        t1=2.*PI*RANFL(mseed)
        s2=sqrt(1.-c2*c2)
        if(prx.eq.0.0.and.pry.eq.0.0) goto 55
        t2=atan2(pry,prx)
 55     if(c1.gt.1)  c1=1.
	  if(c1.lt.-1) c1=-1.
	  s1=sqrt(1.-c1*c1)
        ct1=cos(t1)
        st1=sin(t1)
        ct2=cos(t2)
        st2=sin(t2)

        prz=pr*(c1*c2-s1*s2*ct1)
        ss =c2*s1*ct1+s2*c1 
        prx=pr*(ss*ct2-s1*st1*st2)
        pry=pr*(ss*st2+s1*st1*ct2)

        return
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine elaskaon(prx,pry,prz,srt,i1,i2,sref)
c****** elastic scattering and charge-exchange reaction for kaon-baryon collision
c****** KB --> KB, KB --> K'B'
c****** i1=kaon; i2=N,D
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/ee/mseed,nseed
	  IF(sref.LE.0.) RETURN
	  IF(idr(i2).GE.2) RETURN
	  PI=3.141592654

        TTX=pxms(i1)-px(i2)
        TTY=pyms(i1)-py(i2)
        TTZ=pzms(i1)-pz(i2)
        PKB=SQRT(TTX**2+TTY**2+TTZ**2)
	  xx=RANFL(nseed)

        pr=sqrt(prx*prx+pry*pry+prz*prz)
        c2=prz/pr
        REL=CKELS(srt,PKB)/sref
        IF(xx.GT.REL) RETURN

	  ich=icp(i1)+iso(i2)
	  iex=0
	  IF((idr(i2).LE.0) .AND. (ich.EQ.1)) THEN
	  iex=1
	  ikc=1-icp(i1)
	  ELSEIF((idr(i2).EQ.1) .AND. (ich.LE.2 .and. ich.GE.0))THEN
	  iex=1
	  ikc=1-icp(i1)
	  ENDIF
	  IF((iex.GT.0) .AND. (RANFL(mseed).LT.0.5)) THEN   !CHARGE EXCHARGE
	  icp(i1)=ikc
	  iso(i2)=ich-ikc
	  ENDIF

        c1=1.-2.*RANFL(mseed)
        t1=2.*PI*RANFL(mseed)
        s2=sqrt(1.-c2*c2)
        if(prx.eq.0.0.and.pry.eq.0.0) goto 510
        t2=atan2(pry,prx)
 510    if(c1.gt.1)  c1=1.
	  if(c1.lt.-1) c1=-1.
	  s1=sqrt(1.-c1*c1)
        ct1=cos(t1)
        st1=sin(t1)
        ct2=cos(t2)
        st2=sin(t2)

        prz=pr*(c1*c2-s1*s2*ct1)
        ss =c2*s1*ct1+s2*c1 
        prx=pr*(ss*ct2-s1*st1*st2)
        pry=pr*(ss*st2+s1*st1*ct2)

        return 
        end
c///////////////////////////////////////////////////////////////////////////////
     	subroutine absantikn(prx,pry,prz,srt,i1,i2,sref)
c****** antikaon-baryon collisions for antikaon absorption
c****** KbB --> KbB, KbN --> PY
c****** i1=antikaon; i2=N,D
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/ee/mseed,nseed

	IF(sref.LE.0.)   RETURN
      IF(srt.LE.1.254) RETURN
      IF(idr(i2).GE.1) RETURN

	PI=3.141592654
      DML= 1.115
      DMS= 1.189
	DMPC=0.13957
	DMP0=0.13496

        TTX=pxms(i1)-px(i2)
        TTY=pyms(i1)-py(i2)
        TTZ=pzms(i1)-pz(i2)
        PKB=SQRT(TTX**2+TTY**2+TTZ**2)
	  xx=RANFL(nseed)

        pr=sqrt(prx*prx+pry*pry+prz*prz)
        c2=prz/pr
        REL=CROKB(srt,PKB,2)/sref
        RPL=CROKB(srt,PKB,4)/sref+REL
        RPS=3.*CROKB(srt,PKB,5)/sref+REL+RPL
        IF(xx.LE.REL) THEN
	  GOTO 51 
	  ELSEIF(xx.GT.REL .AND. xx.LE.RPL) THEN
	icpn= icp(i1)+iso(i2)
	idr(i2)= 4
	iso(i2)= 0
	idms(i1)=1
	icp(i1)= icpn
	if(icp(i1).eq.1 .OR. icp(i1).eq.-1) XM1=DMPC
	if(icp(i1).eq.0) XM1=DMP0
	XM2=DML
	  ELSEIF(xx.GT.RPL .AND. xx.LE.RPS) THEN
	icnm= icp(i1)+iso(i2)
	yy=RANFL(mseed)
	if(yy.LT.0.33333) then
	icsg= 1
	elseif(yy.GE.0.33333 .and. yy.LT.0.66667) then
	icsg= 0
	else
	icsg= -1
	endif
	icpn= icnm-icsg
	if(abs(icpn).GT.1) return
	if(icpn.eq.1 .OR. icpn.eq.-1) XM1=DMPC
	if(icpn.eq.0) XM1=DMP0
	idr(i2)= 5
	iso(i2)= icsg
	idms(i1)=1
	icp(i1)= icpn
	XM2=DMS
	  ELSE
	  RETURN
	  ENDIF

        SS=srt*srt
 	  PP2=(SS-(XM1+XM2)**2)*(SS-(XM1-XM2)**2)/(4.*SS)
	  IF(PP2.LE.0.) RETURN
	  pr=sqrt(PP2)
	  emson(i1)=XM1
	  em(i2)=XM2
	  isy(i2)=0

 51     c1=1.-2.*RANFL(mseed)
        t1=2.*PI*RANFL(mseed)
        s2=sqrt(1.-c2*c2)
        if(prx.eq.0.0.and.pry.eq.0.0) goto 55
        t2=atan2(pry,prx)
 55     if(c1.gt.1)  c1=1.
	  if(c1.lt.-1) c1=-1.
	  s1=sqrt(1.-c1*c1)
        ct1=cos(t1)
        st1=sin(t1)
        ct2=cos(t2)
        st2=sin(t2)

        prz=pr*(c1*c2-s1*s2*ct1)
        ss =c2*s1*ct1+s2*c1 
        prx=pr*(ss*ct2-s1*st1*st2)
        pry=pr*(ss*st2+s1*st1*ct2)

        return 
        end
c******* antibaryon-baryon collisions
c///////////////////////////////////////////////////////////////////////////////
     	  subroutine proantinn(srt,i1,i2)
c******* elastic, charge exchange, annihilation, meson and hyperon production
c******* include the final channels of ann., NbarN, NbarDelta, DbarN, YYbar, XiXibar
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/ee/mseed,nseed
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/densi1/denstl(500),dens1(500),velcm(3)
	  common/signbd/xsNbDt(0:200), xsNbNl(0:200)

	  pi= 3.1415926
	  xpi=31.415926
        xminsrts=2.014
	  delta_srts=0.01
        DMN = 0.938
        DML = 1.115
        DMS = 1.189
	  DMD = 1.232
	  DMX = 1.315

	  sanni=0.
        santiNuNu=0.
        santiNuDt=0.
	  sNuantiDt=0.
	  slambar=0.
	  ssigbar=0.
	  sxibar=0.

	 is1=iso(i1)
	 is2=iso(i2)
	 id1=idr(i1)
	 id2=idr(i2)
	 ia1=idanti(i1)
	 ia2=idanti(i2)
	 if(idr(i1).GT.idr(i2)) then
	 is1=iso(i2)
	 is2=iso(i1)
	 id1=idr(i2)
	 id2=idr(i1)
	 ia1=idanti(i2)
	 ia2=idanti(i1)
	 endif
	 if(ia1*ia2 .GT. 0) return   !exclude antibaryon-antibaryon collisions
	  icht= is1+is2
	  idmt= id1+id2

  	  call  SIGANTIBB(srt,i1,i2,stotal,selastic,schex,sanni,sigmeson,
     &                  shyperon,slambar,ssigbar,sxibar)

	  if(stotal.LE.0.) return

        x1=x(i1)
        y1=y(i1)
        z1=z(i1)
        x2=x(i2)
        y2=y(i2)
        z2=z(i2)
        r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
        pz1=pz(i1)
        pz2=pz(i2)
        px1=px(i1)
        px2=px(i2)
        py1=py(i1)
        py2=py(i2)
	  em1=em(i1)
        em2=em(i2)
        e1=sqrt(em(i1)*em(i1)+pz1*pz1+px1*px1+py1*py1) 
        e2=sqrt(em(i2)*em(i2)+pz2*pz2+px2*px2+py2*py2)
        bx=velcm(1)
        by=velcm(2)
        bz=velcm(3)
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
   	  pr2=((srt**2+em(i1)**2-em(i2)**2)/(2.*srt))**2-em(i1)**2
        pri= sqrt(max(0.,pr2))

       if(idmt.eq.0) then          ! nucleon+antinucleon

	  santiNuNu= selastic+schex

	  if(srt .LE. xminsrts) then
	  sigma0= 0.
	  else
        srts_Index=NINT((srt-xminsrts)/delta_srts)
	  if(srts_Index .GT. 200)  srts_Index=200
	  sigma0= xsNbDt(srts_Index)
	  endif
        if(icht.eq.0) then
        santiNuDt= 2./3.*sigma0/pri
        else
        santiNuDt= 4./3.*sigma0/pri
        end if
        sNuantiDt= santiNuDt

       else if(idmt.eq.1) then     ! antinucleon+delta or nucleon+antidelta
	  
        if(ia1.NE.0) then
	  santiNuDt=selastic+schex
	  sNuantiDt=0.
	  else
	  santiNuDt=0.
	  sNuantiDt=selastic+schex
	  endif

        if(abs(icht) .LE. 1) then
	   if(srt .LE. xminsrts) then
	   sigma0= 0.
	   else
         srts_Index=NINT((srt-xminsrts)/delta_srts)
	   if(srts_Index .GT. 200)  srts_Index=200
	   sigma0= xsNbNl(srts_Index)
	   endif
          if((ia2.NE.0) .and. ( is2.eq.-2 .OR. is2.eq.1 ) .OR.
     &       (ia2.EQ.0) .and. ( is2.eq.-1 .OR. is2.eq.2 )) then
          santiNuNu= sigma0/pri
          else if(icht.eq.0) then
          santiNuNu= 2./3.*sigma0/pri
          else
          santiNuNu= sigma0/3./pri
          end if
        else
	  santiNuNu=0.
        end if

       end if

       if(idmt .LE. 1) then
       sigmatot= sanni + santiNuNu + santiNuDt + sNuantiDt + slambar +
     &           ssigbar + sxibar
       else
       sigmatot= sanni + selastic
       end if

       If(sigmatot.LT.0.001) return
	 dint=sigmatot/xpi

        b2=.0001+bx*bx+by*by+bz*bz 
        pb=bx*px1+by*py1+bz*pz1 
        gg=(gam-1.)*pb/b2-gam*e1 
        pzr=gg*bz+pz1  
        pxr=gg*bx+px1 
        pyr=gg*by+py1
        rb=bx*(x1-x2)+by*(y1-y2)+bz*(z1-z2)
        gg=(gam-1.)*rb/b2 
        dx=gg*bx+x1-x2  
        dy=gg*by+y1-y2 
        dz=gg*bz+z1-z2
        dr=sqrt(dx*dx+dy*dy+dz*dz)
        pr=sqrt(pxr*pxr+pyr*pyr+pzr*pzr)
        dzz=(pxr*dx+pyr*dy+pzr*dz)/pr                                           
        c00=dr*dr-dzz*dzz                           !!  feng
	  if(c00.le.0.) c00=0.
        bb=sqrt(c00)
        if(bb.gt.sqrt(dint)) return
        e1r=sqrt(em(i1)*em(i1)+pr*pr)                                        
        e2r=sqrt(em(i2)*em(i2)+pr*pr)                                           
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.                                         
        if(dd.lt.abs(dzz))   return
        c2=pzr/pr

!	Monte-Carlo sampling
        cut=sigmatot*RANFL(nseed)

        if(sanni .GE. cut) then                         ! annihilation
	  iso(i1)=0
	  iso(i2)=0
	  isy(i1)=0
	  isy(i2)=0
	  idr(i1)=0
	  idr(i2)=0
	  idanti(i1)=0
	  idanti(i2)=0

	  return
        end if

       if(idmt.GT.1) then       !only elastic scattering
	 goto 515
	 endif


       if((sanni + santiNuNu) .GE. cut ) then        ! antinucleon+nucleon production

       if(idmt.EQ.0) then      ! antinucleon+nucleon collision:

        if(icht .EQ. 0) then
        scut=santiNuNu*RANFL(nseed)
         if(selastic .GE. scut) then
	   goto 515
         else
         iso(i1)= -1 - is1
         end if
        else
	  goto 515
        end if

       else       ! antinucleon+Delta or antiDelta+nucleon collision:

        if(icht .EQ. 0) then
        xx=RANFL(nseed)
          if(xx .LT. 0.5) then
          iso(i1)= -1
          else
          iso(i1)= 0
          end if
        else if(icht .EQ. -1) then
          iso(i1)= -1
        else if(icht .EQ. 1) then
          iso(i1)= 0
        else
	  return
        end if
	 em(i1)=DMN
	 em(i2)=DMN
       end if

	  iso(i2)= icht-iso(i1)
	  if(abs(iso(i1)).EQ.1)  isy(i1)= -1
	  if(abs(iso(i2)).EQ.1)  isy(i2)= -1
	  if(iso(i1).EQ.0)  isy(i1)= 1
	  if(iso(i2).EQ.0)  isy(i2)= 1
	  idr(i1)= 0
	  idr(i2)= 0
	  idanti(i1)= 1
	  idanti(i2)= 0
	  goto 515
 
       else if((sanni + santiNuNu + santiNuDt) .GE. cut) then  ! antinucleon+Delta production
	          
	  if(idmt.EQ.0) then    ! antinucleon+nucleon collision:
        xx=RANFL(nseed)
         if(icht .EQ. 0) then
          if(xx .LT. 0.5) then
          iso(i1)= -1
          else
          iso(i1)= 0
          end if
         else if(icht .EQ. -1) then
          if(xx .LT. 0.25) then
          iso(i1)= -1
          else
          iso(i1)= 0
          end if
         else if(icht .EQ. 1) then
          if(xx .LT. 0.25) then
          iso(i1)= 0
          else
          iso(i1)= -1
          end if
         end if

        else        ! only antinucleon+delta collision:
	  if(ia2.GT.0)  return  !exclude nucleon+antidelta collision

        if(abs(icht) .EQ. 2) then
	  goto 515
        else
        scut=santiNuDt*RANFL(nseed)
         if(selastic .GE. scut) then
	   goto 515
         else
         iso(i1)= -1 - is1
         end if
        end if

        end if

	  iso(i2)= icht-iso(i1)
	  if(abs(iso(i1)).EQ.1)  isy(i1)= -1
	  if(iso(i1).EQ.0)  isy(i1)= 1
	  idr(i1)= 0
	  idr(i2)= 1
	  idanti(i1)= 1
	  idanti(i2)= 0
	  isy(i2)= 0
	  em(i1)=DMN
	  em(i2)=DMD
 	  goto 515

	 else if((sanni + santiNuNu + santiNuDt+ sNuantiDt) .GE. cut) then   ! nucleon+antiDelta production

        If(idmt.EQ.0) then      ! antinucleon+nucleon collision:
        xx=RANFL(nseed)

        if(icht .EQ. 0 ) then
         if(xx .LT. 0.5) then
           iso(i2)= -1
         else
           iso(i2)= 0
         end if
        else if(icht .EQ. -1) then
         if(xx .LT. 0.25) then
           iso(i2)= -1
         else
           iso(i2)= -2
         end if
        else if(icht .EQ. 1) then
         if(xx .LT. 0.25) then
           iso(i2)= 0
         else
           iso(i2)= 1
         end if
        end if
    	  iso(i1)= icht-iso(i2)

        else       ! only antiDelta+nucleon collision:
	  if(ia1.GT.0)  return  !exclude antinucleon+delta collision

        if( abs(icht) .EQ. 2 ) then
 	  goto 515
        else
        scut=sNuantiDt*RANFL(nseed)
         if(selastic .GE. scut) then
 	   goto 515
         else
          iso(i1)= 1 - is1
         end if
        end if
    	  iso(i2)= icht-iso(i1)

        end if

	  if(abs(iso(i1)).EQ.1)  isy(i1)= -1
	  if(iso(i1).EQ.0)  isy(i1)= 1
   	  idr(i1)= 0
	  idr(i2)= 1
	  idanti(i1)= 0
	  idanti(i2)= 1
	  isy(i2)= 0
	  em(i1)=DMN
	  em(i2)=DMD
 	  goto 515

       else if((sanni+santiNuNu+santiNuDt+sNuantiDt+slambar) .GE. cut) 
     &       then      ! Lambda AntiLambda production
       if(icht .NE. 0)  return

   	  idr(i1)= 4
	  idr(i2)= 4
	  em(i1)= DML
	  em(i2)= DML
	  idanti(i1)= 0
	  idanti(i2)= 1
	  iso(i1)= 0
	  iso(i2)= 0
	  isy(i1)= 0
	  isy(i2)= 0
 	  goto 515

       else if((sanni+santiNuNu+santiNuDt+sNuantiDt+slambar+ssigbar) 
     &         .GE. cut) then    ! Lambda SigmaBar / LambdaBar Sigma production

        if(abs(icht) .GT. 1) return
	  isy(i1)= 0
	  isy(i2)= 0

        xx=RANFL(nseed)
        if(xx .LT. 0.5) then    !LambdaBar Sigma
   	  idr(i1)= 4
	  idr(i2)= 5
	  em(i1)= DML
	  em(i2)= DMS
	  idanti(i1)= 1
	  idanti(i2)= 0
	  iso(i1)= 0
	  iso(i2)= icht
        else
   	  idr(i1)= 4
	  idr(i2)= 5
	  em(i1)= DML
	  em(i2)= DMS
	  idanti(i1)= 0
	  idanti(i2)= 1
	  iso(i1)= 0
	  iso(i2)= icht
        endif
 	  goto 515

        else if((sanni+santiNuNu+santiNuDt+sNuantiDt+slambar+ssigbar
     &           + sxibar) .GE. cut) then   ! Xi+XiBar production

	  isy(i1)= 0
	  isy(i2)= 0

       if( icht .EQ. 0 ) then
	  
        xx=RANFL(nseed)
        if(xx .LT. 0.5) then             ! Xi^- + XiBar^+  
	  iso(i1)= 1
        else                             ! Xi^0 + XiBar^0  
	  iso(i1)= 0
        end if
       else if(icht .EQ. -1) then     ! Xi^- + XiBar^0  
        iso(i1)= 0
       else if(icht .EQ. 1) then      ! Xi^0 + XiBar^+
        iso(i1)= 1
       else
        return
       end if

        iso(i2)= icht - iso(i1)
   	  idr(i1)= 6
	  idr(i2)= 6
	  em(i1)= DMX
	  em(i2)= DMX
	  idanti(i1)= 1
	  idanti(i2)= 0


       else
	 return
	 endif

 515   pr2=((srt**2+em(i1)**2-em(i2)**2)/(2.*srt))**2-em(i1)**2
	 if(pr2 .LT. 0.)  pr2=0.
	 pr=sqrt(pr2)

	 c1=1.-2.*RANFL(mseed)
       t1=2.*pi*RANFL(mseed)
       s2=sqrt(1.-c2*c2)
        if(pxr.eq.0. .or. pyr.eq.0.) goto 516
        t2=atan2(pyr,pxr)
 516    s1=sqrt(1.-c1*c1)
        ct1=cos(t1)                                                             
        st1=sin(t1)                                                             
        ct2=cos(t2)                                                             
        st2=sin(t2)                                                             
        pzr=pr*(c1*c2-s1*s2*ct1)                                                 
        ss=c2*s1*ct1+s2*c1                                                      
        pxr=pr*(ss*ct2-s1*st1*st2)                                               
        pyr=pr*(ss*st2+s1*st1*ct2) 
        eaa=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        ebb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)
        
        CALL LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)
        CALL LORENTZ(-bx,-by,-bz,ebb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)
        px(i1)= pax 
        py(i1)= pay 
        pz(i1)= paz 
        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz

        return 
        end