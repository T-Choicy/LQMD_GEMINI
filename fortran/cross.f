cc	Nucleon-Nucleon elastic and inelatic cross section
c///////////////////////////////////////////////////////////////////////////////
c	nucleon-nucleon cross sections for total, elastic, inelastic scattering
c///////////////////////////////////////////////////////////////////////////////
	FUNCTION CROSStot(is,elab)
	et=elab
	ei=alog10(elab)
	ef=3.+ei
	if(is.eq.1) goto 721
	a=  45.8
	b1= 14.52
	b2=-51.15
	b3= 49.14
	b4=-21.81
	b5= 4.68
	b6=-0.39 
      if(et.le.0.001) then
	CROSStot=2000.
      elseif((et.gt.0.001) .and. (et.le.0.04)) then
      CROSStot=3.414+3.0885/et-0.0011748/et**2
	elseif((et.gt.0.04) .and. (et.le.0.31)) then
      CROSStot=22.-0.011148/et+0.093074/et**2
	elseif((et.gt.0.31) .and. (et.le.0.7)) then
      CROSStot=3.544+53.31*et+0.88737/et
	else
	CROSStot=a+b1*ei+b2*ei**2+b3*ei**3+b4*ei**4+b5*ei**5+b6*ei**6
	endif
	return
 721	a=   0.648
	b1= -0.523
	b2= -0.105
	b3= -0.068
	b4= -0.011
	b5=  0.012
	b6=0.00083
	aa=  41.1
	bb1= 7.04
	bb2=-9.85
	bb3= 2.89
	if(et.le.2.5E-5) then
	CROSStot=20500.
	elseif((et.gt.2.5E-5).and.(et.le.0.4)) then
	af=a+b1*ef+b2*ef**2+b3*ef**3+b4*ef**4+b5*ef**5+b6*ef**6
	CROSStot=(10.**af)*1000.
	else
	CROSStot=aa+bb1*ei+bb2*ei**2+bb3*ei**3
	endif
	return
	END

	FUNCTION CROSSels(is,elab)
	et=elab
	ei=alog10(elab)
	ef=3.+ei
	if(is.ne.0) goto 722
	a=  45.8
	b1= 14.52
	b2=-51.15
	b3= 49.14
	b4=-21.81
	b5= 4.68
	b6=-0.39 
	als=  22.49
	bls1= -7.47
	bls2=-27.77
	bls3=-24.83
	als2=  21.81
	bls11=-17.84
	bls22= 6.95
	bls33=-0.85
      if(et.le.0.001) then
	CROSSels=2000.
      elseif((et.gt.0.001) .and. (et.le.0.04)) then
      CROSSels=3.414+3.0885/et-0.0011748/et**2
	elseif((et.gt.0.04) .and. (et.le.0.3)) then
      CROSSels=22.-0.011148/et+0.093074/et**2
	elseif((et.gt.0.3).and.(et.le.2.)) then
	CROSSels=als+bls1*ei+bls2*ei**2+bls3*ei**3
	else
	CROSSels=als2+bls11*ei+bls22*ei**2+bls33*ei**3
	endif
	if(CROSSels.le.0.) CROSSels=0.
	return
 722	a=   0.648
	b1= -0.523
	b2= -0.105
	b3= -0.068
	b4= -0.011
	b5=  0.012
	b6=0.00083
	als=   24.89
	bls1= -21.63
	bls2=  8.54
	bls3= -1.09
	if(et.le.2.5E-5) then
	CROSSels=20500.
	elseif((et.gt.2.5E-5).and.(et.le.0.3)) then
	af=a+b1*ef+b2*ef**2+b3*ef**3+b4*ef**4+b5*ef**5+b6*ef**6
	CROSSels=(10.**af)*1000.
	else
	CROSSels=als+bls1*ei+bls2*ei**2+bls3*ei**3
	endif
	if(CROSSels.le.0.) CROSSels=0.
	return
	END

	FUNCTION CROSSins(is,elab)
	et=elab
	ei=alog10(elab)
	if(is.eq.1) goto 723
	ain=  22.81 
	bin1= 26.58 
	bin2=-38.62 
	bin3= 25.62 
	bin4= -7.63
	bin5=  0.84
	if(et.le.0.3) then
	CROSSins=0.
	else
	CROSSins= ain+bin1*ei+bin2*ei**2+bin3*ei**3+bin4*ei**4+bin5*ei**5
	endif
	if(CROSSins.le.0.) CROSSins=0.
	return
723	ain=  15.82
	bin1= 26.86
	bin2=-14.7
	bin3= 2.58
	if(et.le.0.3) then
	CROSSins=0.
	else
	CROSSins= ain+bin1*ei+bin2*ei**2+bin3*ei**3
	endif
	if(CROSSins.le.0.) CROSSins=0.

	return
	END

c///////////////////////////////////////////////////////////////////////////////
	    function sigmand(srt)
*     S. Huber, J. Aichelin, Nucl. Phys. A 573(1994)587-625
	sigmand=0.
	s=srt
	if(s.gt.4.5) s=4.5
	if((s.ge.2.014).and.(srt.lt.2.2)) then
	sigmand=1341.*(s-2.014)**2.819
	elseif((s.ge.2.2).and.(s.lt.2.45)) then
	sigmand=18.51-235.2*(s-2.37)**2
	elseif((s.ge.2.45).and.(s.le.4.5)) then
	sigmand=1444.3*s**(-4.957)
	endif

	return
	end
c///////////////////////////////////////////////////////////////////////////////
	    function sigmans(srt,is1,is2)
	sigmans=0.
	s=srt
	if(s.gt.4.5) s=4.5
	IF(is1.eq.is2) THEN
	if((s.ge.2.0158).and.(s.lt.2.35))  then
	sigmans=1.704*(s-2.0158)**2.225
	elseif((s.ge.2.35).and.(s.lt.2.6)) then
	sigmans=1.015*s-2.241
	elseif((s.ge.2.6).and.(s.le.4.5))  then
	sigmans=0.6773-0.1356*abs(s-3.877)**3.014
	endif
	ELSE
	if((s.ge.2.016).and.(s.lt.2.33))  then
	sigmans=1.343*(s-2.016)**2.091
	elseif((s.ge.2.33).and.(s.lt.2.53)) then
	sigmans=0.8031*s-1.752
	elseif((s.ge.2.53).and.(s.lt.3.04)) then
	sigmans=0.3404-0.3272*(s-2.945)**2
	elseif((s.ge.3.04).and.(s.le.4.5))  then
	sigmans=0.5326-0.06395*s
	endif
	ENDIF

	return
	end     
c///////////////////////////////////////////////////////////////////////////////
	    function sigmadd(srt,is1,is2,sigmadd1,sigmadd2)
*	idn=0, nn or pp; 1, np->d(0)+d(+); 2, np->d(-)+d(++)  
	s=srt
	sigmadd=0.
	sigmadd1=0.
	sigmadd2=0.

	if(s.gt.4.5) s=4.5
	IF(is1.eq.is2) THEN
	if((s.ge.2.152).and.(s.lt.2.44))  then
	sigmadd=150.4*(s-2.152)**5.299
	elseif((s.ge.2.44).and.(s.lt.2.63)) then
	sigmadd=2.794*s-6.626
	elseif((s.ge.2.63).and.(s.lt.2.98))  then
	sigmadd=0.9205-2.763*(s-2.898)**2
	elseif((s.ge.2.98).and.(s.lt.3.24))  then
	sigmadd=2.192-0.432*s
	elseif((s.ge.3.24).and.(s.le.4.5))   then
	sigmadd=6.898*s**(-1.417)-0.514
	endif
	ELSE
	if((s.ge.2.152).and.(s.lt.2.412))  then
	sigmadd1=144.8*(s-2.152)**5.083
	elseif((s.ge.2.412).and.(s.lt.2.58)) then
	sigmadd1=2.841*s-6.707
	elseif((s.ge.2.58).and.(s.lt.2.98))  then
	sigmadd1=0.7819-1.777*(s-2.869)**2
	elseif((s.ge.2.98).and.(s.lt.3.26))  then
	sigmadd1=1.894-0.3798*s
	elseif((s.ge.3.26).and.(s.le.4.5))   then
	sigmadd1=8.371*s**(-1.926)-0.2035
	endif

	if((s.ge.2.152).and.(s.lt.2.412))  then
	sigmadd2=239.6*(s-2.152)**5.148
	elseif((s.ge.2.412).and.(s.lt.2.6)) then
	sigmadd2=4.738*s-11.21
	elseif((s.ge.2.6).and.(s.lt.2.98))  then
	sigmadd2=1.395-3.615*(s-2.884)**2
	elseif((s.ge.2.98).and.(s.lt.3.69))  then
	sigmadd2=3.517-0.7197*s
	elseif((s.ge.3.69).and.(s.le.4.5))   then
	sigmadd2=9.358*s**(-1.242)-0.9878
	endif
	sigmadd=sigmadd1+sigmadd2

	ENDIF

	return
	end     
c///////////////////////////////////////////////////////////////////////////////
	    function sigddnn(srt,is1,is2)
	s=srt
	isps=is1+is2
	sigddnn=0.

	if(s.gt.4.5) s=4.5
	if(isps.eq.2) then
	if(is1.eq.2 .or. is2.eq.2) then
	if((s.ge.2.152).and.(s.lt.2.44))  then
	sigddnn=150.4*(s-2.152)**5.299
	elseif((s.ge.2.44).and.(s.lt.2.63)) then
	sigddnn=2.794*s-6.626
	elseif((s.ge.2.63).and.(s.lt.2.98))  then
	sigddnn=0.9205-2.763*(s-2.898)**2
	elseif((s.ge.2.98).and.(s.lt.3.24))  then
	sigddnn=2.192-0.432*s
	elseif((s.ge.3.24).and.(s.le.4.5))   then
	sigddnn=6.898*s**(-1.417)-0.514
	endif
	elseif(is1.eq.1 .or. is2.eq.1) then
	if((s.ge.2.152).and.(s.lt.2.44))  then
	sigddnn=150.4*(s-2.152)**5.299
	elseif((s.ge.2.44).and.(s.lt.2.63)) then
	sigddnn=2.794*s-6.626
	elseif((s.ge.2.63).and.(s.lt.2.98))  then
	sigddnn=0.9205-2.763*(s-2.898)**2
	elseif((s.ge.2.98).and.(s.lt.3.24))  then
	sigddnn=2.192-0.432*s
	elseif((s.ge.3.24).and.(s.le.4.5))   then
	sigddnn=6.898*s**(-1.417)-0.514
	endif
	sigddnn=0.66667*sigddnn
	else
	sigddnn=0.
	endif
	elseif(isps.eq.1) then
	if(is1.eq.0 .or. is2.eq.0) then
	if((s.ge.2.152).and.(s.lt.2.412))  then
	sigddnn=144.8*(s-2.152)**5.083
	elseif((s.ge.2.412).and.(s.lt.2.58)) then
	sigddnn=2.841*s-6.707
	elseif((s.ge.2.58).and.(s.lt.2.98))  then
	sigddnn=0.7819-1.777*(s-2.869)**2
	elseif((s.ge.2.98).and.(s.lt.3.26))  then
	sigddnn=1.894-0.3798*s
	elseif((s.ge.3.26).and.(s.le.4.5))   then
	sigddnn=8.371*s**(-1.926)-0.2035
	endif
	elseif(is1.eq.2 .or. is2.eq.2) then
	if((s.ge.2.152).and.(s.lt.2.412))  then
	sigddnn=239.6*(s-2.152)**5.148
	elseif((s.ge.2.412).and.(s.lt.2.6)) then
	sigddnn=4.738*s-11.21
	elseif((s.ge.2.6).and.(s.lt.2.98))  then
	sigddnn=1.395-3.615*(s-2.884)**2
	elseif((s.ge.2.98).and.(s.lt.3.69))  then
	sigddnn=3.517-0.7197*s
	elseif((s.ge.3.69).and.(s.le.4.5))   then
	sigddnn=9.358*s**(-1.242)-0.9878
	endif
	else
	sigddnn=0.
	endif
	elseif(isps.eq.0) then
	if(is1.eq.1 .or. is2.eq.1) then
	if((s.ge.2.152).and.(s.lt.2.44))  then
	sigddnn=150.4*(s-2.152)**5.299
	elseif((s.ge.2.44).and.(s.lt.2.63)) then
	sigddnn=2.794*s-6.626
	elseif((s.ge.2.63).and.(s.lt.2.98))  then
	sigddnn=0.9205-2.763*(s-2.898)**2
	elseif((s.ge.2.98).and.(s.lt.3.24))  then
	sigddnn=2.192-0.432*s
	elseif((s.ge.3.24).and.(s.le.4.5))   then
	sigddnn=6.898*s**(-1.417)-0.514
	endif
	elseif(is1.eq.0 .or. is2.eq.0) then
	if((s.ge.2.152).and.(s.lt.2.44))  then
	sigddnn=150.4*(s-2.152)**5.299
	elseif((s.ge.2.44).and.(s.lt.2.63)) then
	sigddnn=2.794*s-6.626
	elseif((s.ge.2.63).and.(s.lt.2.98))  then
	sigddnn=0.9205-2.763*(s-2.898)**2
	elseif((s.ge.2.98).and.(s.lt.3.24))  then
	sigddnn=2.192-0.432*s
	elseif((s.ge.3.24).and.(s.le.4.5))   then
	sigddnn=6.898*s**(-1.417)-0.514
	endif
	sigddnn=0.66667*sigddnn
	else
	sigddnn=0.
	endif
	else
	sigddnn=0.
	endif

	return
	end
c///////////////////////////////////////////////////////////////////////////////	
	  function signlk(srt,iso1,iso2,idk,nch1,nch2)
c***** NN --> NLK
c***** srt cm invariant energy, iso1, iso2 input charge numbers of particle 1 and 2
c***** idk -- 0 for K0 and 1 for K+
c***** nch1, nch2 output charge numbers for the reaction 1+2 --> 1'+2', signlk--CS mb
c***** K. Tsushima, A. Sibirtsev, A.W. Thomas, G.Q. Li, PRC 59 (1999) 369.
	s0=6.504
	s2=srt**2
	signlk=0.
	if(s2.le.s0) return
	if(idk.eq.1) goto 11
	if(iso1.eq.0 .and. iso2.eq.0) then   !K^0 production 
	a= 1.879
	b= 2.176
	c= 5.264
	nch1= 0
	nch2= 0
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	a= 2.812
	b= 2.121
	c= 4.893
	nch1= 1
	nch2= 0
	else
	nch1= 0
	nch2= 0
	signlk=0.
	return
	endif
	goto 12

 11	if(iso1.eq.1 .and. iso2.eq.1) then   !K^+ production
	a= 1.879
	b= 2.176
	c= 5.264
	nch1=1
	nch2=0
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	a= 2.812
	b= 2.121
	c= 4.893
	nch1=0
	nch2=0
	else
	nch1= 0
	nch2= 0
	signlk=0.
	return
	endif
 12	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	signlk=a*xb*xc

	  return
	  end 
c///////////////////////////////////////////////////////////////////////////////	
	  function signsk(srt,iso1,iso2,idk,nch1,nch2)
        common/ee/mseed,nseed
c***** NN --> NSK
	s0=6.904
	s2=srt**2
	signsk=0.
	if(s2.le.s0) return
	if(idk.eq.1) goto 11

	if(iso1.eq.0 .and. iso2.eq.0) then   !K^0 production 
	a= 5.321
	b= 2.753
	c= 8.51
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc

	a= 1.466
	b= 2.743
	c= 3.271
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc

	sigt=sig1+sig2
	rsig=sig1/sigt
	if(RANFL(mseed).LE.rsig) then
	nch1= 0
	nch2= 0
	signsk=sig1
	else
	nch1= 1
	nch2= -1
	signsk=sig2
	endif
	return

	elseif(iso1.eq.0 .and. iso2.eq.1) then
	a= 6.31
	b= 2.773
	c= 7.82
	nch1= 1
	nch2= 0
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	a= 11.02
	b= 2.782
	c= 7.674
	nch1= 0
	nch2= 1
	else
	a= 7.079
	b= 2.76
	c= 8.164
	nch1= 1
	nch2= 1
	endif
	goto 12

 11	if(iso1.eq.1 .and. iso2.eq.1) then   !K^+ production
	a= 5.321
	b= 2.753
	c= 8.51
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc

	a= 1.466
	b= 2.743
	c= 3.271
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc

	sigt=sig1+sig2
	rsig=sig1/sigt
	if(RANFL(mseed).LE.rsig) then
	nch1= 1
	nch2= 0
	signsk=sig1
	else
	nch1= 0
	nch2= 1
	signsk=sig2
	endif
	return

	elseif(iso1.eq.1 .and. iso2.eq.0) then
	a= 6.31
	b= 2.773
	c= 7.82
	nch1= 0
	nch2= 0
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	a= 11.02
	b= 2.782
	c= 7.674
	nch1= 1
	nch2= -1
	else
	a= 7.079
	b= 2.76
	c= 8.164
	nch1= 0
	nch2= -1
	endif

 12	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	signsk=a*xb*xc

	  return
	  end 

c///////////////////////////////////////////////////////////////////////////////	
	  function sigdlk(srt,iso1,iso2,idk,nch1,nch2)	
c***** NN --> DLK
	s0=8.085
	s2=srt**2
	sigdlk=0.
	if(s2.le.s0) return
	a= 6.166
	b= 2.842
	c= 1.96
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig=a*xb*xc

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.1) then   !K^0 production 
	nch1= 2
	nch2= 0
	sigdlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sigdlk=sig/3.
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sigdlk=sig/3.
	else
	nch1= 0
	nch2= 0
	sigdlk=0.
	endif
	return

 11	if(iso1.eq.0 .and. iso2.eq.0) then   !K^+ production
	nch1= -1
	nch2= 0
	sigdlk=sig
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sigdlk=sig/3.
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sigdlk=sig/3.
	else
	nch1= 0
	nch2= 0
	sigdlk=0.
	endif

	  return
	  end
c///////////////////////////////////////////////////////////////////////////////
	  function sigdsk(srt,iso1,iso2,idk,nch1,nch2)
        common/ee/mseed,nseed
c***** NN --> DSK
	s0=8.531
	s2=srt**2
	sigdsk=0.
	if(s2.le.s0) return
	a= 10.
	b= 2.874
	c= 2.543
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig=a*xb*xc

	if(idk.eq.1) goto 11
	if(iso1.eq.0 .and. iso2.eq.0) then   !K^0 production
	 if(RANFL(mseed).LE.0.85714) then
	nch1= -1
	nch2= 1
	sigdsk=sig
	 else
	nch1= 0
	nch2= 0
	sigdsk=sig/6.
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sigdsk=sig/6.
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 0
	nch2= 1
	sigdsk=sig/3.
	else
	 if(RANFL(mseed).LE.0.6) then
	nch1= 2
	nch2= 0
	sigdsk=sig/2.
	 else
	nch1= 1
	nch2= 1
	sigdsk=sig/3.
	 endif
	endif
	return

 11	if(iso1.eq.0 .and. iso2.eq.0) then   !K^+ production
	 if(RANFL(mseed).LE.0.6) then
	nch1= -1
	nch2= 0
	sigdsk=sig/2.
	 else
	nch1= 0
	nch2= -1
	sigdsk=sig/3.
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= -1
	sigdsk=sig/3.
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sigdsk=sig/6.
	else
	 if(RANFL(mseed).LE.0.85714) then
	nch1= 2
	nch2= -1
	sigdsk=sig
	 else
	nch1= 1
	nch2= 0
	sigdsk=sig/6.
	 endif
	endif
	
	return
	end
c///////////////////////////////////////////////////////////////////////////////	
	  function sndnlk(srt,iso1,iso2,idk,nch1,nch2)
c***** ND --> NLK
c***** iso1, iso2 are charge numbers of nucleon and delta, respectively.
	s0=6.504
	s2=srt**2
	sndnlk=0.
	if(s2.le.s0) return
	a= 8.337
	b= 2.227
	c= 2.511
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig=a*xb*xc

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.-1) then   !K^0 production 
	nch1= 0
	nch2= 0
	sndnlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sndnlk=sig/3.
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 1
	nch2= 0
	sndnlk=sig/3.
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sndnlk=sig/3.
	else
	nch1= 0
	nch2= 0
	sndnlk=0.
	endif
	return

 11	if(iso1.eq.0 .and. iso2.eq.2) then   !K^+ production
	nch1= 1
	nch2= 0
	sndnlk=sig
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sndnlk=sig/3.
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 0
	nch2= 0
	sndnlk=sig/3.
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sndnlk=sig/3.
	else
	nch1= 0
	nch2= 0
	sndnlk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	  function sndnsk(srt,iso1,iso2,idk,nch1,nch2)
        common/ee/mseed,nseed
c***** ND --> NSK
c***** iso1, iso2 are charge numbers of nucleon and delta, respectively.
	s0=6.904
	s2=srt**2
	sndnsk=0.
	if(s2.le.s0) return
	a= 52.72
	b= 2.799
	c= 6.303
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig=a*xb*xc

	if(idk.eq.1) goto 11
	if(iso1.eq.0 .and. iso2.eq.2) then   !K^0 production 
	nch1= 1
	nch2= 1
	sndnsk=sig
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	nch1= 0
	nch2= 0
	sndnsk=sig/2.
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 1
	sndnsk=sig/3.
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	 if(RANFL(mseed).LE.0.66667) then
	nch1= 0
	nch2= 1
	sndnsk=sig/3.
	 else
	nch1= 1
	nch2= 0
	sndnsk=sig/6.
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	 if(RANFL(mseed).LE.0.66667) then
	nch1= 0
	nch2= 1
	sndnsk=sig/3.
	 else
	nch1= 1
	nch2= 0
	sndnsk=sig/6.
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	sndnsk=sig/6.
	else
	nch1= 0
	nch2= 0
	sndnsk=0.
	endif
	return

 11	if(iso1.eq.1 .and. iso2.eq.-1) then   !K^+ production
	nch1= 0
	nch2= -1
	sndnsk=sig
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	nch1= 1
	nch2= 0
	sndnsk=sig/2.
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	 if(RANFL(mseed).LE.0.66667) then
	nch1= 1
	nch2= -1
	sndnsk=sig/3.
	 else
	nch1= 0
	nch2= 0
	sndnsk=sig/6.
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	 if(RANFL(mseed).LE.0.66667) then
	nch1= 1
	nch2= -1
	sndnsk=sig/3.
	 else
	nch1= 0
	nch2= 0
	sndnsk=sig/6.
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= -1
	sndnsk=sig/3.
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	sndnsk=sig/6.
	else
	nch1= 0
	nch2= 0
	sndnsk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function snddlk(srt,iso1,iso2,idk,nch1,nch2)	
c***** ND --> DLK
c***** iso1, iso2 are charge numbers of nucleon and delta, respectively.
	s0=8.085
	s2=srt**2
	snddlk=0.
	if(s2.le.s0) return
	a= 2.704
	b= 2.303
	c= 5.551
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc      !p+ D++ --> D++ +L+K+

	a= 0.312
	b= 2.11
	c= 2.165
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc      !n+ D+ --> D0 +L+K+

	a= 2.917
	b= 2.35
	c= 6.557
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig3=a*xb*xc      !p+ D+ --> D+ +L+K+

	if(idk.eq.1) goto 11
	if(iso1.eq.0 .and. iso2.eq.2) then   !K^0 production 
	nch1= 2
	nch2= 0
	snddlk=sig1
	elseif(iso1.eq.0 .and. iso2.eq.-1) then
	nch1= -1
	nch2= 0
	snddlk=sig1
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 1
	nch2= 0
	snddlk=sig2
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 2
	nch2= 0
	snddlk=0.75*sig2
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	nch1= 0
	nch2= 0
	snddlk=0.75*sig2
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	snddlk=sig3
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	snddlk=sig3
	else
	nch1= 0
	nch2= 0
	snddlk=0.
	endif
	return

 11	if(iso1.eq.1 .and. iso2.eq.2) then   !K^+ production
	nch1= 2
	nch2= 0
	snddlk=sig1
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	nch1= -1
	nch2= 0
	snddlk=sig1
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 0
	nch2= 0
	snddlk=sig2
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	nch1= 1
	nch2= 0
	snddlk=0.75*sig2
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= -1
	nch2= 0
	snddlk=0.75*sig2
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	snddlk=sig3
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	snddlk=sig3
	else
	nch1= 0
	nch2= 0
	snddlk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function snddsk(srt,iso1,iso2,idk,nch1,nch2)
        common/ee/mseed,nseed
c***** ND --> DSK
c***** iso1, iso2 are charge numbers of nucleon and delta, respectively.
	s0=8.531
	s2=srt**2
	snddsk=0.
	if(s2.le.s0) return
	a= 10.33
	b= 2.743
	c= 8.915
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc      !n+ D++ --> D++ +S- +K+

	a= 2.128
	b= 2.843
	c= 5.986
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc      !p+ D0 --> D+ +S- +K+

	a= 10.57
	b= 2.757
	c= 10.11
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig3=a*xb*xc      !n+ D+ --> D+ +S- +K+

	a= 10.3
	b= 2.748
	c= 9.321
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig4=a*xb*xc      !p+ D++ --> D++ +S0 +K+

	a= 1.112
	b= 2.846
	c= 5.943
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig5=a*xb*xc      !n+ D+ --> D0 +S0 +K+

	a= 10.62
	b= 2.759
	c= 10.2
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig6=a*xb*xc      !p+ D+ --> D+ +S0 +K+

	a= 0.647
	b= 2.83
	c= 3.862
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig7=a*xb*xc      !p+ D+ --> D0 +S+ +K+

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.2) then   !K^0 production 
	nch1= 2
	nch2= 1
	snddsk=sig1
	elseif(iso1.eq.0 .and. iso2.eq.-1) then
	nch1= -1
	nch2= 0
	snddsk=sig4
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	stot=sig1+0.75*sig5+0.75*sig7
	rat1=sig1/stot
	rat2=(sig1+0.75*sig5)/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 1
	 elseif(xx.GT.rat1 .and. xx.LE.rat2) then
	nch1= 0
	nch2= 0
	 else
	nch1= 1
	nch2= -1
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	stot=sig2+sig6+0.75*sig7
	rat1=sig2/stot
	rat2=(sig2+sig6)/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 0
	nch2= 1
	 elseif(xx.GT.rat1 .and. xx.LE.rat2) then
	nch1= 1
	nch2= 0
	 else
	nch1= 2
	nch2= -1
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	stot=0.75*sig2+sig4
	rat1=0.75*sig2/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 1
	nch2= 1
	 else
	nch1= 2
	nch2= 0
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	stot=sig3+0.75*sig5
	rat1=sig3/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 1
	nch2= 1
	 else
	nch1= 2
	nch2= 0
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	stot=sig3+sig5
	rat1=sig3/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 0
	nch2= 1
	 else
	nch1= 1
	nch2= 0
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	stot=0.75*sig2+sig6+sig7
	rat1=0.75*sig2/stot
	rat2=(0.75*sig2+sig6)/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 1
	 elseif(xx.GT.rat1 .and. xx.LE.rat2) then
	nch1= 0
	nch2= 0
	 else
	nch1= 1
	nch2= -1
	 endif
	nch1= 0
	nch2= 0
	snddsk=0.
	endif
	return

 11	if(iso1.eq.0 .and. iso2.eq.-1) then   !K^+ production
	nch1= -1
	nch2= -1
	snddsk=sig1
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	stot=sig3+sig5
	rat1=sig3/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 1
	nch2= -1
	 else
	nch1= 0
	nch2= 0
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	stot=sig3+0.75*sig5
	rat1=sig3/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 0
	nch2= -1
	 else
	nch1= -1
	nch2= 0
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	stot=sig1+0.75*sig5
	rat1=sig1/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= -1
	 else
	nch1= 1
	nch2= 0
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.2) then
	stot=sig4+0.75*sig7
	rat1=sig4/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= 0
	 else
	nch1= 1
	nch2= 1
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	stot=0.75*sig2+sig6+sig7
	rat1=0.75*sig2/stot
	rat2=(0.75*sig2+sig6)/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= -1
	 elseif(xx.GT.rat1 .and. xx.LE.rat2) then
	nch1= 1
	nch2= 0
	 else
	nch1= 0
	nch2= 1
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	stot=sig2+sig6+0.75*sig7
	rat1=sig2/stot
	rat2=(sig2+sig6)/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 1
	nch2= -1
	 elseif(xx.GT.rat1 .and. xx.LE.rat2) then
	nch1= 0
	nch2= 0
	 else
	nch1= -1
	nch2= 1
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	stot=0.75*sig2+sig4
	rat1=0.75*sig2/stot
	xx=RANFL(mseed)
	snddsk=stot
	 if(xx.LE.rat1) then
	nch1= 0
	nch2= -1
	 else
	nch1= -1
	nch2= 0
	 endif
	else
	nch1= 0
	nch2= 0
	snddsk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function cdddlk(srt,iso1,iso2,idk,nch1,nch2)
c***** DD --> DLK

	s0=8.085
	s2=srt**2
	cdddlk=0.
	if(s2.le.s0) return
	a= 1.054
	b= 2.149
	c= 7.969
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc      !D+ + D++ --> D++ +L+K+

	a= 0.881
	b= 2.15
	c= 7.977
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc      !D0 + D++ --> D+ +L+K+

	a= 0.291
	b= 2.148
	c= 7.934
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig3=a*xb*xc      !D0 + D+ --> D0 +L+K+

	if(idk.eq.1) goto 11
	if(iso1.eq.-1 .and. iso2.eq.1) then   !K^0 production 
	nch1= -1
	nch2= 0
	cdddlk=sig1
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 2
	nch2= 0
	cdddlk=sig1/3.
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= 0
	nch2= 0
	cdddlk=sig1/9.
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch1= 0
	nch2= 0
	cdddlk=sig2
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	cdddlk=sig3
	else
	nch1= 0
	nch2= 0
	cdddlk=0.
	endif
	return

 11	if(iso1.eq.1 .and. iso2.eq.2) then   !K^+ production
	nch1= 2
	nch2= 0
	cdddlk=sig1
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch1= -1
	nch2= 0
	cdddlk=sig1/3.
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch1= 1
	nch2= 0
	cdddlk=sig1/9.
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	nch1= 1
	nch2= 0
	cdddlk=sig2
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch1= 0
	nch2= 0
	cdddlk=sig3
	else
	nch1= 0
	nch2= 0
	cdddlk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	  function cdddsk(srt,iso1,iso2,idk,nch1,nch2)
        common/ee/mseed,nseed
c***** DD --> DSK

	s0=8.531
	s2=srt**2
	cdddsk=0.
	if(s2.le.s0) return
	a= 3.532
	b= 2.953
	c= 12.06
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig1=a*xb*xc      !D0 + D++ --> D++ + S- + K+

	a= 7.047
	b= 2.952
	c= 12.05
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig2=a*xb*xc      !D- + D0 --> D- + S- + K+

	a= 2.931
	b= 2.952
	c= 12.03
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig3=a*xb*xc      !D0 + D++ --> D+ + S0 + K+

	a= 5.861
	b= 2.952
	c= 12.04
	xb=(s2/s0-1.)**b
	xc=(s0/s2)**c
	sig4=a*xb*xc      !D- + D+ --> D0 + S- + K+

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.2) then   !K^0 production 
	nch1= 2
	nch2= 1
	cdddsk=sig2
	elseif(iso1.eq.-1 .and. iso2.eq.0) then
	nch1= -1
	nch2= 0
	cdddsk=sig2/2.
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	stot=sig1+sig4
	rat1=sig1/stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= 0
	 else
	nch1= 1
	nch2= 1
	 endif
	elseif(iso1.eq.-1 .and. iso2.eq.2) then
	xx=RANFL(mseed)
	cdddsk=2.*sig1
	 if(xx.LE.0.5) then
	nch1= 1
	nch2= 0
	 else
	nch1= 0
	nch2= 1
	 endif
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	stot=sig1+sig3
	rat1=sig1/stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 1
	 else
	nch1= 0
	nch2= 0
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	stot=sig2/3.+sig2/18.
	rat1=sig2/3./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 1
	 else
	nch1= 0
	nch2= 0
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	stot=sig2/6.+sig2/9.
	rat1=sig2/6./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= 0
	 else
	nch1= 1
	nch2= 1
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	stot=2.*sig3/3.+sig4/6.
	rat1=2.*sig3/3./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 0
	nch2= 1
	 else
	nch1= 1
	nch2= 0
	 endif
	else
	nch1= 0
	nch2= 0
	cdddsk=0.
	endif
	return

 11	if(iso1.eq.-1 .and. iso2.eq.0) then   !K^+ production
	nch1= -1
	nch2= -1
	cdddsk=sig2
	elseif(iso1.eq.1 .and. iso2.eq.2) then
	nch1= 2
	nch2= 0
	cdddsk=sig2/2.
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	stot=sig1+sig3
	rat1=sig1/stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= -1
	 else
	nch1= 1
	nch2= 0
	 endif
	elseif(iso1.eq.-1 .and. iso2.eq.2) then
	xx=RANFL(mseed)
	cdddsk=2.*sig1
	 if(xx.LE.0.5) then
	nch1= 1
	nch2= -1
	 else
	nch1= 0
	nch2= 0
	 endif
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	stot=sig1+sig4
	rat1=sig1/stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 0
	 else
	nch1= 0
	nch2= -1
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	stot=sig2/6.+sig2/9.
	rat1=sig2/6./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= -1
	nch2= 0
	 else
	nch1= 0
	nch2= -1
	 endif
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	stot=sig2/3.+sig2/18.
	rat1=sig2/3./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 2
	nch2= -1
	 else
	nch1= 1
	nch2= 0
	 endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	stot=2.*sig3/3.+sig4/6.
	rat1=2.*sig3/3./stot
	xx=RANFL(mseed)
	cdddsk=stot
	 if(xx.LE.rat1) then
	nch1= 1
	nch2= -1
	 else
	nch1= 0
	nch2= 0
	 endif
	else
	nch1= 0
	nch2= 0
	cdddsk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function spnlk(srt,iso1,iso2,idk)
c***** PN --> LK
	if(idk.EQ.0) then
	s0=1.613
	else
	s0=1.609
	endif
	spnlk=0.
	if(srt.le.s0) return

	if(srt.LE.1.7) then
	sig= 0.9*(srt-s0)/0.091   !mb pion- + p --> L + K0
	else
	sig= 0.09/(srt-1.6)
	endif

	if(idk.eq.1) goto 11
	if(iso1.eq.-1 .and. iso2.eq.1) then
	spnlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	spnlk=sig/2.
	else
	spnlk=0.
	endif
	return

 11	if(iso1.eq.1 .and. iso2.eq.0) then
	spnlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	spnlk=sig/2.
	else
	spnlk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function spnsk(srt,iso1,iso2,idk,nch)
c***** PN --> SK
	s0=1.688
	spnsk=0.
	if(srt.le.s0) return

	sig1= 0.03591*(srt-s0)**0.9541/((srt-1.89)**2+0.01548) + 
     &      0.1594*(srt-s0)**0.01056/((srt-3.)**2+0.9412)           !pion+ + p --> S+ + K+
	sig2= 0.009803*(srt-s0)**0.6021/((srt-1.742)**2+0.006583) + 
     &      0.006521*(srt-s0)**1.4728/((srt-1.94)**2+0.006248)      !pion- + p --> S- + K+
	sig3= 0.05014*(srt-s0)**1.2878/((srt-1.73)**2+0.006455)       !pion+ + n --> S0 + K+
	sig4= 0.003978*(srt-s0)**0.5848/((srt-1.74)**2+0.00667) + 
     &      0.04709*(srt-s0)**2.165/((srt-1.905)**2+0.006358)       !pion0 + p --> S0 + K+

	if(idk.eq.1) goto 11
	if(iso1.eq.-1 .and. iso2.eq.0) then
	nch= -1
	spnsk=sig1
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch= 1
	spnsk=sig2
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch= 1
	spnsk=sig3
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch= 0
	spnsk=sig3
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch= 0
	spnsk=sig4
	else
	spnsk=0.
	endif
	return

 11	if(iso1.eq.1 .and. iso2.eq.1) then
	nch= 1
	spnsk=sig1
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch= -1
	spnsk=sig2
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch= 0
	spnsk=sig3
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch= -1
	spnsk=sig3
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch= 0
	spnsk=sig4
	else
	spnsk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function cpdlk(srt,iso1,iso2,idk)
c***** PD --> LK
c***** K. Tsushima, S.W. Huang, A. Faessler, J. Phys. G 21 (1995) 33-42
	s0=1.613
	cpdlk=0.
	if(srt.le.s0) return

	sig= 0.006545*(srt-s0)**0.7866/((srt-1.72)**2+0.004852)   !mb pion- + D++ --> L + K+

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.-1) then
	cpdlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	cpdlk=0.66667*sig
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	cpdlk=0.33333*sig
	else
	cpdlk=0.
	endif
	return

 11	if(iso1.eq.-1 .and. iso2.eq.2) then
	cpdlk=sig
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	cpdlk=0.66667*sig
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	cpdlk=0.33333*sig
	else
	cpdlk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////	
	  function cpdsk(srt,iso1,iso2,idk,nch)
c***** PD --> SK
	s0=1.688
	cpdsk=0.
	if(srt.le.s0) return

	sig1= 0.004959*(srt-s0)**0.7785/((srt-1.725)**2+0.008147)        !pion- + D++ --> S0 + K+
	sig2= 0.006964*(srt-s0)**0.814/((srt-1.725)**2+0.007713)         !pion0 + D0 --> S- + K+
	sig3= 0.002053*(srt-s0)**0.9853/((srt-1.725)**2+0.005414) + 
     &      0.3179*(srt-s0)**0.9025/((srt-2.675)**2+44.88)             !pion+ + D0 --> S0 + K+
	sig4= 0.01741*(srt-s0)**1.2078/((srt-1.725)**2+0.003777)         !pion+ + D- --> S- + K+

	if(idk.eq.1) goto 11
	if(iso1.eq.1 .and. iso2.eq.-1) then
	nch= 0
	cpdsk=sig1
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch= 1
	cpdsk=0.66667*sig1
	elseif(iso1.eq.0 .and. iso2.eq.-1) then
	nch= -1
	cpdsk=1.5*sig3-0.5*sig4+0.5*sig1
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch= 0
	cpdsk=0.33333*sig4
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch= 1
	cpdsk=sig2
	elseif(iso1.eq.-1 .and. iso2.eq.0) then
	nch= -1
	cpdsk=0.66667*(1.5*sig3-0.5*sig4+0.5*sig1)
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch= 0
	cpdsk=sig3
	elseif(iso1.eq.-1 .and. iso2.eq.2) then
	nch= 1
	cpdsk=sig4
	else
	nch=0
	cpdsk=0.
	endif
	return

 11	if(iso1.eq.-1 .and. iso2.eq.2) then
	nch= 0
	cpdsk=sig1
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch= -1
	cpdsk=0.66667*sig1
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	nch= 1
	cpdsk=1.5*sig3-0.5*sig4+0.5*sig1
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	nch= 0
	cpdsk=0.33333*sig4
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	nch= -1
	cpdsk=sig2
	elseif(iso1.eq.1 .and. iso2.eq.1) then
	nch= 1
	cpdsk=0.66667*(1.5*sig3-0.5*sig4+0.5*sig1)
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	nch= 0
	cpdsk=sig3
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	nch= -1
	cpdsk=sig4
	else
	nch=0
	cpdsk=0.
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	function cnnkkb(srt,iso1,iso2,DMKB0,DMKB,idk,nch1,nch2)
c***** NN --> NNKKB
	cnnkkb= 0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
      DMN= 0.938

	if(idk.eq.1) s0=(2.*DMN+DMK0+DMKB0)**2       !kch=0
	if(idk.eq.2) s0=(2.*DMN+DMK0+DMKB)**2        !kch=-1
	if(idk.eq.3) s0=(2.*DMN+DMK+DMKB0)**2        !kch=1
	if(idk.eq.4) s0=(2.*DMN+DMK+DMKB)**2         !kch=0
	if(s2.le.s0) return

	sig= 1.5*(1.-s0/s2)**3.17*(s0/s2)**1.96

	if(iso1+iso2.eq.3) then
	if(idk.eq.3) then
	cnnkkb=sig
	nch1=1
	nch2=1
	else
	cnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.2) then
	if(idk.eq.1) then 
	cnnkkb=sig
	nch1=1
	nch2=1
	elseif(idk.eq.3) then
	cnnkkb=sig
	nch1=0
	nch2=1
	elseif(idk.eq.4) then
	cnnkkb=sig
	nch1=1
	nch2=1
	else
	cnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.1) then
	if(idk.eq.1) then
	cnnkkb=sig
	nch1=0
	nch2=1
	elseif(idk.eq.2) then
	cnnkkb=sig
	nch1=1
	nch2=1
	elseif(idk.eq.3) then
	cnnkkb=sig
	nch1=0
	nch2=0
	elseif(idk.eq.4) then
	cnnkkb=sig
	nch1=0
	nch2=1
	else
	cnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.0) then
	if(idk.eq.1) then
	cnnkkb=sig
	nch1=0
	nch2=0
	elseif(idk.eq.2) then
	cnnkkb=sig
	nch1=0
	nch2=1
	elseif(idk.eq.4) then
	cnnkkb=sig
	nch1=0
	nch2=0
	else
	cnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-1) then
	if(idk.eq.2) then
	cnnkkb=sig
	nch1=0
	nch2=0
	else
	cnnkkb=0.
	nch1=0
	nch2=0
	endif
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	function cndtkkb(srt,iso1,iso2,DMKB0,DMKB,idk,nch1,nch2)
c***** ND --> N'D'KKB
c***** iso1 -- charge number of nucleon; iso2 -- charge number of delta
c***** nch1 -- charge number of nucleon'; nch2 -- charge number of delta'
      common/ee/mseed,nseed
	cndtkkb=0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
      DMN= 0.938
	DMD= 1.232

	if(idk.eq.1) s0=(DMN+DMD+DMK0+DMKB0)**2       !kch=0
	if(idk.eq.2) s0=(DMN+DMD+DMK0+DMKB)**2        !kch=-1
	if(idk.eq.3) s0=(DMN+DMD+DMK+DMKB0)**2        !kch=1
	if(idk.eq.4) s0=(DMN+DMD+DMK+DMKB)**2         !kch=0
	if(s2.le.s0) return

	sig= 1.5*(1.-s0/s2)**3.17*(s0/s2)**1.96

	if(iso1+iso2.eq.4) then
	if(idk.eq.3) then
	cndtkkb=sig
	nch1=1
	nch2=2
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.3) then
	if(idk.eq.1) then
	cndtkkb=sig
	nch1=1
	nch2=2
	elseif(idk.eq.3) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=1
	 else
	nch1=0
	nch2=2
	 endif
	elseif(idk.eq.4) then
	cndtkkb=sig
	nch1=1
	nch2=2
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.2) then
	if(idk.eq.1) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=1
	 else
	nch1=0
	nch2=2
	 endif
	elseif(idk.eq.2) then
	cndtkkb=sig
	nch1=1
	nch2=2
	elseif(idk.eq.3) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=0
	 else
	nch1=0
	nch2=1
	 endif
	elseif(idk.eq.4) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=1
	 else
	nch1=0
	nch2=2
	 endif
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.1) then
	if(idk.eq.1) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=0
	 else
	nch1=0
	nch2=1
	 endif
	elseif(idk.eq.2) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=1
	 else
	nch1=0
	nch2=2
	 endif
	elseif(idk.eq.3) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=0
	nch2=0
	 else
	nch1=1
	nch2=-1
	 endif
	elseif(idk.eq.4) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=0
	 else
	nch1=0
	nch2=1
	 endif
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.0) then
	if(idk.eq.1) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	elseif(idk.eq.2) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=0
	 else
	nch1=0
	nch2=1
	 endif
	elseif(idk.eq.3) then
	cndtkkb=sig
	nch1=0
	nch2=-1
	elseif(idk.eq.4) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-1) then
	if(idk.eq.1) then
	cndtkkb=sig
	nch1=0
	nch2=-1
	elseif(idk.eq.2) then
	cndtkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	elseif(idk.eq.4) then
	cndtkkb=sig
	nch1=0
	nch2=-1
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-2) then
	if(idk.eq.2) then
	cndtkkb=sig
	nch1=0
	nch2=-1
	else
	cndtkkb=0.
	nch1=0
	nch2=0
	endif
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	function cddnnkkb(srt,iso1,iso2,DMKB0,DMKB,idk,nch1,nch2)
c***** DD --> D'D'KKB
c***** iso1 -- charge number of nucleon; iso2 -- charge number of delta
c***** nch1 -- charge number of nucleon'; nch2 -- charge number of delta'
      common/ee/mseed,nseed
	cddnnkkb=0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
      DMN= 0.938
	DMD= 1.232

	if(idk.eq.1) s0=(DMD+DMD+DMK0+DMKB0)**2       !kch=0
	if(idk.eq.2) s0=(DMD+DMD+DMK0+DMKB)**2        !kch=-1
	if(idk.eq.3) s0=(DMD+DMD+DMK+DMKB0)**2        !kch=1
	if(idk.eq.4) s0=(DMD+DMD+DMK+DMKB)**2         !kch=0
	if(s2.le.s0) return

	sig= 1.5*(1.-s0/s2)**3.17*(s0/s2)**1.96

	if(iso1+iso2.eq.5) then
	if(idk.eq.3) then
	cddnnkkb=sig
	nch1=2
	nch2=2
	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.4) then
	if(idk.eq.1) then
	cddnnkkb=sig
	nch1=2
	nch2=2
	elseif(idk.eq.3) then
	cddnnkkb=sig
	nch1=2
	nch2=1
	elseif(idk.eq.4) then
	cddnnkkb=sig
	nch1=2
	nch2=2
	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.3) then
	if(idk.eq.1) then
	cddnnkkb=sig
	nch1=2
	nch2=1
	elseif(idk.eq.2) then
	cddnnkkb=sig
	nch1=2
	nch2=2
	elseif(idk.eq.3) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=0
	 else
	nch1=1
	nch2=1
	 endif
	elseif(idk.eq.4) then
	cddnnkkb=sig
	nch1=2
	nch2=1
	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.2) then
	if(idk.eq.1) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=0
	 else
	nch1=1
	nch2=1
	 endif
	elseif(idk.eq.2) then
	cddnnkkb=sig
	nch1=2
	nch2=1
	elseif(idk.eq.3) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=-1
	 else
	nch1=1
	nch2=0
	 endif
	elseif(idk.eq.4) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=0
	 else
	nch1=1
	nch2=1
	 endif
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.1) then
	if(idk.eq.1) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=-1
	 else
	nch1=1
	nch2=0
	 endif
	elseif(idk.eq.2) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=0
	 else
	nch1=1
	nch2=1
	 endif
	elseif(idk.eq.3) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	elseif(idk.eq.4) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=-1
	 else
	nch1=1
	nch2=0
	 endif
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.0) then
	if(idk.eq.1) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	elseif(idk.eq.2) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=2
	nch2=-1
	 else
	nch1=1
	nch2=0
	 endif
	elseif(idk.eq.3) then
	cddnnkkb=sig
	nch1=0
	nch2=-1
	elseif(idk.eq.4) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-1) then
	if(idk.eq.1) then
	cddnnkkb=sig
	nch1=0
	nch2=-1
	elseif(idk.eq.2) then
	cddnnkkb=2.*sig
	xx=RANFL(mseed)
	 if(xx.LE.0.5) then
	nch1=1
	nch2=-1
	 else
	nch1=0
	nch2=0
	 endif
	elseif(idk.eq.3) then
	cddnnkkb=sig
	nch1=-1
	nch2=-1
	elseif(idk.eq.4) then
	cddnnkkb=sig
	nch1=0
	nch2=-1
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-2) then
	if(idk.eq.1) then
	cddnnkkb=sig
	nch1=-1
	nch2=-1
	elseif(idk.eq.2) then
	cddnnkkb=sig
	nch1=0
	nch2=-1
	elseif(idk.eq.4) then
	cddnnkkb=sig
	nch1=-1
	nch2=-1
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	elseif(iso1+iso2.eq.-3) then
	if(idk.eq.2) then
	cddnnkkb=sig
	nch1=-1
	nch2=-1
 	else
	cddnnkkb=0.
	nch1=0
	nch2=0
	endif
	endif

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	function cpbkkb(srt,iso1,iso2,DMKB0,DMKB,idk,nch)
c***** PN --> N'KKB
c***** iso1 -- charge number of pion, iso2 -- charge number of nucleon, nch -- charge number of nucleon'
	cpbkkb= 0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
      DMN= 0.938

	if(idk.eq.1) s0=(DMN+DMK0+DMKB0)**2
	if(idk.eq.2) s0=(DMN+DMK0+DMKB)**2
	if(idk.eq.3) s0=(DMN+DMK+DMKB0)**2
	if(idk.eq.4) s0=(DMN+DMK+DMKB)**2
	if(s2.le.s0) return
	sig=1.121*(1.-s0/s2)**1.86*(s0/s2)**2

	if(iso1.eq.1 .and. iso2.eq.1) then
	if(idk.eq.3) then 
	cpbkkb=sig
	nch=1
	else
	cpbkkb=0.
	nch=0
	endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	if(idk.eq.1) then
	cpbkkb=2.*sig
	nch=1
	elseif(idk.eq.3) then
	cpbkkb=sig
	nch=0
	elseif(idk.eq.4) then
	cpbkkb=2.*sig
	nch=1
	else
	cpbkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	if(idk.eq.1) then
	cpbkkb=0.5*sig
	nch=1
	elseif(idk.eq.3) then
	cpbkkb=2.*sig
	nch=0
	elseif(idk.eq.4) then
	cpbkkb=0.5*sig
	nch=1
	else
	cpbkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	if(idk.eq.1) then
	cpbkkb=0.5*sig
	nch=0
	elseif(idk.eq.2) then
	cpbkkb=2.*sig
	nch=1
	elseif(idk.eq.4) then
	cpbkkb=0.5*sig
	nch=0
	else
	cpbkkb=0.
	nch=0
	endif
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	if(idk.eq.1) then
	cpbkkb=2.*sig
	nch=0
	elseif(idk.eq.2) then
	cpbkkb=sig
	nch=1
	elseif(idk.eq.4) then
	cpbkkb=2.*sig
	nch=0
	else
	cpbkkb=0.
	nch=0
	endif
	elseif(iso1.eq.-1 .and. iso2.eq.0) then
	if(idk.eq.2) then
	cpbkkb=sig
	nch=0
	else
	cpbkkb=0.
	nch=0
	endif
	else
	cpbkkb=0.
	nch=0
	endif

	return
	end
c///////////////////////////////////////////////////////////////////////////////
	function cspdkkb(srt,iso1,iso2,DMKB0,DMKB,idk,nch)
c***** PD --> N'KKB
c***** iso1 -- charge number of pion, iso2 -- charge number of delta, nch -- charge number of nucleon'
	cspdkkb=0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
      DMN= 0.938

	if(idk.eq.1) s0=(DMN+DMK0+DMKB0)**2
	if(idk.eq.2) s0=(DMN+DMK0+DMKB)**2
	if(idk.eq.3) s0=(DMN+DMK+DMKB0)**2
	if(idk.eq.4) s0=(DMN+DMK+DMKB)**2
	if(s2.le.s0) return
	sig=1.121*(1.-s0/s2)**1.86*(s0/s2)**2

	if(iso1.eq.1 .and. iso2.eq.1) then
	if(idk.eq.3) then 
	cspdkkb=sig
	nch=1
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.1 .and. iso2.eq.0) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.3) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=1
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.1 .and. iso2.eq.-1) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.2) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=0
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.2) then
	if(idk.eq.3) then 
	cspdkkb=sig
	nch=1
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.1) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.3) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=1
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.0) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.2) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=0
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.0 .and. iso2.eq.-1) then
	if(idk.eq.2) then
	cspdkkb=sig
	nch=0
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.-1 .and. iso2.eq.2) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.3) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=1
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.-1 .and. iso2.eq.1) then
	if(idk.eq.1) then
	cspdkkb=sig
	nch=0
	elseif(idk.eq.2) then
	cspdkkb=sig
	nch=1
	elseif(idk.eq.4) then
	cspdkkb=sig
	nch=0
	else
	cspdkkb=0.
	nch=0
	endif
	elseif(iso1.eq.-1 .and. iso2.eq.0) then
	if(idk.eq.2) then
	cspdkkb=sig
	nch=0
	else
	cspdkkb=0.
	nch=0
	endif
	else
	cspdkkb=0.
	nch=0
	endif

	return
	end
c///////////////////////////////////////////////////////////////////////////////
	function cpnpnkkb(srt,DMKB0,DMKB,idk)
c***** PP --> KKB
	cpnpnkkb=0.
	s2=srt**2
      DMK0=0.498
      DMK= 0.494
	if(idk.eq.1) s0=(DMK0+DMKB0)**2
	if(idk.eq.2) s0=(DMK0+DMKB)**2
	if(idk.eq.3) s0=(DMK+DMKB0)**2
	if(idk.eq.4) s0=(DMK+DMKB)**2
	if(s2.le.s0) return

	cpnpnkkb=2.7*(1.-s0/s2)**0.76   !isospin averaged cross section

	return
	end

	FUNCTION CKELS(SRT,PK)
c **** cross section for kaon-N elastic scattering
	CKELS=0.

	S1=3.+11.5/(1.+exp((PK-1.06)/0.8))
	IF(PK.LT.1.15) S2=18.6*PK**0.86
	IF(PK.GE.1.15 .AND. PK.LT.10.) S2=0.6*PK**2-4.4*PK+25.
	IF(PK.GE.10.) S2=41.
	IF(PK.LT.0.8) S3=0.
	IF(PK.GE.0.8 .AND. PK.LT.1.5) S3=(100.*PK-80.)/7.
	IF(PK.GE.1.5) S3=10.
	CKELS=0.5*(S1+S2-S3)

      RETURN
      END

      FUNCTION CROKB(SRT,PK,IPAR)
c **** cross section for Kb-N interaction in mb
c **** Pk is kaon momentum  in GeV/c
c **** IPAR=1=total, IPAR=2=elastic, IPAR=3=inelastic
c **** ipar=4: k- p --> Lambda pi   from Cugnon et.al.PRC41/4 (1990) 1701
c **** ipar=5: k- p --> Sigma  pi   from Cugnon et.al.PRC41/4 (1990) 1701

      CROKB=0.
      IF(SRT.LE.1.254) RETURN

      S1=22.6*PK**(-1.14)
      IF(PK.GT.0.7)  S1=47.54*PK**(0.94)
      IF(PK.GT.1.1)  S1=69.87*PK**(-3.1)
      IF(PK.GT.1.3)  S1=32.4*PK**(-0.17)
      IF(PK.GT.10.)  S1=22.
      IF(PK.GT.0.03) S2=10.58*PK**(-0.98)
      IF(PK.GT.0.7)  S2=23.*PK**1.2
      IF(PK.GT.1.)   S2=23.*PK**(-2.6)
      IF(PK.GT.1.5)  S2=9.56*PK**(-0.44)
      IF(PK.GT.20.)  S2=2.56
      IF(IPAR.EQ.1)  FSIG=S1
      IF(IPAR.EQ.2)  FSIG=S2
      IF(IPAR.EQ.3) THEN
        FSIG=(S1-S2)
        IF(PK.LT.0.03) FSIG=902.
      ENDIF
      IF(IPAR.EQ.4)THEN
        FSIG=130.
        IF(PK.GT.0.01.AND.PK.LE.0.2) FSIG=12.6*0.2/PK
        IF(PK.GT.0.2.AND.PK.LE.0.9)  FSIG=50.*PK**2-67.*PK+24.
        IF(PK.GT.0.9) FSIG=3.*PK**(-2.6)
c..factor 1.5 from isospin average !
        FSIG=1.5*FSIG
      ENDIF
      IF(IPAR.EQ.5)THEN
c..k- p --> Sigma0 pi0
        SIGMA1=70.
        IF(PK.GT.0.00001.AND.PK.LT.0.2) SIGMA1=10.87*0.2/PK
        IF(PK.GE.0.2) SIGMA1=0.6*PK**(-1.8)
c..k- n --> Sigma0 pi-
        SIGMA2=30.
        IF(PK.GT.0.000001.AND.PK.LT.0.5) SIGMA2=2.954*0.5/PK
        IF(PK.GE.0.5.AND.PK.LT.1.0) SIGMA2=1.2*PK**(-1.3)
        IF(PK.GE.1.0) SIGMA2=1.2*PK**(-2.3)
c..isospin average:
        FSIG=1.5*(SIGMA1+SIGMA2)
      ENDIF
      CROKB=FSIG

      RETURN
      END

	FUNCTION CHYPN(SRT,PK,IPAR)
c **** cross section for hyperon-nucleon interaction in mb
c **** IPAR=1: LN --> LN, IPAR=2: LN --> SN, IPAR=3: SN --> SN, IPAR=4: SN --> LN
	CHYPN=0.

	S1=320.9
	IF(PK.GT.0.1) S1=12.2*PK**(-1.42)
	IF(PK.GT.0.9) S1=14.4*PK**(-0.38)
	S2=0.
	IF(PK.GT.0.636) S2=30.*PK**4.9
	IF(PK.GT.0.8)   S2=5.7*PK**(-2.5)
	IF(IPAR.EQ.1) CHYPN=S1
	IF(IPAR.EQ.2) CHYPN=3.*S2

	SIG1=123.2
	IF(PK.GT.0.15) SIG1=38.*PK**(-0.62)
	SIG2=240.
	IF(PK.GT.0.1)  SIG2=13.5*PK**(-1.25)
	SIG3=240.
	IF(PK.GT.0.1)  SIG3=13.5*PK**(-1.25)
	IF(IPAR.EQ.3)  CHYPN=0.5*(SIG1+SIG2+SIG3)
	SIG4=200.
	IF(PK.GT.0.1)  SIG4=13.2*PK**(-1.18)
	IF(IPAR.EQ.4)  CHYPN=SIG4

      RETURN
      END

	 SUBROUTINE SIGANTIBB(srt,i1,i2,stot,selast,scex,sanni,smeson,
     &                      shyp,slmb,ssgb,sxib)
c***** total and partial cross sections in antinucleon-nucleon/delta and antidelta-nucleon/delta collisions
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        real, parameter :: pCut=1.85     ! momentum above which the PDG parameterization is working (GeV/c)  

	stot=0.
	selast=0.
	scex=0.
	sanni=0.
	smeson=0.
	shyp=0.
	slmb=0.
	ssgb=0.
	sxib=0.
	id1=idanti(i1)
	id2=idanti(i2)
	idsum=id1+id2
	idpul=id1*id2
	icht=iso(i1)+iso(i2)
	idbar=idr(i1)+idr(i2)
	if(idsum.LE.0 .OR. idpul.GT.0)   return    ! antibaryon-baryon collisions

	 if(idbar.EQ.1) then
          I3_Delta=2*icht-1
       if(idsum.GT.0)  I3_Delta=2*icht+1     ! twice I3 projection by Gell-Mann - Nishidjima formula
       end if

       s2=srt**2
       xm1=em(i1)
       xm2=em(i2)

       Etot=(s2 - xm1**2 - xm2**2)/(2.*xm2)           ! total energy of 1-st particle in lab frame
       p=sqrt(Etot**2 - xm1**2)
       v_rel=p/Etot
       pap=0.938*v_rel/sqrt(1.-v_rel**2)              ! momentum of antiproton for the same relative velocity 

	 if(pap .LT. 2.03 ) then
       ! Parameterization from J. Cugnon and J. Vandermeulen, 
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):
       !       selast= 42.3/p**0.54 + 4.3*exp(-(p-1.5)**2)
       selast= 40./pap**0.56 + 5.8*exp(-(pap-1.85)**2)
       else
       selast= 10.2+52.7*pap**(-1.16)+0.125*log(pap)**2-1.28*log(pap)
       end if

       if((idbar.eq.0 .and. icht.eq.0) .OR.
     &    (idbar.eq.1 .and. abs(icht).ne.2)) then
       ! Parameterizations for pbar+p -> nbar+n from J. Cugnon and J. Vandermeulen, 
       ! Annales de Physique (France) 14, 49 (1989);
       ! see also C.B. Dover et al., Prog. Part. Nucl. Phys. 29, 87 (1992):

       if(pap .LE. pCut) then
         if(pap .LT. 0.5) then
           if(pap .LE. 0.1) then
             scex= 0.
           else
             scex= 10.9*(pap-0.1)/pap**1.6
           end if
         else
             scex= 7.1/pap**0.9       ! This is actually not good at p > 3 GeV/c 
                                    ! (above experiment)
         end if
       end if
       else
       scex= 0.
       end if

       if(pap .LE. pCut) then
        if(pap.lt.0.51) then
        if(idbar.eq.0 .and. icht.ne.0 .and. pap.lt.0.382) then
           ! Parameterization for low energy nbar+p annihilation cross section
           ! from T. Armstrong et al., PRD 36, 659 (1987):
           sanni= 41.4 + 29./pap
        else
           ! pbar+p data parameterization by A.L.:            
           sanni= 51.52/pap**0.85 + 0.034/pap**2.94
        end if
        else if(pap.lt.6.34) then                       
        sanni= 88.8/pap**0.4 - 24.2             ! A.L.
        else
        sanni= 38./pap**0.5 + 24./pap**1.1         ! Cugnon
        end if
c        sanni=sanni*( 1. + 4.59*u + 10.6*u**2 + 12.8*u**3 )  in-medium

        if(idbar .LE. 1) then
         if(pap .LE. 0.793) then
           smeson= 0.
         else
           smeson= 30.*(pap-0.793)**1.5/(2.+(pap-0.793)**1.5)
         end if
         if(pap .LE. 1.435) then
           shyp= 0.
         else
           shyp= 3.*(pap-1.435)/(10.+(pap-1.435))
         end if
        else
         smeson= 0.
	   shyp= 0.
        end if
       end if

       if(idbar.eq.0 .and. icht.eq.0) then
          if (s2 .LE. 4.982) then                     !Lambda+LambdaBar Threshold [GeV**2]
             slmb = 0.0
          else if (s2 .LT. 5.071) then
             exc = srt - 2.232 !excess energy
             slmb = 0.0357345*sqrt(exc)+9.00539*exc**1.5                !mb
          else
!             slmb = 0.554417*(s2/4.982-1.)**0.696658*(4.982/s2)**2.33182  !mb
       slmb= max(0.01,0.725872*(s2/4.982-1.)**0.774485*(4.982/s2)**3.35)  !mb
          endif
       else
          slmb = 0.
       endif

       if (idbar.le.1 .and. abs(icht).le.1) then
          if (s2 .le. 5.327) then  !Lambda+Sigma0 Threshold [GeV**2]
             ssgb = 0.
          else
             !ssgb is fitted to Lambda+Sigma0Bar+c.c.
             !ssgb = 2.43292*(s2/5.327-1.)**1.2983*(5.327/s2)**11.1481    !mb
             ssgb = 0.183665*(s2/5.327-1.)**0.436713*(5.327/s2)**1.85006    !mb
          endif

          if(icht .ne. 0) then     ! Multiplicative factors from isospin relations
             if(idbar .eq. 0) then
                ssgb=2.*ssgb
             else if(idbar .eq. 1) then
                if(abs(I3_Delta).eq.3) then
                    ssgb=1.5*ssgb
                else
                    ssgb=0.5*ssgb
                end if
             end if
          end if

          if(srt .LE. 2.63) then        ! Xi+XiBar threshold [GeV]
             sxib=0.
          else
             sxib=0.004  !mb
          end if
          if(idbar .eq. 1) then ! DeltaBar+N or Delta+Nbar collisions --> multiplicative factors from isospin relations
             if(icht .eq. 0) then
                sxib=0.5*sxib
             else if(abs(I3_Delta).eq.1) then
                sxib=0.25*sxib
             else
                sxib=0.75*sxib
             end if
          end if
       else
         ssgb = 0.
	   sxib=0.
       endif

       if(pap .LE. pCut) then
       stot= selast+scex+sanni+smeson+shyp+slmb+ssgb+sxib
       else
      ! Parameterization of total cross sections from PDG, 
      ! L. Montanet et al., PRD 50, 1173 (1994) (see p. 1335);
      ! see also T. Falter et al., PRC 70, 054609 (2004):

       stot=38.4+77.6*pap**(-0.64)+0.26*log(pap)**2-1.2*log(pap)
       end if


       RETURN
       END

  !*****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarN_to_NbarDelta_Integrated 
  ! FUNCTION 
  ! Returns the cross section for pbar n --> nbar delta^- scattering, integrated over the Delta mass.
	  subroutine NbarDelta(sigma,srts)
c*****  use cross section from GiBUU
  !*****************************************************************************
       real ::  srts, dummy2, dummy3, sigma
       real, parameter :: dummy1=1.232

!       call NbarN_to_NbarDelta_dm(dummy2, sigma, dummy3, srts, dummy1)

       end subroutine NbarDelta

  !*****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarDelta_to_NbarN 
  ! FUNCTION
  ! Returns the cross section for nbar delta^- --> pbar n scattering.
  !*****************************************************************************
        subroutine NbDNbN(sigma,srts)
       real ::  srts, dummy1, dummy2, sigma
       real, parameter :: dummy3=1.232

 !      call NbarN_to_NbarDelta_dm(dummy1, dummy2, sigma, srts, dummy3)

       end subroutine NbDNbN

  !*****************************************************************************
  !****s* NbarN_to_NbarDelta/NbarN_to_NbarDelta_dm
  ! NAME
  ! subroutine NbarN_to_NbarDelta_dm(dsdm,dsdm_Integrated,sigma,srts,mdel)
  ! FUNCTION
  ! this subroutine calculates the cross sections for pbar n <--> nbar delta^-
  ! Stores cross sections at first call to a field and then it only returns those field values.
  ! INPUTS
  ! real srts =sqrt(s) in the process
  ! real mdel= mass of delta
  ! RESULT
  ! real :: sigma= total cross section (for reaction with incoming delta) 
  ! real :: dsdm= mass differential cross section as function of mdel and srts 
  !               (for reaction with outgoing delta)  
  ! real :: dsdm_Integrated= mass differential cross section integrated over mass as function 
  !                          of srts (for reaction with outgoing delta) 
  ! NOTES 
  ! Actual cross section is computed by subroutine calculate below.
  !*****************************************************************************
c        subroutine NbarN_to_NbarDelta_dm(dsdm,dsdm_Integrated,sigma,
c     &             srts,mdel)  feng

	 subroutine NbarN_to_NbarDelta_dm

	 common/signbd/xsNbDt(0:200), xsNbNl(0:200)

       real, parameter     :: minimal_mass=1.076
       real, parameter     :: minimal_srts=2.014
c       real, intent(out) :: dsdm, dsdm_Integrated, sigma
c       real, intent(in)  :: srts, mdel
       real  dsdm, dsdm_Integrated, sigma
       real  srts, mdel
       real dsdm_save(0:200,0:150), dsdm_integrated_save(0:200)
       real sigma_save(0:200,0:150)
       integer :: srts_Index, mass_Index
       real :: mass, wurzelS, delta_mass_new 

	  srts=2.3
	  mdel=1.232            !feng

        delta_mass=0.01
	  delta_srts=0.01
	  maxPoints_mass=150
	  maxPoints_srts=200

       dsdm=0.
       dsdm_Integrated=0.
       sigma=0.

!       if(srts.lt.minimal_srts .or. mdel.lt.minimal_mass) return

       Do srts_Index=0, maxPoints_srts
          wurzelS  =minimal_srts   +  Real(srts_Index)   * delta_srts
          Do mass_Index=0, maxPoints_mass
          mass    =  minimal_Mass +  Real(mass_Index)* delta_mass
          call calculate(dsdm,sigma,wurzelS,mass)
          dsdm_save(srts_Index, mass_Index)=dsdm
          sigma_save(srts_Index, mass_Index)=sigma
          End do
        ! Integrate dsDm_save over the mass index
         dsdm_integrated_save(srts_Index) =
     &      Sum(dsdm_save(srts_Index, : )) * delta_mass

       call calculate(dsdm,sigma,wurzelS,mdel)
	 xsNbDt(srts_Index)= dsdm_integrated_save(srts_Index)
	 xsNbNl(srts_Index)= sigma
       End Do
	 
	 return    !feng
c************************************************************
       srts_Index=NINT((srts-minimal_srts)/delta_srts)
       mass_Index=NINT((mdel-minimal_mass)/delta_mass)

       If((srts_Index.gt.maxPoints_srts) .or. 
     &     (mass_Index.gt.maxPoints_mass)) then

        dsdm_integrated=0.

        if (mass_Index .gt. maxPoints_mass) then
        delta_mass_new=delta_mass*10.
        mass_Index=Max(0,  NINT((mdel-minimal_mass)/delta_mass_new))
        If(mass_Index .gt. maxPoints_mass) maxPoints_mass= mass_Index
        else
        delta_mass_new=delta_Mass
        end if

        Do mass_Index=1,maxPoints_mass
        mass    =  minimal_Mass +  Real(mass_Index)* delta_mass_new
        call calculate(dsdm,sigma,srts,mass)
        dsdm_integrated = dsdm_integrated+dsdm * delta_mass_new
        End do

        call calculate(dsdm,sigma,srts,mdel)
       else
       dsdm=dsdm_save(srts_Index, mass_Index)
       dsdm_integrated=dsdm_integrated_save(srts_Index)
       sigma=sigma_save(srts_Index, mass_Index)
       end if

        end subroutine NbarN_to_NbarDelta_dm

   !****************************************************************************
   !****s* NbarN_to_NbarDelta/calculate
   ! NAME
   ! subroutine calculate(dsdm,sigma,srts,mdel)
   ! FUNCTION
   ! Calculate the cross sections for pbar n <--> nbar delta^-
   ! INPUTS
   ! real srts =sqrt(s) in the process
   ! real mdel= mass of delta
   ! RESULT
   ! real :: sigma= total cross section (for reaction with incoming delta) 
   ! real :: dsdm= mass differential cross section as function of mdel and srts 
   !               (for reaction with outgoing delta)  
   !****************************************************************************
        subroutine calculate(dsdm,sigma,srts,mdel)

      real,intent(out) ::  dsdm, sigma
      real,intent(in) :: mdel,srts
      real :: t,qi2,qf2,g,f,dcost,cost!,q,q0
      real :: E1,E2,E3,E4,melem2,dsdodm,dsdo,fac,dmass, mNbarN_to_NbarD2
      integer :: nt, i

	  pi= 3.1415926
        dmass = 1.232
        DMN = 0.938
        DMP = 0.138

       if((srts-2.*DMN-DMP).lt.0. .or. mdel.gt.(srts-DMN)) then
         dsdm  = 0.
         sigma = 0.
         return
       end if

       qi2=srts**2/4.-DMN**2
       qf2=max((srts**2-(mdel+DMN)**2)*(srts**2-(mdel-DMN)**2)/
     &     (4.*srts**2),0.)

       E1 = sqrt(DMN**2 + qi2)
       E2 = E1
       E3 = sqrt(DMN**2+qf2)
       E4 = sqrt(mdel**2+qf2)

       if(mdel .LT. DMN+DMP) then
         f=0.
       else
!        g=FullWidthBaryon(delta,mdel)           
       q=sqrt(max((mdel**2+DMP**2-DMN**2)**2/(4.*mdel**2)-DMP**2,0.))
       q0=sqrt(max((dmass**2+DMP**2-DMN**2)**2/(4.*dmass**2)-DMP**2,0.))
       g=0.118*(q/q0)**3*(0.2**2+q0**2)/(0.2**2+q**2)
       f=1./pi*mdel*g/((mdel**2-dmass**2)**2+(mdel*g)**2)     
       end if

        if(qf2.gt.0.) then
        dcost = min(1.,0.1*DMP**2/(2.*sqrt(qi2*qf2)))
        nt = nint(1./dcost)
        dcost = 1./Real(nt)
        else
        nt = 1
        dcost = 1.
        end if

        dsdm=0.
        sigma=0.
        do i=-nt,nt

        cost=dcost*Real(i)

        t= mdel**2 + DMN**2 - 2.*E4*E2 + 2.*cost*sqrt(qi2*qf2)

!****   Relativistic matrix element^2 :
        melem2 = mNbarN_to_NbarD2(t,DMN,DMN,DMN,mdel)

!****   pbar n -> nbar delta^- :
        dsdodm= 1./(32.*pi*srts**2)*sqrt(qf2)*melem2*2.*mdel*f*0.389

!****   nbar delta^- -> pbar n, factor 1./2. bec. of spin:
        dsdo=   1./(32.*pi*srts**2)*sqrt(qi2)*melem2*0.389/2.

        fac=1.
        if(abs(i).eq.nt) fac=0.5
        dsdm=dsdm+dsdodm/Real(nt)*fac
        sigma=sigma+dsdo/Real(nt)*fac
        end do

        end subroutine calculate

    !***************************************************************************
    !****s* NbarN_to_NbarDelta/mNbarN_to_NbarD2
    ! NAME
    ! real function mNbarN_to_NbarD2 (t, m_1, m_2, m_3, m_4)
    ! FUNCTION 
    ! Matrix element squared and averaged over initial and summed over 
    ! final spins for pbar_1 n_2 -> nbar_3 delta_4^- 
    ! within OPEM with effective nucleon and delta masses.
    ! INPUTS
    !  * t=(p1-p3)^2      -- Mandelstam variable (GeV^2)
    !  * m_1,m_2,m_3,m_4  -- effective masses of particles (GeV)
    !***************************************************************************
        real function mNbarN_to_NbarD2(t, m_1, m_2, m_3, m_4)

      real, intent(in) :: t,m_1,m_2,m_3,m_4
      real :: dmass,qr2t,q2t,zt,ft
      ! Coupling constants from V. Dmitriev et al, NPA 459, 503 (1986):
      real, parameter :: fps=2.202d0,fp=1.008d0
      real, parameter :: lambda= 0.50d0  ! readjusted by A.L.
      real, parameter :: kappa =0.2d0, mPi= 0.138

      dmass = 1.232


      ft=(lambda**2-mpi**2)/(lambda**2-t)       ! From V. Dmitriev et al
      !ft=0.72/(1.+(-t+mpi**2)/(4.73*mpi**2)) + 0.28   ! From R. Armenteros & B. French at low energies

      qr2t=(dmass**2+m_2**2-t)**2/(4.d0*dmass**2)-m_2**2
      q2t=(m_4**2+m_2**2-t)**2/(4.d0*m_4**2)-m_2**2
      zt=(qr2t+kappa**2)/(q2t+kappa**2)


      mNbarN_to_NbarD2 = 2.* (fp*fps*ft**2)**2*zt/mpi**4 /(t-mpi**2)**2
     &                   * (m_1+m_3)**2 * (t-(m_1-m_3)**2)
     &             * (t-(m_4-m_2)**2) * ((m_4+m_2)**2-t)**2/6./m_4**2
  
        end function mNbarN_to_NbarD2