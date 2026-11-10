**********************************************************************************	                    
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE ocu(xx,yy,zz,pxx,pyy,pzz,ocuph)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00
        common/skyrme/al,sal,t0,t1,arang,vy,width
        dimension pp1(6),pp2(6),epp1(6),epp2(6)
        widthp=2.*0.1973289**2/al                                               
        sgma=sqrt(widthp)                                                       
        widthr=al/2.                                                            
        all=sqrt(widthr)                                                        
        ph2=0.045                                                               
        rh2=6.89                                                                
        sum=0.0                                                                 
        do 10 i=1,np                                                            
        pp1(1)=-(ph2-pxx+px(i))/sgma                                            
        pp1(2)=-(ph2-pyy+py(i))/sgma                                            
        pp1(3)=-(ph2-pzz+pz(i))/sgma                                            
        pp1(4)=-(rh2-xx+x(i))/all                                               
        pp1(5)=-(rh2-yy+y(i))/all                                               
        pp1(6)=-(rh2-zz+z(i))/all                                               
        pp2(1)=(ph2+pxx-px(i))/sgma                                             
        pp2(2)=(ph2+pyy-py(i))/sgma                                             
        pp2(3)=(ph2+pzz-pz(i))/sgma                                             
        pp2(4)= (rh2+xx-x(i))/all                                               
        pp2(5)= (rh2+yy-y(i))/all                                               
        pp2(6)= (rh2+zz-z(i))/all                                               
        epp=1.0                                                                 
        do 20 k=1,6                                                             
        aa=pp1(k)                                                               
        if(aa.le.0.) then                                                       
        call erf(-aa,epp1(k))                                                   
        else                                                                    
        call erf(aa,epp1(k))                                                    
        epp1(k)=-epp1(k)                                                        
        endif
        bb=pp2(k)                                                               
        if(bb.le.0.) then                                                       
        call erf(-bb,epp2(k))                                                   
        epp2(k)=-epp2(k)                                                        
        else                                                                    
        call erf(bb,epp2(k))                                                    
        endif                                                                  
        epp=epp*(epp1(k)+epp2(k))                                               
20      continue                                                                
        sum=sum+epp                                                             
10      continue                                                                
        ocuph=sum/2.**6                                                         
        return                                                                  
        end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE rcollbb
c***** baryon/antibaryon-baryon collisions
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	  common/CLWF3/RTIME0
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
        common/deu2/px1,py1,pz1,px2,py2,pz2                                     
	  common/poli/iflagcro,iflagpau,ipoli
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/densi1/denstl(500),dens1(500),velcm(3)
	  common/medium/mdel,mdin,nsym,npop,nkaon
	common/SURF/ VANISH(500) 

      ncoll=0
	ncolln=0
	ncollp=0
      nblock=0 
	nblkn=0
	nblkp=0
	knnp=0

      BMAX0=1.323142    !fm
      SIG0= 55.         !mb
c----------------------
	xpi=31.415926
	pi=3.141592654
	hb=0.197327
	xeff=10./3.*gtau/(1.5*pi**2*rho0*hb**3)**1.66667

c************* medium	effect in nucleon-nucleon collisions *************
        call densitn
c***********************************************************************

        do 21 i1=2,np
		if(VANISH(I1).EQ.1)	cycle
	  x1=x(i1)
        y1=y(i1)
        z1=z(i1)
        ip=i1-1

        do 22 i2=1,ip
		if(VANISH(I2).EQ.1)	cycle
        if(il(i1).eq.i2.and.il(i2).eq.i1) goto 22
        if(iflag1.eq.1.and.id(i2).eq.1.and.id(i1).eq.1)   goto 22
        if(iflag1.eq.1.and.id(i2).eq.51.and.id(i1).eq.51) goto 22

				if(iflag1.eq.3)	then
		ICtrig=0
		do IIII=1,NP-NTA
			if(VANISH(IIII).EQ.1)	cycle
			if(ID(IIII).NE.1) then
			ICtrig=1
				if(RTIME.LT.RTIME0) then
				RTIME0=RTIME
				write(*,*) rtime0
				endif
		goto 995
			endif
		enddo
995		continue
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.51.and.id(i1).eq.51) goto 22                   
				endif

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
        e1=sqrt(em1*em1+pz1*pz1+px1*px1+py1*py1)  
        e2=sqrt(em2*em2+pz2*pz2+px2*px2+py2*py2) 
        bx=(px1+px2)/(e1+e2)
        by=(py1+py2)/(e1+e2)
        bz=(pz1+pz2)/(e1+e2)
	  velcm(1)=bx
	  velcm(2)=by
	  velcm(3)=bz
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
        s2=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2
        if(s2.lt.3.556896) goto 22
 	  srt=sqrt(s2)
	  call  multcoll(x1,y1,z1,px1,py1,pz1,em1,
     &                 x2,y2,z2,px2,py2,pz2,em2,i1,i2,knnp)
	  if(knnp.NE.0) goto 22

	is=1
	if(iso(i1).eq.iso(i2)) is=0

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

c*****	elastic NN scattering
	if(mdel.eq.0) then
	 xmdel=1.
	else
	 cms11= 4.*cmeps*em1*alog(cmeps*pr**2+1.)/(cmeps*pr**2+1.)*
     &      (cmnn*dens1(i1)+cmnp*(denstl(i1)-dens1(i1)))/rho0
	 rff11= 1./(1.+cms11)
	 esf11= rff11*em1
	 cms22= 4.*cmeps*em2*alog(cmeps*pr**2+1.)/(cmeps*pr**2+1.)*
     &      (cmnn*dens1(i2)+cmnp*(denstl(i2)-dens1(i2)))/rho0
	 rff22= 1./(1.+cms22)
	 esf22= rff22*em2
	 rrms= rff11*rff22*(em1+em2)/(esf11+esf22)
	 xmdel=rrms*rrms                           !1.-0.2*(dens1(i1)+dens1(i2))/rho0
	endif

	if(mdin.eq.0) then
	xmdin=1.
	else
	xmdin=xmdel
	endif

c****************get the Elab to determine the minimum approaching distance
	 elab=(s2-(em1+em2)**2)/2./em2  ! GeV
	if(elab.le.0.) goto 22
	 sigels=xmdel*CROSSels(is,elab)

	 idplus=idr(i1)+idr(i2)
	 idplan=idanti(i1)+idanti(i2)

	 if(idplan.EQ.0) then

        if(idplus.EQ.0) then
	  sigref=xmdel*CROSStot(is,elab)
	  dint=sqrt(sigref/xpi)
	  if(sigels.GT.(xmdel*SIG0)) then
	  dint=sqrt(sigels/xpi)
	  sigref=sigels
	  endif
	  else
        call collres(srt,i1,i2,sigels,xmdin)               ! collisions of resonances
	  goto 22
	  endif
	 else
	 goto 22
	 endif


        if(bb.gt.dint) goto 22
        e1r=sqrt(em(i1)*em(i1)+pr*pr) 
        e2r=sqrt(em(i2)*em(i2)+pr*pr) 
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.                                         
        if(dd.lt.abs(dzz))   goto 22

        ncoll=ncoll+1 
	  if(iso(i1).eq.0 .and. iso(i2).eq.0) then
	  ncolln=ncolln+1
	  elseif(iso(i1).eq.1 .and. iso(i2).eq.1) then
	  ncollp=ncollp+1
 	  endif

	idp=0
	ikn=0
        call prodrsn(pxr,pyr,pzr,ppx,ppy,ppz,srt,i1,i2,
     &               sigels,sigref,idp,ims,ikn,xmdin)             !! resonances (delta, N*), strangeness production and S-state pion by NN collisions
	 if(idp.NE.0) then
	 pix=pxms(ims)
	 piy=pyms(ims)
	 piz=pzms(ims)
 	 endif
	 if(ikn.NE.0) goto 22

       if(idp.EQ.0) then
        eaa=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        ebb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)
        call LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,ebb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)
        px(i1)= pax 
        py(i1)= pay 
        pz(i1)= paz 
        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz
	 else
        e11=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        e22=sqrt(em(i2)*em(i2)+ppx*ppx+ppy*ppy+ppz*ppz)
	  ep0=sqrt(emson(ims)*emson(ims)+pix*pix+piy*piy+piz*piz)
        call LORENTZ(-bx,-by,-bz,e11,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,e22,ppx,ppy,ppz,pb0,pbx,pby,pbz)
        call LORENTZ(-bx,-by,-bz,ep0,pix,piy,piz,ep1,pnx,pny,pnz)
        px(i1)= pax
        py(i1)= pay
        pz(i1)= paz

        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz

	  pxms(ims)=pnx
	  pyms(ims)=pny
	  pzms(ims)=pnz
	  goto 22
	 endif

c	pauli blocking for nn scattering
        if(iflagpau.le.0) call paulin(i1)
        if(iflagpau.eq.1) call paur34(i1,i2)
        if(iflagpau.eq.2) call pauover(i1)
	  if(iflagpau.eq.3) call paurho(i1)
	  if(ntag.eq.-1) goto 32
        if(ntag.eq.0)  goto 222

 32 	nblock=nblock+1
	if(iso(i1).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
	if(iso(i2).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
C
        pz(i1)=pz1                                                         
        px(i1)=px1                                                         
        py(i1)=py1                                                         
        pz(i2)=pz2                                                         
        px(i2)=px2                                                         
        py(i2)=py2       
        goto 22  
222     continue 
        if(iflag6.eq.1)call dtest(i1,i2,em(i1),em(i2))
        il(i1)=i2
        il(i2)=i1
        id(i1)=id(i1)+1                                                    
        id(i2)=id(i2)+1                                                    
22   	   continue
21     continue

	  return
        end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE rcollmeson
c***** collisions of meson-meson or meson-baryon
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/bb/np,nzp,ntag,x00
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/densi1/denstl(500),dens1(500),velcm(3)
	common/SURF/ VANISH(500)

	  BMAXP= 0.892    !pion - baryon
	  SIGP = 25.      !inelastic collisions mb
	  BETA = 0.4      !eta/rho/omega - baryon
	  SIET = 5.
        BKAON= 1.38     !kaon - baryon
	  SIGKN= 60.
        BANTK= 2.52     !antikaon - baryon
	  SIGTK= 200.

	 do 21 i1=1,np

	  if(idms(i1).LE.0)      goto 21
	  x1=xms(i1) 
        y1=yms(i1)                                                                
        z1=zms(i1) 

        px1=pxms(i1)
        py1=pyms(i1)
        pz1=pzms(i1) 
	  em1=emson(i1)
	  if(em1 .LT. 0.01)      goto 21
        do 22 i2=1,np            ! meson-baryon collision
		if(VANISH(I2).EQ.1)	cycle
	  if(idanti(i2).NE.0) cycle
        x2=x(i2)
        y2=y(i2)
        z2=z(i2)
        r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
        px1=pxms(i1)                                                              
        px2=px(i2)                                                              
        py1=pyms(i1)                                                              
        py2=py(i2)   
        pz1=pzms(i1) 
        pz2=pz(i2)
        em2=em(i2)
        e1=sqrt(em1*em1+pz1*pz1+px1*px1+py1*py1) 
        e2=sqrt(em2*em2+pz2*pz2+px2*px2+py2*py2) 
        bx=(px1+px2)/(e1+e2) 
        by=(py1+py2)/(e1+e2) 
        bz=(pz1+pz2)/(e1+e2)       
	  velcm(1)=bx
	  velcm(2)=by
	  velcm(3)=bz
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
        s2=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2
	  if(s2.lt.1.1664) goto 22
 	  srt=sqrt(s2)

	if(idms(i1) .EQ. 1) then
	dint=   BMAXP
	sigref= SIGP
	elseif(idms(i1) .EQ. 2) then
	dint=   BETA
	sigref= SIET
	elseif(idms(i1) .EQ. 3) then
	dint=   BKAON
	sigref= SIGKN
	elseif(idms(i1) .EQ. 4) then
	dint=   BANTK
	sigref= SIGTK
	elseif(idms(i1) .EQ. 8) then
	dint=   BETA
	sigref= SIET
	elseif(idms(i1) .EQ. 9) then
	dint=   BETA
	sigref= SIET

	else
	dint=   0.
	sigref= 0.
	endif

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
        if(bb.gt.dint) goto 22
        e1r=sqrt(emson(i1)*emson(i1)+pr*pr)
        e2r=sqrt(em(i2)*em(i2)+pr*pr)
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.
        if(dd.lt.abs(dzz))   goto 22
	
	  if(idms(i1).EQ.1) then
	  kssb=0
	  call prkpnby(srt,i1,i2,sigref,kssb)            ! kaon production by pion-baryon collisions
	  if(kssb.ne.0)    goto 22
	  kkb=0
	  call prpnkkb(srt,i1,i2,sigref,kkb)             ! kaon-antikaon production by pion-baryon collisions
	  if(kkb.ne.0) goto 22
	  call stranex(pxr,pyr,pzr,srt,i1,i2,sigref)	 ! strange exchange reactions
	  if(idms(i1) .EQ. 1)   goto 22

	  elseif(idms(i1).EQ.2) then                     ! eta-nucleon collisions
	  kssb=0
	  call prkpnby(srt,i1,i2,sigref,kssb)            ! kaon production by eta-baryon collisions
	  if(kssb.ne.0) goto 22

	  elseif(idms(i1).EQ.3) then
	  call elaskaon(pxr,pyr,pzr,srt,i1,i2,sigref)    ! elastic scattering of kaon and baryon
	  if(idr(i2).eq.4 .OR. idr(i2).eq.5)
     &  call prkknhy(srt,i1,i2,sigref)                 ! kaon-hyperon collisions
	  elseif(idms(i1).EQ.4) then
	  call absantikn(pxr,pyr,pzr,srt,i1,i2,sigref)   ! absorption of antikaon by baryon
	  else
	  goto 22
 	  endif

        eaa=sqrt(emson(i1)*emson(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        ebb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)
        
        call LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,ebb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)
        pxms(i1)= pax
        pyms(i1)= pay
        pzms(i1)= paz
        px(i2)  = pbx
        py(i2)  = pby
        pz(i2)  = pbz
 22	  continue

         do 25 i2=1,np    ! pion-pion collision
	  if(idms(i2).EQ.0) goto 25
        x2=xms(i2) 
        y2=yms(i2) 
        z2=zms(i2) 
        r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)                                                             
        px1=pxms(i1)                                                              
        px2=pxms(i2)                                                              
        py1=pyms(i1)                                                              
        py2=pyms(i2)   
        pz1=pzms(i1) 
        pz2=pzms(i2)
	  em1=emson(i1)
        em2=emson(i2)
        e1=sqrt(em1*em1+pz1*pz1+px1*px1+py1*py1)
        e2=sqrt(em2*em2+pz2*pz2+px2*px2+py2*py2)
        bx=(px1+px2)/(e1+e2) 
        by=(py1+py2)/(e1+e2) 
        bz=(pz1+pz2)/(e1+e2)       
	  velcm(1)=bx
	  velcm(2)=by
	  velcm(3)=bz
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
        s2=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2
	  if(s2.lt.0.988) goto 25
 	  srt=sqrt(s2)

	if(idms(i1).EQ.1) then
	dint=   BMAXP
	sigref= SIGP
	else
	dint=   0.
	sigref= 0.
	endif

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
        c00=dr*dr-dzz*dzz                 
	  if(c00.le.0.) c00=0.
        bb=sqrt(c00)
        if(bb.gt.dint) goto 25
        e1r=sqrt(emson(i1)*emson(i1)+pr*pr)
        e2r=sqrt(emson(i2)*emson(i2)+pr*pr)
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.                                         
        if(dd.lt.abs(dzz))   goto 25
	
	if(idms(i1).EQ.1.and.idms(i2).EQ.1) call pnpnkkb(srt,i1,i2,sigref)   ! kaon-antikaon production by pion-pion collisions

 25	  continue
 21 	 continue

    	  return
     	  end
c///////////////////////////////////////////////////////////////////////////////
      	SUBROUTINE resdecay
c*****    production of mesons from decays of resonances
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
	  common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &  potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/bb/np,nzp,ntag,x00
	  common/poli/iflagcro,iflagpau,ipoli
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
	common/SURF/ VANISH(500)
      
	  AVMASS = 0.938868
        PIMASS = 0.137265
	  ETMASS =  0.547
	  HB= 0.197327
	  PI= 3.1415926
	  IWD=1      ! 0 for cons. decay width, 1 mom. dep.
 
       do 11 i1=1,np
		if(VANISH(I1).EQ.1)	cycle

	  x1=x(i1)
        y1=y(i1)
        z1=z(i1)
        px1=px(i1)             
        py1=py(i1) 
        pz1=pz(i1)
        em1=em(i1)
	  is1=iso(i1)
	  idrs=idr(i1)
	 if(idrs .LE. 0)  cycle
	 if(idrs .GE. 4)  cycle

	  e1=sqrt(em1**2+px1**2+py1**2+pz1**2)
        bx = px1/e1                                                    
        by = py1/e1                                                    
        bz = pz1/e1
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))

	  DMASS=em1
        AUX = 0.25*(DMASS**2-AVMASS**2-PIMASS**2)**2
     &        -(AVMASS*PIMASS)**2
        if (AUX .GT. 0.) then
        QAVAIL = SQRT(AUX / DMASS**2)
        else
        QAVAIL = 0.0000001
        endif
	  GDN=WIDTHDN(QAVAIL,idrs,IWD)
	  if(GDN.lt.0.00001) then
	  pdecay=0.
	  else
	  tau=gam*HB/GDN
	  pdecay=1.-xexp(-(dt/tau))
	  endif
	  xx=RANFL(nseed)
	  if(xx .GT. pdecay)  goto 11
	  x2=RANFL(mseed)

      icm=0.
 88	icm=icm+1
	if(idms(icm).NE.0) goto 88
      if(icm.GT.np) STOP 'TOO MANY MESONS IN resdecay !!'
	 RRMASS= AVMASS
       WMASS = PIMASS

	 if(idrs .eq. 1) then
	  nidm = 1
	  if(iso(i1).eq.2) then
	  ncms = 1
	  ncnl = 1
	  elseif(iso(i1).eq.1) then
	  if(x2.le.0.3333) then
	  ncms = 1
	  ncnl = 0
	  else
	  ncms = 0
	  ncnl = 1
	  endif
	  elseif(iso(i1).eq.0) then
	  if(x2.le.0.3333) then
	  ncms = -1
	  ncnl = 1
	  else
	  ncms = 0
	  ncnl = 0
	  endif
	  elseif(iso(i1).eq.-1) then
	  ncms = -1
	  ncnl = 0
	  endif 
	 elseif(idrs .eq. 2) then
	  nidm = 1
	  ncms = 0
	  ncnl = iso(i1)
	  if(x2.GT.0.3333) then
	  ncms = 2*iso(i1)-1
	  ncnl = 1-iso(i1)
	  endif
	 else     !decay of N*(1535)
       RNX = RANFL(nseed)
       if(RNX .GT. 0.55) then
	 WMASS = ETMASS
       ncms = 0
       nidm = 2
	 ncnl = iso(i1)
	 else
	 WMASS = PIMASS
	 nidm = 1
	 ncms = 0
	 ncnl = iso(i1)
	 if(x2.GT.0.3333) then
	 ncms = 2*iso(i1)-1
	 ncnl = 1-iso(i1)
       endif
	 endif
	 endif

	  S=em1**2
        PA2 = 1 /(4.*S)*(S**2+(RRMASS)**4+WMASS**4- 2.*S*RRMASS**2-
     &        2.*S*WMASS**2 - 2.*WMASS**2*RRMASS**2)
        if(PA2.LE.0)  return
	  PA = SQRT(PA2)

	  idms(icm)= nidm
	  icp(icm) = ncms
	  idr(i1) = 0
	  iso(i1) = ncnl
	  if(ncnl.eq.0) isy(i1)=1
	  if(abs(ncnl).eq.1) isy(i1)=-1

        c1=1.-2.*RANFL(mseed)
        t1=2.*PI*RANFL(mseed)
        s1=sqrt(1.-c1*c1)
        ct1=cos(t1)
        st1=sin(t1)
        pzr=PA*c1 
        pxr=PA*s1*ct1
        pyr=PA*s1*st1

        ea=sqrt(RRMASS**2+PA2)
        eb=sqrt(WMASS**2+PA2)
        call LORENTZ(-bx,-by,-bz,ea,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,eb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)

        px(i1) = pax
        py(i1) = pay
        pz(i1) = paz
	  em(i1) = RRMASS
	  ep = sqrt(em(i1)**2+px(i1)**2+py(i1)**2+pz(i1)**2)
	  x(i1) = x(i1) + px(i1)*dt/ep
	  y(i1) = y(i1) + py(i1)*dt/ep
	  z(i1) = z(i1) + pz(i1)*dt/ep

        pxms(icm) = pbx
        pyms(icm) = pby
        pzms(icm) = pbz
	  xms(icm)  = x1
	  yms(icm)  = y1
	  zms(icm)  = z1
	  emson(icm)= WMASS

	  ep = sqrt(emson(icm)**2+pxms(icm)**2+pyms(icm)**2+pzms(icm)**2)
	  xms(icm) = xms(icm) + pxms(icm)*dt/ep
	  yms(icm) = yms(icm) + pyms(icm)*dt/ep
	  zms(icm) = zms(icm) + pzms(icm)*dt/ep

        if(iflagpau .eq. 1)  call paur34(i1,0)
	  if(ntag. eq. -1)  then
	  idr(i1) = idrs 
	  iso(i1) = is1
        px(i1)  = px1
        py(i1)  = py1
        pz(i1)  = pz1
        em(i1)  = em1

	  idms(icm) = 0
	  icp(icm)  = 0
        pxms(icm) = 0.
        pyms(icm) = 0.
        pzms(icm) = 0.
	  xms(icm)  = 0.
	  yms(icm)  = 0.
	  zms(icm)  = 0.
	  emson(icm)= 0.
	  endif

 11    continue
 
	  return
	  end

c///////////////////////////////////////////////////////////////////////////////
	  SUBROUTINE mesonabs
c*****    reabsorption of mesons and baryons
c///////////////////////////////////////////////////////////////////////////////
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
	  common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &  potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/bb/np,nzp,ntag,x00
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
	common/SURF/ VANISH(500)

	 xpi=31.415926
       do 21 i1=1,np
		if(VANISH(I1).EQ.1)	cycle
	  x1=x(i1)
        y1=y(i1)
        z1=z(i1)
        px1=px(i1)             
        py1=py(i1) 
        pz1=pz(i1)
        em1=em(i1)
        e1=sqrt(em1*em1+px1*px1+py1*py1+pz1*pz1)

	  idrs=idr(i1)
	  if(idrs .GE. 1)  goto 21

	  do 22 i2=1, np
	  if(idms(i2) .LE. 0) goto 22     ! pion and eta
	  if(idms(i2) .GE. 3) goto 22     ! pion and eta
	  x2=xms(i2)                                                                
        y2=yms(i2)                                                                
        z2=zms(i2)                                                                                                
        px2=pxms(i2)
        py2=pyms(i2)
        pz2=pzms(i2)
        em2=emson(i2)
        e2=sqrt(em2*em2+px2*px2+py2*py2+pz2*pz2)

	  LB1=0
	  LB2=0
	  XMAXD =0.
	  XMAXN1=0.
	  XMAXN2=0.

	if(iso(i1).eq.1)   LB1=1
	if(iso(i1).eq.0)   LB1=2
	if(icp(i2).eq.-1)  LB2=3
	if(icp(i2).eq.0)   LB2=4
	if(icp(i2).eq.1)   LB2=5
	if((LB1.eq.0).or.(LB2.eq.0))  goto 22
      if(LB1*LB2.EQ.5 .OR. LB1*LB2.EQ.6)  XMAXD= 200.   !PION(+)+PROTON-->DELTA++,PION(-)+NEUTRON-->DELTA(-)
      if(LB1*LB2.EQ.3 .OR. LB1*LB2.EQ.10) XMAXD= 66.7   !PION(-)+PROTON-->DELTA0,PION(+)+NEUTRON-->DELTA+
      if(LB1*LB2.EQ.4 .OR. LB1*LB2.EQ.8)  XMAXD=133.3   !PION0+PROTON-->DELTA+,PION0+NEUTRON-->DELTA0
	  if(LB1*LB2.EQ.3 .OR. LB1*LB2.EQ.10)  then
	  XMAXN1= 24.                                     !PION(-)+PROTON-->N*0(1440), PION(+)+NEUTRON-->N*+(1440)
	  XMAXN2= 32.                                     !PION(-)+PROTON-->N*0(1535), PION(+)+NEUTRON-->N*+(1535)
	  endif
	  if(LB1*LB2.EQ.4 .OR. LB1*LB2.EQ.8)   then
	  XMAXN1= 12.                                     !PION0+NEUTRON-->N*0(1440), PION0+PROTON-->N*+(1440)
	  XMAXN2= 16.                                     !PION0+NEUTRON-->N*0(1535), PION0+PROTON-->N*+(1535)
	  endif
	 if(idms(i2) .EQ. 2)  then                        ! ETA+P/N
	 XMAXD=  0.
	 XMAXN1= 0.
	 XMAXN2= 12.
	 endif

        bx=(px1+px2)/(e1+e2)
        by=(py1+py2)/(e1+e2)
        bz=(pz1+pz2)/(e1+e2)
        s2=(e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2
	  if(s2.lt.1.1664)  goto 22
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
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
	  if(pr .LE. 0.)  goto 22                                        
        dzz=(pxr*dx+pyr*dy+pzr*dz)/pr                                           
        c00=dr*dr-dzz*dzz
	  if(c00.le.0.) c00=0.                                       
        bb=sqrt(c00)
	  srt=sqrt(s2)

	  if(idms(i2) .EQ. 1) then
	  crd =XNPION(srt,1,em2,XMAXD)     !D(1232)
	  crn1=XNPION(srt,2,em2,XMAXN1)    !N*(1440)
	  crn2=XNPION(srt,3,em2,XMAXN2)    !N*(1535)
	  elseif(idms(i2) .EQ. 2) then
	  crd= 0.
	  crn1=0.
	  crn2=XNPION(srt,3,em2,XMAXN2)
	  else
	  crd = 0.
	  crn1= 0.
	  crn2= 0.
	  endif
	  xdn = crd+crn1+crn2
	  xdn1=crd+crn1
	  if(xdn .LE. 0.) goto 22
	  dis=sqrt(xdn/xpi)
        if(bb.gt.dis)   goto 22                                                        
        e1r=sqrt(em1*em1+pr*pr)                                           
        e2r=sqrt(em2*em2+pr*pr)                                           
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.
        if(dd.lt.abs(dzz)) goto 22

	  if(idms(i2) .EQ. 1) then
        xx=RANFL(mseed)
	  if(xx.lt.crd/xdn) then
	  idr(i1)=1
	  iso(i1)=iso(i1)+icp(i2)
	  isy(i1)=0

	  elseif((xx.ge.crd/xdn) .and. (xx.lt.xdn1/xdn)) then
	  nc=iso(i1)+icp(i2)
	  if((nc.le.-1).or.(nc.ge.2)) goto 22
	  idr(i1)=2
	  iso(i1)=nc
	  isy(i1)=0

	  else
	  nc=iso(i1)+icp(i2)
	  if((nc.le.-1).or.(nc.ge.2)) goto 22
        idr(i1)= 3
	  iso(i1)= nc
	  isy(i1)= 0

	  endif
	  endif

	  if(idms(i2) .EQ. 2) then
	  nc=iso(i1)+icp(i2)
	  if((nc.le.-1).or.(nc.ge.2)) goto 22
        idr(i1)= 3
	  iso(i1)= nc
	  isy(i1)= 0
	  endif

        call LORENTZ(-bx,-by,-bz,srt,0.,0.,0.,pa0,pax,pay,paz)
        px(i1) = pax
        py(i1) = pay
        pz(i1) = paz
        em(i1) = srt

	  x(i1)  = (x1+x2)/2.
	  y(i1)  = (y1+y2)/2.
	  z(i1)  = (z1+z2)/2.

	  idms(i2)=0
	  icp(i2)= 0
	  xms(i2)= 0.
	  yms(i2)= 0.
	  zms(i2)= 0.
	  pxms(i2)=0.
	  pyms(i2)=0.
	  pzms(i2)=0.
	  emson(i2)=0.

	  goto  21
 22	  continue
21	 continue

	  return
	  end

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      	SUBROUTINE multcoll(x1,y1,z1,px1,py1,pz1,em1,
     &                        x2,y2,z2,px2,py2,pz2,em2,i1,i2,innp)
c*****    multi-collisions, e.g., NNpi-->NN
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	  parameter(pi=3.141592654, hb=0.197327, nsmax=400, delsigs=0.01, 
     &            sigs0=2.014)
	  common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/bb/np,nzp,ntag,x00
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
	  common/densi1/denstl(500),dens1(500),velcm(3)
	  common/poli/iflagcro,iflagpau,ipoli
  	  common/matx/xmatr(4,0:nsmax)
	common/SURF/ VANISH(500)
	
	  innp=0
        rpinn=1.5
	  ssq0 =sigs0*sigs0

	  if(idr(i1).GT.0 .or. idr(i2).GT.0)  return
	  e1=sqrt(em1**2+px1**2+py1**2+pz1**2)
	  is1=iso(i1)
	  iy1=isy(i1)
	  rho1=denstl(i1)
	  if(rho1.LT.1E-2)  return

	  e2=sqrt(em2**2+px2**2+py2**2+pz2**2)
	  is2=iso(i2)
	  iy2=isy(i2)
	  rho2=denstl(i2)
	  if(rho2.LT.1E-2)  return

	  rho=(rho1+rho2)/2.
        rdec=rpinn*(rho0/rho)**(1./3.)
        if(rdec.gt.2.)  rdec=2.
        dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
	  if(dist.GT.rdec) return
	  rh12=rho1*rho2
	  if(is1.EQ.is2) rh12=rh12/2.

	     do 88 ip=1,np
	     if(idms(ip).NE.1)  cycle
           ict=is1+is2+icp(ip)
		 if(ict.lt.0 .OR. ict.gt.2)  cycle
	  xpn=xms(ip)
	  ypn=yms(ip)
	  zpn=zms(ip)
	  pxp=pxms(ip)
        pyp=pyms(ip)
	  pzp=pzms(ip)
	  emp=emson(ip)
	  ipn=idms(ip)
	  ich=icp(ip)
	  if(emp .lt. 0.01)  cycle

        dis1=sqrt((x1-xpn)**2+(y1-ypn)**2+(z1-zpn)**2)
	  if(dis1.GT.rdec)  cycle
        dis2=sqrt((xpn-x2)**2+(ypn-y2)**2+(zpn-z2)**2)
	  if(dis2.GT.rdec)  cycle

           ems=sqrt(emp**2+pxp**2+pyp**2+pzp**2)
	     pxsum=px1+px2+pxp
	     pysum=py1+py2+pyp
	     pzsum=pz1+pz2+pzp
	     emsum=e1+e2+ems
		 bx=pxsum/emsum
		 by=pysum/emsum
		 bz=pzsum/emsum
	     s2=emsum**2-pxsum**2-pysum**2-pzsum**2
	     if(s2 .LE. ssq0)  cycle
           srt3=sqrt(s2)
		 pcmf=sqrt(srt3**2/4.-0.938**2)

           index=min(max(nint((srt3-sigs0)/delsigs),0),nsmax)
	     if(is1.EQ.is2) then
	     if(icp(ip).eq.0) imat=1
           if(abs(icp(ip)).eq.1) imat=2
	     else
           if(icp(ip).eq.0) imat=3
           if(abs(icp(ip)).eq.1) imat=4
		 endif
           gtot=xmatr(imat,index)*pcmf/(4.*pi*srt3)*rh12/
     &          (8.*ems*e1*e2)*hb**4
cc           gtot=gtot*2.57672
           wd=exp(-dt*gtot/hb)

            xrn=RANFL(mseed)
            if(xrn.lt.wd) goto 88
              if(ict.eq.1) then
              xrn=RANFL(nseed)
                 if(xrn.lt.0.5) then
                 iso(i1)=1
                 iso(i2)=0
		       isy(i1)=-1
		       isy(i2)=1
                 else
                 iso(i1)=0
                 iso(i2)=1
		       isy(i1)=1
		       isy(i2)=-1
                 endif
              else
              iso(i1)=nint(ict/2.+0.05)
              iso(i2)=nint(ict/2.+0.05)
			if(iso(i1).eq.1) isy(i1)=-1
			if(iso(i1).eq.0) isy(i1)=1
			if(iso(i2).eq.1) isy(i2)=-1
			if(iso(i2).eq.0) isy(i2)=1
              endif

                     cost=2.*(xrn-0.5)
                     sint=sqrt(max(1.-cost**2,0.))
                     phi=2.*pi*xrn
                     prx=pcmf*sint*cos(phi)
                     pry=pcmf*sint*sin(phi)
                     prz=pcmf*cost
		           emf=sqrt(em2**2+pcmf**2)

	  call LORENTZ(-bx,-by,-bz,emf,prx,pry,prz,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,emf,-prx,-pry,-prz,pb0,pbx,pby,pbz)

        px(i1)= pax
        py(i1)= pay
        pz(i1)= paz
        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz

	  xms(ip)= 0.
	  yms(ip)= 0.
	  zms(ip)= 0.
	  pxms(ip)=0.
	  pyms(ip)=0.
	  pzms(ip)=0.
	  idms(ip)=0.
	  icp(ip)= 0.

c	check pauli blocking
        if(iflagpau.le.0) call paulin(i1)
        if(iflagpau.eq.1) call paur34(i1,i2)
        if(iflagpau.eq.2) call pauover(i1)
	  if(iflagpau.eq.3) call paurho(i1)
	  if(ntag.eq.-1) goto 111
        if(ntag.eq.0)  then
	  innp=1
	  goto 88
	  endif

 111 	  pz(i1)=pz1                                                         
        px(i1)=px1                                                         
        py(i1)=py1   
	  iso(i1)=is1
	  isy(i1)=iy1	                                                        

        pz(i2)=pz2                                                         
        px(i2)=px2                                                         
        py(i2)=py2 
	  iso(i2)=is2
	  isy(i2)=iy2

	  xms(ip)=xpn
	  yms(ip)=ypn
	  zms(ip)=zpn
	  pxms(ip)=pxp
        pyms(ip)=pyp
	  pzms(ip)=pzp
	  emson(ip)=emp
	  idms(ip)=ipn
	  icp(ip)=ich

 88	 continue

        return                                                                
        end
c///////////////////////////////////////////////////////////////////////////////
	        SUBROUTINE rcollv  
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	  common/CLWF3/RTIME0                                             
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
        common/deu2/px1,py1,pz1,px2,py2,pz2                                     
        common/veloc/velx(500),vely(500),velz(500)                              
	  common/poli/iflagcro,iflagpau,ipoli
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/densi1/denstl(500),dens1(500),velcm(3)
	  common/medium/mdel,mdin,nsym,npop,nkaon
        ncoll=0                                                                 
        nblock=0                                                                
	ncolln=0
	ncollp=0
	nblkn=0
	nblkp=0
      xm=1.7
	xpi=31.415926
	pi=3.141592654
	hb=0.197327
	xeff=10./3.*gtau/(1.5*pi**2*rho0*hb**3)**1.66667
c************* medium	effect in nucleon-nucleon collisions *************
        call densitn

        do 21 i1=2,np                                                           
        x1=x(i1)                                                                
        y1=y(i1)                                                                
        z1=z(i1)                                                                
        ip=i1-1                                                                 
        do 22 i2=1,ip                                                           
       if(il(i1).eq.i2.and.il(i2).eq.i1)goto 22                                
       if(iflag1.eq.1.and.id(i2).eq.1.and.id(i1).eq.1)goto 22                   
       if(iflag1.eq.1.and.id(i2).eq.51.and.id(i1).eq.51)goto 22                 

				if(iflag1.eq.3)	then
		ICtrig=0
		do IIII=1,NP-NTA
			if(ID(IIII).NE.1) then
			ICtrig=1
				if(RTIME.LT.RTIME0) then
				RTIME0=RTIME
				write(*,*) rtime0
				endif
		goto 995
			endif
		enddo
995		continue
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.1.and.id(i1).eq.1)  goto 22                   
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.51.and.id(i1).eq.51)goto 22
				endif

c*****	elastic NN scattering
	if(mdel.eq.0) then
	xmdel=1.
	else
	xmdel=1.-0.2*(dens1(i1)+dens1(i2))/rho0
	endif
c*****	inelastic NN scattering
	if(mdin.eq.0) then
	xmdin=1.
	elseif(mdin.eq.1) then
	  pf1=hb*(3.*pi**2*dens1(i1))**0.3333
	  pf2=hb*(3.*pi**2*dens1(i2))**0.3333
	  reff1=1./(1.+xeff*em(i1)*pf1**3+
     &  4.*cdelt*cmeps*em(i1)*alog(cmeps*pf1**2+1.)/(cmeps*pf1**2+1.))
	  reff2=1./(1.+xeff*em(i2)*pf2**3+
     &  4.*cdelt*cmeps*em(i2)*alog(cmeps*pf2**2+1.)/(cmeps*pf2**2+1.))
	xmdin=(reff1*reff2*(em(i1)+em(i2))/(reff1*em(i1)+reff2*em(i2)))**2
	elseif(mdin.eq.2) then
	xeff=amax1(0.,1.-2.*((dens1(i1)+dens1(i2))/rho0-1.))
	xmdin=amin1(1.,xeff)
	endif

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
c                                                                               
c                                                                               
        vx1 = velx(i1)                                                          
        vx2 = velx(i2)                                                          
        vy1 = vely(i1)                                                          
        vy2 = vely(i2)                                                          
        vz1 = velz(i1)                                                          
        vz2 = velz(i2)                                                          
c                                                                               
        e1=sqrt(em(i1)*em(i1)+pz1*pz1+px1*px1+py1*py1)                          
        e2=sqrt(em(i2)*em(i2)+pz2*pz2+px2*px2+py2*py2)                          
c                                                                               
        bx=(px1+px2)/(e1+e2)                                                    
        by=(py1+py2)/(e1+e2)                                                    
        bz=(pz1+pz2)/(e1+e2)                                                    
c                                                                               
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))                                       
        s=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2                     
        if(s.lt.3.556896) go to 22         
 	  srt=sqrt(s)

	is=1
	if(iso(i1).eq.iso(i2)) is=0
c****************get the Elab to determine the minimum approaching distance
	elab=(s-(em1+em2)**2)/2./em2  ! GeV
	if(elab.le.0.) goto 22

	sigels=xmdel*CROSSels(is,elab)
	sigtot=0.
	if(is.eq.0) then
      sigtot=sigels+1.33333*xmdin*sigmand(srt)+xmdin*sigmans(srt,
     &   iso(i1),iso(i2))+1.66667*xmdin*sigmadd(srt,iso(i1),iso(i2),
     &   sdd1,sdd2)
	else
      sigtot=sigels+0.66667*xmdin*sigmand(srt)+xmdin*sigmans(srt,
     &   iso(i1),iso(i2))+xmdin*sigmadd(srt,iso(i1),iso(i2),sdd1,sdd2)
	endif

	dint=sigtot/xpi	  
		                                       
        b2=.0001+bx*bx+by*by+bz*bz                                              
        pb=bx*px1+by*py1+bz*pz1                                                 
        gg=(gam-1.)*pb/b2-gam*e1                                                
c                                                                               
        betav1 = bx*vx1 + by*vy1 + bz*vz1                                       
        betav2 = bx*vx2 + by*vy2 + bz*vz2                                       
        anenn1 = gam * (1. - betav1)                                            
        anenn2 = gam * (1. - betav2)                                            
        gg1    = gam * (gam * betav1/(gam+1.) -1.)                              
        gg2    = gam * (gam * betav2/(gam+1.) -1.)                              
        vx1r   = (vx1 + bx*gg1)/anenn1                                          
        vx2r   = (vx2 + bx*gg2)/anenn2                                          
        vy1r   = (vy1 + by*gg1)/anenn1                                          
        vy2r   = (vy2 + by*gg2)/anenn2                                          
        vz1r   = (vz1 + bz*gg1)/anenn1                                          
        vz2r   = (vz2 + bz*gg2)/anenn2                                          
c                                                                               
        pzr=gg*bz+pz1                                                           
        pxr=gg*bx+px1                                                           
        pyr=gg*by+py1                                                           
c                                                                               
        rb=bx*(x1-x2)+by*(y1-y2)+bz*(z1-z2)                                     
        gg=(gam-1.)*rb/b2                                                       
        dx=gg*bx+x1-x2                                                          
        dy=gg*by+y1-y2                                                          
        dz=gg*bz+z1-z2                                                          
c                                                                               
c                                                                               
c                                                                               
        dr=sqrt(dx*dx+dy*dy+dz*dz)                                              
        pr=sqrt(pxr*pxr+pyr*pyr+pzr*pzr)                                        
c                                                                               
        dzz=(pxr*dx+pyr*dy+pzr*dz)/pr                                           
c                                                                               
        v1parr = (vx1r*pxr + vy1r*pyr + vz1r*pzr)/pr                            
        v2parr = (vx2r*pxr + vy2r*pyr + vz2r*pzr)/pr                            
c                                                                               
        bb=sqrt(abs(dr*dr-dzz*dzz))                                                  
        if(bb.gt.1.32) goto 22                                                 
        n2=n2+1                                                                 
c                                                                               
        dd = (v1parr + v2parr) * gam *dt/2.                                                                                                                   
        if(abs(dd).lt.abs(dzz)) goto 22                                        
c                                                                               
        ncoll=ncoll+1                                                           
	if(iso(i1).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif
	if(iso(i2).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif

	ikn=0
        call prodrsn(pxr,pyr,pzr,ppx,ppy,ppz,srt,i1,i2,
     &               sigels,sigref,idp,icm,ikn,xmdin)             !! resonances (delta, N*), strangeness production and S-state pion by NN collisions
	 if(idp.NE.0) then
	 pix=pxms(icm)
	 piy=pyms(icm)
	 piz=pzms(icm)
	 endif
        call kaon(s,em(i1),em(i2),i1,i2)


       if(idp.EQ.0) then
        eaa=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        ebb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)
        call LORENTZ(-bx,-by,-bz,eaa,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,ebb,-pxr,-pyr,-pzr,pb0,pbx,pby,pbz)
        px(i1)= pax 
        py(i1)= pay 
        pz(i1)= paz 
        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz
	 else
        e11=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)
        e22=sqrt(em(i2)*em(i2)+ppx*ppx+ppy*ppy+ppz*ppz)
	  ep0=sqrt(em(i2)*em(i2)+pix*pix+piy*piy+piz*piz)
        call LORENTZ(-bx,-by,-bz,e11,pxr,pyr,pzr,pa0,pax,pay,paz)
        call LORENTZ(-bx,-by,-bz,e22,ppx,ppy,ppz,pb0,pbx,pby,pbz)
        call LORENTZ(-bx,-by,-bz,ep0,pix,piy,piz,ep1,pnx,pny,pnz)
        px(i1)= pax 
        py(i1)= pay 
        pz(i1)= paz 

        px(i2)= pbx
        py(i2)= pby
        pz(i2)= pbz

	  pxms(icm)=pnx
	  pyms(icm)=pny
	  pzms(icm)=pnz
	  goto 22
	 endif
                                    
        if(iflagpau.le.0) call paulin(i1)       
        if(iflagpau.eq.1) call paur34(i1,i2) 
        if(iflagpau.eq.2) call pauover(i1) 
 	  if(iflagpau.eq.3) call paurho(i1)
c                                                  
        if(ntag.eq.-1)go to 32                                                
c
        if(iflagpau.le.0) call paulin(i2)     
        if(iflagpau.eq.2) call pauover(i2)
	  if(iflagpau.eq.3) call paurho(i2)
clw                                                    
        if(ntag.eq.0)go to 222                                                
32      nblock=nblock+1 
C
	if(iso(i1).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
	if(iso(i2).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
C
        pz(i1)=pz1                                                            
        px(i1)=px1                                                            
        py(i1)=py1                                                            
        pz(i2)=pz2                                                            
        px(i2)=px2                                                            
        py(i2)=py2                                                            
        em(i1)=em1                                                            
        em(i2)=em2                                                            
        goto22                                                                
222     continue                                                              
        if(iflag6.eq.1)call dtest(i1,i2,em1,em2)                              
        il(i1)=i2                                                             
        il(i2)=i1                                                             
        id(i1)=id(i1)+1                                                       
        id(i2)=id(i2)+1                                                       
22      continue
21      continue

        return
        end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE rcollnew
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp                
     +   ,rrt,px(500),py(500),pz(500),em(500),iso(500),fr(500)                
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     +	ncolln,ncollp,nblkn,nblkp
      common/bb/np,nzp,ntag,x00  
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	common/CLWF3/RTIME0                                             
	common/ee/mseed,nseed                                                
      common/dd/dt,ntstep,nru,lswitch,itest,rmax,radj,initi,llt
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
      common/deu2/px1,py1,pz1,px2,py2,pz2                                     
	common/eane/neles,nneles
	common/eorne/jeles,jneles
	common/boite/ann(800),anp(800)
	common/poli/iflagcro,iflagpau,ipoli
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	neles=0
	nneles=0
      ncoll=0                                                                 
      nblock=0                                                                
	ncolln=0
	ncollp=0
	nblkn=0
	nblkp=0
      xm=1.7
      
	  do 21 i1=2,np  	                                                     
        x1=x(i1)                                                                
        y1=y(i1)                                                                
        z1=z(i1)                                                                
        ip=i1-1        
	                                                           
        do 22 i2=1,ip                                                           
        if(il(i1).eq.i2.and.il(i2).eq.i1) goto 22                                
        if(iflag1.eq.1.and.id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.1.and.id(i2).eq.51.and.id(i1).eq.51) goto 22                 

		  if(iflag1.eq.3)	then
		ICtrig=0
		do IIII=1,NP-NTA
			if(id(IIII).NE.1) then
			ICtrig=1 
				if(RTIME.LT.RTIME0) then
				RTIME0=RTIME
				write(*,*) rtime0
				endif
		goto 995
		    endif
		enddo
995		continue
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.51.and.id(i1).eq.51) goto 22                   
		  endif

	is=1
	if(iso(i1).eq.iso(i2)) is=0

        x2=x(i2)                                                                
        y2=y(i2)                                                                
        z2=z(i2)                                                                
        r=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)                       
        if (r.gt.3.4*3.4) goto 22
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
        bx=(px1+px2)/(e1+e2)                                                    
        by=(py1+py2)/(e1+e2)                                                    
        bz=(pz1+pz2)/(e1+e2)                                                    
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))
        s=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2
        if(s.lt.3.556896) goto 22

c****************get the Elab to determine the minimum approaching distance
	elab=(s-(em1+em2)**2)/2./em2
	if(elab.le.0.) goto 22
	it=elab*1000.
	if(it.ge.800) it=800
	if(it.le.1) it=1
	if(is.eq.0) dint=ann(it)
	if(is.eq.1) dint=anp(it)

        b2=.00001+bx*bx+by*by+bz*bz
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
        if(bb.gt.sqrt(dint)) goto 22
        e1r=sqrt(em(i1)*em(i1)+pr*pr)                                           
        e2r=sqrt(em(i2)*em(i2)+pr*pr)                                           
        dd=pr*(1./e1r+1./e2r)*gam*dt/2.                                         
        if(abs(dd).lt.abs(dzz)) go to 22                                        
        ncoll=ncoll+1                                                           
C
	if(iso(i1).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif
	if(iso(i2).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif
C
        srt=sqrt(s)                  
        s=srt                                   
        as=(3.65*(s-1.866))**6                                                  
        a=6.*as/(1.+as)                                                         
        ta=-2*pr*pr
	  xx=RANFL(mseed)                                         
        t1=alog((1-xx)*XEXP(a*ta)+xx)/a                       
        c1=1.-t1/ta                                                             
        t1=6.28*RANFL(mseed) 
        pr=sqrt(pxr*pxr+pyr*pyr+pzr*pzr)
        c2=pzr/pr

        s2=sqrt(1.-c2*c2)                                                       
      if(pxr.eq.0.0.and.pyr.eq.0.0) goto 6 
        t2=atan2(pyr,pxr) 
   6    s1=sqrt(1.-c1*c1)                                                       
        ct1=cos(t1)                                                             
        st1=sin(t1)                                                             
        ct2=cos(t2)                                                             
        st2=sin(t2)                                                             
        pzr=pr*(c1*c2-s1*s2*ct1)
        ss=c2*s1*ct1+s2*c1                                                      
        pxr=pr*(ss*ct2-s1*st1*st2)
        pyr=pr*(ss*st2+s1*st1*ct2)

        pb=bx*pxr+by*pyr+bz*pzr                                                 
        ea=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*pb/b2+gam*ea                                                
        px(i1)=gg*bx+pxr                                                        
        py(i1)=gg*by+pyr                                                        
        pz(i1)=gg*bz+pzr                                                       
        eb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*(-pb)/b2+gam*eb                                         
        px(i2)=gg*bx-pxr                                                        
        py(i2)=gg*by-pyr                                                        
        pz(i2)=gg*bz-pzr
	 
c        call pauli12(i1,i2)
c	if(ntag.eq.-1)  goto 32
c	if(ntag.eq.0)   goto 222
	
	  if(iflagpau.le.0) call paulin(i1)
        if(iflagpau.eq.1) call paur34(i1,i2)
        if(iflagpau.eq.2) call pauover(i1)
        if(iflagpau.eq.3) call paurho(i1)
        if(ntag.eq.-1)  goto 32                                              
c
        if(iflagpau.le.0) call paulin(i2)
        if(iflagpau.eq.2) call pauover(i2)
	  if(iflagpau.eq.3) call paurho(i2)
        if(ntag.eq.0)   goto 222
c                                           
32 	nblock=nblock+1
C
	if(iso(i1).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
	if(iso(i2).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
C
        pz(i1)=pz1                                                         
        px(i1)=px1                                                         
        py(i1)=py1                                                         
        pz(i2)=pz2                                                         
        px(i2)=px2                                                         
        py(i2)=py2                                                         
        em(i1)=em1                                                         
        em(i2)=em2                                                         
        goto 22                                                             
222     continue  
cc feng                                                         
        if(iflag6.eq.1)call dtest(i1,i2,em1,em2)                           
        il(i1)=i2                                                          
        il(i2)=i1                                                          
        id(i1)=id(i1)+1                                                    
	  id(i2)=id(i2)+1                                                    
22   	continue                                                           
21    	continue                                                           

      	return                                                             
      	end                                                                
c///////////////////////////////////////////////////////////////////////////////

        	SUBROUTINE rcollvnew

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp                
     +   ,rrt,px(500),py(500),pz(500),em(500),iso(500),fr(500)                 
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     +	ncolln,ncollp,nblkn,nblkp
      common/bb/np,nzp,ntag,x00  
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	common/CLWF3/RTIME0
	common/ee/mseed,nseed                                                
      common/dd/dt,ntstep,nru,lswitch,itest,rmax,radj,initi,llt
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
      common/deu2/px1,py1,pz1,px2,py2,pz2                                     
	common/eane/neles,nneles
	common/eorne/jeles,jneles
	common/boite/ann(800),anp(800)
	common/poli/iflagcro,iflagpau,ipoli
	common/veloc/velx(500),vely(500),velz(500)
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	neles=0
	nneles=0
      ncoll=0                                                                 
      nblock=0                                                                
	ncolln=0
	ncollp=0
	nblkn=0
	nblkp=0
        xm=1.7                                                                  
        do 21 i1=2,np                                                           
        x1=x(i1)                                                                
        y1=y(i1)                                                                
        z1=z(i1)                                                                
        ip=i1-1                                                                 
        do 22 i2=1,ip                                                           
        if(il(i1).eq.i2.and.il(i2).eq.i1) goto 22                                
        if(iflag1.eq.1.and.id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.1.and.id(i2).eq.51.and.id(i1).eq.51) goto 22                 

				if(iflag1.eq.3)	then
		ICtrig=0
		do IIII=1,NP-NTA
			if(ID(IIII).NE.1) then
			ICtrig=1
				if(RTIME.LT.RTIME0) then
				RTIME0=RTIME
				write(*,*) rtime0
				endif
		goto 995
			endif
		enddo
995		continue
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.51.and.id(i1).eq.51) goto 22                   
				endif
	is=1
	if(iso(i1).eq.iso(i2)) is=0

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
	vx1=velx(i1)	
	vx2=velx(i2)	
	vy1=vely(i1)	
	vy2=vely(i2)	
	vz1=velz(i1)	
	vz2=velz(i2)	
cc
        e1=sqrt(em(i1)*em(i1)+pz1*pz1+px1*px1+py1*py1)                          
        e2=sqrt(em(i2)*em(i2)+pz2*pz2+px2*px2+py2*py2)                          
        bx=(px1+px2)/(e1+e2)                                                    
        by=(py1+py2)/(e1+e2)                                                    
        bz=(pz1+pz2)/(e1+e2)                                                    
        gam=1./sqrt(ABS(1.-bx*bx-by*by-bz*bz))                                       
        s=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2                     
        if(s.lt.3.556896) goto 22
       
	elab=(s-(em1+em2)**2)/2./em2
	if(elab.le.0.) goto 22
	it=elab*1000.
	if(it.ge.800) it=800
	if(it.le.1) it=1
	if(is.eq.0) dint=ann(it)
	if(is.eq.1) dint=anp(it)
        b2=.0001+bx*bx+by*by+bz*bz                                              
        pb=bx*px1+by*py1+bz*pz1                                                 
        gg=(gam-1.)*pb/b2-gam*e1                                                
c
	betav1=bx*vx1+by*vy1+bz*vz1	
	betav2=bx*vx2+by*vy2+bz*vz2
	anenn1=gam*(1.-betav1)	
	anenn2=gam*(1.-betav2)
	gg1=gam*(gam*betav1/(gam+1.)-1.)
	gg2=gam*(gam*betav2/(gam+1.)-1.)
	vx1r=(vx1+bx*gg1)/anenn1
	vx2r=(vx2+bx*gg2)/anenn2
	vy1r=(vy1+by*gg1)/anenn1
	vy2r=(vy2+by*gg2)/anenn2
	vz1r=(vz1+bz*gg1)/anenn1
	vz2r=(vz2+bz*gg2)/anenn2
c
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
c
	v1parr=(vx1r*pxr+vy1r*pyr+vz1r*pzr)/pr
	v2parr=(vx2r*pxr+vy2r*pyr+vz2r*pzr)/pr
c
        bb=sqrt(abs(dr*dr-dzz*dzz))
        if(bb.gt.sqrt(dint)) go to 22
        n2=n2+1
	dd=(v1parr+v2parr)*gam*dt/2.
        if(abs(dd).lt.abs(dzz)) go to 22                                        
        ncoll=ncoll+1                                                           
C
	if(iso(i1).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif
	if(iso(i2).eq.0) then
	ncolln=ncolln+1
	else
	ncollp=ncollp+1
	endif
C
        srt=sqrt(abs(s))                  
        s=srt                                   
        as=(3.65*(s-1.866))**6                                                  
        a=6.*as/(1.+as)                                                         
        ta=-2*pr*pr
	  xx=RANFL(mseed)                                         
        t1=alog((1-xx)*XEXP(a*ta)+xx)/a                       
        c1=1.-t1/ta                                                             
        t1=6.28*RANFL(mseed) 
        pr=sqrt(pxr*pxr+pyr*pyr+pzr*pzr)
        c2=pzr/pr

        s2=sqrt(1.-c2*c2)                                                       
      if(pxr.eq.0.0.and.pyr.eq.0.0) go to 6 
        t2=atan2(pyr,pxr) 
   6    s1=sqrt(1.-c1*c1)                                                       
        ct1=cos(t1)                                                             
        st1=sin(t1)                                                             
        ct2=cos(t2)                                                             
        st2=sin(t2)                                                             
        pzr=pr*(c1*c2-s1*s2*ct1)
        ss=c2*s1*ct1+s2*c1                                                      
        pxr=pr*(ss*ct2-s1*st1*st2)
        pyr=pr*(ss*st2+s1*st1*ct2)

        pb=bx*pxr+by*pyr+bz*pzr                                                 
        ea=sqrt(em(i1)*em(i1)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*pb/b2+gam*ea                                                
        px(i1)=gg*bx+pxr                                                        
        py(i1)=gg*by+pyr                                                        
        pz(i1)=gg*bz+pzr                                                        
        eb=sqrt(em(i2)*em(i2)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*(-pb)/b2+gam*eb                                             
        px(i2)=gg*bx-pxr                                                        
        py(i2)=gg*by-pyr                                                        
        pz(i2)=gg*bz-pzr   
clw
        if(iflagpau.le.0) call paulin(i1)
        if(iflagpau.eq.1) call paur34(i1,i2)
        if(iflagpau.eq.2) call pauover(i1)
	  if(iflagpau.eq.3) call paurho(i1)
c
         if(ntag.eq.-1)go to 32                                             
c
        if(iflagpau.le.0) call paulin(i2)
        if(iflagpau.eq.2) call pauover(i2)
	  if(iflagpau.eq.3) call paurho(i2)
clw
        if(ntag.eq.0)go to 222                                             
32 	nblock=nblock+1
C
	if(iso(i1).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
	if(iso(i2).eq.0) then
	nblkn=nblkn+1
	else
	nblkp=nblkp+1
	endif
C
        pz(i1)=pz1                                                         
        px(i1)=px1                                                         
        py(i1)=py1                                                         
        pz(i2)=pz2                                                         
        px(i2)=px2                                                         
        py(i2)=py2                                                         
        em(i1)=em1                                                         
        em(i2)=em2                                                         
        goto 22                                                             
222     continue                                                           
        if(iflag6.eq.1)call dtest(i1,i2,em1,em2)                           
        il(i1)=i2                                                          
        il(i2)=i1                                                          
        id(i1)=id(i1)+1
	  id(i2)=id(i2)+1                                                    
22   	continue                                                           
21    	continue                                                           
      	return                                                             
      	end
c///////////////////////////////////////////////////////////////////////////////
      SUBROUTINE nrcoll(crsc)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,               
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &  ncolln,ncollp,nblkn,nblkp
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	  common/CLWF3/RTIME0                                            
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
        common/deu2/px1,py1,pz1,px2,py2,pz2                                     
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	  ncoll=0                                                                 
        nblock=0                                                                
        do 21 i1=2,np                                                           
        pz1=pz(i1)                                                              
        px1=px(i1)                                                              
        py1=py(i1)                                                              
        em1=em(i1)                                                              
        ip=i1-1                                                                 
        do 22 i2=1,ip                                                           
        if(il(i1).eq.i2.and.il(i2).eq.i1) goto 22                                
        if(iflag1.eq.1.and.id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.1.and.id(i2).eq.51.and.id(i1).eq.51) goto 22                 

			if(iflag1.eq.3)	then
		ICtrig=0
		do IIII=1,NP-NTA
			if(ID(IIII).NE.1) then
			ICtrig=1
				if(RTIME.LT.RTIME0) then
				RTIME0=RTIME
				write(*,*) rtime0
				endif
		goto 995
			endif
		enddo
995		continue
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.1.and.id(i1).eq.1)   goto 22                   
        if(iflag1.eq.3.AND.ICtrig.EQ.0.and.
     &	  id(i2).eq.51.and.id(i1).eq.51) goto 22                   
				endif

        pz2=pz(i2)                                                              
        px2=px(i2)                                                              
        py2=py(i2)                                                              
        em2=em(i2)                                                              
        cror1=1.142                                                             
        cror=cror1+crsc
        xr=x(i1)-x(i2)                                                          
        yr=y(i1)-y(i2)                                                          
        zr=z(i1)-z(i2)                                                          
        r2=xr*xr+yr*yr+zr*zr                                                    
        if(r2.gt.cror*cror) goto 22                                           
        prx=(px(i1)-px(i2))/2.                                                  
        pry=(py(i1)-py(i2))/2.                                                  
        prz=(pz(i1)-pz(i2))/2.                                                  
        pr=sqrt(prx*prx+pry*pry+prz*prz)                                        
        vr=2*pr/0.938                                                           
        dzz=(prx*xr+pry*yr+prz*zr)/pr                                           
        bb2=r2-dzz*dzz                                                          
        if(bb2.gt.cror1*cror1)  goto 22                                         
        dd=vr*dt/2.                                                             
        if(abs(dd).lt.abs(dzz)) goto 22                                        
        pcz=(pz(i1)+pz(i2))/2.                                                  
        pcx=(px(i1)+px(i2))/2.                                                  
        pcy=(py(i1)+py(i2))/2.                                                  
        ncoll=ncoll+1 
        cth=1.-2.*RANFL(mseed)                                       
        sth=sqrt(1.-cth**2.)
        phi=6.28*RANFL(mseed)                                                    
        przn=pr*cth                                                             
        prxn=pr*sth*cos(phi)                                                    
        pryn=pr*sth*sin(phi)                                                    
        px(i1)=prxn+pcx                                                         
        py(i1)=pryn+pcy                                                         
        pz(i1)=przn+pcz                                                         
        px(i2)=pcx-prxn                                                         
        py(i2)=pcy-pryn                                                         
        pz(i2)=pcz-przn
                                                                                
        call paulin(i1)                   
        if(ntag.eq.-1) goto 32            
        call paulin(i2)                   
        if(ntag.eq.0)  goto 551            
32      nblock=nblock+1                   
        pz(i1)=pcz+prz                    
        px(i1)=pcx+prx                    
        py(i1)=pcy+pry                    
        pz(i2)=pcz-prz                    
        px(i2)=pcx-prx                    
        py(i2)=pcy-pry                    
        goto 22                            
551     continue                          
        il(i1)=i2                                      
        il(i2)=i1                         
        if(iflag6.eq.1) call dtest(i1,i2,em1,em2)       
        id(i1)=id(i1)+1                                
        id(i2)=id(i2)+1                                
22      continue                                       
21      continue                                       
        return                                         
	end        
c///////////////////////////////////////////////////////////////////////////////
	      SUBROUTINE paulin(i)                                            
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/bb/np,nzp,ntag,x00
	  common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/ee/mseed,nseed
        common/paulis/radn,padn,anorm,const,radn2,padn2                       
	  common/poli/iflagcro,iflagpau,ipoli
	  common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     &	frtwo(500)
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/xx/alpha,beta,gamma,symm
	  common/wavel/sl(500)
        dimension con(500)  
	common/SURF/ VANISH(500)
	pi=3.1415926    
	epsim=cmeps*(0.75*rho0*pi**2)**0.66667*0.197327**2	                                                
c
      if(em(i).gt.1.) return                                                 
	if(ipoli.eq.0) then
	ntag=-1
	return
	endif
	if(ipoli.eq.1) then
	ntag=0
	return
	endif
c
        ntag=0                                                                
        sum=0.                                                                
        psum=0.                                                               
        contrs=0.                                                             
        xcm=0.                                                                
        ycm=0.                                                                
        zcm=0.                                                                
        pxcm=0.                                                               
        pycm=0.                                                               
        pzcm=0.                                                               
        nn=0                                                                  
        nn1=0                                                                 
        nn2=0                                                                 
        x1=x(i)                                                               
        px1=px(i)                                                             
        y1=y(i)                                                               
        py1=py(i)                                                             
        z1=z(i)                                                               
        pz1=pz(i)                                                             
        pt=sqrt(px1**2+py1**2+pz1**2)                                           
        rt=sqrt(x1**2+y1**2+z1**2)                                              
        do 1 j=1,np
		if(VANISH(J).EQ.1) cycle                                                           
        con(j)=0.                                                             
        if(j.eq.i) goto 1                                                     
        if(em(j).gt.1.) goto 1                                                  
        dx=x(j)-x1                                                            
        dy=y(j)-y1                                                            
        dz=z(j)-z1                                                            
        dpx=px(j)-px1                                                         
        dpy=py(j)-py1                                                         
        dpz=pz(j)-pz1                                                         
        pr2=dpx*dpx+dpy*dpy+dpz*dpz                                           
        r2=dx*dx+dy*dy+dz*dz                                                  
        if(pr2.gt.padn2) goto 1                                               
        if(r2.gt.radn2)  goto 1                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)  
	al=2.*(sl(i)+sl(j)) 
	sal=(pi*al)**1.5*rho0
      densty = (fr(i)+1.)/sal                                               
      pot=alpha*densty+beta*densty**gamma+gtau*(eta+1.)*densty**eta+
     &    cdelt*densty*alog(epsim*densty**0.66667+1.)**2 
	if(ifdelt.ge.1) then
	dentwo=(frtwo(i)+1.)/salnew
	pot=-.353*dentwo+.303*densty**(7./6.)
	endif
        etot = pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
        ahr=radn-r/2.                                                         
        ahp=padn-pr/2.                                                        
        contr=const*ahr*ahr*(3*radn-ahr)*ahp*ahp*(3*padn-ahp)                 
        con(j)=contr                                                          
        xcm=xcm+x(j)                                                          
        ycm=ycm+y(j)                                                          
        zcm=zcm+z(j)                                                          
        pxcm=pxcm+px(j)                                                       
        pycm=pycm+py(j)                                                       
        pzcm=pzcm+pz(j)                                                       
        if(contr.lt.0.00001) goto 1                                             
        nn=nn+1                                                               
        contrs=contrs+contr                                                   
1       continue                                                              
        if(nn.eq.0) return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np 
		if(VANISH(K).EQ.1)	cycle                                                          
        if(k.eq.i)goto2                                                       
        dx=x(k)-x1                                                            
        dy=y(k)-y1                                                            
        dz=z(k)-z1                                                            
        pdx=px(k)-px1                                                         
        pdy=py(k)-py1                                                         
        pdz=pz(k)-pz1                                                         
        dr=sqrt(dx**2+dy**2+dz**2)                                            
        pdr=sqrt(pdx**2+pdy**2+pdz**2)
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm                                    
        if(xcos.gt.0.707) sum=sum+con(k)                                       
        if(xcos.gt.0.707) nn1=nn1+1                                            
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm                           
        if(pxcos.gt.0.707) psum=psum+con(k)                                    
        if(pxcos.gt.0.707) nn2=nn2+1                                           
2       continue                                                              
        sum=sum*6.83                                                          
        psum=psum*6.83                                                        
        sumf=amax1(sum,psum)                                                  
        if(sumf.lt.contrs) sumf=contrs
        if(etot.gt.0.005)  sumf=contrs
        if(sumf.gt.RANFL(mseed)) ntag=-1
        return                                                                
        end                      
c///////////////////////////////////////////////////////////////////////////////
	      SUBROUTINE pauli12(i1,i2)                                            
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/bb/np,ntag,x00
      common/ee/mseed,nseed
	common/psdc/nntag,ffnn
      
                                                
      nntag=0
      f1=Fn(i1)
	f1=ffnn
      f2=Fn(i2)
	f2=ffnn
      f0=0.622
      fi1=amin1(1.,(f1-f0)/(1.-f0))
      fi2=amin1(1.,(f2-f0)/(1.-f0))
      
      fi12=1.-(1.-fi1)*(1.-fi2)
      if(fi12.ge.RAN(nseed)) nntag=-1     
      
      return
      end	
c///////////////////////////////////////////////////////////////////////////////
	      SUBROUTINE paur34(i1,i2)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/bb/np,nzp,ntag,x00
	  common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/ee/mseed,nseed
        common/paulis/radn,padn,anorm,const,radn2,padn2 
	  common/poli/iflagcro,iflagpau,ipoli
	common/SURF/ VANISH(500)

        ntag = 0    
	
	if(i2. EQ. 0)  goto 22    ! resonance decay into nucleon and meson

	if(em(i1).gt.1. .or. em(i2).gt.1.) then
	ntag=0
	return                                                 
	endif
	
	if(ipoli.eq.0) then
	ntag=-1
	return
	elseif(ipoli.eq.1) then
	ntag=0
	return
	endif
                                                                                                                        
        ctrs1= 0.                                                             
        ctrs2= 0.

        x1 =x(i1)                                                               
        px1=px(i1)                                                             
        y1 =y(i1)                                                               
        py1=py(i1)                                                             
        z1 =z(i1)                                                               
        pz1=pz(i1)

        x2 =x(i2)                                                               
        px2=px(i2)                                                             
        y2 =y(i2)                                                               
        py2=py(i2)                                                             
        z2 =z(i2)                                                               
        pz2=pz(i2)                                                             

        do 1 j=1,np 
		if(VANISH(J).EQ.1) cycle                                                          
        if(j .eq. i1)       goto 1      
	  if(iso(j) .NE. iso(i1))  goto 1
        dx = x(j)-x1                                                           
        dy = y(j)-y1 
        dz = z(j)-z1 
        dpx= px(j)-px1                                                         
        dpy= py(j)-py1                                                         
        dpz= pz(j)-pz1                                                         
        pr2=dpx*dpx+dpy*dpy+dpz*dpz 
        r2=dx*dx+dy*dy+dz*dz                                                  
        if(pr2 .gt. padn2)  goto 1                                               
        if(r2  .gt. radn2)  goto 1                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)   

        ahr=radn-r/2.                                                         
        ahp=padn-pr/2.                                                        
        contr=const*ahr*ahr*(3.*radn-ahr)*ahp*ahp*(3.*padn-ahp)                                                         
                                                       
        if(contr.lt.0.00001) goto 1                                             
        ctrs1 = ctrs1+contr                                                   
 1      continue 
	  fi1=ctrs1

        do 2 j=1,np  
		if(VANISH(J).EQ.1)	cycle                                                         
        if(j .eq. i2)       goto 2                                               
	  if(iso(j) .NE. iso(i2))  goto 2
        dx = x(j)-x2
        dy = y(j)-y2 
        dz = z(j)-z2
        dpx= px(j)-px2                                                         
        dpy= py(j)-py2                                                         
        dpz= pz(j)-pz2
        pr2=dpx*dpx+dpy*dpy+dpz*dpz 
        r2=dx*dx+dy*dy+dz*dz                                                  
        if(pr2 .gt. padn2)  goto 2                                               
        if(r2  .gt. radn2)  goto 2                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)   

        ahr=radn-r/2.                                                         
        ahp=padn-pr/2.                                                        
        contr=const*ahr*ahr*(3.*radn-ahr)*ahp*ahp*(3.*padn-ahp)                                                         
                                                       
        if(contr.lt.0.00001) goto 2 
        ctrs2 = ctrs2+contr                                                   
 2      continue 
	  fi2=ctrs2

	  if(fi1 .gt. 1.)  fi1=1.
	  if(fi2 .gt. 1.)  fi2=1.
	  xx=ran(mseed)      
        fi12=1.-(1.-fi1)*(1.-fi2)
        if(fi12 .gt. RANFL(nseed))  ntag=-1 

	  return

22  	if(em(i1) .gt. 1.) then   
	ntag=0
	return                                                 
	endif
	
	if(ipoli.eq.0) then
	ntag=-1
	return
	elseif(ipoli.eq.1) then
	ntag=0
	return
	endif
                                                                                                                          
        ctrs = 0. 

        x1 =x(i1)                                                               
        px1=px(i1)                                                             
        y1 =y(i1)                                                               
        py1=py(i1)                                                             
        z1 =z(i1)                                                               
        pz1=pz(i1)                                                  

        do 3 j=1,np 
		if(VANISH(J).EQ.1)	cycle                                                          
        if(j .eq. i1)       goto 3      
	  if(iso(j) .NE. iso(i1))  goto 3
        dx = x(j)-x1                                                           
        dy = y(j)-y1 
        dz = z(j)-z1 
        dpx= px(j)-px1 
        dpy= py(j)-py1 
        dpz= pz(j)-pz1  
        pr2=dpx*dpx+dpy*dpy+dpz*dpz 
        r2=dx*dx+dy*dy+dz*dz                                                  
        if(pr2 .gt. padn2)  goto 3                                      
        if(r2  .gt. radn2)  goto 3                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)   

        ahr=radn-r/2.                                                         
        ahp=padn-pr/2.                                                        
        contr=const*ahr*ahr*(3.*radn-ahr)*ahp*ahp*(3.*padn-ahp)                                                         
                                                       
        if(contr.lt.0.00001) goto 3                                           
        ctrs = ctrs + contr                                                   
 3      continue 
	  fi = ctrs
        if(fi .gt. RANFL(nseed))  ntag=-1 

	return
      end     
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
		SUBROUTINE pauover(i)
c/////////////////////////////////////////////////////////////////////////////
	common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	common/bb/np,nzp,ntag,x00
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
	common/ee/mseed,nseed
	common/poli/iflagcro,iflagpau,ipoli
	common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     &	frtwo(500)
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/xx/alpha,beta,gamma,symm
      common/wavel/sl(500)
	common/SURF/ VANISH(500)
	dimension con(500)
	pi=3.1415926
	epsim=cmeps*(0.75*rho0*pi**2)**0.66667*0.197327**2
c
	if(em(i).gt.1.) return
	if(ipoli.eq.0) then
	ntag=-1
	return
	endif
	if(ipoli.eq.1) then
	ntag=0
	return
	endif
c
	ntag=0
        sum=0.                                                                
        psum=0.                                                               
	contrs=0.
        xcm=0.                                                                
        ycm=0.                                                                
        zcm=0.                                                                
        pxcm=0.                                                               
        pycm=0.                                                               
        pzcm=0.                                                               
        nn=0                                                                  
        nn1=0                                                                 
        nn2=0                                                                 
	x1=x(i)
	y1=y(i)
	z1=z(i)
	px1=px(i)
	py1=py(i)
	pz1=pz(i)
c
	if(RANFL(mseed).le.0.5) then
	ispin=-1
	else
	ispin=1
	endif
c	
	do 1 j=1,np
		if(VANISH(J).EQ.1) cycle
c
	if(RANFL(mseed).le.0.5) then
	jspin=-1
	else
	jspin=1
	endif
c
	con(j)=0.                                                             
	if(j.eq.i) goto 1
	if(em(j).gt.1.) goto 1
	if(iso(j).ne.iso(i)) goto 1
	if(jspin.ne.ispin) goto 1
	dx=x(j)-x1
	dy=y(j)-y1
	dz=z(j)-z1
	dpx=px(j)-px1
	dpy=py(j)-py1
	dpz=pz(j)-pz1
	r2=dx*dx+dy*dy+dz*dz
	p2=dpx*dpx+dpy*dpy+dpz*dpz
	al=2.*(sl(i)+sl(j)) 
	sal=(pi*al)**1.5*rho0
      densty = (fr(i)+1.)/sal                                               
      pot=alpha*densty+beta*densty**gamma+gtau*(eta+1.)*densty**eta+
     &    cdelt*densty*alog(epsim*densty**0.66667+1.)**2                                
	if(ifdelt.ge.1) then
	dentwo=(frtwo(i)+1.)/salnew
	pot=-.353*dentwo+.303*densty**(7./6.)
	endif
        etot = pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
C----------------------------------------------------------------------
	contr=XEXP(-r2/sl(j)/4.-p2*sl(j)/0.197/0.197)
C----------------------------------------------------------------------
	con(j)=contr
        xcm=xcm+x(j)                                                          
        ycm=ycm+y(j)                                                          
        zcm=zcm+z(j)                                                          
        pxcm=pxcm+px(j)                                                       
        pycm=pycm+py(j)                                                       
        pzcm=pzcm+pz(j)                                                       
        if(contr.lt.0.00001) goto 1                                             
        nn=nn+1                                                               
	contrs=contrs+contr	
1	continue
c
        if(nn.eq.0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np 
		if(VANISH(K).EQ.1)	cycle                                                          
        if(k.eq.i)goto2                                                       
	if(ISO(K).NE.ISO(I)) goto 2
        dx=x(k)-x1                                                            
        dy=y(k)-y1                                                            
        dz=z(k)-z1                                                            
        pdx=px(k)-px1                                                         
        pdy=py(k)-py1                                                         
        pdz=pz(k)-pz1                                                         
        dr=sqrt(dx**2+dy**2+dz**2)
        pdr=sqrt(pdx**2+pdy**2+pdz**2)
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm
        if(xcos.gt.0.707)sum=sum+con(k)                                       
        if(xcos.gt.0.707)nn1=nn1+1                                            
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm                           
        if(pxcos.gt.0.707)psum=psum+con(k)                                    
        if(pxcos.gt.0.707)nn2=nn2+1                                           
2       continue                                                              
        sum=sum*6.83                                                          
        psum=psum*6.83                                                        
        sumf=amax1(sum,psum)                                                  
        if(sumf.lt.contrs)sumf=contrs                                         
        if(etot.gt.0.005)sumf=contrs 
        if(sumf.gt.RANFL(mseed))ntag=-1                            
        return                                                                
        end                                                                   


c///////////////////////////////////////////////////////////////////////////////
	      SUBROUTINE paurho(i)                                            
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/bb/np,nzp,ntag,x00
	  common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/ee/mseed,nseed
	  common/poli/iflagcro,iflagpau,ipoli
	  common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     &	frtwo(500)
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/xx/alpha,beta,gamma,symm
	  common/wavel/sl(500)
	common/SURF/ VANISH(500)
        dimension con(500),sden(500),rad(500),pad(500) 
	pi=3.1415926
	epsim=cmeps*(0.75*rho0*pi**2)**0.66667*0.197327**2
c
      if(em(i).gt.1.)return                                                 
	if(ipoli.eq.0) then
	ntag=-1
	return
	endif
	if(ipoli.eq.1) then
	ntag=0
	return
	endif
c
        ntag=0                                                                
        sum=0.                                                                
        psum=0.                                                               
        contrs=0.                                                             
        xcm=0.                                                                
        ycm=0.                                                                
        zcm=0.                                                                
        pxcm=0.                                                               
        pycm=0.                                                               
        pzcm=0.                                                               
        nn=0                                                                  
        nn1=0                                                                 
        nn2=0                                                                 
        x1=x(i)                                                               
        px1=px(i)                                                             
        y1=y(i)                                                               
        py1=py(i)                                                             
        z1=z(i)                                                               
        pz1=pz(i)                                                             
        pt=sqrt(px1**2+py1**2+pz1**2)                                           
        rt=sqrt(x1**2+y1**2+z1**2)                                              

	do j=1,500
	sden(j)=0.
	rad(j)=0.
	pad(j)=0.
	enddo

	ccwl=(2.*PI*sl(i))**1.5
	do j=1,np
		if(VANISH(J).EQ.1)	cycle
	   do 96916 k=1,np
		if(VANISH(K).EQ.1)	cycle
	   if(iso(k).ne.iso(j)) goto 96916
	   dis2=(x(k)-x(j))**2.+(y(k)-y(j))**2.+(z(k)-z(j))**2.
	   sden(j)=sden(j)+XEXP(-dis2/2./sl(j))/ccwl
96916	   continue
	enddo

	consta=1.24*(9./PI/PI/32.)**(1./3.)
	do j=1,np
		if(VANISH(J).EQ.1)	cycle
	pad(j)=.197*(3.*PI*PI*sden(j))**(1./3.)
	rad(j)=consta/pad(j)
	enddo

        do 1 j=1,np
		if(VANISH(J).EQ.1)	cycle                                                           
        con(j)=0.                                                             
        if(j.eq.i)go to 1                                                     
        if(em(j).gt.1.)goto1                                                  
	if(iso(j).ne.iso(i)) goto 1
        dx=x(j)-x1                                                            
        dy=y(j)-y1                                                            
        dz=z(j)-z1                                                            
        dpx=px(j)-px1                                                         
        dpy=py(j)-py1                                                         
        dpz=pz(j)-pz1                                                         
        pr2=dpx*dpx+dpy*dpy+dpz*dpz                                           
        r2=dx*dx+dy*dy+dz*dz
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)
	al=2.*(sl(i)+sl(j)) 
	sal=(pi*al)**1.5*rho0
      densty = (fr(i)+1.)/sal                                               
      pot=alpha*densty+beta*densty**gamma+gtau*(eta+1.)*densty**eta+
     &    cdelt*densty*alog(epsim*densty**0.66667+1.)**2              
	if(ifdelt.ge.1) then
	dentwo=(frtwo(i)+1.)/salnew
	pot=-.353*dentwo+.303*densty**(7./6.)
	endif
        etot = pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
	
	if(r.ge.(rad(j)+rad(i)).or.pr.ge.(pad(j)+pad(j))) goto 1
	if(r.le.abs(rad(j)-rad(i))) then
	rrmin=amin1(rad(j),rad(i))
	phaser=4./3.*PI*(rrmin)**3.
	endif
	if(r.gt.abs(rad(j)-rad(i)).and.r.lt.(rad(j)+rad(j))) then
	rh1=(rad(j)**2.-(r-rad(i))**2.)/2./r
	rh2=(rad(i)**2.-(r-rad(j))**2.)/2./r
	phaser=PI/3.*(rh1**2.*(3.*rad(i)-rh1)+rh2**2.*(3.*rad(j)-rh2))
	endif

	if(pr.le.abs(pad(j)-pad(i))) then
	ppmin=amin1(pad(j),pad(i))
	phasep=4./3.*PI*(ppmin)**3.
	endif
	if(pr.gt.abs(pad(j)-pad(i)).and.pr.lt.(pad(j)+pad(i))) then
	ph1=(pad(j)**2.-(pr-pad(i))**2.)/2./pr
	ph2=(pad(i)**2.-(pr-pad(j))**2.)/2./pr
	phasep=PI/3.*(ph1**2.*(3.*pad(i)-ph1)+ph2**2.*(3.*pad(j)-ph2))
	endif

	phasei=1.24**3./2.
	contr=phaser*phasep/phasei
        con(j)=contr                                                          
        xcm=xcm+x(j)                                                          
        ycm=ycm+y(j)                                                          
        zcm=zcm+z(j)                                                          
        pxcm=pxcm+px(j)                                                       
        pycm=pycm+py(j)                                                       
        pzcm=pzcm+pz(j)                                                       
        if(contr.lt.0.00001)goto1                                             
        nn=nn+1                                                               
        contrs=contrs+contr                                                   
1       continue                                                              
        if(nn.eq.0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np 
	  if(VANISH(K).EQ.1)	cycle
        if(k.eq.i)goto2                                                       
	if(iso(k).ne.iso(i)) goto 2
        dx=x(k)-x1 
	dy=y(k)-y1
        dz=z(k)-z1 
        pdx=px(k)-px1 
        pdy=py(k)-py1
        pdz=pz(k)-pz1 
        dr=sqrt(dx**2+dy**2+dz**2)
        pdr=sqrt(pdx**2+pdy**2+pdz**2) 
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm
        if(xcos.gt.0.707) sum=sum+con(k) 
	  if(xcos.gt.0.707) nn1=nn1+1
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm  
        if(pxcos.gt.0.707) psum=psum+con(k)  
        if(pxcos.gt.0.707) nn2=nn2+1 
2       continue 
        sum=sum*6.83  
        psum=psum*6.83  
        sumf=amax1(sum,psum)
        if(sumf.lt.contrs) sumf=contrs                                         
        if(etot.gt.0.005)  sumf=contrs  
	if(sumf.gt.RANFL(mseed)) ntag=-1
        return                                                                
        end
c///////////////////////////////////////////////////////////////////////////////
	      	SUBROUTINE initp
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	common/paulis/radn,padn,anorm,const,radn2,padn2                       
	common/poli/iflagcro,iflagpau,ipoli
        pi=3.14159                                                            
        aux=(4.*pi/3.)**2                                                     
        hc=1.24            
C
   	if(iflagpau.eq.1) then
      const=4./9.*pi*pi*4./(hc)**3./2.
 	radn=3.367
      padn=hc/radn/aux**(1./3.)/2.**(1./3.)             
	else 
	const=4./9.*pi*pi*4./(hc)**3.
	radn=3.
	padn=hc/radn/aux**.333/4.**.333
	endif
C
        radn2=(2.*radn)**2                                                    
        padn2=(2.*padn)**2     
	                                                
233     format(1x,'pauli radii   r = ',f6.2,'   p = ',f6.2)                   
        return                                                                
        end             
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE init 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/rp0/xi(500),yi(500),zi(500),pxi(500),pyi(500),pzi(500),
     &  emi(500),isoi(500),isyi(500),idi(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/bb/np,nzp,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/gg/nzta,nzpr                                                     
        common/ee/mseed,nseed
        common/skyrme/al,sal,t0,t1,arang,vy,width                                                                            
        common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	  common/max/mden
	  common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
        dimension dvdpx(500),dvdpy(500),dvdpz(500)                              
	  real mden
        hc=1.24                                                        
        vpaul=v01                                                               
        hb=0.197327                                                     
        widthp=p01                                                              
        widthq=q01                                                              
        wp2=2.*widthp**2                                                        
        wq2=2.*widthq**2                                                        
        anorm=vpaul*(hb/widthp/widthq)**3                                       
        x0=x00/2.                                                               
        pav=0.                                                                  
        ddt=0.                                                                  
        gami=1.                                                                 
	ntstep=0
        do 1 i=1,nzpr                                                           
13      r11=rrp*RANFL(mseed)
        r12=rrp*RANFL(mseed)                                                     
        r13=rrp*RANFL(mseed)                                                  
        r=amax1(r11,r12,r13)                                                    
        cth=1.-2.*RANFL(mseed)                                                
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*RANFL(mseed)                                                 
        z(i)=(r*cth)/gami-z0                                            
        x(i)=r*sth*cos(phi)+x0
        y(i)=r*sth*sin(phi)                                                     
        iso(i)=1                                                                
	isy(i)=-1
        id(i)=1                                                                 
        do 12 j=1,i-1                                                           
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
        if(d.lt.1.5)go to 13
12      continue                                                                
1       continue                                                                
        do 11 i=nzpr+1,npr                                                      
155     r11=rrp*RANFL(mseed)            
        r12=rrp*RANFL(mseed)            
        r13=rrp*RANFL(mseed)            
        r=amax1(r11,r12,r13)                                                    
        cth=1.-2.*RANFL(mseed)          
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*RANFL(mseed)           
        z(i)=(r*cth)/gami-z0                                       
        x(i)=r*sth*cos(phi)+x0
        id(i)=1                                                                 
        y(i)=r*sth*sin(phi)                                                     
        iso(i)=0                                                                
	isy(i)=1
        do 14 j=1,i-1                                                           
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
        if(d.lt.1.5)go to 155                    
14      continue                                                                
11      continue                                                                
450     format(1x,'npr=',i6,'  imiss=',i6)                                      
        do 111 i=npr+1,npr+nzta                                                 
19      r11=rrt*RANFL(mseed)
        r12=rrt*RANFL(mseed)       
        r13=rrt*RANFL(mseed)       
        r=amax1(r11,r12,r13)                                                    
        cth=1.-2.*RANFL(mseed)     
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*RANFL(mseed)
        z(i)=(r*cth)/gami+z0                                         
        x(i)=r*sth*cos(phi)-x0
        iso(i)=1                                                                
	isy(i)=-1
        id(i)=51                                                                
        y(i)=r*sth*sin(phi)                                                     
        do 18 j=npr+1,i-1                                                       
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
        if(d.lt.1.5)go to 19
18      continue                                                                
111     continue                                                                
        do 1111 i=npr+nzta+1,np                                                 
17      r11=rrt*RANFL(mseed)        
        r12=rrt*RANFL(mseed)
        r13=rrt*RANFL(mseed) 
        r=amax1(r11,r12,r13)                                                   
        cth=1.-2.*RANFL(mseed)
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*RANFL(mseed)
        z(i)=(r*cth)/gami+z0                                       
        x(i)=r*sth*cos(phi)-x0
        iso(i)=0                                                                
	  isy(i)=1
        id(i)=51                                                                
        y(i)=r*sth*sin(phi)                                                     
        do 16 j=npr+1,i-1                                                       
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
        if(d.lt.1.5)goto17                                                      
16      continue                                                                
1111    continue                                                                
451     format(1x,'nta = ',i6,'  imiss = ',i6)                                                                                              
9100    format(1x,'initial density')                                            
9200    format(1x,10f9.4)                                                       

	p0z=pzpr
	p0x=pxpr

        if(npr.ne.0)dpp=(0.008/float(npr)*3.4)**0.666
        if(nta.ne.0)dpt=(0.008/float(nta)*3.4)**0.666
        if(npr.ne.0.and.npr.lt.30)dpp=(0.008/float(npr)*4.3)**0.666
        if(nta.ne.0.and.nta.lt.30)dpt=(0.008/float(nta)*4.3)**0.666
        dpc=dpp                                                                 
        atz=0.                                                                  
        apz=0.                                                                  
        bpz=0.                                                                  
        btz=0.                                                                  
        atx=0.                                                                  
        apx=0.                                                                  
        bpx=0.                                                                  
        btx=0.                                                                  
        apy=0.                                                                  
        aty=0.                                                                  
        bpy=0.                                                                  
        bty=0.                                                                  
        phmax=0.280**2                                                          
        do 36 i=1,np                                                            
        dens1=densit(x(i),y(i),z(i))                                                         
        ddt=ddt+dens1/float(np)                                              
        if(densty.gt.0.16)densty=0.16                                           
        aux=(3.*densty)/(16.*3.1415)                                            
        if(aux.lt.0)aux=0.                                                      
        pfl=hc*aux**.3333                                                       
352     pr=(pfl-0.02)*RANFL(mseed)**.3333 
        if(pr.lt.0)pr=0.                 
        cth=1.-2.*RANFL(mseed)            
        sth=sqrt(1.-cth*cth)                                                    
        phi=6.28318*RANFL(mseed)          
        if(i.ne.1+npr) goto 35                                                  

	p0z=-pzta
	p0x=-pxta

        phmax=0.28**2                                                           
        dpc=dpt                                                                 
35      pz(i)=p0z+pr*cth                                                        
        px(i)=p0x+pr*sth*cos(phi)
        py(i)=pr*sth*sin(phi)                                                   
        do 357 j=1,i-1                                                          
        dp=(px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2                   
        d=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2                          
        if(dp*d.lt.phmax)goto 352                                               
357     continue                                                                
3666    pav=pav+pr/float(np)                                                    
        em(i)=0.938                                                             
        if(id(i).eq.1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        endif                                                                   
        if(id(i).eq.51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        endif                                                                   
   36   continue                                                                
452     format(1x,'imiss phassp=',i6)                                           
453     format(1x,'average density= ',f10.4,' average fermi momentum            
     +  = ',f10.4)                                                              
        do15 i=1,np                                                             
        if(npr.eq.0.or.nta.eq.0) then
		z0=0.
		x0=0.
		pzpr=0.
		pxpr=0.
		pzta=0.
		pxta=0.

        endif 
        if(id(i).eq.1) then                                                     
                      pz(i)=pz(i)-bpz+pzpr
                      z(i)=z(i)-apz-z0
                      px(i)=px(i)-bpx+pxpr                                    
                      x(i)=x(i)-apx+x0
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
	if(iflag1.lt.0)then
	xi(i)=x(i)
	yi(i)=y(i)
	zi(i)=z(i)
	pxi(i)=px(i)
	pyi(i)=py(i)
	pzi(i)=pz(i)
	emi(i)=em(i)
	isoi(i)=iso(i)
	isyi(i)=isy(i)
	idi(i)=id(i)
	endif
        endif                                                                  
        if(id(i).eq.51) then                                                    
                       pz(i)=pz(i)-btz-pzta
                       z(i)=z(i)-atz+z0
                       px(i)=px(i)-btx-pxta                                   
                       x(i)=x(i)-atx-x0
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
	if(iflag1.lt.0)then
	xi(i)=x(i)
	yi(i)=y(i)
	zi(i)=z(i)
	pxi(i)=px(i)
	pyi(i)=py(i)
	pzi(i)=pz(i)
	emi(i)=em(i)
	isoi(i)=iso(i)
	isyi(i)=isy(i)
	idi(i)=id(i)
        endif 
        endif 
77777	format(1x,7f9.4,2x,2i3,2x,i3)	
88888	format(1x,2f9.4)
15      continue                                                                
345     format(1x,'average excitation energy = ', f10.4)                        
	if((iflag1.lt.0.and.initi.eq.1).or.iflag1.gt.0)then
	do 98721 i=1,np
	write(71,87731)z(i),x(i)
87731	format(1x,2f9.4)	
98721	continue
	endif
        dt=dt/2.                                                                
        if(iflag5.eq.0)goto765                                                  
        factor=4*0.00157*500.                                                   
        do 8010 i1=1,np                                                         
        dvdpx(i1)=0.                                                            
        dvdpy(i1)=0.                                                            
        dvdpz(i1)=0.                                                            
8010    continue                                                                
        do 8011 i1=1,np                                                         
        do 8012 i2=i1+1,np                                                      
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     &  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        vdens=XEXP(-dr2/al)                                                      
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr/sal                                     
        dvdp1=-anorm/widthp**2*XEXP(-dr2/wq2)*XEXP(-p2/wp2)                       
        dvdpx(i1)=dvdpx(i1)+(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i1)=dvdpy(i1)+(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i1)=dvdpz(i1)+(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
        dvdpx(i2)=dvdpx(i2)-(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i2)=dvdpy(i2)-(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i2)=dvdpz(i2)-(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
8012    continue                                                                
8011    continue                                                                
c                                                                               
c                                                                               
        do 8013 i=1,np                                                          
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1+dvdpx(i)*dt                                       
        y(i)=y(i)+dt*py(i)/e1+dvdpy(i)*dt                                       
        z(i)=z(i)+dt*pz(i)/e1+dvdpz(i)*dt                                       
8013    continue                                                                
        goto 766                                                                
                                                                                
765     do 8014 i=1,np                                                          
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1
8014    continue
766     dt=dt*2.                                                                
        return                                                                  
        end             	                                                                                                                                                                       
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE init2 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	character*50 targ, proj

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,               
     & rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     & idr(500),xms(500),yms(500),zms(500),pxms(500),
     & pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp,nanti
      common/bb/np,nzp,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/skyrme/al,sal,t0,t1,arang,vy,width 
      common/gg/nzta,nzpr,nmint
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/ee/mseed,nseed
	common/tt/prtime(8000),nrmin,nrmax,rtime,perc
	common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
	common/data/targ, proj
      integer target_orient_method
      common /target_orient/ target_orient_method
      real orient_matrix(3,3)
      dimension dvdpx(500),dvdpy(500),dvdpz(500) 
		
        hc=  1.24 
        vpaul=v01 
        hb=0.197327  
        widthp=p01 
        widthq=q01                                                              
        wp2=2.*widthp**2                                                        
        wq2=2.*widthq**2                                                        
        anorm=vpaul*(hb/widthp/widthq)**3                                       
        p0=ppo                                                                  
        atz=0.
        apz=0.
        bpz=0.                                                                  
        btz=0.                                                                  
        atx=0.                                                                  
        apx=0.                                                                  
        bpx=0.                                                                  
        btx=0.                                                                  
        apy=0.                                                                  
        aty=0.                                                                  
        bpy=0.                                                                  
        bty=0.

	 z00=2.*z0
	 zp=z00*float(nta)/float(np)                      !fengzhq
	 zt=z00*float(npr)/float(np)

	if(npr.le.nmint) then
	  do i=1,npr
			rr0=1.14*float(npr)**0.3333*RANFL(mseed)**0.3333
			cth=1.-2.*RANFL(mseed)
			sth=sqrt(1.-cth**2)
			phi=6.28*RANFL(mseed)
			z(i)=rr0*cth                                         
			x(i)=rr0*sth*cos(phi)
			y(i)=rr0*sth*sin(phi)
			if(nanti .EQ. 0) then
			if(i.gt.0.and.i.le.nzpr)   iso(i)=1
			else
			if(i.gt.0.and.i.le.nzpr)   iso(i)=-1
			endif
			if(i.gt.nzpr.and.i.le.npr) iso(i)=0
			em(i)=0.938
			pr0=0.262*RANFL(mseed)**0.3333
			pz(i)=pr0*cth                                                        
			px(i)=pr0*sth*cos(phi)
			py(i)=pr0*sth*sin(phi)

			bpz=bpz+pz(i)/npr                                                       
			apz=apz+z(i)/npr                                                        
			bpx=bpx+px(i)/npr                                                       
			apx=apx+x(i)/npr 
			bpy=bpy+py(i)/npr 
			apy=apy+y(i)/npr 
	  enddo

	else

	  open(30,file=proj,form='formatted',status='old')

        do 10 i=1,npr
          read(30,1000) x(i),y(i),z(i),px(i),py(i),pz(i),em(i),iso(i),
     &      isy(i),id(i)
			bpz=bpz+pz(i)/npr                                                       
			apz=apz+z(i)/npr                                                        
			bpx=bpx+px(i)/npr                                                       
			apx=apx+x(i)/npr                                                        
			bpy=bpy+py(i)/npr                                                       
			apy=apy+y(i)/npr                                                        
 10     continue                                                                
	  close(30)
	endif

        do 11 i=1,npr                                                           
        pz(i)=pz(i)-bpz
        z(i)=z(i)-apz
        px(i)=px(i)-bpx
        x(i)=x(i)-apx
        py(i)=py(i)-bpy                                                         
        y(i)=y(i)-apy    
        id(i)=1     
	  idr(i)=0
	  idanti(i)=0
	  if(nanti .NE. 0)	 idanti(i)=1              ! for antinucleon
	if(i.gt.0.and.i.le.nzpr)   isy(i)=-1
	if(i.gt.nzpr.and.i.le.npr) isy(i)=1
  11    continue
        
	p0=-tpo                                                                 
	if(nta.ne.0.and.nzta.ne.0)then
	open(10,file=targ,form='formatted',status='old')
      read(10,*)
        do 20 i=npr+1,np                                                        
        read(10,1000) x(i),y(i),z(i),px(i),py(i),pz(i),em(i),iso(i),
     &  isy(i),id(i)
        btz=btz+pz(i)/nta                                                       
        atz=atz+z(i)/nta                                                        
        btx=btx+px(i)/nta                                                       
        atx=atx+x(i)/nta                                                        
        bty=bty+py(i)/nta                                                       
        aty=aty+y(i)/nta                                                        
  20    continue  


	                                                              
	close(10)

        do 21  i=npr+1,np                                                       
        pz(i)=pz(i)-btz
        z(i)=z(i)-atz                                                    
        px(i)=px(i)-btx
        x(i)=x(i)-atx                                                        
        py(i)=py(i)-bty                                                         
        y(i)=y(i)-aty                                                           
        id(i)=51
	  idr(i)=0
	  idanti(i)=0
	if(i.gt.npr.and.i.le.npr+nzta) isy(i)=-1
	if(i.gt.npr+nzta.and.i.le.np)  isy(i)=1
 21     continue
	endif


      if(iflag1.gt.0.and.((target_orient_method.eq.1.and.
     &   nru.gt.nrmin).or.(target_orient_method.eq.2.and.nru.gt.1)))then
	cth=1.-2.*RANFL(mseed)
	sth=sqrt(1.-cth**2.)
	phi=6.28*RANFL(mseed)
	do 1 i=1,npr
	xx1=x(i)
	yy1=y(i)
	zz1=z(i)
	x(i)=xx1*cos(phi)+(yy1*cth+zz1*sth)*sin(phi)
	y(i)=-xx1*sin(phi)+(yy1*cth+zz1*sth)*cos(phi)
	z(i)=-yy1*sth+zz1*cth
	px1=px(i)
	py1=py(i)
	pz1=pz(i)
	px(i)=px1*cos(phi)+(py1*cth+pz1*sth)*sin(phi)
	py(i)=-px1*sin(phi)+(py1*cth+pz1*sth)*cos(phi)
	pz(i)=-py1*sth+pz1*cth
1	continue
	print*,'rotate proj. two euler angles!'
	!write(9,*) 'rotate proj. two euler angles!'
	endif

        do 3 i=1,npr 
        pz(i)=pz(i)+pzpr
        z(i)=z(i)-zp
        px(i)=px(i)+pxpr
        x(i)=x(i)+x0
        py(i)=py(i)
        y(i)=y(i)
	  id(i)=1
 3	  continue



      if(iflag1.gt.0.and.((target_orient_method.eq.1.and.
     &   nru.gt.nrmin).or.(target_orient_method.eq.2.and.nru.gt.1)))then
      if(target_orient_method.eq.2) then
        print*,'rotate targ. Haar SO(3)!'
        call target_orientation_matrix(mseed,orient_matrix)
        do i=npr+1,np
          call target_orientation_apply(orient_matrix,x(i),y(i),z(i))
          call target_orientation_apply(orient_matrix,px(i),py(i),pz(i))
        enddo
      else
	cth=1.-2.*RANFL(mseed)
	sth=sqrt(1.-cth**2.)
	phi=6.28*RANFL(mseed)
	do 10004 i=npr+1,np
	xx1=x(i)
	yy1=y(i)
	zz1=z(i)
	x(i)=xx1*cos(phi)+(yy1*cth+zz1*sth)*sin(phi)
	y(i)=-xx1*sin(phi)+(yy1*cth+zz1*sth)*cos(phi)
	z(i)=-yy1*sth+zz1*cth
	px1=px(i)
	py1=py(i)
	pz1=pz(i)
	px(i)=px1*cos(phi)+(py1*cth+pz1*sth)*sin(phi)
	py(i)=-px1*sin(phi)+(py1*cth+pz1*sth)*cos(phi)
	pz(i)=-py1*sth+pz1*cth
10004	continue
	print*,'rotate targ. two euler angles!'
      endif
	!write(9,*) 'rotate targ. two euler angles!' 
	endif


C C 	print*, X0,ZT,PXTA,PZTA
C C 	PAUSE

 

        do 4  i=npr+1,np	
        pz(i)=pz(i)-pzta
        z(i)=z(i)+zt                                                    
        px(i)=px(i)-pxta
        x(i)=x(i)-x0 
        py(i)=py(i)
        y(i)=y(i)
        id(i)=51                                                            
 4	  continue




	do i=1,np
	idr(i)=0
	xms(i)=0.
	yms(i)=0.
	zms(i)=0.
	pxms(i)=0.
	pyms(i)=0.
	pzms(i)=0.
	emson(i)=0.
	idms(i)=0
	icp(i)=0
	enddo


	dt=dt/2. 

        if(iflag5.eq.0) goto 765
        factor=4*0.00157*500.                                                   
        do 8010 i1=1,np                                                         
        dvdpx(i1)=0.                                                            
        dvdpy(i1)=0.                                                            
        dvdpz(i1)=0.                                                            
8010    continue                                                                
        do 8011 i1=1,np                                                         
        do 8012 i2=i1+1,np                                                      
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     &  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        vdens=XEXP(-dr2/al)                                                      
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr/sal                                     
        dvdp1=-anorm/widthp**2*XEXP(-dr2/wq2)*XEXP(-p2/wp2)                       
        dvdpx(i1)=dvdpx(i1)+(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i1)=dvdpy(i1)+(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i1)=dvdpz(i1)+(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
        dvdpx(i2)=dvdpx(i2)-(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i2)=dvdpy(i2)-(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i2)=dvdpz(i2)-(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
8012    continue                                                                
8011    continue                                                                
        do 8013 i=1,np                                                          
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1+dvdpx(i)*dt                                       
        y(i)=y(i)+dt*py(i)/e1+dvdpy(i)*dt                                       
        z(i)=z(i)+dt*pz(i)/e1+dvdpz(i)*dt                                     
8013    continue 
        goto 766 

        do 8014 i=1,np
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1                                                   
8014    continue
765	continue !!HG CHENG
766     dt=dt*2.

1000    format(1x,7f9.4,2x,2i3,2x,i3) 
10000   format(1x,i20,1x,i4,4(f8.2,1x))                                         
10001   format(1x,4(i3,1x))                                                     
        return                                                                  
        end                                                                        

c///////////////////////////////////////////////////////////////////////////////
      function densit(x1,y1,z1)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &    rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,nzp,ntag,x00                                                                         
      common/wavel/sl(500)
	common/SURF/ VANISH(500)


    	pi=3.141592654
	dens=0. 

        do 10 i=1,np
		if(VANISH(I).EQ.1)	cycle 
	  width=2.*sl(i) 
        dr2=(x1-x(i))**2.+(y1-y(i))**2.+(z1-z(i))**2. 
        dd=XEXP(-dr2/width)/(pi*width)**1.5 
        dens=dens+dd 
10	  continue 
	  densit=dens 
        return                                                     
        end
	      
c///////////////////////////////////////////////////////////////////////////////
      function denstn(x1,y1,z1)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &    rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,nzp,ntag,x00                                                                         
      common/wavel/sl(500)
	common/SURF/ VANISH(500)


    	pi=3.141592654
	dens=0. 

        do 10 i=1,np 
		if(VANISH(I).EQ.1)	cycle
	  width=2.*sl(i) 
        dr2=(x1-x(i))**2.+(y1-y(i))**2.+(z1-z(i))**2. 
        dd=XEXP(-dr2/width)/(pi*width)**1.5 
        if(isy(i). EQ. 1)  dens=dens+dd 
10	  continue 
	  denstn=dens 
        return                                                     
        end

c///////////////////////////////////////////////////////////////////////////////
      function denstp(x1,y1,z1)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &    rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,nzp,ntag,x00                                                                         
      common/wavel/sl(500)
	common/SURF/ VANISH(500)


    	pi=3.141592654
	dens=0. 

        do 10 i=1,np
		if(VANISH(I).EQ.1)	cycle 
	  width=2.*sl(i) 
        dr2=(x1-x(i))**2.+(y1-y(i))**2.+(z1-z(i))**2. 
        dd=XEXP(-dr2/width)/(pi*width)**1.5 
        if(isy(i). EQ. -1)  dens=dens+dd 
10	  continue 
	  denstp=dens 
        return                                                     
        end
c///////////////////////////////////////////////////////////////////////////////
      		SUBROUTINE densitn                                          
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00  
        common/densi1/denstl(500),dens1(500),velcm(3)
	  common/wavel/sl(500) 
	common/SURF/ VANISH(500) 
	  pi=3.141592654                                                         
        do 10 i=1,np 
	       if(VANISH(I).EQ.1)	cycle
	  denstl(i)=0.
        dens1(i)=0.                                               
10      continue                                                   
        do 20 i=1,np 
			if(VANISH(I).EQ.1)	cycle                                              
        do 30 j=1,np  
			if(VANISH(J).EQ.1)	cycle
	  width=sl(i)+sl(j)
        dr2=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2  
        dd=xexp(-dr2/width) 
        denstl(i)=denstl(i)+dd/(pi*width)**1.5
	  if(iso(i).eq.iso(j)) dens1(i)=dens1(i)+dd/(pi*width)**1.5
30   	  continue                                                                     
20   	  continue  
                                                 
        return                                                     
        end  	                                                     
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE densim(densm)  
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        dimension densm(250)                                                  
        widthm=.1973289**2*2./al                                              
        anorm=1./(3.141592*widthm)**1.5/np                                    
        widthr=al/2.                                                          
        anorm2=1./1.24**3                                                     
        do 10 i=1,np                                                          
10      densm(i)=0.0                                                          
        do 20 i=1,np                                                          
        do 30 j=i+1,np                                                        
        dp2=(px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2                
        dd=XEXP(-dp2/widthm)                                                   
        densm(i)=densm(i)+dd                                                  
        densm(j)=densm(j)+dd                                                  
30      continue                                                              
        densm(i)=anorm*(1.+densm(i))                                          
20      continue                                                              
        return                                                                
        end                                                                   

c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE dendgn 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp
        common/bb/np,nzp,ntag,x00
        common/densi/err(302),dens(50,50,50)                          
        dimension ex(5),ey(5),ez(5)                                   
        do 4 j=1,50                                                   
        do 4 jj=1,50                                                  
        do 4 jjj=1,50                                                 
4       dens(j,jj,jjj)=0.                                             
        do 1 k1=1,np                                                  
        lx=aint(x(k1)+25)                                             
        ly=aint(y(k1)+25)                                             
        lz=aint(z(k1)+25)                                             
        if(lx.lt.3.or.lx.gt.47)goto1                                  
        if(ly.lt.3.or.ly.gt.47)goto1                                  
        if(lz.lt.3.or.lz.gt.47)goto1                                  
        dx=abs(x(k1))-aint(abs(x(k1)))                                
        idx=anint(dx*100)                                             
        iddx=anint((1.-dx)*100)                                       
        ex(1)=err(idx+201)-err(idx+101)                               
        ex(2)=err(idx+101)-err(idx+1)                                 
        ex(3)=err(idx+1)+err(iddx+1)                                  
        ex(4)=err(iddx+101)-err(iddx+1)                               
        ex(5)=err(iddx+201)-err(iddx+101)                             
        dy=abs(y(k1))-aint(abs(y(k1)))                                
        idy=anint(dy*100)                                             
        iddy=anint((1.-dy)*100)                                       
        ey(1)=err(idy+201)-err(idy+101)                               
        ey(2)=err(idy+101)-err(idy+1)                                 
        ey(3)=err(idy+1)+err(iddy+1)                                  
        ey(4)=err(iddy+101)-err(iddy+1)                               
        ey(5)=err(iddy+201)-err(iddy+101)                             
        dz=abs(z(k1))-aint(abs(z(k1)))                                
        idz=anint(dz*100)                                             
        iddz=anint((1.-dz)*100)                                       
        ez(1)=err(idz+201)-err(idz+101)                               
        ez(2)=err(idz+101)-err(idz+1)                                 
        ez(3)=err(idz+1)+err(iddz+1)                                  
        ez(4)=err(iddz+101)-err(iddz+1)                               
        ez(5)=err(iddz+201)-err(iddz+101)                             
        do 3 mx=1,5                                                   
        do 3 my=1,5                                                   
        do 3 mz=1,5                                                   
        mmx=lx-3+mx                                                   
        mmy=ly-3+my                                                   
        mmz=lz-3+mz                                                   
3       dens(mmx,mmy,mmz)=dens(mmx,mmy,mmz)+ex(mx)*ey(my)*ez(mz)      
1       continue                                                      
        return                                                        
        end                                                           

c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE prot
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	character idate*9,itime*8        
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     & rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     & idr(500),xms(500),yms(500),zms(500),pxms(500),
     & pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
	common/rpi/rmst(500),pmst(500),avet(500)
      common/ee/mseed,nseed
      common/bb/np,nzp,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                         
      common/kaot/pktot,nkn                                                                
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
      common/deu2/px1,py1,pz1,px2,py2,pz2                                     
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
	common/cmod/coeft,coefp,epp,ecm,ptcm
	
        hb=0.197327
        pp=0.                                                                   
        tp=0.                                                                   
        pr=0.                                                                   
        tr=0.
	  rr=0.                                                                          
        bpz=0.                                                                  
        btz=0.                                    
        bpx=0.                                                                  
        btx=0.                                                                         
        bpy=0.                                                                  
        bty=0.   
	ndel=0
	ndel1=0
	ndel2=0
	ndel3=0
	ndel4=0
	nstr=0
	nstr1=0
	nstr2=0
	npion=0
	npion1=0
	npion2=0
	npion3=0
        do 136 i=1,np                                                           
        if(id(i).lt.50) then                                                  
                   bpz=bpz+pz(i)/npr                                                  
                   bpx=bpx+px(i)/npr                                                           
                   bpy=bpy+py(i)/npr                                           
        endif 
        if(id(i).gt.50) then                                                    
                   btz=btz+pz(i)/nta                                            
                   btx=btx+px(i)/nta                                              
                   bty=bty+py(i)/nta                                                
        endif 
136     continue    
	  call center(apx,apy,apz,1)
	  call center(atx,aty,atz,2)     
	                                                         
	  bcmx=(apx*npr+atx*nta)/float(np)                                        
        bcmy=(apy*npr+aty*nta)/float(np)                                        
        bcmz=(apz*npr+atz*nta)/float(np) 
	  rpp= apx**2.+apy**2.+apz**2.
	  rtt= atx**2.+aty**2.+atz**2.   
	  rpt=(apx-atx)**2.+(apy-aty)**2.+(apz-atz)**2.
	  rpp=sqrt(rpp)
	  rtt=sqrt(rtt)
	  rpt=sqrt(rpt)
        do 137 i=1,np
        if (id(i) .lt. 50.and.nta.eq.0.and.epp.eq.0.) then
               pr = pr + (x(i)-apx)**2.+(y(i)-apy)**2.+(z(i)-apz)**2.        
               pp = pp + (px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.     
        endif 
        if (id(i) .gt. 50.and.npr.eq.0.and.epp.eq.0.) then
               tr = tr + (x(i)-atx)**2.+(y(i)-aty)**2.+(z(i)-atz)**2.        
               tp = tp + (px(i)-btx)**2.+(py(i)-bty)**2.+(pz(i)-btz)**2.     
        endif
        rr=rr+(x(i)-bcmx)**2.+(y(i)-bcmy)**2.+(z(i)-bcmz)**2. 
	
	if(idr(i).eq.1) then
	ndel=ndel+1
	if(iso(i).eq.2) ndel1=ndel1+1
	if(iso(i).eq.1) ndel2=ndel2+1
	if(iso(i).eq.0) ndel3=ndel3+1
	if(iso(i).eq.-1)ndel4=ndel4+1
	elseif(idr(i).eq.2) then
	nstr=nstr+1
	if(iso(i).eq.1) nstr1=nstr1+1
	if(iso(i).eq.0) nstr2=nstr2+1
	endif
	if(idms(i).eq.1) then
	npion=npion+1
	if(icp(i).eq.1) npion1=npion1+1
	if(icp(i).eq.0) npion2=npion2+1
	if(icp(i).eq.-1)npion3=npion3+1
	endif
137     continue
*	print*, ntstep, npion1, npion2, npion3
	llt=ntstep/kfile+1

	if(npr.gt.0.and.nta.eq.0.and.epp.eq.0.)then
	pr=pr/float(npr)
      pr= sqrt(pr)
	rmst(llt)=pr
      prtime(ntstep) = prtime(ntstep) + pr                                    
	endif

	if(nta.gt.0.and.npr.eq.0.and.epp.eq.0.)then
	tr=tr/float(nta)
      tr=sqrt(tr)
	rmst(llt)=tr
	endif

	if(nta.gt.0.and.npr.gt.0.and.epp.gt.0) then
	rr=rr/float(np) 
	rr=sqrt(rr)
	rmst(llt)=rr
	endif
                                                  
        return                                                                  
        end                                                                     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        SUBROUTINE collnum  
C//////////////////////////////////////////////////////////////////////////

        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,
     &	ncolln,ncollp,nblkn,nblkp
        common/bb/np,nzp,ntag,x00                                                   
	common/gg/nzta,nzpr
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
	common/cmod/coeft,coefp,epp,ecm,ptcm
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	common/poli/iflagcro,iflagpau,ipoli
	common/totcol/ipuli3,iplin3,iplip3
	dimension ipuli1(0:3000),ipuli2(0:3000),ipulit(0:3000),
     +	ratio1(0:3000),
     +	iplin1(0:3000),iplin2(0:3000),iplin(0:3000),ration(0:3000),
     +	iplip1(0:3000),iplip2(0:3000),iplip(0:3000),ratiop(0:3000)

	frun=float(nrmax-nrmin+1)
	fnp=float(np)
	fnz=float(nzta+nzpr)
	fnn=fnp-fnz

	if(nru.eq.1.and.rtime.eq.0.) then
c       write(55,*) ' \\E= ',epp,' MeV/u'
c       write(55,*) ' \\b= ',x00,' fm'
	if(iflagcro.le.0) then
c	write(55,*) ' \\Cugnon''s N-N cross section !!'
	else if(iflagcro.ge.1) then
c	write(55,*) ' \\Experimental N-N cross section !!'
	endif 
       if(ipoli.eq.0.or.ipoli.eq.1) then
	if(ipoli.eq.0) then
c       write(55,*) ' \\All collisions are blocked !!'
	endif
	if(ipoli.eq.1) then
c      write(55,*) ' \\No Pauli blocking !!'
	endif
	goto 96322
       endif
c       if(iflagpau.le.0) write(55,*) ' \\Old Pauli blocking !!'
c       if(iflagpau.eq.1) write(55,*) ' \\R34 Pauli blocking !!'
c       if(iflagpau.eq.2) write(55,*) ' \\Over Pauli blocking !!'
c       if(iflagpau.eq.3) write(55,*) ' \\Rho Pauli blocking !!'
96322  continue
c       write(55,*) ' \\Time(fm/c) Col.Numb Blk.Numb RealColl BC.Ratio
c     + rcolltot'
	endif
C--------------------------------------------------------------------
	iplin1(int(rtime))=iplin1(int(rtime))+ncolln
	iplin2(int(rtime))=iplin2(int(rtime))+nblkn
	iplin(int(rtime))=iplin(int(rtime))+iplin3
	iplip1(int(rtime))=iplip1(int(rtime))+ncollp
	iplip2(int(rtime))=iplip2(int(rtime))+nblkp
	iplip(int(rtime))=iplip(int(rtime))+iplip3
	ipuli1(int(rtime))=ipuli1(int(rtime))+ncoll
	ipuli2(int(rtime))=ipuli2(int(rtime))+nblock
	ipulit(int(rtime))=ipulit(int(rtime))+ipuli3
C
	if(nru.eq.nrmax) then
	if(ipuli1(int(rtime)).eq.0) then
	ratio1(int(rtime))=1.
	else
	ratio1(int(rtime))=float(ipuli2(int(rtime)))/
     +	float(ipuli1(int(rtime)))
	endif
C
	if(iplin1(int(rtime)).eq.0) then
	ration(int(rtime))=1.
	else
	ration(int(rtime))=float(iplin2(int(rtime)))/
     +	float(iplin1(int(rtime)))
	endif
C
	if(iplip1(int(rtime)).eq.0) then
	ratiop(int(rtime))=1.
	else
	ratiop(int(rtime))=float(iplip2(int(rtime)))/
     + 	float(iplip1(int(rtime)))
	endif
C
c	write(55,96515) rtime,ipuli1(int(rtime))/frun/fnp/dt,
c     +ipuli2(int(rtime))/frun/fnp/dt,
c     +float(ipuli1(int(rtime))-ipuli2(int(rtime)))/frun/fnp/dt,
c     +ratio1(int(rtime)),
c     +float(ipulit(int(rtime)))/frun/fnp/dt
C
c	write(56,96515) rtime,iplin1(int(rtime))/frun/fnn/dt,
c     +iplin2(int(rtime))/frun/fnn/dt,
c     +float(iplin1(int(rtime))-iplin2(int(rtime)))/frun/fnn/dt,
c     +ration(int(rtime)),
c     +float(iplin(int(rtime)))/frun/fnn/dt
C
c	write(57,96515) rtime,iplip1(int(rtime))/frun/fnz/dt,
c     +iplip2(int(rtime))/frun/fnz/dt,
c     +float(iplip1(int(rtime))-iplip2(int(rtime)))/frun/fnz/dt,
c     +ratiop(int(rtime)),
c     +float(iplip(int(rtime)))/frun/fnz/dt
	endif
!96515	format(1x,f10.2,2(1x,f8.4),1x,f8.4,2x,f7.5,2x,f9.4)
	return
	end
c///////////////////////////////////////////////////////////////////////////////
      		SUBROUTINE updt 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &  ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        dimension f1x(500),f1y(500),f1z(500),t2(500),dr2(500) 
	  common/SURF/ VANISH(500)                  
c                                                                               
        do 1 i=1,np  
		if(VANISH(I).EQ.1)	cycle                                                           
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1                                                   
        f1x(i)=0.                                                               
        f1y(i)=0.                                                               
        f1z(i)=0.                                                               
1       fr(i)=0.                                                                
c                                                                               
        do 2 i1=1,np 
		if(VANISH(I1).EQ.1)	cycle                                                           
        do 21 i2=i1+1,np  
		if(VANISH(I2).EQ.1)	cycle                                                      
        dr2(i2)=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))         
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        t2(i2)=XEXP(-dr2(i2)/al)                                                 
21      continue                                                                
        do 22 i2=i1+1,np 
		if(VANISH(I2).EQ.1)	cycle                                                       
        f1x(i1)=f1x(i1)+(x(i1)-x(i2))*t2(i2)                                    
        f1x(i2)=f1x(i2)+(x(i2)-x(i1))*t2(i2)                                    
        f1y(i1)=f1y(i1)+(y(i1)-y(i2))*t2(i2)                                    
        f1y(i2)=f1y(i2)+(y(i2)-y(i1))*t2(i2)                                    
        f1z(i1)=f1z(i1)+(z(i1)-z(i2))*t2(i2)                                    
22      f1z(i2)=f1z(i2)+(z(i2)-z(i1))*t2(i2)                                    
        sum=0.                                                                  
        do 23 i2=i1+1,np 
	  if(VANISH(I2).EQ.1)	cycle                                                       
        sum=sum+t2(i2)                                                          
        fr(i2)=fr(i2)+t2(i2)                                                    
23      continue                                                                
        fr(i1)=fr(i1)+sum                                                       
2       continue                                                                
        do 3 i=1,np
	  if(VANISH(I).EQ.1)	cycle                                                             
        st=pota/sal/al*(1.-potb*(fr(i)/sal)**potc)                              
        fx=st*f1x(i)                                                            
        fy=st*f1y(i)                                                            
        fz=st*f1z(i)                                                            
        px(i)=px(i)+dt*fx                                                       
        py(i)=py(i)+dt*fy                                                       
3       pz(i)=pz(i)+dt*fz                                                       
        return                                                                  
        end 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        	SUBROUTINE updtyc  
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,               
     &        potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/coul/cy(5000),cyu(5000),cypion(5000)                                          
	  common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +	frtwo(500)
	common/SURF/ VANISH(500)

        dimension f1x(500),f1y(500),f1z(500),t2(500),t2new(500) 
        dimension fcx(500),fcy(500),fcz(500),u(500)                             
c       dimension fhx(500),fhy(500),fhz(500),upau(i2)                           
        dimension fsx(500),fsy(500),fsz(500)
c	dimension s(500)
        do 1 i=1,np         
		if(VANISH(I).EQ.1)	cycle
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1                                                   
        f1x(i)=0.                                                               
        f1y(i)=0.                                                               
        f1z(i)=0.                                                               
        fcx(i)=0.                                                               
        fcy(i)=0.                                                               
        fcz(i)=0.                                                               
	fsx(i)=0.
	fsy(i)=0.
	fsz(i)=0.
	frtwo(i)=0.
1       fr(i)=0.                                                                
c                                                                               
        do 2 i1=1,np 
	  if(VANISH(I1).EQ.1)	cycle                                                           
        do 21 i2=i1+1,np 
	  if(VANISH(I2).EQ.1)	cycle                                                       
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        t2(i2)=xexp(-dr2/al)                                                     
	if(ifdelt.ge.1) t2new(i2)=xexp(-dr2/alnew)
        idr22=dr2*25.+1.                                                        
        idr2=min0(idr22,100000)
        nccc=ih(idr2)                                                           
        u(i2)=iso(i1)*iso(i2)*cy(nccc)+cyu(nccc)
c       upau(i2)=cyu(nccc)*xexp(-p2**2/width)                                    
21      continue                                                                
        do 22 i2=i1+1,np 
	  if(VANISH(I2).EQ.1)	cycle                                                       
        t3=(x(i1)-x(i2))*u(i2)                                                  
        t4=(y(i1)-y(i2))*u(i2)                                                  
        t5=(z(i1)-z(i2))*u(i2)                                                  
	sx=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)
	sy=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)
	sz=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)

	if(ifdelt.ge.1) then
	wx=(x(i1)-x(i2))*t2new(i2)
	wy=(y(i1)-y(i2))*t2new(i2)
	wz=(z(i1)-z(i2))*t2new(i2)
	frtwo(i1)=frtwo(i1)+t2new(i2)
	frtwo(i2)=frtwo(i2)+t2new(i2)
	else
        wx=(x(i1)-x(i2))*t2(i2)                                                 
        wy=(y(i1)-y(i2))*t2(i2)                                                 
        wz=(z(i1)-z(i2))*t2(i2)                                                 
	endif
        fr(i1)=fr(i1)+t2(i2)                                                    
        fr(i2)=fr(i2)+t2(i2)                                                    
        f1x(i1)=f1x(i1)+wx                                                      
        f1x(i2)=f1x(i2)-wx                                                      
        f1y(i1)=f1y(i1)+wy                                                      
        f1y(i2)=f1y(i2)-wy                                                      
        f1z(i1)=f1z(i1)+wz                                                      
        f1z(i2)=f1z(i2)-wz                                                      
        fcx(i1)=fcx(i1)+t3                                                      
        fcx(i2)=fcx(i2)-t3                                                      
        fcy(i1)=fcy(i1)+t4                                                      
        fcy(i2)=fcy(i2)-t4                                                      
        fcz(i1)=fcz(i1)+t5                                                      
        fcz(i2)=fcz(i2)-t5                                                      
	fsx(i1)=fsx(i1)+sx
	fsx(i2)=fsx(i2)-sx
	fsy(i1)=fsy(i1)+sy
	fsy(i2)=fsy(i2)-sy
	fsz(i1)=fsz(i1)+sz
	fsz(i2)=fsz(i2)-sz
22      continue                                                                
2       continue                                                                
c                                                                               
        do 258 i=1,np 
	  if(VANISH(I).EQ.1)	cycle                                                          
        f3x=0.                                                                  
        f3y=0.                                                                  
        f3z=0.                                                                  
        do 11 i1=1,np 
	  if(VANISH(I1).EQ.1)	cycle                                                          
        dr2=(x(i1)-x(i))*(x(i1)-x(i))+(y(i1)-y(i))*(y(i1)-y(i))                 
     +  +(z(i1)-z(i))*(z(i1)-z(i))                                              
        if(i1.eq.i) goto 11                                                     
        f3x=f3x+(x(i)-x(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3y=f3y+(y(i)-y(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3z=f3z+(z(i)-z(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
11      continue                                                                
        f3x=potd*pote*f3x                                                       
        f3y=potd*pote*f3y                                                       
        f3z=potd*pote*f3z                                                       
	if(ifdelt.ge.1) then
	st=potdnew*(1.-potenew*(frtwo(i)+1.)**potc)
	else
        st=potd*(1.-pote*(fr(i)+1)**potc)                                       
	endif
        fx=f1x(i)*st+fcx(i)-f3x+fsx(i)*pots
        fy=f1y(i)*st+fcy(i)-f3y+fsy(i)*pots
        fz=f1z(i)*st+fcz(i)-f3z+fsz(i)*pots
        px(i)=px(i)+dt*fx                                                       
        py(i)=py(i)+dt*fy                                                       
        pz(i)=pz(i)+dt*fz                                                       
258     continue                                                                
        return                                                                  
        end 
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE updtv
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	  common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &        potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/bb/np,nzp,ntag,x00  
	  common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                                  
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/skyrme/al,sal,t0,t1,arang,vy,width                          
        common/coul/cy(5000),cyu(5000),cypion(5000)                                          
        common/veloc/velx(500),vely(500),velz(500)                              
        common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
        dimension f1x(500),f1y(500),f1z(500),t2(500)                            
        dimension fcx(500),fcy(500),fcz(500),u(500)                             
        dimension fvx(500),fvy(500),fvz(500),v(500)                             
        dimension dvdpx(500),dvdpy(500),dvdpz(500)                              
        dimension fhx(500),fhy(500),fhz(500),upau(500)                          
        dimension fgx(500),fgy(500),fgz(500),gra(500)                           
        dimension fsx(500),fsy(500),fsz(500)
cc
	   wigth=2.*width
	   fnorm=(3.141596*wigth)**1.5*2.   
	   potgr=2./(fnorm*width)**2./0.165
        factor=4*0.00157*500.                                                   
        vpaul=v01                                                               
        hb=0.197327                                                    
        widthp=p01                                                              
        widthq=q01                                                              
        wp2=2.*widthp**2                                                        
        wq2=2.*widthq**2                                                        
        anorm=vpaul*(hb/widthp/widthq)**3                                       
        do 1821 i1=1,np                                                         
        dvdpx(i1)=0.                                                            
        dvdpy(i1)=0.                                                            
        dvdpz(i1)=0.                                                            
1821    continue                                                                
        do 1822 i1=1,np                                                         
        do 1823 i2=i1+1,np                                                      
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     +  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        vdens=XEXP(-dr2/al)/sal                                                  
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr                                         
        dvdp1=-anorm/widthp**2*XEXP(-dr2/wq2)*XEXP(-p2/wp2)                       
        dvdpx(i1)=dvdpx(i1)+(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i1)=dvdpy(i1)+(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i1)=dvdpz(i1)+(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
        dvdpx(i2)=dvdpx(i2)-(px(i1)-px(i2))*(dvdp+dvdp1)                        
        dvdpy(i2)=dvdpy(i2)-(py(i1)-py(i2))*(dvdp+dvdp1)                        
        dvdpz(i2)=dvdpz(i2)-(pz(i1)-pz(i2))*(dvdp+dvdp1)                        
1823    continue                                                                
1822    continue                                                                
c                                                                               
c                                                                               
        do 1 i=1,np                                                              
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
c                                                                               
        velx(i) = px(i)/e1 + dvdpx(i)                                           
        vely(i) = py(i)/e1 + dvdpy(i)                                           
        velz(i) = pz(i)/e1 + dvdpz(i)                                           
c                                                                               
c                                                                               
        x(i) = x(i) + velx(i)*dt                                                
        y(i) = y(i) + vely(i)*dt                                                
        z(i) = z(i) + velz(i)*dt                                                
c                                                                               
c ende                                                                          
c                                                                               
        f1x(i)=0.                                                               
        f1y(i)=0.                                                               
        f1z(i)=0.                                                               
        fvx(i)=0.                                                               
        fvy(i)=0.                                                               
        fvz(i)=0.                                                               
        fcx(i)=0.                                                               
        fcy(i)=0.                                                               
        fcz(i)=0.                                                               
        fr(i)=0.                                                                
        fhx(i)=0.                                                               
        fhy(i)=0.                                                               
        fhz(i)=0.                                                               
        fgx(i)=0.                                                               
        fgy(i)=0.                                                               
        fgz(i)=0.                                                               
	fsx(i)=0.
	fsy(i)=0.
	fsz(i)=0.
1	 continue                                                                              
        do 2 i1=1,np                                                            
        do 21 i2=i1+1,np                                                        
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     +  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        t2(i2)=XEXP(-dr2/al)                                                 
        v(i2)=t2(i2)*(alog(p2*500.+1.))**2.                                     
        idr22=dr2*40.+1.5                                                       
        idr2=min0(idr22,100000)                                                  
        nccc=ih(idr2)                                                           
        u(i2)=iso(i1)*iso(i2)*cy(nccc)+cyu(nccc)
        upau(i2)=anorm/widthq**2*XEXP(-dr2/wq2)*XEXP(-p2/wp2)
        gra(i2)=2.5-dr2/al                                                                                             
21      continue                                                                
        do 22 i2=i1+1,np                                                        
        t3=(x(i1)-x(i2))*u(i2)                                                  
        t4=(y(i1)-y(i2))*u(i2)                                                  
        t5=(z(i1)-z(i2))*u(i2)                                                  
	sx=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)
	sy=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)
	sz=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)
        v3=(x(i1)-x(i2))*v(i2)                                                  
        v4=(y(i1)-y(i2))*v(i2)                                                  
        v5=(z(i1)-z(i2))*v(i2)                                                  
        wx=(x(i1)-x(i2))*t2(i2)                                                 
        wy=(y(i1)-y(i2))*t2(i2)                                                 
        wz=(z(i1)-z(i2))*t2(i2)                                                 
ccc
        grx=wx*gra(i2)                                                          
        gry=wy*gra(i2)                                                          
        grz=wz*gra(i2)                                                          
ccc
        hx=(x(i1)-x(i2))*upau(i2)                                               
        hy=(y(i1)-y(i2))*upau(i2)                                               
        hz=(z(i1)-z(i2))*upau(i2)   
cc                           
ccc  
	  fr(i1)=fr(i1)+t2(i2)                                                    
        fr(i2)=fr(i2)+t2(i2)                                                    
        f1x(i1)=f1x(i1)+wx                                                      
        f1x(i2)=f1x(i2)-wx                                                      
        f1y(i1)=f1y(i1)+wy                                                      
        f1y(i2)=f1y(i2)-wy                                                      
        f1z(i1)=f1z(i1)+wz                                                      
        f1z(i2)=f1z(i2)-wz                                                      
        fcx(i1)=fcx(i1)+t3                                                      
        fcx(i2)=fcx(i2)-t3                                                      
        fcy(i1)=fcy(i1)+t4                                                      
        fcy(i2)=fcy(i2)-t4                                                      
        fcz(i1)=fcz(i1)+t5                                                      
        fcz(i2)=fcz(i2)-t5                                                      
        fvx(i1)=fvx(i1)+v3                                                      
        fvx(i2)=fvx(i2)-v3                                                      
        fvy(i1)=fvy(i1)+v4                                                      
        fvy(i2)=fvy(i2)-v4                                                      
        fvz(i1)=fvz(i1)+v5                                                      
        fvz(i2)=fvz(i2)-v5                                                      
        fhx(i1)=fhx(i1)+hx                                                      
        fhx(i2)=fhx(i2)-hx                                                      
        fhy(i1)=fhy(i1)+hy                                                      
        fhy(i2)=fhy(i2)-hy                                                      
        fhz(i1)=fhz(i1)+hz                                                      
        fhz(i2)=fhz(i2)-hz                                                      
        fgx(i1)=fgx(i1)+grx                                                     
        fgx(i2)=fgx(i2)-grx                                                     
        fgy(i1)=fgy(i1)+gry                                                     
        fgy(i2)=fgy(i2)-gry                                                     
        fgz(i1)=fgz(i1)+grz                                                     
        fgz(i2)=fgz(i2)-grz                                                     
	fsx(i1)=fsx(i1)+sx
	fsx(i2)=fsx(i2)-sx
	fsy(i1)=fsy(i1)+sy
	fsy(i2)=fsy(i2)-sy
	fsz(i1)=fsz(i1)+sz
	fsz(i2)=fsz(i2)-sz
22      continue                                                        
2       continue                                                                
c                                                                               
        do 258 i=1,np                                                           
        f3x=0.                                                                  
        f3y=0.                                                                  
        f3z=0.                                                                  
        do 11 i1=1,np                                                           
        dr2=(x(i1)-x(i))*(x(i1)-x(i))+(y(i1)-y(i))*(y(i1)-y(i))                 
     &  +(z(i1)-z(i))*(z(i1)-z(i))                                              
        if(i1.eq.i) goto 11                                                     
        f3x=f3x+(x(i)-x(i1))*XEXP(-dr2/al)*(fr(i1)+1)**potc                      
        f3y=f3y+(y(i)-y(i1))*XEXP(-dr2/al)*(fr(i1)+1)**potc                      
        f3z=f3z+(z(i)-z(i1))*XEXP(-dr2/al)*(fr(i1)+1)**potc                      
11      continue                                                                
        f3x=potd*pote*f3x                                                       
        f3y=potd*pote*f3y                                                       
        f3z=potd*pote*f3z                                                       
        st=potd*(1.-pote*(fr(i)+1)**potc)                                       
        fx=f1x(i)*st+fcx(i)+fhx(i)+fvx(i)*potf+
     &		fgx(i)*potg-f3x+fsx(i)*pots
        fy=f1y(i)*st+fcy(i)+fhy(i)+fvy(i)*potf+
     &		fgy(i)*potg-f3y+fsy(i)*pots
        fz=f1z(i)*st+fcz(i)+fhz(i)+fvz(i)*potf+
     &		fgz(i)*potg-f3z+fsz(i)*pots
        px(i)=px(i)+dt*fx                                                       
        py(i)=py(i)+dt*fy                                                       
        pz(i)=pz(i)+dt*fz                                         
 258    continue
        return
        end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE updtim  
c*****  mean-field potentials for nucleon, resonances (delta, N*), hyperons (lambda, sigma) 
c*****  and anti-nucleon, anti-delta, anti-hyperons
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &    rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),pyms(500),
     &       pzms(500),emson(500),idms(500),icp(500),idanti(500)
	common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
	common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/bb/np,nzp,ntag,x00  
      common/densi1/denstl(500),dens1(500),velcm(3)
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
	common/xx/alpha,beta,gamma,symm                          
      common/coul/cy(5000),cyu(5000),cypion(5000)                                          
      common/veloc/velx(500),vely(500),velz(500)
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc,numit
      common/CLWF3/RTIME0
	common/wavel/sl(500)
      common/medium/mdel,mdin,nsym,npop,nkaon

      dimension f1x(500),f1y(500),f1z(500),t2(500),t2new(500) 
      dimension fcx(500),fcy(500),fcz(500),f1xc(500),f1yc(500),
     &          f1zc(500),u(500),frp(500)           
	dimension fvx(500),fvy(500),fvz(500),vm(500)                        
      dimension fsx(500),fsy(500),fsz(500),fgrx(500),fgry(500),fgrz(500)
     &      ,fsx1(500),fsy1(500),fsz1(500),fsx2(500),fsy2(500),fsz2(500)
	dimension fsgx(500),fsgy(500),fsgz(500)
	dimension dvdpx(500),dvdpy(500),dvdpz(500),ppp(500)
	dimension fhx(500),fhy(500),fhz(500),upau(500)
	common/SURF/ VANISH(500)

     
	pi=3.141592654
	pi3=pi**3
	hbar=0.197327
	Pf1=hbar*(1.5*pi**2*rho0)**0.3333
	xld= 1.5*Pf1
	xld2=xld**2

        widthp=p01 
        widthq=q01 
        wp2=2.*widthp**2 
        wq2=2.*widthq**2  

        do 1821 i=1,np
		if(VANISH(I).EQ.1)	cycle  
        dvdpx(i)=0.
        dvdpy(i)=0.
	  dvdpz(i)=0.
1821    continue

        call densitn

        do 1822 i1=1,np
		if(VANISH(I1).EQ.1)	cycle
        e1=sqrt(em(i1)*em(i1)+px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
	  pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
	  pf=hbar*(3.*pi**2*dens1(i1))**0.333333
	  if(pf.le.0.) then
	  qrf=0.
	  else 
	  qrf=pr/pf
	  endif

        do 1823 i2=i1+1,np
		if(VANISH(I2).EQ.1)	cycle
	  if((idanti(i1).GT.0) .OR. (idanti(i2).GT.0)) cycle
	  if((idr(i1).GT.3) .OR. (idr(i2).GT.3)) cycle

	  al=2.*(sl(i1)+sl(i2))
	  sal=(pi*al)**1.5*rho0
	  xmm=sqrt(sl(i1)+sl(i2))/hbar
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &  (z(i1)-z(i2))*(z(i1)-z(i2))
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))
        vdens=xexp(-dr2/al)/sal
	  
	  if(momp.eq.0) then                 !without mom. dep.
	  dvdp=0.
	  elseif(momp.eq.1) then             !with mom. by Aichelin
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif
        pbr= p2*cmeps+1.
        dvdp=4.*cfm*cmeps*vdens*alog(pbr)/pbr
	  else                               !with mom. and mass splitting
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif

!!	 Gaussian integration
	  cfmom=-2.*cfm*xmm*vdens/pi3
	  clmd=xld2*xmm**2
	  xk1=xmm*px(i1)
	  yk1=xmm*py(i1)
	  zk1=xmm*pz(i1)
	  xk2=xmm*px(i2)
	  yk2=xmm*py(i2)
	  zk2=xmm*pz(i2)
cc	  call xmoms(clmd,xk1,yk1,zk1,xk2,yk2,zk2,sumx,sumy,sumz)    !Gaussian integration
	  call xpart(1,np,clmd,xk1,yk1,zk1,xk2,yk2,zk2,sumx,sumy,sumz)
	  dvdpx(i1)=dvdpx(i1)+cfmom*sumx
	  dvdpy(i1)=dvdpy(i1)+cfmom*sumy
	  dvdpz(i1)=dvdpz(i1)+cfmom*sumz
cc	  call xmoms(clmd,xk2,yk2,zk2,xk1,yk1,zk1,sumx,sumy,sumz)
	  call xpart(1,np,clmd,xk2,yk2,zk2,xk1,yk1,zk1,sumx,sumy,sumz)
	  dvdpx(i2)=dvdpx(i2)+cfmom*sumx
	  dvdpy(i2)=dvdpy(i2)+cfmom*sumy
	  dvdpz(i2)=dvdpz(i2)+cfmom*sumz
 	  goto 1823
	  endif

	  if(iso(i1).eq.iso(i2)) then         !Pauli potential
        dvdp1=-vpb*vdens/wp2*xexp(-p2/wp2)
	  dvdp2= vpc*e01*vdens/(pf*(1.+xexp(e01*(1.-qrf**2))))**2
	  else
	  dvdp1=0.
	  dvdp2=0.
	  endif

        dvdpx(i1)=dvdpx(i1)+(px(i1)-px(i2))*(dvdp+dvdp1)+px(i1)*dvdp2                        
        dvdpy(i1)=dvdpy(i1)+(py(i1)-py(i2))*(dvdp+dvdp1)+py(i1)*dvdp2                        
        dvdpz(i1)=dvdpz(i1)+(pz(i1)-pz(i2))*(dvdp+dvdp1)+pz(i1)*dvdp2                        
        dvdpx(i2)=dvdpx(i2)-(px(i1)-px(i2))*(dvdp+dvdp1)+px(i2)*dvdp2                        
        dvdpy(i2)=dvdpy(i2)-(py(i1)-py(i2))*(dvdp+dvdp1)+py(i2)*dvdp2                        
        dvdpz(i2)=dvdpz(i2)-(pz(i1)-pz(i2))*(dvdp+dvdp1)+pz(i2)*dvdp2                       
1823    continue
	  if((idanti(i1).LE.0) .AND. (idr(i1).LE.3)) then
        velx(i1) = px(i1)/e1 + dvdpx(i1)
        vely(i1) = py(i1)/e1 + dvdpy(i1)
        velz(i1) = pz(i1)/e1 + dvdpz(i1)
	  else
	   if(nkaon.EQ.0) then
	   esg= 0.
	   else
	  call selfener(idanti(i1),idr(i1),denstl(i1),esg,eom,dvsg,dvom)
	  if((em(i1)+esg) .LE. 0.) esg= -em(i1)
	   endif

	  e1=sqrt((em(i1)+esg)**2 + px(i1)**2 + py(i1)**2 + pz(i1)**2)
        velx(i1) = px(i1)/e1
        vely(i1) = py(i1)/e1
        velz(i1) = pz(i1)/e1
	  endif
1822    continue

        do 1 i=1,np
		if(VANISH(I).EQ.1)	cycle
        x(i) = x(i) + velx(i)*dt
        y(i) = y(i) + vely(i)*dt
        z(i) = z(i) + velz(i)*dt
c
      f1x(i)=0.                                                               
      f1y(i)=0.                                                                
      f1z(i)=0.                                                               
      fcx(i)=0.                                                               
      fcy(i)=0.                                                               
      fcz(i)=0.                                                               
	fsx(i)=0.
	fsy(i)=0.
	fsz(i)=0.
	fsx1(i)=0.
	fsy1(i)=0.
	fsz1(i)=0.
	fsx2(i)=0.
	fsy2(i)=0.
	fsz2(i)=0.
	fsgx(i)=0.
	fsgy(i)=0.
	fsgz(i)=0.
      fr(i)=0.    
      fvx(i)=0.                                                               
      fvy(i)=0.                                                               
      fvz(i)=0.	                                                                   
      fhx(i)=0.                                                               
      fhy(i)=0.                                                               
      fhz(i)=0.                                                               
	fgrx(i)=0. 
	fgry(i)=0.
	fgrz(i)=0.
	f1xc(i)=0.
	f1yc(i)=0.
	f1zc(i)=0.
	frp(i)=0.
 1	continue

        do 2 i1=1,np
		if(VANISH(I1).EQ.1)	cycle
	  pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
	  pf=hbar*(3.*pi**2*dens1(i1))**0.333333
	  if(pf.le.0.) then
	  qrf=0.
	  else 
	  qrf=pr/pf
	  endif	

        do 21 i2=i1+1,np
		if(VANISH(I2).EQ.1)	cycle
	  al=2.*(sl(i1)+sl(i2))
	  sal=(pi*al)**1.5*rho0
	  xmm=sqrt(sl(i1)+sl(i2))/hbar
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))
     &  +(z(i1)-z(i2))*(z(i1)-z(i2))
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))
	  t2(i2)=XEXP(-dr2/al)/sal
        idr22=dr2*40.+1. 
        idr2=min0(idr22,100000)
        nccc=ih(idr2)
        u(i2)=iso(i1)*iso(i2)*cy(nccc)

cc	 mom. dep. interaction
	  if(momp.eq.0) then
	  vm(i2)=0.
	  elseif(momp.eq.1) then
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif
        vm(i2)=cfm*t2(i2)*(alog(p2*cmeps+1.))**2
	  else
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif

!!	 Gaussian integration
	  clmd=xld2*xmm**2
	  xk1=xmm*px(i1)
	  yk1=xmm*py(i1)
	  zk1=xmm*pz(i1)
	  xk2=xmm*px(i2)
	  yk2=xmm*py(i2)
	  zk2=xmm*pz(i2)
cc	  vm(i2)=cfm*t2(i2)/pi3*sumom(clmd,xk1,yk1,zk1,xk2,yk2,zk2)    !Gaussian integration
	  vm(i2)=cfm*t2(i2)/pi3*xmsum(1,np,clmd,xk1,yk1,zk1,xk2,yk2,zk2)
	  endif
cc	 Pauli potential
	  if(iso(i1).eq.iso(i2)) then
        upau(i2)=vpa*(1./al+1./wq2)*xexp(-dr2/wq2)+vpb/al*xexp(-p2/wp2)+
     &           vpc/al/(1.+xexp(e01*(1.-qrf**2)))
	  else
        upau(i2)=0.
	  endif	
 21    continue

        do 22 i2=i1+1,np
		if(VANISH(I2).EQ.1)	cycle
	  al=2.*(sl(i1)+sl(i2))
	  dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &  (z(i1)-z(i2))*(z(i1)-z(i2))

        t3=(x(i1)-x(i2))*u(i2)
        t4=(y(i1)-y(i2))*u(i2) 
        t5=(z(i1)-z(i2))*u(i2)

        fcx(i1)=fcx(i1)+t3 
        fcx(i2)=fcx(i2)-t3 
        fcy(i1)=fcy(i1)+t4 
        fcy(i2)=fcy(i2)-t4 
        fcz(i1)=fcz(i1)+t5 
        fcz(i2)=fcz(i2)-t5 

	  if((idanti(i1).GT.0) .OR. (idanti(i2).GT.0)) cycle     !exclude nuclear mean-field potential
	  if((idr(i1).GT.3) .OR. (idr(i2).GT.3)) cycle

	  v3=2.*(x(i1)-x(i2))*vm(i2)/al
        v4=2.*(y(i1)-y(i2))*vm(i2)/al
        v5=2.*(z(i1)-z(i2))*vm(i2)/al

	  hx=(x(i1)-x(i2))*upau(i2)*t2(i2)
        hy=(y(i1)-y(i2))*upau(i2)*t2(i2)
        hz=(z(i1)-z(i2))*upau(i2)*t2(i2)

 	  grdi=(2.5-dr2/al)/al*16.        !Skyrme  surface

	  sx1=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)/al
	  sy1=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)/al
	  sz1=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)/al
	  sx2=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)**2/al
	  sy2=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)**2/al
	  sz2=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)**2/al

	  sgdx=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)/al*grdi  !surface symmetry
	  sgdy=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)/al*grdi
	  sgdz=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)/al*grdi

	  fgradx=(x(i1)-x(i2))*t2(i2)/al*grdi
	  fgrady=(y(i1)-y(i2))*t2(i2)/al*grdi
	  fgradz=(z(i1)-z(i2))*t2(i2)/al*grdi

        wx=(x(i1)-x(i2))*t2(i2)/al                                                 
        wy=(y(i1)-y(i2))*t2(i2)/al                                                 
        wz=(z(i1)-z(i2))*t2(i2)/al                                                 

        fr(i1)=fr(i1)+t2(i2)                                                    
        fr(i2)=fr(i2)+t2(i2)
	  frp(i1)=frp(i1)+iso(i1)*iso(i2)*t2(i2)*rho0
	  frp(i2)=frp(i2)+iso(i1)*iso(i2)*t2(i2)*rho0
		                                                      
        f1x(i1)=f1x(i1)+wx                                                      
        f1x(i2)=f1x(i2)-wx                                                      
        f1y(i1)=f1y(i1)+wy                                                      
        f1y(i2)=f1y(i2)-wy                                                      
        f1z(i1)=f1z(i1)+wz                                                      
        f1z(i2)=f1z(i2)-wz 
	  fvx(i1)=fvx(i1)+v3                                                      
        fvx(i2)=fvx(i2)-v3                                                      
        fvy(i1)=fvy(i1)+v4                                                      
        fvy(i2)=fvy(i2)-v4                                                      
        fvz(i1)=fvz(i1)+v5
        fvz(i2)=fvz(i2)-v5
	  fhx(i1)=fhx(i1)+hx                                                      
        fhx(i2)=fhx(i2)-hx                                                      
        fhy(i1)=fhy(i1)+hy                                                      
        fhy(i2)=fhy(i2)-hy                                                      
        fhz(i1)=fhz(i1)+hz                                                      
        fhz(i2)=fhz(i2)-hz 
	  f1xc(i1)=f1xc(i1)+iso(i1)*iso(i2)*wx   !Coulomb exchange                                                   
        f1xc(i2)=f1xc(i2)-iso(i1)*iso(i2)*wx                                                      
        f1yc(i1)=f1yc(i1)+iso(i1)*iso(i2)*wy                                                      
        f1yc(i2)=f1yc(i2)-iso(i1)*iso(i2)*wy                                                      
        f1zc(i1)=f1zc(i1)+iso(i1)*iso(i2)*wz                                                      
        f1zc(i2)=f1zc(i2)-iso(i1)*iso(i2)*wz
	  fsx1(i1)=fsx1(i1)+sx1
	  fsx1(i2)=fsx1(i2)-sx1
	  fsy1(i1)=fsy1(i1)+sy1
	  fsy1(i2)=fsy1(i2)-sy1
	  fsz1(i1)=fsz1(i1)+sz1
	  fsz1(i2)=fsz1(i2)-sz1
	  fsx2(i1)=fsx2(i1)+sx2
	  fsx2(i2)=fsx2(i2)-sx2
	  fsy2(i1)=fsy2(i1)+sy2
	  fsy2(i2)=fsy2(i2)-sy2
	  fsz2(i1)=fsz2(i1)+sz2
	  fsz2(i2)=fsz2(i2)-sz2
	  fsgx(i1)=fsgx(i1)+sgdx
	  fsgx(i2)=fsgx(i2)-sgdx
	  fsgy(i1)=fsgy(i1)+sgdy
	  fsgy(i2)=fsgy(i2)-sgdy
	  fsgz(i1)=fsgz(i1)+sgdz
	  fsgz(i2)=fsgz(i2)-sgdz
	  fgrx(i1)=fgrx(i1)+fgradx
        fgrx(i2)=fgrx(i2)-fgradx
	  fgry(i1)=fgry(i1)+fgrady
        fgry(i2)=fgry(i2)-fgrady
	  fgrz(i1)=fgrz(i1)+fgradz
        fgrz(i2)=fgrz(i2)-fgradz
22      continue 
2       continue 
c         
        potc=gamma-1.
	  emass=eta-1.

        do 258 i=1,np

		if(VANISH(I).EQ.1)	cycle
	  fx=0.
	  fy=0.
	  fz=0. 
	  f3x=0.                                                                  
        f3y=0.                                                                  
        f3z=0.
	  ftx=0.                                                                  
        fty=0.                                                                  
        ftz=0.
	  f3xc=0.                                                                  
        f3yc=0.                                                                  
        f3zc=0.
	  fsx13=0.
	  fsy13=0.
        fsz13=0.
	  fsx23=0.
	  fsy23=0.
        fsz23=0.
	  fsx33=0.
	  fsy33=0.
        fsz33=0.
	  fsx43=0.
	  fsy43=0.
        fsz43=0.

        do 11 i1=1,np
		if(VANISH(I1).EQ.1)	cycle
	  if(idanti(i1).GT.0)  cycle     !exclude nuclear mean-field potential
	  if(idr(i1).GT.3)  cycle

	  al=2.*(sl(i1)+sl(i))
	  sal=(pi*al)**1.5*rho0    
        dr2=(x(i1)-x(i))*(x(i1)-x(i))+(y(i1)-y(i))*(y(i1)-y(i))
     &  +(z(i1)-z(i))*(z(i1)-z(i)) 
        if(i1.eq.i) goto 11
        f3x=f3x+(x(i)-x(i1))*XEXP(-dr2/al)/sal/al*fr(i1)**potc                       
        f3y=f3y+(y(i)-y(i1))*XEXP(-dr2/al)/sal/al*fr(i1)**potc                      
        f3z=f3z+(z(i)-z(i1))*XEXP(-dr2/al)/sal/al*fr(i1)**potc 
	  ftx=ftx+(x(i)-x(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass          !effective mass                   
        fty=fty+(y(i)-y(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass                      
        ftz=ftz+(z(i)-z(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass

c	  f3xc=f3xc+(x(i)-x(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al       !Coulomb exchange
c     &  /(frp(i1)+0.0001)**0.6666
c	  f3yc=f3yc+(y(i)-y(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al
c     &  /(frp(i1)+0.0001)**0.6666
c	  f3zc=f3zc+(z(i)-z(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al
c     &  /(frp(i1)+0.0001)**0.6666

	  if(nsym.eq.2) then
	  fsx13=fsx13+(x(i)-x(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(xmu-1.)
	  fsy13=fsy13+(y(i)-y(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(xmu-1.)
	  fsz13=fsz13+(z(i)-z(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(xmu-1.)
	  fsx23=fsx23+(x(i)-x(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(xmu-2.)
	  fsy23=fsy23+(y(i)-y(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(xmu-2.)
	  fsz23=fsz23+(z(i)-z(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(xmu-2.)
 	  elseif(nsym.eq.3) then
	  fsx13=fsx13+(x(i)-x(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamsy-1.)
	  fsy13=fsy13+(y(i)-y(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamsy-1.)
	  fsz13=fsz13+(z(i)-z(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamsy-1.)
	  fsx23=fsx23+(x(i)-x(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamsy-2.)
	  fsy23=fsy23+(y(i)-y(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamsy-2.)
	  fsz23=fsz23+(z(i)-z(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamsy-2.)
	  if(csym.eq.0.) goto 11
	  fsx33=fsx33+(x(i)-x(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(eta-1.)
	  fsy33=fsy33+(y(i)-y(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(eta-1.)
	  fsz33=fsz33+(z(i)-z(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(eta-1.)
	  fsx43=fsx43+(x(i)-x(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(eta-2.)
	  fsy43=fsy43+(y(i)-y(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(eta-2.)
	  fsz43=fsz43+(z(i)-z(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(eta-2.)
	  endif
11      continue

        f3x=2.*beta*gamma/(gamma+1.)*f3x                                                       
        f3y=2.*beta*gamma/(gamma+1.)*f3y                                                       
        f3z=2.*beta*gamma/(gamma+1.)*f3z  
	  ftx=2.*gtau*eta*ftx
	  fty=2.*gtau*eta*fty
	  ftz=2.*gtau*eta*ftz

 	  f1xc(i)= 0.  !Cex/3.*2.*f1xc(i)/(frp(i)+0.001)**0.6666*rho0
	  f1yc(i)= 0.  !Cex/3.*2.*f1yc(i)/(frp(i)+0.001)**0.6666*rho0
	  f1zc(i)= 0.  !Cex/3.*2.*f1zc(i)/(frp(i)+0.001)**0.6666*rho0
	  f3xc=    0.  !Cex/3.*2.*f3xc*rho0
	  f3yc=    0.  !Cex/3.*2.*f3yc*rho0
	  f3zc=    0.  !Cex/3.*2.*f3zc*rho0

	  if(nsym.eq.2) then
	  fs3x=symm*(fsx13+(xmu-1.)*fsx23)
	  fs3y=symm*(fsy13+(xmu-1.)*fsy23)
	  fs3z=symm*(fsz13+(xmu-1.)*fsz23)
	  elseif(nsym.eq.3) then
	  fs3x=2.*bsym*(fsx13+(gamsy-1.)*fsx23)+
     &       2.*csym*(fsx33+(eta-1.)*fsx43)
	  fs3y=2.*bsym*(fsy13+(gamsy-1.)*fsy23)+
     &       2.*csym*(fsy33+(eta-1.)*fsy43)
	  fs3z=2.*bsym*(fsz13+(gamsy-1.)*fsz23)+
     &       2.*csym*(fsz33+(eta-1.)*fsz43)
	  else
	  fs3x=0.
	  fs3y=0.
	  fs3z=0.
	  endif
      st=2.*alpha+2.*beta*gamma/(gamma+1.)*fr(i)**potc+
     &   2.*gtau*eta*fr(i)**emass

	  if(nsym.eq.1) then
	  fsx(i)=2.*symm*fsx1(i)
	  fsy(i)=2.*symm*fsy1(i)
	  fsz(i)=2.*symm*fsz1(i)
	  elseif(nsym.eq.2) then
	  fsx(i)=symm*fsx1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsx2(i)*(fr(i)+0.0001)**(xmu-2.)
	  fsy(i)=symm*fsy1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsy2(i)*(fr(i)+0.0001)**(xmu-2.)
	  fsz(i)=symm*fsz1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsz2(i)*(fr(i)+0.0001)**(xmu-2.)
 	  elseif(nsym.eq.3) then
	  fsx(i)=4.*asym*fsx1(i)+2.*fsx1(i)*(bsym*(fr(i)+0.0001)**
     &       (gamsy-1.)+csym*(fr(i)+0.0001)**(eta-1.))+
     &       2.*fsx2(i)*(bsym*(gamsy-1.)*(fr(i)+0.0001)**(gamsy-2.)+
     &       csym*(eta-1.)*(fr(i)+0.0001)**(eta-2.))
	  fsy(i)=4.*asym*fsy1(i)+2.*fsy1(i)*(bsym*(fr(i)+0.0001)**
     &       (gamsy-1.)+csym*(fr(i)+0.0001)**(eta-1.))+
     &       2.*fsy2(i)*(bsym*(gamsy-1.)*(fr(i)+0.0001)**(gamsy-2.)+
     &       csym*(eta-1.)*(fr(i)+0.0001)**(eta-2.))
	  fsz(i)=4.*asym*fsz1(i)+2.*fsz1(i)*(bsym*(fr(i)+0.0001)**
     &       (gamsy-1.)+csym*(fr(i)+0.0001)**(eta-1.))+
     &       2.*fsz2(i)*(bsym*(gamsy-1.)*(fr(i)+0.0001)**(gamsy-2.)+
     &       csym*(eta-1.)*(fr(i)+0.0001)**(eta-2.))
	  else
	  fsx(i)=0.
	  fsy(i)=0.
	  fsz(i)=0.
	  endif

	  if((idr(i).LE.3) .AND. (idanti(i).LE.0)) then
        fx=f1x(i)*st+fcx(i)+f3x+fsx(i)+0.5*gsur*fgrx(i)+fs3x
        fy=f1y(i)*st+fcy(i)+f3y+fsy(i)+0.5*gsur*fgry(i)+fs3y
        fz=f1z(i)*st+fcz(i)+f3z+fsz(i)+0.5*gsur*fgrz(i)+fs3z
	  fx=fx+f1xc(i)+f3xc+ftx+fvx(i)+fhx(i)+0.5*gsym*fsgx(i)
	  fy=fy+f1yc(i)+f3yc+fty+fvy(i)+fhy(i)+0.5*gsym*fsgy(i)
	  fz=fz+f1zc(i)+f3zc+ftz+fvz(i)+fhz(i)+0.5*gsym*fsgz(i)
	  else
	   if(nkaon.EQ.0) then
	   dvf= 0.
	   else
	  call selfener(idanti(i),idr(i),denstl(i),esg,eom,dvsg,dvom)
	  if((em(i)+esg) .LE. 0.) esg= -em(i)
	  e1=sqrt((em(i)+esg)**2 + px(i)**2 + py(i)**2 + pz(i)**2)
	  dvf= dvsg*(em(i)+esg)/e1 + dvom
	   endif

	  fx=fcx(i) + 2.*rho0*dvf*f1x(i)
	  fy=fcy(i) + 2.*rho0*dvf*f1y(i)
	  fz=fcz(i) + 2.*rho0*dvf*f1z(i)
	  endif

        px(i)=px(i)+dt*fx
        py(i)=py(i)+dt*fy
        pz(i)=pz(i)+dt*fz
 258    continue


        return
        end 
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE udtmsn
c***** Coulomb potential and in-medium dispersion of mesons
c///////////////////////////////////////////////////////////////////////////////
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
	  common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &  potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
        common/bb/np,nzp,ntag,x00
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/ee/mseed,nseed
	  common/coul/cy(5000),cyu(5000),cypion(5000)
	  common/wavel/sl(500)   
        common/medium/mdel,mdin,nsym,npop,nkaon
	common/SURF/ VANISH(500)

	 pi   = 3.1415926
	 hb   = 0.197327
	 AVMS = 0.938868
       PIMS = 0.137265
	 DLMS = 1.232
	 CUTF = 7.*PIMS
	 CPND = 2.
	 GPR  = 0.6             !Migdal parameter
	 gsp  = 2.
	 ss   = 0.

	 do 2 i1=1,np
		
	 if(idms(i1)  .eq. 0)      goto 2
	 if(emson(i1) .lt. 0.01)   goto 2

	  fcmx= 0.
	  fcmy= 0.
	  fcmz= 0.
	  fmdx= 0.
	  fmdy= 0.
	  fmdz= 0.

	  fopx= 0.
	  fopy= 0.
	  fopz= 0.
	  fisx= 0.
	  fisy= 0.
	  fisz= 0.

	  velx=0.
	  vely=0.
	  velz=0.

	  xp=xms(i1)
	  yp=yms(i1)
	  zp=zms(i1)
	  rho = densit(xp,yp,zp)
	  rhon= denstn(xp,yp,zp)
	  rhop= denstp(xp,yp,zp)
	  pabs=pxms(i1)**2+pyms(i1)**2+pzms(i1)**2
	  pab=sqrt(pabs)

	  DRH =  rho+rho/10.
        dens= denscal(DRH)
	  dvsc= 10.*(dens-denscal(rho))/rho     ! value[0,1]
	  rtis= rho**(gsp-1.)/rho0**gsp

	 if(npop .le. 0) then
	 apn  = 0.
	 ciso = 0.
	 elseif(npop .eq. 1) then
	 apn  = 0.154
	 ciso = 0.036 
	 else
	 apn= 0.
	 ciso = 0.036 
	 endif

	  if(idms(i1).eq.1) then

	 if(npop .LT. 2)  then
        xx=xexp(-apn*rho/rho0)
	  xm0=emson(i1)*(1.+6.5*(1.-xx**10))
	  xp0=sqrt((1.-xx)**2*emson(i1)**2+2.*xm0*emson(i1)*(1.-xx))
	  if(pab.eq.0.) goto 2

	  cm0=65.*emson(i1)*apn*xx**10
	  if(xp0.eq.0.) then
	  cp0=0.
	  else
	  cp0=apn*xx*((1.-xx)*emson(i1)**2+xm0*emson(i1))/xp0
	  endif

	  xp=(pab-xp0)/pab
	  ep=sqrt(xm0**2+(pab-xp0)**2)

	  cmd=xm0/ep-xm0/sqrt(xp0**2+xm0**2)
	  cpd=(xp0-pab)/ep-xp0/sqrt(xp0**2+xm0**2)
	  velx = xp*pxms(i1)/ep
	  vely = xp*pyms(i1)/ep
	  velz = xp*pzms(i1)/ep
	  goto  222

	 else
C----------------------- EVALUATING PION ENERGY FROM DELTA-HOLE MODEL -----
	  omdl = sqrt(DLMS**2+pabs) - AVMS
	  ckaf = 0.888889*GPR*hb**3*(CPND/PIMS)**2*omdl*rho*
     &         exp(-2.*pabs/CUTF**2)                                   !GeV^2
	  aopm = pabs*ckaf*(1.-GPR)/GPR + (omdl**2-PIMS**2)*ckaf +
     &         0.25*(ckaf+PIMS**2+pabs-omdl**2)**2                     !GeV^4
	  if(aopm .LT. 0.)  aopm = 0.
	  cpp = 0.5*(ckaf+PIMS**2+pabs-omdl**2)
	  ydl = cpp + sqrt(aopm) + omdl**2
	  ypn = cpp - sqrt(aopm) + omdl**2
	  if(ydl .LT. 0.)  ydl=0.
	  if(ypn .LT. 0.)  ypn=0.
	  ompnlk = sqrt(ypn)
	  omdllk = sqrt(ydl)

	  if(ompnlk .EQ. omdl)  then
	  ss=1.
	  else
	  cksi    = -0.888889*hb**3*(CPND/PIMS)**2*rho*omdl/
     &             (omdl**2-ompnlk**2)*exp(-2.*pabs/CUTF**2)
	  ss=1.+pabs*ckaf/GPR/(1.-GPR*cksi)**2/(omdl**2-ompnlk**2)**2
	  ss=1./ss
	  endif

	  dadp = (1.-2.*pabs/CUTF**2)*ckaf*(1.-GPR)/GPR - 
     &          2.*(omdl**2-PIMS**2)*ckaf/CUTF**2 +
     &         (1.-2.*ckaf/CUTF**2)*(ckaf+PIMS**2+pabs-omdl**2)/2.
	  if(aopm .LE. 0.) then
	  dypn = (1.-2.*ckaf/CUTF**2)
	  dydl = (1.-2.*ckaf/CUTF**2)
	  else
	  dypn = (1.-2.*ckaf/CUTF**2) - dadp/sqrt(aopm)
	  dydl = (1.-2.*ckaf/CUTF**2) + dadp/sqrt(aopm)
	  endif
	  dompn = (omdl/sqrt(DLMS**2+pabs) + 0.5*dypn)/ompnlk
	  domdl = (omdl/sqrt(DLMS**2+pabs) + 0.5*dydl)/omdllk

	  domdp = ss*dompn + (1.-ss)*domdl

	  velx  = pxms(i1)*domdp
	  vely  = pyms(i1)*domdp
	  velz  = pzms(i1)*domdp
	  goto  222
	  endif

	  elseif(idms(i1).eq.2) then
	   if(nkaon.EQ.0) then
	   ak= 0.
	   bk= 0.
	   else
	   ak= 0.249         !GeV^2fm^3  SigmaEtaN= 280 MeV, Veta= -94. MeV (rho0)
	   bk= 1.8           !fm^3
	   endif
	  ometa=(emson(i1)**2-ak*denscal(rho))/(1.+bk*denscal(rho)) + pabs
	  if(ometa.LE.0.)  goto 2
	  ep=sqrt(ometa)

	  elseif(idms(i1).eq.3) then
	  is= 2*icp(i1)-1
	  if(nkaon.EQ.0) then
	  ak= 0.
	  bk= 0.
	  else
	  ak= 0.18         !GeV^2fm^3,  SigmaKN= 202.62MeV, Vk= 27.8 MeV (rho0)
	  bk= 0.333        !GeVfm^3
	  ck= 0.0298       !isospin splitting
	  dk= 0.111
	  endif
	  omkn=emson(i1)**2+pabs-ak*denscal(rho)-is*ck*(denscal(rhon)-
     &	   denscal(rhop))+(bk*rho+is*dk*(rhon-rhop))**2
	  if(omkn.LE.0.)  goto 2
	  ep=sqrt(omkn)

	  elseif(idms(i1).eq.4) then
	  is= 2*icp(i1)+1
	  if(nkaon.EQ.0) then
	  ak= 0.
	  bk= 0.
	  else
	  ak= 0.31         !SigmaKN=348.9 MeV, Vanti-K= -100.3 MeV
	  bk= 0.333
	  ck= 0.0298
	  dk= 0.111
	  endif
	  omak=emson(i1)**2+pabs-ak*denscal(rho)-is*ck*(denscal(rhon)-
     &	   denscal(rhop))+(bk*rho+is*dk*(rhon-rhop))**2
	  if(omak.LE.0.)  goto 2
	  ep=sqrt(omak)

	  else
	  ep=sqrt(emson(i1)**2+pabs)
	  endif

	  velx=pxms(i1)/ep
	  vely=pyms(i1)/ep
	  velz=pzms(i1)/ep

 222	  xms(i1) = xms(i1) + dt*velx
	  yms(i1) = yms(i1) + dt*vely
	  zms(i1) = zms(i1) + dt*velz
c***** evolutions of mesons





        do 21 i2=1,np
		if(VANISH(I2).EQ.1)	cycle
c*****  interaction of mesons and baryons
	  al = 2.*sl(i2)
	  sal= (pi*al)**1.5
        dr2= (xms(i1)-x(i2))*(xms(i1)-x(i2))+(yms(i1)-y(i2))*
     &       (yms(i1)-y(i2))+(zms(i1)-z(i2))*(zms(i1)-z(i2))
	  roij=XEXP(-dr2/al)/sal
        idr22=dr2*40.+1.
        idr2= min0(idr22,100000)
        nccc= ih(idr2)
        upion=icp(i1)*iso(i2)*cypion(nccc)

        fcmx= fcmx + (xms(i1)-x(i2))*upion
        fcmy= fcmy + (yms(i1)-y(i2))*upion
        fcmz= fcmz + (zms(i1)-z(i2))*upion
	  fopx= fopx + (xms(i1)-x(i2))*roij/al
	  fopy= fopy + (yms(i1)-y(i2))*roij/al
	  fopz= fopz + (zms(i1)-z(i2))*roij/al
	  fisx= fisx + isy(i2)*(xms(i1)-x(i2))*roij/al
	  fisy= fisy + isy(i2)*(yms(i1)-y(i2))*roij/al
	  fisz= fisz + isy(i2)*(zms(i1)-z(i2))*roij/al

 21	  continue

	 if(idms(i1).eq.1) then
	   if(npop .LT. 2) then
	 fmdx= 2.*(cm0*cmd+cp0*cpd)*fopx/rho0 - 2.*ciso*icp(i1)*fisx*rtis
	 fmdy= 2.*(cm0*cmd+cp0*cpd)*fopy/rho0 - 2.*ciso*icp(i1)*fisy*rtis
	 fmdz= 2.*(cm0*cmd+cp0*cpd)*fopz/rho0 - 2.*ciso*icp(i1)*fisz*rtis
	   else
	  omdl = sqrt(DLMS**2+pabs) - AVMS
	  ckaf = 0.888889*GPR*hb**3*(CPND/PIMS)**2*omdl*
     &         exp(-2.*pabs/CUTF**2)                                         !GeV^2
	  aopm = pabs*ckaf*rho*(1.-GPR)/GPR + (omdl**2-PIMS**2)*ckaf*rho +
     &         0.25*(ckaf*rho+PIMS**2+pabs-omdl**2)**2                       !GeV^4
	  if(aopm .LT. 0.)  aopm = 0.
	  cpp = 0.5*(ckaf*rho+PIMS**2+pabs-omdl**2)
	  ydl = cpp + sqrt(aopm) + omdl**2
	  ypn = cpp - sqrt(aopm) + omdl**2
	  if(ydl .LT. 0.)  ydl=0.
	  if(ypn .LT. 0.)  ypn=0.
	  ompnlk = sqrt(ypn)
	  omdllk = sqrt(ydl)

	  if(ompnlk .EQ. omdl)  then
	  ss=1.
	  else
	  cksi    = -0.888889*hb**3*(CPND/PIMS)**2*rho*omdl/
     &             (omdl**2-ompnlk**2)*exp(-2.*pabs/CUTF**2)
	  ss=1.+pabs*ckaf*rho/GPR/(1.-GPR*cksi)**2/(omdl**2-ompnlk**2)**2
	  ss=1./ss
	  endif

	  dadp = 0.5*pabs*ckaf*(1.-GPR)/GPR + 0.5*(omdl**2-PIMS**2)*ckaf +
     &         0.25*ckaf*(ckaf*rho+PIMS**2+pabs-omdl**2)
	  if(aopm .LE. 0.) then
	  dypn = 0.5*ckaf
	  dydl = 0.5*ckaf
	  else
	  dypn = 0.5*ckaf - dadp/sqrt(aopm)
	  dydl = 0.5*ckaf + dadp/sqrt(aopm)
	  endif

	  domdr = ss*dypn/ompnlk + (1.-ss)*dydl/omdllk

	 fmdx= domdr*fopx - 2.*ciso*icp(i1)*fisx*rtis
	 fmdy= domdr*fopy - 2.*ciso*icp(i1)*fisy*rtis
	 fmdz= domdr*fopz - 2.*ciso*icp(i1)*fisz*rtis

	   endif

	 elseif(idms(i1).eq.2) then
	 fmdx= -ak*dvsc*fopx/(1.+bk*denscal(rho))/ep - bk*dvsc*fopx*
     &        (emson(i1)**2-ak*denscal(rho))/(1.+bk*denscal(rho))**2/ep
	 fmdy= -ak*dvsc*fopy/(1.+bk*denscal(rho))/ep - bk*dvsc*fopy*
     &        (emson(i1)**2-ak*denscal(rho))/(1.+bk*denscal(rho))**2/ep
	 fmdz= -ak*dvsc*fopz/(1.+bk*denscal(rho))/ep - bk*dvsc*fopz*
     &        (emson(i1)**2-ak*denscal(rho))/(1.+bk*denscal(rho))**2/ep

	 elseif(idms(i1).eq.3) then
! derivative on baryon density
	 fmdx= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopx+is*dk*fisx)/ep+
     &       2.*bk*fopx + 2.*is*dk*fisx
	 fmdy= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopy+is*dk*fisy)/ep+
     &       2.*bk*fopy + 2.*is*dk*fisy
	 fmdz= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopz+is*dk*fisz)/ep+
     &       2.*bk*fopz + 2.*is*dk*fisz
! derivative on scalar density
	 fmdx= fmdx - (ak*dvsc*fopx + is*ck*dvsc*fisx)/ep
	 fmdy= fmdy - (ak*dvsc*fopy + is*ck*dvsc*fisy)/ep
	 fmdz= fmdz - (ak*dvsc*fopz + is*ck*dvsc*fisz)/ep

	   if(nkaon .EQ. 1) then            !Lorentz-force 
	 fmdx= fmdx - 2.*bk*velx**2*fopx
	 fmdy= fmdy - 2.*bk*vely**2*fopy
	 fmdz= fmdz - 2.*bk*velz**2*fopz
	   endif

	 elseif(idms(i1).eq.4) then
! derivative on baryon density
	 fmdx= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopx+is*dk*fisx)/ep-
     &       2.*bk*fopx - 2.*is*dk*fisx
	 fmdy= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopy+is*dk*fisy)/ep-
     &       2.*bk*fopy - 2.*is*dk*fisy
	 fmdz= 2.*(bk*rho+is*dk*(rhon-rhop))*(bk*fopz+is*dk*fisz)/ep-
     &       2.*bk*fopz - 2.*is*dk*fisz
! derivative on scalar density
	 fmdx= fmdx - (ak*dvsc*fopx+is*ck*dvsc*fisx)/ep
	 fmdy= fmdy - (ak*dvsc*fopy+is*ck*dvsc*fisy)/ep
	 fmdz= fmdz - (ak*dvsc*fopz+is*ck*dvsc*fisz)/ep

	   if(nkaon .EQ. 1) then
	 fmdx= fmdx + 2.*bk*velx**2*fopx
	 fmdy= fmdy + 2.*bk*vely**2*fopy
	 fmdz= fmdz + 2.*bk*velz**2*fopz
	   endif

	 else
	 fmdx= 0.
	 fmdy= 0.
	 fmdz= 0.
	 endif

	 fx=fcmx + fmdx
	 fy=fcmy + fmdy
	 fz=fcmz + fmdz

	 pxms(i1) = pxms(i1)+dt*fx
	 pyms(i1) = pyms(i1)+dt*fy
	 pzms(i1) = pzms(i1)+dt*fz

 2	 continue

        return 
        end 


c///////////////////////////////////////////////////////////////////////////////         
	  function  denscal(rho)            !Scalar density as function of baryon density
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
	 eps=1.E-8
	 hb= 0.197327
	 pi= 3.1415926
	 pf= hb*(1.5*pi**2*rho)**0.33333
	 ex= 0.938/(1.+0.52*rho/rho0)
	 cs= pi**2*hb**3/2.
	 call valgss(0.,pf,eps,ex,vds)
 	 denscal= vds/cs
        return
        end 
c///////////////////////////////////////////////////////////////////////////////         
	SUBROUTINE valgss(A,B,EPS,EFF,G)
	real C(5),T(5)
	real A,B,F,G,S,P,H,AA,BB,W,X,Q
	DATA T/-0.90618,-0.53847,0.,0.53847,0.90618/
	DATA C/0.23693,0.47863,0.56889,0.47863,0.23693/

	M=1
	S=(B-A)*0.01
	P=0.
	H=(B-A)/M
	G=0.
	 do l=1,M
	AA=A+(l-1)*H
	BB=A+l*H
	W=0.
	  do j=1,5
           x=((BB-AA)*T(j)+(BB+AA))/2.
           W=W+xfun(x,EFF)*C(j)
	  enddo
	G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))

 111	P=G
	M=M+1
	H=(B-A)/M
      G=0.
	 do l=1,M
      AA=A+(l-1)*H
      BB=A+l*H
      W=0.
	  do j=1,5
	     x=((BB-AA)*T(j)+(BB+AA))/2.
	     W=W+xfun(x,EFF)*C(j)
	  enddo
	  G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))
	if((Q.GE.EPS).AND.(ABS(H).GT.ABS(S))) goto 111

 	 return
	 end

       function xfun(xx,eff)
	 xfun=xx**2*eff/sqrt(xx**2+eff**2)
       return
       end 
c///////////////////////////////////////////////////////////////////////////////
	  SUBROUTINE selfener(iant,iren,rho,esig,eomg,dvsg,dvom)

	 XMSCA= 0.508194     !GeV
	 XMVEC= 0.782501
	 XMRHO= 0.763

	 GSN= 80.82
	 GVN= 155.
	 GRN= 20.017
	 XAT= 0.25

	 esig= 0.
	 eomg= 0.
	 dvsg= 0.
	 dvom= 0.

	   if(iren .LE. 3) then
	   GS2= GSN
	   GV2= GVN
	   GR2= GRN
	   elseif(iren .GT. 3 .and. iren .LE. 5) then   !lamda and sigma
	   GS2= 0.66667*GSN
	   GV2= 0.66667*GVN
	   elseif(iren .GT. 5 .and. iren .LE. 6) then   !Xi
	   GS2= 0.33333*GSN
	   GV2= 0.33333*GVN
	   else
	   return
	   endif

         if((iant.GT.0) .and. (iren.LE.3)) then   !anti-particles
	   GS2= XAT*GS2
	   GV2= XAT*GV2
	   GR2= XAT*GR2
	   endif

	 DRH= rho+rho/10.

	 call  GFINT(denscal(rho), XMSCA, GS2, PSCA)
	 call  GFINT(rho, XMVEC, GV2, PVEC)

	 call  GFINT(denscal(DRH), XMSCA, GS2, PSCAFN)
	 call  GFINT(DRH, XMVEC, GV2, PVECFN)

	 esig=-PSCA
	 dvsg= 10.*(PSCA-PSCAFN)/rho
	 if(iant .LE. 0) then
	 eomg= PVEC
	 dvom= 10.*(PVECFN-PVEC)/rho
	 else
	 eomg=-PVEC
	 dvom= 10.*(PVEC-PVECFN)/rho
	 endif

	 return
	 end

      SUBROUTINE GFINT(RHOI,XMASS,CUPCON,VPOT)
C
C This SUBROUTINE does the Green function integration
C
      REAL MESH,NINE,MASSFM
      REAL ARRAY(500),ITILD(500),HALF(500),I0
      REAL EMASS(2000)
C
      DATA ZERO,ONE,TWO,FOUR,FIVE,SIX,NINE,FIFTN,SIXTN/0.0D0,1.0D0,
     &     2.0D0,4.0D0,5.0D0,6.0D0,9.0D0,15.0D0,16.0D0/
C
C
	 MESH= 0.05
	 NPT= 400
       MASSFM= XMASS/0.197327
	  
	 do I=1,NPT
	 ARRAY(I)=RHOI*FLOAT(NPT-I+1)/FLOAT(NPT)
	 enddo

C
C Set up exponential array
C
       MPT=2*NPT
       do 5 I=1,MPT
         X=(I-1)*MESH/TWO
         EMASS(I)=XEXP(MASSFM*X)
 5     continue
C
C Interpolate (four point here)
C
       HALF(2)=(FIVE*ARRAY(1)+FIFTN*ARRAY(2)-
     &         FIVE*ARRAY(3)+ARRAY(4))/SIXTN
       HALF(NPT)=(ARRAY(NPT-3)-FIVE*ARRAY(NPT-2)+
     &           FIFTN*ARRAY(NPT-1)+FIVE*ARRAY(NPT))/SIXTN
C
      NSTOP=NPT-1
       do 50 I=3,NSTOP
         HALF(I)=(NINE*(ARRAY(I)+ARRAY(I-1))-
     &           (ARRAY(I+1)+ARRAY(I-2)))/SIXTN
 50    continue
C
C Calculate ITILD integral
C
         ITILD(NPT)=ZERO
         do 200 I=1,NSTOP
            J=NPT-I+1
            X=(J-1)*MESH
            XTEMP=X-MESH/TWO
C
            TEMP=XTEMP/EMASS(2*J-2)*HALF(J)*FOUR
     &           + (X-MESH)/EMASS(2*J-3)*ARRAY(J-1)
     &           + X/EMASS(2*J-1)*ARRAY(J)
C
            ITILD(J-1)=ITILD(J)+TEMP*MESH/SIX
 200     continue

         I0= ITILD(1)
         VPOT= I0*CUPCON*HBARC

 	  return
	  end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE coulg1
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/coul/cy(5000),cyu(5000),cypion(5000)                                          
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &        potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/gg/nzta,nzpr                                                                                            
c                                                                               
        a1=al/4.                                                                
        gam=1./arang                                                            
        ct0r=0.005                                                              
        ceff=float(nzpr+nzta)/float(npr+nta)                                    
        c1=0.5*XEXP(a1*gam*gam)                                                  
        c2=1./sqrt(a1*3.1415927)                                                
        c3=0.5*vy                                                               
        c4=sqrt(al)                                                             
        c5=0.5*0.00144*ceff*ceff                                                
c                                                                               
        t0=t0-ct0r/gam/gam/gam*vy                                               
c                                                                               
        cy(1)=0.                                                                
        cy(5000)=0.                                                             
        do 1 j=2,4999                                                           
         dr=(float(j)-1.)/100.                                                  
         dr2=dr*dr                                                              
         t1=dr/c4                                                               
         t2=XEXP(-t1*t1)                                                         
         call erf(t1,t3)                                                        
         t4=t3/dr-1.128/c4*t2                                                   
         t5=c5/dr2*t4                                                           
         t11=2.*a1*gam-dr                                                       
         t12=2.*a1*gam+dr                                                       
         if(t11.lt.0.) then                                                     
         call erf(-t11/c4,t13)                                                  
         t13=-t13                                                               
         else                                                                   
         call erf(t11/c4,t13)                                                   
         endif                                                                 
         call erf(t12/c4,t14)                                                   
         t15=1.-t13                                                             
         t16=1.-t14                                                             
         t18=XEXP(-gam*dr)                                                       
         t19=XEXP(gam*dr)                                                        
         t20=t18*t15-t19*t16                                                    
         t21=t18*(gam*t15-c2*XEXP(-t11*t11/a1))                                  
         t22=t19*(gam*t16-c2*XEXP(-t12*t12/a1))                                  
         t23=c3*c1/dr2/gam*(1./dr*t20+t21+t22)                                  
         cy(j)=t5+t23                                                           
c        fr=XEXP(-dr2/al)                                                        
c        st=pota/sal/al*(1.-potb*(fr/sal)**potc)*dr                             
c        t551=t5*dr                                                             
c        t231=t23*dr                                                            
c        t999=t321+t551+st                                                                               
1       continue                                                                
        return                                                                  
        end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE erf(x,er) 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        p=.3275911                                                              
        a1=.254829592                                                           
        a2=-.284496736                                                          
        a3=1.421413741                                                          
        a4=-1.453152027                                                         
        a5=1.061405429                                                          
        t=1./(1+p*x)                                                            
        er=1-(a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5)*XEXP(-x**2)                  
        return                                                                  
        end                                                                     
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE kaon(s,em1,em2,i1,i2)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common/kaot/pktot,nkn                                                 
        data amk,aml,amn/ 0.495,1.115,0.938 /                                 
        pmax=(s-(amk+aml+amn)**2)*(s-(aml+amn-amk)**2)/4./s                   
        pk=0.0                                                                
        if (pmax.lt.0.0) return                                               
        pmax=sqrt(pmax)/amk/40.                                               
        if (em1.lt.1..and.em2.lt.1.) pk=0.072*pmax                            
        if (em1.lt.1..and.em2.gt.1.) pk=0.054*pmax                            
        if (em1.gt.1..and.em2.lt.1.) pk=0.054*pmax                            
        if (em1.gt.1..and.em2.gt.1.) pk=0.036*pmax                            
        pktot=pktot+pk                                                        
        nkn=nkn+1                                                             
	  I1=I1
	  I2=I2
	  return                                                                
        end                                                                   
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE minspa    
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        dimension i(500,500),it(500),itt(500),ifrag(500)                        
        do 8015 k=1,500                                                         
        it(k)=0                                                                 
        itt(k)=0                                                                
        ifrag(k)=0                                                              
        ipy=0                                                                   
        do 8016 kk=1,500                                                        
        i(k,kk)=0                                                               
8016    continue                                                                
8015    continue                                                                
        do 8017 m=1,np                                                          
        mmm=0                                                                   
        do12 mm=1,np                                                            
        if(m.eq.mm)goto12                                                       
        dr=sqrt((x(m)-x(mm))**2+(y(m)-y(mm))**2+(z(m)-z(mm))**2)                
        if(dr.lt.3.)mmm=mmm+1                                                   
        if(dr.lt.3.)i(m,mmm)=mm                                                 
        if(dr.lt.3.)it(m)=mmm                                                   
12      continue                                                                
        if(it(m).eq.0)ipy=ipy+1                                                 
124     format(1x,22i5)                                                         
8017    continue                                                                
        do 1 k=1,np                                                             
        if(it(k).eq.0)goto 1                                                    
        m1=0                                                                    
4       m1=m1+1                                                                 
        if(i(k,m1).eq.0)goto11                                                  
        do 2 mm=1,it(i(k,m1))                                                   
        ic=0                                                                    
        do 3 m=1,it(k)                                                          
        if(i(i(k,m1),mm).eq.i(k,m))ic=1                                         
        if(i(i(k,m1),mm).eq.k)ic=1                                              
3       continue                                                                
        if(ic.eq.0)it(k)=it(k)+1                                                
        if(ic.eq.0)i(k,it(k))=i(i(k,m1),mm)                                     
        i(i(k,m1),mm)=0                                                         
        it(i(k,m1))=0                                                           
2       continue                                                                
        goto 4                                                                  
11      itt(k)=it(k)+1                                                          
        ifrag(itt(k))=ifrag(itt(k))+1                                           
1       continue                                                                
        do 8018 k=1,np                                                          
        iev(k)=1                                                                
8018    continue                                                                
        do 51 k=1,np                                                            
        ifram=3                                                                 
c       ifram= smallest cluster considered as fragment                          
        if(itt(k).lt.ifram) goto 51                                             
        iev(k)=0                                                                
        do 8019 l=1,it(k)                                                       
        iev(i(k,l))=0                                                           
8019    continue                                                                
51      continue                                                                
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE clallnew 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/deu3/dg(50,10,2),dtg(2),dtg1(2),
     &  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
	do 3 l=1,50
	do 3 ll=1,10
	do 3 ideu=1,2
	dcm(l,ll,ideu)=0.
3	dg(l,ll,ideu)=0.
	do 4 ideu=1,2
	dtgt(ideu)=0.
4	dtg1t(ideu)=0.
ccc 
        return                                                                  
        end                                                                     
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE clall 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
         common/deu3/dg(50,10,2),dtg(2),dtg1(2),
     &  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
        common/bb/np,nzp,ntag,x00                                                   
        do 1 k=1,500
        x(k)=0.                                                                 
        y(k)=0.                                                                 
        z(k)=0.                                                                 
        px(k)=0.                                                                
        py(k)=0.                                                                
        pz(k)=0.                                                                
        em(k)=0.                                                                
        iso(k)=0                                                                
	isy(k)=0
        id(k)=0                                                                 
        il(k)=0                                                                 
        fr(k)=0.                                                                
1       continue                                                                
        do 2 k=np+1,500 
        iev(k)=0                                                                
2       continue                                                                
        return                                                                  
        end                                                                     

	function ranf(m,n)
	m=973*m+1
	m=m-m/1594323*1594323
	n=769*n+1
	n=n-n/2097152*2097152
	ranf=6.27225e-7*m+4.76837e-7*n
	if(ranf.gt.1.0) ranf=ranf-1.0
	return
	end

      FUNCTION RANFL(IDUM)
      implicit real (m)
      PARAMETER (MBIG=4000000.00,MSEED=1618033.00,MZ=0.,FAC=2.5e-7)
!      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      real MA(55)
      common/randumb/inext,inextp,ma     ! needed on HP
      DATA IFF /0/

       RANFL=ran(IDUM)
	 return    !windows random number

      if(IDUM.LT.0.OR.IFF.EQ.0)then
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        do 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          if(MK.LT.MZ) MK=MK+MBIG
          MJ=MA(II)
11      continue
        do 13 K=1,4
          do 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            if(MA(I).LT.MZ) MA(I)=MA(I)+MBIG
12        continue
13      continue
        INEXT=0
        INEXTP=31
        IDUM=1
      endif
      INEXT=INEXT+1
      if(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      if(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      if(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RANFL=MJ*FAC

	return
      end
c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE dtest(i1,i2,em1,em2)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
        common/deu2/px1,py1,pz1,px2,py2,pz2    
	  common/SURF/ VANISH(500)
                                 
        common/deu3/dg(50,10,2),dtg(2),dtg1(2),                                 
     1  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
        fiso=3./8.                                                              
        fa=8.*3./4.                                                             
        do 2 k=1,np 
		if(VANISH(K).EQ.1)	cycle
        im1=1                                                                   
        im2=1                                                                   
        im1o=1                                                                  
        im2o=1                                                                  
        if(k.eq.i1)goto 2                                                       
        if(k.eq.i2)goto 2                                                       
        xr_proj=((x(i1)-x(k))**2+(y(i1)-y(k))**2+(z(i1)-z(k))**2)                   
        xr_targ=((x(i2)-x(k))**2+(y(i2)-y(k))**2+(z(i2)-z(k))**2)                   
        if(xr_proj.gt.70)im1=0                                                      
        if(xr_targ.gt.70)im2=0                                                      
        if(xr_proj.gt.70)im1o=0                                                     
        if(xr_targ.gt.70)im2o=0                                                     
        pxr1=((px(i1)-px(k))**2+(py(i1)-py(k))**2+(pz(i1)-pz(k))**2)            
        pxr2=((px(i2)-px(k))**2+(py(i2)-py(k))**2+(pz(i2)-pz(k))**2)            
        if(pxr1.gt.2)im1=0                                                      
        if(pxr2.gt.2)im2=0                                                      
        pxs1=((px(i1)+px(k))**2+(py(i1)+py(k))**2+(pz(i1)+pz(k))**2)            
        pxs2=((px(i2)+px(k))**2+(py(i2)+py(k))**2+(pz(i2)+pz(k))**2)            
        pxr1o=((px1-px(k))**2+(py1-py(k))**2+(pz1-pz(k))**2)                    
        pxr2o=((px2-px(k))**2+(py2-py(k))**2+(pz2-pz(k))**2)                    
        if(pxr1o.gt.2)im1o=0                                                    
        if(pxr2o.gt.2)im2o=0                                                    
        pxs1o=((px1+px(k))**2+(py1+py(k))**2+(pz1+pz(k))**2)                    
        pxs2o=((px2+px(k))**2+(py2+py(k))**2+(pz2+pz(k))**2)                    
        prob1=fd1*fiso*XEXP(-xr_proj*fdr*im1)*XEXP(-pxr1*fdp*im1)*im1                 
        prob2=fd1*fiso*XEXP(-xr_targ*fdr*im2)*XEXP(-pxr2*fdp*im2)*im2                 
        prob1o=fd1*fiso*XEXP(-xr_proj*fdr*im1)*XEXP(-pxr1o*fdp*im1o)*
     &  im1o             
        prob2o=fd1*fiso*XEXP(-xr_targ*fdr*im2)*XEXP(-pxr2o*fdp*im2o)*
     &  im2o
        iecm1=int((sqrt(pxs1+3.519)-1.876)/0.05)+1                              
        iecm2=int((sqrt(pxs2+3.519)-1.876)/0.05)+1                              
        iecm1o=int((sqrt(pxs1o+3.519)-1.876)/0.05)+1                            
        iecm2o=int((sqrt(pxs2o+3.519)-1.876)/0.05)+1                            
        acm1=(pz(i1)+pz(k))/sqrt(pxs1)                                          
        acm2=(pz(i2)+pz(k))/sqrt(pxs2)                                          
        acm1o=(pz1+pz(k))/sqrt(pxs1o)                                           
        acm2o=(pz2+pz(k))/sqrt(pxs2o)                                           
        if(acm1.gt.0.9999)acm1=0.9999                                           
        if(acm2.gt.0.9999)acm2=0.9999                                           
        if(acm1o.gt.0.9999)acm1o=0.9999                                         
        if(acm2o.gt.0.9999)acm2o=0.9999                                         
        acm1=acos(acm1)*180./3.1415                                             
        acm2=acos(acm2)*180./3.1415                                             
        acm1o=acos(acm1o)*180./3.1415                                           
        acm2o=acos(acm2o)*180./3.1415                                           
        icm1=10                                                                 
        icm2=10                                                                 
        icm1o=10                                                                
        icm2o=10                                                                
        icm1=int(acm1/20.)+1                                                    
        icm2=int(acm2/20.)+1                                                    
        icm1o=int(acm1o/20.)+1                                                  
        icm2o=int(acm2o/20.)+1                                                  
        if(iso(i1)+iso(k).ne.1) goto 212                                         
        if(iev(i1).eq.0.or.iev(k).eq.0) goto 212                                 
        if(iecm1.lt.1.or.iecm1.gt.50) goto 211                                    
        dcm(iecm1,icm1,ideu)=dcm(iecm1,icm1,ideu)+prob1                         
        dtg(ideu)=dtg(ideu)+prob1                                               
211     if(iecm1o.lt.1.or.iecm1o.gt.50) goto 212                                  
        if(id(i1).eq.1.or.id(i1).eq.51) goto 212                                  
        dcm(iecm1o,icm1o,ideu)=dcm(iecm1o,icm1o,ideu)-prob1o                    
        dtg(ideu)=dtg(ideu)-prob1o                                              
212     if(iso(i2)+iso(k).ne.1)goto 214                                         
        if(iev(i2).eq.0.or.iev(k).eq.0) goto 214                                  
        if(iecm2.lt.1.or.iecm2.gt.50) goto 213                                    
        dcm(iecm2,icm2,ideu)=dcm(iecm2,icm2,ideu)+prob2                         
        dtg(ideu)=dtg(ideu)+prob2                                               
213     if(iecm2o.lt.1.or.iecm2o.gt.50) goto 214                                  
        if(id(i2).eq.1.or.id(i2).eq.51) goto 214                                  
        dcm(iecm2o,icm2o,ideu)=dcm(iecm2o,icm2o,ideu)-prob2o                    
        dtg(ideu)=dtg(ideu)-prob2o                                              
214     continue                                                                
        t1=pz(i1)+pz(k)                                                         
        e1=sqrt(pxs1+1.876**2)                                                  
        t2=e1                                                                   
        e2=sqrt(pxs2+1.876**2)                                                  
        e1o=sqrt(pxs1o+1.876**2)                                                
        e2o=sqrt(pxs2o+1.876**2)                                                
        pb1=bzz*(pz(i1)+pz(k))                                                  
        pb2=bzz*(pz(i2)+pz(k))                                                  
        pb1o=bzz*(pz1+pz(k))                                                    
        pb2o=bzz*(pz2+pz(k))                                                    
        b2=bzz**2                                                               
        gg1=(gam3-1.)*pb1/b2+gam3*e1                                            
        gg2=(gam3-1.)*pb2/b2+gam3*e2                                            
        gg1o=(gam3-1.)*pb1o/b2+gam3*e1o                                         
        gg2o=(gam3-1.)*pb2o/b2+gam3*e2o                                         
        pl1=gg1*bzz+(pz(i1)+pz(k))                                              
        pl2=gg2*bzz+(pz(i2)+pz(k))                                              
        pl1o=gg1o*bzz+(pz1+pz(k))                                               
        pl2o=gg2o*bzz+(pz2+pz(k))                                               
        p1=sqrt(pl1**2+(px(i1)+px(k))**2+(py(i1)+py(k))**2)                     
        p2=sqrt(pl2**2+(px(i2)+px(k))**2+(py(i2)+py(k))**2)                     
        e1=sqrt(p1**2+1.876**2)-1.876                                           
        e2=sqrt(p2**2+1.876**2)-1.876                                           
        p1o=sqrt(pl1o**2+(px1+px(k))**2+(py1+py(k))**2)                         
        p2o=sqrt(pl2o**2+(px2+px(k))**2+(py2+py(k))**2)                         
        e1o=sqrt(p1o**2+1.876**2)-1.876                                         
        e2o=sqrt(p2o**2+1.876**2)-1.876                                         
        wp1=pl1/p1                                                              
        if(wp1.gt.0.9999)wp1=0.9999                                             
        wp2=pl2/p2                                                              
        if(wp2.gt.0.9999)wp2=0.9999                                             
        wp3=pl1o/p1o                                                            
        if(wp3.gt.0.9999)wp3=0.9999                                             
        wp4=pl2o/p2o                                                            
        if(wp4.gt.0.9999)wp4=0.9999                                             
        t1=acos(wp1)*180/3.1415                                                 
        t2=acos(wp2)*180/3.1415                                                 
        t1o=acos(wp3)*180/3.1415                                                
        t2o=acos(wp4)*180/3.1415                                                
        ibit1=10                                                                
        ibit2=10                                                                
        ibit1o=10                                                               
        ibit2o=10                                                               
        ibit1=int(t1/20.)+1                                                     
        ibit2=int(t2/20.)+1                                                     
        ibit1o=int(t1o/20.)+1                                                   
        ibit2o=int(t2o/20.)+1                                                   
        n1=int(p1/0.1)+1                                                        
        n2=int(p2/0.1)+1                                                        
        n1o=int(p1o/0.1)+1                                                      
        n2o=int(p2o/0.1)+1                                                      
        if(iev(i1).eq.0.or.iev(k).eq.0)goto20                                   
        if(iso(i1)+iso(k).ne.1)goto20                                           
        if(n1.lt.1.or.n1.gt.50)goto10                                           
        dg(n1,ibit1,ideu)=dg(n1,ibit1,ideu)+prob1                               
        dtg1(ideu)=dtg1(ideu)+prob1                                             
10      if(n1o.lt.1.or.n1o.gt.50)goto20                                         
        if(id(i1).eq.1.or.id(i1).eq.51)goto 20                                  
        dg(n1o,ibit1o,ideu)=dg(n1o,ibit1o,ideu)-prob1o                          
20      if(iso(i2)+iso(k).ne.1)goto 2                                           
        if(iev(i2).eq.0.or.iev(k).eq.0) goto 2                                    
        if(n2.lt.1.or.n2.gt.50)goto30                                           
        dg(n2,ibit2,ideu)=dg(n2,ibit2,ideu)+prob2                               
        dtg1(ideu)=dtg1(ideu)+prob2                                             
30      if(n2o.lt.1.or.n2o.gt.50) goto 2                                          
        if(id(i2).eq.1.or.id(i2).eq.51)goto 2                                   
        dg(n2o,ibit2o,ideu)=dg(n2o,ibit2o,ideu)-prob2o                          
2       continue                                                                
	  EM1=EM1
	  EM2=EM2
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////
         	SUBROUTINE rdensity
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
	  common/gg/nzta,nzpr
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	  common/rdens/rtd(200)
        atz=0.                                                                  
        apz=0.                                                                  
        bpz=0.                                                                  
        btz=0.                                                                  
        atx=0.                                                                  
        apx=0.                                                                  
        bpx=0.                                                                  
        btx=0.                                                                  
        apy=0.                                                                  
        aty=0.                                                                  
        bpy=0.                                                                  
        bty=0.                                                                  
        do 36 i=1,np                                                            
        em(i)=0.938                                                             
        if(id(i).eq.1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        endif 
        if(id(i).eq.51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        endif
   36   continue                                                                
        do 15 i=1,np
        if(id(i).eq.1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        endif                                                                  
        if(id(i).eq.51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
	endif
15      continue                                                                
	do 5 i=1,30
5	rtd(i)=0.
	a1=2.*width
	pi=3.141592654
	coff=1./(a1*pi)**1.5

	do 10 radius=0.,rrp+rrt+4.5,.5
	number=int((radius+.5)/.5)
	if(radius.eq.0.)then
	do 20 nkl=1,np
	r2=x(nkl)**2.+y(nkl)**2.+z(nkl)**2.
	rtd(number)=rtd(number)+coff*XEXP(-r2/a1)
20	continue
	else
	texp=0.
	do 30 ki1=1,np

	r2=x(ki1)*x(ki1)+y(ki1)*y(ki1)+z(ki1)*z(ki1)
		if(sqrt(r2).eq.0)then
		rtd(number)=rtd(number)+coff*XEXP(-radius*radius/a1)
		else
	exp12=(XEXP(sqrt(r2)*radius/width)-XEXP(-sqrt(r2)*
     &	radius/width))/sqrt(r2)
	texp=exp12*XEXP(-r2/a1)*width/2./radius*XEXP(-radius*radius/a1)
	rtd(number)=rtd(number)+coff*texp
		endif
30	continue  	
	endif
10	continue
	write(11,*)'time=',rtime
	limit=int((rrp+rrt+4.5)/.5)
	do 40 i=1,limit
	radius=float(i-1)*.5
	write(11,100)radius,rtd(i)
100	format(1x,2f9.4)
40	continue
	return
	end
c///////////////////////////////////////////////////////////////////////////////
         	SUBROUTINE rdeninteg
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
 	  common/gg/nzta,nzpr
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	  common/rint/rinteg(200)
	  common/kkk/number,ilj
	  common radius
        atz=0.                                                                  
        apz=0.                                                                  
        bpz=0.                                                                  
        btz=0.                                                                  
        atx=0.                                                                  
        apx=0.                                                                  
        bpx=0.                                                                  
        btx=0.                                                                  
        apy=0.                                                                  
        aty=0.                                                                  
        bpy=0.                                                                  
        bty=0.                                                                  
        do 36 i=1,np                                                            
        em(i)=0.938                                                             
        if(id(i).eq.1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        endif 
        if(id(i).eq.51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        endif
   36   continue                                                                
        do 15 i=1,np
        if(id(i).eq.1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        endif                                                                  
        if(id(i).eq.51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
	endif
15      continue                                                                
	do 5 i=1,30
5	rinteg(i)=0.
	a1=2.*width
	pi=3.141592654
	do 10 radius=0.,rrp+rrt+4.5,.5
	number=int((radius+.5)/.5)
	if(radius.eq.0.)then
	factor=1./(a1*pi)**1.5/0.16
	do 20 nkl=1,np
	r2=x(nkl)**2.+y(nkl)**2.+z(nkl)**2.
	rinteg(number)=rinteg(number)+factor*XEXP(-r2/a1)
20	continue
	go to 10
	endif
	eps=.0001
	factor=1./4./pi/(a1*pi)**1.5*XEXP(-radius**2./a1)	
	call simps2(0.,pi,eps,kflag,n,sum)
	rinteg(number)=rinteg(number)+sum*factor
10	continue
	write(95,*)'time=',rtime
	limit=int((rrp+rrt+4.5)/.5)
	do 30 i=1,limit
	radius=float(i-1)*.5
	write(95,100)radius,rinteg(i)
100	format(1x,2f9.4)
30	continue
	return
	end

c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE simps2(a,b,eps,kflag,n,sum)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	n2=1
	kflag=0
	c=abs(a)+abs(b)
	h2=.5*(b-a)
	call simps1(a,eps,kflag,n1,ss1)
	call simps1(b,eps,kflag,n3,ss2)
	n=n1+n3+2
	tb1=h2*(ss1+ss2)
4	xli=a-h2
	tb2=.5*tb1
	do 5 j=1,n2
	xli=xli+2.*h2
	call simps1(xli,eps,kflag,n1,s)
	n=n+n1+1
5	tb2=tb2+h2*s
	sum=.3333333333*(4.*tb2-tb1)
	n2=n2+n2
	if(n2.lt.8)go to 6
	if(abs(sum-sum0).le.eps*(abs(sum)+1.)) return
6	sum0=sum
	tb1=tb2
	h2=.5*h2
	if(c+h2.ne.c)go to 4
	n=-n
	return
	end

c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE simps1(xli,eps,kflag,n1,s)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	n1=1
	call fls(xli,y1,y2)
	h1=.5*(y2-y1)
	e=abs(y1)+abs(y2)
	t1=h1*(f(xli,y1)+f(xli,y2))
1	yli=y1-h1
	t2=.5*t1
	do 2 i=1,n1
	yli=yli+2.*h1
2	t2=t2+h1*f(xli,yli)
	s=.3333333333*(4.*t2-t1)
	n1=n1+n1
	if(n1.lt.8)go to 3
	if(abs(s-s0).le.eps*(abs(s)+1.))return
3	s0=s
	t1=t2
	h1=.5*h1
	if(e+h1.ne.e)go to 1
	if(y1.eq.y2)return
	kflag=kflag+1
	return
	end

c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE fls(xli,y1,y2)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	pi=3.141592654
	y1=0.
	y2=2.*pi
	XLI=XLI
	return
	end

c///////////////////////////////////////////////////////////////////////////////
		function f(xli,yli)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
	  common/gg/nzta,nzpr
	  common/kkk/number,ilj
	  common radius
	ccc=0.
	pi=3.141592654
	a1=2.*width
	do 40 ilj=1,np
	ri2=x(ilj)**2.+y(ilj)**2.+z(ilj)**2.
	aaa=XEXP(-ri2/2./width)
	bbb=x(ilj)*sin(xli)*cos(yli)+y(ilj)*sin(xli)*sin(yli)+
     &z(ilj)*cos(xli)
	ddd=XEXP(radius/width*bbb)
	ccc=ccc+ddd*aaa*sin(xli)
40	continue
	f=ccc/0.16                                                      
	return
	end

c///////////////////////////////////////////////////////////////////////////////
          	SUBROUTINE axyzdensity
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
	  common/gg/nzta,nzpr
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	  common/xyz/xtd(200),ytd(200),ztd(200),axyz(200)
	  real xxx
        atz=0.                                                                  
        apz=0.                                                                  
        bpz=0.                                                                  
        btz=0.                                                                  
        atx=0.                                                                  
        apx=0.                                                                  
        bpx=0.                                                                  
        btx=0.                                                                  
        apy=0.                                                                  
        aty=0.                                                                  
        bpy=0.                                                                  
        bty=0.                                                                  
        do 36 i=1,np                                                            
        if(id(i).eq.1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        endif
        if(id(i).eq.51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        endif 
36      continue                                                                
        do 15 i=1,np 
        if(id(i).eq.1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        endif                                                                  
        if(id(i).eq.51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
	endif
15      continue                                                                
	a1=al/2.
	pi=3.141592654
	factor=1./(a1*pi)**1.5
	limit=int((rrt+rrp+4.5)/0.5)
	do 5 i=1,limit
	xtd(i)=0.
	ytd(i)=0.
	ztd(i)=0.
	axyz(i)=0.
5	continue
	do 10 xxx=0.,rrp+rrt+4.5,0.5
	i=int((xxx+0.5)/0.5)
	do 20 j=1,np
	xi2=(xxx-x(j))**2.+y(j)**2.+z(j)**2.
	a=xi2/a1
	aexp=XEXP(-a)
	aexp=factor*aexp
	xtd(i)=xtd(i)+aexp
	yi2=(xxx-y(j))**2.+z(j)**2.+x(j)**2.
	b=yi2/a1
	bexp=XEXP(-b)
	bexp=factor*bexp
	ytd(i)=ytd(i)+bexp
	zi2=(xxx-z(j))**2.+y(j)**2.+x(j)**2.
	c=zi2/a1
	cexp=XEXP(-c)
	cexp=factor*cexp
	ztd(i)=ztd(i)+cexp
20	continue
10	continue
	write(21,*)'time=',rtime
	do 100 i=1,limit
	axyz(i)=(xtd(i)+ytd(i)+ztd(i))/3.
	dist=float(i-1)*0.5
	write(21,200)dist,axyz(i)
100	continue
200	format(1x,2f9.4)
	return
	end

c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE maxdensit
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,nzp,ntag,x00                                                                              
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	common/max/mden
	dimension dens1(500)
	real mden
	mden=0.  	
        witdh=al/2                                                    
        anorm=(3.141596*witdh)**1.5                                   
        do 10 i=1,np                                                  
        dens1(i)=0.0                                                  
10      continue                                                      
        do 20 i=1,np                                                  
        do 30 j=i+1,np                                                
        dr2=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2              
        dd=XEXP(-dr2/witdh)                                            
        dens1(i)=dens1(i)+dd                                          
        dens1(j)=dens1(j)+dd                                          
30      continue                                                      
        dens1(i)=(1.+dens1(i))/anorm                                  
	if(dens1(i).gt.mden)mden=dens1(i)
20      continue                                                      
	mden=mden/0.16
	write(110,55555)rtime,mden
55555	format(1x,2f9.4)
        return                                                        
        end                                                           

c///////////////////////////////////////////////////////////////////////////////
       	SUBROUTINE coulgn
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        double precision xdble                                                  
        common/coul/cy(5000),cyu(5000),cypion(5000)                                          
        common/atcy/sco(5000),syu(5000)                                       
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/bb/np,nzp,ntag,x00                                                   
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &        potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	  common/gg/nzta,nzpr
        a1=al/4.                                                                
        gam=1./arang                                                            
        c1 = .5                                                                 
        c2 = 1./sqrt(a1*3.1415927)                                                
        c3 = vy                                                               
        c4=sqrt(al)                                                             
        c5 = .00144
        cy(1)=0.
        cy(5000)=0.
        cyu(1)=0.                                                               
        cyu(5000)=0.                                                            
        cypion(1)=0.                                                               
        cypion(5000)=0.                                                            
	sco(1)=0.
	sco(5000)=0.
	syu(1)=0.
	syu(5000)=0.
        do 1 j=2,4999 
        dr=(float(j)-1.)/100.                                                  
        dr2=dr*dr                                                              
        t1=dr/c4                                                               
        t2=XEXP(-t1*t1)                                                         
        xdble = dble( t1 )                                                   
        t3 = 1. - erfc ( xdble )                                              
        t4=t3/dr-1.128/c4*t2                                                   
        t5=c5/dr2*t4 
	sco(j)=c5/dr*t3
	cy(j)=t5
	cypion(j)=c5/dr**3
        t11=2.*a1*gam-dr                                                     
        t12=2.*a1*gam+dr                                                       
        xdble = dble( t11/c4 )                                                 
        t15 = erfc( xdble )                                                   
        xdble = dble( t12/c4 )                                                 
        t16 = erfc( xdble )                                                   
        t18=XEXP(-gam*dr)                                                       
        t19=XEXP(gam*dr)                                                        
        t20= XEXP(a1*gam*gam)*(t18*t15-t19*t16)                               
        t21=t18*( gam*XEXP(a1*gam*gam)*t15 -                                    
     &      c2*XEXP(-t11*t11/ al+a1*gam*gam ))                            
        t22=t19*( gam*XEXP(a1*gam*gam)*t16 -                                    
     &      c2*XEXP(-t12*t12/ al+a1*gam*gam ))                            
        t23=c3*c1/dr2*(1./dr*t20+t21+t22)                                    
	  syu(j)=c3*c1/dr*t20
        cyu(j)=t23                                                           
1       continue

        return 
        end 
c///////////////////////////////////////////////////////////////////////////////
	function erfc (x)                                   
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      double precision s15adf                                                   
      double precision x,xdexp                                        
      integer                          ifail                                    
      double precision                 t, xhi, xlo, y                           
      intrinsic                        abs, exp                                 
      data xlo/-6.5d0/                                                          
      data xhi/9.5d0/                                                           
      ifail = 0                                                                 
      if (x.ge.xhi) go to 20                                                    
      if (x.le.xlo) go to 40                                                    
      t = 1.0d0 - 7.5d0/(abs(x)+3.75d0)                                         
      y = (((((((((((((((+3.328130055126039d-10                                 
     *    *t-5.718639670776992d-10)*t-4.066088879757269d-9)                     
     *    *t+7.532536116142436d-9)*t+3.026547320064576d-8)                      
     *    *t-7.043998994397452d-8)*t-1.822565715362025d-7)                      
     *    *t+6.575825478226343d-7)*t+7.478317101785790d-7)                      
     *    *t-6.182369348098529d-6)*t+3.584014089915968d-6)                      
     *    *t+4.789838226695987d-5)*t-1.524627476123466d-4)                      
     *    *t-2.553523453642242d-5)*t+1.802962431316418d-3)                      
     *    *t-8.220621168415435d-3)*t+2.414322397093253d-2                       
      y = (((((y*t-5.480232669380236d-2)*t+1.026043120322792d-1)                
     *    *t-1.635718955239687d-1)*t+2.260080669166197d-1)                      
     *    *t-2.734219314954260d-1)*t + 1.455897212750385d-1                     
      s15adf = XDEXP(-x*x)*y                                                      
      if (x.lt.0.0d0) s15adf = 2.0d0 - s15adf                                   
      erfc = sngl (s15adf)
      return                                                                    
  20  s15adf = 0.0d0                                                            
      erfc = sngl (s15adf)
      return                                                                    
  40  s15adf = 2.0d0                                                            
      erfc = sngl (s15adf)
      return                                                                    
      end                                                                      	

c///////////////////////////////////////////////////////////////////////////////
        	SUBROUTINE avene(avpot,aveki,Esh,idn)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,               
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	common/rpi/rmst(500),pmst(500),avet(500)
      common/bb/np,nzp,ntag,x00                                                   
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
	common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                       
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
	common/gg/nzta,nzpr
      common/densi1/denstl(500),dens1(500),velcm(3)
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/coul/cy(5000),cyu(5000),cypion(5000)                                          
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
      common/atcy/sco(5000),syu(5000)
	common/xx/alpha,beta,gamma,symm
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01                                     
	common/wavel/sl(500)  
	common/shell/ish,ds,r0p,r0t,r0c,Eshpp,Eshpn,Eshtp,Eshtn,Eshcp,
     &             Eshcn,Ri0
      common/medium/mdel,mdin,nsym,npop,nkaon
	dimension t2(500),t2new(500),frp(500),frtwop(500),
     &   frt(500),frtwot(500)
     
      nzp=nzpr+nzta

      if(idn.eq.1)then
      imin=1
      imax=npr
      nnuc=npr
      endif
      if(idn.eq.2)then
      imin=npr+1
      imax=np
      nnuc=nta
      endif
      if(idn.eq.3)then
      imin=1
      imax=np
      nnuc=np
      endif 
cc
	tkine=0.
	ucou=0.
	uyuk=0.
	utwo=0.
	uthr=0.
	upaul=0.
	umomd=0.
	usym=0.
	ussr=0.
	ugradi=0.
	Ekin=0.
	utau=0.
	couex=0.
	Eshp=0.
	Eshn=0.
	Esh=0.
cc
	pi=3.141592654
	pi3=pi**3
      hb=0.197327
	Ekfac=12.4515      !MeV
	ccf=3.*pi**2./2.
	Pf1=hb*(1.5*pi**2*rho0)**0.3333
	xld= 1.5*Pf1
	xld2=xld**2

	Ri=Ri0
	Rf=r0c-r0t
	shp=1.+1.6228*ds/r0p
	sht=1.+1.6228*ds/r0t
	shc=1.+1.6228*ds/r0c
cc
    
	if(iflag5.eq.1)then
      factor=cdelt
	else
	factor=0.
	endif

        widthp=p01                                                              
        widthq=q01                                                              
        wp2=2.*widthp**2                                                        
        wq2=2.*widthq**2                                                        

	if(ish.eq.0)then	  
	OD3=1./3.
	TD3=2./3.
	CB=0.0058
	CS=0.26

	ip=int(npr-nzpr+0.0001)
	jp=nzpr
	Eshp=CB*((FS(ip)+GS(jp))/(float(npr)/2.)**TD3-CS*float(npr)**OD3)
	Eshp=Eshp/npr*shp
	ip=int(nta-nzta+0.0001)
	jp=nzta
	Esht=CB*((FS(ip)+GS(jp))/(float(nta)/2.)**TD3-CS*float(nta)**OD3)
	Esht=Esht/nta*sht
	ip=int(np-nzpr-nzta+0.0001)
	jp=nzpr+nzta
	Eshc=CB*((FS(ip)+GS(jp))/(float(np)/2.)**TD3-CS*float(np)**OD3)
	Eshc=Eshc/np*shc
	else
	Eshpp=Eshpp/nzpr*shp
	Eshpn=Eshpn/(npr-nzpr)*shp
	Eshtp=Eshtp/nzta*sht
	Eshtn=Eshtn/(nta-nzta)*sht
	Eshcp=Eshcp/nzp*shc
	Eshcn=Eshcn/(np-nzp)*shc
	endif

	call center(apx,apy,apz,1)
	call center(atx,aty,atz,2) 
	Rpt=(apx-atx)**2+(apy-aty)**2+(apz-atz)**2
	Rpt=sqrt(Rpt)
	if(idn.ne.3) Rpt=Ri
	r01=r0p*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))
	r02=r0t*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))

	if(ish.eq.0)then
	Esh1=Eshp*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))
	Esh2=Esht*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))
	else
	Esh1p=Eshpp*SWITCH(Ri,Rf,Rpt)+Eshcp*(1.-SWITCH(Ri,Rf,Rpt))
	Esh1n=Eshpn*SWITCH(Ri,Rf,Rpt)+Eshcn*(1.-SWITCH(Ri,Rf,Rpt))	
	Esh2p=Eshtp*SWITCH(Ri,Rf,Rpt)+Eshcp*(1.-SWITCH(Ri,Rf,Rpt))
	Esh2n=Eshtn*SWITCH(Ri,Rf,Rpt)+Eshcn*(1.-SWITCH(Ri,Rf,Rpt))		
	endif

	do i=imin,imax
	fr(i)=0.
	frp(i)=0.
	enddo

      call densitn                                                           
	do 22 i1=imin,imax
	  pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        e1=sqrt(em(i1)*em(i1)+px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
	  pf=hb*(3.*pi**2*dens1(i1))**0.3333
	  if(pf.le.0.) then
	  qrf=0.
	  else 
	  qrf=pr/pf
	  endif
        do 23 i2=i1+1,imax  
        al=2.*(sl(i1)+sl(i2))
   	  sal=(pi*al)**1.5*rho0        
	  xmm=sqrt(sl(i1)+sl(i2))/hb	                                
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &  (z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        pbr=p2*cmeps+1.
        t2(i2)=xexp(-dr2/al)/sal
	                                                       
        idr22=dr2*40.+1.5                                                       
        idr2=min0(idr22,100000)
      nccc=ih(idr2)                                                          
	ucou=ucou+iso(i1)*iso(i2)*sco(nccc)
	ugradi=ugradi+gsur*(6./al-4.*dr2/al**2.)*t2(i2)
c	uyuk=uyuk+syu(nccc)   fengzhq
	if(iso(i1).eq.iso(i2))then
	upaul=upaul+vpa*xexp(-dr2/wq2)*t2(i2)+vpb*xexp(-p2/wp2)*t2(i2)+
     &      vpc*t2(i2)/(1.+xexp(e01*(1.-qrf**2)))
	endif
	  if(momp.eq.0) then
	umomd=0.
	  elseif(momp.eq.1) then
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif
	umomd=umomd+cfm*alog(pbr)*alog(pbr)*t2(i2)/2.
	  else
	   if(iso(i1).eq.iso(i2)) then
	   cfm=cmnn
	   else
	   cfm=cmnp
	   endif
	  clmd=xld2*xmm**2
	  xk1=xmm*px(i1)
	  yk1=xmm*py(i1)
	  zk1=xmm*pz(i1)
	  xk2=xmm*px(i2)
	  yk2=xmm*py(i2)
	  zk2=xmm*pz(i2)
cc	  xss=sumom(clmd,xk1,yk1,zk1,xk2,yk2,zk2)
	  xss=xmsum(imin,imax,clmd,xk1,yk1,zk1,xk2,yk2,zk2)
	umomd=umomd+cfm*t2(i2)/pi3*xss/2.
	  endif
23      continue
c
        do 21 i2=i1+1,imax                                                       
        fr(i1)=fr(i1)+t2(i2)                                                   
        fr(i2)=fr(i2)+t2(i2)  
	  
	  frp(i1)=frp(i1)+iso(i1)*iso(i2)*t2(i2)*rho0
	  frp(i2)=frp(i2)+iso(i1)*iso(i2)*t2(i2)*rho0
                                                
21      continue                     
22	continue 
	   	                          
        do 258 i=imin,imax
        skine=(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))/2./em(i) 
	tkine=tkine+skine
	Ekin=Ekin+Ekfac*(ccf*fr(i)*rho0)**0.6667
	
	if(i.le.npr)then
	r=(x(i)-apx)**2+(y(i)-apy)**2+(z(i)-apz)**2
	r=sqrt(r)         
	 if(ish.eq.0)then
      Esh=Esh+Esh1/(1.+xexp((r-r01)/ds))
	 else
	  if(iso(i).eq.1)then
	Eshp=Eshp+Esh1p/(1.+xexp((r-r01)/ds))
	  else
	Eshn=Eshn+Esh1n/(1.+xexp((r-r01)/ds))
	  endif
	Esh=Eshp+Eshn
	 endif
	else
	r=(x(i)-atx)**2+(y(i)-aty)**2+(z(i)-atz)**2
	r=sqrt(r)         
	 if(ish.eq.0)then
      Esh=Esh+Esh2/(1.+xexp((r-r02)/ds))
	 else
	  if(iso(i).eq.1)then
	Eshp=Eshp+Esh2p/(1.+xexp((r-r02)/ds))
	  else
	Eshn=Eshn+Esh2n/(1.+xexp((r-r02)/ds))
	  endif
	Esh=Eshp+Eshn
	 endif
	endif
258     continue 

	do 220 i1=imin,imax
	do 230 i2=imin,imax
	  al=2.*(sl(i1)+sl(i2))  
	  sal=(pi*al)**1.5*rho0
	 if(i1.eq.i2) goto 230
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &      (z(i1)-z(i2))*(z(i1)-z(i2))                                            
	ex12=xexp(-dr2/al)
	ussr=ussr+gsym*(6./al-4.*dr2/al**2.)*isy(i1)*isy(i2)*ex12/2./sal
	if(nsym.eq.1) then
	usym=usym+symm*isy(i1)*isy(i2)*ex12/2./sal
	elseif(nsym.eq.2) then
	usym=usym+symm*isy(i1)*isy(i2)*ex12/2./sal*fr(i1)**(xmu-1.)
	elseif(nsym.eq.3) then
	usym=usym+asym*isy(i1)*isy(i2)*ex12/sal+bsym*isy(i1)*isy(i2)*
     &     ex12/sal*fr(i1)**(gamsy-1.)+csym*isy(i1)*isy(i2)*
     &     ex12/sal*fr(i1)**(eta-1.)
	else
	usym=0.
	endif
230	continue
220	continue
	do 1 i=imin,imax
	utwo=utwo+1./2.*alpha*fr(i)
	uthr=uthr+1./(1.+gamma)*beta*fr(i)**gamma
	utau=utau+gtau*fr(i)**eta
	couex=couex+Cex*frp(i)**0.3333
1	continue
	usky=utwo+uthr

	tene=tkine+usky+ucou+usym+ussr+ugradi+utau+couex+umomd+upaul+Esh
	aene= tene*1000./float(nnuc)

	avpot=(usky+ucou+usym+ugradi+utau+couex+umomd+upaul)*1000.
	aveki=Ekin
	Esh=Esh*1000.

        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE search
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	dimension nc(300,300)
	real logx
	do 135 i=1,np
	do 136 j=1,np
	nc(i,j)=0
136	continue
135	continue
      hb=0.197327                                                    
	hb2=hb/2.
	ncount=0
        do 137 i=1,np                                                           
	k=i+1
	do 138 j=k,np
               drx=x(i)-x(j)
	dpx=px(i)-px(j)
		logx=abs(drx*dpx)      
	if(logx.lt.hb2)then
	nc(i,j)=1
	ncount=ncount+1
        endif 
138	continue
137     continue                                                                
	write(25,18476)rtime,ncount
	do 139 i=1,np
	write(25,*)(nc(i,j),j=1,np)
139	continue
18476	format(1x,f9.4,i5)
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE coulomb_modify
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/bb/np,nzp,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	common/gg/nzta,nzpr
	common/cmod/coeft,coefp,epp,ecm,ptcm
	coeft=float(npr)/float(npr+nta)
	coefp=float(nta)/float(npr+nta)
	angm=x00*sqrt(2.*0.938*nta*npr/(npr+nta)*ecm)/0.197327 ! 非相对论的, hbarc = 0.197327 (fm GeV)
	call coulomb_second
	return
	end
c///////////////////////////////////////////////////////////////////////////////
		SUBROUTINE coulomb_second
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/bb/np,nzp,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr
      common/cmod/coeft,coefp,epp,ecm,ptcm
	 pi=3.1415926

	redm=float(nta*npr)/float(nta+npr)
	redm=redm*0.938
	ech12=0.00144*float(nzpr*nzta)
	orbl=x00*sqrt(2.*redm*ecm)
	epsi=1.+2.*ecm*orbl*orbl/redm/ech12**2
	epsi=sqrt(epsi)
	dista=ptcm
	epsi1=1./epsi
	if(epsi1.ge.0.9999)  epsi1=0.9999
	epsi2=1.+orbl**2/redm/ech12/dista
	epsi3=epsi2*epsi1
	if(epsi3.ge.0.9999)  epsi3=0.9999
	thei=pi-acos(epsi1)+acos(epsi3)

	ekin = ecm-ech12/dista-orbl**2./(2.*redm*dista**2)
	rdot = 2./redm*ekin*0.938*0.938
	if(rdot.le.0.) then
	rdot=0.
	else
	rdot=sqrt(rdot)
	endif

	if(x00.eq.0.) then
	thdot=0.
	else
	thdot=sqrt(abs(epsi**2.-epsi2**2.))
	thdot=rdot*orbl**2/redm/ech12/dista**2/thdot
	endif
	zdot=rdot*cos(thei)-dista*sin(thei)*thdot
	xdot=rdot*sin(thei)+dista*cos(thei)*thdot
	zdot=abs(zdot)
	xdot=abs(xdot)
	pzta=zdot*coeft
	pzpr=zdot*coefp
	pxta=xdot*coeft
	pxpr=xdot*coefp
	x0i=dista*sin(thei)
	x0=x0i/2.
	z0i=dista*cos(thei)
	z0=abs(z0i/2.)

	return
	end
c///////////////////////////////////////////////////////////////////////////////
	SUBROUTINE total_energy
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,nzp,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
	  common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp
	  common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc
	tote=0.
	do 1 i=1,np
        tote =tote+sqrt(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i)+gam3*em(i)
     &	*gam3*em(i))
1	continue
	if(iprin.eq.1)then
	print*,'total energy=',tote
	endif
	return
	end                                                  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		function xexp(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if(x.le.-87.) then
	xexp=exp(-87.)
	return
	endif
	
      if(x.ge.87.) then
	xexp=exp(87.)
	return
	endif

	if(x.gt.-87..and.x.lt.87.) then
	xexp=exp(x)
	return
	endif

	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		function xdexp(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	double precision x,xdexp

	if(x.le.-700.D0) then
	xdexp=dexp(-700.D0)
	return
	endif
	
      if(x.ge.700.D0) then
	xdexp=dexp(700.D0)
	return
	endif

	if(x.gt.-700.D0.and.x.lt.700.D0) then
	xdexp=dexp(x)
	return
	endif
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE width0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
	common/bb/np,nzp,ntag,x00  
	common/CW/c1,c2
	common/wavel/sl(500) 
	  
	do 1 i=1,npr
	sl(i)=(c1*float(npr)**0.3333+c2)**2.
1	continue

	do 2 i=npr+1,np
	sl(i)=(c1*float(nta)**0.3333+c2)**2. 
2	continue

	return
	end
cc    fengzhq
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE updtsh
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

	common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     & rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	common/bb/np,nzp,ntag,x00 
	common/gg/nzta,nzpr   
	common/xx/alpha,beta,gamma,symm     
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile                                        
	common/shell/ish,ds,r0p,r0t,r0c,Eshpp,Eshpn,Eshtp,Eshtn,Eshcp,
     &             Eshcn,Ri0
      common/ths/th1,th2,th3
	common/SURF/ VANISH(500)

!TESTMF

	Ri=Ri0
	Rf=r0c-r0t
	shp=1.+1.6228*ds/r0p
	sht=1.+1.6228*ds/r0t
	shc=1.+1.6228*ds/r0c
	nzp=nzpr+nzta

	if(ish.eq.0)then
	OD3=1./3.
	TD3=2./3.
	CB=0.0058
	CS=0.26

	ip=int(npr-nzpr+0.0001)
	jp=nzpr
	Sp=CB*((FS(ip)+GS(jp))/(float(npr)/2.)**TD3-CS*float(npr)**OD3)
	Eshp=Sp*(1.-2.*th1)*xexp(-th1)
	Eshp=Eshp/npr*shp
	ip=int(nta-nzta+0.0001)
	jp=nzta
	St=CB*((FS(ip)+GS(jp))/(float(nta)/2.)**TD3-CS*float(nta)**OD3)
	Esht=St*(1.-2.*th2)*xexp(-th2)
	Esht=Esht/nta*sht
	ip=int(np-nzpr-nzta+0.0001)
	jp=nzpr+nzta
	Sc=CB*((FS(ip)+GS(jp))/(float(np)/2.)**TD3-CS*float(np)**OD3)
	Eshc=Sc*(1.-2.*th3)*xexp(-th3)
	Eshc=Eshc/np*shc
	else
	Eshpp=Eshpp/nzpr*shp
	Eshpn=Eshpn/(npr-nzpr)*shp
	Eshtp=Eshtp/nzta*sht
	Eshtn=Eshtn/(nta-nzta)*sht
	Eshcp=Eshcp/nzp*shc
	Eshcn=Eshcn/(np-nzp)*shc
	endif

	call center(apx,apy,apz,1)
	call center(atx,aty,atz,2) 
	Rpt=(apx-atx)**2+(apy-aty)**2+(apz-atz)**2
	Rpt=sqrt(Rpt)
	r01=r0p*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))
	r02=r0t*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))

	if(ish.eq.0)then
	Esh1=Eshp*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))
	Esh2=Esht*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))
	else
	Esh1p=Eshpp*SWITCH(Ri,Rf,Rpt)+Eshcp*(1.-SWITCH(Ri,Rf,Rpt))
	Esh1n=Eshpn*SWITCH(Ri,Rf,Rpt)+Eshcn*(1.-SWITCH(Ri,Rf,Rpt))	
	Esh2p=Eshtp*SWITCH(Ri,Rf,Rpt)+Eshcp*(1.-SWITCH(Ri,Rf,Rpt))
	Esh2n=Eshtn*SWITCH(Ri,Rf,Rpt)+Eshcn*(1.-SWITCH(Ri,Rf,Rpt))		
	endif

	do 2 i=1,np
		if(VANISH(I).EQ.1)	cycle
	if(i.le.npr)then
	r2=(x(i)-apx)**2+(y(i)-apy)**2+(z(i)-apz)**2
	r=sqrt(r2)
	if(r.eq.0.) goto 2

	 if(ish.eq.0)then
	fx=Esh1/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(x(i)-apx)/r
	fy=Esh1/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(y(i)-apy)/r
	fz=Esh1/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(z(i)-apz)/r
	 else
	  if(iso(i).eq.1)then
	fx=Esh1p/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(x(i)-apx)/r
	fy=Esh1p/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(y(i)-apy)/r
	fz=Esh1p/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(z(i)-apz)/r
	  else
	fx=Esh1n/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(x(i)-apx)/r
	fy=Esh1n/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(y(i)-apy)/r
	fz=Esh1n/ds*xexp((r-r01)/ds)/(1.+xexp((r-r01)/ds))**2*(z(i)-apz)/r
	  endif
	 endif
	else
	r2=(x(i)-atx)**2+(y(i)-aty)**2+(z(i)-atz)**2
	r=sqrt(r2)
	if(r.eq.0.) goto 2
	
	 if(ish.eq.0)then
	fx=Esh2/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(x(i)-atx)/r
	fy=Esh2/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(y(i)-aty)/r
	fz=Esh2/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(z(i)-atz)/r
	 else
	  if(iso(i).eq.1)then
	fx=Esh2p/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(x(i)-atx)/r
	fy=Esh2p/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(y(i)-aty)/r
	fz=Esh2p/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(z(i)-atz)/r
	  else
	fx=Esh2n/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(x(i)-atx)/r
	fy=Esh2n/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(y(i)-aty)/r
	fz=Esh2n/ds*xexp((r-r02)/ds)/(1.+xexp((r-r02)/ds))**2*(z(i)-atz)/r
	  endif
	 endif
	endif

	px(i)=px(i)+fx*dt
	py(i)=py(i)+fy*dt
	pz(i)=pz(i)+fz*dt

2	continue	



	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
	SUBROUTINE constr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500) 
	common/bb/np,nzp,ntag,x00
	common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
	common/ee/mseed,nseed
	common/poli/iflagcro,iflagpau,ipoli
	common/wavel/sl(500) 
	common/psdc/nntag,ffnn

      pi=3.14159265
	hb=0.197327
	xpi=31.4159265

	do 1 i=1,np
	widthr=sl(i)
	widthp=hb**2./widthr/4. 
	widthsr=3.*sqrt(widthr)
	widthsm=3.*sqrt(widthp)
	f1=Fn(i)
	f1=ffnn

	if(f1.gt.1.) then
        do 2 j=i+1,np

	  if(iso(i).eq.iso(j)) then
	xxr=x(i)-x(j)
	yyr=y(i)-y(j)
	zzr=z(i)-z(j)
	ppxr=px(i)-px(j)
	ppyr=py(i)-py(j)
	ppzr=pz(i)-pz(j)
	rrs=sqrt(xxr**2.+yyr**2.+zzr**2.)
	pps=sqrt(ppxr**2.+ppyr**2.+ppzr**2.)

	  if(rrs.le.widthsr.and.pps.le.widthsm) then
	  x1=x(i)                                                                
        y1=y(i)                                                                
        z1=z(i)
	  x2=x(j)                                                                
        y2=y(j)                                                                
        z2=z(j)
	  pz1=pz(i)                                                              
        pz2=pz(j)                                                              
        px1=px(i)                                                              
        px2=px(j)                                                              
        py1=py(i)                                                              
        py2=py(j)                                                              
        em1=em(i)                                                              
        em2=em(j)
	  e1=sqrt(em1*em1+pz1*pz1+px1*px1+py1*py1)                          
        e2=sqrt(em2*em2+pz2*pz2+px2*px2+py2*py2) 
        bx=(px1+px2)/(e1+e2)                                                    
        by=(py1+py2)/(e1+e2)                                                    
        bz=(pz1+pz2)/(e1+e2) 
	  gam=1./sqrt(1.-bx*bx-by*by-bz*bz) 
	  s=(e1+e2)**2-(pz1+pz2)**2-(px1+px2)**2-(py1+py2)**2
	  if(s.lt.3.556896) goto 2	   
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

	  srt=sqrt(s)                  
        s=srt                                   
        as=(3.65*(s-1.866))**6                                               
        a=6.*as/(1.+as)                                                         
        ta=-2*pr*pr
	  xx=RAn(mseed)                                  
        t1=alog((1-xx)*XEXP(a*ta)+xx)/a                    
        c1=1.-t1/ta
        t1=2.*pi*RAn(mseed)
        c2=pzr/pr
        s2=sqrt(1.-c2*c2)                                                       
        if(pxr.eq.0.0.and.pyr.eq.0.0) goto 3
        t2=atan2(pyr,pxr)
3       s1=sqrt(1.-c1*c1)                                                      
        ct1=cos(t1)                                                             
        st1=sin(t1)                                                             
        ct2=cos(t2)                                                             
        st2=sin(t2) 	                                                              
        pzr=pr*(c1*c2-s1*s2*ct1)
        ss=c2*s1*ct1+s2*c1                                                      
        pxr=pr*(ss*ct2-s1*st1*st2)
        pyr=pr*(ss*st2+s1*st1*ct2)

        pb=bx*pxr+by*pyr+bz*pzr                                              
        ea=sqrt(em(i)*em(i)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*pb/b2+gam*ea   
  	  px(i)=gg*bx+pxr                                                        
        py(i)=gg*by+pyr                                                        
        pz(i)=gg*bz+pzr  
	  eb=sqrt(em(j)*em(j)+pxr*pxr+pyr*pyr+pzr*pzr)                          
        gg=(gam-1.)*(-pb)/b2+gam*eb                                             
        px(j)=gg*bx-pxr                                                        
        py(j)=gg*by-pyr                                                        
        pz(j)=gg*bz-pzr

	call pauli12(i,j)
        if(nntag.eq.-1) goto 4  
	                                             
        if(nntag.eq.0)  goto 2 
	
c	call paur34(i,j)
c        if(ntag.eq.0)  goto 2
	
4	  pz(i)=pz1                                                              
        pz(j)=pz2                                                              
        px(i)=px1                                                              
        px(j)=px2                                                              
        py(i)=py1                                                              
        py(j)=py2                                                              
        em(i)=em1                                                              
        em(j)=em2

         endif
         endif
2	continue
	endif	

1	continue	 

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    SUBROUTINE neck
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	common/bb/np,nzp,ntag,x00
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc,numit
	common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt,kfile
	common/prneck/ratio_np,nprogp,nprogn,ntargp,ntargn,rrn
	common/cmod/coeft,coefp,epp,ecm,ptcm
	common/diss/xldis,Ediss

	open(184,file='diarho.dat',status='unknown')

	pi=3.141592654
	rhodef=0.02
	r0p=1.45*float(npr)**0.3333
	r0t=1.55*float(nta)**0.3333
	rstep=0.2	 
	nump=int(r0p/rstep)
	numt=int(r0t/rstep)

      mprogp=0.
	mprogn=0.
	mtargp=0.
	mtargn=0.
	nprogp=0.
	nprogn=0.
	ntargp=0.
	ntargn=0.
	sum=0.
	sum_p=0.
	ratio_np=0.

	sum_pp=0.
	sum_pn=0.
	sum_tp=0.
	sum_tn=0.
	flow_pph=0.
	flow_pnh=0.
	flow_tph=0.
	flow_tnh=0.
	flow_ppv=0.
	flow_pnv=0.
	flow_tpv=0.
	flow_tnv=0.

	call center(xc,yc,zc,0)
	call center(xcprog,ycprog,zcprog,1)
	call center(xctarg,yctarg,zctarg,2)
	r12=(xctarg-xcprog)**2.+(yctarg-ycprog)**2.+(zctarg-zcprog)**2.
	r12=sqrt(r12)
	rneck=r0p+r0t+5.0

c----------dissipation of kinetic energy and angular momentum of relative motion
	x12=(xctarg-xcprog)**2.+(yctarg-ycprog)**2.
	x12=sqrt(x12)
	ARED=float(npr*nta)/float(np)
      xlamda=197.32/SQRT(2.*931.477736*ARED*ecm)
	xldis=x12/xlamda
	Ediss=0.
	do i=1,np
	Ediss=Ediss+(2.*pz(i)**2-px(i)**2-py(i)**2)/2./em(i)
	enddo
	return
c-------------------------------------------------------------------------------
	ntime=rtime
	if(ntime.ge.400) then
	nump=100
	numt=100
	endif
	if(mod(ntime,100).eq.0.and.nru.eq.1) then
	write(184,*) 'The time is :', rtime

	  do i1=-nump,nump
	   zzp=rstep*i1+zcprog
	   do i2=-nump,nump
	   xxp=rstep*i2+xcprog
	   rho0p=densit(xxp,ycprog,zzp)
	write(184,*) zzp,xxp,rho0p
	   enddo 
	  enddo

	  do i1=-numt,numt
	   zzt=rstep*i1+zctarg
	   do i2=-numt,numt
	   xxt=rstep*i2+xctarg
	   rho0t=densit(xxt,yctarg,zzt)
	write(184,*) zzt,xxt,rho0t
	   enddo 
	  enddo
	endif

	if(r12.gt.rneck) return
	dxc=xcprog-xctarg
	dzc=zcprog-zctarg
	if(dxc.eq.0.or.x00.eq.0)then
	theta0=pi/2.
	else
	ac=dzc/dxc
	theta0=atan(ac)
	endif	
	fxcprog=frxif(xcprog,zcprog,xc,zc,theta0)
	fxctarg=frxif(xctarg,zctarg,xc,zc,theta0)
	rhotemp=densit(xcprog,ycprog,zcprog)
      do xx=fxcprog,fxctarg,0.05
	  call frx(xx,xc,zc,theta0,xi,zi)
	  rho=densit(xi,yc,zi)          
	  if(rho.lt.rhodef)return
	  if(rho.lt.rhotemp)then
	  rhotemp=rho
	  xneck=xx
	  endif
	enddo

	rs=0.05
	do i1=1,500
	xt=rs*i1
	rho0=densit(xt,yc,xneck)
	if(rho0.lt.rhodef)exit
	rrn1=xt
	enddo	
	do i2=1,500
	xt=-rs*i2
	rho0=densit(xt,yc,xneck)
	if(rho0.lt.rhodef)exit
	rrn2=xt
	enddo
	rrn=abs(rrn1-rrn2)/2.

	thetamax=theta0
	do 10 i=1,np
      xk=frxif(x(i),z(i),xc,zc,thetamax)-xneck
	zk=frzif(x(i),z(i),xc,zc,thetamax)
	if(abs(zk).le.6.0)then
	if(xk.lt.0.and.id(i).lt.50.and.iso(i).eq.1)mprogp=mprogp+1
	if(xk.lt.0.and.id(i).lt.50.and.iso(i).eq.0)mprogn=mprogn+1
	if(xk.gt.0.and.id(i).gt.50.and.iso(i).eq.1)mtargp=mtargp+1
	if(xk.gt.0.and.id(i).gt.50.and.iso(i).eq.0)mtargn=mtargn+1
	if(xk.gt.0.and.id(i).lt.50.and.iso(i).eq.1)nprogp=nprogp+1
	if(xk.gt.0.and.id(i).lt.50.and.iso(i).eq.0)nprogn=nprogn+1
	if(xk.lt.0.and.id(i).gt.50.and.iso(i).eq.1)ntargp=ntargp+1
	if(xk.lt.0.and.id(i).gt.50.and.iso(i).eq.0)ntargn=ntargn+1

	if(abs(xk).le.2.0)then
	sum=sum+1.0
	if(iso(i).eq.1)sum_p=sum_p+1.0

	ph=frxif(px(i),pz(i),0.0,0.0,thetamax)          ! transverse flow at neck region
	if(id(i).lt.50.and.iso(i).eq.1)flow_pph=flow_pph+ph
	if(id(i).lt.50.and.iso(i).eq.0)flow_pnh=flow_pnh+ph
	if(id(i).gt.50.and.iso(i).eq.1)flow_tph=flow_tph+ph
	if(id(i).gt.50.and.iso(i).eq.0)flow_tnh=flow_tnh+ph

	pzk=frzif(px(i),pz(i),0.0,0.0,thetamax)         ! radial flow at neck region
	pv=(zk*pzk+y(i)*py(i))/sqrt(zk*zk+y(i)*y(i))

	if(id(i).lt.50.and.iso(i).eq.1)flow_ppv=flow_ppv+pv
	if(id(i).lt.50.and.iso(i).eq.0)flow_pnv=flow_pnv+pv
	if(id(i).gt.50.and.iso(i).eq.1)flow_tpv=flow_tpv+pv
	if(id(i).gt.50.and.iso(i).eq.0)flow_tnv=flow_tnv+pv

	if(id(i).lt.50.and.iso(i).eq.1)sum_pp=sum_pp+1
	if(id(i).lt.50.and.iso(i).eq.0)sum_pn=sum_pn+1
	if(id(i).gt.50.and.iso(i).eq.1)sum_tp=sum_tp+1
	if(id(i).gt.50.and.iso(i).eq.0)sum_tn=sum_tn+1

	endif
	endif
10	continue

	if(sum_p.ne.0) ratio_np=(sum-sum_p)/sum_p

	if(sum_pp.ne.0)flow_pph=flow_pph/sum_pp
	if(sum_pn.ne.0)flow_pnh=flow_pnh/sum_pn
	if(sum_tp.ne.0)flow_tph=flow_tph/sum_tp
	if(sum_tn.ne.0)flow_tnh=flow_tnh/sum_tn

	if(sum_pp.ne.0)flow_ppv=flow_ppv/sum_pp
	if(sum_pn.ne.0)flow_pnv=flow_pnv/sum_pn
	if(sum_tp.ne.0)flow_tpv=flow_tpv/sum_tp
	if(sum_tn.ne.0)flow_tnv=flow_tnv/sum_tn

	return
	end
cc   fengzhq
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE center(xc,yc,zc,icase)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp                
     +   ,rrt,px(500),py(500),pz(500),em(500),iso(500),fr(500)
	common/bb/np,nzp,ntag,x00
	common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
	dimension IDkey(500)
	common/SURF/ VANISH(500)

	if(icase.eq.0)then
	imin=1
	imax=np
	mass=np
	endif
	if(icase.eq.1)then
	imin=1
	imax=npr
	mass=npr
	endif
	if(icase.eq.2)then
	imin=npr+1
	imax=np
	mass=nta
	endif

!      subtract single particles
        numout=0
	  do 10 i=imin,imax
		if(VANISH(I).EQ.1)	cycle
	  IDkey(i)=0
	  key=0
	  do 20 j=imin,imax
		if(VANISH(J).EQ.1)	cycle
	  if(i.eq.j)cycle
	  dr=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
	  if(dr.le.3.5)then
	  key=1
	  goto 15
	  endif
20	  continue
15      if(key.eq.0)then
	  IDkey(i)=1           ! one particle escape
	  numout=numout+1
	  endif
10	  continue

	    xc=0.
		yc=0.
		zc=0.
		do 30 i=imin,imax
		if(VANISH(I).EQ.1)	cycle
	    if(IDkey(i).eq.1)cycle
		xc=xc+x(i)/float(mass-numout)
		yc=yc+y(i)/float(mass-numout)
		zc=zc+z(i)/float(mass-numout) 
30	    continue 

	return
	end


cc
	SUBROUTINE  frx(xx,xo,zo,theta,xi,zi)      
      ! rotate xx to origin coodinate by (xo,zo) by anti-clock-theta
       xi=xx*cos(theta)+xo
       zi=xx*sin(theta)+zo
	 return
	 end

      function   frxif(xi,zi,xo,zo,theta)
!      rotate (xi,zi) to final coodinate x by clock-theta degree 
      frxif=(xi-xo)*cos(theta)+(zi-zo)*sin(theta)
      return
	end

      function   frzif(xi,zi,xo,zo,theta)
!      rotate (xi,zi) to final coodinate z by clock-theta degree 
      frzif=-(xi-xo)*sin(theta)+(zi-zo)*cos(theta)
      return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE cluster(NNFUS,NFPRO,NFNEU)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
	common/bb/np,nzp,ntag,x00 
	common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
	common/tt/prtime(8000),nrmin,nrmax,nt,rtime,perc,numit  
	DIMENSION IC(500),XJR(500,2)
	LOGICAL   PROTON 

	RR0=3.5
	PP0=0.2 
	RCR=RR0**2.
	PCR=PP0**2.
!!    Coalescence model is used to construct cluster
CCCCCCCCCCCCCCCCCC* Begin construct cluster *CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	do I=1,NP
	IC(I)=0.
	XJR(I,1)=0.
	XJR(I,2)=0.
	enddo

        III=0                                                          
        do 10 I1=1,NP
            if(IC(I1).EQ.0)then 
		  III=III+1 
		  IC(I1)=III                                             
            endif
	
         do 20 I2=I1+1,NP
	xI1=0.
	yI1=0.
	zI1=0.
	Icc1=0                             
            DRX=X(I1)-X(I2)                                                
            DRY=Y(I1)-Y(I2)                                                
            DRZ=Z(I1)-Z(I2)                                                
            DRR=DRX*DRX+DRY*DRY+DRZ*DRZ                                    
            if(DRR.GT.RCR) goto 20                                        
            DPX=PX(I1)-PX(I2)                                              
            DPY=PY(I1)-PY(I2)                                              
            DPZ=PZ(I1)-PZ(I2)                                              
            DPP=DPX*DPX+DPY*DPY+DPZ*DPZ                                    
            if(DPP.GT.PCR) goto 20                                        
         if(IC(I2).EQ.0)then
	do K1=1,NP
	if(IC(K1).EQ.IC(I1)) then
	xI1=xI1+X(K1)
	yI1=yI1+Y(K1)
	zI1=zI1+Z(K1)
	Icc1=Icc1+1
	endif
	enddo
	xI1=xI1/float(Icc1)
	yI1=yI1/float(Icc1)
	zI1=zI1/float(Icc1)
	xxI1=0.
	yyI1=0.
	zzI1=0.
	do K1=1,NP
	if(IC(K1).EQ.IC(I1)) then
	xxI1=xxI1+(xI1-X(K1))*(xI1-X(K1))
	yyI1=yyI1+(yI1-Y(K1))*(yI1-Y(K1))
	zzI1=zzI1+(zI1-Z(K1))*(zI1-Z(K1))
	endif
	enddo
	crms1=sqrt((xxI1+yyI1+zzI1)/Icc1)

		if(Icc1.ge.2) then
		DDx=xI1-X(I2)	
		DDy=yI1-Y(I2)	
		DDz=zI1-Z(I2)
		DD=sqrt(DDx*DDx+DDy*DDy+DDz*DDz)
		if(DD.GT.(CRMS1+RR0)) goto 20
		endif

	IC(I2)=IC(I1)      
	                                            
         else

              ICK=MAX0(IC(I1),IC(I2))                                        
              ICS=MIN0(IC(I1),IC(I2))                                        
               do J=1,NP
                if(IC(J).EQ.ICK)IC(J)=ICS                                      
               enddo
         endif                                                          
20     continue                                                       
10     continue                                                       
CCCCCCCCCCCCCCCCCC* end construct cluster *CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
        do 30 I=1,NP
        IS=IC(I)                                                       
        PROTON=ISO(I).EQ.1                                          
        if(PROTON)then                                                 
          XJR(IS,1)=XJR(IS,1)+1.
        else                                                           
          XJR(IS,2)=XJR(IS,2)+1.
        endif                                                          
30     continue

	ICMAX=1
	do I=1,NP
	if(IC(I).GT.ICMAX) ICMAX=IC(I)
	enddo

	NFPRO=0
	NFNEU=0
	NNFUS=0
	do I=1,ICMAX
	NN=XJR(I,1)+XJR(I,2)
	  if(NN.NE.0) then
	if(XJR(I,1).EQ.1.)    NFPRO=NFPRO+1
	if(XJR(I,2).EQ.1.)    NFNEU=NFNEU+1
	if(NN.GE.(NP-4))      NNFUS=1

cc	if(rtime.ge.test)     AAN(NN)=AAN(NN)+1.  
!!      AAN(NN)=AAN(NN)+1. 
	  endif
	enddo

	return
	end

	FUNCTION FS(N)
      DIMENSION M(10)
      M(1)=0
	M(2)=2
	M(3)=8
	M(4)=20
	M(5)=28
	M(6)=50
	M(7)=82
	M(8)=126
	M(9)=184
	FD3=5./3.
	do 1 I=2,9
1	if(N.LT.M(I)) goto 2
2	QN=0.6*(M(I)**FD3-M(I-1)**FD3)/(M(I)-M(I-1))
	FS=QN*(N-M(I-1))-0.6*(N**FD3-M(I-1)**FD3)
      return
	end

      FUNCTION GS(K)
	DIMENSION L(10)
      L(1)=0
	L(2)=2
	L(3)=8
	L(4)=20
	L(5)=28
	L(6)=50
	L(7)=82
	L(8)=114
	L(9)=126
	FD3=5./3.
	do 3 I=2,9
3	if(K.LT.L(I)) goto 4
4	QN=0.6*(L(I)**FD3-L(I-1)**FD3)/(L(I)-L(I-1))
	GS=QN*(K-L(I-1))-0.6*(K**FD3-L(I-1)**FD3)
      return
	end

	FUNCTION SWITCH(Ri,Rf,R)
	C0=0.
	C1=0.
	C2=0.
	C3= 10.
	C4=-15.
	C5= 6.
	X0=(R-Rf)/(Ri-Rf)
	SWITCH=C0+C1*X0+C2*X0**2+C3*X0**3+C4*X0**4+C5*X0**5
	return
	end

	FUNCTION Fn(m)
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)                                 
      common/bb/np,nzp,ntag,x00   
      common/wavel/sl(500) 
	common/psdc/nntag,ffnn

      double precision  xxr,yyr,zzr,ppxr,ppyr,ppzr,xxr1,yyr1,zzr1,ppxr1,
     +ppyr1,ppzr1,xxr2,yyr2,zzr2,ppxr2,ppyr2,ppzr2,occupr,occupp
      real Fn ! 设置双精度时会处bug，CHEN SM
     	 
	ffnn=0.
        pi=3.14159265
	hb=0.197327
	widthr=sl(m)
	widthp=hb**2./widthr/4.
	hcr=sqrt(4.*pi*widthr)
	hcp=sqrt(4.*pi*widthp) 
	widthsr=3.*sqrt(widthr)
	widthsm=3.*sqrt(widthp)
	Fn=0.d0

	do 1 j=1,np

	if(iso(m).eq.iso(j))then
	xxr=x(m)-x(j)
	yyr=y(m)-y(j)
	zzr=z(m)-z(j)
	ppxr=px(m)-px(j)
	ppyr=py(m)-py(j)
	ppzr=pz(m)-pz(j)
	xxr1=(hcr/2.+xxr)/sqrt(2.*widthr)
	xxr2=(hcr/2.-xxr)/sqrt(2.*widthr)
	yyr1=(hcr/2.+yyr)/sqrt(2.*widthr)
	yyr2=(hcr/2.-yyr)/sqrt(2.*widthr)
      zzr1=(hcr/2.+zzr)/sqrt(2.*widthr)
	zzr2=(hcr/2.-zzr)/sqrt(2.*widthr)
	ppxr1=(hcp/2.+ppxr)/sqrt(2.*widthp)
	ppxr2=(hcp/2.-ppxr)/sqrt(2.*widthp)
	ppyr1=(hcp/2.+ppyr)/sqrt(2.*widthp)
	ppyr2=(hcp/2.-ppyr)/sqrt(2.*widthp)
	ppzr1=(hcp/2.+ppzr)/sqrt(2.*widthp)
	ppzr2=(hcp/2.-ppzr)/sqrt(2.*widthp)

	occupr=(2.-erfc(xxr1)-erfc(xxr2))*(2.-erfc(yyr1)-erfc(yyr2))*
     +(2.-erfc(zzr1)-erfc(zzr2))
	occupp=(2.-erfc(ppxr1)-erfc(ppxr2))*(2.-erfc(ppyr1)-erfc(ppyr2))*
     +(2.-erfc(ppzr1)-erfc(ppzr2))
	
	if(m.eq.j) then
	Fn=Fn+occupr*occupp/64.
	else
	Fn=Fn+0.5*occupr*occupp/64.
	endif
	endif
1     continue

	ffnn=fn
	return
	end
C***********************************************************************                                                                      
      SUBROUTINE crossnew                                                 
C***********************************************************************     
      common/boite/ann(800),anp(800)

      xpi=31.415926
      do 105 i=1,800                                                    
      e=float(i)                                                       
      if(i.lt.25)then                                                   
      ann(i)=80.6e0                                                     
      anp(i)=381.0e0                                                    
      goto 100                                                         
      endif                                                             
      if(i.lt.40)then                                                   
      ann(i)=5.3107e0+3088.5e0/e-1174.8e0/e/e                           
      anp(i)=6.9466e0+9069.2e0/e-5057.4e0/e/e                           
      goto 100                                                          
      endif                                                             
      if(i.lt.310)then                                                  
      ann(i)=22.429e0-11.148e0/e+93074.e0/e/e                           
      anp(i)=27.147e0+1802.0e0/e+239380.e0/e/e                          
      goto 100                                                          
      endif                                                             
      if(i.le.800)then                                                  
      ann(i)=3.5475e0+.05331e0*e +887.37e0/e                            
      anp(i)=34.5e0                                                     
      goto 100                                                          
      endif                                                             
  100 ann(i)=ann(i)/xpi                                                 
      anp(i)=anp(i)/xpi                                                 
  105 continue

      return
      end
c///////////////////////////////////////////////////////////////////////////////
	function fdb(srt)
	eps=1.E-6
      AVMASS=0.938868
      PIMASS=0.137265
	xms=AVMASS+PIMASS
	xlow=xms**2
	xupp=(srt-AVMASS)**2
	if(xupp.le.xlow) then
	fdb=0.0001
	return
	else
	call floren(xlow,xupp,eps,fdb)
	endif

	return
	end     
cc
	SUBROUTINE floren(A,B,EPS,G)
	real C(5),T(5)
	real A,B,F,G,S,P,H,AA,BB,W,X,Q
	DATA T/-0.90618,-0.53847,0.,0.53847,0.90618/
	DATA C/0.23693,0.47863,0.56889,0.47863,0.23693/

	M=1
	S=(B-A)*0.01
	P=0.
	H=(B-A)/M
	G=0.
	 do l=1,M
	AA=A+(l-1)*H
	BB=A+l*H
	W=0.
	  do j=1,5
           x=((BB-AA)*T(j)+(BB+AA))/2.
           W=W+gams(x)*C(j)
	  enddo
	G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))

 710	P=G
	M=M+1
	H=(B-A)/M
      G=0.
	 do l=1,M
      AA=A+(l-1)*H
      BB=A+l*H
      W=0.
	  do j=1,5
	     x=((BB-AA)*T(j)+(BB+AA))/2.
	     W=W+gams(x)*C(j)
	  enddo
	  G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))
	if((Q.GE.EPS).AND.(ABS(H).GT.ABS(S))) goto 710

	return
	end

      function gams(DMASS)
      PI=3.1415926
      AVMASS=0.938868
      PIMASS=0.137265
      xms=AVMASS+PIMASS
      AUX = 0.25*(DMASS**2-AVMASS**2+PIMASS**2)**2
     &      -(DMASS*PIMASS)**2
        if (AUX .GT. 0.) then
        QV = SQRT(AUX / DMASS**2)
        else
        QV = 0.0000001
        endif
	gams=xms*WIDTHDN(QV,1,1)/((DMASS**2-xms**2)**2+
     &     xms**2*WIDTHDN(QV,1,1)*WIDTHDN(QV,1,1))/PI

	return
	end

c///////////////////////////////////////////////////////////////////////////////
	function  fpw(srt)
	eps=1.E-6
      AVMASS=0.938
      PIMASS=0.139
	xlow =AVMASS+PIMASS
	xupp =srt-AVMASS
	if(xupp.le.xlow) then
	fdb=0.0001
	return
	else
	call fintdel(xlow,xupp,eps,fpw,srt)
	endif

	return
	end     
cc
	SUBROUTINE fintdel(A,B,EPS,G,srt)
	real C(5),T(5)
	real A,B,F,G,S,P,H,AA,BB,W,X,Q
	DATA T/-0.90618,-0.53847,0.,0.53847,0.90618/
	DATA C/0.23693,0.47863,0.56889,0.47863,0.23693/

	M=1
	S=(B-A)*0.01
	P=0.
	H=(B-A)/M
	G=0.
	 do l=1,M
	AA=A+(l-1)*H
	BB=A+l*H
	W=0.
	  do j=1,5
           x=((BB-AA)*T(j)+(BB+AA))/2.
           W=W+gdel(x,srt)*C(j)
	  enddo
	G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))

 710	P=G
	M=M+1
	H=(B-A)/M
      G=0.
	 do l=1,M
      AA=A+(l-1)*H
      BB=A+l*H
      W=0.
	  do j=1,5
	     x=((BB-AA)*T(j)+(BB+AA))/2.
	     W=W+gdel(x,srt)*C(j)
	  enddo
	  G=G+W
	 enddo
	G=G*H/2.
	Q=ABS(G-P)/(1.0+ABS(G))
	if((Q.GE.EPS).AND.(ABS(H).GT.ABS(S))) goto 710

	return
	end

      function gdel(DMASS,srt)
       PI=3.1415926
       AVMASS=0.938
       PIMASS=0.139
       GDWD = 0.115
	 DLMAS = 1.232
	 xdm=4.*GDWD*DLMAS**2/((DMASS**2-DLMAS**2)**2+GDWD**2*DLMAS**2)
  	 pr2=((srt**2+DMASS**2-AVMASS**2)/(2.*srt))**2-DMASS**2

        if (pr2 .GT. 0.) then
        QV = SQRT(pr2)
        else
        QV = 0.0000001
        endif
	  gdel=xdm*QV*DMASS/PI/2.


	return
	end
c///////////////////////////////////////////////////////////////////////////////
* PURPOSE : CALCULATE THE MESON+NUCLEON CROSS SECTION ACCORDING TO 
! THE BREIT-WIGNER FORMULA/(p**2)
* VARIABLE : LA = 1 FOR DELTA(1232) RESONANCE
*            LA = 2 FOR N*(1440) RESONANCE
*            LA = 3 FOR N*(1535) RESONANCE

        REAL FUNCTION XNPION(DM,LA,AVPI,XMAX)
c///////////////////////////////////////////////////////////////////////////////
      PARAMETER(AMN=0.939457,AMP=0.93828)
        AVMASS=0.5*(AMN+AMP)
	  DMRES= 1.232
	  if(LA.EQ.2) DMRES= 1.44
	  if(LA.EQ.3) DMRES= 1.535
        PSTAR2=((DM**2-AVMASS**2+AVPI**2)/(2.*DM))**2-AVPI**2
	  PMESN2=((DMRES**2-AVMASS**2+AVPI**2)/(2.*DMRES))**2-AVPI**2

	  if(PSTAR2.LE.0. .or. PMESN2.LE.0.) then
        XNPION=0.000001
	  return
	  endif

	  QM=sqrt(PSTAR2)
        GAM=WIDTHDN(QM,LA,1)
        F1=0.25*GAM**2/(0.25*GAM**2+(DM-DMRES)**2)
        XNPION=F1*(PMESN2/PSTAR2)*XMAX

        return
        end

        FUNCTION PWIDTH(DMASS)
        AVMASS=0.938868
        PIMASS=0.137265
        AUX = 0.25*(DMASS**2-AVMASS**2+PIMASS**2)**2
     &        -(DMASS*PIMASS)**2
        if (AUX .GT. 0.) then
        QAVAIL = SQRT(AUX / DMASS**2)
        else
        QAVAIL = 0.0000001
        endif
        PWIDTH = 0.47 * QAVAIL**3 /
     &           (PIMASS**2 * (1.+0.6*(QAVAIL/PIMASS)**2))
        return
        end

	  FUNCTION WIDTHDN(Q,ID,N)
      PARAMETER (PSCAL1=.238,PSCAL2=.318,P0REF=.227,PIMASS=0.137265)
      empi2=1./PIMASS**2
	if(ID.EQ.1) then
	GAM= 0.12
	A1 = 22.48
	A2 = 39.69
	A3 = 0.04
      elseif(ID.EQ.2) then
	GAM= 0.2
	A1 = 17.22
	A2 = 39.69
	A3 = 0.09
      elseif(ID.EQ.3) then
	GAM= 0.15
	A1 = 17.22
	A2 = 39.69
	A3 = 0.09
	endif

	  WIDTHDN=0.
	  if(N.EQ.0) then
	WIDTHDN=GAM
	  elseif(N.EQ.1) then
      WIDTHDN = A1*Q**3*GAM/(1.+A2*Q**2)/(A3+Q**2)
	  elseif(N.EQ.2) then
C      Kitazoe et al. PLB 166 (1986) 35  taken from Li-Bauer MSUCL759/91        
      WIDTHDN = 0.47*(Q**3)/(1.+ 0.6*(Q**2 *empi2))*empi2
	  elseif(N.EQ.3) then
      WIDTHDN=GAM*(Q**3/(1.+(Q/PSCAL1)**2+(Q/PSCAL2)**4))
     @ /(P0REF**3/(1.+(P0REF/PSCAL1)**2+(P0REF/PSCAL2)**4))   
	  endif

	  return
	  end	   
C*********************************************************************
C     Lorentz Transformation with velocity (BX,BY,BZ)
C     IN: (A0,AX,AY,AZ) --> OUT: (AP0,APX,APY,APZ)
      SUBROUTINE LORENTZ(BX,BY,BZ,A0,AX,AY,AZ,AP0,APX,APY,APZ)
      REAL*8 BS,G,BFACTOR,BS2
      BS=dble(BX)**2+dble(BY)**2+dble(BZ)**2
      BS2=1.D0-BS
       if(BS2.LE.1D-10) BS2=1.D-10
      G=1.D0/DSQRT(BS2)
      if(BS.LT.0.00000001) then
        BFACTOR=0.D0
      else
        BFACTOR=(G-1.D0)/BS*(BX*AX+BY*AY+BZ*AZ)
      endif
      APX=AX+BFACTOR*BX-G*A0*BX
      APY=AY+BFACTOR*BY-G*A0*BY
      APZ=AZ+BFACTOR*BZ-G*A0*BZ
      AP0=G*(A0-(BX*AX+BY*AY+BZ*AZ))
      return
      end
C************************************************************************
      SUBROUTINE EVENT(P2,DMK,PX,PY,PZ,P0)
      common/ee/mseed,nseed
	PI=3.141592654
c***** generate kaon event in NN cms
      PX=0.0
      PY=0.0
      PZ=0.0
      P0=DMK
      if(P2.LE.0.0) return
      PS=SQRT(P2)
C     do SOME MONTE CARLO
 10   PP=RANFL(nseed)*PS
      F=(1.-PP/PS)*(PP/PS)**2
      if(RANFL(nseed).GT.F)                             goto 10
C     now accept P-value
C     generate angles
      COST=2.*RANFL(nseed)-1.
      SINT=SQRT(1.-COST**2)
      PHI=RANFL(nseed)*2.*PI
      PX=PP*SINT*COS(PHI)
      PY=PP*SINT*SIN(PHI)
      PZ=PP*COST
      P0=SQRT(DMK**2+PP**2)

      return
      end
C***************************************************************
      SUBROUTINE TWOPARTDECAY(PERR,SS,DM1,DM2,PX,PY,PZ,P01,P02)
C A-> a+b decay routine
      common/ee/mseed,nseed
	PI=3.141592654
      PERR=1

        PMA2=((SS-(DM1+DM2)**2)*(SS-(DM1-DM2)**2))/4./SS
        if(PMA2.LE.0.0)                 return
        PMAX=SQRT(PMA2)

        PHI=2.*PI*RANFL(nseed)
        CTETA=2.*RANFL(nseed)-1.
        STETA=SQRT(1.-CTETA**2)
        XX=STETA*COS(PHI)
        YY=STETA*SIN(PHI)
        ZZ=CTETA

        PX=XX*PMAX
        PY=YY*PMAX
        PZ=ZZ*PMAX
      P01=SQRT(DM1**2+PMA2)
      P02=SQRT(DM2**2+PMA2)
      PERR=0

      return
      end
C************************************************************
	SUBROUTINE xmoms(xlm,x1,y1,z1,x2,y2,z2,sx,sy,sz)
	dimension  U(3),W(3),KS(6),C(6),X(6),KK(2,6),DD(4,7)

cc	 Gaussian integration points U(3) and weights W(3)
        U(1)=-0.7746
	  U(2)= 0.
	  U(3)=-u(1)
	  W(1)= 0.555556
	  W(2)= 0.888889
	  W(3)= w(1)
        LL=2
	  KS(1)=LL
	  KS(2)=LL
	  KS(3)=LL
	  KS(4)=LL
	  KS(5)=LL
	  KS(6)=LL

	 M=1
	 N=6
	 DD(1,7)=1.
	 DD(2,7)=1.
	 DD(3,7)=1.
	 DD(4,7)=1.
 50    do 10 J=M,N
     	 A=flow(J,X,x1,y1,z1,x2,y2,z2)
	 B=KS(J)
	 DD(1,J)=0.5*(fupp(J,X,x1,y1,z1,x2,y2,z2)-A)/B
	 C(J)=DD(1,J)+A
	 X(J)=DD(1,J)*U(1)+C(J)
	 DD(2,J)=0.
	 DD(3,J)=0.
	 DD(4,J)=0.
	 KK(1,J)=1
 10    KK(2,J)=1	
       j=N
 20    fx=fnsx(xlm,j,x,x1,y1,z1,x2,y2,z2)
	 fy=fnsy(xlm,j,x,x1,y1,z1,x2,y2,z2)
	 fz=fnsz(xlm,j,x,x1,y1,z1,x2,y2,z2)
	 kb=KK(1,j)
	 DD(2,j)=DD(2,j+1)*DD(1,j+1)*fx*W(kb)+DD(2,j)
	 DD(3,j)=DD(3,j+1)*DD(1,j+1)*fy*W(kb)+DD(3,j)
	 DD(4,j)=DD(4,j+1)*DD(1,j+1)*fz*W(kb)+DD(4,j)
	 KK(1,j)=KK(1,j)+1
	 if(KK(1,j).LE.3)     goto 30
	 if(KK(2,j).LT.KS(j)) goto 40
	 j=j-1
	 if(j.NE.0) goto 20
	 sx=DD(2,1)*DD(1,1)
	 sy=DD(3,1)*DD(1,1)
	 sz=DD(4,1)*DD(1,1)
	return
 40	 KK(2,j)=KK(2,j)+1
	 C(j)=C(j)+DD(1,j)*2.
	 KK(1,j)=1
 30	 kc=KK(1,j)
	 x(j)=dd(1,j)*U(kc)+C(j)
	 if(j.eq.N) goto 20
	 M=j+1
	 goto 50

	return
	end

	function sumom(xlm,x1,y1,z1,x2,y2,z2)
	dimension  U(3),W(3),KS(6),C(6),X(6),KK(2,6),DD(2,7)

cc	 Gaussian integration points U(3) and weights W(3)
        U(1)=-0.7746
	  U(2)= 0.
	  U(3)=-u(1)
	  W(1)= 0.555556
	  W(2)= 0.888889
	  W(3)= w(1)
        LL=2
	  KS(1)=LL
	  KS(2)=LL
	  KS(3)=LL
	  KS(4)=LL
	  KS(5)=LL
	  KS(6)=LL

	 M=1
	 N=6
	 DD(1,7)=1.
	 DD(2,7)=1.
 50    do 10 J=M,N
     	 A=flow(J,x,x1,y1,z1,x2,y2,z2)
	 B=KS(J)
	 DD(1,J)=0.5*(fupp(J,x,x1,y1,z1,x2,y2,z2)-A)/B
	 C(J)=DD(1,J)+A
	 X(J)=DD(1,J)*U(1)+C(J)
	 DD(2,J)=0.
	 KK(1,J)=1
 10    KK(2,J)=1	
       j=N
 20    F=funsm(xlm,j,x,x1,y1,z1,x2,y2,z2)
	 kb=KK(1,j)
	 DD(2,j)=DD(2,j+1)*DD(1,j+1)*F*W(kb)+DD(2,j)
	 KK(1,j)=KK(1,j)+1
	 if(KK(1,j).LE.3)     goto 30
	 if(KK(2,j).LT.KS(j)) goto 40
	 j=j-1
	 if(j.NE.0) goto 20
	 sumom=DD(2,1)*DD(1,1)
	return
 40	 KK(2,j)=KK(2,j)+1
	 C(j)=C(j)+DD(1,j)*2.
	 KK(1,j)=1
 30	 kc=KK(1,j)
	 x(j)=dd(1,j)*U(kc)+C(j)
	 if(j.eq.N) goto 20
	 M=j+1
	 goto 50

	return
	end

	function flow(I,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	cc=2.
	goto(1,2,3,4,5,6), I
 1	 flow=x1-cc
	return
 2     flow=y1-cc
	return
 3     flow=z1-cc
	return
 4     flow=x2-cc
	return
 5     flow=y2-cc
	return
 6     flow=z2-cc
	return
	end

	function fupp(J,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	cc=2.
	goto(1,2,3,4,5,6), J
 1	 fupp=x1+cc
	return
 2     fupp=y1+cc
	return
 3     fupp=z1+cc
	return
 4     fupp=x2+cc
	return
 5     fupp=y2+cc
	return
 6     fupp=z2+cc
	return
	end

	function funsm(xlm,j,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	 cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
	 cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
	 cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
	 cf= xexp(-cx1-cx2)/(1.+cx3/xlm)
      if(j.ne.6) goto 40
       funsm=cf
      return
 40	if(j.ne.5) goto 30
       funsm=1.
	return
 30	if(j.ne.4)goto 20
       funsm=1.
	return
 20   if(j.ne.3)goto 10
       funsm=1.
	return
 10   if(j.ne.2)goto 5
       funsm=1.
      return
 5	 funsm=1.
	end
	
	function fnsx(xlm,j,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	 cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
	 cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
	 cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
	 cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if(j.ne.6) goto 40
       fnsx=cf
      return
 40	if(j.ne.5) goto 30
       fnsx=1.
	return
 30	if(j.ne.4)goto 20
       fnsx=1.
	return
 20   if(j.ne.3)goto 10
       fnsx=1.
	return
 10   if(j.ne.2)goto 5
       fnsx=1.
      return
 5	 fnsx=x1-x(1)
	end		

	function fnsy(xlm,j,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	 cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
	 cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
	 cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
	 cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if(j.ne.6) goto 40
       fnsy=cf
      return
 40	if(j.ne.5) goto 30
       fnsy=1.
	return
 30	if(j.ne.4)goto 20
       fnsy=1.
	return
 20   if(j.ne.3)goto 10
       fnsy=1.
	return
 10   if(j.ne.2)goto 5
       fnsy=y1-x(2)
      return
 5	 fnsy=1.
	end

	function fnsz(xlm,j,x,x1,y1,z1,x2,y2,z2)
	dimension x(6)
	 cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
	 cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
	 cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
	 cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if(j.ne.6) goto 40
       fnsz=cf
      return
 40	if(j.ne.5) goto 30
       fnsz=1.
	return
 30	if(j.ne.4)goto 20
       fnsz=1.
	return
 20   if(j.ne.3)goto 10
       fnsz=z1-x(3)
	return
 10   if(j.ne.2)goto 5
       fnsz=1.
      return
 5	 fnsz=1.
	end

	function xmsum(imin,imax,xlm,x1,y1,z1,x2,y2,z2)
       common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),
     &       pyms(500),pzms(500),emson(500),idms(500),icp(500)
	 common/wavel/sl(500) 
       hb=0.197327
	 xmsum=0.
	 do 11 i1=imin,imax
	    do 22 i2=imin,imax
	  if(i2.EQ.i1) goto 22
	  xmm=sqrt(sl(i1)+sl(i2))/hb
	  xk1=xmm*px(i1)
	  yk1=xmm*py(i1)
	  zk1=xmm*pz(i1)
	  xk2=xmm*px(i2)
	  yk2=xmm*py(i2)
	  zk2=xmm*pz(i2)
	  x11=(x1-xk1)*(x1-xk1)+(y1-yk1)*(y1-yk1)+(z1-zk1)*(z1-zk1)
	  x22=(x2-xk2)*(x2-xk2)+(y2-yk2)*(y2-yk2)+(z2-zk2)*(z2-zk2)
        pk2=(xk1-xk2)*(xk1-xk2)+(yk1-yk2)*(yk1-yk2)+(zk1-zk2)*(zk1-zk2)
	  xmsum=xmsum+xexp(-x11-x22)/(1.+pk2/xlm)
 22	   continue
 11	 continue

	return
	end

	SUBROUTINE xpart(imin,imax,xlm,x1,y1,z1,x2,y2,z2,sx,sy,sz)
       common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),
     &       pyms(500),pzms(500),emson(500),idms(500),icp(500)
	 common/wavel/sl(500) 
	common/SURF/ VANISH(500)

       hb=0.197327
	 sx=0.
	 sy=0.
	 sz=0.
	 do 11 i1=imin,imax
	if(VANISH(I1).EQ.1)	cycle
	    do 22 i2=imin,imax
			if(VANISH(I2).EQ.1)	cycle
	  if(i2.EQ.i1) goto 22
	  xmm=sqrt(sl(i1)+sl(i2))/hb
	  xk1=xmm*px(i1)
	  yk1=xmm*py(i1)
	  zk1=xmm*pz(i1)
	  xk2=xmm*px(i2)
	  yk2=xmm*py(i2)
	  zk2=xmm*pz(i2)
	  x11=(x1-xk1)*(x1-xk1)+(y1-yk1)*(y1-yk1)+(z1-zk1)*(z1-zk1)
	  x22=(x2-xk2)*(x2-xk2)+(y2-yk2)*(y2-yk2)+(z2-zk2)*(z2-zk2)
        pk2=(xk1-xk2)*(xk1-xk2)+(yk1-yk2)*(yk1-yk2)+(zk1-zk2)*(zk1-zk2)
	  cff=xexp(-x11-x22)/(1.+pk2/xlm)
	  sx=sx+cff*(x1-xk1)
	  sy=sy+cff*(y1-yk1)
	  sz=sz+cff*(z1-zk1)
 22	   continue
 11	 continue

	return
	end
c****************************************************************************	
	  SUBROUTINE matrixs
	 parameter (rmass=0.938,pmass=0.138,pi=3.1416)
	 parameter(hb=0.197327,nsmax=400,delsigs=0.01,nres=30,sigs0=2.014,
     &           delmassdi=0.01, dimimass0=1.076,srtseta=2.51724)
       parameter(nmassmax=400,etafac1=5.,etafac2=2.5)
	 parameter (dsrts=0.01,isrtm=200,nsres=21,numchs=7)
  	 common/matx/xmatr(4,0:nsmax)

      real sigtot,mdel,dsdm,sigma,ps(4),maxcontr,integral
      real sigmasp(4),pini,fac,pss
      real s12(0:nsmax),s12n(0:nsmax),s32(0:nsmax),s32n(0:nsmax),
     &     sigma1pi,sigmapir(4,0:nsmax),sppom(0:nsmax)
      real s120(0:nsmax),s12n0(0:nsmax),s320(0:nsmax),s32n0(0:nsmax)
      real sigmapir0(4,0:nsmax),dimixsec(0:nsmax),sigmasr(2,0:isrtm)
	real dmdimi(2,0:nsmax,0:nsmax),nrmassmax(2,2:nres+1,0:nsmax),
     &     rrmassmax(2,1,0:nsmax),sigmastpi(4,0:nsmax)

      do i=0,nsmax
         s12(i)=0.
         s32(i)=0.
         s12n(i)=0.
         s32n(i)=0.
         s120(i)=0.
         s12n0(i)=0.
         s320(i)=0.
         s32n0(i)=0.
      enddo

*Delta cross sections
      do i=0,nsmax
         srts=float(i)*delsigs+sigs0
         sigtot=0.
         sigmaom=0.
         do j=0,nmassmax
            mdel=dimimass0+float(j)*delmassdi
            
            call dimidm(dsdm,sigma,srts,mdel)

            if(i.eq.1) then
               dsdm=0.
               sigma=0.
            endif
*dsdm for p p -> n Delta++
            dmdimi(1,i,j)=dsdm
*sigma for n Delta++ -> p p
            dmdimi(2,i,j)=sigma
            sigtot=sigtot+dsdm*delmassdi
         enddo
         dimixsec(i)=sigtot

         pini=sqrt(srts**2/4.-rmass**2)
         s32(i)=sigtot/3./pini
         s32n(i)=s32(i)
c         if(irres(1).eq.0) then
c            s320(i)=s32(i)
c            s32n0(i)=s32n(i)
c         endif
      enddo

*N N -> N N pi
      do i=0,nsmax
*index:
*1: p p -> p p pi0
*2: p n -> n n pi+
*3: p n -> p n pi0
*4: p p -> p n pi+

*resonance contributions
         sigmapir(1,i)=2./3.*s32(i)+1./3.*s12(i)
         sigmapir(4,i)=10./3.*s32(i)+2./3.*s12(i)
         sigmapir(2,i)=1./3.*s32n(i)+1./3.*s12n(i)
         sigmapir(3,i)=4./3.*s32n(i)+1./3.*s12n(i)
*resonance contributions of resonances that are not propagated
         sigmapir0(1,i)=2./3.*s320(i)+1./3.*s120(i)
         sigmapir0(4,i)=10./3.*s320(i)+2./3.*s120(i)
         sigmapir0(2,i)=1./3.*s32n0(i)+1./3.*s12n0(i)
         sigmapir0(3,i)=4./3.*s32n0(i)+1./3.*s12n0(i)
*
*background
         srts=float(i)*delsigs+sigs0
         call sstatepi(srts,sigmasp)
*this routine gives the background contribution if all resonances are
*included, therefore missing resonance contributions have to be added

         do k=1,4
            sigmastpi(k,i)=sigmasp(k)+sigmapir0(k,i)
         enddo

*detailed balance: pi N N -> N N
         call bops3(pss,srts,rmass,rmass,pmass)
*pss=\int dm12 dm23 /s (!)
         if(pss.gt.0) then
            do k=1,4
*factor for identical particles
               fac=1.
               if(k.eq.2) fac=2.
               if(k.eq.4) fac=0.5
               pini=sqrt((srts**2-4.*rmass**2)/4.)
               xmatr(k,i)=fac*64.*(2.*pi)**3*pini*srts*
     &              sigmastpi(k,i)/pss/(pi*10.)
            enddo
         else
*if potentials are also included for the calculation of production
*cross section the following should be changed
            do k=1,4
               xmatr(k,i)=0
            enddo
         endif
      enddo

	return 
	end

      SUBROUTINE dimidm(dsdm,sigma,srts,mdel)
c*************   the cross sections for N N <-> Delta N
      implicit none
      real  pi,mdel,mn,srts,m12,t,pi2,pf2,lamda,mpi,
     &     fps,gp,dmass,gdelt,qr2,c,q2,g,f,cost,u,m1,m2,
     &     pt,pu,ft,fu,dsdodm,z,q2t,q2u,qr2t,qr2u,zu,zt,fp,dsdo,
     &     sigma,fac,ratio(10)
      real dsdm
      integer i
      integer nt

      parameter(mn=0.938,pi=3.141592654,mpi=0.14)
      parameter(fps=2.202,fp=1.008,dmass=1.232,gdelt=0.12)
      parameter(lamda= 0.63)
      parameter(nt=100)

      if(srts-2.*mn-mpi.lt.0) then
         dsdm = 0.0
         sigma = 0.0
         return
      endif
      gp=fp*2.*mn/mpi
      qr2=(dmass**2-(mn-mpi)**2)*(dmass**2-(mn+mpi)**2)/
     &     (4*dmass**2)
      c =0.3

      pi2=srts**2/4.-mn**2
      
   
      if(mdel.gt.srts-mn) then
         dsdm=0.
         sigma=0.
         return
      endif
      if(mdel.lt.mn+mpi) then
         f=0.
      else
         q2=(mdel**2-(mn-mpi)**2)*(mdel**2-(mn+mpi)**2)
     &        /(4*mdel**2)

         z=(qr2+c**2)/(q2+c**2)
         z=z**2
         g=sqrt((q2/qr2)**3)*dmass/mdel*z*gdelt

c         call manley(mdel,2,0,0,ratio,g,0,200)

         f=1./pi*mdel*g/((mdel**2-dmass**2)**2+(mdel*g)**2)     

      endif
      pf2=max((srts**2-(mdel+mn)**2)*(srts**2-(mdel-mn)**2)
     &     /(4.*srts**2),0.)

      dsdm=0.
      sigma=0.
      do i=-nt,nt

         cost=float(i)/float(nt)

         t=mdel**2+mn**2-2*sqrt(mdel**2+pf2)*sqrt(mn**2+pi2)
     &        +2*cost*sqrt(pi2*pf2)
         u=2.*mn**2-2*sqrt(mn**2+pf2)*sqrt(mn**2+pi2)
     &        -2*cost*sqrt(pi2*pf2)

         qr2t=(dmass**2+mn**2-t)**2/(4.*dmass**2)-mn**2
         q2t=(mdel**2+mn**2-t)**2/(4.*mdel**2)-mn**2
         zt=(qr2t+c**2)/(q2t+c**2)
         qr2u=(dmass**2+mn**2-u)**2/(4.*dmass**2)-mn**2
         q2u=(mdel**2+mn**2-u)**2/(4.*mdel**2)-mn**2
         zu=(qr2u+c**2)/(q2u+c**2)


         m1=t*(t-(mdel-mn)**2)*((mdel+mn)**2-t)**2/(3.*mdel**2)
         m2=u*(u-(mdel-mn)**2)*((mdel+mn)**2-u)**2/(3.*mdel**2)
         m12=1./(2.*mdel**2)*((t*u+(mdel**2-mn**2)*(t+u)-mdel**4+
     &        mn**4)*
     &        (t*u+mn*(mn+mdel)*(mdel**2-mn**2))-1./3.*(t*u-
     &        (mdel+mn)**2*(t+u)+(mn+mdel)**4)*(t*u-mn*(mdel-mn)*
     &        (mdel**2-mn**2)))
       
         ft=(lamda**2-mpi**2)/(lamda**2-t)
         fu=(lamda**2-mpi**2)/(lamda**2-u)
         pt=1./(t-mpi**2)
         pu=1./(u-mpi**2)

*p p -> n Delta++
         dsdodm=1./(32.*pi*srts**2)*sqrt(pf2)*(fps*gp/mpi)**2*
     &        (ft**4*zt*pt**2*m1+fu**4*zu*pu**2*m2+
     &        ft**2*fu**2*sqrt(zu*zt)*pt*pu*m12)*
     &        2.*mdel*f*0.389
*n Delta++ -> p p, factor 1./4. bec. of spin and identical particles
         dsdo=1./(32.*pi*srts**2)*sqrt(pi2)*(fps*gp/mpi)**2*
     &        (ft**4*zt*pt**2*m1+fu**4*zu*pu**2*m2+
     &        ft**2*fu**2*sqrt(zu*zt)*pt*pu*m12)*
     &        0.389/4.
         fac=1.
         if(abs(i).eq.nt) fac=0.5
         dsdm=dsdm+dsdodm/float(nt)*fac
         sigma=sigma+dsdo/float(nt)*fac
      enddo
      return
      end
***********************************************************************
      SUBROUTINE dimido(dsdo,srts,mdel,cost)
c*************   the cross sections for N N <-> Delta N
      implicit none
      real  pi,mdel,mn,srts,m12,t,pi2,pf2,lamda,mpi,
     &     fps,gp,dmass,c,q2,cost,u,m1,m2,gdelt,
     &     pt,pu,ft,fu,dsdodm,z,q2t,q2u,qr2t,qr2u,zu,zt,fp,dsdo,
     &     fac
      integer i
      integer nt

      parameter(mn=0.938,pi=3.141592654,mpi=0.14)
      parameter(fps=2.202,fp=1.008,dmass=1.232,gdelt=0.12)
      parameter(lamda= 0.63)
      parameter(nt=100)

      if(srts-2.*mn-mpi.lt.0) then
         dsdo=0.0
         return
      endif
      gp=fp*2.*mn/mpi
      c =0.3

      pi2=srts**2/4.-mn**2
      
   
      if(mdel.gt.srts-mn) then
         dsdo=0.
         return
      endif
      pf2=max((srts**2-(mdel+mn)**2)*(srts**2-(mdel-mn)**2)
     &     /(4.*srts**2),0.)



      t=mdel**2+mn**2-2*sqrt(mdel**2+pf2)*sqrt(mn**2+pi2)
     &     +2*cost*sqrt(pi2*pf2)
      u=2.*mn**2-2*sqrt(mn**2+pf2)*sqrt(mn**2+pi2)
     &     -2*cost*sqrt(pi2*pf2)
      
      qr2t=(dmass**2+mn**2-t)**2/(4.*dmass**2)-mn**2
      q2t=(mdel**2+mn**2-t)**2/(4.*mdel**2)-mn**2
      zt=(qr2t+c**2)/(q2t+c**2)
      qr2u=(dmass**2+mn**2-u)**2/(4.*dmass**2)-mn**2
      q2u=(mdel**2+mn**2-u)**2/(4.*mdel**2)-mn**2
      zu=(qr2u+c**2)/(q2u+c**2)


      m1=t*(t-(mdel-mn)**2)*((mdel+mn)**2-t)**2/(3.*mdel**2)
      m2=u*(u-(mdel-mn)**2)*((mdel+mn)**2-u)**2/(3.*mdel**2)
      m12=1./(2.*mdel**2)*((t*u+(mdel**2-mn**2)*(t+u)-mdel**4+
     &     mn**4)*
     &     (t*u+mn*(mn+mdel)*(mdel**2-mn**2))-1./3.*(t*u-
     &     (mdel+mn)**2*(t+u)+(mn+mdel)**4)*(t*u-mn*(mdel-mn)*
     &     (mdel**2-mn**2)))
      
      ft=(lamda**2-mpi**2)/(lamda**2-t)
      fu=(lamda**2-mpi**2)/(lamda**2-u)
      pt=1./(t-mpi**2)
      pu=1./(u-mpi**2)

*n Delta++ -> p p, factor 1./4. bec. of spin and identical particles
      if(pf2.gt.0) then
         dsdo=1./(32.*pi*srts**2)*sqrt(pi2/pf2)*(fps*gp/mpi)**2*
     &        (ft**4*zt*pt**2*m1+fu**4*zu*pu**2*m2+
     &        ft**2*fu**2*sqrt(zu*zt)*pt*pu*m12)*
     &        0.389/4.
      else
         dsdo=0.
      endif
      return
      end

c**********************************************************************
      SUBROUTINE bops3(ps,srts,mass1,mass2,mass3)
      implicit none
      real m12,e3star,e1star,m13min,m13max,m13,srts,pcm,
     &     pi,ps,mass1,mass2,mass3
      parameter(pi=3.141592654)
      integer l,i,j,nm12,nm13,idfile
      nm12=100
cc    Int. von  von m12 [(m1+m2)**2,(srts-m3)**2]

        ps=0
      if(srts.LE.mass1+mass2+mass3)   return

      do i=1,nm12
         m12=(float(i)-0.5)/float(nm12)*((srts-mass3)**2-
     &        (mass1+mass2)**2)+(mass1+mass2)**2
*Int. von m13
        m12=max(m12,0.)
        e3star=(srts**2-m12-mass3**2)/2./sqrt(m12)
        e1star=(m12+mass1**2-mass2**2)/2./sqrt(m12)
        m13min=(e1star+e3star)**2-(sqrt(max(e1star**2-mass1**2,0.))+
     &       sqrt(max(e3star**2-mass3**2,0.)))**2
        m13max=(e1star+e3star)**2-(sqrt(max(e1star**2-mass1**2,0.))-
     &       sqrt(max(e3star**2-mass3**2,0.)))**2
*        do j=1,nm13
*          m13=(float(j)-0.5)/float(nm13)*(m13max-m13min)+m13min
          
          ps=ps+1./srts**2*(m13max-m13min)*((srts-mass3)**2-
     &         (mass1+mass2)**2)/float(nm12)
*        enddo
      enddo

	return
      end

******************************************************************************
      SUBROUTINE sstatepi(srts,sigma)
cc	 calculation of N N -> N N pi background

      implicit none
      real, intent(in) :: srts
      real, intent(out), dimension(1:4) :: sigma

      real xfit, s0charged,par1(11),par2(7),par3(7),par4(10)
      real s0uncharged
      integer i,k 

*index:
*1: p p -> p p pi0
*2: p n -> n n pi+
*3: p n -> p n pi0
*4: p p -> p n pi+
 
* parameters of the fit:

      data (par1(i), i=1,11)/ 0.00421147, 0.987631, -71.7949, 2073.54, 
     & -24143.3, 130013., -319134., -16.1583, 306499., 0.154771,
     & 0.000823241/
      
      data (par2(i), i=1,7) / 0.214637, -0.264214, 114.609, -835.136,
     & 1615.66, 18.485, -54.2815/

      data (par3(i), i=1,7)/ 10.6202,-178.313, 958.58, -40.171, 
     & 76.1355, 336.747, -425.802/

      data (par4(i), i=1,10)/ 2.50696,-38.0018, 186.721, 14.5268, 
     &  -17.8739, 245.919,-319.183, 0.111811, 7.78858, -0.752894/

      s0charged = 2.018
      s0uncharged=2.012

** sigma for PiZero in final state      
      xfit = srts-s0uncharged
      if(xfit.gt.0.0) then
       sigma(1)= 
     &  exp(par1(8)*xfit)*(xfit*par1(1)+xfit**2*par1(2)+xfit**3*par1(3)+
     &    xfit**4*par1(4)+xfit**5*par1(5)+xfit**6*par1(6)+
     &    xfit**7*par1(7)+xfit**8*par1(9))/(par1(11)+(xfit-par1(10))**2)    
     
       sigma(3)=
     &  exp(par3(6)*xfit+par3(7)*xfit**1.2)*
     &   (par3(1)+xfit*par3(2)+xfit**2*par3(3))/
     &   (par3(5)+(xfit-par3(4))**2)    
 
      else
            sigma(1)=0.
            sigma(3)=0.
      endif

** sigma for charged Pion in Final state
      xfit = srts-s0charged
      if(xfit.gt.0.0) then
       sigma(2)=
     &  exp(par2(6)*xfit+par2(7)*xfit**2)*
     &    (xfit*par2(1)+xfit**2*par2(2)+xfit**3*par2(3)+
     &    xfit**4*par2(4)+xfit**5*par2(5))    
       
       sigma(4)=
     &  exp(par4(6)*xfit+par4(7)*xfit**1.2)*
     &  (par4(1)+xfit*par4(2)+xfit**2*par4(3)+
     &    xfit**3*par4(4)+xfit**4*par4(5))/
     &  (par4(9)*xfit**par4(10)+(xfit-par4(8))**2)    
      
      else
           sigma(2)=0.
           sigma(4)=0.
      endif

      return          
      end
!*****************************************************************


      SUBROUTINE sstatepiOld(srts,sigma)
*SUBROUTINE for calculation of N N -> N N pi
      implicit none
      real x(5,3) , xfit, s0,srts,sigma(4) 
      integer i,k 
*index:
*1: p p -> p p pi0
*2: p n -> n n pi+
*3: p n -> p n pi0
*4: p p -> p n pi+
       
      data (x(i,1), i=1,5) / 61.32865, 6.1777919, 1.520738,3.481812,
     &     2.503206/ 
      data (x(i,2), i=1,5) / 24.94817,1.928099,3.300971,2.4472713e-03,
     &     0.8496530/ 
      data (x(i,3), i=1,5) / 7.25069,2.306925, 0.883237, 
     &                     3.641924, 6.6390312e-05 /

      s0 = 2.015 
      xfit = (srts-s0)*5
      if(xfit.gt.0.0) then
         sigma(1)= 
     &        x(1,1)* xfit**x(2,1)*exp(-(x(3,1)*xfit**x(4,1)
     &        +x(5,1)*xfit))
         sigma(2)=
     &        x(1,2)* xfit**x(2,2)*exp(-(x(3,2)*xfit**x(4,2)
     &        +x(5,2)*xfit))
         sigma(3)=
     &        x(1,3)* xfit**x(2,3)*exp(-(x(3,3)*xfit**x(4,3)
     &        +x(5,3)*xfit))
         sigma(4)=2*sigma(1)
      else
         do k=1,4
            sigma(k)=0.
         enddo
      endif
      return          
      end  
***************************************************************
       function sigmega(srts)
c****  cross section p p -> p p omega in mb
	 parameter (rmass=0.938,xomeg=0.782)
      real srts,sigmega,sqrts0,s0,s,x

      sqrts0=2.*rmass+xomeg
      s0=sqrts0**2
      s=srts**2
      if(s.gt.s0) then
         x=s/s0      
         sigmega=2.5*(x-1.)**1.47*x**(-1.11)
      else
         sigmega=0.
      endif

      return
      end


	!获得此核子的总能量(MeV)
	FUNCTION GET_TOTAL_ENERGY(I_L)

	PARAMETER (PI=3.14159265)
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
      common/bb/np,nzp,ntag,x00 
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf,
     &      potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp 
      common/xx/alpha,beta,gamma,symm                                                 
	common/wavel/sl(500)   
	common/SURF/ VANISH(500)

	REAL	KIN_E, GET_TOTAL_ENERGY
	INTEGER, INTENT(IN) :: I_L

	KIN_E = 0.
	GET_TOTAL_ENERGY = 0.

	!密度计算与库伦交换准备
	AL = 2.*(2*(0.08*float(NP)**0.3333+0.92)**2.)
	SAL = (PI*AL)**1.5*RHO0
      cex=-3./4.*0.00144*(3/pi)**0.3333

	!求fr, frp及复杂能项
	FR1 = 0.
	FRP1 = 0.
	UCOU = 0.
	USSR = 0.
	UGRD = 0.
	USYM = 0.
	UMOD = 0.	

	do I = 1, NP

		if(VANISH(I).EQ.1)	cycle

		if(I.EQ.I_L)	cycle
		
		DR = (X(I) - X(I_L))**2 + (Y(I) - Y(I_L))**2 + 
     &		(Z(I) - Z(I_L))**2
		T2 = XEXP(-DR/AL)/SAL
		FR1 = FR1 + T2
		FRP1 = FRP1 + ISO(I_L)*ISO(I)*T2*RHO0

		!库伦能
		IDR22 = DR*40. + 1.5
		IDR2 = MIN0(IDR22, 100000)
		NCCC = IH(IDR2)
		UCOU = UCOU + ISO(I_L)*ISO(I)*SCOU(NCCC,AL)*8.		

		!表面能对称
		USSR = USSR + gsym*(6./al-4.*dr/al**2.)*isy(I_L)*isy(I)*t2

		!表面能
		UGRD = UGRD + gsur*(6./al-4.*dr/al**2.)*t2

		!对称能
		USYM = USYM + symm*isy(I_L)*isy(I)*t2

		!动量相关
		PR = (PX(I) - PX(I_L))**2 + (PY(I) - PY(I_L))**2 + 
     &		(PZ(I) - PZ(I_L))**2
		PBR = PR*CMEPS + 1.
		if(ISO(I_L).EQ.ISO(I)) then
			CFM = CMNN
		else
			CFM = CMNP
		endif
		if(MOMP.EQ.0)	CFM = 0.
		UMOD = UMOD + cfm*alog(pbr)*alog(pbr)*t2/2.
	
	enddo

	!求正比于fr或frp的能项
	UTAU = 2.*GTAU*FR1**ETA
	UTWO = 2.*0.5*ALPHA*FR1
	UTHR = 2.*1./(1. + gamma)*beta*fr1**gamma
	COUEX = 2.*CEX*FRP1**0.3333

	!动能
	KIN_E = PX(I_L)**2 + PY(I_L)**2 + PZ(I_L)**2 + 0.938**2
	KIN_E = SQRT(KIN_E) - 0.938
C 	KIN_E = (PX(I_L)**2 + PY(I_L)**2 + PZ(I_L)**2)/2./EM(I_L)

	!返回总能量
	GET_TOTAL_ENERGY = KIN_E + UCOU + UGRD + USSR +USYM +UMOD
     &					+ UTWO +UTHR +COUEX +UTAU
	GET_TOTAL_ENERGY = GET_TOTAL_ENERGY*1000.

	return
	end


	function scou(ncc,al)

	  if(ncc .EQ. 1)    then
	  scou=0.
	  return
	  endif
	  if(ncc .EQ. 5000) then
	  scou=0.
	  return
	  endif

        ceff= 0.5
        c4= sqrt(al)
        c5= 0.5*0.00144*ceff*ceff

        dr=(float(ncc)-1.)/100.
        t1=dr/c4 
        call  erf(t1,t3)                                                        
	  scou= c5/dr*t3

	return
      end

c///////////////////////////////////////////////////////////////////////////////
c     （1）现在TREE(0:STORE_MAX,2:SURF_MAX_A)是用来装各级结团的ID的
c     （2）现在CLUST_STOR(0:STORE_MAX,SURF_MAX_A)是用来装各结团的核子的
	SUBROUTINE	CONSTR_LOOP
c///////////////////////////////////////////////////////////////////////////////
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
	INTEGER, PARAMETER :: PICK_MAX = 233
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),
     &			TOT_CLUST_NUM, He4_ERROR
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),pyms(500),
     &       pzms(500),emson(500),idms(500),icp(500),idanti(500)
	common/SURF_LOOP/NUM_RANK(SURF_MAX_A)
	INTEGER	EXTRA(PICK_MAX), PICKED, IDD, JMAX, RANK1, JJJ
	INTEGER	TEMP_NEW(STORE_MAX), IDM(STORE_MAX), CURR_N, THIS
	INTEGER	COMA(STORE_MAX), COMB(STORE_MAX), COMC(STORE_MAX)
	INTEGER	COMD(STORE_MAX), IN_CLUST(1:(SURF_MAX_A-1))
	
      ! 主循环：按核子数逐级构造
	do	I = 1, SURF_MAX_A - 1

          ! 阶段1：收集候选核子（抓新核子）
          CURR_N = 0 ! 当前结团核子总数
		IDM(1:STORE_MAX) = 0
		JMAX = INT(NUM_RANK(I)) ! 上一轮的最大候选核子数
          
          ! 对某一核子数的各结团循环
		do	J = 1, JMAX 
			IDD = TREE(J, I) ! 含有I+1个核子的第J种结团的ID
			IN_CLUST(1:SURF_MAX_A - 1) = 
     &		CLUST_STOR(IDD, 1:SURF_MAX_A - 1) ! IID母结团的各核子ID
			call PICK_NUCLEON(EXTRA, IN_CLUST, PICKED)
              
              ! 记录找到的候选核子及其母簇信息
			do	K = CURR_N + 1, CURR_N + PICKED
                  if(K.GT.STORE_MAX)exit ! 这里直接截断，有待检查，CHEN SM
				IDM(K) = IDD ! 记录母簇ID
				THIS = K - CURR_N
				TEMP_NEW(K) = EXTRA(THIS) ! 更新当前候选核子ID，
                                            ! 注意某种结团的新候选核子ID不会覆盖前一种结团的候选核子ID
              enddo
              
			CURR_N = CURR_N + PICKED
              if(CURR_N.GT.STORE_MAX) CURR_N = STORE_MAX ! 防止数组越界，有待检查，CHEN SM
          enddo
		
          if(CURR_N.EQ.0)EXIT ! 没有找到候选核子，终止构造
		
		! 出来之后直接在这个地方sort them
          ! 阶段2：去重处理 - 使用多项式哈希方法
		COMA(1:STORE_MAX) = 0
		COMB(1:STORE_MAX) = 0
		COMC(1:STORE_MAX) = 0
		COMD(1:STORE_MAX) = 0
		do	J = 1, CURR_N
              ! 计算母结团核子的贡献
			do	K = 1, I
				!COMA(J) = COMA(J) + CLUST_STOR(IDM(J), I)
				!COMB(J) = COMB(J) + CLUST_STOR(IDM(J), I)**2
				!COMC(J) = COMC(J) + CLUST_STOR(IDM(J), I)**3
				!COMD(J) = COMD(J) + CLUST_STOR(IDM(J), I)**4
                  COMA(J) = COMA(J) + CLUST_STOR(IDM(J), K)
				COMB(J) = COMB(J) + CLUST_STOR(IDM(J), K)**2
				COMC(J) = COMC(J) + CLUST_STOR(IDM(J), K)**3
				COMD(J) = COMD(J) + CLUST_STOR(IDM(J), K)**4
              enddo
              ! 加上候选核子的贡献，得到每个候选核子的哈希值，哈希值相同的属于同一结团
			COMA(J) = COMA(J) + TEMP_NEW(J)
			COMB(J) = COMB(J) + TEMP_NEW(J)**2
			COMC(J) = COMC(J) + TEMP_NEW(J)**3
			COMD(J) = COMD(J) + TEMP_NEW(J)**4
          enddo     
          ! 如果四个幂和都相同，说明是相同的核子集合，即属于同一结团       
		do	J = 1, CURR_N - 1
			do	K = J + 1, CURR_N
				if(TEMP_NEW(K).EQ.-1)cycle
				if(COMA(J).EQ.COMA(K).AND.COMB(J).EQ.COMB(K)
     &			.AND.COMC(J).EQ.COMC(K).AND.COMD(J).EQ.COMD(K))then
					TEMP_NEW(K) = -1 ! 标记为重复
				endif
			enddo
		enddo ! 最终每种结团只保留一个

		! 阶段3：储存结团
		NUM_RANK(I + 1) = 0 ! 存储某一核子数的结团种类数，NUM_RANK(2)就是核子数为2的结团种类数
		do	J = 1, CURR_N
			if(TEMP_NEW(J).EQ.-1)cycle ! 跳过重复项
              ! 排除5Be（4个质子+1个中子）
			if(I.EQ.4)then
				JJJ = 0
				do JJ = 1, 4
					if(ISY(INT(CLUST_STOR(IDM(J), JJ)))
     &				    .EQ.-1) JJJ = JJJ + 1
				enddo
				if(ISY(TEMP_NEW(J)).EQ.-1) JJJ = JJJ + 1
				if(JJJ.EQ.4)cycle
			endif

			NUM_RANK(I + 1) = NUM_RANK(I + 1) + 1 ! NUM_RANK(I+1) = 1, CURR_N
              
              if(TOT_CLUST_NUM.lt.STORE_MAX)then
			    TOT_CLUST_NUM = TOT_CLUST_NUM + 1
			    THIS = NINT(TOT_CLUST_NUM)
              endif
              
			RANK1 = NUM_RANK(I + 1)
			TREE(RANK1, I + 1) = THIS
			CLUST_STOR(THIS, 1:I) = CLUST_STOR(IDM(J), 1:I)
			! 完整性检查
              do	FKY = 1, I
                  ! CLUST_STOR(THIS, I).EQ.TEMP_NEW(J) 避免重复核子 CHEN SM
				if(CLUST_STOR(THIS, FKY).EQ.0.OR.
     &                CLUST_STOR(THIS, I).EQ.TEMP_NEW(J))then
					TOT_CLUST_NUM = 0
					NUM_RANK(1:SURF_MAX_A) = 0.
					print*, 'FAULT !'
					goto 1131
                  endif
			enddo
			CLUST_STOR(THIS, I + 1) = TEMP_NEW(J)
              CLUST_STOR(THIS, I + 1) = TEMP_NEW(J)
		enddo
	enddo
1131	continue

	END SUBROUTINE CONSTR_LOOP

c///////////////////////////////////////////////////////////////////////////////
c     輸入leading核子的時候，INPUT_CLUST_ID = 0, TREE(0) = (0, 0, 0)
c     （也就是2-4下標對應的元素都為0）CLUST_STOR(0,1:4) = (I_L, 0, 0, 0)（祗有1號元素不為0）
	RECURSIVE SUBROUTINE CONSTR_AND_TUNNEL(INPUT_CLUST_ID)
c///////////////////////////////////////////////////////////////////////////////
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
	INTEGER, PARAMETER :: PICK_MAX = 233 ! 最大候选核子数，上限应该是总核子数-1
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),
     &			TOT_CLUST_NUM, He4_ERROR
	common/SURF/TOO_MANY

	INTEGER, INTENT(IN) :: INPUT_CLUST_ID
	INTEGER EXTRA(PICK_MAX), PICKED, RANK, COUNT_PICK, TOT_CLUST_NUM1
	INTEGER	IN_CLUST(SURF_MAX_A - 1), I_TEMP

	if(CLUST_STOR(INPUT_CLUST_ID, 1).EQ.0)	
     &    print*, 'EMPTY INPUT TO CLUST_CONSTR!'

	! 计算输入簇的核子数（秩）
	RANK = 0
	do I = 1, SURF_MAX_A
		if(CLUST_STOR(INPUT_CLUST_ID, I).NE.0)	RANK = RANK + 1
      enddo
      
      ! 检查是否已达到最大核子数
	if(RANK.EQ.SURF_MAX_A)then
		print*, 'SURF_MAX_A INPUT TO CONSTR_AND_TUNNEL'
		He4_ERROR = 1.
		return
	endif

	! 輸入為leading核子（单个核子）的時候1-1000號cluster的記憶全部清零
	! CLUSTER總數TOT_CLUST_NUM = 0，祗有0號CLUSTER
	if(RANK.EQ.1)then
		do I = 1, STORE_MAX
			TREE(I,2:SURF_MAX_A) = 0 ! 这里的SURF_MAX_A本来为4
			CLUST_STOR(I,1:SURF_MAX_A) = 0
		enddo
		TOT_CLUST_NUM = 0 ! 重置簇计数器
		TOO_MANY = 0      ! 重置溢出标志
	endif

	if(NINT(TOO_MANY).EQ.1)return ! 如果已溢出，直接返回

	EXTRA(1:PICK_MAX) = 0 ! 存储候选核子
	IN_CLUST(1:(SURF_MAX_A - 1)) = ! 准备输入簇信息，用于调用PICK_NUCLEON
     &    CLUST_STOR(INPUT_CLUST_ID, 1:(SURF_MAX_A - 1)) ! CLUST_STOR(0,1)=I_L
	call PICK_NUCLEON(EXTRA, IN_CLUST, PICKED) ! 获取符合条件的候选核子

	TOT_CLUST_NUM1 = NINT(TOT_CLUST_NUM)

	COUNT_PICK = 1
	if(PICKED.GE.1)then
		do I = TOT_CLUST_NUM1 + 1, TOT_CLUST_NUM1 + PICKED
			! 检查团簇存储空间是否已满
              if(I.EQ.(STORE_MAX - 1))then
				TOO_MANY = 1.
				I_TEMP = I - 1
				EXIT
			endif
			! 儲存結團結構
			do J = 1, RANK	! SURF_MAX_A - 1
				if(CLUST_STOR(INPUT_CLUST_ID, J).EQ.0)cycle
				CLUST_STOR(I, J) = CLUST_STOR(INPUT_CLUST_ID, J)
			enddo
			CLUST_STOR(I, RANK + 1) = EXTRA(COUNT_PICK)
			! 儲存結團方位
			do J = 1, RANK
				if(RANK.EQ.1)EXIT
				if(J.EQ.1)cycle
				TREE(I, J) = TREE(INPUT_CLUST_ID, J)
			enddo
			TREE(I, RANK + 1) = COUNT_PICK
			COUNT_PICK = COUNT_PICK + 1
		enddo
	else ! 沒有抓到核子
		return
	endif

	TOT_CLUST_NUM = TOT_CLUST_NUM + PICKED

	if(NINT(TOO_MANY).EQ.1)then
		TOT_CLUST_NUM = I_TEMP
		return
	endif

	! 已經是獲取的已經是SURF_MAX_A個核子
	if((RANK + 1).EQ.SURF_MAX_A)then
		return
	endif

	do I = TOT_CLUST_NUM1 + 1, TOT_CLUST_NUM1 + PICKED
		call CONSTR_AND_TUNNEL(I)
		if(He4_ERROR.NE.0)EXIT
	enddo
	
	return
	END	SUBROUTINE	CONSTR_AND_TUNNEL



c///////////////////////////////////////////////////////////////////////////////
c     輸入的IN_CLUST(3)中，前A_IN_CL個元素不為零，後面全部設零
	SUBROUTINE	PICK_NUCLEON(PICK_NUCLEON1, IN_CLUST, PICKED)
c///////////////////////////////////////////////////////////////////////////////
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
	INTEGER, PARAMETER :: PICK_MAX = 233
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
      common/bb/np,nzp,ntag,x00 
	common/SURF/ VANISH(500), H0_SURF(5)
	common/SURF1/CLUST_LEAST_DIST, XCC_TARG, YCC_TARG, ZCC_TARG
     &,R0_SURF, D_SURF, COL_PATH, NDATA
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)


	
	INTEGER, INTENT(OUT) :: PICK_NUCLEON1(PICK_MAX)	!為了到時候便於對比，核子編號應該從小到大
													!對於d, t, 3He，第四元素表示質子數
	INTEGER, INTENT(OUT) :: PICKED
	INTEGER, INTENT(IN)	 :: IN_CLUST(SURF_MAX_A - 1)

	INTEGER	INIT_COUNT(2), A_IN_CL, CASE_NO
	REAL	H0_P, H0_N, R_IM, P_IM, CLUST_DIST_X
	REAL	CLUST_DIST_Y, CLUST_DIST_Z, CLUST_DIST
	REAL	FACT

	INIT_COUNT(1:2) = 0
	CASE_NO = 0
	H0_P = 0.
	H0_N = 0.
	PICKED = 0

	PICK_NUCLEON1(1:PICK_MAX) = 0

	A_IN_CL = 0
	XCC = 0.
	YCC = 0.
	ZCC = 0.
	PXCC = 0.
	PYCC = 0.
	PZCC = 0.
	do I = 1, SURF_MAX_A - 1

		if(IN_CLUST(I).EQ.0)cycle
		A_IN_CL = A_IN_CL + 1 ! 结团核子数
		if(ISY(IN_CLUST(I)).EQ.-1)then
			INIT_COUNT(1) = INIT_COUNT(1) + 1 ! 初始質子數
		elseif(ISY(IN_CLUST(I)).EQ.1)then
			INIT_COUNT(2) = INIT_COUNT(2) + 1 ! 初始中子數
		endif
		XCC = XCC + X(IN_CLUST(I))
		YCC = YCC + Y(IN_CLUST(I))
		ZCC = ZCC + Z(IN_CLUST(I))
		PXCC = PXCC + PX(IN_CLUST(I))
		PYCC = PYCC + PY(IN_CLUST(I))
		PZCC = PZCC + PZ(IN_CLUST(I))

      enddo
      
      ! 计算质心坐标和动量
	XCC = XCC/FLOAT(A_IN_CL)
	YCC = YCC/FLOAT(A_IN_CL)
	ZCC = ZCC/FLOAT(A_IN_CL)
	PXCC = PXCC/FLOAT(A_IN_CL)
	PYCC = PYCC/FLOAT(A_IN_CL)
	PZCC = PZCC/FLOAT(A_IN_CL)
	
	! 基于 Hui-Gan Cheng and Zhao-Qing Feng 2021 ChinesePhys.C 45 084107 的模型
	! 设置相空间约束参数 H0
	if(INT(INCL4_6).EQ.0)then
          ! INIT_COUNT(1)质子数、INIT_COUNT(2)中子数、A_IN_CL核子数
		if(INIT_COUNT(1).EQ.1.AND.A_IN_CL.EQ.1)then
			CASE_NO = 1		! p
		elseif(INIT_COUNT(1).EQ.0.AND.A_IN_CL.EQ.1)then
			CASE_NO = 2		! n
		elseif(INIT_COUNT(1).EQ.1.AND.A_IN_CL.EQ.2)then
			CASE_NO = 3		! d
		elseif(INIT_COUNT(1).EQ.1.AND.A_IN_CL.EQ.3)then
			CASE_NO = 4		! t
		elseif(INIT_COUNT(1).EQ.2.AND.A_IN_CL.EQ.3)then
			CASE_NO = 5		! 3He
		endif

		SELECT CASE(CASE_NO)
			CASE(1)
				H0_N = H0_SURF(1)
			CASE(2)
				H0_P = H0_SURF(1)
			CASE(3)
				H0_N = H0_SURF(2)
				H0_P = H0_SURF(3)
			CASE(4)
				H0_P = H0_SURF(4)
			CASE(5)
				H0_N = H0_SURF(5)
		end	SELECT 
	endif

      ! 主循环：遍历所有核子，寻找符合条件的候选核子
	do I = 1, NP
		CLUST_ENERGY = 0. ! 结团能量（此处未使用）
		
          ! 排除条件
          do J = 1, SURF_MAX_A - 1 ! 母结团中的核子
              if(IN_CLUST(J).EQ.I)goto 6660 
          enddo

		if(IDR(I).GT.0)     cycle  ! 其他重子
		if(VANISH(I).EQ.1)  cycle  ! 已发射核子
          
		if(INT(INCL4_6).EQ.0.AND.A_IN_CL.LE.3)then ! ORIGINAL
			if(A_IN_CL.EQ.1)then ! 在輸入為leading核子的時候衹產生d
				if((ISY(IN_CLUST(1))*ISY(I)).EQ.1)cycle
			elseif(CASE_NO.EQ.4.AND.ISY(I).EQ.1)then ! 輸入t核
				H0_N = 300. 
			elseif(CASE_NO.EQ.5.AND.ISY(I).EQ.-1)then ! 輸入3He核
				H0_P = 300.
			endif
		elseif(A_IN_CL.LE.3)then ! INCL4.6 模型
			if(A_IN_CL.EQ.1)then
				if((ISY(IN_CLUST(1))*ISY(I)).EQ.1)cycle
				H0_N = 424.
                  H0_P = H0_N
			elseif(A_IN_CL.EQ.2)then
				H0_N = 300.
                  H0_P = H0_N
			elseif(A_IN_CL.EQ.3)then
				H0_N = 300.
                  H0_P = H0_N
			endif
		else ! 对于 A_IN_CL >= 4 的大核簇
			if(ISY(I).EQ.-1)then
				if(INIT_COUNT(2).EQ.1)then
					if(INIT_COUNT(1).EQ.5)cycle
				endif
			else
				if(INIT_COUNT(2).EQ.6)then
					if(INIT_COUNT(1).EQ.1)cycle
				endif
			endif
			H0_N = H1_SURF*((A_IN_CL + 1)/5.)**0.333
			H0_P = H0_N
          endif
          
          ! 条件1：检查核子与靶核的距离约束
		CLUST_DIST = (X(I) - XCC_TARG)**2 +
     &				 (Y(I) - YCC_TARG)**2 +
     &				 (Z(I) - ZCC_TARG)**2
		CLUST_DIST = SQRT(CLUST_DIST)
		CLUST_DIST = CLUST_DIST - R0_SURF - D_SURF - COL_PATH
		if(CLUST_DIST.GT.0.)cycle ! 距离太远，跳过

		! 条件2：检查与初始核簇的相对距离
		R_IM = (X(I) - XCC)**2 + (Y(I) - YCC)**2 +
     &			(Z(I) - ZCC)**2
		R_IM = SQRT(R_IM)
		if(R_IM.LT.1.)cycle ! 距离太近，跳过
          
          ! 条件3：检查相对动量（相空间约束的核心）
		FACT = A_IN_CL + 1
		if(NDATA.EQ.1)then
			PXXCC = PXCC*A_IN_CL/FACT ! 合并后的平均每核子动量
			PYYCC = PYCC*A_IN_CL/FACT
			PZZCC = PZCC*A_IN_CL/FACT
			P_IM = (PX(I)*A_IN_CL/FACT - PXXCC)**2 + ! PX(I)*A_IN_CL/FACT 为质心系动量
     &			   (PY(I)*A_IN_CL/FACT - PYYCC)**2 +
     &			   (PZ(I)*A_IN_CL/FACT - PZZCC)**2
		else
			P_IM = (PX(I) - PXCC)**2 + 
     &               (PY(I) - PYCC)**2 +
     &               (PZ(I) - PZCC)**2
		endif
		P_IM = SQRT(P_IM)
		P_IM = 1000.*P_IM ! (MeV)
		if(INT(INCL4_6).EQ.0)then
			if(ISY(I).EQ.-1)then ! 質子
				if((R_IM*P_IM).GE.H0_P)cycle
			elseif(ISY(I).EQ.1)then	! 中子
				if((R_IM*P_IM).GE.H0_N)cycle
			else
				cycle
			endif
		else		
			if((R_IM*P_IM).GE.H0_N)cycle ! 相空间约束条件
          endif
          
          ! 条件4：检查新形成的核簇与靶核的最小距离
		CLUST_DIST_X = (X(I) + XCC)/(A_IN_CL+1)
		CLUST_DIST_Y = (Y(I) + YCC)/(A_IN_CL+1)
		CLUST_DIST_Z = (Z(I) + ZCC)/(A_IN_CL+1)

		CLUST_DIST = (CLUST_DIST_X - XCC_TARG)**2 +
     &				 (CLUST_DIST_Y - YCC_TARG)**2 +
     &				 (CLUST_DIST_Z - ZCC_TARG)**2
		CLUST_DIST = SQRT(CLUST_DIST)
		if(CLUST_DIST.LT.CLUST_LEAST_DIST)cycle ! 距离太小，跳过
          
          ! 通过所有条件，记录候选核子
		PICKED = PICKED + 1
		PICK_NUCLEON1(PICKED) = I ! 记录核子编号
6660  continue
	enddo

	return
	END SUBROUTINE PICK_NUCLEON



	FUNCTION	EMISSION_TEST(Z_C, A_C, INPUT_ID, TOT_E)
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)              
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),TOT_CLUST_NUM, 
     &			He4_ERROR
	common/TUNNEL/ Z_RESID, A_RESID
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)
      
      INTEGER, INTENT(IN) :: Z_C, A_C
	INTEGER, INTENT(IN)	:: INPUT_ID
	INTEGER	RESULT_TEST, IIK
	REAL	BIN_E, R_CAL, II
	REAL, INTENT(OUT)	::	TOT_E
	INTEGER	PROTON_NUM, NEUTRON_NUM
	INTEGER	Z_CLUST, A_CLUST
	REAL(KIND = 8)	E_TOT

	RESULT_TEST = 1

	TOT_E = 0.
      ! 计算当前结团质子数、中子数、核子数、总能量
	do I = 1, SURF_MAX_A
          IIK = CLUST_STOR(INPUT_ID, I)
		if(IIK.EQ.0)cycle
		TOT_E = TOT_E + GET_TOTAL_ENERGY(IIK)
	enddo
	NEUTRON_NUM = A_C - Z_C
	BIN_E = 7.289*Z_C + 8.071*(NEUTRON_NUM) -
     &    READ_XBIN(Z_C, A_C) ! 结合能
	TOT_E = TOT_E + BIN_E
     
	if(TOT_E.LE.0.)then
		EMISSION_TEST = 0
		return
	endif

	E_TOT = TOT_E
	EMISSION_TEST = 0
	call TUNNEL1(Z_C, A_C, E_TOT, RESULT_TEST) ! Z_RESID, A_RESID异常，不应为0
	EMISSION_TEST = RESULT_TEST

	return
	end


	SUBROUTINE	TUNNEL1(Z_C, A_C, E_TOT, TUNN)
	REAL(KIND = 8), PARAMETER :: EE = 1.44, H_BAR = 6.58211889, 
     &	M_0 = 938., C = 2.99792458, PI = 3.141592654, 
     &    D_F = 2.4/4./LOG(3.)
	COMPLEX(KIND = 8), PARAMETER :: CJ = (0., 1.)
      INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
      
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)              
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),TOT_CLUST_NUM, 
     &			He4_ERROR
	common/TUNNEL/ Z_RESID, A_RESID
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)
      
	INTEGER	N, I, I_ARR, FLAG_DV, I_MOD, I_MAX
	INTEGER, INTENT(IN) :: Z_C, A_C
	REAL(KIND = 8), INTENT(IN) :: E_TOT
	REAL(KIND = 8)	V_BOT, T, Z_RR, Z_RRR, RHO_C, R_R, XX, RM
	REAL(KIND = 8)	V_C, V_WS, V_TOT, DV_C, DV_WS, EXP_WS
	REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: V(:), AA(:)
	COMPLEX(KIND = 8)	A1, B1, A1_, B1_, BETA, BETA_, E, E_, B, B_
	COMPLEX(KIND = 8)	R, R1
	REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: TEMP_V1(:),
     &				TEMP_A1(:), TEMP_V2(:), TEMP_A2(:)
	INTEGER			V_MOD(20), E_ID, J, SORT3, SORT2(0:20)
	REAL(KIND = 8)	SORT(0:20), SORT1, RAN0
	LOGICAL	FLAG_PRINT
	INTEGER, INTENT(OUT) :: TUNN	

	open(1, FILE = 'Sqr_pot.dat', STATUS = 'UNKNOWN') ! output

	I_MAX = 0	!為了防止死循環，規定每個循環語句最多執行1萬次；

	if(Z_C.EQ.0)	then
		goto 45682
		if(E_TOT.LE.0.001)	then
			T = 0.
			goto	12432
		endif
	endif
!1: 划分势垒---------------------------------------------------------------------------------------
	if(E_TOT.LE.1.)	then
		T = 0.
		goto	12432
	endif

	I_MAX = 0
	XX   = 0.
	Z_RR = 0.
	R_R  = 1.15*A_RESID**0.333333
	do
		if(I_MAX.GE.10000)	EXIT
		if(XX.GE.2.*R_R)	EXIT
		RHO_C = D_F
		RHO_C = (XX - R_R)/RHO_C
		RHO_C = 1./(1. + EXP(RHO_C))
		Z_RR = Z_RR + 0.1*RHO_C*XX**2
		XX = XX + 0.1
		I_MAX = I_MAX + 1
	enddo
	
	XX    = 0.1
	RM    = 0.1
	Z_RRR = 0.
	I = 1
	I_ARR = 1
	FLAG_DV = 0
	I_MOD = 0
	FLAG_PRINT = .TRUE.
	ALLOCATE(TEMP_V1(100))
	ALLOCATE(TEMP_A1(100))
	I_MAX = 0

	do

		if(I_MAX.GE.10000)	EXIT
		if(FLAG_DV.EQ.1.AND.V_TOT.LT.(E_TOT/10.))	EXIT

		if(XX.LT.2.*R_R)	then
			RHO_C = D_F
			RHO_C = (XX - R_R)/RHO_C
			RHO_C = 1./(1. + EXP(RHO_C))
			Z_RRR = Z_RRR + 0.1*RHO_C*XX**2
			
			V_WS  = (XX - R_R)/0.5
			EXP_WS = EXP(V_WS)
			V_WS  = 1. + EXP_WS
			V_WS  = -50.*A_C/V_WS
		endif

		V_C  = EE*Z_C*(Z_RESID*Z_RRR/Z_RR)/XX
		if(XX.GE.2.*R_R)	RHO_C = 0.
		DV_C = -Z_RRR/XX**2 + XX*RHO_C
		DV_C = DV_C*Z_C*Z_RESID*EE/Z_RR

		DV_WS = EXP_WS/(1. + EXP_WS)**2.
		if(XX.GE.2.*R_R)	DV_WS = 0.
		DV_WS = DV_WS*100.*A_C

		V_TOT = V_C + V_WS
		
		DV_C = DV_C + DV_WS
		if(DV_C.LT.0.)	then
			FLAG_DV = 1		!过了库伦位垒
			if(FLAG_PRINT)	then
C 				print*, 'Coulomb Barrier =', V_TOT
				FLAG_PRINT = .FALSE.
			endif
		endif

		if(NINT(10*XX).EQ.NINT(10*RM))	then

			if(I_ARR.EQ.1)		then
				TEMP_A1(I) = 0.1*NINT(10./(ABS(DV_C) + 1.D-3)**0.15)
				if(FLAG_DV.EQ.0.AND.TEMP_A1(I).GE.2.) TEMP_A1(I) = 2.	!库伦位垒之前应该精细
				TEMP_V1(I) = V_TOT + 0.5*DV_C*TEMP_A1(I)
				if(ABS(TEMP_V1(I) - E_TOT).LT.0.005)	then
					I_MOD = I_MOD + 1
					if(I_MOD.GT.20)	print*, 'V_MOD EXCEEDED!'
					V_MOD(I_MOD) = I
				endif
				WRITE(1,*)	I, RM, TEMP_V1(I) ! 检测是否异常，TEMP_V1异常，可能源于Z_R=0
				RM = RM + TEMP_A1(I)
			elseif(I_ARR.EQ.-1)	then
				TEMP_A2(I) = 0.1*NINT(10./(ABS(DV_C) + 1.D-3)**0.15)
				if(FLAG_DV.EQ.0.AND.TEMP_A2(I).GE.2.) TEMP_A2(I) = 2.	!库伦位垒之前应该精细
				TEMP_V2(I) = V_TOT + 0.5*DV_C*TEMP_A2(I)
				if(ABS(TEMP_V2(I) - E_TOT).LT.0.005)	then
					I_MOD = I_MOD + 1
					if(I_MOD.GT.20)	print*, 'V_MOD EXCEEDED!'
					V_MOD(I_MOD) = I
				endif
				WRITE(1,*) I, RM, TEMP_V2(I)
				RM = RM + TEMP_A2(I)
			endif

			if(MOD(I,100).EQ.0)	then
				I_ARR = -1.*I_ARR
				if(I_ARR.EQ.1)	then
					ALLOCATE(TEMP_V1(I + 100))
					ALLOCATE(TEMP_A1(I + 100))
					TEMP_V1(1:I) = TEMP_V2(1:I)
					TEMP_A1(1:I) = TEMP_A2(1:I)
					DEALLOCATE(TEMP_V2)
					DEALLOCATE(TEMP_A2)
				else
					ALLOCATE(TEMP_V2(I + 100))
					ALLOCATE(TEMP_A2(I + 100))
					TEMP_V2(1:I) = TEMP_V1(1:I)
					TEMP_A2(1:I) = TEMP_A1(1:I)
					DEALLOCATE(TEMP_V1)
					DEALLOCATE(TEMP_A1)
				endif
			endif

			I = I + 1

		endif

		XX   = XX + 0.1

		I_MAX = I_MAX + 1

	enddo

	N = I - 1

	ALLOCATE(V(N))
	ALLOCATE(AA(N))

	if(I_ARR.EQ.1)	then
		V(1:N)  = TEMP_V1(1:N)
		AA(1:N) = TEMP_A1(1:N)
		DEALLOCATE(TEMP_V1)
		DEALLOCATE(TEMP_A1)
	else
		V(1:N)  = TEMP_V2(1:N)
		AA(1:N) = TEMP_A2(1:N)
		DEALLOCATE(TEMP_V2)
		DEALLOCATE(TEMP_A2)	
	endif

	I_MAX = 0
	do	I = 1, I_MOD
		if(I_MAX.GE.10000)	EXIT
		if(V(V_MOD(I)).LE.E_TOT)	V(V_MOD(I)) = V(V_MOD(I)) - 0.01
		if(V(V_MOD(I)).GT.E_TOT)	V(V_MOD(I)) = V(V_MOD(I)) + 0.01
		I_MAX = I_MAX + 1
	enddo
	
45682	continue
	!出射粒子是中子
	if(Z_C.EQ.0)	then
		ALLOCATE(V(2))
		ALLOCATE(AA(2))
		N     = 2

		V(1)  = 0.
		AA(1) = 1.
		V(2)  = 0.
		AA(2) = 1.
	endif
	
!2: 求穿透几率-------------------------------------------------------------------------------------

	V_BOT = -50.*A_C

	BETA_ = CJ*SQRT(2*A_C*M_0*E_TOT)/(H_BAR*C*10)
	BETA  = ABS(E_TOT - V(N))
	BETA  = SQRT(2*A_C*M_0*BETA)/(H_BAR*C*10)

	if(E_TOT.GT.V(N))	then
		BETA = BETA*CJ
	endif

	!第N个势垒的A和B
	A1_ = 1. + BETA_/BETA
	B1_ = 1. - BETA_/BETA

	I_MAX = 0
	BETA_ = BETA
	I = N - 1
	do

		if(I_MAX.GE.10000)	EXIT
		BETA = ABS(E_TOT - V(I))
		BETA = SQRT(2*A_C*M_0*BETA)/(H_BAR*C*10)
		if(E_TOT.GT.V(I))	BETA = CJ*BETA
		
		E  = EXP(-AA(I+1)*BETA_)
		E_ = EXP( AA(I+1)*BETA_)
		B  = +1.*BETA_*EXP(-AA(I+1)*BETA_)
		B_ = -1.*BETA_*EXP( AA(I+1)*BETA_)

		A1 = (A1_*(E + B/BETA) + B1_*(E_ + B_/BETA))/2.
		B1 = (A1_*(E - B/BETA) + B1_*(E_ - B_/BETA))/2.

		A1_ = A1
		B1_ = B1
		BETA_ =BETA

		I = I - 1
		if(I.EQ.0)	EXIT	!已经求出了I = 1
		I_MAX = I_MAX + 1
	enddo


	BETA = E_TOT - V_BOT
	BETA = SQRT(2*A_C*M_0*BETA)/(H_BAR*C*10.)	!BETA_0
	BETA = CJ*BETA

	E  = EXP(-AA(1)*BETA_)
	E_ = EXP( AA(1)*BETA_)
	B  = +1.*BETA_*EXP(-AA(1)*BETA_)
	B_ = -1.*BETA_*EXP( AA(1)*BETA_)

	R  = A1_*(E*BETA - B) + B1_*(E_*BETA - B_)
	R1 = A1_*(E*BETA + B) + B1_*(E_*BETA + B_)
	R  = R/R1
					
	DEALLOCATE(V)
	DEALLOCATE(AA)

	T = 1. - (ABS(R))**2
	if(T.LT.0.)	T = 0.

12432	continue

	call RANDOM_NUMBER(RAN0)	
	TUNN = 0
	if(RAN0.LE.T) TUNN = 1

C 	print*, E_TOT, RAN0, T

	return
	end

	!這個子程序的作用是：在構造結團的時候因為是遞歸構造的，最後同一級的結團上
	!可能有一些結團是由完全相同的核子構成的，衹不過這些結團的構造經歷了不同的
	!路徑，所以有必要將相同的結團化成一個；
	SUBROUTINE	SORT_CANDID(RANK, CLUST_NUM, NUM_CANDID)
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,             
     &       px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &       idr(500),xms(500),yms(500),zms(500),pxms(500),pyms(500),
     &       pzms(500),emson(500),idms(500),icp(500),idanti(500)
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),TOT_CLUST_NUM
	common/EMIT_CANDID/EMIT_CANDID(STORE_MAX), 
     &		CLUST_BIN(2:SURF_MAX_A, STORE_MAX)
	INTEGER, INTENT(IN) :: RANK
	INTEGER, INTENT(IN) :: CLUST_NUM
	INTEGER, INTENT(OUT):: NUM_CANDID
	INTEGER	I, J, K, SUM1, SUM2, K1, K2, K3, JJJ


	do I = 1, CLUST_NUM - 1
		if(CLUST_BIN(RANK, I).EQ.0.)	cycle
		do J = I + 1, CLUST_NUM
			K1 = 0
			K2 = 0
			K3 = 0
			if(INT(CLUST_BIN(RANK, J)).EQ.0)	cycle
			do K = 2, RANK
				K1 = K1 + CLUST_STOR(INT(CLUST_BIN(RANK, I)), K) -
     &					CLUST_STOR(INT(CLUST_BIN(RANK, J)), K)
				K2 = K2 + CLUST_STOR(INT(CLUST_BIN(RANK, I)), K)**2 - 
     &					CLUST_STOR(INT(CLUST_BIN(RANK, J)), K)**2
				K3 = K3 + CLUST_STOR(INT(CLUST_BIN(RANK, I)), K)**3 - 
     &					CLUST_STOR(INT(CLUST_BIN(RANK, J)), K)**3
			enddo
			if(K1.EQ.0.AND.K2.EQ.0.AND.K3.EQ.0)	then	!I和J所裝的結團的核子構成 完全相同的話，則取消J
				CLUST_BIN(RANK, J) = 0.
			endif
		enddo
	enddo

	J = 0
	do	I = 1, CLUST_NUM

		if(CLUST_BIN(RANK, I).NE.0.)	then

			if(RANK.EQ.5)					then	!排除掉5Be
				JJJ = 0
				do	JJ = 1, RANK
					if(ISY(CLUST_STOR(INT(CLUST_BIN(RANK, I)), JJ)).
     &				EQ.-1)	JJJ = JJJ + 1
				enddo
				if(JJJ.EQ.4)	cycle
			endif

			J = J + 1
			EMIT_CANDID(J) = CLUST_BIN(RANK, I)

		endif

	enddo

	NUM_CANDID = J

	return
	end

	SUBROUTINE	GET_CLUST_KIN(CLUST_ID, E_CL1, MOM, AGM)
	INTEGER, PARAMETER :: surf_max_a = 4
	INTEGER, PARAMETER :: STORE_MAX = 1000
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,        
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500),
     &  idr(500),xms(500),yms(500),zms(500),pxms(500),
     &  pyms(500),pzms(500),emson(500),idms(500),icp(500),idanti(500)
	common/SURF/VANISH(500),H0_SURF(5),TREE(0:STORE_MAX,1:SURF_MAX_A),
     &			CLUST_STOR(0:STORE_MAX,SURF_MAX_A),TOT_CLUST_NUM, 
     &			He4_ERROR

	REAL, PARAMETER		:: ANGMCOEF=5.067731
	REAL, INTENT(OUT)	:: AGM(3), MOM(3)
	REAL, INTENT(IN)	:: E_CL1
	REAL	E_CL
	INTEGER, INTENT(IN)	:: CLUST_ID
	INTEGER	RANK, I, J, I00
	INTEGER	C(SURF_MAX_A)
	REAL	XC, YC, ZC
	
	C(1:SURF_MAX_A) = 0
	RANK = 1
	J = CLUST_ID

	do	I = 2, SURF_MAX_A
		I00 = SURF_MAX_A - I + 2
		if(CLUST_STOR(J,I00).NE.0)then
			RANK = I00
			EXIT
		endif
	enddo

	do I = 1, SURF_MAX_A
		C(I) = CLUST_STOR(J, I)
	enddo

	XC = 0.
	YC = 0.
	ZC = 0.
	MOM(1:3) = 0.
	do I = 1, RANK
		MOM(1) = MOM(1) + PX(C(I))
		MOM(2) = MOM(2) + PY(C(I))
		MOM(3) = MOM(3) + PZ(C(I))

		XC = XC + X(C(I))
		YC = YC + Y(C(I))
		ZC = ZC + Z(C(I))
	enddo
	XC = XC/FLOAT(RANK)
	YC = YC/FLOAT(RANK)
	ZC = ZC/FLOAT(RANK)

	AGM(1) = YC*MOM(3) - ZC*MOM(2)		!C.M. FRAME
	AGM(2) = ZC*MOM(1) - XC*MOM(3)
	AGM(3) = XC*MOM(2) - YC*MOM(1)
	AGM(1) = AGM(1)*ANGMCOEF
	AGM(2) = AGM(2)*ANGMCOEF
	AGM(3) = AGM(3)*ANGMCOEF

	E_CL = E_CL1
	E_CL = 0.001*E_CL
	E_CL = 2*0.938*RANK*E_CL
	E_CL = SQRT(E_CL)
	E_CL = E_CL/SQRT(MOM(1)**2 + MOM(2)**2 + MOM(3)**2)

	do I = 1, 3
		MOM(I) = MOM(I)*E_CL
	enddo

	return
	end


	SUBROUTINE	CLUSTER_STRONG_DECAY(Z, A, M_PRIM, 
     &			Z_LAST, A_LAST, M_LAST, ST_P, ST_M, NP)
	IMPLICIT NONE
	INTEGER, PARAMETER :: MAX_EMIT = 4

	INTEGER, INTENT(IN) :: Z, A
	REAL, INTENT(IN) ::	M_PRIM(3)
	REAL, INTENT(OUT) :: M_LAST(3), ST_M(MAX_EMIT,3)
	INTEGER, INTENT(OUT) ::	A_LAST, Z_LAST, ST_P(MAX_EMIT), NP
	REAL	V(3), E_M, M_3BOD(3, 3), E_R, E_0
	INTEGER	DEC, ZZ, AA
	LOGICAL I


	if((A - Z).EQ.1)	then
		if(Z.EQ.4)	print*, '5BE INPUT!!'
		if(Z.EQ.6)	print*, '7C INPUT!!'
		if(Z.EQ.7)	print*, '8N INPUT!!'
	elseif((A - Z).EQ.7)then
		print*, '8H INPUT!!'
	endif

	E_M = SQRT((A*0.938)**2 + M_PRIM(1)**2 + M_PRIM(2)**2 + 
     &		M_PRIM(3)**2)
	V(1:3) = M_PRIM(1:3)/E_M

	!0	中子
	!1	質子
	!2	alpha

	Z_LAST = Z
	A_LAST = A
	ZZ = Z
	AA = A
	M_LAST(1:3) = M_PRIM(1:3)


	I = .FALSE.
	NP = 0
	do

	I = .TRUE.

	call STRONG_DECAY(ZZ, AA, DEC, M_3BOD, A_LAST, Z_LAST)

	if	  (DEC.EQ.0)	then
		EXIT
	elseif(DEC.EQ.1.OR.DEC.EQ.3)	then
		NP = NP + 1
		if(DEC.EQ.1)	then
			ST_P(NP) = 1
		else
			ST_P(NP) = 0
		endif
		E_R = SQRT(0.938**2 + M_3BOD(2, 1)**2 + M_3BOD(2, 2)**2 +
     &			M_3BOD(2, 3)**2)
		call LORENTZ(-V(1),-V(2),-V(3),E_R,M_3BOD(2, 1),M_3BOD(2, 2)
     &	,M_3BOD(2, 3),E_0,ST_M(NP, 1),ST_M(NP, 2),ST_M(NP, 3))
	elseif(DEC.EQ.2.OR.DEC.EQ.4)	then
		NP = NP + 1
		if(DEC.EQ.2)	then
			ST_P(NP) = 1
		else
			ST_P(NP) = 0
		endif
		E_R = SQRT(0.938**2 + M_3BOD(2, 1)**2 + M_3BOD(2, 2)**2 +
     &			M_3BOD(2, 3)**2)
		call LORENTZ(-V(1),-V(2),-V(3),E_R,M_3BOD(2, 1),M_3BOD(2, 2)
     &	,M_3BOD(2, 3),E_0,ST_M(NP, 1),ST_M(NP, 2),ST_M(NP, 3))
		NP = NP + 1
		if(DEC.EQ.2)	then
			ST_P(NP) = 1
		else
			ST_P(NP) = 0
		endif
		E_R = SQRT(0.938**2 + M_3BOD(3, 1)**2 + M_3BOD(3, 2)**2 +
     &			M_3BOD(3, 3)**2)
		call LORENTZ(-V(1),-V(2),-V(3),E_R,M_3BOD(3, 1),M_3BOD(3, 2)
     &	,M_3BOD(3, 3),E_0,ST_M(NP, 1),ST_M(NP, 2),ST_M(NP, 3))
	elseif(DEC.EQ.5)	then
		NP = NP + 1
		ST_P(NP) = 2
		E_R = SQRT((4*0.938)**2 + M_3BOD(2, 1)**2 + M_3BOD(2, 2)**2 +
     &			M_3BOD(2, 3)**2)
		call LORENTZ(-V(1),-V(2),-V(3),E_R,M_3BOD(2, 1),M_3BOD(2, 2)
     &	,M_3BOD(2, 3),E_0,ST_M(NP, 1),ST_M(NP, 2),ST_M(NP, 3))
	endif

	if(DEC.NE.0)	then

		ZZ = Z_LAST
		AA = A_LAST
		E_R = SQRT((A_LAST*0.938)**2 + M_3BOD(1, 1)**2 + 
     &			M_3BOD(1, 2)**2 + M_3BOD(1, 3)**2)
		call LORENTZ(-V(1),-V(2),-V(3),E_R,M_3BOD(1, 1),M_3BOD(1, 2)
     &		,M_3BOD(1, 3),E_0,M_LAST(1),M_LAST(2),M_LAST(3))

		E_M = (A_LAST*0.938)**2 + M_LAST(1)**2 + M_LAST(2)**2 
     &			+ M_LAST(3)**2
		E_M = SQRT(E_M)
		V(1:3) = M_LAST(1:3)/E_M

	endif

	enddo


	end


	SUBROUTINE	STRONG_DECAY(Z, A, DEC_MOD, M_3BOD, A1, Z1)

	REAL, PARAMETER :: PI = 3.14159
	REAL	Q, RAN
	REAL	DELTA, XBIN(0:1000,0:1000)
	INTEGER, INTENT(IN)	::  Z, A
	INTEGER SWITCH, I, IZ, IA
	REAL	I1, Be, RR
	REAL	P, THETA, PHI
	CHARACTER*4 ELS
	REAL, INTENT(OUT) :: M_3BOD(3, 3)
	REAL	P1X, P1Y, P1Z
	REAL	PHYPX, PHYPY, PHYPZ, EPP2, EPP1
	REAL	PPX, PPY, PPZ, EP2, EP1
	REAL	PXX, PYY, PZZ, PA, PA2, S
	REAL	XX, YY, ZZ, BETX, BETY, BETZ
	REAL	PK0, PKX, PKY, PKZ
	DOUBLE PRECISION EM0, EM1, EM2, EM3, PPP
	INTEGER	NN1, ZZ1
	INTEGER,INTENT(OUT) ::	DEC_MOD, A1, Z1
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)

	NN1 = A - Z
	ZZ1 = Z

	if(NN1.EQ.1)	then
		if(ZZ1.EQ.3)	then
			DEC_MOD = 1
		elseif(ZZ1.EQ.5)then
			DEC_MOD = 2
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return
		endif
	elseif(NN1.EQ.2)then

		if(ZZ1.EQ.3)	then
			DEC_MOD = 1
		elseif(ZZ1.EQ.4)then
			DEC_MOD = 2
		elseif(ZZ1.EQ.5)then
			DEC_MOD = 1
		elseif(ZZ1.EQ.6)then
			DEC_MOD = 2
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return			
		endif
	elseif(NN1.EQ.3)then
		if(ZZ1.EQ.2)	then
			DEC_MOD = 3
		elseif(ZZ1.EQ.1)then
			DEC_MOD = 3
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return
		endif
	elseif(NN1.EQ.4)then
		if(ZZ1.EQ.4)	then
			DEC_MOD = 5
		elseif(ZZ1.EQ.1)then
			DEC_MOD = 4
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return			
		endif
	elseif(NN1.EQ.5)then
		if(ZZ1.EQ.2)	then
			DEC_MOD = 3
		elseif(ZZ1.EQ.1)then
			DEC_MOD = 3
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return
		endif
	elseif(NN1.EQ.6)then
		if(ZZ1.EQ.1)	then
			DEC_MOD = 4
		else
			DEC_MOD = 0
			A1 = A
			Z1 = Z
			return
		endif
	else
		DEC_MOD = 0
		A1 = A
		Z1 = Z
		return
	endif


	if(DEC_MOD.EQ.1)	then	!p
		Q = 7.289*(Z - 1) + 8.071*(A - Z) - READ_XBIN(Z - 1, A - 1)
		Q = Q - (7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A))
		if(Q.LE.0.)	Q = 0.000001
		p = 2*Q*0.938*(A - 1)/FLOAT(A)/1000.
		P = SQRT(P)
		A1 = A - 1
		Z1 = Z - 1
		goto	1
	elseif(DEC_MOD.EQ.2)then	!2p
		Q = 7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A)
		EM0 = Z*0.9382719980 + (A - Z)*0.939565330 - Q*0.001
		Q = 7.289*(Z - 2) + 8.071*(A - Z) - READ_XBIN(Z - 2, A - 2)
		EM1 = (Z- 2)*0.9382719980 + (A - Z)*0.939565330 - 
     &			Q*0.001
		EM2 = 0.9382719980
		EM3 = 0.9382719980
		A1 = A - 2
		Z1 = Z - 2
		goto	2
	elseif(DEC_MOD.EQ.3)then	!n
		Q = 7.289*Z + 8.071*(A - Z - 1) - READ_XBIN(Z, A - 1)
		Q = Q - (7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A))

		if(Q.LE.0.)	Q = 0.000001
		p = 2*Q*0.938*(A - 1)/FLOAT(A)/1000.
		P = SQRT(P)
		A1 = A - 1
		Z1 = Z
		goto	1
	elseif(DEC_MOD.EQ.4)then	!2n
		Q = 7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A)
		EM0 = Z*0.9382719980 + (A - Z)*0.939565330 - 
     &			0.001*Q
		Q = 7.289*Z + 8.071*(A - Z - 2) - READ_XBIN(Z, A - 2)
		EM1 = Z*0.9382719980 + (A - Z - 2)*0.939565330 - 
     &			0.001*Q
		EM2 = 0.939565330
		EM3 = 0.939565330
		A1 = A - 2
		Z1 = Z
		goto	2
	elseif(DEC_MOD.EQ.5)then	!α
		Q = 7.289*(Z - 2) + 8.071*(A - Z - 2) -READ_XBIN(Z - 2, A - 4)
		Q = Q + 28.295
		Q = Q - (7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A))

		if(Q.LE.0.)	Q = 0.000001
		p = 2*Q*0.938*(A - 4)*4./FLOAT(A)/1000.
		P = SQRT(P)
		A1 = A - 4
		Z1 = Z - 2
		goto	1
	endif


1	continue	!两体衰变
		call RANDOM_SEED()
		call RANDOM_NUMBER(RAN)
		THETA = ACOS(1. - 2*RAN)
		call RANDOM_SEED()
		call RANDOM_NUMBER(RAN)
		PHI = 2*PI*RAN

		M_3BOD(1, 1) = P*SIN(THETA)*COS(PHI)
		M_3BOD(1, 2) = P*SIN(THETA)*SIN(PHI)
		M_3BOD(1, 3) = P*COS(THETA)
		M_3BOD(2, 1) = -M_3BOD(1, 1)
		M_3BOD(2, 2) = -M_3BOD(1, 2)
		M_3BOD(2, 3) = -M_3BOD(1, 3)
		M_3BOD(3, 1) = 0.
		M_3BOD(3, 2) = 0.
		M_3BOD(3, 3) = 0.
		
	goto	3
2	continue	!三体衰变
	
	!冯老师的方式
C  		PPP=(em0*em0 - (em1 + em2 + em3)**2)*(em0*em0-
C      &      (em1 + em2 - em3)**2)/(4.*em0*em0)
C 
C 		if(PPP.LE.0..OR.(EM0 - EM1 - EM2 - EM3).LE.0.)	PPP = 0.0001
C 
C 		CALL	EVENT(PPP, em3, PKX, PKY, PKZ, PK0)
C 		print*, PKX, PKY, PKZ, PK0
C 		PAUSE
C 		BETX=PKX/(em0-PK0)
C 		BETY=PKY/(em0-PK0)
C 		BETZ=PKZ/(em0-PK0)
C 
C 35		call RANDOM_SEED()
C 		call RANDOM_NUMBER(RAN)
C  		XX    = 1. - 2.*RAN
C 		call RANDOM_SEED()
C 		call RANDOM_NUMBER(RAN)
C 		YY    = 1. - 2.*RAN
C 		call RANDOM_SEED()
C 		call RANDOM_NUMBER(RAN)
C 		ZZ    = 1. - 2.*RAN
C 		RR    = SQRT( XX**2 + YY**2 + ZZ**2 )
C 		if((RR .LT. 0.001) .OR. (RR .GT. 1.))             goto 35
C 
C 		S=(em0 - PK0)**2 - PKX**2 - PKY**2 - PKZ**2
C 
C 		PA2 = (S**2 + em1**4 + em2**4 - 2.*S*em1**2 -
C      &      2.*S*em2**2 - 2.*em1**2*em2**2)/(4.*S)
C 		if(PA2.LE.0.)	PA2 = 0.00000001
C 		PA = SQRT(PA2)
C 		PXX = PA * XX / RR
C 		PYY = PA * YY / RR
C 		PZZ = PA * ZZ / RR
C 		EP1 = SQRT(em1**2 + PA2)
C 		EP2 = SQRT(em2**2 + PA2)
C 
C 		call LORENTZ(BETX,BETY,BETZ,EP1,PXX,PYY,PZZ,
C      &     EPP1,PPX,PPY,PPZ)
C 		call LORENTZ(BETX,BETY,BETZ,EP2,-PXX,-PYY,-PZZ,
C      &     EPP2,PHYPX,PHYPY,PHYPZ)

		M_3BOD(3, 1) = 0.
		M_3BOD(3, 2) = 0.
		M_3BOD(3, 3) = 0.

		M_3BOD(1, 1) = 0.
		M_3BOD(1, 2) = 0.
		M_3BOD(1, 3) = 0.

		M_3BOD(2, 1) = 0.
		M_3BOD(2, 2) = 0.
		M_3BOD(2, 3) = 0.

	!Dalitz图的方法：

3	continue

	return
	end

	FUNCTION	GET_ACC_EKIN(Z, A, M)
	common/SURF_EXT/H1_SURF, INCL4_6, READ_XBIN(0:1000,0:1000)
	DOUBLE PRECISION MASS
	REAL, INTENT(IN) ::	M(3)
	INTEGER, INTENT(IN) ::	Z, A

	MASS = 7.289*Z + 8.071*(A - Z) - READ_XBIN(Z, A)
	if(A.EQ.1)	MASS = 0.	!输入为质子或中子的时候不用结合能
	MASS = Z*0.9382719980 + (A - Z)*0.939565330 - MASS*0.001
	GET_ACC_EKIN = MASS**2 + M(1)**2 + M(2)**2 + M(3)**2
	GET_ACC_EKIN = SQRT(GET_ACC_EKIN)
	GET_ACC_EKIN = (GET_ACC_EKIN - MASS)*1000.
	if(GET_ACC_EKIN.LE.0.)	GET_ACC_EKIN = 0.

	return
	end

      subroutine target_orientation_matrix(mseed,r)
      implicit none
      integer mseed
      real r(3,3),ct,st,phi,psi,cp,sp,cs,ss,ranfl
      real, parameter :: twopi=6.2831853071795864769
      ct=1.-2.*ranfl(mseed)
      st=sqrt(max(0.,1.-ct*ct))
      phi=twopi*ranfl(mseed)
      psi=twopi*ranfl(mseed)
      cp=cos(phi)
      sp=sin(phi)
      cs=cos(psi)
      ss=sin(psi)
      r(1,1)=cp*cs+sp*ct*ss
      r(1,2)=-cp*ss+sp*ct*cs
      r(1,3)=sp*st
      r(2,1)=-sp*cs+cp*ct*ss
      r(2,2)=sp*ss+cp*ct*cs
      r(2,3)=cp*st
      r(3,1)=-st*ss
      r(3,2)=-st*cs
      r(3,3)=ct
      end

      subroutine target_orientation_apply(r,x,y,z)
      implicit none
      real r(3,3),x,y,z,a,b,c
      a=x
      b=y
      c=z
      x=r(1,1)*a+r(1,2)*b+r(1,3)*c
      y=r(2,1)*a+r(2,2)*b+r(2,3)*c
      z=r(3,1)*a+r(3,2)*b+r(3,3)*c
      end
