C*****************************************************************************C
C****************** Initialization code - LQMDini.for for ********************C
C************* Lanzhou 'Quantum' Molecular Dynamics (LQMD) Model *************C
C*****************************************************************************C
C*****************************************************************************C
      program LQMDini
C*****************************************************************************C
      character*50 outfile,aini,bini,iniout,rmsout,
     &pmsout,rdeout,aveout,NEWAINI,NEWBINI,INDAT,DENinn,DENinp,dtiniout
      COMMON/WORK1/ Mnumn,Mnump,Rmaxn,Rmaxp,RMS0N,RMS0P
      common/cmclw/apx,apy,apz,atx,aty,atz
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/rp0/xi(5000),yi(5000),zi(5000),pxi(5000),pyi(5000),
     &pzi(5000),emi(5000),isoi(5000),isyi(5000),idi(5000)
      common/rpi/rmst(5000),pmst(5000),pmsn(5000),pmsz(5000),avet(5000)
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1, 
     &ncolln,ncollp,nblkn,nblkp
      common/pro/nblk1,nbt,nct,pions                                          
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/ee/mseed,nseed
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo,rpt                                 
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/gg/nzta,nzpr,nsym
      common/densi1/denstl(500),dens1(500)                                                
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp
      common/deu2/px1,py1,pz1,px2,py2,pz2                                     
      common/deu3/dg(50,10,2),dtg(2),dtg1(2),
     &dcm(50,10,2),ideu,dtgt(2),dtg1t(2)
      common/kaot/pktot,nkn 
      common/xx/alpha,beta,gamma,symm
      common/aven/tene
      common/encomps/tkine_comp,upaul_comp,upot_comp
      common/tt/prtime(0:8000),nruns,nt,rtime,perc
      common/veloc/velx(500),vely(500),velz(500)                              
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/cmod/coeft,coefp,epp,ecm,ptcm
      common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/data/aini(200),bini(200)
      common/poli/iflagcro,iflagpau,ipoli
      common/totcol/ipuli3,iplin3,iplip3
      common/boite/ann(800),anp(800)
      common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +frtwo(500)
      common/riclw/rrrr(500)
      common/derot/alxyz,idelrot
      common/rdens/rtd(200),ptd(200)
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/bind/progEb,aene,pfc,laesw
      common/CW/c1,c2      
      common/wavel/sl(500) 
      common/deform/beta2
      common/shell/ish,r0,ds,Ri0
      COMMON/INICLW/DISTNN,DISTNP,PHMANN,PHMANP
      COMMON/frictn/gamma_fric_P, tcool, gamma_fric_R
      common/raden/irad,rcount(100),r_tkine(100),r_usky(100),
     +r_ucou(100),r_usym(100),r_ussr(100),r_ugradi(100),
     +r_utau(100),r_couex(100),r_umomd(100),r_upau(100),
     +r_Esh(100),r_tene(100)

      integer :: nstable

      dimension dtot(10),dtotcm(10)                                           
      dimension m1(100)
      character idate*9,itime*8 
      character date1*10,time1*8,zone1*5  
      integer   values(8)     

                                                                                 
************************************************************************
CCC**************** 1 The begin of input             ****************CCC
************************************************************************
      open(2003,file='Input file/LQMDin.in',status='old')
      read (2003,*) npr, nzpr, nta, nzta                ! number of nucleons and protons of projectile and target
      read (2003,*) progEb, beta2                       ! binding energy of projectile and deformation parameter
      read (2003,*) alpha, beta, gamma                  ! Skyrme parameters
      read (2003,*) gsur, gtau, eta, gsym               ! gradient, symmetry, and effective mass parameters
      read (2003,*) symm, rho0, c1, c2                  ! symmetry energy, saturation density, and width parameters
      read (2003,*) asym, bsym, csym, gamsy             ! symmetry energy parameters
      read (2003,*) vpa, vpb, vpc, q01, p01, e01        ! Pauli potential parameters
      read (2003,*) ish, rmax, r0, ds                   ! shell correction parameters
      read (2003,*) initi, nruns                        ! initialization flag and number of runs
      read (2003,*) nt, dt                              ! number of time steps and time step size
      read (2003,*) Ri0, rdt                            ! initial distance and distance step for static fusion barrier calculation
      read (2003,*) cdelt, cmeps, cmom, xsym, momp      !
      read (2003,*) iflag1,iflag2,iflag3,iflag4,iflag5,iflag6
                                                        !iflag1 = -1  ! 1 for initialization; 0 for no; -1 for initialization and output
                                                        !iflag2 = 1   ! 1 for Column (updtyc/updtv); 0 for no (updt)
                                                        !iflag3 = 3   ! 1 for hard equation of state; 
                                                        !             ! 2 for soft equation of state;
                                                        !             ! 3 for soft with momentum dependent potential (updtim+updtsh+constr)
                                                        !iflag4 = 0   ! 1 for static fusion barrier calculation; 0 for no
                                                        !iflag5 = 0   ! 1 for mom. dep. force, 0 for no
                                                        !iflag6 = 0   ! 1 for deuteron production, 0 for no
      read (2003,*) nsym, xmu                           ! symmetry energy parameters
      read (2003,*) gamma_fric_P, gamma_fric_R, tcool   ! friction parameters
      read (2003,*) irad

      if (momp == 0) cdelt = 0.
      if (momp == 1) then
        cmnn = cdelt*(1. + xsym)
        cmnp = cdelt*(1. - xsym)
      else
        cmnn = cmom*(1. + xsym)
        cmnp = cmom*(1. - xsym)
      end if

      v01 = 0.0     ! Pauli potential
      x00min = 0.
      x00max = 0.      
      epp = 0.

      if (initi.ne.0) then
      print*, 'please input nta, nzta'
      read (*,*) nta, nzta
      endif

      vy1 = -0.0024
      arang = 0.5
      ibin = 1
      ird = 0    ! 1 for radial density
      crsc = 0
c     kfile = 10
      kfile = 1       
      ifile = 8
      gradi = 1.
      
      width = (c1*float(npr)**0.3333 + c2)**2.  ! important
      print*, "width=", width       ! take fixed Gaussian wave packet width
      np = nta + npr   
      call width0                   ! Gaussian wave packet width
      
      iprin = 0
      
      if (initi == 0) then
        if (iflag1 < 0.and.epp == 0.and.(npr == 0.or.nta == 0)) then
!          print 901 
!901       format(t2,'input initial outfile name:',$)
          read (2003,2006) iniout ! CHG
!          print 902
!902       format(t2,'input rms outfile name:',$)
          read (2003,2006) rmsout
        end if
        if (iflag1 < 0.and.ibin == 1.and.epp == 0) then
!          print 905
!905       format(t2,'input average binding energy outfile name:',$)
          read (2003,2006) aveout
        end if
        if (iflag1 < 0.and.epp == 0.and.(npr == 0.or.nta == 0)) then
!          print 907
!907       format(t2,'input dtini outfile name:',$)
          read (2003,2006) dtiniout
        end if
      else     
        if (npr > 0.and.nzpr > 0) then 
          print 520, nruns
520       format(t2,'please input aini,number is',1x,i4,':',$)
          read (*,2006) NEWAINI
        end if
        if (epp == 0.and.nta > 0.and.nzta > 0) then 
          print 522, nruns
522       format(t2,'please input bini,number is',1x,i4,':',$)
          read (*,2006) NEWBINI
        end if
        aini(nruns) = NEWAINI
        bini(nruns) = NEWBINI
      end if
2006  format(a50)

      iflagcro = 0
      iflagpau = 1
      ipoli = 2
      ifdelt = 0
      fdelta = 0.5

      idelrot = 0     ! 1  fengzhq
      DISTNN = 1.8    ! 1.5
      DISTNP = 1.0    ! 0.5
      PHMANN = 0.0784
      PHMANP = PHMANN/100.

      open(UNIT=55,FILE='Input file/rmsnp.dat',STATUS='old')
      read (55,*) RMS0N, RMS0P

      Mnumn = 80
      Rmaxn = 8.0
      Mnump = 80
      Rmaxp = 8.0
      deninn = 'Input file/nden.dat'
      deninp = 'Input file/pden.dat'

      write (*,*)
      write (*,*)'The end of all input!!!'
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\ The end of all input!!! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!/////////////////////////////////////////////////////////////////////////////
      pi=3.1415927
      al =4.0*width      
      sal =(al*pi)**1.5*rho0
      alnew=4.*(width+fdelta**2./4.)
      salnew=(alnew*pi)**1.5*.165
      t0 =       -0.124                                                       
      t1 =        0.0705                                                      
      vy=vy1  
                        
      rrp=r0*npr**0.3333
      rrt=r0*nta**0.3333
      RMS0N=RMS0N*SQRT(1.+5.*BETA2**2./4./PI)
      RMS0P=RMS0P*SQRT(1.+5.*BETA2**2./4./PI)
                
        alp=2./al                                                               
        bet=.346                                                                
        fdr=alp*bet/(alp+2*bet)                                                 
        fd1=8*(alp*bet)**1.5/(alp+2*bet)**1.5/(2*alp+bet)**1.5                  
        fdp=1./4./(2*alp+bet)/.197/.197                                         
        xfdr=sqrt(3./2./fdr)                                                    
        xfdp=sqrt(3./2./fdp)                                                    
1237    format(///,1x,'fd1 = ',f10.4,' d rms = ',f10.4,'  fm            
     &  d rms p  =',f10.4,'  gev ')  
                                                
        call initp 
        do 664 j=1,99999                                       
        dr2=float(j)/40.                                                        
        dr=sqrt(dr2)                                                            
        ir=int(dr*100.+0.5)                                                     
        ih(j+1)=ir                                                              
664     continue                                                                
        ih(1)=1                                                                 
1006    format(1x,'v0=',f6.3,'q0=',f6.3,'p0=',f6.3,'crsc=',f6.3,
     &        'vs=',f6.3)
111     format(1x,'sal= ',f6.2,' al= ',f6.2,' vy= ',f6.3,                        
     &    'arang= ',f6.2)                                                       
      iwrite=abs(kfile/dt) + 0.5
      kwrite=abs(kfile/dt)+0.5

       if (iflag3  ==  1.and.iflag5  ==  0) then       !HM for QMD or IQMD                    
             alpha = -0.124                                                     
             beta  =  0.0705                                                    
             gamma =  2.     
       else if (iflag3  ==  1.and.iflag5  ==  1) then  !Mom. Dep. HM for QMD or IQMD          
             alpha = -0.1292                                                    
             beta  =  0.0594                                                    
             gamma =  2.09                                                      
       if (v01 == 0.) goto 11111                                                
             alpha= -0.1276                                                      
             beta=   0.080                                                   
             gamma=  2.                                                          
             gradi=  0.                                                          
11111  continue
                                                            
       else if (iflag3  ==  2.and.iflag5  ==  0) then  !SM for QMD or IQMD                    
             alpha = -0.356                                                     
             beta  =  0.303                                                     
             gamma =  7./6.                                                     
       else if (iflag3  ==  2.and.iflag5  ==  1) then  !Mom. Dep. SM for QMD or IQMD         
             alpha = -0.3901                                                    
             beta  =  0.3203                                                    
             gamma =  1.14                                                      
       if (v01 == 0.) goto 11112                                                 
             alpha= -0.332                                                       
             beta=   0.274                                                         
             gamma=  7./6.                                                        
             gradi=  0.                                                          
11112  continue                                                                                                                            
       end if
c                                                                         
        pota =2.*alpha                                                          
        potb=beta/(-alpha)*gamma/(gamma+1.)                                     
        potc = gamma - 1.                                                       
        potd=pota/sal/al                                                        
        pote=potb/(sal)**potc                                                   
        potf=cdelt/sal*2./al                                                  
        potg=gradi/sal*0.16/(al/4)**2                                           
        pots=2.*symm/sal/al
        potdnew=pota/salnew/alnew
        potenew=potb/salnew**potc
        Cex=-3./4.*0.00144*(3/pi)**0.3333
c-----------***********************-------------
        if (initi == 1) then
      open(400,file='Uint-R.dat',status='unknown')
      z0=Ri0/2.  
        call INIT2
        if (iflag2 == 1) call coulgn 

1122      continue
      call avene(avpot,aveki,Eshp,1)
      EEp=avpot+aveki+Eshp
      call avene(avpot,aveki,Esht,2)
      EEt=avpot+aveki+Esht
      call avene(avpot,aveki,Eshpt,3)      
      Ept=avpot+aveki+Eshpt
      write (400,*) rpt, Ept-EEp-EEt 
      z0=z0-rdt
      if (z0 > 0.5) goto 1122
      stop 'The calculation of the static fusion barrier is ok !'
        end if  
c-----------***********************-------------

      if (epp > 0.and.npr > 0.and.nzpr > 0
     &  .and.nta > 0.and.nzta > 0) then
      iopen=0
        go to 19103
19922   print 29933,ifile
29933   format(/t2,'unit=',i3,1x,'open error 1 !'/)
        iopen=iopen+1
        if (iopen > 100)stop
19103   open(ifile,file=outfile,err=19922,form='formatted',status='new')

        write (ifile,1120) nt,x00min,x00max,epp,dt,al
        write (ifile,1130) nta,npr,nzta,nzpr                                     

        close(ifile)
      end if

1120  format(/1x,i4,5(f8.2,1x))                                    
1130  format(1x,4(i3,1x))                                                     
112   format(1x,'steps=',i4,' bmin=',f6.2,' bmax=',f6.2,
     &       'energy=',f6.2,'  dt=  ',f6.3)
113   format(1x,'n target = ',i3,'  n proj = ',i3,                            
     &       '  z target = ',i3,'  z proj = ',i3)                                    
114   format(1x,'number of repetitions = ',i3)                                
1251  format(1x,' velocity dependent forces')                                 
115   format(1x,'soft e.o.s')                                                 
1151  format(1x,'stiff e.o.s')                                                
116   format(1x,'no coulomb')                                                 
117   format(1x,'coulomb')                                                    
118   format(1x,'mean field only')                                            
119   format(1x,'reduced collisions')                                         
120   format(1x,'all collisions')                                             
121   format(1x,'relativistic dynamics')                                      
122   format(1x,'non relaticvistic dynamics')                                 
c                                                                               
      if (iflag2 == 1) call coulgn       
!!    call coulg1                                      
c                                                                               
      if (npr == 0.or.nta == 0.or.epp == 0.) then
        gam3=1.
        goto 707  
      end if
      e1=0.938+0.001*epp                                                      
      pz1=sqrt(e1*e1-0.938**2)                                                
      bcm=pz1/(e1+float(nta)/float(npr)*0.938)                                
      bzz=bcm                                                                 
      gamcm=1./sqrt(1-bcm**2)                                                 
      gam3=gamcm                                                              
      pb=bcm*pz1                                                              
      gg=(gamcm-1)*pb/bcm/bcm-gamcm*e1                                        
      ppo=gg*bcm+pz1                                                          
      tpo=float(npr)/float(nta)*ppo                                           
      if (iflag6  ==  0) goto 707                                              
707   continue  
      do 4723 ijk = 1, nt                                                      
4723    prtime(ijk) = 0                                                      
        if (iflag1 < 0.and.epp == 0.and.(npr == 0.or.nta == 0)) then
          open(50,file=iniout,form='formatted',status='unknown')
          open(60,file=rmsout,form='formatted',status='unknown')
          open(80,file=dtiniout,form='formatted',status='unknown')
        end if
        if (iflag1 < 0.and.ibin == 1.and.epp == 0) then
          open(70,file=aveout,form='formatted',status='unknown')
          open(90,file='encomps.dat',form='formatted',status='unknown')
        end if
        call clallnew

        nstable = 0
        do 556 nru=1,nruns
        write (*,*) 'Running NO: ',nru                                                   
cc      if (nru == 1)          write (*,*) 'Running NO: ',nru
cc      if (mod(nru,10) == 0)  write (*,*) 'Running NO: ',nru
      mseed=28+nru*15
      nseed=48+nru*25
c      call date_and_time(date1,time1,zone1,values)
c        mseed=values(5)*31+values(6)*95+values(7)*89439
c        nseed=values(5)*21+values(6)*41+values(7)*13974
      if (iflag1 < 0.and.epp == 0.and.(npr == 0.or.nta == 0)) then
        do 101 i = 1, nt/iwrite+1
          rmst(i)=0.      
          pmst(i)=0.
          pmsN(i)=0.
          pmsZ(i)=0.
          avet(i)=0.
          xi(i)=0.
          yi(i)=0.
          zi(i)=0.
          pxi(i)=0.
          pyi(i)=0.
          pzi(i)=0.
          emi(i)=0.
          isoi(i)=0
          isyi(i)=0
          idi(i)=0
101     continue
      end if
      if (nta > 0.and.nzta > 0.and.npr > 0.and.nzpr > 0.and.epp
     &       >= 0.) then
c        open(ifile,file=outfile,form='formatted',
c     &             status='old',access='append')

        if (x00min >= x00max) then
          x00 = x00max
        else
          x00 = x00min+(x00max-x00min)*sqrt(ran(mseed))
        end if
        write (ifile,19201)x00
19201   format(1x,f9.4)
        call coulomb_modify
        close(ifile)
      end if
      
      ideut = 1
      if (iflag6 == 1) ideut = 2
      do 5561 ideu = 1, ideut
        do klk = 1, 2
          dtg(klk) = 0
          dtg1(klk) = 0
        end do
        if (ideu == 2) goto 5562
        do kkl = 1, np
          iev(kkl) = 1
        end do
        goto 5563
5562    call minspa
5563    continue   

      ntstep=0
      rtime=0.                                                 

c  -------------************---------------
        if (initi == 0) call init0
c  -------------************---------------
c     fengzhq
cc      if (NTA <= 0) THEN
cc      open(14,FILE=DENinP,FORM='FORMATTED',STATUS='OLD')
cc      open(15,FILE=DENinN,FORM='FORMATTED',STATUS='OLD')
cc      END IF
cc      CLOSE(14)
cc      CLOSE(15)
cc        END IF                                                
      if (iflag1 < 0.and.epp == 0..and.ird == 1
     +   .and.(nta == 0.or.npr == 0)) then
      call rpdens
      end if

      if (nt > 0) then
      perc=(float((nru-1)*nt)+0.5)/float(nruns*nt)*100.
      else if (nt == 0) then
      perc=100.
      end if
      ncoll=0
      ncolln=0
      ncollp=0
      nblock=0
      nblkn=0
      nblkp=0
      ipuli3=0
      iplin3=0
      iplip3=0

c     fengzhq
cc      call prot
                                 
        ncoll1=0                                                                
        nblk1=0  
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                                                                 
        do 26  k=1,nt                                                        
        ntstep=k                                                           
      perc=float((nru-1)*nt+k)/float(nruns*nt)*100.
      if (dt > 0.)rtime=ntstep*dt
      if (dt < 0.)rtime=-ntstep*dt
      do 779 i= 1, np
        write (80,'(f10.4,1x,3f12.6,i3)') rtime, x(i), y(i), z(i), 
     &    iso(i)
779   continue

       if (iflag2  ==  0) call updt                                             
       if (iflag2  ==  1.and.iflag5 == 0.and.iflag3.ne.3) call updtyc
       if (iflag2  ==  1.and.iflag5 == 1.and.iflag3.ne.3) call updtv
       if (iflag3  ==  3) call  updtim
       if (iflag3  ==  3) call  updtsh
       if (iflag3  ==  3) call  constr
                                    
      if (iflag1 < 0.and.epp == 0..and.ird == 1.and.
     &      mod(k,kwrite) == 0.and.(nta == 0.or.npr == 0)) then
      call rpdens
      end if

      if ((iflag1 < 0.and.initi == 1).or.iflag1 > 0) then
      if (mod(k,kwrite) == 0) then
      junit=71+ntstep/kwrite
      write (junit,*)'\\ time=',rtime
      do 9 i=1,np
      write (junit,911)z(i),x(i)
9      continue
911      format(1x,2f9.4)
      end if
      end if

          if (mod(k,iwrite) == 0) then      
C
      ipuli3=ipuli3+ncoll-nblock
      iplin3=iplin3+ncolln-nblkn
      iplip3=iplip3+ncollp-nblkp
C
      call prot
            write (60,'(f10.4,1x,f9.4)') rtime, rmst(llt)

      if (epp == 0..and.ibin == 1) then
      call avene(avpot,aveki,Esh,1)
            write (70,'(f10.4,1x,f9.4)') rtime, aene
            write (90,'(f10.4,4(1x,f9.4))') rtime, aene,
     &  tkine_comp, upaul_comp, upot_comp
      end if
        end if

      if (epp == 0..and.lswitch == 0) goto 54321
9955    continue                                                                
          if (k == nt) then
        if (iflag1 < 0.and.(nta == 0.or.npr == 0).and.epp == 0.) then
      open(4,file='0result.inf',status='unknown',access='append')
      write (4,*)'Running NO: ',nru,' is stable !!'
      write (4,*)'            The eorrors of C. M. (fm): ' 
      if (npr > 0) then
        write (4,*)'  Xcm= ',apx,'  Ycm= ',apy,'  Zcm= ',apz
      else
        write (4,*)'  Xcm= ',atx,'  Ycm= ',aty,'  Zcm= ',atz
      end if
      write (4,*)'     The R(i)s of protons:       '
      if (npr > 0) then
        do i=1,nzpr
          write (4,*) rrrr(i)
        end do
      else
        do i=1,nzta
          write (4,*) rrrr(i)
        end do
      end if
      write (4,*)'     The R(i)s of neutrons:       '
      if (npr > 0) then
        do i=nzpr+1,npr
          write (4,*) rrrr(i)
        end do
      else
        do i=nzta+1,nta
          write (4,*) rrrr(i)
        end do
      end if
      close (4)
      do i=1,npr+nta
        if (I == 1) write (50,*) '    Running NO: ',nru
        write (50,77777)xi(i),yi(i),zi(i),pxi(i),pyi(i),pzi(i),emi(i),
     &               isoi(i),isyi(i),idi(i)
      end do
      nstable = nstable + 1
        end if
        end if
77777      format(1x,7f9.4,2x,2i3,2x,i3)      

26      continue      
                                                          
54321   pionst=pionst+pions                                                     
334     format(1x,'run # = ',i3,' # of pions all runs = ',f6.1)                 
        dtgt(ideu)=dtgt(ideu)+dtg(ideu)                                         
        dtg1t(ideu)=dtg1t(ideu)+dtg1(ideu)                                      
3334     format(1x,'red   # of deut all runs = ',f6.2,' # of prod               
     1  deuteroons = ',f10.4)                                                   
3335     format(1x,'full   # of deut all runs = ',f6.2,' # of prod              
     1  deuterons = ',f10.4)                                                    
5561    continue                                                                
        if (nstable  >=  3.and.epp  ==  0.) goto 557

556     continue                                                              
557   continue
c  -----------------------------------------------------------------------------
      if (iflag1 < 0.and.npr == 0.and.epp == 0.) then
      close(50)
      close(60)
            close(80)
      close(90)
      else if (iflag1 < 0.and.nta == 0.and.epp == 0) then
           close(40)
      close(61)
            close(80)
      close(91)
      end if
      if (iflag1 < 0.and.ibin == 1.and.epp == 0) close(70)
        do 8009 ideu=1,2                                                        
1289    format(////1x,' full number of deuterons ',////)                        
1299    format(////1x,' reduced number of deuterons ',////)                     
        do 8003 ll=1,9                                                          
        dtot(ll)=0.                                                             
        dtotcm(ll)=0                                                            
8003    continue                                                                
1346    format(1x,' deuteron lab spectra ')                                     
        do 8004 l=1,50                                                          
        do 8005 ll=1,9                                                          
        dtot(ll)=dtot(ll)+dg(l,ll,ideu)                                         
8005    continue                                                                
1234    format(1x,9f7.3)                                                        
8004    continue                                                                
1347    format(1x,' deuteron lab angular distribution ')                        
1345    format(1x,' deuteron cm spectra ')                                      
        do  8006 l=1,20                                                         
        do 8007 ll=1,9                                                          
        dtotcm(ll)=dtotcm(ll)+dcm(l,ll,ideu)                                    
8007    continue                                                                
8006    continue                                                                
        cmd=0.                                                                  
        do 8008 ll=1,9                                                          
        cmd=cmd+dtotcm(ll)                                                      
8008    continue                                                                
1547    format(1x,' deuteron cm angular distribution ')                         
1348    format(1x,' tot# = ',f10.4,'  prod #= ',f10.4)                          
8009    continue                                                                
        nx=ncoll1-nblk1                                                         
1160      format(1x,3i5,f6.1,f9.4)
c       call date(idate)
c       call time(itime) 
      print 124,idate,itime
124      format(/t2,'the program end at ',a9,1x,a8/)

        if (irad == 1) call avene(avpot,aveki,Esh,3)

        stop 'the program end normally, ok!'
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

              subroutine ocu(xx,yy,zz,pxx,pyy,pzz,ocuph)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
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
        if (aa <= 0.) then                                                       
        call erf(-aa,epp1(k))                                                   
        else                                                                    
        call erf(aa,epp1(k))                                                    
        epp1(k)=-epp1(k)                                                        
        end if
        bb=pp2(k)                                                               
        if (bb <= 0.) then                                                       
        call erf(-bb,epp2(k))                                                   
        epp2(k)=-epp2(k)                                                        
        else                                                                    
        call erf(bb,epp2(k))                                                    
        end if                                                                  
        epp=epp*(epp1(k)+epp2(k))                                               
20      continue                                                                
        sum=sum+epp                                                             
10      continue                                                                
        ocuph=sum/2.**6                                                         
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

              subroutine nrcoll(crsc)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1, 
     +      ncolln,ncollp,nblkn,nblkp
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/ee/mseed,nseed
        common/deu2/px1,py1,pz1,px2,py2,pz2                                     
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
      ncoll=0                                                                 
        nblock=0                                                                
c find nn collisions involving particle i1                                      
        do 21 i1=2,np                                                           
        pz1=pz(i1)                                                              
        px1=px(i1)                                                              
        py1=py(i1)                                                              
        em1=em(i1)                                                              
        ip=i1-1                                                                 
c  find nn collisions between i1 and i2                                         
        do 22 i2=1,ip                                                           
        if (il(i1) == i2.and.il(i2) == i1)goto 22                                
        if (iflag1 == 1.and.id(i2) == 1.and.id(i1) == 1)goto22                   
        if (iflag1 == 1.and.id(i2) == 51.and.id(i1) == 51)goto22                 
        pz2=pz(i2)                                                              
        px2=px(i2)                                                              
        py2=py(i2)                                                              
        em2=em(i2)                                                              
        cror1=1.142                                                             
        cror=cror1+crsc
c    cross section =crsc input 
        xr=x(i1)-x(i2)                                                          
        yr=y(i1)-y(i2)                                                          
        zr=z(i1)-z(i2)                                                          
        r2=xr*xr+yr*yr+zr*zr                                                    
        if (r2 > cror*cror) go to 22                                           
        prx=(px(i1)-px(i2))/2.                                                  
        pry=(py(i1)-py(i2))/2.                                                  
        prz=(pz(i1)-pz(i2))/2.                                                  
        pr=sqrt(prx*prx+pry*pry+prz*prz)                                        
        vr=2*pr/0.938                                                           
c  cross section 41mb  1.142 fm                                                 
        dzz=(prx*xr+pry*yr+prz*zr)/pr                                           
        bb2=r2-dzz*dzz                                                          
        if (bb2 > cror1*cror1) go to 22                                         
        dd=vr*dt/2.                                                             
        if (abs(dd) < abs(dzz)) go to 22                                        
c a collision has occurred                                                      
        pcz=(pz(i1)+pz(i2))/2.                                                  
        pcx=(px(i1)+px(i2))/2.                                                  
        pcy=(py(i1)+py(i2))/2.                                                  
        ncoll=ncoll+1 
        cth=1.-2.*ran(mseed) 
        phi=6.28*ran(mseed) 
        sth=sqrt(1.-cth**2.)
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
        if (ntag == -1)go to 32            
        call paulin(i2)                   
        if (ntag == 0)go to 551            
32      nblock=nblock+1                   
        pz(i1)=pcz+prz                    
        px(i1)=pcx+prx                    
        py(i1)=pcy+pry                    
        pz(i2)=pcz-prz                    
        px(i2)=pcx-prx                    
        py(i2)=pcy-pry                    
        goto22                            
551     continue                          
        il(i1)=i2                                      
        il(i2)=i1                         
        if (iflag6 == 1)call dtest(i1,i2,em1,em2)       
        id(i1)=id(i1)+1                                
        id(i2)=id(i2)+1                                
22      continue                                       
21      continue                                       
        return                                         
      end                                            
c///////////////////////////////////////////////////////////////////////////////

            subroutine paulin(i)                                            

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/bb/np,ntag,x00                                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/ee/mseed,nseed
        common/paulis/radn,padn,anorm,const,radn2,padn2                       
      common/poli/iflagcro,iflagpau,ipoli
      common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +      frtwo(500)
        dimension con(500)                                                    
        if (em(i) > 1.)return                                                 
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
        con(j)=0.                                                             
        if (j == i)go to 1                                                     
        if (em(j) > 1.)goto1                                                  
        dx=x(j)-x1                                                            
        dy=y(j)-y1                                                            
        dz=z(j)-z1                                                            
        dpx=px(j)-px1                                                         
        dpy=py(j)-py1                                                         
        dpz=pz(j)-pz1                                                         
        pr2=dpx*dpx+dpy*dpy+dpz*dpz                                           
        r2=dx*dx+dy*dy+dz*dz                                                  
        if (pr2 > padn2) goto 1                                               
        if (r2 > radn2) goto 1                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)                                                            
        densty = (fr(i)+1.)/sal                                               
        pot=- .353*densty+.303*densty**(7./6.)                                
      if (ifdelt >= 1) then
      dentwo=(frtwo(i)+1.)/salnew
        pot=- .353*dentwo+.303*densty**(7./6.)
      end if                                
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
        if (contr < 0.00001)goto 1                                             
        nn=nn+1                                                               
        contrs=contrs+contr                                                   
1       continue                                                              
        if (nn == 0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np                                                           
        if (k == i)goto2                                                       
        dx=x(k)-x1                                                            
        dy=y(k)-y1                                                            
        dz=z(k)-z1                                                            
        pdx=px(k)-px1                                                         
        pdy=py(k)-py1                                                         
        pdz=pz(k)-pz1                                                         
        dr=sqrt(dx**2+dy**2+dz**2)                                            
        pdr=sqrt(pdx**2+pdy**2+pdz**2)                                        
c
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm                                    
        if (xcos > 0.707)sum=sum+con(k)                                       
        if (xcos > 0.707)nn1=nn1+1                                            
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm                           
        if (pxcos > 0.707)psum=psum+con(k)                                    
        if (pxcos > 0.707)nn2=nn2+1                                           
2       continue                                                              
        sum=sum*6.83                                                          
        psum=psum*6.83                                                        
        sumf=amax1(sum,psum)                                                  
        if (sumf < contrs)sumf=contrs                                         
        if (etot > 0.005)sumf=contrs                                          
        if (sumf > ran(mseed))ntag=-1        
      if (ipoli == 0) ntag=-1
      if (ipoli == 1) ntag=0
      if (ipoli < 0.or.ipoli > 1) goto 3
3      continue
        return                                                                
        end    
c///////////////////////////////////////////////////////////////////////////////

            subroutine pauli12(i1,i2)                                            

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/bb/np,ntag,x00
      common/ee/mseed,nseed      
                                                
      ntag=0
      f1=Fn(i1)
      f2=Fn(i2)
      f0=0.622
      fi1=amin1(1.,(f1-f0)/(1.-f0))
      fi2=amin1(1.,(f2-f0)/(1.-f0))
      
      fi12=1.-(1.-fi1)*(1.-fi2)
      if (fi12 >= ran(nseed)) ntag=-1     
      
      return
      end      
c///////////////////////////////////////////////////////////////////////////////

            subroutine paur34(i)                                            

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/bb/np,ntag,x00                                                 
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/ee/mseed,nseed
      common/paulis/radn,padn,anorm,const,radn2,padn2                       
      common/poli/iflagcro,iflagpau,ipoli
      common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +      frtwo(500)
      dimension con(500) 
                                                           
        if (em(i) > 1.) return                                                 
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
        con(j)=0.                                                             
        if (j == i)           goto 1                                                     
        if (em(j) > 1.)      goto 1                                                  
      if (iso(j).ne.iso(i)) goto 1
        dx=x(j)-x1                                                            
        dy=y(j)-y1                                                            
        dz=z(j)-z1                                                            
        dpx=px(j)-px1                                                         
        dpy=py(j)-py1                                                         
        dpz=pz(j)-pz1                                                         
        pr2=dpx*dpx+dpy*dpy+dpz*dpz                                           
        r2=dx*dx+dy*dy+dz*dz                                                  
        if (pr2 > padn2) goto 1                                               
        if (r2 > radn2)  goto 1                                                
        pr=sqrt(pr2)                                                          
        r=sqrt(r2)                                                            
        densty= (fr(i)+1.)/sal                                               
        pot=-.353*densty+.303*densty**(7./6.)                                
      if (ifdelt >= 1) then
      dentwo=(frtwo(i)+1.)/salnew
        pot=-.353*dentwo+.303*densty**(7./6.)
      end if                                
        etot=pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
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
        if (contr < 0.00001) goto 1                                             
        nn=nn+1                                                               
        contrs=contrs+contr                                                   
1       continue                                                              
        if (nn == 0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np                                                           
        if (k == i) goto 2                                                       
      if (ISO(K).NE.ISO(I)) GOTO 2
        dx=x(k)-x1                                                            
        dy=y(k)-y1                                                            
        dz=z(k)-z1                                                            
        pdx=px(k)-px1                                                         
        pdy=py(k)-py1                                                         
        pdz=pz(k)-pz1                                                         
        dr=sqrt(dx**2+dy**2+dz**2)                                            
        pdr=sqrt(pdx**2+pdy**2+pdz**2)                                        
c      print*,'dr=',dr,'rcm=',rcm
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm                                    
        if (xcos > 0.707)sum=sum+con(k)                                       
        if (xcos > 0.707)nn1=nn1+1                                            
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm                           
        if (pxcos > 0.707)psum=psum+con(k)                                    
        if (pxcos > 0.707)nn2=nn2+1                                           
2       continue                                                              
        sum=sum*6.83                                                          
        psum=psum*6.83                                                        
        sumf=amax1(sum,psum)                                                  
        if (sumf < contrs)sumf=contrs                                         
        if (etot > 0.005) sumf=contrs 
        if (sumf > ran(mseed))ntag=-1
      if (ipoli == 0) ntag=-1
      if (ipoli == 1) ntag=0
      if (ipoli < 0.or.ipoli > 1) goto 3
3      continue
        
      return                                                                
        end                                                                   
c///////////////////////////////////////////////////////////////////////////////

                  subroutine initp

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common/paulis/radn,padn,anorm,const,radn2,padn2                       
      common/poli/iflagcro,iflagpau,ipoli
        pi=3.14159                                                            
        aux=(4.*pi/3.)**2                                                     
        hc=1.24            
C
         if (iflagpau == 1) then
        const=4./9.*pi*pi*4./(hc)**3./2.
       radn=3.367
        padn=hc/radn/aux**(1./3.)/2.**(1./3.)             
      else 
      const=4./9.*pi*pi*4./(hc)**3.
      radn=3.
      padn=hc/radn/aux**.333/4.**.333
      end if
C
        radn2=(2.*radn)**2                                                    
        padn2=(2.*padn)**2                                                    
233     format(1x,'pauli radii   r = ',f6.2,'   p = ',f6.2)                   
        return                                                                
        end                                                                   

clw\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

            subroutine pauover(i)

clw/////////////////////////////////////////////////////////////////////////////
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     +  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,ntag,x00
      common/skyrme/al,sal,t0,t1,arang,vy,width
      common/ee/mseed,nseed
      common/poli/iflagcro,iflagpau,ipoli
      common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +      frtwo(500)
      dimension con(500)
c
      if (em(i) > 1.) return
      ntag=0
      contrs=0.
      x1=x(i)
      y1=y(i)
      z1=z(i)
      px1=px(i)
      py1=py(i)
      pz1=pz(i)
c
      do 1 j=1,np
      if (j == i) goto 1
      if (em(j) > 1.) goto 1
      if (iso(j).ne.iso(i)) goto 1
      dx=x(j)-x1
      dy=y(j)-y1
      dz=z(j)-z1
      dpx=px(j)-px1
      dpy=py(j)-py1
      dpz=pz(j)-pz1
      r2=dx*dx+dy*dy+dz*dz
      p2=dpx*dpx+dpy*dpy+dpz*dpz
        densty = (fr(i)+1.)/sal                                               
        pot=- .353*densty+.303*densty**(7./6.)                                
      if (ifdelt >= 1) then
      dentwo=(frtwo(i)+1.)/salnew
        pot=- .353*dentwo+.303*densty**(7./6.)
      end if                                
        etot = pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
C----------------------------------------------------------------------
      contr=xexp(-r2/width/4.-p2*width/0.197/0.197)
C----------------------------------------------------------------------
      con(j)=contr
        xcm=xcm+x(j)                                                          
        ycm=ycm+y(j)                                                          
        zcm=zcm+z(j)                                                          
        pxcm=pxcm+px(j)                                                       
        pycm=pycm+py(j)                                                       
        pzcm=pzcm+pz(j)                                                       
        if (contr < 0.00001)goto1                                             
        nn=nn+1                                                               
      contrs=contrs+contr
1      continue
        if (nn == 0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np                                                           
        if (k == i)goto2                                                       
      if (ISO(K).NE.ISO(I)) GOTO2
        dx=x(k)-x1                                                            
        dy=y(k)-y1                                                            
        dz=z(k)-z1                                                            
        pdx=px(k)-px1                                                         
        pdy=py(k)-py1                                                         
        pdz=pz(k)-pz1                                                         
        dr=sqrt(dx**2+dy**2+dz**2)                                            
        pdr=sqrt(pdx**2+pdy**2+pdz**2)                                        
c      print*,'dr=',dr,'rcm=',rcm
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm                                    
        if (xcos > 0.707)sum=sum+con(k)                                       
        if (xcos > 0.707)nn1=nn1+1                                            
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm                           
        if (pxcos > 0.707)psum=psum+con(k)                                    
        if (pxcos > 0.707)nn2=nn2+1                                           
2       continue                                                              
        sum=sum*6.83                                                          
        psum=psum*6.83                                                        
        sumf=amax1(sum,psum)                                                  
        if (sumf < contrs)sumf=contrs                                         
        if (etot > 0.005)sumf=contrs 
        if (sumf > ran(mseed))ntag=-1                            
      if (ipoli == 0) ntag=-1
      if (ipoli == 1) ntag=0
      if (ipoli < 0.or.ipoli > 1) goto 3
3      continue
        return                                                                
        end                                                                   


c///////////////////////////////////////////////////////////////////////////////

            subroutine paurho(i)                                            

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/bb/np,ntag,x00 
        common/skyrme/al,sal,t0,t1,arang,vy,width 
        common/ee/mseed,nseed
      common/poli/iflagcro,iflagpau,ipoli
      common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +      frtwo(500)
        dimension con(500),sden(500),rad(500),pad(500) 
        if (em(i) > 1.)return                                                 
        PI=3.1415927
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
      end do

      ccwl=(2.*PI*width)**1.5
      do j=1,np
         do 96916 k=1,np
         if (iso(k).ne.iso(j)) goto 96916
         dis2=(x(k)-x(j))**2.+(y(k)-y(j))**2.+(z(k)-z(j))**2.
         sden(j)=sden(j)+xexp(-dis2/2./width)/ccwl
96916         continue
      end do

      consta=1.24*(9./PI/PI/32.)**(1./3.)
      do j=1,np
      pad(j)=.197*(3.*PI*PI*sden(j))**(1./3.)
      rad(j)=consta/pad(j)
      end do

        do 1 j=1,np                                                           
        con(j)=0.                                                             
        if (j == i)go to 1                                                     
        if (em(j) > 1.) goto 1                                                  
      if (iso(j).ne.iso(i)) goto 1
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
      densty=(fr(i)+1.)/sal
        pot=- .353*densty+.303*densty**(7./6.)
      if (ifdelt >= 1) then
      dentwo=(frtwo(i)+1.)/salnew
        pot=- .353*dentwo+.303*densty**(7./6.)
      end if                                
        etot = pot + sqrt(px1*px1+py1*py1+pz1*pz1+em(i)*em(i))-em(i)          
      
      if (r >= (rad(j)+rad(i)).or.pr >= (pad(j)+pad(j))) goto 1
      if (r <= abs(rad(j)-rad(i))) then
      rrmin=amin1(rad(j),rad(i))
      phaser=4./3.*PI*(rrmin)**3.
      end if
      if (r > abs(rad(j)-rad(i)).and.r < (rad(j)+rad(j))) then
      rh1=(rad(j)**2.-(r-rad(i))**2.)/2./r
      rh2=(rad(i)**2.-(r-rad(j))**2.)/2./r
      phaser=PI/3.*(rh1**2.*(3.*rad(i)-rh1)+rh2**2.*(3.*rad(j)-rh2))
      end if

      if (pr <= abs(pad(j)-pad(i))) then
      ppmin=amin1(pad(j),pad(i))
      phasep=4./3.*PI*(ppmin)**3.
      end if
      if (pr > abs(pad(j)-pad(i)).and.pr < (pad(j)+pad(i))) then
      ph1=(pad(j)**2.-(pr-pad(i))**2.)/2./pr
      ph2=(pad(i)**2.-(pr-pad(j))**2.)/2./pr
      phasep=PI/3.*(ph1**2.*(3.*pad(i)-ph1)+ph2**2.*(3.*pad(j)-ph2))
      end if

      phasei=1.24**3./2.
      contr=phaser*phasep/phasei
        con(j)=contr                                                          
        xcm=xcm+x(j)                                                          
        ycm=ycm+y(j)                                                          
        zcm=zcm+z(j)                                                          
        pxcm=pxcm+px(j)                                                       
        pycm=pycm+py(j)                                                       
        pzcm=pzcm+pz(j)                                                       
        if (contr < 0.00001) goto 1                                             
        nn=nn+1                                                               
        contrs=contrs+contr                                                   
1       continue                                                              
        if (nn == 0)return                                                     
        xcm=xcm/nn-x1                                                         
        ycm=ycm/nn-y1                                                         
        zcm=zcm/nn-z1                                                         
        rcm=sqrt(xcm*xcm+ycm*ycm+zcm*zcm)                                     
        pxcm=pxcm/nn-px1                                                      
        pycm=pycm/nn-py1                                                      
        pzcm=pzcm/nn-pz1                                                      
        prcm=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm)                              
        do 2 k=1,np                                                           
        if (k == i) goto 2                                                       
      if (iso(k).ne.iso(i)) goto 2
        dx=x(k)-x1 
      dy=y(k)-y1
        dz=z(k)-z1 
        pdx=px(k)-px1 
        pdy=py(k)-py1
        pdz=pz(k)-pz1 
        dr=sqrt(dx**2+dy**2+dz**2)
        pdr=sqrt(pdx**2+pdy**2+pdz**2) 
c      print*,'dr=',dr,'rcm=',rcm
        xcos=(dx*xcm+dy*ycm+dz*zcm)/dr/rcm  
        if (xcos > 0.707)sum=sum+con(k) 
      if (xcos > 0.707)nn1=nn1+1
        pxcos=(pdx*pxcm+pdy*pycm+pdz*pzcm)/pdr/prcm  
        if (pxcos > 0.707)psum=psum+con(k)  
        if (pxcos > 0.707)nn2=nn2+1 
2       continue 
        sum=sum*6.83  
        psum=psum*6.83  
        sumf=amax1(sum,psum)
        if (sumf < contrs)sumf=contrs                                         
        if (etot > 0.005)sumf=contrs  
      if (sumf > ran(mseed))ntag=-1
      if (ipoli == 0) ntag=-1
      if (ipoli == 1) ntag=0
      if (ipoli < 0.or.ipoli > 1) goto 3
3      continue
        return                                                                
        end                                                                                                                                  
c///////////////////////////////////////////////////////////////////////////////

              subroutine init2 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      character*50 aini,bini,targ,proj

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1, 
     +      ncolln,ncollp,nblkn,nblkp
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/gg/nzta,nzpr,nsym                                                     
        common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
        common/ee/mseed,nseed
        common/tt/prtime(0:8000),nruns,nt,rtime,perc
        common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
        common/data/aini(200),bini(200)
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp

        dimension dvdpx(500),dvdpy(500),dvdpz(500)                              
        hc=1.24                                                        
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
      targ=bini(nruns)
      proj=aini(nruns)
cc      print*,targ,proj
      if (npr.ne.0.and.nzpr.ne.0) then

      open(30,file=proj,form='formatted',status='old')
        do 10 i=1,npr                                                         
        read (30,1000)x(i),y(i),z(i),px(i),py(i),pz(i),em(i),iso(i),
     &isy(i),id(i)
        bpz=bpz+pz(i)/npr                                                       
        apz=apz+z(i)/npr                                                        
        bpx=bpx+px(i)/npr                                                       
        apx=apx+x(i)/npr                                                        
        bpy=bpy+py(i)/npr                                                       
        apy=apy+y(i)/npr                                                        
 10     continue                                                                

      close(30)

        do 11 i=1,npr                                                           
        pz(i)=pz(i)-bpz
        z(i)=z(i)-apz
        px(i)=px(i)-bpx
        x(i)=x(i)-apx
        py(i)=py(i)-bpy                                                         
        y(i)=y(i)-apy                                                           
        id(i)=1                                                                 
      if (i > 0.and.i <= nzpr)isy(i)=-1
      if (i > nzpr.and.i <= npr)isy(i)=1      
  11    continue                                                                
      end if

        p0=-tpo                                                                 
                                                    
      if (nta.ne.0.and.nzta.ne.0) then

      open(10,file=targ,form='formatted',status='old')

        do 20 i=npr+1,np                                                        
        read (10,1000)x(i),y(i),z(i),px(i),py(i),pz(i),em(i),iso(i),
     &isy(i),id(i)
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
      if (i > npr.and.i <= npr+nzta)isy(i)=-1
      if (i > npr+nzta.and.i <= np)isy(i)=1
21      continue                                                                
      end if

      goto 767  !fengzhq
C      do 91410 in=1,nru-1
C      if (aini(nru) == aini(in)) then
C      iangle=1
C      go to 91411
C      else
C      iangle=0
C      end if
C91410      continue
C91411      if (iflag1 < 0.and.iangle == 1) then
      cth=0.   !1.-2.*ran(mseed)
      sth=1.   !sqrt(1.-cth**2.)
      phi=6.28*ran(mseed)
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
1      continue
      print*,'rotate proj. two euler angles!'
C      end if
C        do 3 i=1,npr                                                           
C        pz(i)=pz(i)+pzpr
C        z(i)=z(i)-z0                                                    
C        px(i)=px(i)+pxpr
C        x(i)=x(i)+x0                                                        
C        py(i)=py(i)
C        y(i)=y(i)
C        id(i)=1
C3      continue
C      do 81410 ln=1,nru-1
C      if (bini(nru) == bini(ln)) then
C      langle=1
C      go to 81411
C      else
C      langle=0
C      end if
C81410      continue
C81411      if (iflag1 < 0.and.langle == 1) then
      cth=0.   !1.-2.*ran(mseed)
      sth=1.   !sqrt(1.-cth**2.)
      phi=6.28*ran(mseed)
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
10004      continue
      print*,'rotate targ. two euler angles!'
C      end if
C        do 4  i=npr+1,np                                                       
C        pz(i)=pz(i)-pzta
C        z(i)=z(i)+z0                                                    
C        px(i)=px(i)-pxta
C        x(i)=x(i)-x0                                                        
C        py(i)=py(i)
C        y(i)=y(i)
C        id(i)=51                                                                
C4      continue
c
        goto 767    ! fengzhq                                                                         
      if ((iflag1 < 0.and.initi == 1).or.iflag1 > 0) then
      do 3269 i=1,np
      write (71,911)z(i),x(i)
911      format(1x,2f9.4)
3269      continue      
      end if
      dt=dt/2.                                                                
c      
c       added 21.1.87                                                           
c                                                                               
        if (iflag5 == 0) goto 765                                                  
c       print*,'im in updatev'                                                  
        factor=4.*cdelt*cmeps                                                 
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
        vdens=xexp(-dr2/al)                                                      
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr/sal                                     
        dvdp1=-anorm/widthp**2*xexp(-dr2/wq2)*xexp(-p2/wp2)                       
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
                                                                                
765     do 8014 i=1,np                                                          
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1                                                   
8014    continue                                                                
766     dt=dt*2.    
767      continue                                                            
1000    format(1x,7f9.4,2x,2i3,2x,i3) 
10000   format(1x,i20,1x,i4,4(f8.2,1x))                                         
10001   format(1x,4(i3,1x))                                                     
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

              subroutine init1

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/rp0/xi(5000),yi(5000),zi(5000),pxi(5000),pyi(5000),
     @ pzi(5000),emi(5000),isoi(5000),isyi(5000),idi(5000)
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1, 
     +ncolln,ncollp,nblkn,nblkp
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr,nsym                                                     
      common/ee/mseed,nseed
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/densi1/denstl(500),dens1(500)                                                
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/tt/prtime(0:8000),nruns,nt,rtime,perc
      common/max/mden
      common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/riclw/rrrr(500)
      common/derot/alxyz,idelrot
      common/cmclw/apx,apy,apz,atx,aty,atz
      COMMON/INICLW/DISTNN,DISTNP,PHMANN,PHMANP
      COMMON/WORK1/ Mnumn,Mnump,Rmaxn,Rmaxp,RMS0N,RMS0P
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/bind/progEb,aene,pfc,laesw
      common/deform/beta2
      dimension dvdpx(500),dvdpy(500),dvdpz(500)                              
      dimension tx(500),ty(500),tz(500),rr(500)
      real mden

      if (NTA > 0) then
        DISTNN=1.5
        DISTNP=0.5
        PHMANN=0.28**2                                                           
        PHMANP=0.28**2/100.
      end if    
                                                            
!98513 continue
!LW9778      continue
!98513if (NTA <= 0.AND.NPR > 0) THEN
!       CALL RADRI
!     END IF
98513 continue
      IFSPSP=0
      IFSPSN=0
      IFMOM=0
      lloss=1
      iclwfg=0

c     distance between surfaces = 5 fm                                        
      hc = 1.24                                                        
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
      iimis=0
c96111   do 1 i=1,nzpr
c13      cth=1.-2.*ran(mseed)                        
c        sth=sqrt(1.-cth**2)                                                     
c        phi=6.28*ran(mseed)                
c        z(i)=rrrr(i)*cth/gami-z0                                          
c        x(i)=rrrr(i)*sth*cos(phi)+x0
c        y(i)=rrrr(i)*sth*sin(phi)
c        iso(i)=1                                                                
c        isy(i)=-1
c        id(i)=1                                                           
c        do 12 j=1,i-1                                                           
c        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
c      if (ISO(I) == ISO(J)) THEN
c        if (d < DISTNN) THEN
c        IFSPSP=IFSPSP+1
c        if (IFSPSP > 1000000) GOTO 98513
c        go to 13    
c        END IF      
c      ELSE                
c        if (d < DISTNP)  THEN  
c        IFSPSP=IFSPSP+1
c        if (IFSPSP > 1000000) GOTO 98513
c        go to 13
c        END IF
c      END IF
c12      continue                                                                
c1       continue                                                                
c
c        do 11 i=nzpr+1,npr
c155     cth=1.-2.*ran(mseed)          
c        sth=sqrt(1.-cth**2)                                                     
c        phi=6.28*ran(mseed)           
c        z(i)=rrrr(i)*cth/gami-z0                                       
c        x(i)=rrrr(i)*sth*cos(phi)+x0
c        y(i)=rrrr(i)*sth*sin(phi)
c        id(i)=1 
c        iso(i)=0                                                                
c        isy(i)=1
c        do 14 j=1,i-1                                                           
c        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
c      if (ISO(I) == ISO(J)) THEN
c        if (d < DISTNN) THEN    
c        IFSPSN=IFSPSN+1
c        if (IFSPSN > 1000000) GOTO 98513
c        go to 155
c        END IF
c      ELSE                
c        if (d < DISTNP) THEN
c        IFSPSN=IFSPSN+1
c        if (IFSPSN > 1000000) GOTO 98513
c        go to 155    
c        END IF
c      END IF
c14      continue                                                                
c11      continue    
            
c       sampling using rectangular shape distribution   fengzhq
cc-------------------------------------------------------------->
      pi=3.1415927

      r1=1.28*float(npr)**0.3333-0.76+0.8*float(npr)**(-0.3333)
      beta20=beta2*sqrt(5./(4.*pi))

      ! Projectile proton sampling
      do
        rms = 0.
        do i1 = 1, nzpr
          rr(i1) = r1 * ran(nseed)**0.3333
          rms = rms + rr(i1) / float(nzpr)
        end do
        ! Check after sampling, exit loop if condition is met
        if (abs(rms - rms0p) <= 0.3) exit
      end do
      ! Projectile neutron sampling
      do
        rms = 0.
        do i2 = nzpr + 1, npr
          rr(i2) = r1 * ran(nseed)**0.3333
          rms = rms + rr(i2) / float(npr - nzpr)
        end do 
        if (abs(rms - rms0n) <= 0.3) exit
      end do

        do 1 i=1,nzpr
13      cth=1.-2.*ran(mseed)                        
        sth=sqrt(1.-cth**2)
        phi=6.28*ran(mseed)
        rr(i)=rr(i)*(1.+beta20*(3.*cth**2.-1.)/2.)                
        z(i)=rr(i)*cth                                 
        x(i)=rr(i)*sth*cos(phi)
        y(i)=rr(i)*sth*sin(phi)
        id(i)=1   
        iso(i)=1
        isy(i)=-1
        do 12 j=1,i-1
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
      if (ISO(I) == ISO(J)) THEN
        if (d < DISTNN) then
        IFSPSP=IFSPSP+1
        if (IFSPSP > 100000) goto 98513
        goto 13
        endif          
      ELSE                
        if (d < DISTNP) then
        IFSPSP=IFSPSP+1
        if (IFSPSP > 100000) goto 98513
        goto 13
        endif
      END IF
12      continue
1       continue

        do 11 i=nzpr+1,npr
155     cth=1.-2.*ran(mseed)          
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*ran(mseed)   
        rr(i)=rr(i)*(1.+beta20*(3.*cth**2.-1.)/2.)         
        z(i)=rr(i)*cth                                      
        x(i)=rr(i)*sth*cos(phi)
        y(i)=rr(i)*sth*sin(phi)
        id(i)=1 
        iso(i)=0                                                                
        isy(i)=1
        do 14 j=1,i-1                                                           
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
      if (ISO(I) == ISO(J)) THEN
        if (d < DISTNN) then
        IFSPSN=IFSPSN+1
        if (IFSPSN > 100000) goto 98513
        goto 155
        endif
      ELSE                
        if (d < DISTNP) then
        IFSPSN=IFSPSN+1
        if (IFSPSN > 100000) goto 98513      
        goto 155
        endif    
      END IF
14      continue                                                                
11      continue                                

96112   do 111 i=npr+1,npr+nzta  
19      r11=(rrt+0.2*rrt)*ran(mseed)
        r12=(rrt+0.2*rrt)*ran(mseed)       
        r13=(rrt+0.2*rrt)*ran(mseed)       
        r=amax1(r11,r12,r13)                                                    
        cth=1.-2.*ran(mseed)     
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*ran(mseed)
        z(i)=(r*cth)/gami+z0                                         
        x(i)=r*sth*cos(phi)-x0
        iso(i)=1                                                                
        isy(i)=-1
        id(i)=51                                                                
        y(i)=r*sth*sin(phi)                                                     
        do 18 j=npr+1,i-1                                                       
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
      if (ISO(I) == ISO(J)) THEN
        if (d < DISTNN)goto 19    
      ELSE                
        if (d < DISTNP)goto 19    
      END IF
18      continue                                                                
111     continue                                                                

        do 1111 i=npr+nzta+1,np                                                 
17      r11=(rrt+0.2*rrt)*ran(mseed)        
        r12=(rrt+0.2*rrt)*ran(mseed)
        r13=(rrt+0.2*rrt)*ran(mseed) 
        r=amax1(r11,r12,r13)                                                   
        cth=1.-2.*ran(mseed)
        sth=sqrt(1.-cth**2)                                                     
        phi=6.28*ran(mseed)
        z(i)=(r*cth)/gami+z0                                       
        x(i)=r*sth*cos(phi)-x0
        iso(i)=0                                                                
        isy(i)=1
        id(i)=51                                                                
        y(i)=r*sth*sin(phi)                                                     
        do 16 j=npr+1,i-1                                                       
        d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)                    
      if (ISO(I) == ISO(J)) THEN
        if (d < DISTNN)goto 17
      ELSE
        if (d < DISTNP)goto 17
      END IF                                                
16      continue                                                                
1111    continue                                                                
C
      if (nta > 0) then
      do ii=npr+1,np
      tx(ii)=0.
      ty(ii)=0.
      tz(ii)=0.
      txx=0.
      tyy=0.
      tzz=0.
      end do
      do ii=npr+1,np
      txx=txx+x(ii)
      tyy=tyy+y(ii)
      tzz=tzz+z(ii)
      end do
      do ii=npr+1,np
      tx(ii)=x(ii)-txx/float(nta)
      ty(ii)=y(ii)-tyy/float(nta)
      tz(ii)=z(ii)-tzz/float(nta)
      end do
      dddd=0.0
      do i=npr+1,np
      dddd=dddd+tx(i)*tx(i)+ty(i)*ty(i)+tz(i)*tz(i)
      end do
      ddddd=sqrt(dddd/float(nta))
      if (abs(ddddd-rrt) > 0.05*rrt) then
      iimis=iimis+1
C
      if (iimis > 1000) stop 'Reinput parameters RRT0 to initialize T'
C
      goto 96112
      end if
      end if
C
c      sampling in momentum space

      call densitn  

      p0z=0.
      p0x=0.                                                           
c      p0z=pzpr
c      p0x=pxpr
c        if (npr.ne.0)dpp=(0.008/float(npr)*3.4)**0.666
c        if (nta.ne.0)dpt=(0.008/float(nta)*3.4)**0.666
c        if (npr.ne.0.and.npr < 30)dpp=(0.008/float(npr)*4.3)**0.666
c        if (nta.ne.0.and.nta < 30)dpt=(0.008/float(nta)*4.3)**0.666
c      dpc=dpp                                                                 
c      ALFA=FLOAT((NP-2*(NZPR+NZTA)))/FLOAT(NP)

351      continue
      do 36 i=1,np
        densty=dens1(i)                                                      
        ddt=ddt+dens1(i)/float(np)                                          
c       aux=(3.*densty)/(16.*3.14159)    
        aux=(3.*densty)/(8.*pi)                                                 
        if (aux < 0)aux=0.                                                      
        pfl=hc*aux**.3333 
c352   if (i <= nzpr.OR.(I > NPR.AND.I <= (NPR+NZTA))) THEN
cc      pr=(1-alfa/3.)*(pfl-pfc)*ran(mseed)**.3333 
cc      else
cc      pr=(1+alfa/3.)*(pfl-pfc)*ran(mseed)**.3333 
cc      pr=(pfl-pfc)*ran(mseed)**.3333
cc      END IF
352      pr=(pfl-pfc)*ran(mseed)**.3333

        if (pr < 0)pr=0.                 
        cth=1.-2.*ran(mseed)            
        sth=sqrt(1.-cth*cth)                                                    
        phi=2.*pi*ran(mseed)          
cc        if (i.ne.1+npr) goto 35                                                  
cc        p0z=-pzta
cc        p0x=-pxta
cc        dpc=dpt                                                                 
35      pz(i)=p0z+pr*cth                                                        
        px(i)=p0x+pr*sth*cos(phi)
        py(i)=pr*sth*sin(phi)  
                                                   
        do 357 j=1,i-1                                                          
        dp=(px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2                   
        d=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2                          
      if (ISO(I) == ISO(J)) THEN
        if (dp*d < PHMANN) THEN
        IFMOM=IFMOM+1
        if (IFMOM >= 1000000.AND.NPR > 0)      GOTO 98513
        goto 352
        END IF
      ELSE
        if (dp*d < PHMANP) then 
        IFMOM=IFMOM+1
        if (IFMOM >= 1000000.AND.NPR > 0) GOTO 98513
        goto 352
        END IF
      END IF                                 
357     continue   
                                                             
        pav=pav+pr/float(np)                                                    
        em(i)=0.938                                                        
36         continue                                                                

      if (iclwfg == 1) then
      alxyz=10.
      goto 97710
clw      goto 96111
      end if


      if (idelrot <= 0) goto 97710

      call delrot
      if (lloss > 100) then
      write (*,*)'Lxyz=',alxyz
      write (*,*) 'The rotation cannot be gotten rid of in this
     & running because Lloss>100 !!'
      write (*,*)
      iclwfg=1
      goto 98513
      end if

CLW      do i=1,npr
      do i=1,NP
        do j=1,i-1 

        dp=(px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2                   
        d=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2                          

      if (ISO(I) == ISO(J)) THEN
       if (dp*d < PHMANN) then
       lloss=lloss+1
       goto 351
       end if
      ELSE
       if (dp*d < PHMANP) then
       lloss=lloss+1
       goto 351
       end if
      END IF

        end do
      end do
C*

97710      continue


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
      DO I=1,NP
        if (id(i) == 1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        endif                                                                   
        if (id(i) == 51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        endif 
      END DO

clwnew       if (nru == 1)write (54,453)ddt,pav                                        
453     format(1x,'average density= ',f10.4,' average fermi momentum            
     +  = ',f10.4)                                                              

        do 15 i=1,np                                                             
        if (npr == 0.or.nta == 0) then
            z0=0.
            x0=0.
            pzpr=0.
            pxpr=0.
            pzta=0.
            pxta=0.
        end if 

        if (id(i) == 1) then                                                     
                      pz(i)=pz(i)-bpz+pzpr
                      z(i)=z(i)-apz-z0
                      px(i)=px(i)-bpx+pxpr                                    
                      x(i)=x(i)-apx+x0
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
      if (iflag1 < 0) then
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
      end if
        end if                                                                  

        if (id(i) == 51) then                                                    
                       pz(i)=pz(i)-btz-pzta
                       z(i)=z(i)-atz+z0
                       px(i)=px(i)-btx-pxta                                   
                       x(i)=x(i)-atx-x0
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
      if (iflag1 < 0) then
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
      end if 
        end if 
77777      format(1x,7f9.4,2x,2i3,2x,i3)      

15      continue
                                                                
345     format(1x,'average excitation energy = ', f10.4)                        
      if ((iflag1 < 0.and.initi == 1).or.iflag1 > 0) then
      write (71,*)'time=',rtime
      do 98721 i=1,np
      write (71,87731)z(i),x(i)
87731      format(1x,2f9.4)      
98721      continue
      end if
c      
      dt=dt/2.                                                                       
c      added 21.1.87                                                           
        if (iflag5 == 0)goto 765                                                  
c       print*,'im in updatev'                                                  
        factor=4.*cdelt*cmeps                                                   
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
        vdens=xexp(-dr2/al)                                                      
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr/sal                                     
        dvdp1=-anorm/widthp**2*xexp(-dr2/wq2)*xexp(-p2/wp2)                       
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
c       end of changes                                                          
                                                                                
765     do 8014 i=1,np                                                          
        e1=sqrt(em(i)*em(i)+px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))                
        x(i)=x(i)+dt*px(i)/e1                                                   
        y(i)=y(i)+dt*py(i)/e1                                                   
        z(i)=z(i)+dt*pz(i)/e1                                                   
8014    continue                                                      
766     dt=dt*2.
                                                              
cc      write (*,*)'Lloss=',lloss
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////
              subroutine init0 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  
      common/bind/progEb,aene,pfc,laesw
      common/skyrme/al,sal,t0,t1,arang,vy,width      

        integer iter_count
        hb=0.197327
      pfc0=0.92*sqrt(-2.*log(0.5)*hb**2./4/width)   !0.025
      pfc =pfc0
        iter_count = 0

10        iter_count = iter_count + 1
        if (iter_count  >  10000) then
              print *, 'FATAL: init0 >10000 iters, aborting'
              stop
        endif
        call init1      
        call prot 

        laesw=0
        call avene(avpot,aveki,Esh,1)    
        Eb=aene
        Ebmin=progEb-0.2
        Ebmax=progEb+0.2

        if (Eb < Ebmin) then
        pfc=pfc-0.002
        goto 10
        endif
        if (Eb > Ebmax) then
        pfc=pfc+0.002
        goto 10
        endif

        laesw=1
        call avene(avpot,aveki,Esh,1)

      return
      end
clw\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

                  subroutine delrot

clw\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp 
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/derot/alxyz,idelrot

            NP=NPR+NTA
            FNP=FLOAT(NP)
          xxr1=0.
          yyr1=0.
          zzr1=0.
          xxp1=0.
          yyp1=0.
          zzp1=0.
          xxl1=0.
          yyl1=0.
          zzl1=0.
c
          rixx=0.
          riyy=0.
          rizz=0.
          rixy=0.
          rixz=0.
          riyz=0.
c
         do icm1=1,NP
          xxr1=xxr1+ X(icm1)
          yyr1=yyr1+ Y(icm1)
          zzr1=zzr1+ Z(icm1)
          xxp1=xxp1+PX(icm1)
          yyp1=yyp1+PY(icm1)
          zzp1=zzp1+PZ(icm1)
          xxl1=xxl1+ Y(icm1)*Pz(icm1)-Z(icm1)*Py(icm1)
          yyl1=yyl1+ Z(icm1)*Px(icm1)-X(icm1)*Pz(icm1)
          zzl1=zzl1+ X(icm1)*Py(icm1)-Y(icm1)*Px(icm1)
c
          rixx=rixx +Y(icm1)**2.+Z(icm1)**2.
          riyy=riyy +Z(icm1)**2.+X(icm1)**2.
          rizz=rizz +X(icm1)**2.+Y(icm1)**2.
          rixy=rixy -X(icm1)*Y(icm1)
          rixz=rixz -X(icm1)*Z(icm1)
          riyz=riyz -Y(icm1)*Z(icm1)
         end do

         xxr1=xxr1/FNP
         yyr1=yyr1/FNP
         zzr1=zzr1/FNP
         xxp1=xxp1/FNP
         yyp1=yyp1/FNP
         zzp1=zzp1/FNP
         xxl1=xxl1/FNP -(yyr1*zzp1-zzr1*yyp1)
         yyl1=yyl1/FNP -(zzr1*xxp1-xxr1*zzp1)
         zzl1=zzl1/FNP -(xxr1*yyp1-yyr1*xxp1)
c
         rixx=rixx/FNP -yyr1**2.-zzr1**2.
         riyy=riyy/FNP -zzr1**2.-xxr1**2.
         rizz=rizz/FNP -xxr1**2.-yyr1**2.
         rixy=rixy/FNP +xxr1*yyr1 
         rixz=rixz/FNP +xxr1*zzr1
         riyz=riyz/FNP +yyr1*zzr1
c
         DET=rixx*riyy*rizz + 2.*rixy*rixz*riyz
     +        -rixx*riyz**2. -riyy*rixz**2. -rizz*rixy**2.
         WX=(xxl1 * ( RIYY*RIZZ-RIYZ**2   )+
     +         yyl1 * ( RIXZ*RIYZ-RIZZ*RIXY )+
     +         zzl1 * ( RIXY*RIYZ-RIYY*RIXZ ) )/DET
         WY=(xxl1 * ( RIXZ*RIYZ-RIXY*RIZZ )+
     +         yyl1 * ( RIXX*RIZZ-RIXZ**2   )+
     +         zzl1 * ( RIXY*RIXZ-RIXX*RIYZ ) )/DET
         WZ=(xxl1 * ( RIXY*RIYZ-RIYY*RIXZ )+
     +         yyl1 * ( RIXZ*RIXY-RIXX*RIYZ )+
     +         zzl1 * ( RIXX*RIYY-RIXY**2   ) )/DET
         xxl1=xxl1*FNP
         yyl1=yyl1*FNP
         zzl1=zzl1*FNP

         do iran=1,NP
         X(iran)= X(iran)-xxr1
         Y(iran)= Y(iran)-yyr1
         Z(iran)= Z(iran)-zzr1
         PX(iran)=PX(iran)-xxp1-(WY*Z(iran)-WZ*Y(iran))
         PY(iran)=PY(iran)-yyp1-(WZ*X(iran)-WX*Z(iran))
         PZ(iran)=PZ(iran)-zzp1-(WX*Y(iran)-WY*X(iran))
         enddo
c
      xxl1=0.
      yyl1=0.
      zzl1=0.
      do icm1=1,NP
          xxl1=xxl1+ Y(icm1)*Pz(icm1)-Z(icm1)*Py(icm1)
          yyl1=yyl1+ Z(icm1)*Px(icm1)-X(icm1)*Pz(icm1)
          zzl1=zzl1+ X(icm1)*Py(icm1)-Y(icm1)*Px(icm1)
      end do

      xxl1=xxl1/fNP/0.197
      yyl1=yyl1/fNP/0.197
      zzl1=zzl1/fNP/0.197
      alxyz=sqrt(xxl1*xxl1+yyl1*yyl1+zzl1*zzl1)

      RETURN
      END

c///////////////////////////////////////////////////////////////////////////////

                  subroutine densitn                                          

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                      
        common/densi1/denstl(500),dens1(500)                                   
        common/wavel/sl(500)  
        pi=3.1415927                                                         
        do 10 i=1,np   
        denstl(i)=0.
        dens1(i)=0.0                                               
10      continue                                                   
        do 20 i=1,np                                               
        do 30 j=1,np   
        width=sl(i)+sl(j)                                           
        dr2=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2           
        dd=xexp(-dr2/width)                                         
        denstl(i)=denstl(i)+dd/(pi*width)**1.5
        if (iso(i) == iso(j)) dens1(i)=dens1(i)+dd/(pi*width)**1.5
30           continue                                                                     
20           continue                                                   
        return                                                     
        end                                                        

c///////////////////////////////////////////////////////////////////////////////

              subroutine densim(densm)  
 
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        dimension densm(500)        
                                                  
        widthm=.1973289**2*2./al                                          
!        anorm=1./(3.1415927*widthm)**1.5/np                                                                                         
        anorm=1./(3.1415927*widthm)**1.5/np**2    !fengzhq
        do 10 i=1,np                                                          
10      densm(i)=0.0                                                          
        do 20 i=1,np                                                          
        do 30 j=i+1,np                                                     
        dp2=(px(i)-px(j))**2+(py(i)-py(j))**2+(pz(i)-pz(j))**2                
        dd=xexp(-dp2/widthm)                                                   
        densm(i)=densm(i)+dd                                                  
        densm(j)=densm(j)+dd                                                                                             
30      continue                                                              
        densm(i)=anorm*(1.+densm(i))                                                      
20      continue                                                              
        return                                                                
        end                                                                   

c///////////////////////////////////////////////////////////////////////////////

              subroutine dendgn 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1,
     +      ncolln,ncollp,nblkn,nblkp
        common/bb/np,ntag,x00
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
        if (lx < 3.or.lx > 47)goto1                                  
        if (ly < 3.or.ly > 47)goto1                                  
        if (lz < 3.or.lz > 47)goto1                                  
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

              subroutine prot

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      character idate*9,itime*8        
c       common/densi/err(302),dens(50,50,50)                                    
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     +  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/rpi/rmst(5000),pmst(5000),pmsn(5000),pmsz(5000),avet(5000)
      common/aa/iflag1,iflag2,iflag3,iflag4,ncoll,nblock,ncoll1, 
     +      ncolln,ncollp,nblkn,nblkp
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr,nsym
      common/wavel/sl(500)      
      common/poli/iflagcro,iflagpau,ipoli                     
      common/tt/prtime(0:8000),nruns,nt,rtime,perc             
      common/cmod/coeft,coefp,epp,ecm,ptcm
      COMMON /WORK1/ Mnumn,Mnump,Rmaxn,Rmaxp,RMS0N,RMS0P
      common/riclw/rrrr(500)
      dimension ipuli1(0:300),ipuli2(0:300),ipulit(0:300),
     +      ratio1(0:300),
     +      iplin1(0:300),iplin2(0:300),iplin(0:300),ration(0:300),
     +      iplip1(0:300),iplip2(0:300),iplip(0:300),ratiop(0:300)
c      real trms,prms
        hb=0.197327   
                                                         
        pp=0.                                                                   
        ppn=0.
        ppz=0.
        tp=0.                                                                   
        tpn=0.
        tpz=0.
        pr=0.                                                                   
        tr=0.                                                                   
        atz=0                                                                   
        apz=0                                                                   
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
        do 136 i=1,np                                                           
        if (id(i) < 50) then                                                    
                   bpz=bpz+pz(i)/npr                                            
                   apz=apz+z(i)/npr                                             
                   bpx=bpx+px(i)/npr                                            
                   apx=apx+x(i)/npr                                             
                   bpy=bpy+py(i)/npr                                            
                   apy=apy+y(i)/npr                                             
        end if 
        if (id(i) > 50) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        end if 
136     continue                                                                
        bcmx=(apx*npr+atx*nta)/float(np)                                        
        bcmy=(apy*npr+aty*nta)/float(np)                                        
        bcmz=(apz*npr+atz*nta)/float(np)                                          
        do 137 i=1,np                                                          
        if (id(i)  <  50.and.nta == 0.and.epp == 0.) then
               pr = pr + (x(i)-apx)**2.+(y(i)-apy)**2.+(z(i)-apz)**2.        
               pp = pp + (px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.     
         if (iso(i) == 0) then
           ppn=ppn+(px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.
         ELSE
           ppZ=ppZ+(px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.
         END IF
        end if 
        if (id(i)  >  50.and.npr == 0.and.epp == 0.) then
               tr = tr + (x(i)-atx)**2.+(y(i)-aty)**2.+(z(i)-atz)**2.        
               tp = tp + (px(i)-btx)**2.+(py(i)-bty)**2.+(pz(i)-btz)**2.     
         if (iso(i) == 0) then
           tpn=tpn+(px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.
         ELSE
           tpZ=tpZ+(px(i)-bpx)**2.+(py(i)-bpy)**2.+(pz(i)-bpz)**2.
         END IF

        end if               
137     continue                                                                

      llt=ntstep/int(1./dt+0.0001)+1
                                                
      if (npr > 0.and.nta == 0.and.epp == 0.) then
      pr=pr/float(npr)+2.*sl(npr)
cc      pr=pr/float(npr)
      pr = sqrt(pr)
      rmst(llt)=pr
      rms=pr  
      prtime(ntstep) = prtime(ntstep) + pr                                  
      end if

      if (nta > 0.and.npr == 0.and.epp == 0.) then
        tr=tr/float(nta)
        tr = sqrt(tr)
      rmst(llt)=tr
      end if
      if (npr > 0.and.nta == 0.and.epp == 0.) then
      pp=pp/float(npr)
      ppN=ppN/float(npr-nzpr)
      ppZ=ppZ/float(nzpr)
        pp=sqrt(pp)
        ppN=sqrt(ppN)
        ppZ=sqrt(ppZ)
      pmst(llt)=pp
      pmsN(llt)=ppN
      pmsZ(llt)=ppZ
      end if
      if (nta > 0.and.npr == 0.and.epp == 0.) then
      tp=tp/float(nta)
      tpN=tpN/float(nta-nzta)
      tpZ=tpZ/float(nzta)
        tp = sqrt(tp)
        tpN=sqrt(tpN)
        tpZ=sqrt(tpZ)
      pmst(llt)=tp
      pmsN(llt)=tpN
      pmsZ(llt)=tpZ
      end if
cc*************************************************************
      if (iflag1 < 0.and.epp == 0.and.rtime == 0.)
     +radj=(rms0n*(npr-nzpr)+rms0p*nzpr)/npr
      if (iflag1 < 0.and.epp == 0.and.abs(rms-radj) > rmax) then
      lswitch=0
      else
      lswitch=1
      end if
      do i=1,np
        rrrr(i)=sqrt((x(i)-bcmx)**2+(y(i)-bcmy)**2+(z(i)-bcmz)**2)
      end do
      goto 26                                      
      if (iflag1 < 0.and.iprin == 1.and.npr > 0.and.nta == 0) then
      write (*,99)
      write (*,98)pr,pp,perc,idate,itime
      end if
      if (iflag1 < 0.and.iprin == 1.and.nta > 0.and.npr == 0) then
      write (*,98)tr,tp,perc,idate,itime
      end if
      if (epp > 0.and.nta > 0.and.npr > 0.and.iprin == 1) then
      write (*,98)rtime,rms,perc,idate,itime
      end if
98      format(5x,f9.4,6x,f9.4,5x,f9.4,'%',7x,a9,6x,a8)
99      format(6x,'rms radius',5x,'pms radius',5x,'completed',8x,'date',
     &11x,'time')
62      format(1x,2f9.4)
222     format(////1x,'  ****  run # = ',i3,' timestep # = ',                   
     &  i5,' real time = ',f10.4,' fm/c   *****')                               
145     format(1x,'pz=',f8.3,' ppz=',f8.3,' tz=',f8.3,' tpz=',f8.3)             
156     format(1x,'py=',f8.3,' ppy=',f8.3,' ty=',f8.3,' tpy=',f8.3)             
165     format(1x,'px=',f8.3,' ppx=',f8.3,' tx=',f8.3,' tpx=',f8.3)             
166     format(1x,'center of mass,x=,y=,z=',3f8.3,'rms=',1f8.3)                 
181     format(1x,'root mean square values')                                    
cc                                                                              
      GOTO 98512
      if (nru == 1.and.rtime == 0.) then
       write (55,*) ' \\E= ',epp,' MeV/u'
       write (55,*) ' \\b= ',x00,' fm'
      if (iflagcro <= 0) then
      write (55,*) ' \\Cugnon''s N-N cross section !!'
      else if (iflagcro >= 1) then
      write (55,*) ' \\Experimental N-N cross section !!'
      end if 
      if (ipoli == 0.or.ipoli == 1) then
      if (ipoli == 0) then
      write (55,*) ' \\All collisions are blocked !!'
      end if
      if (ipoli == 1) then
      write (55,*) ' \\No Pauli blocking !!'
      end if
      GOTO 96322
      end if

       if (iflagpau <= 0) write (55,*) ' \\Old Pauli blocking !!'
       if (iflagpau == 1) write (55,*) ' \\R34 Pauli blocking !!'
       if (iflagpau == 2) write (55,*) ' \\Over Pauli blocking !!'
       if (iflagpau == 3) write (55,*) ' \\Rho Pauli blocking !!'
96322  continue
       write (55,*) ' \\Time(fm/c) Col.Numb Blk.Numb RealColl BC.Ratio
     + rcolltot'
      end if
C--------------------------------------------------------------------
      iplin1(int(rtime))=iplin1(int(rtime))+ncolln
      iplin2(int(rtime))=iplin2(int(rtime))+nblkn
      iplin(int(rtime))=iplin(int(rtime))+iplin3
      iplip1(int(rtime))=iplip1(int(rtime))+ncollp
      iplip2(int(rtime))=iplip2(int(rtime))+nblkp
      iplip(int(rtime))=iplip(int(rtime))+iplip3
96515      format(1x,f10.2,2(1x,i8),1x,f8.4,2x,f7.5,2x,f9.4)
      ipuli1(int(rtime))=ipuli1(int(rtime))+ncoll
      ipuli2(int(rtime))=ipuli2(int(rtime))+nblock
      ipulit(int(rtime))=ipulit(int(rtime))+ipuli3
C
      if (nru == nruns) then
      if (ipuli1(int(rtime)) == 0) then
      ratio1(int(rtime))=1.
      else
      ratio1(int(rtime))=float(ipuli2(int(rtime)))/
     +      float(ipuli1(int(rtime)))
      end if
C
      if (iplin1(int(rtime)) == 0) then
      ration(int(rtime))=1.
      else
      ration(int(rtime))=float(iplin2(int(rtime)))/
     +      float(iplin1(int(rtime)))
      end if
C
      if (iplip1(int(rtime)) == 0) then
      ratiop(int(rtime))=1.
      else
      ratiop(int(rtime))=float(iplip2(int(rtime)))/
     +       float(iplip1(int(rtime)))
      end if
C
      write (55,96515) rtime,ipuli1(int(rtime)),ipuli2(int(rtime)),
     +float(ipuli1(int(rtime))-ipuli2(int(rtime)))/float(nru),
     +ratio1(int(rtime)),
     +float(ipulit(int(rtime)))/float(nru)
C
      write (56,96515) rtime,iplin1(int(rtime)),iplin2(int(rtime)),
     +float(iplin1(int(rtime))-iplin2(int(rtime)))/float(nru),
     +ration(int(rtime)),
     +float(iplin(int(rtime)))/float(nru)
C
      write (57,96515) rtime,iplip1(int(rtime)),iplip2(int(rtime)),
     +float(iplip1(int(rtime))-iplip2(int(rtime)))/float(nru),
     +ratiop(int(rtime)),
     +float(iplip(int(rtime)))/float(nru)
      END if

98512      CONTINUE                                                    
cc                                                                              
        if (iflag4 == 2) goto 26                                                  
        pions=0.                                                                
        do 950 i=1,np                                                           
        if (abs(em(i)-0.938) > .01)pions=pions+1                                
950     continue                                                                                         
c       write (6,9991)dtg(ideu),dtg1(ideu)                                       
9991    format(1x,'no of d this run = ',f10.4,'  no of prod d this              
     1  run = ',f10.4)   
                                                                                                  
26      continue       
                                                         
        return                                                                  
        end                                                                     
c///////////////////////////////////////////////////////////////////////////////

                  subroutine updt 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        dimension f1x(500),f1y(500),f1z(500),t2(500),dr2(500)                   
c                                                                               
        do 1 i=1,np                                                             
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
        do 21 i2=i1+1,np                                                        
        dr2(i2)=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))         
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        t2(i2)=xexp(-dr2(i2)/al)                                                 
21      continue                                                                
        do 22 i2=i1+1,np                                                        
        f1x(i1)=f1x(i1)+(x(i1)-x(i2))*t2(i2)                                    
        f1x(i2)=f1x(i2)+(x(i2)-x(i1))*t2(i2)                                    
        f1y(i1)=f1y(i1)+(y(i1)-y(i2))*t2(i2)                                    
        f1y(i2)=f1y(i2)+(y(i2)-y(i1))*t2(i2)                                    
        f1z(i1)=f1z(i1)+(z(i1)-z(i2))*t2(i2)                                    
22      f1z(i2)=f1z(i2)+(z(i2)-z(i1))*t2(i2)                                    
        sum=0.                                                                  
        do 23 i2=i1+1,np                                                        
        sum=sum+t2(i2)                                                          
        fr(i2)=fr(i2)+t2(i2)                                                    
23      continue                                                                
        fr(i1)=fr(i1)+sum                                                       
2       continue                                                                
        do 3 i=1,np                                                             
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

              subroutine updtyc  

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/pot/pota,potb,potc,ih(10000),potd,pote,iflag5,potf               
     &        ,potg,pots,xmu,asym,bsym,csym,gamsy
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/coul/cy(2000),cyu(2000)                                          
        common/fdelt/ifdelt,fdelta,alnew,salnew,potdnew,potenew,
     +      frtwo(500)
        dimension f1x(500),f1y(500),f1z(500),t2(500),t2new(500) 
        dimension fcx(500),fcy(500),fcz(500),u(500)                             
c       dimension fhx(500),fhy(500),fhz(500),upau(i2)                           
        dimension fsx(500),fsy(500),fsz(500)
c      dimension s(500)
        do 1 i=1,np                                                             
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
        do 21 i2=i1+1,np                                                        
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        t2(i2)=xexp(-dr2/al)                                                     
      if (ifdelt >= 1) t2new(i2)=xexp(-dr2/alnew)
        idr22=dr2*25.+1.                                                        
        idr2=min0(idr22,10000)                                                  
        nccc=ih(idr2)                                                           
        u(i2)=iso(i1)*iso(i2)*cy(nccc)+cyu(nccc)
c       upau(i2)=cyu(nccc)*xexp(-p2**2/width)                                    
21      continue                                                                
        do 22 i2=i1+1,np                                                        
        t3=(x(i1)-x(i2))*u(i2)                                                  
        t4=(y(i1)-y(i2))*u(i2)                                                  
        t5=(z(i1)-z(i2))*u(i2)                                                  
      sx=(x(i1)-x(i2))*isy(i1)*isy(i2)*t2(i2)
      sy=(y(i1)-y(i2))*isy(i1)*isy(i2)*t2(i2)
      sz=(z(i1)-z(i2))*isy(i1)*isy(i2)*t2(i2)

      if (ifdelt >= 1) then
      wx=(x(i1)-x(i2))*t2new(i2)
      wy=(y(i1)-y(i2))*t2new(i2)
      wz=(z(i1)-z(i2))*t2new(i2)
      frtwo(i1)=frtwo(i1)+t2new(i2)
      frtwo(i2)=frtwo(i2)+t2new(i2)
      else
        wx=(x(i1)-x(i2))*t2(i2)                                                 
        wy=(y(i1)-y(i2))*t2(i2)                                                 
        wz=(z(i1)-z(i2))*t2(i2)                                                 
      end if
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
        f3x=0.                                                                  
        f3y=0.                                                                  
        f3z=0.                                                                  
        do 11 i1=1,np                                                           
        dr2=(x(i1)-x(i))*(x(i1)-x(i))+(y(i1)-y(i))*(y(i1)-y(i))                 
     +  +(z(i1)-z(i))*(z(i1)-z(i))                                              
        if (i1 == i) goto 11                                                     
        f3x=f3x+(x(i)-x(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3y=f3y+(y(i)-y(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3z=f3z+(z(i)-z(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
11      continue                                                                
        f3x=potd*pote*f3x                                                       
        f3y=potd*pote*f3y                                                       
        f3z=potd*pote*f3z                                                       
      if (ifdelt >= 1) then
      st=potdnew*(1.-potenew*(frtwo(i)+1.)**potc)
      else
        st=potd*(1.-pote*(fr(i)+1)**potc)                                       
      end if
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

              subroutine updtv

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,xmu,asym,bsym,csym,gamsy
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/coul/cy(5000),cyu(5000)                                          
        common/veloc/velx(500),vely(500),velz(500)                              
        common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
         common/tt/prtime(0:8000),nruns,nt,rtime,perc
        common/densi1/denstl(500),dens1(500)
        dimension f1x(500),f1y(500),f1z(500),t2(500)                            
        dimension fcx(500),fcy(500),fcz(500),u(500)                             
        dimension fvx(500),fvy(500),fvz(500),v(500)                             
        dimension dvdpx(500),dvdpy(500),dvdpz(500),ppp(500)                              
        dimension fhx(500),fhy(500),fhz(500),upau(500)                          
        dimension fgx(500),fgy(500),fgz(500),gra(500)                           
        dimension fsx(500),fsy(500),fsz(500)

         wigth=2.*width
         fnorm=(3.141596*wigth)**1.5
         potgr=1./(fnorm*wigth)**2./0.165          
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
        vdens=xexp(-dr2/al)/sal                                                  
        pbr=p2*500.+1                                                           
        dvdp=factor*vdens*alog(pbr)/pbr                                         
        dvdp1=-anorm/widthp**2*xexp(-dr2/wq2)*xexp(-p2/wp2)                       
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
c                                                                               
        velx(i) = px(i)/e1 + dvdpx(i)                                           
        vely(i) = py(i)/e1 + dvdpy(i)                                           
        velz(i) = pz(i)/e1 + dvdpz(i)                                           
c                                                                               
        x(i) = x(i) + velx(i)*dt                                                
        y(i) = y(i) + vely(i)*dt                                                
        z(i) = z(i) + velz(i)*dt                                             
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
1       continue                                                                
c                                                                               
        do 2 i1=1,np                                                            
        do 21 i2=i1+1,np                                                        
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     +  +(z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     +  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        t2(i2)=xexp(-dr2/al)                                                     
        v(i2)=t2(i2)*(alog(p2*500.+1.))**2.                                     
        idr22=dr2*40.+1.5                                                       
        idr2=min0(idr22,100000)                                                  
        nccc=ih(idr2)                                                           
        u(i2)=iso(i1)*iso(i2)*cy(nccc)+cyu(nccc)
c       upau(i2)=anorm*cyu(nccc)*xexp(-p2**2/width)                              
        upau(i2)=anorm/widthq**2*xexp(-dr2/wq2)*xexp(-p2/wp2)  
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
ccc
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
     +  +(z(i1)-z(i))*(z(i1)-z(i))                                              
        if (i1 == i) goto 11                                                     
        f3x=f3x+(x(i)-x(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3y=f3y+(y(i)-y(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
        f3z=f3z+(z(i)-z(i1))*xexp(-dr2/al)*(fr(i1)+1)**potc                      
11      continue                                                                
        f3x=potd*pote*f3x                                                       
        f3y=potd*pote*f3y                                                       
        f3z=potd*pote*f3z                                                       
        st=potd*(1.-pote*(fr(i)+1)**potc)                                       
        fx=f1x(i)*st+fcx(i)+fhx(i)+fvx(i)*potf+
     &            fgx(i)*potg-f3x+fsx(i)*pots
        fy=f1y(i)*st+fcy(i)+fhy(i)+fvy(i)*potf+
     &            fgy(i)*potg-f3y+fsy(i)*pots
        fz=f1z(i)*st+fcz(i)+fhz(i)+fvz(i)*potf+
     &            fgz(i)*potg-f3z+fsz(i)*pots
        px(i)=px(i)+dt*fx                                                       
        py(i)=py(i)+dt*fy                                                       
        pz(i)=pz(i)+dt*fz                                            

258     continue                                                                
        return                                                                  
        end                                                                     
c///////////////////////////////////////////////////////////////////////////////

              subroutine updtim

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,          
     &rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt                            
      common/coul/cy(5000),cyu(5000)                                         
      common/rpi/rmst(5000),pmst(5000),pmsn(5000),pmsz(5000),avet(5000)
      common/densi1/denstl(500),dens1(500)
      common/gg/nzta,nzpr,nsym
      common/xx/alpha,beta,gamma,symm
      common/veloc/velx(500),vely(500),velz(500)
      common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01
      common/wavel/sl(500)
      common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
      common/frictn/gamma_fric_P, tcool, gamma_fric_R
      dimension f1x(500),f1y(500),f1z(500),t2(500),t2new(500) 
      dimension fcx(500),fcy(500),fcz(500),u(500)
      dimension fvx(500),fvy(500),fvz(500),vm(500) 
      dimension fsx(500),fsy(500),fsz(500),fgrx(500),fgry(500),fgrz(500)
     &      ,fsx1(500),fsy1(500),fsz1(500),fsx2(500),fsy2(500),fsz2(500)
      dimension fsgx(500),fsgy(500),fsgz(500)
      dimension dvdpx(500),dvdpy(500),dvdpz(500),ppp(500)
      dimension fhx(500),fhy(500),fhz(500),upau(500)
      dimension f1xc(500),f1yc(500),f1zc(500),frp(500)
cc
      pi=3.1415927
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
        dvdpx(i)=0.                                                            
        dvdpy(i)=0.                                                            
        dvdpz(i)=0.
1821    continue
 
        call densitn 
        do 1822 i1=1,np
        
        e1=sqrt(em(i1)*em(i1)+px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        pf=hbar*(3.*pi**2*dens1(i1))**0.3333
        if (pf <= 0.) then
        qrf=0.
        else 
        qrf=pr/pf
        endif                                                        
                                                                 
        do 1823 i2=i1+1,np
        al=2.*(sl(i1)+sl(i2))
        sal=(pi*al)**1.5*rho0
        xmm=sqrt(sl(i1)+sl(i2))/hbar
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &  (z(i1)-z(i2))*(z(i1)-z(i2))                                            
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                 
        vdens=xexp(-dr2/al)/sal

        if (momp == 0) then                 !without mom. dep.
        dvdp=0.
        elseif (momp == 1) then             !with mom. by Aichelin
          if (iso(i1) == iso(i2)) then
            cfm=cmnn
          else
            cfm=cmnp
          endif
          pbr= p2*cmeps+1.
          dvdp=4.*cfm*cmeps*vdens*alog(pbr)/pbr
        else                               !with mom. and mass splitting
          if (iso(i1) == iso(i2)) then
            cfm=cmnn
          else
            cfm=cmnp
          endif

          cfmom=-2.*cfm*xmm*vdens/pi3
          clmd=xld2*xmm**2
          xk1=xmm*px(i1)
          yk1=xmm*py(i1)
          zk1=xmm*pz(i1)
          xk2=xmm*px(i2)
          yk2=xmm*py(i2)
          zk2=xmm*pz(i2)
cc        call xmoms(clmd,xk1,yk1,zk1,xk2,yk2,zk2,sumx,sumy,sumz)     !Gaussian integration
          call xpart(1,np,clmd,xk1,yk1,zk1,xk2,yk2,zk2,sumx,sumy,sumz)
          dvdpx(i1)=dvdpx(i1)+cfmom*sumx
          dvdpy(i1)=dvdpy(i1)+cfmom*sumy
          dvdpz(i1)=dvdpz(i1)+cfmom*sumz
cc        call xmoms(clmd,xk2,yk2,zk2,xk1,yk1,zk1,sumx,sumy,sumz)     !Gaussian integration
          call xpart(1,np,clmd,xk2,yk2,zk2,xk1,yk1,zk1,sumx,sumy,sumz)
          dvdpx(i2)=dvdpx(i2)+cfmom*sumx
          dvdpy(i2)=dvdpy(i2)+cfmom*sumy
          dvdpz(i2)=dvdpz(i2)+cfmom*sumz
          goto 1823
        endif

        if (iso(i1) == iso(i2)) then
          dvdp1=-vpb*vdens/wp2*xexp(-p2/wp2)
          dvdp2=vpc*e01*vdens/(pf*(1.+xexp(e01*(1.-qrf**2))))**2
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
        velx(i1) = px(i1)/e1 + dvdpx(i1)                                          
        vely(i1) = py(i1)/e1 + dvdpy(i1)                                           
        velz(i1) = pz(i1)/e1 + dvdpz(i1)  
1822    continue

      do 1 i=1,np                                                                                                                                                         
        x(i) = x(i) + velx(i)*dt                                                
        y(i) = y(i) + vely(i)*dt                                                
        z(i) = z(i) + velz(i)*dt
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
        fhx(i)=0.                                                               
        fhy(i)=0.                                                               
        fhz(i)=0.                                                               
        fr(i)=0.
        frp(i)=0.
        f1xc(i)=0.
        f1yc(i)=0.
        f1zc(i)=0.
        fgrx(i)=0. 
        fgry(i)=0.
        fgrz(i)=0.
1     continue
c                                                                             
        do 2 i1=1,np                                                    

        pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        pf=hbar*(3.*pi**2*dens1(i1))**0.3333
        if (pf <= 0.) then
        qrf=0.
        else 
        qrf=pr/pf
        endif                                                        

        do 21 i2=i1+1,np
        al=2.*(sl(i1)+sl(i2))
        sal=(pi*al)**1.5*rho0    
        xmm=sqrt(sl(i1)+sl(i2))/hbar
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))             
     &  +(z(i1)-z(i2))*(z(i1)-z(i2))
        p2=(px(i1)-px(i2))*(px(i1)-px(i2))+(py(i1)-py(i2))*(py(i1)-             
     &  py(i2))+(pz(i1)-pz(i2))*(pz(i1)-pz(i2))                                            
        t2(i2)=xexp(-dr2/al)/sal
        idr22=dr2*40.+1.
        idr2=min0(idr22,100000)                                                  
        nccc=ih(idr2)                                                        
        u(i2)=iso(i1)*iso(i2)*cy(nccc)      !! cyu(nccc)
cc       mom. dep. interaction
        if (momp == 0) then
        vm(i2)=0.
        elseif (momp == 1) then
         if (iso(i1) == iso(i2)) then
         cfm=cmnn
         else
         cfm=cmnp
         endif
        vm(i2)=cfm*t2(i2)*(alog(p2*cmeps+1.))**2
        else
         if (iso(i1) == iso(i2)) then
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
cc        vm(i2)=cfm*t2(i2)/pi3*sumom(clmd,xk1,yk1,zk1,xk2,yk2,zk2)       !Gaussian integration
        vm(i2)=cfm*t2(i2)/pi3*xmsum(1,np,clmd,xk1,yk1,zk1,xk2,yk2,zk2)
        endif

        if (iso(i1) == iso(i2)) then                                         
        upau(i2)=vpa*(1./al+1./wq2)*xexp(-dr2/wq2)+vpb/al*xexp(-p2/wp2)+
     &           vpc/al/(1.+xexp(e01*(1.-qrf**2)))
        else
        upau(i2)=0.
        endif
21      continue

        do 22 i2=i1+1,np     
        al=2.*(sl(i1)+sl(i2))
        dr2=(x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &  (z(i1)-z(i2))*(z(i1)-z(i2))

        v3=2.*(x(i1)-x(i2))*vm(i2)/al                                                  
        v4=2.*(y(i1)-y(i2))*vm(i2)/al                                                  
        v5=2.*(z(i1)-z(i2))*vm(i2)/al

        t3=(x(i1)-x(i2))*u(i2)                                                  
        t4=(y(i1)-y(i2))*u(i2)                                                  
        t5=(z(i1)-z(i2))*u(i2)

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
cc 
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
        f1xc(i1)=f1xc(i1)+iso(i1)*iso(i2)*wx                                                      
        f1xc(i2)=f1xc(i2)-iso(i1)*iso(i2)*wx                                                      
        f1yc(i1)=f1yc(i1)+iso(i1)*iso(i2)*wy                                                      
        f1yc(i2)=f1yc(i2)-iso(i1)*iso(i2)*wy                                                      
        f1zc(i1)=f1zc(i1)+iso(i1)*iso(i2)*wz                                                      
        f1zc(i2)=f1zc(i2)-iso(i1)*iso(i2)*wz                                                   
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

        potc=gamma-1.
        emass=eta-1.
                                                                                  
        do 258 i=1,np                                                           
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
        al=2.*(sl(i1)+sl(i)) 
        sal=(pi*al)**1.5*rho0                                                        
        dr2=(x(i1)-x(i))*(x(i1)-x(i))+(y(i1)-y(i))*(y(i1)-y(i))                 
     &  +(z(i1)-z(i))*(z(i1)-z(i))                                              
        if (i1 == i) goto 11                                                     
        f3x=f3x+(x(i)-x(i1))*xexp(-dr2/al)/sal/al*fr(i1)**potc                      
        f3y=f3y+(y(i)-y(i1))*xexp(-dr2/al)/sal/al*fr(i1)**potc                      
        f3z=f3z+(z(i)-z(i1))*xexp(-dr2/al)/sal/al*fr(i1)**potc  
        ftx=ftx+(x(i)-x(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass          !effective mass                   
        fty=fty+(y(i)-y(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass                      
        ftz=ftz+(z(i)-z(i1))*xexp(-dr2/al)/sal/al*fr(i1)**emass
        f3xc=f3xc+(x(i)-x(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al       !Coulomb exchange
     &  /(frp(i1)+0.001)**0.6666   
        f3yc=f3yc+(y(i)-y(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al
     &  /(frp(i1)+0.001)**0.6666
        f3zc=f3zc+(z(i)-z(i1))*iso(i)*iso(i1)*xexp(-dr2/al)/sal/al
     &  /(frp(i1)+0.001)**0.6666     
        if (nsym == 2) then
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
         elseif (nsym == 3) then
        fsx13=fsx13+(x(i)-x(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamma-1.)
        fsy13=fsy13+(y(i)-y(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamma-1.)
        fsz13=fsz13+(z(i)-z(i1))*isy(i)*isy(i1)*XEXP(-dr2/al)/sal/al*
     &        (fr(i1)+0.0001)**(gamma-1.)
        fsx23=fsx23+(x(i)-x(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamma-2.)
        fsy23=fsy23+(y(i)-y(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamma-2.)
        fsz23=fsz23+(z(i)-z(i1))*isy(i)*isy(i1)*(XEXP(-dr2/al)/sal)**2
     &        /al*(fr(i1)+0.0001)**(gamma-2.)
        if (csym == 0.) goto 11
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

      f1xc(i)=Cex/3.*2.*f1xc(i)/(frp(i)+0.001)**0.6666*rho0
      f1yc(i)=Cex/3.*2.*f1yc(i)/(frp(i)+0.001)**0.6666*rho0
      f1zc(i)=Cex/3.*2.*f1zc(i)/(frp(i)+0.001)**0.6666*rho0
      f3xc=Cex/3.*2.*f3xc*rho0
      f3yc=Cex/3.*2.*f3yc*rho0
      f3zc=Cex/3.*2.*f3zc*rho0

        if (nsym == 2) then
        fs3x=symm*(fsx13+(xmu-1.)*fsx23)
        fs3y=symm*(fsy13+(xmu-1.)*fsy23)
        fs3z=symm*(fsz13+(xmu-1.)*fsz23)
        elseif (nsym == 3) then
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

        if (nsym == 1) then
        fsx(i)=2.*symm*fsx1(i)
        fsy(i)=2.*symm*fsy1(i)
        fsz(i)=2.*symm*fsz1(i)
        elseif (nsym == 2) then
        fsx(i)=symm*fsx1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsx2(i)*(fr(i)+0.0001)**(xmu-2.)
        fsy(i)=symm*fsy1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsy2(i)*(fr(i)+0.0001)**(xmu-2.)
        fsz(i)=symm*fsz1(i)*(fr(i)+0.0001)**(xmu-1.)+
     &         symm*(xmu-1.)*fsz2(i)*(fr(i)+0.0001)**(xmu-2.)
         elseif (nsym == 3) then
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

c        if (momp == 0) then
c        fsmx=0.
c        fsmy=0.
c        fsmz=0.
c        elseif (momp == 1) then
c        pfn=hbar*(3.*pi**2*dens1(i))**0.333333
c        pfp=hbar*(3.*pi**2*(denstl(i)-dens1(i)))**0.333333
c        xkpp=fmlog(pfp,pfp,cmeps)*cmnn
c        xknn=fmlog(pfn,pfn,cmeps)*cmnn
c        xkpn=fmlog(pfp,pfn,cmeps)*cmnp
c        cmmy=-0.5*(0.5*(xkpp+xknn)-xkpn)
c
c        fsmx=2.*cmmy*fsx1(i)
c        fsmy=2.*cmmy*fsy1(i)
c        fsmz=2.*cmmy*fsz1(i)
c        else
c        pfn=hbar*(3.*pi**2*dens1(i))**0.333333
c        pfp=hbar*(3.*pi**2*(denstl(i)-dens1(i)))**0.333333
c        xkpp=fmnp0(pfp,pfp,xld)*cmnn
c        xknn=fmnp0(pfn,pfn,xld)*cmnn
c        xkpn=fmnp0(pfp,pfn,xld)*cmnp
c        cmmy=-0.375*(0.5*(xkpp+xknn)-xkpn)
c
c        fsmx=2.*cmmy*fsx1(i)
c        fsmy=2.*cmmy*fsy1(i)
c        fsmz=2.*cmmy*fsz1(i)
c        endif

        fx=f1x(i)*st+fcx(i)+f3x+fsx(i)+0.5*gsur*fgrx(i)+fs3x
        fy=f1y(i)*st+fcy(i)+f3y+fsy(i)+0.5*gsur*fgry(i)+fs3y
        fz=f1z(i)*st+fcz(i)+f3z+fsz(i)+0.5*gsur*fgrz(i)+fs3z
        fx=fx+ftx+f1xc(i)+f3xc+fvx(i)+fhx(i)+0.5*gsym*fsgx(i)  !+fsmx
        fy=fy+fty+f1yc(i)+f3yc+fvy(i)+fhy(i)+0.5*gsym*fsgy(i)  !+fsmy
        fz=fz+ftz+f1zc(i)+f3zc+fvz(i)+fhz(i)+0.5*gsym*fsgz(i)  !+fsmz
        tnow = ntstep*dt
        if (tnow  <=  tcool.and.gamma_fric_P  >  0.0) then
          damp = gamma_fric_P * (1.0 - tnow/tcool)
          px(i) = px(i) + dt*fx - damp*px(i)*dt
          if (gamma_fric_R  >  0.0) then
            damp_R = gamma_fric_R * (1.0 - tnow/tcool)
            x(i) = x(i) - damp_R*x(i)*dt
            y(i) = y(i) - damp_R*y(i)*dt
            z(i) = z(i) - damp_R*z(i)*dt
          end if
          py(i) = py(i) + dt*fy - damp*py(i)*dt
          pz(i) = pz(i) + dt*fz - damp*pz(i)*dt
        else
          px(i) = px(i) + dt*fx
          py(i) = py(i) + dt*fy
          pz(i) = pz(i) + dt*fz
        end if
 258    continue
                                                           
        return                                                                  
        end    
c///////////////////////////////////////////////////////////////////////////////

              subroutine coulg1

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/coul/cy(5000),cyu(5000)                                          
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/gg/nzta,nzpr,nsym 
        common/atcy/sco(5000),syu(5000)                                                    
c                                                                               
        a1=al/4.                                                                
        gam=1./arang                                                            
        ct0r=0.005                                                              
        ceff=float(nzpr+nzta)/float(npr+nta)                                    
        c1=0.5*xexp(a1*gam*gam)                                                  
        c2=1./sqrt(a1*3.1415927)                                                
        c3=0.5*vy                                                               
        c4=sqrt(al)                                                             
        c5=0.5*0.00144*ceff*ceff 
        sco(1)=0.
        sco(5000)=0.                                               
c                                                                               
        t0=t0-ct0r/gam/gam/gam*vy                                               
c                                                                               
        cy(1)=0.                                                                
        cy(5000)=0.                                                             
        do 1 j=2,4999                                                           
         dr=(float(j)-1.)/100.                                                  
         dr2=dr*dr                                                              
         t1=dr/c4                                                               
         t2=xexp(-t1*t1)                                                         
         call erf(t1,t3)                                                        
         t4=t3/dr-1.128/c4*t2                                                   
         t5=c5/dr2*t4  
         sco(j)=c5/dr*t3                                                         
         t11=2.*a1*gam-dr                                                       
         t12=2.*a1*gam+dr                                                       
         if (t11 < 0.) then                                                     
         call erf(-t11/c4,t13)                                                  
         t13=-t13                                                               
         else                                                                   
         call erf(t11/c4,t13)                                                   
         end if                                                                 
         call erf(t12/c4,t14)                                                   
         t15=1.-t13                                                             
         t16=1.-t14                                                             
         t18=xexp(-gam*dr)                                                       
         t19=xexp(gam*dr)                                                        
         t20=t18*t15-t19*t16                                                    
         t21=t18*(gam*t15-c2*xexp(-t11*t11/a1))                                  
         t22=t19*(gam*t16-c2*xexp(-t12*t12/a1))                                  
         t23=c3*c1/dr2/gam*(1./dr*t20+t21+t22)                                  
         cy(j)=t5+t23                                                           
                                                                                
1       continue                                                                
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

              subroutine erf(x,er) 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        p=.3275911                                                              
        a1=.254829592                                                           
        a2=-.284496736                                                          
        a3=1.421413741                                                          
        a4=-1.453152027                                                         
        a5=1.061405429                                                          
        t=1./(1+p*x)                                                            
        er=1-(a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5)*xexp(-x**2)                  
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

              subroutine kaon(s,em1,em2,i1,i2)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/kaot/pktot,nkn                                                 
        data amk,aml,amn/ 0.495,1.115,0.938 /                                 
        pmax=(s-(amk+aml+amn)**2)*(s-(aml+amn-amk)**2)/4./s                   
        pk=0.0                                                                
        if (pmax < 0.0) return                                               
        pmax=sqrt(pmax)/amk/40.                                               
        if (em1 < 1..and.em2 < 1.) pk=0.072*pmax                            
        if (em1 < 1..and.em2 > 1.) pk=0.054*pmax                            
        if (em1 > 1..and.em2 < 1.) pk=0.054*pmax                            
        if (em1 > 1..and.em2 > 1.) pk=0.036*pmax                            
        pktot=pktot+pk                                                        
        nkn=nkn+1                                                             
        I1=I1
        I2=I2
        return                                                                
        end                                                                   

c///////////////////////////////////////////////////////////////////////////////

      subroutine minspa
      ! Minimum Spanning Tree (MST)
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,rrt,
     &  px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,ntag,x00
      dimension i(500,500),it(500),itt(500),ifrag(500)

      do k = 1, 500
        it(k) = 0
        itt(k) = 0
        ifrag(k) = 0
        ipy = 0
        do kk = 1, 500
          i(k, kk) = 0
        end do
      end do

      do m = 1, np
        mmm = 0
        do mm = 1, np
          if (m == mm) cycle
          dr = sqrt((x(m) - x(mm))**2 + (y(m) - y(mm))**2
     &        + (z(m) - z(mm))**2)
          if (dr < 3.) mmm = mmm + 1
          if (dr < 3.) i(m, mmm) = mm
          if (dr < 3.) it(m) = mmm
        end do
        if (it(m) == 0) ipy = ipy + 1
      end do

      do k = 1, np
        if (it(k) == 0) cycle     
        m1 = 0
        do while (i(k, m1 + 1) /= 0)
          m1 = m1 + 1
          ! 遍历当前邻居的子邻居
          do mm = 1, it(i(k, m1))
            ic = 0
            do m = 1, it(k)
              if (i(i(k, m1), mm) == i(k, m)) ic = 1
              if (i(i(k, m1), mm) == k) ic = 1
            end do
            ! 发现新核子，吸纳进当前碎片
            if (ic == 0) then
              it(k) = it(k) + 1
              i(k, it(k)) = i(i(k, m1), mm)
            end if         
            ! 清空已被融合的邻居的标记，防止死循环
            i(i(k, m1), mm) = 0
            it(i(k, m1)) = 0
          end do
        end do
        itt(k) = it(k) + 1
        ifrag(itt(k)) = ifrag(itt(k)) + 1
      end do

      do k = 1, np
        iev(k) = 1
      end do

      do k = 1, np
        ifram = 3
!       ifram= smallest cluster considered as fragment
        if (itt(k) < ifram) cycle
        iev(k) = 0
        do l = 1, it(k)
          iev(i(k, l)) = 0
        end do
      end do

      return
      end

c///////////////////////////////////////////////////////////////////////////////

              subroutine clallnew 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/deu3/dg(50,10,2),dtg(2),dtg1(2),
     &  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
ccc      nov. 27 1992 added
      do 3 l=1,50
      do 3 ll=1,10
      do 3 ideu=1,2
      dcm(l,ll,ideu)=0.
3      dg(l,ll,ideu)=0.
      do 4 ideu=1,2
      dtgt(ideu)=0.
4      dtg1t(ideu)=0.
ccc 
        return                                                                  
        end                                                                     
c///////////////////////////////////////////////////////////////////////////////

              subroutine clall 

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
         common/deu3/dg(50,10,2),dtg(2),dtg1(2),
     &  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
        common/bb/np,ntag,x00                                                   
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
      if (ranf > 1.0)ranf=ranf-1.0
      return
      end

c///////////////////////////////////////////////////////////////////////////////

              subroutine dtest(i1,i2,em1,em2)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp                                 
        common/deu2/px1,py1,pz1,px2,py2,pz2                                     
        common/deu3/dg(50,10,2),dtg(2),dtg1(2),                                 
     1  dcm(50,10,2),ideu,dtgt(2),dtg1t(2)                                      
        fiso=3./8.                                                              
        fa=8.*3./4.                                                          
        do 2 k=1,np                                                             
        im1=1                                                                   
        im2=1                                                                   
        im1o=1                                                                  
        im2o=1                                                                  
        if (k == i1)goto 2                                                       
        if (k == i2)goto 2                                                       
        xr1=((x(i1)-x(k))**2+(y(i1)-y(k))**2+(z(i1)-z(k))**2)                   
        xr2=((x(i2)-x(k))**2+(y(i2)-y(k))**2+(z(i2)-z(k))**2)                   
        if (xr1 > 70)im1=0                                                      
        if (xr2 > 70)im2=0                                                      
        if (xr1 > 70)im1o=0                                                     
        if (xr2 > 70)im2o=0                                                     
        pxr1=((px(i1)-px(k))**2+(py(i1)-py(k))**2+(pz(i1)-pz(k))**2)            
        pxr2=((px(i2)-px(k))**2+(py(i2)-py(k))**2+(pz(i2)-pz(k))**2)           
        if (pxr1 > 2)im1=0                                                      
        if (pxr2 > 2)im2=0                                                      
        pxs1=((px(i1)+px(k))**2+(py(i1)+py(k))**2+(pz(i1)+pz(k))**2)            
        pxs2=((px(i2)+px(k))**2+(py(i2)+py(k))**2+(pz(i2)+pz(k))**2)            
        pxr1o=((px1-px(k))**2+(py1-py(k))**2+(pz1-pz(k))**2)                    
        pxr2o=((px2-px(k))**2+(py2-py(k))**2+(pz2-pz(k))**2)                    
        if (pxr1o > 2)im1o=0                                                    
        if (pxr2o > 2)im2o=0                                                    
        pxs1o=((px1+px(k))**2+(py1+py(k))**2+(pz1+pz(k))**2)                    
        pxs2o=((px2+px(k))**2+(py2+py(k))**2+(pz2+pz(k))**2)                    
        prob1=fd1*fiso*xexp(-xr1*fdr*im1)*xexp(-pxr1*fdp*im1)*im1                 
        prob2=fd1*fiso*xexp(-xr2*fdr*im2)*xexp(-pxr2*fdp*im2)*im2                 
        prob1o=fd1*fiso*xexp(-xr1*fdr*im1)*xexp(-pxr1o*fdp*im1o)*im1o             
        prob2o=fd1*fiso*xexp(-xr2*fdr*im2)*xexp(-pxr2o*fdp*im2o)*im2o             
        iecm1=int((sqrt(pxs1+3.519)-1.876)/0.05)+1                              
        iecm2=int((sqrt(pxs2+3.519)-1.876)/0.05)+1                              
        iecm1o=int((sqrt(pxs1o+3.519)-1.876)/0.05)+1                            
        iecm2o=int((sqrt(pxs2o+3.519)-1.876)/0.05)+1                            
        acm1=(pz(i1)+pz(k))/sqrt(pxs1)                                          
        acm2=(pz(i2)+pz(k))/sqrt(pxs2)                                          
        acm1o=(pz1+pz(k))/sqrt(pxs1o)                                           
        acm2o=(pz2+pz(k))/sqrt(pxs2o)                                           
        if (acm1 > 0.9999)acm1=0.9999                                           
        if (acm2 > 0.9999)acm2=0.9999                                           
        if (acm1o > 0.9999)acm1o=0.9999                                         
        if (acm2o > 0.9999)acm2o=0.9999                                         
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
        if (iso(i1)+iso(k).ne.1)goto 212                                         
        if (iev(i1) == 0.or.iev(k) == 0)goto 212                                 
        if (iecm1 < 1.or.iecm1 > 50)goto211                                    
        dcm(iecm1,icm1,ideu)=dcm(iecm1,icm1,ideu)+prob1                         
        dtg(ideu)=dtg(ideu)+prob1                                               
211     if (iecm1o < 1.or.iecm1o > 50)goto212                                  
        if (id(i1) == 1.or.id(i1) == 51)goto212                                  
        dcm(iecm1o,icm1o,ideu)=dcm(iecm1o,icm1o,ideu)-prob1o                    
        dtg(ideu)=dtg(ideu)-prob1o                                              
212     if (iso(i2)+iso(k).ne.1)goto 214                                         
        if (iev(i2) == 0.or.iev(k) == 0)goto214                                  
        if (iecm2 < 1.or.iecm2 > 50)goto213                                    
        dcm(iecm2,icm2,ideu)=dcm(iecm2,icm2,ideu)+prob2                         
        dtg(ideu)=dtg(ideu)+prob2                                               
213     if (iecm2o < 1.or.iecm2o > 50)goto214                                  
        if (id(i2) == 1.or.id(i2) == 51)goto214                                  
        dcm(iecm2o,icm2o,ideu)=dcm(iecm2o,icm2o,ideu)-prob2o                    
        dtg(ideu)=dtg(ideu)-prob2o                                              
214     continue                                                                
        t1=pz(i1)+pz(k)                                                         
        e1=sqrt(pxs1+1.876**2)                                                  
        t2=e1                                                                   
        e2=sqrt(pxs2+1.876**2)                                                  
        e1o=sqrt(pxs1o+1.876**2)                                                
        e2o=sqrt(pxs2o+1.876**2)                                                
c       transformation to lab system                                            
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
c       print*,gam3,bzz,pz(i1),pz(k),t1,t2,pl1,e1                               
        e2=sqrt(p2**2+1.876**2)-1.876                                           
        p1o=sqrt(pl1o**2+(px1+px(k))**2+(py1+py(k))**2)                         
        p2o=sqrt(pl2o**2+(px2+px(k))**2+(py2+py(k))**2)                         
        e1o=sqrt(p1o**2+1.876**2)-1.876                                         
        e2o=sqrt(p2o**2+1.876**2)-1.876                                         
        wp1=pl1/p1                                                              
        if (wp1 > 0.9999)wp1=0.9999                                             
        wp2=pl2/p2                                                              
        if (wp2 > 0.9999)wp2=0.9999                                             
        wp3=pl1o/p1o                                                            
        if (wp3 > 0.9999)wp3=0.9999                                             
        wp4=pl2o/p2o                                                            
        if (wp4 > 0.9999)wp4=0.9999                                             
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
        if (iev(i1) == 0.or.iev(k) == 0)goto20                                   
        if (iso(i1)+iso(k).ne.1)goto20                                           

        if (n1 < 1.or.n1 > 50)goto10                                           
        dg(n1,ibit1,ideu)=dg(n1,ibit1,ideu)+prob1                               
        dtg1(ideu)=dtg1(ideu)+prob1                                             
10      if (n1o < 1.or.n1o > 50)goto20                                         

        if (id(i1) == 1.or.id(i1) == 51)goto 20                                  
        dg(n1o,ibit1o,ideu)=dg(n1o,ibit1o,ideu)-prob1o                          
20      if (iso(i2)+iso(k).ne.1)goto 2                                           
        if (iev(i2) == 0.or.iev(k) == 0)goto2                                    
        if (n2 < 1.or.n2 > 50)goto30                                           
        dg(n2,ibit2,ideu)=dg(n2,ibit2,ideu)+prob2                               
        dtg1(ideu)=dtg1(ideu)+prob2                                             
30      if (n2o < 1.or.n2o > 50)goto2                                          
        if (id(i2) == 1.or.id(i2) == 51)goto 2                                   
        dg(n2o,ibit2o,ideu)=dg(n2o,ibit2o,ideu)-prob2o                          
2       continue                                                                
        EM1=EM1
        EM2=EM2
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

               subroutine rpdens

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
c      calculation of radial density distribution 
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/gg/nzta,nzpr,nsym
        common/tt/prtime(0:8000),nruns,nt,rtime,perc
        common/rdens/rtd(200),ptd(200)
c      common/kkk/number,ilj
c      common radius
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
        if (id(i) == 1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        end if 
        if (id(i) == 51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        end if
   36   continue                                                                
        do 15 i=1,np
        if (id(i) == 1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        end if                                                                  
        if (id(i) == 51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
      end if
15      continue 
                                                               
      do 5 i=1,100
      rtd(i)=0.
      ptd(i)=0.
5      continue

      a1=2.*width
      a2=.1973289**2/a1
      pi=3.1415927
      coff1=1./(a1*pi)**1.5
      coff2=1./(pi*a2)**1.5/np**2

c      factor=width/2./radius/(a1*pi)**1.5*xexp(-radius**2./a1)
      do 10 radius=0.,rrp+5.,0.2
      number=int((radius+0.25)/0.2)
      if (radius == 0.) then
      do 20 nkl=1,np
      r2=x(nkl)**2.+y(nkl)**2.+z(nkl)**2.
      rtd(number)=rtd(number)+coff1*xexp(-r2/a1)
20      continue
      
      else

      texp=0.
      do 30 ki1=1,np
      r2=x(ki1)*x(ki1)+y(ki1)*y(ki1)+z(ki1)*z(ki1)
            if (sqrt(r2) == 0.) then
            rtd(number)=rtd(number)+coff1*xexp(-radius*radius/a1)
            else
      exp12=(xexp(sqrt(r2)*radius/width)-xexp(-sqrt(r2)*
     &radius/width))/sqrt(r2)
      texp=exp12*xexp(-r2/a1)*width/2./radius*xexp(-radius*radius/a1)
      rtd(number)=rtd(number)+coff1*texp
            end if
30      continue        
      end if
10      continue

      do 40 prd=0.,0.4,0.005
      number=int((prd+0.0055)/0.005)
      if (prd == 0.) then
      do 50 nkl=1,np
      pr2=px(nkl)**2.+py(nkl)**2.+pz(nkl)**2.
      ptd(number)=ptd(number)+coff2*xexp(-pr2/a2)
50      continue
      
      else
      
      texp=0.
      do 60 ki2=1,np
      pr2=px(ki2)*px(ki2)+py(ki2)*py(ki2)+pz(ki2)*pz(ki2)
            if (sqrt(pr2) == 0.) then
            ptd(number)=ptd(number)+coff2*xexp(-prd*prd/a2)
            else
      exp12=(xexp(sqrt(pr2)*prd/a2*2.)-xexp(-sqrt(pr2)*
     &prd/a2*2.))/sqrt(pr2)
      texp=exp12*xexp(-pr2/a2)*a2/4./prd*xexp(-prd*prd/a2)
      ptd(number)=ptd(number)+coff2*texp
            end if
60      continue        
      end if
40      continue

      return
      end
c///////////////////////////////////////////////////////////////////////////////

               subroutine rdeninteg

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

c      calculation of radial density distribution by integrated to cubic angle
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,ntag,x00                                                   
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/gg/nzta,nzpr,nsym
      common/tt/prtime(0:8000),nruns,nt,rtime,perc
c      common/rdens/td(200)
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
        if (id(i) == 1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        end if 
        if (id(i) == 51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        end if
   36   continue                                                                
        do 15 i=1,np
        if (id(i) == 1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        end if                                                                  
        if (id(i) == 51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
      end if
15      continue                                                                
      do 5 i=1,30
5      rinteg(i)=0.
      a1=2.*width
      pi=3.1415927
      do 10 radius=0.,rrp+rrt+4.5,.5
      number=int((radius+.5)/.5)
      if (radius == 0.) then
      factor=1./(a1*pi)**1.5/0.16
      do 20 nkl=1,np
      r2=x(nkl)**2.+y(nkl)**2.+z(nkl)**2.
      rinteg(number)=rinteg(number)+factor*xexp(-r2/a1)
20      continue
      go to 10
      end if
      eps=.0001
      factor=1./4./pi/(a1*pi)**1.5*xexp(-radius**2./a1)      
      call simps2(0.,pi,eps,kflag,n,sum)
      rinteg(number)=rinteg(number)+sum*factor
10      continue
      write (95,*)'time=',rtime
      limit=int((rrp+rrt+4.5)/.5)
      do 30 i=1,limit
      radius=float(i-1)*.5
      write (95,100)radius,rinteg(i)
100      format(1x,2f9.4)
30      continue
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            subroutine simps2(a,b,eps,kflag,n,sum)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

c      fixed xli values to calculate values of integral
      n2=1
      kflag=0
      c=abs(a)+abs(b)
      h2=.5*(b-a)
      call simps1(a,eps,kflag,n1,ss1)
      call simps1(b,eps,kflag,n3,ss2)
      n=n1+n3+2
      tb1=h2*(ss1+ss2)
4      xli=a-h2
      tb2=.5*tb1
      do 5 j=1,n2
      xli=xli+2.*h2
      call simps1(xli,eps,kflag,n1,s)
      n=n+n1+1
5      tb2=tb2+h2*s
      sum=.3333333333*(4.*tb2-tb1)
      n2=n2+n2
      if (n2 < 8)go to 6
      if (abs(sum-sum0) <= eps*(abs(sum)+1.)) return
6      sum0=sum
      tb1=tb2
      h2=.5*h2
      if (c+h2.ne.c)go to 4
      n=-n
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            subroutine simps1(xli,eps,kflag,n1,s)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      n1=1
      call fls(xli,y1,y2)
      h1=.5*(y2-y1)
      e=abs(y1)+abs(y2)
      t1=h1*(f(xli,y1)+f(xli,y2))
1      yli=y1-h1
      t2=.5*t1
      do 2 i=1,n1
      yli=yli+2.*h1
2      t2=t2+h1*f(xli,yli)
      s=.3333333333*(4.*t2-t1)
      n1=n1+n1
      if (n1 < 8)go to 3
      if (abs(s-s0) <= eps*(abs(s)+1.))return
3      s0=s
      t1=t2
      h1=.5*h1
      if (e+h1.ne.e)go to 1
      if (y1 == y2)return
      kflag=kflag+1
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            subroutine fls(xli,y1,y2)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      pi=3.1415927
      y1=0.
      y2=2.*pi
      XLI=XLI
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            function f(xli,yli)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/gg/nzta,nzpr,nsym
      common/kkk/number,ilj
      common radius
      ccc=0.
      pi=3.1415927
      a1=2.*width
      do 40 ilj=1,np
      ri2=x(ilj)**2.+y(ilj)**2.+z(ilj)**2.
      aaa=xexp(-ri2/2./width)
      bbb=x(ilj)*sin(xli)*cos(yli)+y(ilj)*sin(xli)*sin(yli)+
     &      z(ilj)*cos(xli)
      ddd=xexp(radius/width*bbb)
      ccc=ccc+ddd*aaa*sin(xli)
40      continue
      f=ccc/0.16                                                      
      return
      end

c///////////////////////////////////////////////////////////////////////////////

                subroutine axyzdensity

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/gg/nzta,nzpr,nsym
        common/tt/prtime(0:8000),nruns,nt,rtime,perc
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
        if (id(i) == 1) then                                                     
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        end if
        if (id(i) == 51) then                                                    
                       btz=btz+pz(i)/nta                                        
                       atz=atz+z(i)/nta                                         
                       btx=btx+px(i)/nta                                        
                       atx=atx+x(i)/nta                                         
                       bty=bty+py(i)/nta                                        
                       aty=aty+y(i)/nta                                         
        end if 
36      continue                                                                
        do 15 i=1,np 
        if (id(i) == 1) then                                                     
                      pz(i)=pz(i)-bpz
                      z(i)=z(i)-apz
                      px(i)=px(i)-bpx                                           
                      x(i)=x(i)-apx
                      py(i)=py(i)-bpy                                           
                      y(i)=y(i)-apy                                             
        end if                                                                  
        if (id(i) == 51) then                                                    
                       pz(i)=pz(i)-btz
                       z(i)=z(i)-atz
                       px(i)=px(i)-btx                                          
                       x(i)=x(i)-atx
                       py(i)=py(i)-bty                                          
                       y(i)=y(i)-aty                                            
      end if
15      continue                                                                
      a1=al/2.
      pi=3.1415927
      factor=1./(a1*pi)**1.5
      limit=int((rrt+rrp+4.5)/0.5)
      do 5 i=1,limit
      xtd(i)=0.
      ytd(i)=0.
      ztd(i)=0.
      axyz(i)=0.
5      continue
      do 10 xxx=0.,rrp+rrt+4.5,0.5
      i=int((xxx+0.5)/0.5)
      do 20 j=1,np
      xi2=(xxx-x(j))**2.+y(j)**2.+z(j)**2.
      a=xi2/a1
      aexp=xexp(-a)
      aexp=factor*aexp
      xtd(i)=xtd(i)+aexp
      yi2=(xxx-y(j))**2.+z(j)**2.+x(j)**2.
      b=yi2/a1
      bexp=xexp(-b)
      bexp=factor*bexp
      ytd(i)=ytd(i)+bexp
      zi2=(xxx-z(j))**2.+y(j)**2.+x(j)**2.
      c=zi2/a1
      cexp=xexp(-c)
      cexp=factor*cexp
      ztd(i)=ztd(i)+cexp
20      continue
10      continue
      write (21,*)'time=',rtime
      do 100 i=1,limit
      axyz(i)=(xtd(i)+ytd(i)+ztd(i))/3.
      dist=float(i-1)*0.5
      write (21,200)dist,axyz(i)
100      continue
200      format(1x,2f9.4)
      return
      end

c///////////////////////////////////////////////////////////////////////////////

              subroutine maxdensit

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                         
        common/densi1/denstl(500),dens1(500)                                      
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/tt/prtime(0:8000),nruns,nt,rtime,perc
        common/max/mden
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
        dd=xexp(-dr2/witdh)                                            
        dens1(i)=dens1(i)+dd                                          
        dens1(j)=dens1(j)+dd                                          
30      continue                                                      
        dens1(i)=(1.+dens1(i))/anorm                                  
      if (dens1(i) > mden)mden=dens1(i)
20      continue                                                      
      mden=mden/0.16
      write (110,55555)rtime,mden
55555      format(1x,2f9.4)
        return                                                        
        end                                                           

c///////////////////////////////////////////////////////////////////////////////

             subroutine coulgn 
! precomputes radial lookup tables for Coulomb and Yukawa potentials (with derivatives) 
! to enable fast interpolation during simulation, 
! avoiding repeated evaluation of error functions and exponentials.
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        double precision xdble                                                  
        common/coul/cy(5000),cyu(5000)                                          
        common/atcy/sco(5000),syu(5000)                                       
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
        common/bb/np,ntag,x00                                                   
        common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
        common/gg/nzta,nzpr,nsym
        a1=al/4.                                                               
        gam=1./arang                                                            
c        ceff=float(nzta+nzpr)/float(np)                                    
c umgeschrieben                                                                 
        c1 = .5                                                                 
        c2=1./sqrt(a1*3.1415927)                                                
        c3 = vy                                                                 
        c4=sqrt(al)                                                             
        c5 = .00144                                                 
        cy(1)=0.                                                                
        cy(5000)=0.                                                             
        cyu(1)=0.                                                               
        cyu(5000)=0.                                                            
      sco(1)=0.
      sco(5000)=0.
      syu(1)=0.
      syu(5000)=0.
        do 1 j=2,4999 
        dr=(float(j)-1.)/100.                                                  
        dr2=dr*dr                                                              
        t1=dr/c4                                                               
        t2=xexp(-t1*t1)                                                         
        xdble = dble( t1 )                                                     
        t3 = 1. - erfc ( xdble )                                              
        t4=t3/dr-1.128/c4*t2                                                   
        t5=c5/dr2*t4                                                           
      sco(j)=c5/dr*t3
      cy(j)=t5
        t11=2.*a1*gam-dr                                                       
        t12=2.*a1*gam+dr                                                       
        xdble = dble( t11/c4 )                                                 
        t15 = erfc( xdble )                                                   
        xdble = dble( t12/c4 )                                                 
        t16 = erfc( xdble )                                                   
        t18=xexp(-gam*dr)                                                       
        t19=xexp(gam*dr)                                                        
        t20= xexp(a1*gam*gam) * (t18*t15-t19*t16)                               
        t21=t18*( gam*xexp(a1*gam*gam)*t15 -                                    
     &      c2*xexp(-t11*t11/ al+a1*gam*gam ))                            
        t22=t19*( gam*xexp(a1*gam*gam)*t16 -                                    
     &      c2*xexp(-t12*t12/ al+a1*gam*gam ))                            
        t23=c3*c1/dr2 * (1./dr*t20+t21+t22)                                    
        syu(j)=c3*c1/dr*t20
        cyu(j)=t23                                                         
1       continue                                                                
        return                                                                  
        end                                                                     
c                                                                               

c///////////////////////////////////////////////////////////////////////////////

      function erfc (x)                                   

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      double precision s15adf                                                   
      double precision                 x                                        
      integer                          ifail                                    
      double precision                 t, xhi, xlo, y                           
      intrinsic                        abs, exp                                 
      data xlo/-6.5d0/                                                          
      data xhi/9.5d0/                                                           
      ifail = 0                                                                 
      if (x >= xhi) go to 20                                                    
      if (x <= xlo) go to 40                                                    
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
      s15adf = Dexp(-x*x)*y                                                      
      if (x < 0.0d0) s15adf = 2.0d0 - s15adf                                   
      erfc = sngl (s15adf)
      return                                                                    
   20 s15adf = 0.0d0                                                            
      erfc = sngl (s15adf)
      return                                                                    
   40 s15adf = 2.0d0                                                            
      erfc = sngl (s15adf)
      return                                                                    
      end                                                                            

c///////////////////////////////////////////////////////////////////////////////

              subroutine avene(avpot,aveki,Esh,idn)

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

       common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &     rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
       common/rpi/rmst(5000),pmst(5000),pmsn(5000),pmsz(5000),avet(5000)
       common/bb/np,ntag,x00 
       common/pot/pota,potb,potc,ih(100000),potd,pote,iflag5,potf               
     &        ,potg,pots,gsym,gtau,eta,Cex,xmu,asym,bsym,csym,gamsy
       common/aaa/npr,nta,gami,hc,nesc,tpo,ppo,rpt                                 
       common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp,pp                                 
       common/gg/nzta,nzpr,nsym
       common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt 
       common/pmom/pxta,pxpr,pzta,pzpr,x0,z0   
       common/coul/cy(5000),cyu(5000)                                          
       common/tt/prtime(0:8000),nruns,nt,rtime,perc
       common/atcy/sco(5000),syu(5000)
       common/xx/alpha,beta,gamma,symm
       common/vpaupa/v01,vpa,vpb,vpc,q01,p01,e01                                         
       common/densi1/denstl(500),dens1(500)
        common/grad/gsur,rho0,cdelt,cmeps,cmnn,cmnp,xsym,momp
       common/wavel/sl(500)
       common/encomps/tkine_comp,upaul_comp,upot_comp
       common/bind/progEb,aene,pfc,laesw
       common/shell/ish,r0,ds,Ri0
      common/raden/irad,rcount(100),r_tkine(100),r_usky(100),
     +r_ucou(100),r_usym(100),r_ussr(100),r_ugradi(100),
     +r_utau(100),r_couex(100),r_umomd(100),r_upau(100),
     +r_Esh(100),r_tene(100)
       real frme
       dimension t2(500),t2new(500),frp(500),rdist(500)

      ibinr(rr)=min0(100,max0(1,int(rr/0.2)+1))

      if (idn == 1) then
      imin=1
      imax=npr
      nnuc=npr
      end if
      if (idn == 2) then
      imin=npr+1
      imax=np
      nnuc=nta
      end if
      if (idn == 3) then
      imin=1
      imax=np
      nnuc=np
      end if
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
      Esh=0.
cc
      pi=3.1415927
      pi3=pi**3
      hb=0.197327
      Ekfac=12.4515   !MeV
      ccf=3.*pi**2./2.
      Pf1=hb*(1.5*pi**2*rho0)**0.3333
      xld= 1.5*Pf1
      xld2=xld**2

      r0p=r0*float(npr)**0.3333
      r0t=r0*float(nta)**0.3333
      r0c=r0*float(np)**0.3333
      Ri=Ri0
      Rf=r0c-r0t
      shp=1.+1.6228*ds/r0p
      sht=1.+1.6228*ds/r0t
      shc=1.+1.6228*ds/r0c
cc
    
      if (iflag5 == 1) then
      factor=cdelt                                                   
      else
      factor=0.
      end if

        widthp=p01                                                              
        widthq=q01                                                              
        wp2=2.*widthp**2
        wq2=2.*widthq**2                      
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

      OD3=1./3.
      TD3=2./3.
      CB=0.0058
      CS=0.26        
      ip=int(npr-nzpr+0.0001)
      jp=nzpr
      Eshp=CB*((FS(ip)+GS(jp))/(float(npr)/2.)**TD3-CS*float(npr)**OD3)
      Eshp=Eshp/npr*shp  
      it=int(nta-nzta+0.0001)
      jt=nzta
      Esht=CB*((FS(it)+GS(jt))/(float(nta)/2.)**TD3-CS*float(nta)**OD3)
      Esht=Esht/nta*sht 
      ic=int(np-nzpr-nzta+0.0001)
      jc=nzpr+nzta
      Eshc=CB*((FS(ic)+GS(jc))/(float(np)/2.)**TD3-CS*float(np)**OD3)
      Eshc=Eshc/np*shc           
                                                                         
        do 36 i=imin,imax                                                            
        fr(i)=0.
        frp(i)=0.
        em(i)=0.938                                                    
        if (id(i) == 1) then                                                    
                      bpz=bpz+pz(i)/npr                                         
                      apz=apz+z(i)/npr                                          
                      bpx=bpx+px(i)/npr                                         
                      apx=apx+x(i)/npr                                          
                      bpy=bpy+py(i)/npr                                         
                      apy=apy+y(i)/npr                                          
        end if
        if (id(i) == 51) then                                                    
                      btz=btz+pz(i)/nta                                        
                      atz=atz+z(i)/nta                                         
                      btx=btx+px(i)/nta                                        
                      atx=atx+x(i)/nta                                         
                      aty=bty+py(i)/nta                                        
                      aty=aty+y(i)/nta                                         
        end if 
   36   continue  

        Rpt=(apx-atx)**2+(apy-aty)**2+(apz-atz)**2
        Rpt=sqrt(Rpt)
      if (idn.ne.3) Rpt=Ri
      r01=r0p*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))
      r02=r0t*SWITCH(Ri,Rf,Rpt)+r0c*(1.-SWITCH(Ri,Rf,Rpt))
      Esh1=Eshp*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))
      Esh2=Esht*SWITCH(Ri,Rf,Rpt)+Eshc*(1.-SWITCH(Ri,Rf,Rpt))

        do 136 i=imin,imax 
        if (id(i) == 1) then                                                   
                     pz(i)=pz(i)-bpz
                     z(i)=z(i)-apz-z0
                     px(i)=px(i)-bpx                                           
                     x(i)=x(i)-apx
                     py(i)=py(i)-bpy                                           
                     y(i)=y(i)-apy                             
        end if                                                                
        if (id(i) == 51) then                                                    
                     pz(i)=pz(i)-btz
                     z(i)=z(i)-atz+z0
                     px(i)=px(i)-btx                                          
                     x(i)=x(i)-atx
                     py(i)=py(i)-bty                                          
                     y(i)=y(i)-aty                                        
        end if
136     continue
      if (irad == 1.and.idn == 3) then
        xcms=0.
        ycms=0.
        zcms=0.
        do 137 i=imin,imax
          xcms=xcms+x(i)
          ycms=ycms+y(i)
          zcms=zcms+z(i)
137     continue
        xcms=xcms/float(np)
        ycms=ycms/float(np)
        zcms=zcms/float(np)
        do 138 i=imin,imax
          rdist(i)=sqrt((x(i)-xcms)**2+(y(i)-ycms)**2+
     &                  (z(i)-zcms)**2)
138     continue
        do 139 i=1,100
          rcount(i)=0.
          r_tkine(i)=0.
          r_usky(i)=0.
          r_ucou(i)=0.
          r_usym(i)=0.
          r_ussr(i)=0.
          r_ugradi(i)=0.
          r_utau(i)=0.
          r_couex(i)=0.
          r_umomd(i)=0.
          r_upau(i)=0.
          r_Esh(i)=0.
          r_tene(i)=0.
139     continue
      end if
cc
      call densitn                                                           
      do 22 i1=imin,imax
        pr=sqrt(px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        e1=sqrt(em(i1)*em(i1)+px(i1)*px(i1)+py(i1)*py(i1)+pz(i1)*pz(i1))
        pf=hb*(3.*pi**2*dens1(i1))**0.3333
        if (pf <= 0.) then
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
      cucon=iso(i1)*iso(i2)*sco(nccc)
      ucou=ucou+cucon
      ugcon=gsur*(6./al-4.*dr2/al**2.)*t2(i2)
      ugradi=ugradi+ugcon
c      uyuk=uyuk+syu(nccc)
      if (iso(i1) == iso(i2)) then
      upcon=vpa*xexp(-dr2/wq2)*t2(i2)+vpb*xexp(-p2/wp2)*t2(i2)+
     &      vpc*t2(i2)/(1.+xexp(e01*(1.-qrf**2)))
      upaul=upaul+upcon
      endif
        if (momp == 0) then
      umomd=0.
        elseif (momp == 1) then
         if (iso(i1) == iso(i2)) then
         cfm=cmnn
         else
         cfm=cmnp
         endif
      umcon=cfm*alog(pbr)*alog(pbr)*t2(i2)/2.
      umomd=umomd+umcon
        else
         if (iso(i1) == iso(i2)) then
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
cc        xss=sumom(clmd,xk1,yk1,zk1,xk2,yk2,zk2)     !Gaussian integration
        xss=xmsum(imin,imax,clmd,xk1,yk1,zk1,xk2,yk2,zk2)
      umcon=cfm*t2(i2)/pi3*xss/2.
      umomd=umomd+umcon
        endif
        if (irad == 1.and.idn == 3) then
          ib1=ibinr(rdist(i1))
          ib2=ibinr(rdist(i2))
          r_ucou(ib1)=r_ucou(ib1)+cucon/2.
          r_ucou(ib2)=r_ucou(ib2)+cucon/2.
          r_ugradi(ib1)=r_ugradi(ib1)+ugcon/2.
          r_ugradi(ib2)=r_ugradi(ib2)+ugcon/2.
          if (iso(i1) == iso(i2)) then
            r_upau(ib1)=r_upau(ib1)+upcon/2.
            r_upau(ib2)=r_upau(ib2)+upcon/2.
            if (momp.ne.0) then
              r_umomd(ib1)=r_umomd(ib1)+umcon/2.
              r_umomd(ib2)=r_umomd(ib2)+umcon/2.
            end if
          end if
        end if

23      continue                                       
c
        do 21 i2=i1+1,imax                                                       
        fr(i1)=fr(i1)+t2(i2)                                                   
        fr(i2)=fr(i2)+t2(i2)  
        
        frp(i1)=frp(i1)+iso(i1)*iso(i2)*t2(i2)*rho0
        frp(i2)=frp(i2)+iso(i1)*iso(i2)*t2(i2)*rho0
                                                
21      continue                     
22      continue 
                                         
        do 258 i=imin,imax
      skine=(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))/2./em(i)
      tkine=tkine+skine
      Ekin=Ekin+Ekfac*(ccf*fr(i)*rho0)**0.66667

      if (i <= npr) then
      r=(x(i)-apx)**2+(y(i)-apy)**2+(z(i)-apz)**2
      r=sqrt(r)
      eshcon=Esh1/(1.+xexp((r-r01)/ds))
      Esh=Esh+eshcon
      else
      r=(x(i)-atx)**2+(y(i)-aty)**2+(z(i)-atz)**2
      r=sqrt(r)
      eshcon=Esh2/(1.+xexp((r-r02)/ds))
      Esh=Esh+eshcon
      endif
        if (irad == 1.and.idn == 3) then
          ib=ibinr(rdist(i))
          rcount(ib)=rcount(ib)+1.
          r_tkine(ib)=r_tkine(ib)+skine
          r_Esh(ib)=r_Esh(ib)+eshcon
        end if
258     continue

      do 220 i1 = imin, imax
        do 230 i2 = i1+1, imax
          al = 2.*(sl(i1)+sl(i2))  
          sal = (pi*al)**1.5*rho0
          if (i1 == i2) goto 230
        
          dr2 = (x(i1)-x(i2))*(x(i1)-x(i2))+(y(i1)-y(i2))*(y(i1)-y(i2))+
     &      (z(i1)-z(i2))*(z(i1)-z(i2))                         
          ex12 = xexp(-dr2/al)
        
          usscon=gsym*(6./al-4.*dr2/al**2.)*isy(i1)*isy(i2)*ex12/2./sal
          ussr = ussr + usscon
        
          if (nsym == 1) then
            usycon = symm*isy(i1)*isy(i2)*ex12/2./sal
            usym = usym + usycon
          elseif (nsym == 2) then
            frme = 0.5*(fr(i1) + fr(i2))
            usycon = symm*isy(i1)*isy(i2)*ex12/2./sal*frme**(xmu-1.)
            usym = usym + usycon
          elseif (nsym == 3) then
            frme = 0.5*(fr(i1) + fr(i2))
            usycon = asym*isy(i1)*isy(i2)*ex12/2./sal  
     &             + bsym*isy(i1)*isy(i2)*ex12/2./sal*frme**(gamsy-1.)
     &             + csym*isy(i1)*isy(i2)*ex12/2./sal*frme**(eta-1.)
            usym = usym + usycon
          else
            usym = 0.
            usycon = 0.
          end if
          if (irad == 1.and.idn == 3) then
            ib1 = ibinr(rdist(i1))
            ib2 = ibinr(rdist(i2))
            r_ussr(ib1) = r_ussr(ib1) + usscon
            r_ussr(ib2) = r_ussr(ib2) + usscon
            r_usym(ib1) = r_usym(ib1) + usycon
            r_usym(ib2) = r_usym(ib2) + usycon
          end if
230     continue
220   continue

      do 1 i = imin, imax
        utwcon=1./2.*alpha*fr(i)
        uthcon=1./(1.+gamma)*beta*fr(i)**gamma
        utacon=gtau*fr(i)**eta
        cexcon=Cex*frp(i)**0.3333
        utwo=utwo+utwcon
        uthr=uthr+uthcon
        utau=utau+utacon
        couex=couex+cexcon
        if (irad == 1.and.idn == 3) then
          ib=ibinr(rdist(i))
          r_usky(ib)=r_usky(ib)+utwcon+uthcon
          r_utau(ib)=r_utau(ib)+utacon
          r_couex(ib)=r_couex(ib)+cexcon
        end if
1      continue
      usky=utwo+uthr

      tene=tkine+usky+ucou+usym+ussr+ugradi+utau+couex+umomd+upaul+Esh
      aene= tene*1000./float(nnuc)
      avpot=(usky+ucou+usym+ugradi+utau+couex+umomd+upaul)*1000.
      aveki=Ekin

	if (mod(ntstep,100) == 0.and.laesw == 1) then
	open(3,file='007E.out',status='unknown')
	if (rtime == 0.) write (3,*)'     Running NO: ', nru
      write (3,*)'time        =',rtime,                  ' fm/c'
	write (3,*)'K total           =',tkine*1000./float(nnuc), ' MeV/u'
	write (3,*)'V Coulomb         =',ucou*1000./float(nnuc),  ' MeV/u'
      write (3,*)'V Coulomb exchange=',couex*1000./float(nnuc), ' MeV/u'
	write (3,*)'V Yukawa          =',uyuk*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V isospin         =',usym*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V surface symmetry=',ussr*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V pauli           =',upaul*1000./float(nnuc), ' MeV/u'
	write (3,*)'V mom. dep.       =',umomd*1000./float(nnuc), ' MeV/u'
	write (3,*)'V skyrme          =',usky*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V skyrme (2-body) =',utwo*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V skyrme (3-body) =',uthr*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V gradient        =',ugradi*1000./float(nnuc),' MeV/u'
	write (3,*)'V effective mass  =',utau*1000./float(nnuc),  ' MeV/u'
	write (3,*)'V shell correction=',Esh*1000./float(nnuc),   ' MeV/u'
	write (3,*)'E average         =',aene,                    ' MeV/u'
	write (3,*)
	end if
      if (llt > 0) avet(llt)=aene
      tkine_comp=tkine*1000./float(nnuc)
      upaul_comp=upaul*1000./float(nnuc)
      upot_comp=(tene-tkine)*1000./float(nnuc)

      if (irad == 1.and.idn == 3) then
        open(33,file='rad_ener.dat',status='unknown')
        write (33,909) 'r(fm)','Ncount','tkine','usky','ucou',
     &  'usym','utau','umomd','upaul','ugradi','ussr','couex',
     &  'Esh','total'
        do 200 i=1,100
          rmid=(float(i)-0.5)*0.2
          if (rcount(i) > 0.5) then
            r_tkine(i)=r_tkine(i)*1000.
            r_usky(i)=r_usky(i)*1000.
            r_ucou(i)=r_ucou(i)*1000.
            r_usym(i)=r_usym(i)*1000.
            r_ussr(i)=r_ussr(i)*1000.
            r_ugradi(i)=r_ugradi(i)*1000.
            r_utau(i)=r_utau(i)*1000.
            r_couex(i)=r_couex(i)*1000.
            r_umomd(i)=r_umomd(i)*1000.
            r_upau(i)=r_upau(i)*1000.
            r_Esh(i)=r_Esh(i)*1000.
            rtotal=r_tkine(i)+r_usky(i)+r_ucou(i)+r_usym(i)+
     &             r_ussr(i)+r_ugradi(i)+r_utau(i)+r_couex(i)+
     &             r_umomd(i)+r_upau(i)+r_Esh(i)
            write (33,910) rmid,int(rcount(i)),
     &        r_tkine(i)/rcount(i), r_usky(i)/rcount(i),
     &        r_ucou(i)/rcount(i), r_usym(i)/rcount(i),
     &        r_utau(i)/rcount(i), r_umomd(i)/rcount(i),
     &        r_upau(i)/rcount(i), r_ugradi(i)/rcount(i),
     &        r_ussr(i)/rcount(i), r_couex(i)/rcount(i),
     &        r_Esh(i)/rcount(i), rtotal/rcount(i)
          end if
200     continue
        close(33)
909     format(14(1x,a12))
910     format(f12.4,i12,13(1x,f12.6))
      end if

      return
      end                                                                     

c///////////////////////////////////////////////////////////////////////////////

            subroutine search

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
c       common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
c       common/skyrme/al,sal,arang,vy,width
      common/tt/prtime(0:8000),nruns,nt,rtime,perc
c      common/gg/nzta,nzpr
      dimension nc(300,300)
      real logx
      do 135 i=1,np
      do 136 j=1,np
      nc(i,j)=0
136      continue
135      continue      
        hb=0.197327                                                     
      hb2=hb/2.
      ncount=0
        do 137 i=1,np                                                           
      k=i+1
      do 138 j=k,np
               drx=x(i)-x(j)
      dpx=px(i)-px(j)
            logx=abs(drx*dpx)      
      if (logx < hb2) then
      nc(i,j)=1
      ncount=ncount+1
        end if 
138      continue
137     continue                                                                
      write (25,18476)rtime,ncount
      do 139 i=1,np
      write (25,*)(nc(i,j),j=1,np)
139      continue
18476      format(1x,f9.4,i5)
        return                                                                  
        end                                                                     

c///////////////////////////////////////////////////////////////////////////////

            subroutine coulomb_modify

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr,nsym
      common/cmod/coeft,coefp,epp,ecm,ptcm
      coeft=float(npr)/float(npr+nta)
      coefp=float(nta)/float(npr+nta)
      ecm=float(npr*nta)/float(nta+npr)*epp
      angm=(938./9)**.5*.1/6.5821*x00*sqrt(2.0*nta*npr/(npr+nta)*ecm)
      call coulomb_second
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            subroutine coulomb_first

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common/bb/np,ntag,x00                                                   
      common/pmom/pxta,pxpr,pzta,pzpr,thei,dptc
      common/cmod/coeft,coefp,epp,ecm,ptcm
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr,nsym
      common/rp1t/ xx(6),yy(6),zz(6)
      angm=(938./9)**.5*.1/6.5821*x00*sqrt(2.0*nta*npr/(npr+nta)*ecm)
      xx(1)=.72*(nzta*nzpr)/ecm+sqrt(.5187*(nzta*nzta*nzpr
     +      *nzpr)/(ecm*ecm)+20.89*angm*angm*(nta+npr)/(nta*npr*ecm))      
      xx(2)=atan(.6346*10.*angm*sqrt((nta+npr)*ecm/(nta*npr)
     +      /(nzta*nzpr)))
      xx(3)=0.0
      xx(4)=634.9*(nta+npr)/(nta*npr*xx(1)*xx(1))*angm
      angc=57.3*xx(2)
      csa=180.-2.*angc
      p2=13900.*(nta+npr)*(nzta*nzpr)/(nta*npr*xx(1)*xx(1))
50      d1=p2
      xx(5)=xx(1)*xx(4)*xx(4)+d1
      xx(6)=-2.*xx(3)*xx(4)/xx(1)
      do 3 j=1,2
      yy(j)=xx(j)-dt*xx(j+2)
      yy(j+2)=xx(j+2)-dt*xx(j+4)
3      continue
      p2=13900.*(nta+npr)*nzta*nzpr/(nta*npr*yy(1)*yy(1))
51      d1=p2
      yy(5)=yy(1)*yy(4)*yy(4)+d1
      yy(6)=-2.*yy(3)*yy(4)/yy(1)
      do 4 j=1,2
      zz(j)=xx(j)-.5*dt*(xx(j+2)+yy(j+2))
      zz(j+2)=xx(j+2)-.5*dt*(xx(j+4)+yy(j+4))
4      continue
      p2=13900.*(nta+npr)*nzta*nzpr/(nta*npr*zz(1)*zz(1))
      d1=p2
      zz(5)=zz(1)*zz(4)*zz(4)+d1
      zz(6)=-2.*zz(3)*zz(4)/zz(1)
      if (abs(zz(1)-yy(1)) < 0.001)go to 52
      do 5 j=1,2
      yy(j)=zz(j)
      yy(j+2)=zz(j+2)
5      continue
      go to 51
52      if (zz(1) >= ptcm)go to 53
      do 6 j=1,2
      xx(j)=zz(j)
      xx(j+2)=zz(j+2)
6      continue
      go to 50
53      do 7 j=1,2
      xx(j)=zz(j)
      xx(j+2)=zz(j+2)
7      continue
      px=312.67*(xx(3)*sin(xx(2))+xx(1)*xx(4)*cos(xx(2)))/dt
      pz=312.67*(xx(3)*cos(xx(2))+xx(1)*xx(4)*sin(xx(2)))/dt
      pxta=px*coeft*.001
      pxpr=px*coefp*.001
      pzta=pz*coeft*.001
      pzpr=pz*coefp*.001
      dptc=xx(1)      
      thei=xx(2)
      x0i=abs(dptc*sin(thei))
      x0=x0i/2.
      z0i=abs(dptc*cos(thei))
      z0=z0i/2.
      print*,dptc,thei,pxta,pxpr,pzta,pzpr,x0,z0
100      format(1x,6f9.4)
      return
      end

c///////////////////////////////////////////////////////////////////////////////

            subroutine coulomb_second

c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      common/pmom/pxta,pxpr,pzta,pzpr,x0,z0
      common/bb/np,ntag,x00                                                   
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/gg/nzta,nzpr,nsym
      common/cmod/coeft,coefp,epp,ecm,ptcm
      pi=3.1415927
      redm=float(nta*npr)/float(nta+npr)
      redm=redm*938.
      ech12=1.44*float(nzpr*nzta)
      orbl=x00*sqrt(2.0*redm*ecm)
      epsi=1.0+2.0*ecm*orbl*orbl/redm/ech12**2.0
      epsi=sqrt(epsi)
      dista=ptcm
      epsi1=1.0/epsi
      epsi2=1.0+orbl**2./redm/ech12/dista
      epsi3=epsi2*epsi1
      thei=pi-acos(epsi1)+acos(epsi3)
      rdot=2./redm*(ecm-ech12/dista-orbl**2./(2.*redm*dista**2))
      if (rdot <= 0.) then
      rdot=0.
      else
      rdot=sqrt(rdot)
      endif
      if (x00 == 0) then
      thdot=0.
      else
      thdot=sqrt(epsi**2.-epsi2**2.)
      thdot=rdot*orbl**2./redm/ech12/dista**2./thdot
      end if
      zdot=rdot*cos(thei)-dista*sin(thei)*thdot
      xdot=rdot*sin(thei)+dista*cos(thei)*thdot
      zdot=abs(zdot)*0.938
      xdot=abs(xdot)*0.938
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
            
            subroutine total_energy
            
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
        common/bb/np,ntag,x00                                                   
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
        common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/deu1/iflag6,bzz,gam3,fd1,fdr,fdp
      common/tt/prtime(0:8000),nruns,nt,rtime,perc
      tote=0.
      do 1 i=1,np
        tote =tote+sqrt(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i)+gam3*em(i)
     +      *gam3*em(i))
1      continue
      if (iprin == 1) then
      print*,'total energy=',tote
      end if
      return
      end


CLW*********************************************************************** 
      
      subroutine crossnew
                                                       
CLW***********************************************************************     
      common/boite/ann(800),anp(800)   
                    
      xpi=3.1415965e0*10.e0                                              
      do 105 i=1,800                                                    
      e=float(i)                           
      if (i < 25) then                                                   
      ann(i)=80.6e0                                                     
      anp(i)=381.0e0                                                    
      goto 100                                                         
      endif                                                             
      if (i < 40) then                                                   
      ann(i)=5.3107e0+3088.5e0/e-1174.8e0/e/e                           
      anp(i)=6.9466e0+9069.2e0/e-5057.4e0/e/e                           
      goto 100                                                          
      endif                                                             
      if (i < 310) then                                                  
      ann(i)=22.429e0-11.148e0/e+93074.e0/e/e                           
      anp(i)=27.147e0+1802.0e0/e+239380.e0/e/e                          
      goto 100                                                          
      endif                                                             
      if (i <= 800) then                                                  
      ann(i)=3.5475e0+.05331e0*e +887.37e0/e                            
      anp(i)=34.5e0                                                     
      goto 100                                                          
      endif                                                             
  100 ann(i)=ann(i)/xpi                                                 
      anp(i)=anp(i)/xpi                                                 
  105 continue                                                          

      return                                                            
      end                                                               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function xexp(x)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (x <= -87.) then
      xexp=exp(-87.)
      return
      end if
      
      if (x >= 87.) then
      xexp=exp(87.)
      return
      end if

      if (x > -87..and.x < 87.) then
      xexp=exp(x)
      return
      end if

      end

clwcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE RADRI
clwcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PARAMETER (MG=500,PI=3.14159)
      EXTERNAL rpdens
      common/ee/mseed,ird                                                
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo                                 
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/gg/nzta,nzpr
      common/riclw/rrrr(500)
      COMMON /WORK2/ RI(MG),PD(MG)
      COMMON /WORK1/ Mnumn,Mnump,Rmaxn,Rmaxp,RMS0N,RMS0P


      if (NPR <= 0) THEN
      stop 'Mass number of projectile nucleus should be positive'
      END IF
      NZ=NZPR
      NN=NPR-NZPR

      Imax=1
      DO I=1,MnumP
      read (14,*) RI(I),PD(I)
      PD(I)=4.*PI*RI(I)*RI(I)*PD(I)
       if (I > 1) THEN
         if (PD(I) >= PD(I-1))Imax=I
       END IF
      END DO
      REWIND(14)
c
9712      do 20 I=1,NZ
2        rand1=ran(mseed)
        rand2=ran(mseed)
      Rtest=RI(MnumP-1)*RAND1
       DO II=1,MnumP-1
        if (RI(II) <= Rtest.AND.RI(II+1) > Rtest) Itest=II
       END DO
      SLOP=(PD(Itest+1)-PD(Itest))/(RI(Itest+1)-RI(Itest))
      PDtest=PD(Itest)+SLOP*(Rtest-RI(Itest))
      RANtest=PDtest/PD(Imax)
      if (RAND2 >= RANtest) GOTO 2
      RRRR(I)=Rtest
      if (Rtest > Rmaxp) GOTO 2
  20       CONTINUE

      rmst=0.
      do i=1,NZ
      rmst=rmst+rrrr(i)*rrrr(i)
      end do
      rmst=sqrt(rmst/float(NZ))
      if (ABS(RMST-RMS0P) > 0.3) GOTO 9712


      Imax=1
      DO I=1,MnumN
      read (15,*) RI(I),PD(I)
      PD(I)=4.*PI*RI(I)*RI(I)*PD(I)
       if (I > 1) THEN
         if (PD(I) >= PD(I-1))Imax=I
       END IF
      END DO
      REWIND(15)
c
9713      do 220 I=NZ+1,NPR
22        rand1=ran(mseed)
        rand2=ran(mseed)
      Rtest=RI(MnumN-1)*RAND1
       DO II=1,MnumN-1
        if (RI(II) <= Rtest.AND.RI(II+1) > Rtest) Itest=II
       END DO
      SLOP=(PD(Itest+1)-PD(Itest))/(RI(Itest+1)-RI(Itest))
      PDtest=PD(Itest)+SLOP*(Rtest-RI(Itest))
      RANtest=PDtest/PD(Imax)
      if (RAND2 >= RANtest) GOTO 22
      RRRR(I)=Rtest
      if (Rtest > RmaxN) GOTO 22
  220       CONTINUE

      rmst=0.
      do i=NZ+1,NPR
      rmst=rmst+rrrr(i)*rrrr(i)
      end do
      rmst=sqrt(rmst/float(NN))
      if (ABS(RMST-RMS0N) > 0.3) GOTO 9713

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine width0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo,rpt
      common/bb/np,ntag,x00  
      common/CW/c1,c2
      common/wavel/sl(500) 
        
      do i = 1, npr
        sl(i)=(c1*float(npr)**0.3333 + c2)**2.
      end do

      do i = npr + 1, np
        sl(i) = (c1*float(nta)**0.3333 + c2)**2. 
      end do

      return
      end
cc    fengzhq
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
      subroutine constr
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,          
     &  rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)    
        common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/bb/np,ntag,x00   
        common/tt/prtime(0:8000),nruns,nt,rtime,perc
      common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt
      common/ee/mseed,nseed
      common/poli/iflagcro,iflagpau,ipoli
      
      pi=3.1415927
      hb=0.197327
      widthr=width
      widthp=hb**2./widthr/4. 
      widthsr=3.*sqrt(widthr)
      widthsm=3.*sqrt(widthp)

      do 1 i=1,np
      f1=Fn(i)
      if (f1 > 1.) then
        do 2 j=i+1,np
        if (iso(i) == iso(j)) then
      xxr=x(i)-x(j)
      yyr=y(i)-y(j)
      zzr=z(i)-z(j)
      ppxr=px(i)-px(j)
      ppyr=py(i)-py(j)
      ppzr=pz(i)-pz(j)
      rrs=sqrt(xxr**2.+yyr**2.+zzr**2.)
      pps=sqrt(ppxr**2.+ppyr**2.+ppzr**2.)
        if (rrs <= widthsr.and.pps <= widthsm) then
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
      if (s < 3.556896) goto 2
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
        xx=ran(mseed)                                  
        t1=alog((1-xx)*XEXP(a*ta)+xx)/a                    
        c1=1.-t1/ta  
        t1=2.*pi*ran(mseed)                                                                   
        c2=pzr/pr
        s2=sqrt(1.-c2*c2)                                                       
        if (pxr == 0.0.and.pyr == 0.0) goto 3
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
        if (ntag == -1) goto 4                                             
        if (ntag == 0)  goto 2 

c        call paur34(i)
c        if (ntag == -1) goto 4                                             
c        call paur34(j)
c        if (ntag == 0)  goto 2
      
4        pz(i)=pz1                                                              
        pz(j)=pz2                                                              
        px(i)=px1                                                              
        px(j)=px2                                                              
        py(i)=py1                                                              
        py(j)=py2                                                              
        em(i)=em1                                                              
        em(j)=em2

         endif
         endif
2      continue
      endif

1      continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine updtsh
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/aaa/npr,nta,gami,hc,nesc,tpo,ppo
        common x(500),y(500),z(500),id(500),il(500),iev(500),rrp               
     &,rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/bb/np,ntag,x00 
      common/gg/nzta,nzpr   
      common/xx/alpha,beta,gamma,symm     
        common/dd/dt,ntstep,nru,lswitch,iprin,rmax,radj,initi,llt      
        common/deform/beta2                                  
      common/shell/ish,r0,ds,Ri0
 
      r0p=r0*float(npr)**0.3333
      shc=1.+1.6228*ds/r0p
      dr2=3.*beta2**2/(8.*3.1415927)
      dis=0.32*r0
      th2=dr2/dis**2
      
      xpc=0.
      ypc=0.
      zpc=0.

      if (ish == 0) then
      OD3=1./3.
      TD3=2./3.
      CB=0.0058
      CS=0.26

      ip=int(npr-nzpr+0.0001)
      jp=nzpr
      S=CB*((FS(ip)+GS(jp))/(float(npr)/2.)**TD3-CS*float(npr)**OD3)
      Esh=S*(1.-2.*th2)*xexp(-th2)
      Esh=Esh/npr*shc
      else
      if (nru == 1.and.ntstep == 1) then
      print*, 'please input Eshp and Eshn'
      read (*,*) Eshp, Eshn
      Eshp=Eshp/nzpr*shc
      Eshn=Eshn/(npr-nzpr)*shc
      endif
      endif

      do 1 i=1,npr
      xpc=xpc+x(i)/npr
      ypc=ypc+y(i)/npr
      zpc=zpc+z(i)/npr
1      continue

      do 2 i=1,npr
      r2=(x(i)-xpc)**2+(y(i)-ypc)**2+(z(i)-zpc)**2
      r=sqrt(r2)
      if (r == 0.) goto 2

      if (ish == 0) then
      fx=Esh/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(x(i)-xpc)/r
      fy=Esh/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(y(i)-ypc)/r
      fz=Esh/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(z(i)-zpc)/r
      else
       if (iso(i) == 1) then
      fx=Eshp/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(x(i)-xpc)/r
      fy=Eshp/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(y(i)-ypc)/r
      fz=Eshp/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(z(i)-zpc)/r
       else
      fx=Eshn/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(x(i)-xpc)/r
      fy=Eshn/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(y(i)-ypc)/r
      fz=Eshn/ds*xexp((r-r0p)/ds)/(1.+xexp((r-r0p)/ds))**2*(z(i)-zpc)/r
       endif
      endif

      px(i)=px(i)+fx*dt
      py(i)=py(i)+fy*dt
      pz(i)=pz(i)+fz*dt
      
2      continue      

      return
      end
c
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
      DO 1 I=2,9
1      if (N < M(I)) GOTO 2
2      QN=0.6*(M(I)**FD3-M(I-1)**FD3)/(M(I)-M(I-1))
      FS=QN*(N-M(I-1))-0.6*(N**FD3-M(I-1)**FD3)
      RETURN
      END

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
      DO 3 I=2,9
3      if (K < L(I)) GOTO 4
4      QN=0.6*(L(I)**FD3-L(I-1)**FD3)/(L(I)-L(I-1))
      GS=QN*(K-L(I-1))-0.6*(K**FD3-L(I-1)**FD3)
      RETURN
      END

      FUNCTION SWITCH(Ri,Rf,R)
      C0=0.
      C1=0.
      C2=0.
      C3= 10.
      C4=-15.
      C5= 6.
      X0=(R-Rf)/(Ri-Rf)
      SWITCH=C0+C1*X0+C2*X0**2+C3*X0**3+C4*X0**4+C5*X0**5
      RETURN
      END

      FUNCTION Fn(m)
      common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,          
     &rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
      common/skyrme/al,sal,t0,t1,arang,vy,width                            
      common/bb/np,ntag,x00 
      double precision  xxr,yyr,zzr,ppxr,ppyr,ppzr,xxr1,yyr1,zzr1,ppxr1,
     &ppyr1,ppzr1,xxr2,yyr2,zzr2,ppxr2,ppyr2,ppzr2,occupr,occupp
      real Fn
      
      pi=3.1415927
      hb=0.197327
      widthr=width
      widthp=hb**2./widthr/4.
      hcr=sqrt(4.*pi*widthr)
      hcp=sqrt(4.*pi*widthp) 
      widthsr=3.*sqrt(widthr)
      widthsm=3.*sqrt(widthp)
      Fn=0.d0

      do 1 j=1,np
      if (iso(m) == iso(j)) then
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
      if (m == j) then
      Fn=Fn+occupr*occupp/64.
      else
      Fn=Fn+0.5*occupr*occupp/64.
      endif
      endif
1     continue

      return
      end

C************************************************************
      FUNCTION fmnp0(pfn,pfp,xlm)
      qf=(pfn+pfp)/2.
      x=qf/xlm
      A=atan(2.*x)
      B=alog(1.+4.*x**2)
      C=(pfn**3)*(pfp**3)
      if (C <= 1.E-12) THEN
      fmnp0=0.
      RETURN
      ENDIF
      fmnp0=4.*qf*(1.-0.5*A/x)*(pfn**3+pfp**3)-(3.*(pfn**2+pfp**2)+
     &    0.5*xlm**2)*qf**2+qf**4+(0.75*(pfn**2+pfp**2)*xlm**2+
     &    xlm**4/8.-0.375*(pfn**2-pfp**2)**2)*B
      fmnp0=fmnp0*xlm**2/C

      RETURN
      END
C************************************************************
      FUNCTION fmlog(pfn,pfp,eps)
      qf=(pfn+pfp)/2.
      C=(pfn**3)*(pfp**3)
      if (C <= 1.E-12) THEN
      fmlog=0.
      RETURN
      ENDIF
      fmlog=4.*(pfn**3+pfp**3)*qf**3-4.5*(pfn**2+pfp**2)*qf**4+2.*qf**6-
     &      1.125*(pfn**2-pfp**2)**2*qf**2
      fmlog=fmlog*alog(1.+eps*qf**2)*alog(1.+eps*qf**2)/C

      RETURN
      END
C************************************************************
      subroutine xmoms(xlm,x1,y1,z1,x2,y2,z2,sx,sy,sz)
      dimension  U(3),W(3),KS(6),C(6),X(6),KK(2,6),DD(4,7)

cc       Gaussian integration points U(3) and weights W(3)
        U(1)=-0.7746
        U(2)= 0.
        U(3)=-u(1)
        W(1)= 0.555556
        W(2)= 0.888889
        W(3)= w(1)
        LL=1
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
       if (KK(1,j) <= 3)     goto 30
       if (KK(2,j) < KS(j)) goto 40
       j=j-1
       if (j.NE.0) goto 20
       sx=DD(2,1)*DD(1,1)
       sy=DD(3,1)*DD(1,1)
       sz=DD(4,1)*DD(1,1)
      return
 40       KK(2,j)=KK(2,j)+1
       C(j)=C(j)+DD(1,j)*2.
       KK(1,j)=1
 30       kc=KK(1,j)
       x(j)=dd(1,j)*U(kc)+C(j)
       if (j == N) goto 20
       M=j+1
       goto 50

      return
      end

      function sumom(xlm,x1,y1,z1,x2,y2,z2)
      dimension  U(3),W(3),KS(6),C(6),X(6),KK(2,6),DD(2,7)

cc       Gaussian integration points U(3) and weights W(3)
        U(1)=-0.7746
        U(2)= 0.
        U(3)=-u(1)
        W(1)= 0.555556
        W(2)= 0.888889
        W(3)= w(1)
        LL=1
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
       if (KK(1,j) <= 3)     goto 30
       if (KK(2,j) < KS(j)) goto 40
       j=j-1
       if (j.NE.0) goto 20
       sumom=DD(2,1)*DD(1,1)
      return
 40       KK(2,j)=KK(2,j)+1
       C(j)=C(j)+DD(1,j)*2.
       KK(1,j)=1
 30       kc=KK(1,j)
       x(j)=dd(1,j)*U(kc)+C(j)
       if (j == N) goto 20
       M=j+1
       goto 50

      return
      end

      function flow(I,x,x1,y1,z1,x2,y2,z2)
      dimension x(6)
      cc=1.5
      goto(1,2,3,4,5,6), I
 1       flow=x1-cc
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
      cc=1.5
      goto(1,2,3,4,5,6), J
 1       fupp=x1+cc
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
      if (j.ne.6) goto 40
       funsm=cf
      return
 40      if (j.ne.5) goto 30
       funsm=1.
      return
 30      if (j.ne.4)goto 20
       funsm=1.
      return
 20   if (j.ne.3)goto 10
       funsm=1.
      return
 10   if (j.ne.2)goto 5
       funsm=1.
      return
 5       funsm=1.
      end
      
      function fnsx(xlm,j,x,x1,y1,z1,x2,y2,z2)
      dimension x(6)
       cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
       cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
       cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
       cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if (j.ne.6) goto 40
       fnsx=cf
      return
 40      if (j.ne.5) goto 30
       fnsx=1.
      return
 30      if (j.ne.4)goto 20
       fnsx=1.
      return
 20   if (j.ne.3)goto 10
       fnsx=1.
      return
 10   if (j.ne.2)goto 5
       fnsx=1.
      return
 5       fnsx=x1-x(1)
      end            

      function fnsy(xlm,j,x,x1,y1,z1,x2,y2,z2)
      dimension x(6)
       cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
       cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
       cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
       cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if (j.ne.6) goto 40
       fnsy=cf
      return
 40      if (j.ne.5) goto 30
       fnsy=1.
      return
 30      if (j.ne.4)goto 20
       fnsy=1.
      return
 20   if (j.ne.3)goto 10
       fnsy=1.
      return
 10   if (j.ne.2)goto 5
       fnsy=y1-x(2)
      return
 5       fnsy=1.
      end

      function fnsz(xlm,j,x,x1,y1,z1,x2,y2,z2)
      dimension x(6)
       cx1=(x1-x(1))**2+(y1-x(2))**2+(z1-x(3))**2
       cx2=(x2-x(4))**2+(y2-x(5))**2+(z2-x(6))**2
       cx3=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
       cf=xexp(-cx1-cx2)/(1.+cx3/xlm)
      if (j.ne.6) goto 40
       fnsz=cf
      return
 40      if (j.ne.5) goto 30
       fnsz=1.
      return
 30      if (j.ne.4)goto 20
       fnsz=1.
      return
 20   if (j.ne.3)goto 10
       fnsz=z1-x(3)
      return
 10   if (j.ne.2)goto 5
       fnsz=1.
      return
 5       fnsz=1.
      end

      function xmsum(imin,imax,xlm,x1,y1,z1,x2,y2,z2)
       common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &     rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
       common/wavel/sl(500)
       hb=0.197327
       xmsum=0.
       do 11 i1=imin,imax
        xk1=xmm*px(i1)
        yk1=xmm*py(i1)
        zk1=xmm*pz(i1)
          do 22 i2=imin,imax
        if (i2 == i1) goto 22
        xmm=sqrt(sl(i1)+sl(i2))/hb
        xk2=xmm*px(i2)
        yk2=xmm*py(i2)
        zk2=xmm*pz(i2)
        x11=(x1-xk1)*(x1-xk1)+(y1-yk1)*(y1-yk1)+(z1-zk1)*(z1-zk1)
        x22=(x2-xk2)*(x2-xk2)+(y2-yk2)*(y2-yk2)+(z2-zk2)*(z2-zk2)
        pk2=(xk1-xk2)*(xk1-xk2)+(yk1-yk2)*(yk1-yk2)+(zk1-zk2)*(zk1-zk2)
        if (x11+x22 > 30.) then
        cff=0.
        else
        cff=xexp(-x11-x22)/(1.+pk2/xlm)
        endif
        xmsum=xmsum+cff
 22         continue
 11       continue

      return
      end

      subroutine xpart(imin,imax,xlm,x1,y1,z1,x2,y2,z2,sx,sy,sz)
       common x(500),y(500),z(500),id(500),il(500),iev(500),rrp,
     &     rrt,px(500),py(500),pz(500),em(500),iso(500),isy(500),fr(500)
       common/wavel/sl(500)
       hb=0.197327
       sx=0.
       sy=0.
       sz=0.
          do 11 i1=imin,imax
        sum=0.
        xk1=xmm*px(i1)
        yk1=xmm*py(i1)
        zk1=xmm*pz(i1)
          do 22 i2=imin,imax
        if (i2 == i1) goto 22
        xmm=sqrt(sl(i1)+sl(i2))/hb
        xk2=xmm*px(i2)
        yk2=xmm*py(i2)
        zk2=xmm*pz(i2)
        x11=(x1-xk1)*(x1-xk1)+(y1-yk1)*(y1-yk1)+(z1-zk1)*(z1-zk1)
        x22=(x2-xk2)*(x2-xk2)+(y2-yk2)*(y2-yk2)+(z2-zk2)*(z2-zk2)
        pk2=(xk1-xk2)*(xk1-xk2)+(yk1-yk2)*(yk1-yk2)+(zk1-zk2)*(zk1-zk2)
        if (x11+x22 > 30.) then
        cff=0.
        else
        cff=xexp(-x11-x22)/(1.+pk2/xlm)
        endif
        sum=sum+cff
 22         continue
         sx=sx+sum*(x1-xk1)
        sy=sy+sum*(y1-yk1)
        sz=sz+sum*(z1-zk1)
 11       continue

      return
      end