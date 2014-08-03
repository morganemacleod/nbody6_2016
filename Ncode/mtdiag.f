c
c diagnostic module for Nbody6/Nbody4
c

      subroutine wrapper(densrad)

c      implicit none
C____________________________________
      INCLUDE 'common6.h'

C____________________________________
      INCLUDE 'mtdiag.h'
      
      integer mti,mtj
      real*8 densrad(3)
      
      do mti=1,Ntot
         do mtj=1,3
            mtx(mtj,mti)=x(mtj,mti)
            mtvx(mtj,mti)=xdot(mtj,mti)
         enddo
         mtm(mti) = body(mti)
      enddo

      mtNTOT=NTOT
      mtNS=N-Ifirst
      mtNB=NTOT-N
      mtIFIRST=Ifirst
      mtN=N
      mtNZERO = NZERO
      mttime=TIME
      mtdr(1)=densrad(1)
      mtdr(2)=densrad(2)
      mtdr(3)=densrad(3)

      mtZMBAR = ZMBAR

      ! MMadd
      write(*,*) "ETAR = ", ETAR, "ETAI = ", ETAI, "ETAU = ", ETAU  
      write(*,*) "BODY1 = ",BODY1
c$$$      write(*,*) "IN WRAPPER", mtNTOT,mtNS,mtNB,mtIFIRST,mtN
c$$$      write(*,*) "IN WRAPPER", NTOT,IFIRST,N
      call diagnostic()

C_______________________________
c      do I=1,N
c         if(body(i).gt.0.01) write(98,*) time,body(i),radius(i)
c      enddo


      return
      end



      subroutine diagnostic()

      implicit none

C____________________________________
      INCLUDE 'mtdiag.h'

      integer i,j,jj,L
  
      integer ibh,flag
      double precision mbh

      double precision xbh(3),vxbh(3)

      double precision R2(mtNMAX)
      real*8 R2cm(mtNMAX)
      integer nsort,jlist(mtNMAX)

      double precision mass,mtot,lagr(11),lagm(11)

      double precision rho0

      double precision eps
      parameter(eps=1.d-11)

      double precision pi
      parameter(pi=3.14159271)

      integer iest
      parameter(iest=8)

C_______________________
C density grid
      integer Ngrid,Nadv,icell(mtNMAX)
      double precision r(mtNMAX),rho(mtNMAX),r2old
      double precision SR(mtNMAX),ST(mtNMAX),alpha(mtNMAX),KR,KT
      double precision V2,Vr,Vr2,Vt2,th(mtNMAX),ph(mtNMAX),aux

      double precision rc,rho2,rc2,mcs,mcb
      integer NCS,NCB,js,jb

      double precision mtots,mtotb,mtotbound
      double precision lagrS(11),lagrB(11),maS,maB


      double precision logmass
      integer logcounter

      double precision lowerM, upperM
      lowerM = 0.2  !modified for stellar evolution 
      upperM = 0.8  !modified for stellar evolution

C_______________________________________________
C lagrangian radii definitions
      lagm(1)=0.01
      lagm(2)=0.02
      lagm(3)=0.05
      lagm(4)=0.10
      lagm(5)=0.20
      lagm(6)=0.50
      lagm(7)=0.80
      lagm(8)=0.90
      lagm(9)=0.95
      lagm(10)=0.98
      lagm(11)=1.d0
     
C____________________________________
C find bh position as most massive body among singles bodies
C avoids binaries fromed ith the BH
      mbh=-99999999.
      ibh=-1
      do i=1,mtN
         if(mtm(i).gt.mbh) then
            ibh=i
            mbh=mtm(i)
         else
         endif
      enddo
           
      flag=0
      if(mbh.le.0.005) flag = 100
      do i=1,mtn
         if(mtm(i).eq.mbh) then
            flag=flag+1
         else
         endif
      enddo

      if(flag.gt.1) then
         write(*,*) "WARNING: are you sure that there is a BH?", flag
         write(*,*) "WARNING: USING CENTER OF DENSITY"
         do j=1,3
            xbh(j)=mtdr(j)
            vxbh(j)=0
         enddo
         mbh=1000000.
         goto 23
      else
      endif
C____________________________________
C save bh position
      do j=1,3
         xbh(j)=mtx(j,ibh)
         vxbh(j)=mtvx(j,ibh)
      enddo

      write(*,*) "Time: ", mttime
      write(*,*) "BH found at: ", mttime,ibh, mbh,xbh,vxbh
c$$$      write(*,*) "mtPARAMS: ", mtNTOT,mtNS,mtNB,mtIFIRST,mtN

 23   continue

C_____________________________________
C recenters snapshot around bh position
      do i=1,mtNTOT
         r2(i)=0
c$$$         write(*,*) r2(i),mtx(1,i),mtx(2,i),mtx(3,i)

         do j=1,3
            mtx(j,i)=mtx(j,i)-xbh(j)
            mtvx(j,i)=mtvx(j,i)-vxbh(j)
            r2(i)=r2(i)+mtx(j,i)**2
         enddo
         
       enddo

C______________________________________
C     Sort square distances of CMs with respect to the bh.
      Nsort=mtNTOT-mtIFIRST+1
      
c$$$      write(*,*) "sorting: ",NSORT,mtiFIRST,R2(1),R2(10)

      do i=1,Nsort
         r2cm(i)=r2(mtiFIRST-1+i)
         Jlist(i)=mtiFIRST-1+i
      enddo
      CALL SORT1(Nsort,R2cm,JList)

C______________________________________________________
C check1 OK
C      goto 111


C_____________________________________
C total mass minus BH.
      mtot=0
      mtots=0.
      mtotb=0.
      mtotbound=0.
      do i=2,Nsort
         if(mtm(Jlist(i)).lt.mbh) then
            mtot=mtot+mtm(Jlist(i))
            if(Jlist(i).le.mtN) then
               mtots=mtots+mtm(Jlist(i))
            else
               mtotb=mtotb+mtm(Jlist(i))
            endif
         else
            mtot=mtot+(mtm(Jlist(i))-mbh) !only binary-bhmass added
            mtotbound=mtotbound+mtm(Jlist(i))-mbh
            write(*,*) "Star bound to BH ",mtm(Jlist(i))-mbh
         endif
c         write(78,*) i,r2cm(i),Jlist(i),mtot,mtm(Jlist(i))
      enddo

      write(*,*) "total mass: ",mtot,mtots,mtotb,mtotbound



c$$$C__________________________
c$$$C write positions and velocities of the closest particles to the bh
c$$$      
c$$$      write(72,*) "time = ", mttime, sqrt(r2cm(5))
c$$$      do i=1,mtNTOT
c$$$         if(r2(i).le.r2cm(5)) then
c$$$            write(72,787) i,mtm(i),sqrt(r2(i)),mtx(1,i),mtx(2,i)
c$$$     &           ,mtx(3,i),mtvx(1,i),mtvx(2,i),mtvx(3,i)
c$$$         else
c$$$         endif
c$$$      enddo
c$$$ 787  format(I6,1X,E14.6,1X,7E16.8)
C______________________________________
C     compute lagrangian radii (system minus BH)
      j=1
      js=1
      jb=1
      maS=0.
      maB=0.
      mass=0.

      do i=2,Nsort              !bh is excluded **
         
         L=Jlist(i)
         if((L.lt.mtIFIRST).or.(L.gt.mtNTOT)) then
            write(*,*) "ERROR IN LINKED LIST!!!"
            call flush()
         else
         endif

         if((j.lt.1).or.(j.gt.11)) then
c            write(*,*) "ERROR IN J VALUE!!!", J,i,L,mass
c            call flush()
            j=11
         else
         endif
         if((js.lt.1).or.(js.gt.11)) then
c            write(*,*) "ERROR IN JS VALUE!!!", Js,i,L,maS,mtotS
c            call flush()
            js=11
         else
         endif
         if((jb.lt.1).or.(jb.gt.11)) then
c            if(abs(maB-mtotB).ge.eps)
c     &           write(*,*) "ERROR IN JB VALUE!!!", Jb,i,L,maB,mtotB
c            call flush()
            jb=11
         else
         endif

C_____________________________________
C check for exclusion of binaries with bh. 
         if(mtm(L).lt.mbh) then
            !total lagr radii
            mass=mass+mtm(L)
            if(mass.ge.(lagm(j)*mtot-eps)) then      
               lagr(j)=sqrt(r2cm(i))
               j=j+1
            else
            endif
            !lagr radii bin/singles
            if(L.le.mtN) then
               maS=maS+mtm(L)
            else
               maB=maB+mtm(L)
            endif

            if(maS.ge.(lagm(jS)*mtotS-eps)) then      
               lagrS(jS)=sqrt(r2cm(i))
               js=js+1
            else
            endif
            if(maB.ge.(lagm(jB)*mtotB-eps)) then      
               lagrB(jB)=sqrt(r2cm(i))
               jB=jB+1
            else
            endif
            
         else
            mass=mass+mtm(L)-mbh
            if(mass.ge.(lagm(j)*mtot-eps)) then      
               lagr(j)=sqrt(r2cm(i))
               j=j+1
            else
            endif
         endif
         
         
      enddo
C_______________________________________
      write(*,44) lagr
 44   format("LRtot = ",11E15.6)
      write(*,43) lagrS
 43   format("LRsing = ",11E15.6)
      write(*,42) lagrB
 42   format("LRbin = ",11E15.6)


C__________________________
C check 3 
      goto 17 


C______________________________________
C compute central density (Casertano & Hut) 
      mass=0
      do i=2,iest+1
         if(i.ge.nsort) then
            mass=mass+mtm(Jlist(i))
            goto 11
         else
            mass=mass+mtm(Jlist(i))
         endif
      enddo
 11   rho0=3.*mass/(4.*pi*r2cm(i)**1.5)

      write(*,*) "central density ", rho0


 17   continue

C___________________________________________
C compute spherical coordinates of the particles
      do i=2,Nsort
         L=Jlist(i)
         aux=sqrt(r2cm(i))
         TH(L)=ACOS(mtx(3,L)/aux)
         PH(L)=ATAN2(mtx(2,L),mtx(1,L))
         IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI 

c$$$         write(*,*) i,L,TH(L),mtx(3,L),aux

      enddo


C______________________________________
C compute density/anisotropy over a grid with sqrt(N) spacing
      Ngrid=sqrt(float(Nsort))!*2.
      if(Ngrid.le.iest) Ngrid=iest


c$$$      do i=1,mtn
c$$$         write(75,*) i,mtvx(1,i) ,mtvx(2,i) ,mtvx(3,i) 
c$$$      enddo

C______________________________________
C init vel dis grid
      do j=1,mtNMAX
         SR(J)=0.
         ST(J)=0.
      enddo

C_______________________________________      
      Nadv=(Nsort-1)/Ngrid

      KR=0.
      KT=0.
      j=1
      jj=0
      mass=0
      r2old=0
      do i=2,Nsort
         
         jj=jj+1
         L=Jlist(i)
         icell(i)=j

         mass=mass+mtm(L)
         
         V2=mtVX(1,L)**2+mtVx(2,L)**2+mtVx(3,L)**2 !V*V
         VR=(mtVx(1,L)*DCOS(PH(L))+mtVx(2,L)*DSIN(PH(L)))
     &        *DSIN(TH(L))+mtVx(3,L)*DCOS(TH(L)) 
         VR2=VR*VR           
         VT2=V2-VR2             !TANGENTIAL VEL ^2

         SR(j)=SR(j)+VR2*mtm(L)
         ST(j)=ST(j)+VT2*mtm(L)
         KR=KR+VR2*mtm(L)*0.5
         KT=KT+VT2*mtm(L)*0.5
         
c$$$         write(*,*) L,V2,mtvx(1,L),mtvx(2,L),mtvx(3,L),PH(L),TH(L)

         if(jj.eq.Nadv) then
            r(j)=sqrt(r2cm(i))
            rho(j)= 3*mass/(4.*pi*(r2cm(i)**1.5-r2old**1.5))
            SR(J)=SR(j)/mass
            ST(J)=ST(j)/mass
            mass=0
            j=j+1
            jj=0
            r2old=r2cm(i)
         else
         endif
      enddo
      
C___________________________________________
C scrive OUTPUT
      write(*,*) "Kinetic energy: KR,KT,2KR/KT = " ,KR,KT,2*KR/KT

      alpha(1)=2-st(1)/sr(1)
      write(76,*) "-1 TIME = ", mttime
      write(76,*) 0.5*r(1),rho(1),sr(1),st(1),alpha(1)
      do j=2,Ngrid
         alpha(j)=2-st(j)/sr(j)
         write(76,*) 0.5*(r(j)+r(j-1)),rho(j),sr(j),st(j),alpha(j)
      enddo

 45   format(5E15.6)
C-------------------------------------------- 



C____________________________________________
C compute core radius as density radius
      rc=0.
      rho2=0.
      do i=2,Nsort
         L=Jlist(i)
         if(mtm(L).le.mbh) then
            rc=rc+sqrt(r2cm(i))*rho(icell(i))*mtm(L)!modified mt 03ago05 
            rho2=rho2+rho(icell(i))*mtm(L)
         else
         endif
      enddo
      rc=rc/rho2

      write(*,*) 'CORE RADIUS (CH84) +hmradius=', rc, lagr(6) 

C___________________________________________________
C compute number of singles and binaries in the core

      mcs=0.
      mcb=0.
      NCS=0.
      NCB=0.

      rc2=rc*rc
      do i=2,Nsort
         L=Jlist(i)
         if(r2cm(i).le.rc2) then
            if(L.le.mtN) then
               NCS=NCS+1
               mcs=mcs+mtm(L)
            else
               NCB=NCB+1
               if(mtm(L).le.mbh) then
                  mcb=mcb+mtm(L)
               else
               endif
            endif
         else
            goto 112
         endif
      enddo
 112  continue
      
      write(*,*) 'NCS, NCB,mcs,mcb = ',NCS,NCB,mcs,mcb
      


C__________________________________________________________________
C compute number of singles and binaries in twice the core radius

      mcs=0.
      mcb=0.
      NCS=0.
      NCB=0.

      rc2=rc*rc*4               !
      do i=2,Nsort
         L=Jlist(i)
         if(r2cm(i).le.rc2) then
            if(L.le.mtN) then
               NCS=NCS+1
               mcs=mcs+mtm(L)
            else
               NCB=NCB+1
               if(mtm(L).le.mbh) then
                  mcb=mcb+mtm(L)
               else
               endif
            endif
         else
            goto 114
         endif
      enddo
 114  continue
      
      write(*,*) 'NCS2, NCB2,mcs2,mcb2 = ',NCS,NCB,mcs,mcb
 

ccccccccccccccccccccccccccccc
c New experimental addition to print out a mass spectrum

c__________________
      do i=1,mtNMASS
         mtMASSvalue(i)=0.1*i
         mtMASSdis(i)=0
      enddo

      do L=1,mtN
         if((mtm(L).le.mbh).and.(L.le.mtN)) then
           j = (mtm(L)*mtNZERO)/0.1+1
           if(j.gt.mtNMASS) j=mtNMASS   !done to fix possible run outs
           mtMASSdis(j)=mtMASSdis(j)+1
           !write(*,*) "deeeeeeeeeb: ", j, mtm(L)*mtNZERO
	 else
         endif
      enddo
      write(*,*) "mass df: ", mtMASSdis

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C computes the sum of the log of the masses for ML estimation of power law
      logmass = 0
      logcounter = 0
      do L=1,mtN
       if(((mtm(L)*mtZMBAR).ge.lowerM).and.((mtm(L)*mtZMBAR).le.upperM))
     & then
            logmass = logmass + log(mtm(L)*mtZMBAR)
            logcounter = logcounter+1
        else
        endif
      enddo
      write(*,*) "NLogMass_SumLogMass_in ", logcounter, 
     &         logmass, lowerM, upperM,  " modified for stellar ev"
     &   ," mtZMBAR = ", mtZMBAR

cccccccccccccccccccccccccccccc
     

 111  continue
      return
      end
