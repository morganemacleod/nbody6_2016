C______________________________________________
C mt 06 May 2005: saves a snapshot on the disk

      subroutine snap_to_disk()

      
C____________________________
      INCLUDE 'common6.h'

      integer L

      integer Nbinaries         !number of binaries Ntot-N
      
C_____________________________________
C set Nbinaries
      Nbinaries = Ntot - N

C_____________________________________
C write to disk
      write(56,1022) time,Ntot,N 

 1022 format("time ",F14.4,1X,I12,1X,I12)
 1023 format(E15.8,1X,E17.10,1X,E17.10,1X,E17.10,
     &     1X,E17.10,1X,E17.10,1X,E17.10,1X,I7)

      do L=1,NTOT
         write(56,1023) body(L),X(1,L),X(2,L),X(3,L),
     &        XDOT(1,L),XDOT(2,L),XDOT(3,L),name(L)
      enddo

      end

C__________________________________________
C subroutine to save on disk only at programmed time intervals


      subroutine save_snap()

C____________________________
      INCLUDE 'common6.h'

C___________________________ table for saving snapshot
      INCLUDE 'save.h'

      integer iii

      if((time.gt.table(ltable)).and.ltable.lt.ltablemax) then
 
         !loop to avoid multiple writings associated with restart
         do iii=0,10000
            ltable = ltable+1
            if((time.lt.table(ltable)).or.(ltable.ge.ltablemax)) then
               goto 543
            else
            endif
         enddo
 543     continue
         
         call snap_to_disk()
      else
      endif

      end

C__________________________________________
C subroutine to load programmed time intervals


      subroutine load_snap()

C____________________________
      INCLUDE 'common6.h'

C___________________________ table for saving snapshot
      INCLUDE 'save.h'

      L = 1
      open(7,file='snaptable.dat')
      open(56,file='snapdata.dat')
 5    continue                  !read loop
      read (7,*,end=10) table(L)
      L = L + 1
      go to 5
 10   continue
      Ltablemax = L - 1
      Ltable=1

      close(7)

      end

      
