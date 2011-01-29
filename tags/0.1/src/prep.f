************************************************************************
      subroutine grid
      include 'mult.inc'
C Compute grid metrics.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       n1 = nop(j,k,1)
       n2 = nop(j,k,2)
       n3 = nop(j,k,3)
       n4 = nop(j,k,4)
       xc(j,k) = 0.25D0*(x(n1) + x(n2) + x(n3) + x(n4))
       yc(j,k) = 0.25D0*(y(n1) + y(n2) + y(n3) + y(n4))
       zc0(j,k) = 0.25D0*(z0(n1) + z0(n2) + z0(n3) + z0(n4)) 
       zc(j,k) = zc0(j,k)
       zcini(j,k) = 0.25D0*(zini(n1) + zini(n2) + zini(n3) + zini(n4)) 
       dxdxi = 0.5D0*(-x(n1) + x(n2) + x(n3) - x(n4))
       dxdeta = 0.5D0*(-x(n1) - x(n2) + x(n3) + x(n4))
       dydxi = 0.5D0*(-y(n1) + y(n2) + y(n3) - y(n4))
       dydeta = 0.5D0*(-y(n1) - y(n2) + y(n3) + y(n4))
	 area(j,k) = dxdxi*dydeta - dxdeta*dydxi
       if(area(j,k) .le. 0.D0)then
           write(*,*) 'area error in cell ',j,k
                 pause 

           stop
       endif
	 dxi(j,k,1) = dydeta/area(j,k)
	 deta(j,k,1) = -dydxi/area(j,k)
	 dxi(j,k,2) = -dxdeta/area(j,k)
	 deta(j,k,2) = dxdxi/area(j,k)
       do i=1,inum(j,k)
        if(ipos(j,k,i) .eq. 1) then
          area(j,k-1) = area(j,k)
	    dxi(j,k-1,1) = dxi(j,k,1)
	    deta(j,k-1,1) = deta(j,k,1)
	    dxi(j,k-1,2) = dxi(j,k,2)
	    deta(j,k-1,2) = deta(j,k,2)
         elseif(ipos(j,k,i) .eq. 2) then
          area(j+1,k) = area(j,k)
	    dxi(j+1,k,1) = dxi(j,k,1)
	    deta(j+1,k,1) = deta(j,k,1)
	    dxi(j+1,k,2) = dxi(j,k,2)
	    deta(j+1,k,2) = deta(j,k,2)
         elseif(ipos(j,k,i) .eq. 3) then
          area(j,k+1) = area(j,k)
	    dxi(j,k+1,1) = dxi(j,k,1)
	    deta(j,k+1,1) = deta(j,k,1)
	    dxi(j,k+1,2) = dxi(j,k,2)
	    deta(j,k+1,2) = deta(j,k,2)
         elseif(ipos(j,k,i) .eq. 4) then
          area(j-1,k) = area(j,k)
	    dxi(j-1,k,1) = dxi(j,k,1)
	    deta(j-1,k,1) = deta(j,k,1)
	    dxi(j-1,k,2) = dxi(j,k,2)
	    deta(j-1,k,2) = deta(j,k,2)
        endif
       enddo
      enddo
      enddo
C Compute face angles.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       dx = x(nop(j,k,2)) - x(nop(j,k,1))
       dy = y(nop(j,k,2)) - y(nop(j,k,1))     ! Horizontal face.
       ds(j,k,2) = dsqrt(dx*dx + dy*dy)
       sn(j,k,2) = dx/ds(j,k,2)
       cn(j,k,2) = -dy/ds(j,k,2)
       dx = x(nop(j,k,4)) - x(nop(j,k,1))
       dy = y(nop(j,k,4)) - y(nop(j,k,1))     ! Vertical face.
       ds(j,k,1) = dsqrt(dx*dx + dy*dy)
       sn(j,k,1) = -dx/ds(j,k,1)
       cn(j,k,1) = dy/ds(j,k,1)
       do i=1,inum(j,k)
         if(ipos(j,k,i) .eq. 3) then
          dx = x(nop(j,k,3)) - x(nop(j,k,4))
          dy = y(nop(j,k,3)) - y(nop(j,k,4))  ! Top (boundary) faces.
          ds(j,k+1,2) = dsqrt(dx*dx + dy*dy)
          sn(j,k+1,2) = dx/ds(j,k+1,2)
          cn(j,k+1,2) = -dy/ds(j,k+1,2)
         elseif(ipos(j,k,i) .eq. 2) then
          dx = x(nop(j,k,3)) - x(nop(j,k,2))
          dy = y(nop(j,k,3)) - y(nop(j,k,2))  ! Right (boundary) faces.
          ds(j+1,k,1) = dsqrt(dx*dx + dy*dy)
          sn(j+1,k,1) = -dx/ds(j+1,k,1)
          cn(j+1,k,1) = dy/ds(j+1,k,1)
         endif
       enddo
      enddo
      enddo
C Set some things in ghost cells.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       do i=1, inum(j,k)
         call findbc(i,j,k,jj,kk,j2,k2,jl,kl,id)
         area(jj,kk) = area(j,k)
         xc(jj,kk) = 2.D0*xc(j,k) - xc(j2,k2)
         yc(jj,kk) = 2.D0*yc(j,k) - yc(j2,k2)
	 enddo
      enddo
	enddo

      return
      end 
************************************************************************      
      subroutine input
      include 'mult.inc'
      character*72 dum, filename(100)        
!C Read run data from file 'mult.dat'.
      read(3,'(a72)') dum
      read(3,*) grav, nt, tmax, xsplit, epsh, beta, ilim
      dt = tmax/dfloat(nt)
      read(3,'(a72)') dum
      read(3,*) cdrag, vnu
      read(3,'(a72)') dum
      read(3,*) nsed, sigma
      ndim = nsed + 3
      read(3,'(a72)') dum
      d50   = 0.D0
	    dum2 = 0.d0
     
      do i=1,nsed
        read(3,*) dsed(i), csf(i), pow(i), den(i), p(i)
        d50 = d50 + p(i)*dsed(i)
	dum2= dum2+p(i)*log(dsed(i)) 	
        den(i) = grav*den(i)
      enddo
!c      dsg = 10.d0**dum2 
	write(*,*)d50
	dum1 = 0.d0
      do i=1,nsed
	  dum1 = dum1 + (log(dsed(i))-dum2)**2*p(i)
	enddo

       dum3 =2.d0**(-dsqrt(dum1))
       cappa = 1.D0 - 0.288D0*sigma*dum3
	write(*,*)cappa
 
	 read(3,'(a72)') dum	
      read(3,*) ew, es
      read(3,'(a72)') dum
      read(3,*) fcor, cbeddum
      read(3,'(a72)') dum
      read(3,*) nbcell, nplt, nprt, iout
      read(3,'(a72)') dum
      read(3,*) istart, ibedvol, ierode, ifront, icheck	
      read(3,'(a72)') dum
      write(*,*) 'got initial data input'
       
      do ii=1,nbcell
         read(3,*) j, k, inum(j,k),
     &   (itype(j,k,i),i=1,inum(j,k)),(ipos(j,k,i),i=1,inum(j,k))
      enddo
      write(*,*) 'got boundary data'   
      
        
       read(3,'(a72)') dum
      read(3,*) ncol
      if(ncol+1 .gt. nx) then
        write(*,*) 'ncol too big-increase parameter nx to', ncol+1
        stop
      endif
      read(3,'(a72)') dum
      do j=1,ncol       
         read(3,*) idum,kbeg(j),kend(j)
         if(kend(j)+1 .gt. ny) then
           write(*,*) 'must increase parameter ny to ', kend(j)+1
           stop
         endif
      enddo
      if(nsed+3 .gt. numsed) then
        write(*,*) 'nsed too big-increase parameter numsed to', nsed+3
        stop
      endif
      write(*,*) 'got column data'
         
      read(3,'(a72)') dum
      read(3,*) h0l, u0l, v0l, (c0l(i),i=1,nsed)
      read(3,'(a72)') dum
      read(3,*) h0r, u0r, v0r, (c0r(i),i=1,nsed)
      read(3,'(a72)') dum
      read(3,*) ndir
      read(3,'(a72)') dum
      do i=1,ndir
         read(3,*) j, k, fix(j,k,1), fix(j,k,2), fix(j,k,3),
     &             (fix(j,k,ii), ii=4, ndim)      
     
       enddo
      
      read(3,'(a72)') dum
      read(3,*) nmon
      read(3,'(a72)') dum
      do i=1,nmon
      read(3,*) idum, xmon(idum), ymon(idum)
      enddo
     
      read(3,'(a72)') dum
      do i=1,nmon
	 read(3,'(a72)') filename(i)
  	 idev = 13 + i
	 open(idev,file=filename(i))
      enddo
C Read grid data from file 'coords'.
      read(2,*) np,ne  
      pause     
      if(np .gt. nn) then
         write(*,*) 'too many grid points-must increase nn to', np
         stop
      endif
      do i=1,np
        read(2,*) x(i),y(i),z0(i) 
!        z0(i)=z0(i)*0.3048          
      enddo
        
      write(*,*) 'got grid coords'
         do j=1,ncol
      do k=kbeg(j),kend(j)
       read(2,*) (nop(j,k,i),i=1,4)
      enddo
      enddo           
         
        
      write(*,*) 'got connectivity array'     
      
      
!      do i=1,np
!        read(22,*) x(i),y(i),zini(i)       
!      enddo
!
!      do j=1,ncol
!      do k=kbeg(j),kend(j)
!       do i = 1, nsed
!	 read(24,*) dzt_old(j,k,i)
!	 enddo
!	enddo
!	enddo
      
       
      write(*,*) ' '     
      write(*,*) 'data input was successful.'
      write(*,*) 'Enter a 1 to input sediment fall velocity'
      read (*,*) ifall
      if(ifall .eq. 1) then
        do i=1,nsed
         write(*,*) 'enter fall velocity in [m/sec] for sediment #',i
         read(*,*) vfall(i)
        enddo
      endif
      write(*,*) 'Enter a 1 to input bed porosity'
      read (*,*) ibed
      if(ibed .eq. 1) then
         do i=1,nsed
          write(*,*) 'enter bed porosity for sediment #',i
          read(*,*) xlam(i)
         enddo
      endif
      call grid  
      write(*,*)'grid set-up is complete'                          
C Set initial conditions.
      if(istart .eq. 1)then
        write(*,*) 'Enter the restart time'
        read(*,*) t
        dt = (tmax - t)/nt
        do j=1,ncol
        do k=kbeg(j),kend(j)
         read(11,*) jdum,kdum,h(j,k),u(j,k),v(j,k),
     &   (c(j,k,i),i=1,nsed),(dzt(j,k,i),i=1,nsed) 
        enddo
        enddo 
        rewind 11       
        do j=1,ncol
        do k=kbeg(j),kend(j)
        do i=1,nsed
          zc(j,k) = zc0(j,k) + dzt(j,k,i)
        enddo
        enddo
        enddo 
       else 
C Cold start.
        t = 0.D0 
        do j=1,ncol
        do k=kbeg(j),kend(j)
          if(xc(j,k) .le. xsplit)then
           h(j,k) = h0l
           u(j,k) = u0l
           v(j,k) = v0l
           do i=1,nsed
            c(j,k,i) = c0l(i)
           enddo
          else
           h(j,k) = h0r
           u(j,k) = u0r
           v(j,k) = v0r
           do i=1,nsed
            c(j,k,i) = c0r(i)
           enddo
          endif
        enddo
        enddo     
      endif 
C Initialize conservative variables, apply initial BCs.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       q(j,k,1) = h(j,k)
       q(j,k,2) = h(j,k)*u(j,k)
       q(j,k,3) = h(j,k)*v(j,k)
       do i=4,ndim
        q(j,k,i) = h(j,k)*c(j,k,i-3)
       enddo
       ct(j,k) = 0.d0
       do i=1,nsed
         ct(j,k) = ct(j,k) + den(i)*c(j,k,i)
       enddo
      enddo
      enddo 
      write(*,*) 'initial conditions set'
	if(iout .eq. 1) then	
       if(nsed .eq. 1) then
	 write(8,*) 'VARIABLES ="X","Y","H","U","V","Z","C","CT",
     & "Ri","Net-elev","d50"'
       elseif(nsed .eq. 2) then
	 write(8,*) 'VARIABLES ="X","Y","H","U","V","Z","C1","C2",
     &"CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 3) then
	 write(8,*) 'VARIABLES ="X","Y","H","U","V","Z","C1","C2",
     &"C3","CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 4) then
	 write(8,*) 'VARIABLES="X","Y","H","U","V","Z","C1","C2",
     &"C3","C4","CT","Ri","Netelev","d50"'
       endif
	elseif(iout .eq. 2) then
       if(nsed .eq. 1) then
	 write(8,*) 'VARIABLES ="X","Y","Z","U","V","H","C","CT",
     & "Ri","Net-elev","d50"'
      elseif(nsed .eq. 2) then
	 write(8,*) 'VARIABLES ="X","Y","Z","U","V","H","C1","C2", 
     &"CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 3) then
	 write(8,*) 'VARIABLES ="X","Y","Z","U","V","H","C1","C2",
     &"C3","CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 4) then
	 write(8,*) 'VARIABLES="X","Y","Z","U","V","H","C1","C2",
     &"C3","C4","CT","Ri","Netelev","d50"'
       endif
	elseif(iout .eq. 3) then
       if(nsed .eq. 1) then
	 write(8,*) 'VARIABLES ="X","Y","H+Z","U","V","H","Z","C", 
     &"CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 2) then
	 write(8,*) 'VARIABLES ="X","Y","H+Z","U","V","H","Z","C1",
     &"C2","CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 3) then
	 write(8,*) 'VARIABLES ="X","Y","H+Z","U","V","H","Z","C1",
     &"C2","C3","CT","Ri","Net-elev","d50"'
       elseif(nsed .eq. 4) then
	 write(8,*) 'VARIABLES="X","Y","H+Z","U","V","H","Z","C1",
     &"C2","C3","C4","CT","Ri","Netelev","d50"'
       endif
	endif
      call interp
      call output
      itest(1) = 1

      return
      end
************************************************************************
      subroutine sedmnt
      include 'mult.inc'

C Determine the sediment fall velocity.
      write(*,*) ' '
      write(*,*) '************** SEDIMENT/BED PROPERTIES **************'
      write(*,*) ' '
      write(*,*) ' DIAMETER    FALL VELOCITY    REYNOLDS #    POROSITY '
      do i=1,nsed
       if(dsed(i) .eq. 0.D0) then
        vfall(i) = 0.D0 
        fun(i) = 0.D0
       elseif(ifall .ne. 1) then 
        dstar = den(i)*dsed(i)*dsed(i)*dsed(i)/(vnu*vnu)
        dum = dlog10(dstar)
        dum2 = 1.D0 - csf(i)
        dum3 = 1.D0 + (3.5D0 - pow(i))/2.5D0
        r1 =-3.76715D0 + 1.92944D0*dum - 0.09815D0*dum*dum -
     &       0.00575D0*dum*dum*dum + 0.00056D0*dum*dum*dum*dum
        r2 = dlog10(1.D0 - dum2/0.85D0) - dum2**(2.3)*tanh(dum - 4.6D0)
     &     + 0.3D0*(0.5D0 - csf(i))*dum2*dum2*(dum - 4.6D0)
        r3 = (0.65D0 - csf(i)/2.83D0*tanh(dum - 4.6D0))**dum3 
        wstar = r3*10.D0**(r1 + r2)
        vfall(i) = (den(i)*vnu*wstar)**(1.D0/3.D0)
         
C Compute the particle Reynold's number and the entrainment function.
        rp = dsed(i)*dsqrt(den(i)*dsed(i))/vnu
       if(rp .ge. 3.5D0) then
         fun(i) = rp**(0.6d0)*(dsed(i)/d50)**(0.2d0)
        else
         fun(i) = 0.586D0*rp**(1.23d0)*(dsed(i)/d50)**(0.2d0)
        endif
       endif
C Compute the bed porosity.
        if(dsed(i) .le. 1.d-8)then
          xlam(i) = 0.D0
        elseif(ibed .ne. 1) then
          ws = 2.D0 - 0.229D0/(100.D0*dsed(i))**(0.21D0)
          xlam(i) = 1.D0 - ws/(den(i) + 1.D0)
        endif
        write(*,100) dsed(i),vfall(i),rp,xlam(i) 
      enddo
      write(*,*) ' '
      write(*,*)'*****************************************************'
      write(*,*) ' '
      pause 'hit enter to continue'

      return
  100 format(' ',f7.6,3f15.6)
  101 format(' ',f7.6,8x,f15.6)
      end
