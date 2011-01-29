************************************************************************      
      subroutine bedvol
      include 'mult.inc'
C This subroutine computes the mass of settled sediment per unit area 
C of bed in gm/cm^2.
      write(97,*)'Zone'
      do j=1,ncol
      do k=kbeg(j),kend(j)
       bden(j,k) = 0.d0
       do i=1,nsed
         bden(j,k) = bden(j,k) + dzt(j,k,i)*
     &               100.d0*(1.D0 - xlam(i))*(den(i)/grav + 1.D0)      
     
       enddo
       if(k.eq.mm)write(97,*)xc(j,k),bden(j,k),zc0(j,k)
      enddo
      enddo

      return
      end 
************************************************************************
      subroutine findbc(i,j,k,jj,kk,j2,k2,jl,kl,id)
      include 'mult.inc'

      if(ipos(j,k,i) .eq. 1) then  
         jj = j
         kk = k - 1
	   j2 = j
	   k2 = k + 1
	   id = 2
         jl = j
         kl = k  
       elseif(ipos(j,k,i) .eq. 2) then
         jj = j + 1
         kk = k 
	   j2 = j - 1
	   k2 = k
	   id = 1
         jl = j + 1 
         kl = k  
       elseif(ipos(j,k,i) .eq. 3) then 
         jj = j
         kk = k + 1
	   j2 = j
	   k2 = k - 1
	   id = 2
         jl = j
         kl = k + 1 
       elseif(ipos(j,k,i) .eq. 4) then 
         jj = j - 1
         kk = k
	   j2 = j + 1
	   k2 = k 
	   id = 1
         jl = j
         kl = k 
      endif

      return
      end 
************************************************************************
      subroutine interp
      include 'mult.inc'
      dimension cdum(nsed)

	do j=1,ncol
	do k=kbeg(j), kend(j)
C Extrapolate to corner cells.
	 ibc1 = 0
	 ibc2 = 0
	 ibc3 = 0
	 ibc4 = 0
	 do i=1, inum(j,k)
	  if(ipos(j,k,i) .eq. 1) ibc1 = 1
	  if(ipos(j,k,i) .eq. 2) ibc2 = 1
	  if(ipos(j,k,i) .eq. 3) ibc3 = 1
	  if(ipos(j,k,i) .eq. 4) ibc4 = 1	  
	 enddo
       if(ibc1 .eq. 1 .and. ibc4 .eq. 1)then
         xc(j-1,k-1) = 2.D0*xc(j,k) - xc(j+1,k+1)
         yc(j-1,k-1) = 2.D0*yc(j,k) - yc(j+1,k+1)
         zc(j-1,k-1) = 2.D0*zc(j,k) - zc(j+1,k+1)
         h(j-1,k-1) = dmax1(0.d0, 2.D0*h(j,k) - h(j+1,k+1))
         u(j-1,k-1) = 2.D0*u(j,k) - u(j+1,k+1)
         v(j-1,k-1) = 2.D0*v(j,k) - v(j+1,k+1)
         do i=1,nsed
           c(j-1,k-1,i) = dmax1(0.d0, 2.D0*c(j,k,i) - c(j+1,k+1,i))
         enddo
       elseif(ibc1 .eq. 1 .and. ibc2 .eq. 1)then
         xc(j+1,k-1) = 2.D0*xc(j,k) - xc(j-1,k+1)
         yc(j+1,k-1) = 2.D0*yc(j,k) - yc(j-1,k+1)
         zc(j+1,k-1) = 2.D0*zc(j,k) - zc(j-1,k+1)
         h(j+1,k-1) = dmax1(0.d0,2.D0*h(j,k) - h(j-1,k+1))
         u(j+1,k-1) = 2.D0*u(j,k) - u(j-1,k+1)
         v(j+1,k-1) = 2.D0*v(j,k) - v(j-1,k+1)
         do i=1,nsed
           c(j+1,k-1,i) = dmax1(0.d0, 2.D0*c(j,k,i) - c(j-1,k+1,i))
         enddo
       elseif(ibc2 .eq. 1 .and. ibc3 .eq. 1)then
         xc(j+1,k+1) = 2.D0*xc(j,k) - xc(j-1,k-1)
         yc(j+1,k+1) = 2.D0*yc(j,k) - yc(j-1,k-1)
         zc(j+1,k+1) = 2.D0*zc(j,k) - zc(j-1,k-1)
         h(j+1,k+1) = dmax1(0.d0,2.D0*h(j,k) - h(j-1,k-1))
         u(j+1,k+1) = 2.D0*u(j,k) - u(j-1,k-1)
         v(j+1,k+1) = 2.D0*v(j,k) - v(j-1,k-1)
         do i=1,nsed
           c(j+1,k+1,i) = dmax1(0.d0, 2.D0*c(j,k,i) - c(j-1,k-1,i))
         enddo
       elseif(ibc3 .eq. 1 .and. ibc4 .eq. 1)then
         xc(j-1,k+1) = 2.D0*xc(j,k) - xc(j+1,k-1)
         yc(j-1,k+1) = 2.D0*yc(j,k) - yc(j+1,k-1)
         zc(j-1,k+1) = 2.D0*zc(j,k) - zc(j+1,k-1)
         h(j-1,k+1) = dmax1(0.d0,2.D0*h(j,k) - h(j+1,k-1))
         u(j-1,k+1) = 2.D0*u(j,k) - u(j+1,k-1)
         v(j-1,k+1) = 2.D0*v(j,k) - v(j+1,k-1)
         do i=1,nsed
           c(j-1,k+1,i) = dmax1(0.d0, 2.D0*c(j,k,i) - c(j+1,k-1,i))
         enddo
       endif
       call bconds(j,k,h,u,v,c)
      enddo
      enddo
      do n=1,np
       call sample(1,n,x(n),y(n),hn(n),un(n),vn(n),cdum,z(n))
       do i=1,nsed
         cnn(n,i) = cdum(i)
       enddo
      enddo
C Flag partially wet cells according to wettest face.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       iwet(j,k) = 0
	 n1 = nop(j,k,1)
	 n2 = nop(j,k,2)
	 n3 = nop(j,k,3)
	 n4 = nop(j,k,4)
       hsum = 0.25d0*(hn(n1) + hn(n2) + hn(n3) + hn(n4))
       hprod = hn(n1)*hn(n2)*hn(n3)*hn(n4)
       if(hsum .gt. 0.d0 .and. hprod .eq. 0.d0) then
         h1 = hn(n1) + hn(n2)
         h2 = hn(n2) + hn(n3)
         h3 = hn(n3) + hn(n4)
         h4 = hn(n4) + hn(n1)
         hmax1 = dmax1(h2,h3,h4)
         hmax2 = dmax1(h1,h3,h4)
         hmax3 = dmax1(h2,h1,h4)
         hmax4 = dmax1(h2,h3,h1)
         if(h1 .gt. hmax1) then
            iwet(j,k) = 1
           elseif(h2 .gt. hmax2) then
            iwet(j,k) = 2
           elseif(h3 .gt. hmax3) then
            iwet(j,k) = 3 
           elseif(h4 .gt. hmax4) then
            iwet(j,k) = 4
           else
            iwet(j,k) = 5
         endif
       endif
	enddo
	enddo
              
      return
      end  
************************************************************************
      function grad(i,beta,d1,d2)
      implicit real*8(a-h,o-z)

      grad = 0.d0
      if(i .eq. 1) then
!C Double minmod.
       if(d1*d2 .gt. 0.D0) then
         grad = fmin1(0.5D0*(d1 + d2), 2.D0*d1, 2.D0*d2)  
       endif
      elseif(i .eq. 2) then
!C Superbee.
       if(d1*d2 > 0.D0) then
         grad = fmin2(fmax2(d1,d2), beta*fmin2(d1,d2)) 
       endif
      endif

      return     
      end
************************************************************************      
      subroutine output 
      include 'mult.inc'

!C Write output to file 'case2D.plt'.
	if(iout .eq. 1) then	
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
      do n=1,np
       csum = 0.d0
        ct1 = 0.d0 
       do i=1,nsed
         csum = csum + cnn(n,i)
         ct1= ct1+ den(i)*cnn(n,i)
       enddo  
       dum2 = (un(n)*un(n))+(vn(n)*vn(n))
       if(dum2.gt.0.d0)then
       Rich(n) = hn(n)*ct1/dum2
       else
       Rich(n) = 0.d0
       endif
       write(8,200)x(n),y(n),hn(n),un(n),vn(n),z(n),(cnn(n,i),i=1,nsed),
     &             csum,Rich(n),(z(n)-z0(n)),pp(n)     
!       write(8,200)x(n),y(n),hn(n),un(n),vn(n),z(n),(cnn(n,i),i=1,nsed),
!     &             csum,Rich(n),(z(n)-zini(n))
      enddo
	elseif(iout .eq. 2) then
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
       do n=1,np
       csum = 0.d0
        ct1 = 0.d0 
       do i=1,nsed
         csum = csum + cnn(n,i)
         ct1= ct1+ den(i)*cnn(n,i)
       enddo  
       dum2 = (un(n)*un(n))+(vn(n)*vn(n))
       if(dum2.gt.0.d0)then
       Rich(n) = hn(n)*ct1/dum2
       else
       Rich(n) = 0.d0
       endif
       write(8,200)x(n),y(n),z(n),un(n),vn(n),hn(n),(cnn(n,i),i=1,nsed),
     &             csum,Rich(n),(z(n)-z0(n)),pp(n)      
!       write(8,200)x(n),y(n),z(n),un(n),vn(n),hn(n),(cnn(n,i),i=1,nsed),
!     &             csum,Rich(n),(z(n)-zini(n))
      enddo
	elseif(iout .eq. 3) then
      write(8,*)'ZONE N= ',np,' E=',ne,' F=FEPOINT', ' ET=QUADRILATERAL'
      do n=1,np
       csum = 0.d0
        ct1 = 0.d0 
       do i=1,nsed
         csum = csum + cnn(n,i)
         ct1= ct1+ den(i)*cnn(n,i)
       enddo  
       dum2 = (un(n)*un(n))+(vn(n)*vn(n))
       if(dum2.gt.0.d0)then
       Rich(n) = hn(n)*ct1/dum2
       else
       Rich(n) = 0.d0
       endif 
        write(8,200) x(n),y(n),hn(n)+z(n),un(n),vn(n),hn(n),z(n),
     &           (cnn(n,i),i=1,nsed), csum,Rich(n),z(n)-z0(n),pp(n)    
!       write(8,200) x(n),y(n),hn(n)+z(n),un(n),vn(n),hn(n),z(n),
!     &             (cnn(n,i),i=1,nsed), csum,Rich(n),(z(n)-zini(n))
      enddo
	endif

	
      do j=1,ncol
      do k=kbeg(j),kend(j)
	 write(8,*) (nop(j,k,i), i=1,4)
	enddo
	enddo

      return
 200  format(' ',15e15.5)
      end   
************************************************************************
      subroutine sample(iflag,i,xdum,ydum,hdum,udum,vdum,cdum,zdum)
	                                     
      include 'mult.inc'
      dimension cdum(nsed), id(4) 
C This subroutine interpolates h, u, and v at 
C user selected points and prints them to output files.
C Compute interpolation weights first time through only.
      if(itest(iflag) .eq. 0) then
C Bracket sampling location.
       do j=1,ncol
       do k=kbeg(j),kend(j)
        n1=nop(j,k,1) ; n2=nop(j,k,2) ; n3=nop(j,k,3) ; n4=nop(j,k,4)
        dx1 = x(n2) - x(n1) ; dy1 = y(n2) - y(n1)
        dx2 = x(n3) - x(n2) ; dy2 = y(n3) - y(n2)
        dx3 = x(n3) - x(n4) ; dy3 = y(n3) - y(n4)
        dx4 = x(n4) - x(n1) ; dy4 = y(n4) - y(n1)
        f1 = dx1*(ydum - y(n1)) - dy1*(xdum - x(n1)) 
        f2 = dx2*(ydum - y(n2)) - dy2*(xdum - x(n2))
        f3 = dx3*(ydum - y(n4)) - dy3*(xdum - x(n4))
        f4 = dx4*(ydum - y(n1)) - dy4*(xdum - x(n1))
        if(f1 >= 0.d0 .and. f2 >= 0.d0 .and. 
     &     f3 <= 0.d0 .and. f4 <= 0.d0) then
           dx1 = 0.5d0*(x(n2) + x(n3) - x(n1) - x(n4))
           dy1 = 0.5d0*(y(n2) + y(n3) - y(n1) - y(n4))
           f1 = (ydum - 0.5d0*(y(n1) + y(n4)))*dx1 - 
     &          (xdum - 0.5d0*(x(n1) + x(n4)))*dy1
           dx2 = 0.5d0*(x(n4) + x(n3) - x(n1) - x(n2))
           dy2 = 0.5d0*(y(n4) + y(n3) - y(n1) - y(n2))
           f2 = (ydum - 0.5d0*(y(n1) + y(n2)))*dx2 - 
     &          (xdum - 0.5d0*(x(n1) + x(n2)))*dy2
           if(f1 >= 0.d0) then
              k1(i,iflag) = k+1
             else
              k1(i,iflag) = k
           endif 
           if(f2 >= 0.d0) then
              j1(i,iflag) = j
             else
              j1(i,iflag) = j+1
           endif 
           goto 10
        endif
	 enddo
	 enddo
C Compute interpolation distances.
10     d1 = dsqrt((xdum - xc(j1(i,iflag),k1(i,iflag)))**2 + 
     &            (ydum - yc(j1(i,iflag),k1(i,iflag)))**2)
       d2 = dsqrt((xdum - xc(j1(i,iflag)-1,k1(i,iflag)))**2 + 
     &            (ydum - yc(j1(i,iflag)-1,k1(i,iflag)))**2)
       d3 = dsqrt((xdum - xc(j1(i,iflag)-1,k1(i,iflag)-1))**2 + 
     &            (ydum - yc(j1(i,iflag)-1,k1(i,iflag)-1))**2)
       d4 = dsqrt((xdum - xc(j1(i,iflag),k1(i,iflag)-1))**2 + 
     &            (ydum - yc(j1(i,iflag),k1(i,iflag)-1))**2)
C Compute weights as distance inverses.
	  if(d1 .gt. 0.D0) then
	    w1(i,iflag) = d1**(-3.5D0)
	   else
	    w1(i,iflag) = 1.D0
	    w2(i,iflag) = 0.D0
	    w3(i,iflag) = 0.D0
	    w4(i,iflag) = 0.D0
	    goto 4
	  endif
	  if(d2 .gt. 0.D0) then
	    w2(i,iflag) = d2**(-3.5D0)
	   else
	    w1(i,iflag) = 0.D0
	    w2(i,iflag) = 1.D0
	    w3(i,iflag) = 0.D0
	    w4(i,iflag) = 0.D0
	    goto 4
	  endif
	  if(d3 .gt. 0.D0) then
	    w3(i,iflag) = d3**(-3.5D0)
	   else
	    w1(i,iflag) = 0.D0
	    w2(i,iflag) = 0.D0
	    w3(i,iflag) = 1.D0
	    w4(i,iflag) = 0.D0
	    goto 4
	  endif
	  if(d4 .gt. 0.D0) then
	    w4(i,iflag) = d4**(-3.5D0)
	   else
	    w1(i,iflag) = 0.D0
	    w2(i,iflag) = 0.D0
	    w3(i,iflag) = 0.D0
	    w4(i,iflag) = 1.D0
	    goto 4
	  endif
      endif
C Interpolate data.
  4	do ii=1,4
        id(ii) = 1
      enddo
      sumw = w1(i,iflag)+w2(i,iflag)+w3(i,iflag)+w4(i,iflag)
	zdum = (w1(i,iflag)*zc(j1(i,iflag),k1(i,iflag)) + 
     &        w2(i,iflag)*zc(j1(i,iflag)-1,k1(i,iflag)) + 
     &        w3(i,iflag)*zc(j1(i,iflag)-1,k1(i,iflag)-1) + 
     &        w4(i,iflag)*zc(j1(i,iflag),k1(i,iflag)-1))/sumw
      if(h(j1(i,iflag),k1(i,iflag)) .le. epsh) then
        sumw = sumw - w1(i,iflag)
        id(1) = 0
      endif
      if(h(j1(i,iflag)-1,k1(i,iflag)) .le. epsh) then
        sumw = sumw - w2(i,iflag)
        id(2) = 0
      endif
      if(h(j1(i,iflag)-1,k1(i,iflag)-1) .le. epsh) then
        sumw = sumw - w3(i,iflag) 
        id(3) = 0
      endif
      if(h(j1(i,iflag),k1(i,iflag)-1) .le. epsh) then
        sumw = sumw - w4(i,iflag)
        id(4) = 0
      endif
      if(sumw <= 0.d0) sumw = 1.d-10
	hdum = (id(1)*w1(i,iflag)*h(j1(i,iflag),k1(i,iflag)) + 
     &        id(2)*w2(i,iflag)*h(j1(i,iflag)-1,k1(i,iflag)) + 
     &        id(3)*w3(i,iflag)*h(j1(i,iflag)-1,k1(i,iflag)-1) + 
     &        id(4)*w4(i,iflag)*h(j1(i,iflag),k1(i,iflag)-1))/sumw
      if(hdum .lt. 0.d0) hdum = 0.d0
	udum = (id(1)*w1(i,iflag)*u(j1(i,iflag),k1(i,iflag)) + 
     &        id(2)*w2(i,iflag)*u(j1(i,iflag)-1,k1(i,iflag)) + 
     &        id(3)*w3(i,iflag)*u(j1(i,iflag)-1,k1(i,iflag)-1) + 
     &        id(4)*w4(i,iflag)*u(j1(i,iflag),k1(i,iflag)-1))/sumw
	vdum = (id(1)*w1(i,iflag)*v(j1(i,iflag),k1(i,iflag)) + 
     &        id(2)*w2(i,iflag)*v(j1(i,iflag)-1,k1(i,iflag)) + 
     &        id(3)*w3(i,iflag)*v(j1(i,iflag)-1,k1(i,iflag)-1) + 
     &        id(4)*w4(i,iflag)*v(j1(i,iflag),k1(i,iflag)-1))/sumw
      do n=1,nsed
	cdum(n) =(id(1)*w1(i,iflag)*c(j1(i,iflag), k1(i,iflag), n) + 
     &          id(2)*w2(i,iflag)*c(j1(i,iflag)-1, k1(i,iflag), n) + 
     &          id(3)*w3(i,iflag)*c(j1(i,iflag)-1, k1(i,iflag)-1, n) + 
     &          id(4)*w4(i,iflag)*c(j1(i,iflag), k1(i,iflag)-1, n))/sumw
      if(cdum(n) .lt. 0.d0) cdum(n) = 0.d0
      enddo
     
      if(iflag .eq. 2) then
	  idev = 13 + i
	  write(idev,99) i, t, hdum, udum, vdum, (cdum(n), n=1,nsed), 
     &                 hdum + zdum
      endif

      return
 99   format(' ',i5,15e15.5)
	end
************************************************************************
      real*8 function f1(h,up)
      implicit real*8(a-h,o-z)
      f1 = h*up
      return
      end
************************************************************************
      real*8 function f2(h,udum,up,dum,sz)
      implicit real*8(a-h,o-z)
      f2 = h*udum*up + 0.5d0*sz*h*dum
      return
      end
************************************************************************
      real*8 function f3(h,c,up)
      implicit real*8(a-h,o-z)
      f3 = h*c*up
      return
      end
************************************************************************
      real*8 function fmin1(a,b,c)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.0D0) then
         fmin1 = -dmin1(dabs(a),dabs(b),dabs(c)) 
        else
         fmin1 = dmin1(a,b,c) 
      endif 
      return
      end
************************************************************************
      real*8 function fmin2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.0D0) then
          fmin2 = -dmin1(dabs(a),dabs(b))
         else
          fmin2 = dmin1(a,b)
      endif 
      return
      end
************************************************************************
      real*8 function fmax2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.0D0) then
          fmax2 = -dmax1(dabs(a),dabs(b))
         else
          fmax2 = dmax1(a,b) 
      endif
      return
      end
