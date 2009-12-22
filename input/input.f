
         program input
         implicit real(a-h,o-z)
!         include 'data.inc'
          
         call Initial_condition
         stop
         end program input
         
!         **********************************************     
         
         Subroutine Initial_condition
!         Include 'data.inc'      
          
         dimension itype(1000),ipos(1000),ncell(1001,2)
         dimension x(1001,2),y(1001,2),z(1001,2)

       open(unit=1,file='mult_1D.dat')
       open(unit=2,file='coords_1D.dat')
       
          mx=100
          my=100
          fix_h=4.0
          fix_u=1.3
          fix_v=0.0
          fix_c1=0.012
          fix_c2=0.0
          fix_c3=0.0
          L=10000.0
          W=10000.0
          n1=1
          n2=2
          nw=1
          no=0
          nf=3
          slope=0.001

       write(1,*)'gravity  nt  tmax   xsplit epsh   beta ilim'
       write(1,*)'9.81d0   3600   600.d0   0.5d0  1.d-5  2.d0   1'
       write(1,*)'drag   kin. vis. of clean water'
       write(1,*)'0.005d0      0.000001d0'
       write(1,*)'number of sediments       std. dev. of mixture'
       write(1,*)'1                           1.d0'
       write(1,*)'sed. dia CSF  Powround. den. ratio  vol. frac'
       write(1,*)'0.00003d0      1.d0   6.d0     1.65d0         1.d0'
       !write(1,*)'0.00004d0      1.d0   6.d0     1.65d0         0.25d0'
       !write(1,*)'0.00006d0      1.d0   6.d0     1.65d0         0.25d0'
       write(1,*)'fluid ent. constant  sediment ent. (1.0=yes)'
      write(1,*)'1.d0                 1.d0'
      write(1,*)'Coriolis constant      near bed conc. factor'
      write(1,*)'0.d0                 1.d0'
      write(1,*)'number of boundary cells  nplt   nprt  iout'
      write(1,*)'396       1800              10          3'
      write(1,*)'istart ibedvol  ierode ifront icheck '
      write(1,*)'0        1         0     1       0'

      write(1,*)'j   k    inum   itype   ipos'

       k=1
       j=1
       itype(1)=nf
       itype(2)=nw
       itype(3)=nw
       inum=3
       ipos(1)=4
       ipos(2)=1
       ipos(3)=3
       write(1,99)j,k,inum,(itype(i),i=1,inum),(ipos(i),i=1,inum)
       do j = 2,(mx-1)
       k=1
      inum=2
      itype(1)=nw
      itype(2)=nw
      ipos(1)=1
      ipos(2)=3
       write(1,100)j,k,inum,(itype(i),i=1,inum),(ipos(i),i=1,inum)
       enddo

        k=1
        j=mx
       itype(1)=no
       itype(2)=nw
       itype(3)=nw
       inum=3
       ipos(1)=2
       ipos(2)=1
       ipos(3)=3
       write(1,99)j,k,inum,(itype(i),i=1,inum),(ipos(i),i=1,inum)

   99 format (1x,9i6)
  100 format (1x,5i6,i12,i6)

       write(1,*)'ncol'
      ncol=mx
       write(1,97)ncol
       write(1,*)'j   nstart  nend'
       nstart=1
      nend=my
       do j= 1,mx
       write(1,98)j,nstart,nend
       enddo
   97 format(i5)
   98 format(1x,3i9)

      write(1,*)'h0l   u0l  v0l  c0l(1)  c0l(2)......c0l(nsed)'
      write(1,*)'0.0   0.0  0.0   0.000  0.000  0.000' 
      write(1,*)'h0r   u0r  v0r  c0r(1)  c0r(2)......c0r(nsed)'
      write(1,*)'0.0   0.0  0.0   0.000  0.000  0.000' 
      write(1,*)'number of fixed bc cells, ndir'
      ndir=n2-n1
      write(1,*)ndir
      write(1,*)'j k fix h fix u fix v fix c(1)fix c(2).....fix c(nsed)'
      k=1
      write(1,96)1,k,fix_h,fix_u,fix_v,fix_c1
   96 format (1x,2i5,4f9.4) 

      write(1,*)'nmon'
      write(1,*)'0'
      write(1,*)'i x y'
      write(1,*)
      write(1,*)
      write(1,*)

      ngrid=(mx+1)*(my+1)
      ndiv=mx*my
      write(2,*)ngrid,ndiv
      x=0.0
      z=0.0
      dx=L/float(mx)
      dy=W/float(my)
      do j=1,(mx+1)
      y=0.0
      do k= 1,(my+1)
      z(j,k)=z(1,k)+(j-1)*dx*(-slope)
      write(2,*)x(j,k),y(j,k),z(j,k)
      y=y+dy
       enddo
       x=x+dx
      enddo

      ncell(1,1)=1
      do j=1,(mx+1)
      do k=2,(my+1)
      ncell(j,k)=ncell(j,k-1)+1
      enddo
      if(j.lt.(mx+1))then
      ncell(j+1,1)=ncell(j,my+1)+1
      endif
      enddo
      do j=1,mx
      do k= 1,my
      write(2,*)ncell(j,k),ncell(j+1,k),ncell(j+1,k+1),ncell(j,k+1)
       enddo
      enddo

      return
      end subroutine Initial_condition























