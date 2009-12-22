!************************************************************************
!c      Artificial viscosity is added
!c      Along  Shore current is added
!************************************************************************


      program mult
      include 'mult.inc'
      dimension cdum(numsed)
      open(2,file = 'coords_2D.dat')
      open(3,file = 'mult_2D.dat')
      open(8,file = 'Case2D.plt')
      open(11,file = 'restart.out')
!c     open(99,file='diag.out')
      open(97,file='bedvol.out')
      write(97,*) 'VARIABLES = "X","Bed deposit density", "Z0"'
      

!C Read input, set initial conditions, set up grid.
      
      call input 

      open(15,file='xy_Plot.out')
      if(nsed.eq.1)then
      write(15,*)'VARIABLES="X","H","U","Ri","C","CT","net-Z","Z"'
      endif
      if(nsed.eq.2)then
      write(15,*)'VARIABLES="X","H","U","Ri","C1","C2","CT","net-Z","Z"'
      endif
      if(nsed.eq.3)then
      write(15,*)'VARIABLES="X","H","U","Ri","C1","C2","C3","CT",'
      write(15,*)'          "net-Z","Z"'
      endif
      if(nsed.eq.4)then
      write(15,*)'VARIABLES="X","H","U","Ri","C1","C2","C3","C4","CT",'
      write(15,*)'           "net-Z","Z"'
      endif
      
      
      if(ifront .eq. 1) open(13,file = 'front.out')
!C Compute sediment properties.            
      call sedmnt     
       call xy_plot
!C Begin time loop. 
      iplt = 0
      iprt = 0
      do it=1,nt
       t = t + dt
        amax = 0.D0        
!C Compute predictor.
       call predict
!C Compute fluxes.
       do j=1,ncol
       do k=kbeg(j),kend(j)
c        if(k==1)write(99,*) j,zc0(j,k),zc(j,k)
        call bconds(j,k,hp,up,vp,cp)
        h1 = hn(nop(j,k,1))
        h4 = hn(nop(j,k,4))
        dh4 = h4 - h1
        call fluxes(j-1,j,k,k,1,dh4)         ! vert. faces.
        h2 = hn(nop(j,k,2))
        dh1 = h2 - h1              
        call fluxes(j,j,k-1,k,2,dh1)         ! hori. faces.
        do i=1,inum(j,k)
          if(ipos(j,k,i) .eq. 3) then
            h3 = hn(nop(j,k,3))
            dh3 = h4 - h3
            call fluxes(j,j,k,k+1,2,dh3)        ! top (boundary) faces. 
           elseif(ipos(j,k,i) .eq. 2) then
            h3 = hn(nop(j,k,3))
            dh2 = h2 - h3
            call fluxes(j,j+1,k,k,1,dh2)        ! right (boundary) faces
          endif
        enddo
       enddo
       enddo
C Update solution.
       
       do j=1,ncol
       do k=kbeg(j),kend(j)
        call source(j,k,dend)          
        cel = dsqrt(ct(j,k)*h(j,k))
        dx =  deta(j,k,2)*area(j,k) 
	  dy =  dxi(j,k,1)*area(j,k)	  
        asum = (dabs(u(j,k)) + cel)/dx + (dabs(v(j,k)) + cel)/dy  
        if(asum .gt. amax) amax = asum
        do l=1,ndim
         q(j,k,l) = q(j,k,l) + dt/area(j,k)*(f(j,k,l,2)*ds(j,k,2)
     &            + f(j,k,l,1)*ds(j,k,1) - f(j+1,k,l,1)*ds(j+1,k,1)
     &            - f(j,k+1,l,2)*ds(j,k+1,2)) + dt*qsource(j,k,l)   
        enddo
        
        q(j,k,2) = q(j,k,2)/dend
        q(j,k,3) = q(j,k,3)/dend
       enddo
       enddo
       
       call store                  
       call bathym           
       call ARTVIS  

C Interpolate cell-centered quantities to corner nodes.
       call interp  
       iplt = iplt + 1
       iprt = iprt + 1
       if(iplt .eq. nplt) then
          iplt = 0      

          call output
          call xy_plot

C Compute bed deposit density.
          if(ibedvol .eq. 1) call bedvol 
  
       endif

       if(iprt .eq. nprt) then
          iprt = 0
       if(icheck .eq. 1) then
         vol = 0.D0
         do j=1,ncol
         do k=kbeg(j),kend(j)
            vol = vol + h(j,k)*area(j,k)
         enddo
         enddo
       endif
       if(ifront .eq. 1) then
         jf = 0
         do j=1,ncol
          if(h(j,25) .gt. epsh) then 
            jf = j
          endif
         enddo
         if(jf .ne. 0) then
          dum = ct(jf,25)*h(jf,25)
          if(dum .gt. 0.d0) fr = u(jf,25)/dsqrt(dum)
        write(13,'(6e15.5)') t,xc(jf,25)+0.5d0*deta(jf,25,2)*area(jf,25)
     &                       ,h(jf,25),u(jf,25),ct(jf,25),fr 
          endif
         endif
         do i=1,nmon
	     call sample(2,i,xmon(i),ymon(i),hn(i),un(i),vn(i),cdum,z(i))
         enddo
         itest(2) = 1
       endif
       write(*,200) t, amax*dt  
	   
	     
      enddo
C End of time loop.  Write to restart file.
      do j=1,ncol
      do k=kbeg(j),kend(j)
      do i = 1, nsed
      enddo
       write(11,'(2i5,15e15.5)') j,k,h(j,k),u(j,k),v(j,k),
     &  (c(j,k,i),i=1,nsed),(dzt(j,k,i),i=1,nsed)
      enddo
      enddo
      
c     final bathymetry + coordinates
c      call  coord

      stop
 200  format(' time is ',e10.3,'  max CFL # is ',f7.2)
      end
************************************************************************
      subroutine bathym
      include 'mult.inc'

      do j=1,ncol
      do k=kbeg(j),kend(j)
       zc(j,k) = zc0(j,k)
C Compute dz for each grain size for this time step.
       do i=1,nsed
        dzt(j,k,i) = dzt(j,k,i) + dt*0.5D0*
     &           (qsource(j,k,i+3) + qold(j,k,i))/(xlam(i) - 1.D0)     
      
C Compute dzt for each grain size.
        if(dzt(j,k,i) .lt. 0.D0 .and. ierode .eq. 0) dzt(j,k,i) = 0.D0
C Compute new bed elevation.
        zc(j,k) = zc(j,k) + dzt(j,k,i)
       enddo
      enddo
      enddo

      return
      end
************************************************************************
      subroutine bconds(j,k,hdum,udum,vdum,cdum)
      include 'mult.inc'
      dimension hdum(0:nx,0:ny), udum(0:nx,0:ny), vdum(0:nx,0:ny), 
     &          cdum(0:nx,0:ny,numsed)
C Loop over all boundary faces in the cell.
      do i=1,inum(j,k)
       call findbc(i,j,k,jj,kk,j2,k2,jl,kl,id)
       zc(jj,kk) = 2.d0*zc(j,k) - zc(j2,k2)
C Open boundary.
       if(itype(j,k,i) .eq. 0)then
         hdum(jj,kk) = hdum(j,k) 
         udum(jj,kk) = udum(j,k)
         vdum(jj,kk) = vdum(j,k) 
         do ii=1,nsed
           cdum(jj,kk,ii) = cdum(j,k,ii)
         enddo
C Solid wall.
       elseif(itype(j,k,i) .eq. 1)then
        dh(jj,kk,id) = -dh(j,k,id) 
        du(jj,kk,id) = du(j,k,id)*(cn(jl,kl,id)*cn(jl,kl,id) - 
     &                 sn(jl,kl,id)*sn(jl,kl,id))	+
     &                 2.d0*dv(j,k,id)*sn(jl,kl,id)*cn(jl,kl,id)
        dv(jj,kk,id) = dv(j,k,id)*(sn(jl,kl,id)*sn(jl,kl,id) - 
     &                 cn(jl,kl,id)*cn(jl,kl,id))	+
     &                 2.d0*du(j,k,id)*sn(jl,kl,id)*cn(jl,kl,id)
        hdum(jj,kk) = hdum(j,k) 
        udum(jj,kk) = udum(j,k)*(sn(jl,kl,id)*sn(jl,kl,id) -
     &                cn(jl,kl,id)*cn(jl,kl,id)) -
     &                2.D0*vdum(j,k)*sn(jl,kl,id)*cn(jl,kl,id)
        vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,id)*cn(jl,kl,id) -
     &                sn(jl,kl,id)*sn(jl,kl,id)) -
     &                2.D0*udum(j,k)*sn(jl,kl,id)*cn(jl,kl,id)
        do ii=1,nsed
         dc(jj,kk,id,ii) = -dc(j,k,id,ii)
         cdum(jj,kk,ii) = cdum(j,k,ii)
        enddo
C Fixed velocity (subcritical).
       elseif(itype(j,k,i) .eq. 2) then
        hdum(jj,kk) = hdum(j,k)
        udum(jj,kk) = fix(j,k,2)
        vdum(jj,kk) = fix(j,k,3)
        do ii=1,nsed                
         cdum(jj,kk,ii) = cdum(j,k,ii) 
        enddo  
C Fixed depth, velocity (supercritical).             
       elseif(itype(j,k,i) .eq. 3) then
        hdum(jj,kk) = fix(j,k,1)
        udum(jj,kk) = fix(j,k,2)
        vdum(jj,kk) = fix(j,k,3) 
        do ii=1,nsed  
          cdum(jj,kk,ii) = fix(j,k,ii+3) 
        enddo 
C Fixed depth (subcritical).                            
       elseif(itype(j,k,i) .eq. 4) then
        hdum(jj,kk) = fix(j,k,1)
        udum(jj,kk) = udum(j,k)
        vdum(jj,kk) = vdum(j,k)
        do ii=1,nsed
          cdum(jj,kk,ii) = fix(j,k,ii+3)
        enddo
C Fixed flowrate (subcritical).             
       elseif(itype(j,k,i) .eq. 5) then
        hdum(jj,kk) = hdum(j,k)
        udum(jj,kk) = fix(j,k,2)/hdum(jj,kk)
        vdum(jj,kk) = fix(j,k,3)/hdum(jj,kk) 
        do ii=1,nsed                
          cdum(jj,kk,ii)  =  fix(j,k,ii+3)/hdum(jj,kk)
        enddo                               
       endif
      enddo

      return
      end   
************************************************************************
      subroutine fluxes(jl,jr,kl,kr,i1,dhface)
      include 'mult.inc'
      dimension  fdum(numsed),cl(numsed),cr(numsed)
C MUSCL reconstruction at cell face.
      hl = dmax1(0.d0, hp(jl,kl) + 0.5D0*dh(jl,kl,i1))
      hr = dmax1(0.d0, hp(jr,kr) - 0.5D0*dh(jr,kr,i1))
      clt = 0.D0 ; crt = 0.D0 ; ctot = 0.D0
      do i=1,nsed
       cl(i) = dmax1(0.d0, cp(jl,kl,i) + 0.5D0*dc(jl,kl,i1,i))
       cr(i) = dmax1(0.d0, cp(jr,kr,i) - 0.5D0*dc(jr,kr,i1,i))
       ctot = ctot + 0.5D0*(cl(i) + cr(i))
       clt = clt + den(i)*cl(i)
       crt = crt + den(i)*cr(i)
      enddo
      szl = hl*clt
      szr = hr*crt
C Zero fluxes for dry face.
      if(szl + szr .le. 0.d0) then
        do i=1,ndim
          f(jr,kr,i,i1) = 0.d0
        enddo
C Normal face.
       else
        zetal = hp(jl,kl) + zc(jl,kl)
        zetar = hp(jr,kr) + zc(jr,kr)
        if(hl .gt. epsh .and. hr .le. epsh)then
           if(zetal .le. zc(jr,kr)) then
             hr = hl
             do i=1,nsed  
               cr(i) = cl(i)
             enddo
           endif
          elseif(hr .gt. epsh .and. hl .le. epsh) then
           if(zetar .le. zc(jl,kl)) then
             hl = hr
             do i=1,nsed  
               cl(i) = cr(i)
             enddo
           endif
        endif
        ul = up(jl,kl) + 0.5D0*du(jl,kl,i1)
        vl = vp(jl,kl) + 0.5D0*dv(jl,kl,i1)
        ur = up(jr,kr) - 0.5D0*du(jr,kr,i1)
        vr = vp(jr,kr) - 0.5D0*dv(jr,kr,i1)
        call solver(hl,hr,ul,ur,vl,vr,cl,cr,szl,szr,fdum,sn(jr,kr,i1),
     &              cn(jr,kr,i1),ctot)
        hdum = (clt + crt)*dhface**2/48.d0
        f(jr,kr,1,i1) = fdum(1)
        f(jr,kr,2,i1) = fdum(2) + hdum*cn(jr,kr,i1)
        f(jr,kr,3,i1) = fdum(3) + hdum*sn(jr,kr,i1)
        do i=4,ndim                
          f(jr,kr,i,i1) = fdum(i)
        enddo
      endif 
 
      return
      end
************************************************************************
      subroutine predict
      include 'mult.inc'
      dimension dct(2), q4(numsed)
C Compute cell average gradients.
      if(ilim .gt. 0) then
       do j=1,ncol
       do k=kbeg(j),kend(j)
        n1 = nop(j,k,1)
        n2 = nop(j,k,2)
        n3 = nop(j,k,3)
        n4 = nop(j,k,4)
        call bconds(j,k,h,u,v,c)
        if(h(j,k) .gt. epsh) then
        
C Depth.
	  dzb(j,k,1) = 0.5d0*(z(n2) + z(n3) - z(n1) - z(n4))
	  dzb(j,k,2) = 0.5d0*(z(n3) + z(n4) - z(n1) - z(n2))
        zeta = h(j,k) + zc(j,k)
        zetal = h(j-1,k) + zc(j-1,k)
        zetar = h(j+1,k) + zc(j+1,k)
        if(h(j-1,k) .le. epsh) then 
           zetal = zeta
         elseif(h(j+1,k) .le. epsh) then
           zetar = zeta
        endif
        d1 = zeta - zetal
        d2 = zetar - zeta
        dfs1 = grad(ilim, beta, d1, d2)
        dh(j,k,1) = dfs1 - dzb(j,k,1)	
        zetal = h(j,k-1) + zc(j,k-1)
        zetar = h(j,k+1) + zc(j,k+1)
        if(h(j,k-1) .le. epsh) then 
           zetal = zeta
         elseif(h(j,k+1) .le. epsh) then
           zetar = zeta
        endif
        d1 = zeta - zetal
        d2 = zetar - zeta
        dfs2 = grad(ilim, beta, d1, d2)
        dh(j,k,2) = dfs2 - dzb(j,k,2)	
C U velocity.
        d1 = u(j,k) - u(j-1,k)
        d2 = u(j+1,k) - u(j,k)
        du(j,k,1) = grad(ilim, beta, d1, d2)
        d1 = u(j,k) - u(j,k-1)
        d2 = u(j,k+1) - u(j,k)
        du(j,k,2) = grad(ilim, beta, d1, d2)
C V velocity.
        d1 = v(j,k) - v(j-1,k)
        d2 = v(j+1,k) - v(j,k)
        dv(j,k,1) = grad(ilim, beta, d1, d2)
        d1 = v(j,k) - v(j,k-1)
        d2 = v(j,k+1) - v(j,k)
        dv(j,k,2) = grad(ilim, beta, d1, d2)
C Sediment concentration.
        dct(1) = 0.d0 ; dct(2) = 0.d0
        do i=1,nsed
         d1 = c(j,k,i) - c(j-1,k,i)
         d2 = c(j+1,k,i) - c(j,k,i)
         dc(j,k,1,i) = grad(ilim, beta, d1, d2)
         dct(1) = dct(1) + den(i)*dc(j,k,1,i)
         d1 = c(j,k,i) - c(j,k-1,i)
         d2 = c(j,k+1,i) - c(j,k,i)
         dc(j,k,2,i) = grad(ilim, beta, d1, d2)
         dct(2) = dct(2) + den(i)*dc(j,k,2,i)
        enddo
	  uxi = u(j,k)*dxi(j,k,1) + (v(j,k)+ALSC)*dxi(j,k,2)
	  ueta = u(j,k)*deta(j,k,1) + (v(j,k)+ALSC)*deta(j,k,2)
C Predictor.
        hp(j,k) = dmax1(0.d0, h(j,k) - 0.5D0*dt*(
     &    uxi*dh(j,k,1) + h(j,k)*(dxi(j,k,1)*du(j,k,1) + 
     &                            dxi(j,k,2)*dv(j,k,1)) +   
     &    ueta*dh(j,k,2) + h(j,k)*(deta(j,k,1)*du(j,k,2) +
     &                             deta(j,k,2)*dv(j,k,2)))) 
     

        up(j,k) = u(j,k) - 0.5D0*dt*(uxi*du(j,k,1) + ueta*du(j,k,2) +
     &    0.5d0*h(j,k)*(dxi(j,k,1)*dct(1) + deta(j,k,1)*dct(2)) +
     &         ct(j,k)*(dxi(j,k,1)*dfs1 + deta(j,k,1)*dfs2)) 
        vp(j,k) = v(j,k) - 0.5D0*dt*(uxi*dv(j,k,1) + ueta*dv(j,k,2) +
     &    0.5d0*h(j,k)*(dxi(j,k,2)*dct(1) + deta(j,k,2)*dct(2)) +
     &    ct(j,k)*(dxi(j,k,2)*dfs1 + deta(j,k,2)*dfs2))   
        do i=1,nsed
         cp(j,k,i) = dmax1(0.d0, c(j,k,i) - 0.5D0*dt*(
     &              uxi*dc(j,k,1,i) + ueta*dc(j,k,2,i)))
        
        enddo
        else
         hp(j,k) = 0.d0
         up(j,k) = 0.d0
         vp(j,k) = 0.d0
         do i=1,nsed
          cp(j,k,i) = 0.d0
         enddo
        endif
        enddo
        enddo
       else
       do j=1,ncol
       do k=kbeg(j),kend(j)
        call bconds(j,k,h,u,v,c)
        hp(j,k) = h(j,k)  
        up(j,k) = u(j,k) 
        vp(j,k) = v(j,k) 
        do i=1,nsed
         cp(j,k,i) = c(j,k,i) 
        enddo
       enddo
       enddo
      endif

      return
      end
************************************************************************
      subroutine solver(hl,hr,ul,ur,vl,vr,cl,cr,szl,szr,f,
     &                  sinn,coss,ct)
      parameter (nn=10)
      implicit real*8(a-h,o-z)
      parameter (ALSC = 0.0d0)
      common/sol1/ epsh, dt, t
      common/sol2/ ndim
      dimension cl(nn),cr(nn),dc(nn),chat(nn)
      dimension ws(nn),e(nn,nn),a(nn),astar(nn),eps(nn),f(nn),suma(nn)
C Compute averages at cell face.
      hhat  = dsqrt(hr*hl)
      duml  = dsqrt(hl)
      dumr  = dsqrt(hr)
      uhat  = (duml*ul + dumr*ur)/(duml + dumr)
      vhat  = (duml*vl + dumr*vr)/(duml + dumr)
      dct   = 0.D0
      do i=1,ndim-3
       chat(i) = (duml*cl(i) + dumr*cr(i))/(duml + dumr)
       dc(i) = cr(i) - cl(i)
       dct = dct + dc(i)
      enddo      
      ahat = dsqrt(0.5D0*(szl + szr))
c      uperp = uhat*coss + vhat*sinn
       uperp = uhat*coss + ((duml*(vl+ALSC) + dumr*(vr+ALSC))
     &            /(duml + dumr))*sinn


C Compute approximate wave strengths.
      dh = hr - hl
      du = ur - ul
      dv = vr - vl
      dupar = -du*sinn + dv*coss
      duperp = du*coss + dv*sinn
      ws(1) = 0.5D0*(dh - hhat*duperp/ahat + 0.5D0*hhat*dct/ct)
      ws(2) = hhat*dupar     
      ws(3) = 0.5D0*(dh + hhat*duperp/ahat + 0.5D0*hhat*dct/ct)
      do i=4,ndim
        ws(i) = hhat/ct*(0.5D0*chat(i-3)*dc(i-3)/ct - dc(i-3)) 
      enddo
C 1st eigenvector.
      e(1,1) = 1.D0
      e(2,1) = uhat - ahat*coss
      e(3,1) = vhat - ahat*sinn
      do i=4,ndim
        e(i,1) = chat(i-3)
      enddo
C 2nd eigenvector.
      e(1,2) = 0.D0
      e(2,2) = -sinn 
      e(3,2) = coss
      do i=4,ndim
        e(i,2) = 0.D0
      enddo
C 3rd eigenvector.
      e(1,3) = 1.D0
      e(2,3) = uhat + ahat*coss
      e(3,3) = vhat + ahat*sinn
      do i=4,ndim
        e(i,3) = chat(i-3)
      enddo
C Remaining eigenvectors.
      do j=4,ndim
       e(1,j) = 1.D0
       e(2,j) = uhat
       e(3,j) = vhat
       do i=4,ndim
        if(i .eq. j)then
          e(i,j) = -ct
         else
          e(i,j) = 0.D0
        endif
       enddo
      enddo
C Compute eigenvalues.
      a(1) = uperp - ahat
      a(2) = uperp
      a(3) = uperp + ahat
      do i=4,ndim
       a(i) = uperp 
      enddo
C Entropy fix.
      al = dsqrt(szl)
      ar = dsqrt(szr)
      uperpl = ul*coss + (vl+ALSC)*sinn
      uperpr = ur*coss + (vr+ALSC)*sinn
      al1 = uperpl - al
      ar1 = uperpr - ar
      al3 = uperpl + al
      ar3 = uperpr + ar
      eps(1) = dmax1(0.D0, 4.D0*(ar1 - al1))
      eps(2) = 0.D0
      eps(3) = dmax1(0.D0, 4.D0*(ar3 - al3))
      do i=4,ndim
         eps(i) = 0.D0
      enddo
      do k=1,ndim   
       if(dabs(a(k)) .lt. 0.5D0*eps(k)) then
         astar(k) = a(k)*a(k)/eps(k) + 0.25D0*eps(k)
        else
         astar(k) = dabs(a(k))
       endif
      enddo
C Compute flux increments.
      do i=1,ndim
      suma(i) = 0.D0
      do k=1,ndim
       suma(i) = suma(i) + astar(k)*ws(k)*e(i,k)
      enddo
      enddo
C Add flux to appropriate cell faces.
      f(1) = 0.5D0*(f1(hl,uperpl) + f1(hr,uperpr) - suma(1))
      f(2) = 0.5D0*(f2(hl,ul,uperpl,coss,szl) 
     &            + f2(hr,ur,uperpr,coss,szr) - suma(2))
      f(3) = 0.5D0*(f2(hl,vl,uperpl,sinn,szl)   
     &            + f2(hr,vr,uperpr,sinn,szr) - suma(3))
      do i=4,ndim
       f(i) = 0.5D0*(f3(hl,cl(i-3),uperpl)
     &             + f3(hr,cr(i-3),uperpr) - suma(i))
      enddo

      return
      end
************************************************************************
      subroutine source(j,k,dend)
      include 'mult.inc'

C Check if cell is wet.
      if(h(j,k) .gt. epsh) then
C Compute bed slope for each cell.
       n1 = nop(j,k,1) 
       n2 = nop(j,k,2)
       n3 = nop(j,k,3)
       n4 = nop(j,k,4)
       sx = ((z(n2)-z(n4))*(y(n3)-y(n1))-(z(n3)-z(n1))*(y(n2)-y(n4)))
     &     /(2.D0*area(j,k))
       sy = ((z(n3)-z(n1))*(x(n2)-x(n4))-(z(n2)-z(n4))*(x(n3)-x(n1)))
     &     /(2.D0*area(j,k))
     
       d1 = z(n4)*z(n4)+z(n4)*z(n1)-z(n1)*z(n2)-z(n2)*z(n2)
       d2 = z(n1)*z(n1)+z(n1)*z(n2)-z(n2)*z(n3)-z(n3)*z(n3) 
       d3 = z(n2)*z(n2)+z(n2)*z(n3)-z(n3)*z(n4)-z(n4)*z(n4)
       d4 = z(n3)*z(n3)+z(n3)*z(n4)-z(n4)*z(n1)-z(n1)*z(n1)       
       sx2 = (y(n1)*d1+y(n2)*d2+y(n3)*d3+y(n4)*d4)/(6.d0*area(j,k))
       sy2 =-(x(n1)*d1+x(n2)*d2+x(n3)*d3+x(n4)*d4)/(6.d0*area(j,k))
C Compute shear velocity for each cell.
        vmag  = dsqrt(u(j,k)*u(j,k) + v(j,k)*v(j,k))
        ustar = cdrag*u(j,k)*vmag
        vstar = cdrag*(v(j,k)+ALSC)*vmag
      
C Compute quiescent fluid entrainment.
         if(vmag .gt. 1.d-10) then
          ri= h(j,k)*ct(j,k)/(vmag*vmag) 

!       dum2 = 0.075D0*vmag/dsqrt(1.D0 + 718.D0*ri**(2.4D0))

!      dum2= 0.072*(dabs(zc0(j+1,k)-zc0(j,k))/(xc(j+1,k)-xc(j,k)))*vmag 
c       dum2= 0.072*slope*vmag 
         dum2 = 0.00153*vmag/(0.0204+ri)

         else
        dum2 = 0.D0
        endif
               

        zeta = h(j,k) + zc(j,k)        
        dend = 1.d0 + cdrag*vmag*dt/h(j,k)

c        dend = 1.d0 + cdrag*1.43d0*vmag*dt/h(j,k)
c        dend = 1.d0
        
        qsource(j,k,1) =  ew*dum2           
        qsource(j,k,2) = -ct(j,k)*(zeta*sx - sx2)  + 
     &                   fcor*h(j,k)*v(j,k)
        qsource(j,k,3) = -ct(j,k)*(zeta*sy - sy2)  - 
     &                   fcor*h(j,k)*u(j,k)
     
c        qsource(j,k,2) = -ct(j,k)*(zeta*sx - sx2) -(ustar*1.43) + 
c     &                   fcor*h(j,k)*v(j,k)
c        qsource(j,k,3) = -ct(j,k)*(zeta*sy - sy2) -(vstar*1.43) - 
c     &                   fcor*h(j,k)*u(j,k)
     
     
C Compute sediment entrainment/deposition.
        do i=1,nsed
c       cbed = cbeddum*c(j,k,i)
        cbed = cbeddum*c(j,k,i)*2.04
       

c For non-uniform sediment-------modified by sadia
c          cbed = cbeddum*c(j,k,i)*(0.4d0*((dsed(i)/dsg)**1.64)+1.64d0)
c --------------------------------------------             
         if(vfall(i) .gt. 1.d-10) then
	 
	    zu = dsqrt(dabs(ustar) + dabs(vstar))/vfall(i)*fun(i)
          ent = 1.3d-7*(cappa*zu)**5/(1.D0 + 4.3333d-7*(cappa*zu)**5)

c        zu = dsqrt(dsqrt(den(i)*dsed(i)*dsed(i))/vnu)
c     &                     *(dsqrt(abs(ustar))/vfall(i))
c		if(zu.le.5.0)ent = 0.0
c          if(zu.gt.5.0.and.zu.le.13.2)ent = 3.d-12*zu**10*(1.0-(5.0/zu))
c          if(zu.gt.13.2)ent = 0.3

		 
          flux = es*p(i)*ent - cbed 
		 
		if(flux.gt.0.d0.and.dzt(j,k,i).lt.0.d0.and. ierode .eq. 0)then
          flux = 0.d0  
          endif
         
          
!	if(dzt(j,k,i).lt.(-0.5d0))then
!	   flux = 0.d0
!	endif
!
!	if((dzt(j,k,i)+dzt_old(j,k,i)).lt.(-0.5d0))then
!	   flux = 0.d0
!	endif

       
           qsource(j,k,i+3) = vfall(i)*flux


!          ent = (cdrag*1053.d0*vmag**2 - 0.2d0)/(3.5d0*86400.d0)
!          qsource(j,k,i+3) =  ent*es*p(i) - vfall(i)*cbed
          
          endif             
         
        enddo     
        
       else
C Dry cell.
         do i=1,ndim
           qsource(j,k,i) = 0.d0
         enddo   
         dend = 1.d0     
      endif

      return
      end
************************************************************************      
      subroutine store 
      include 'mult.inc'

C Update cells.
      do j=1,ncol
      do k=kbeg(j),kend(j)
        h(j,k) = dmax1(0.d0, q(j,k,1))
        
C Update flow variables in fully wet cells.
        if(h(j,k) .gt. epsh .and. iwet(j,k) == 0) then       
          u(j,k) = q(j,k,2)/h(j,k)
          v(j,k) = q(j,k,3)/h(j,k)
          ct(j,k) = 0.d0
          do i=1,nsed
           c(j,k,i) = dmax1(0.d0, q(j,k,i+3)/h(j,k))
           ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
          enddo
         elseif(h(j,k) .le. epsh .and. iwet(j,k) .eq. 0) then
          u(j,k) = 0.d0
          v(j,k) = 0.d0
          ct(j,k) = 0.d0
          do i=1,nsed
           c(j,k,i) = 0.d0
          enddo
        endif
        do i=4,ndim
          qold(j,k,i-3) = qsource(j,k,i)    
        enddo                              
      enddo                              
      enddo    
        
       
      do j=1,ncol
      do k=kbeg(j),kend(j)
        if(iwet(j,k) .gt. 0) then
         if(iwet(j,k) .eq. 1 .and. iwet(j,k-1) .eq. 0) then        
             u(j,k) = u(j,k-1)
             v(j,k) = v(j,k-1)
             ct(j,k) = 0.d0
             do i=1,nsed
              c(j,k,i) = c(j,k-1,i)
              ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
             enddo
           elseif(iwet(j,k) .eq. 2 .and. iwet(j+1,k) .eq. 0) then
            u(j,k) = u(j+1,k)
             v(j,k) = v(j+1,k)            
             ct(j,k) = 0.d0
             do i=1,nsed
              c(j,k,i) = c(j+1,k,i)
              ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
             enddo
           elseif(iwet(j,k) .eq. 3 .and. iwet(j,k+1) .eq. 0) then
            u(j,k) = u(j,k+1)
             v(j,k) = v(j,k+1)             
             ct(j,k) = 0.d0
             do i=1,nsed
              c(j,k,i) = c(j,k+1,i)
              ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
             enddo
          elseif(iwet(j,k) .eq. 4 .and. iwet(j-1,k) .eq. 0) then        
             u(j,k) = u(j-1,k) 
             v(j,k) = v(j-1,k) 
             ct(j,k) = 0.d0
             do i=1,nsed
              c(j,k,i) = c(j-1,k,i)
              ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
             enddo
           else
             u(j,k) = 0.d0
             v(j,k) = 0.d0
             ct(j,k) = 0.d0
             do i=1,nsed
              c(j,k,i) = 0.d0
             enddo
         endif
         q(j,k,2) = h(j,k)*u(j,k)
         q(j,k,3) = h(j,k)*v(j,k)
         do i=1,nsed
           q(j,k,i+3) = h(j,k)*c(j,k,i)
         enddo
        endif
      enddo                              
      enddo

      return
      end 
      
!************************************************************************      


      SUBROUTINE ARTVIS	
      INCLUDE 'mult.inc'
    
          
      do j=1,ncol
      do k=kbeg(j),kend(j)
    
      if(h(j,k).gt.epsh)then
      
      if (j.eq.1)then
      vneu_xi(j,k)=dabs(h(j+1,k)-h(j,k))/(dabs(h(j+1,k)) + dabs(h(j,k)))
      endif

       if (j.eq.ncol)then
      vneu_xi(j,k) = abs(h(j,k)-h(j-1,k))/(dabs(h(j,k))+ dabs(h(j-1,k)))
       endif

      if(j.gt.1.and. j.lt.ncol)then
      vneu_xi(j,k) = dabs(h(j+1,k)-2*h(j,k) + h(j-1,k))/ 
     &                (dabs(h(j+1,k)) + dabs(h(j,k)) + dabs(h(j-1,k)))
      endif
      
      endif

      enddo
      enddo
!      
!      
      do j=1,ncol
      do k=kbeg(j),kend(j)
      
      if(h(j,k).gt.epsh)then
      
      if (k.eq.kbeg(j))then
      vneu_eta(j,k) =dabs(h(j,k+1)-h(j,k))/(dabs(h(j,k+1))+dabs(h(j,k)))
      endif

       if (k.eq.kend(j))then
      vneu_eta(j,k)=dabs(h(j,k)-h(j,k-1))/(dabs(h(j,k))+ dabs(h(j,k-1)))
       endif

      if(k.gt.kbeg(j).and. k.lt.kend(j))then
      vneu_eta(j,k) =dabs(h(j,k+1)-2*h(j,k) + h(j,k-1))/ 
     &                (dabs(h(j,k+1)) +dabs(h(j,k)) +dabs(h(j,k-1)))
      endif
      
      endif
      
      enddo
      enddo
!      
!      
       do j=1,ncol
       do k=kbeg(j),kend(j)
       
       if(h(j,k).gt.epsh)then
      
       coeff_xi_L(j,k)=Art_vis*dmax1(vneu_xi(j-1,k),vneu_xi(j,k))     
       coeff_xi_R(j,k)=Art_vis*dmax1(vneu_xi(j,k),vneu_xi(j+1,k))
       coeff_eta_L(j,k)=Art_vis*dmax1(vneu_eta(j,k-1),vneu_eta(j,k))
       coeff_eta_R(j,k)=Art_vis*dmax1(vneu_eta(j,k),vneu_eta(j,k+1))  
       
       endif
       
       enddo
       enddo
!      
!       
       do j=1,ncol
       do k=kbeg(j),kend(j)      
	 
	 											   
       if(h(j,k).gt.epsh)then
       h(j,k) = h(j,k)+ (coeff_xi_R(j,k)*(h(j+1,k)-h(j,k))- 
     &                   coeff_xi_L(j,k)*(h(j,k)-h(j-1,k)))+
     &                   (coeff_eta_R(j,k)*(h(j,k+1)-h(j,k))- 
     &                    coeff_eta_L(j,k)*(h(j,k)-h(j,k-1)))
!
!       u(j,k) = u(j,k)+ (coeff_xi_R(j,k)*(u(j+1,k)-u(j,k))- 
!     &                    coeff_xi_L(j,k)*(u(j,k)-u(j-1,k)))+
!     &                           (coeff_eta_R(j,k)*(u(j,k+1)-u(j,k))- 
!     &                           coeff_eta_L(j,k)*(u(j,k)-u(j,k-1)))
!     
!     
!       v(j,k) = v(j,k)+ (coeff_xi_R(j,k)*(v(j+1,k)-v(j,k))- 
!     &                    coeff_xi_L(j,k)*(v(j,k)-v(j-1,k)))+
!     &                           (coeff_eta_R(j,k)*(v(j,k+1)-v(j,k))- 
!     &                           coeff_eta_L(j,k)*(v(j,k)-v(j,k-1)))
!     
!     
!      do i=1,nsed
!      
!      c(j,k,i) = c(j,k,i)+ (coeff_xi_R(j,k)*(c(j+1,k,i)-c(j,k,i))- 
!     &                    coeff_xi_L(j,k)*(c(j,k,i)-c(j-1,k,i)))+
!     &                         (coeff_eta_R(j,k)*(c(j,k+1,i)-c(j,k,i))- 
!     &                           coeff_eta_L(j,k)*(c(j,k,i)-c(j,k-1,i)))
!             
!              ct(j,k)  = ct(j,k) + den(i)*c(j,k,i)
!             enddo
             
!          ct(j,k) = ct(j,k)+ (coeff_xi_R(j,k)*(ct(j+1,k)-ct(j,k))- 
!     &                    coeff_xi_L(j,k)*(ct(j,k)-ct(j-1,k)))+
!     &                         (coeff_eta_R(j,k)*(ct(j,k+1)-ct(j,k))- 
!     &                           coeff_eta_L(j,k)*(ct(j,k)-ct(j,k-1)))
!     
!               zc(j,k) = zc(j,k)+ (coeff_xi_R(j,k)*(zc(j+1,k)-zc(j,k))- 
!     &                    coeff_xi_L(j,k)*(zc(j,k)-zc(j-1,k)))+
!     &                           (coeff_eta_R(j,k)*(zc(j,k+1)-zc(j,k))- 
!     &                           coeff_eta_L(j,k)*(zc(j,k)-zc(j,k-1)))
!     
!     
     
      
      
      endif
      q(j,k,1) = h(j,k)
      q(j,k,2) = h(j,k)*u(j,k)
      q(j,k,3) = h(j,k)*v(j,k)
      do i = 1, nsed
      q(j,k,i+3) = h(j,k)*c(j,k,i)
      enddo
        
      enddo
      enddo 

       RETURN
       END	


!	***	***	***	***	***	***	***	***	*** 

 !      SUBROUTINE coord	
 !      INCLUDE 'mult.inc'
 !      
 !      write(16,*)np,ne 
 !       do i=1,np
 !      write(16,*) x(i),y(i),z(i)
 !    enddo
      
 !     do j=1,ncol
 !     do k=kbeg(j),kend(j)
!	 write(16,*) (nop(j,k,i), i=1,4)
!	enddo
!	enddo
      
 !     do j=1,ncol
 !     do k=kbeg(j),kend(j)
 !      do i= 1, nsed
!	 write(23,*) dzt(j,k,i)
!	 enddo
!	enddo
!	enddo
	
!	 RETURN
!       END
       
!	***	***	***	***	***	***	***	***	*** 

       SUBROUTINE xy_plot	
       INCLUDE 'mult.inc'

c      *************** along x-dir *****************************  
       
      write(15,*)'Zone'
      
	 znet(0,mm)= (zc(0,mm)-zcini(0,mm))  
	 csum=0.0 
	 do i = 1,nsed       
        csum = csum + c(0,mm,i)
        enddo
	 Ri_n(0,mm)= h(0,mm)*1.65d0*9.81d0*csum/(u(0,mm)**2)	
	write(15,*)xc(0,mm),h(0,mm),u(0,mm),Ri_n(0,mm),
     &           (c(0,mm,i),i=1,nsed),csum,znet(0,mm),zc(0,mm)
       do j =1,ncol
	 znet(j,mm)= (zc(j,mm)-zcini(j,mm))
	dum = dsqrt((u(j,mm)**2)+(v(j,mm)**2))
	csum=0.0 
	 do i = 1,nsed       
        csum = csum + c(j,mm,i)
        enddo
	if(dum.gt.0.d0)then	
	 Ri_n(j,mm)= h(j,mm)*1.65d0*9.81d0*csum/(dum*dum)
	endif
	 if(dum.le.0.d0) Ri_n(j,mm)=0.d0	  
      write(15,*)xc(j,mm),h(j,mm),u(j,mm),Ri_n(j,mm),
     &           (c(j,mm,i),i=1,nsed),csum,znet(j,mm),zc(j,mm)
       enddo
    
     
!  55   Format(4f15.4,3f15.10) 
	
	 RETURN
       END
       
!	***	***	***	***	***	***	***	***	*** 
     	

 