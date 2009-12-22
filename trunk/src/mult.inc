************************************************************************
      parameter (nn=46700,nx=1001,ny=301,numsed=10,mm=1)
************************************************************************
      implicit real*8(a-h,o-z)
      parameter (ALSC = 0.0d0,slope=0.002,Art_vis=0.1d0,Le=0.4d0)
       
      common/ints/ ibed, ibedvol, icheck, ierode, ifall, ifront, ilim, 
     &             iout, istart, itest(2), nbcell, ncol, ne, nmon, np,    
     &             nplt, nprt, nsed, nt   
      common/reals/ beta,cappa,cbeddum,cdrag, d50, es,ew,fcor,grav,ct1, 
     &          h0l,h0r,rp,tmax,u0l,u0r,vnu,v0l,v0r,xsplit,ri0,c0,dsg              
      common/a0/ f(0:nx,0:ny,numsed,2), q(0:nx,0:ny,numsed),
     &           qsource(0:nx,0:ny,numsed), qold(0:nx,0:ny,numsed),
     &           fix(0:nx,0:ny,numsed)
      common/a1/ h(0:nx,0:ny), u(0:nx,0:ny), v(0:nx,0:ny),
     &       hn(nn), un(nn), vn(nn), znet(0:nx,0:ny),Ri_n(0:nx,0:ny),
     &           hp(0:nx,0:ny), up(0:nx,0:ny), vp(0:nx,0:ny),
     &           dh(0:nx,0:ny,2), du(0:nx,0:ny,2), dv(0:nx,0:ny,2),
     &           c(0:nx,0:ny,numsed), cp(0:nx,0:ny,numsed),
     &           dc(0:nx,0:ny,2,numsed), cnn(nn,numsed), ct(0:nx,0:ny),
     &           zini(nn),zcini(0:nx,0:ny),dzt_old(0:nx,0:ny,numsed)
      common/h/  itype(nx,ny,4), ipos(nx,ny,4),	kbeg(nx), kend(nx),
     &           inum(nx,ny), nop(nx,ny,4), iwet(0:nx,0:ny)
      common/i/  x(nn), y(nn), z(nn), xc(0:nx,0:ny), yc(0:nx,0:ny), 
     &           zc(0:nx,0:ny), dzt(0:nx,0:ny,numsed),
     &           z0(nn), zc0(0:nx,0:ny), dz(nn,numsed), bden(0:nx,0:ny),
     &           dzb(nx,ny,2),Rich(nn)
      common/k/  sn(nx,ny,2), cn(nx,ny,2), area(0:nx,0:ny),
     &           ds(nx,ny,2), dxi(0:nx,0:ny,2), deta(0:nx,0:ny,2)
      common/seds/ dsed(numsed), pow(numsed), csf(numsed), qq(numsed),
     &            vfall(numsed), fun(numsed), xlam(numsed), den(numsed),
     &            p(numsed), c0l(numsed), c0r(numsed),FF(numsed) 
      common/sol1/ epsh, dt, t
      common/sol2/ ndim
      common/datar/ w1(nn,2), w2(nn,2), w3(nn,2), w4(nn,2), xmon(100), 
     &             ymon(100),pstar(numsed),pmap(0:nx,0:ny),pp(nn)
      common/datai/ j1(nn,2), k1(nn,2)
      common/art/ vneu_eta(0:nx,0:ny), vneu_xi(0:nx,0:ny),
     &            coeff_xi_L(0:nx,0:ny), coeff_xi_R(0:nx,0:ny),
     &            coeff_eta_L(0:nx,0:ny), coeff_eta_R(0:nx,0:ny)
   
