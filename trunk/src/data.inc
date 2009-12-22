integer i,j,k,inum,itype,ipos,nstart,nend,ncell,m,ncol,mx,my
integer n1,n2,nw,no,nf,ndir,ngrid,ndiv
real x,y,z,fix_h,fix_u,fix_v,dx,dy,L,W,slope
real fix_c
real fix_c1,fix_c2,fix_c3,fix_c4

parameter (mz=100,my=1,fix_h=4.0,fix_u=1.3,fix_v=0.0)
parameter (fix_c=0.012)
parameter (L=1000,W=0.001,n1=1,n2=1,nw=1,no=0,nf=3,slope=0.01)

dimension itype(mx),ipos(mx)
dimension x(mx+1,my+1),y(mx+1,my+1),z(mx+1,my+1),ncell(mx+1,my+1)

common/data_1/ inum,itype,ipos,nstart,nend,x,y,z,ncell,m,ncol,dx,dy,ngrid,ndiv