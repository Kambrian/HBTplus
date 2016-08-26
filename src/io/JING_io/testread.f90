program read
implicit none
INTEGER*8 ::np,ips,i
integer*8,parameter:: ng=2048,ipsp=24,iseed=-1000,ipos=1000,nu=0
real*4,parameter:: omega0=0.258,lambda0=0.742,delta=0.8,alpha=1.0
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale,zip1,bl
real*4,allocatable:: pos(:,:),vel(:,:)
real*4::Hratio,time,vunit,sqa
integer*8, allocatable:: id(:)

open(10,file='id6610.5000.01',form='unformatted',status='old')
! read(10) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
! write(*,*) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
! write(*,*) 2048**3
! allocate(pos(3,np))
allocate(id(1000))
read(10) (id(i),i=1,1000)
write(*,*) id(1), id(102), id(120), id(121)
! read(10) pos
! write(*,*) minval(pos),maxval(pos)
close(10)
end program