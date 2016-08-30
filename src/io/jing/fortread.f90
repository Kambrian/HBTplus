function alloc_file_unit(reset)
integer alloc_file_unit
integer,save::fileno=10
integer, intent(in)::reset
if(reset.eq.1) then
  fileno=10
else
  fileno=fileno+1
  if(fileno.eq.100) fileno=103 !skip 100-102, reserved for stdin,out,err
endif
alloc_file_unit=fileno  
end function

subroutine open_fortran_file(filename,fileno,bigendian_flag,error_stat)
implicit none
character(1024) filename
integer*4:: fileno,bigendian_flag,error_stat
integer alloc_file_unit
!$omp critical (fortran_open) !this is necessary for parallel fortran read. no effect is outside parallel region.
fileno=alloc_file_unit(0)
if(bigendian_flag.NE.0) then
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat,convert='big_endian')
else
open(unit=fileno,FILE=filename,form='unformatted',status='old',iostat=error_stat)
endif
!$omp end critical (fortran_open)
end subroutine

subroutine close_fortran_file(fileno)
implicit none
integer*4 fileno
!$omp critical (fortran_open)
close(fileno)
!$omp end critical (fortran_open)
end subroutine 

subroutine skip_fortran_record(fileno)
integer*4 fileno
read(fileno)
end subroutine skip_fortran_record

subroutine read_fortran_record1(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno
integer*1 arr(arr_len)
read(fileno) arr
end subroutine 

subroutine read_fortran_record2(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno
integer*2 arr(arr_len)
read(fileno) arr
end subroutine 

subroutine read_fortran_record4(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno
integer*4 arr(arr_len)
read(fileno) arr
end subroutine 

subroutine read_fortran_record8(arr,arr_len,fileno)
implicit none
integer*8 arr_len
integer*4 fileno
integer*8 arr(arr_len)
read(fileno) arr
end subroutine 

subroutine read_group_header_int4(b,ngrp,fileno)
implicit none
real*4 b
integer*4 ngrp,fileno
read(fileno) b,ngrp
end subroutine

subroutine read_group_header_int8(b,ngrp,fileno)
implicit none
real*4 b
integer*8 ngrp
integer*4 fileno
read(fileno) b,ngrp
end subroutine

subroutine read_part_header_int4(np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale,fileno)
implicit none
INTEGER*4 ::np,ips,fileno
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
read(fileno) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
end subroutine read_part_header_int4

subroutine read_part_header_int8(np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale,fileno)
implicit none
INTEGER*8 ::np,ips
integer*4 fileno
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
read(fileno) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
end subroutine read_part_header_int8

subroutine read_part_arr_imajor(part_arr,np,fileno)
!arr stored as arr(n,3) on disk, i.e., ID varies first, then dimension
implicit none
integer*8 np
integer*4 fileno
real*4,intent(inout):: part_arr(3,np)
real*4,allocatable:: parr(:,:)
allocate(parr(np,3))
read(fileno) parr
!write(*,*) parr(1,1:3),parr(np,1:3)
part_arr=transpose(parr)
!write(*,*) part_arr(1:3,1),part_arr(1:3,np)
deallocate(parr)
!read(fileno) ((part_arr(i,j),j=1,np),i=1,3)
end subroutine 

subroutine read_part_arr_xmajor(part_arr,np,fileno)
!arr stored as arr(3,n) on disk, i.e., dimension varies first, then particle 
implicit none
integer*8 np
integer*4 fileno
real*4,intent(inout):: part_arr(3,np)
read(fileno) part_arr
end subroutine

subroutine read_pos_file(filename,bigendian_flag,error_stat,header,posarr,np)
implicit none
TYPE:: headstruct
sequence
INTEGER*4 ::np,ips
real*4::ztp,omgt,lbdt,boxsize,xscale,vscale
real*4::hz,vunit,time,mass(2)
END TYPE headstruct
type(headstruct),intent(OUT)::header
character (len=1024) filename
integer*4:: bigendian_flag,error_stat,i,j,np
real*4,intent(out)::posarr(3,np)
real*4,allocatable::parr(:,:),pparr(:,:)
if(bigendian_flag.NE.0) then
open(unit=10,FILE=filename,form='unformatted',status='old',iostat=error_stat,convert='big_endian')
else
open(unit=10,FILE=filename,form='unformatted',status='old',iostat=error_stat)
endif
write(*,*) size(posarr,dim=1),size(posarr,dim=2)
read(10) header%np,header%ips,header%ztp,header%omgt,header%lbdt,header%boxsize,header%xscale,header%vscale
write(*,*) header%np,header%ips,header%ztp,header%omgt,header%lbdt,header%boxsize,header%xscale,header%vscale
allocate(parr(np,3))
read(10) parr
write(*,*) parr(1,1:3),parr(np,1:3)
write(*,*) parr(1:3,1),parr(np-2:np,3)
posarr=transpose(parr)
deallocate(parr)
write(*,*) posarr(1:3,1),posarr(1:3,np)
write(*,*) posarr(1,1:3),posarr(3,np-2:np)
!read(10) ((posarr(i,j),j=1,np),i=1,3)
close(10)
end subroutine read_pos_file
