      subroutine intreadxvid(npmax,x,v,np,inf,infv,iior)
!     npmax: for memory alloc, equal to np
!     np: total n?
!     x,v: pos,vel
!     inf,infv: pos and vel input file
!     iior: unused. for backward compatibility
!     what precision? real*8 and int*8?
      real x(3,npmax),v(3,npmax),lambdat
      integer id(npmax) !integer*8
      character inf1*100,inf2*3,inf*(*),infv*(*)

      nrwproc=10 !parallel processes

!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nrwproc)
!$OMP&PRIVATE(nunit,i,j,inf1,inf2)
!$OMP DO SCHEDULE(STATIC,1)
      do k=1,nrwproc !read in parallel

         if(k.eq.1)then
            write(inf2,'(i2.2)')k !2digits int
            inf1=inf(1:lblnk(inf))//'.'//inf2(1:lblnk(inf2)) !lblnk:index of blank char. construct ($inf.%02d)
            write(*,*)inf,inf1,inf2
            nunit=9+k !file id
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale !header
            write(*,*)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
            read(nunit)((x(j,i),j=1,3),i=1,np/4) !read 1/4
            close(nunit)
         else if(k.eq.2)then
            write(inf2,'(i2.2)')k
            inf1=inf(1:lblnk(inf))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((x(j,i),j=1,3),i=1+npmax/4,npmax/2)!read 2/4, no header
            close(nunit)
         else if(k.eq.3)then
            write(inf2,'(i2.2)')k
            inf1=inf(1:lblnk(inf))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((x(j,i),j=1,3),i=1+npmax/2,npmax/4*3)!read 3/4
            close(nunit)
         else if(k.eq.4)then
            write(inf2,'(i2.2)')k
            inf1=inf(1:lblnk(inf))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((x(j,i),j=1,3),i=1+npmax/4*3,npmax) !read 4/4
            close(nunit)
         else if(k.eq.5)then
            write(inf2,'(i2.2)')k-4
            inf1=infv(1:lblnk(infv))//'.'//inf2(1:lblnk(inf2)) !vel file
            nunit=9+k
            write(*,*)inf1
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale !header
            write(*,*)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
            read(nunit)((v(j,i),j=1,3),i=1,npmax/4) 
            close(nunit)
         else if(k.eq.6)then
            write(inf2,'(i2.2)')k-4
            inf1=infv(1:lblnk(infv))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((v(j,i),j=1,3),i=1+npmax/4,npmax/2)
            close(nunit)
         else if(k.eq.7)then
            write(inf2,'(i2.2)')k-4
            inf1=infv(1:lblnk(infv))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((v(j,i),j=1,3),i=1+npmax/2,npmax/4*3)
            close(nunit)
         else if(k.eq.8)then
            write(inf2,'(i2.2)')k-4
            inf1=infv(1:lblnk(infv))//'.'//inf2(1:lblnk(inf2))
            nunit=9+k
            open(nunit,file=inf1,status='old',form='unformatted')
            read(nunit)((v(j,i),j=1,3),i=1+npmax/4*3,npmax)
            close(nunit)
         else if(k.eq.9)then
            write(rvout,'(''id'',i4.4,''.'',i4.4''.'',i2.2)')irun,nu1,k
     $           -8
            outfile=ndir(1:lblnk(ndir))//rvout
            nunit=9+k
            open(nunit,file=outfile,status='old',form='unformatted') !no header
            read(nunit)(id(i),i=1,np/2)
            close(nunit)
         else if(k.eq.10)then
            write(rvout,'(''id'',i4.4,''.'',i4.4''.'',i2.2)')irun,nu1,k
     $           -8
            outfile=ndir(1:lblnk(ndir))//rvout
            nunit=9+k
            open(nunit,file=outfile,status='old',form='unformatted') !no header
            read(nunit)(id(i),i=1+np/2,np)
            close(nunit)
         endif
      enddo
!$OMP END DO
!$OMP END PARALLEL
      end
