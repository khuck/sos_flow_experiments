!----------------------------------------------------------------
! (c) Copyright, 2001 by the Regents of the University of California.
! FFTclass: Fourier function class in Math Function module of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LANL, 3/7/01
! Description: This class defines the 3d FFT transformation subject to
!              open or periodic conditions, Fourier Sine transformation,
!              Complex-Complex, Complex-Real, and Real-Complex FFT.
! Comments:
!----------------------------------------------------------------
      module FFTclass
      use Timerclass
      use Transposeclass
      interface fftcrlocal_FFT
        module procedure fftcrlocal1_FFT,fftcrlocal2_FFT
      end interface
      interface fftrclocal_FFT
        module procedure fftrclocal1_FFT,fftrclocal2_FFT
      end interface
      contains
!----------------------------------------------------------------
! FFT for 3D open boundary conditions. 
! The original computational domain is doubled in each dimension
! to apply the FFT for the new domain.
        ! 3_D FFT.
        subroutine fft3d1_FFT(nx,ny,nz,nsizez,nsizey,nsizexy,&
                nsizeyz,ksign,scalex,scaley,scalez,x,xstable,xrtable,&
            ystable,yrtable,nprocrow,commrow,nproccol,commcol,comm2d,&
            myidx,myidy,xout)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer,intent(in) :: nprocrow,commrow,nproccol,commcol,comm2d
        integer,intent(in) :: ksign,myidx,myidy
        double precision, intent(in) :: scalex,scaley,scalez
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
        double precision, dimension(nx/2,nsizey,nsizez), intent(in) &
                            :: x
        double complex, dimension(nz,nsizexy,nsizeyz), intent(out) &
                            :: xout
        double precision, dimension(nx,nsizey) :: tmp1
        double complex, dimension(nx/2+1,nsizey) :: tmp10
        double complex, dimension(ny,nsizexy) :: tmp2
        double complex, dimension(nz,nsizexy) :: tmp3
        integer :: i,j,k
        double precision :: t0
        integer :: ierr,nxx
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0

        call starttime_Timer(t0)

        nxx = nx/2 + 1
        !FFTs along x dimensions: could be a lot of cache miss.
        allocate(x0(nx/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, nx/2
              tmp1(i,j) = x(i,j,k)
            enddo
            do i = nx/2+1, nx
              tmp1(i,j) = 0.0
            enddo
          enddo

          ! FFTs along x dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,nx,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, nx/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(ny/2,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
! yrtable needs to be changed.
! transpose between the x and y dimension.
        call trans3d_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
                     ystable,yrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(ny,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, ny/2
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = ny/2+1,ny
              tmp2(i,j) = (0.0,0.0)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

          do j = 1, nsizexy 
            do i = 1, ny
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(nz/2,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
! xrtable needs to be changed.
! transpose between serial ny (stored in i index) and z.
        call trans3d3_TRANSP(ny,nsizexy,nsizez,nsizeyz,x0,x1,nprocrow,&
                      xstable,xrtable,commrow,myidx,nz/2)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, nz/2
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = nz/2+1,nz
              tmp3(i,j) = (0.0,0.0)
            enddo
          enddo

          !FFT along Z.
          call fftlocal_FFT(ksign,scalez,tmp3,nz,nsizexy)

          do j = 1, nsizexy
            do i = 1, nz
              xout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine fft3d1_FFT

        ! 3_D inverse FFT for open BCs.
        subroutine invfft3d1_FFT(nz,ny,nx,nsizexy,nsizeyz,nsizey,&
                nsizez,ksign,scalex,scaley,scalez,x,xstable,xrtable,&
            ystable,yrtable,nprocrow,commrow,nproccol,commcol,comm2d,&
            myidx,myidy,xout)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer,intent(in) :: nprocrow,commrow,nproccol,commcol,comm2d
        integer,intent(in) :: ksign,myidx,myidy
        double precision, intent(in) :: scalex,scaley,scalez
        integer,dimension(0:nprocrow-1),intent(in) :: xstable,xrtable
        integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
        double complex, dimension(nz,nsizexy,nsizeyz), intent(inout) &
                            :: x
        double precision, dimension(nx/2,nsizey,nsizez), intent(out) &
                            :: xout
        double complex, dimension(nz,nsizexy) :: tmp1
        double complex, dimension(ny,nsizexy) :: tmp2
        double complex, dimension(nx/2+1,nsizey) :: tmp3
        double precision, dimension(nx,nsizey) :: tmp30
        integer :: i,j,k
        double precision :: t0
        integer :: ierr,nxx
        double complex, allocatable, dimension(:,:,:) :: x0
        double complex, allocatable, dimension(:,:,:) :: x1

        call starttime_Timer(t0)

        nxx = nx/2 + 1
        allocate(x0(nz/2,nsizexy,nsizeyz))
        !FFTs along y and z dimensions: could be a lot of cache miss.
        do k = 1, nsizeyz
          do j = 1, nsizexy 
            do i = 1, nz
              tmp1(i,j) = x(i,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftlocal_FFT(ksign,scalez,tmp1,nz,nsizexy)

          do j = 1, nsizexy 
            do i = 1, nz/2
              x0(i,j,k) = tmp1(i,j)
!              x0(i,j,k) = tmp1(i,j)*scalez
            enddo
          enddo
        enddo

        allocate(x1(ny,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
! xrtable needs to be changed.
        call trans3d3_TRANSP(nz/2,nsizexy,nsizeyz,nsizez,x0,x1,nprocrow,&
                      xstable,xrtable,commrow,myidx,ny)

        deallocate(x0)
        allocate(x0(ny/2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, ny
              tmp2(i,j) = x1(i,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

          do j = 1, nsizexy 
            do i = 1, ny/2
              x0(i,j,k) = tmp2(i,j) 
!              x0(i,j,k) = tmp2(i,j)*scaley 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(nxx,nsizey,nsizez))
        call trans3d_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
                     ystable,yrtable,commcol,nsizez)
!        call MPI_BARRIER(commcol,ierr)
        deallocate(x0)

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nxx
              tmp3(i,j) = x1(i,j,k) 
            enddo
          enddo

          call fftcrlocal_FFT(ksign,scalex,tmp3,nx,nsizey,tmp30)

          do j = 1, nsizey
            do i = 1, nx/2
              xout(i,j,k) = tmp30(i,j)
!              xout(i,j,k) = tmp30(i,j)*scalex*2
            enddo
          enddo
        enddo

        deallocate(x1)

        t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

        return
        end subroutine invfft3d1_FFT

!----------------------------------------------------------------
      !used to find the first derivative after FFT.
      ! Subroutine to perform 1D FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftlocal0_FFT(ksign,scale,x,ny,nsizex)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(inout) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           do i = 1, ny
             tempi(i) = x(i,j)
           enddo
           call four1(tempi,ny/2,ksign)
           do i = 1, ny
             x(i,j) = tempi(i)*scale
           enddo
      end do

      end subroutine fftlocal0_FFT

      ! Subroutine to perform 1D FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftlocal_FFT(ksign,scale,x,ny,nsizex)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double complex, dimension(ny,nsizex), intent(inout) :: x
      real*8, dimension(2*ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           do i = 1, ny
             tempi(2*i-1) = real(x(i,j))
             tempi(2*i) = aimag(x(i,j))
           enddo
           call four1(tempi,ny,ksign)
           do i = 1, ny
             x(i,j) = cmplx(tempi(2*i-1),tempi(2*i))*scale
             !x(i,j) = cmplx(tempi(2*i-1),tempi(2*i))
           enddo
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftlocal_FFT

      ! Subroutine to perform 1D real to complex
      ! FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftrclocal1_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(in) :: x
      double complex, dimension(ny/2+1,nsizex), intent(out) :: y
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
         tempi = x(:,j)
         call realft(tempi,ny,ksign)
         y(1,j) = cmplx(tempi(1),0.0)*scale
         y(ny/2+1,j) = cmplx(tempi(2),0.0)*scale
         do i = 2,ny/2
           y(i,j) = cmplx(tempi(2*i-1),tempi(2*i))*scale
         enddo
         !y(1,j) = cmplx(tempi(1),0.0)
         !y(ny/2+1,j) = cmplx(tempi(2),0.0)
         !do i = 2,ny/2
         !  y(i,j) = cmplx(tempi(2*i-1),tempi(2*i))
         !enddo
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftrclocal1_FFT

      subroutine fftrclocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(in) :: x
      double precision, dimension(ny,nsizex), intent(out) :: y
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
         tempi = x(:,j)
         call realft(tempi,ny,ksign)
         y(:,j) = tempi*scale
         !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftrclocal2_FFT

      ! Subroutine to perform 1D complex to real 
      ! FFT along y in 2D array. Here
      ! y is local to each processor.
      subroutine fftcrlocal1_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(out) :: y
      double complex, dimension(ny/2+1,nsizex), intent(in) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           tempi(1) = real(x(1,j))
           tempi(2) = real(x(ny/2+1,j))
           do i = 2, ny/2
             tempi(2*i-1) = real(x(i,j))
             tempi(2*i) = aimag(x(i,j))
           enddo
           call realft(tempi,ny,ksign)
           y(:,j) = tempi*scale*2
           !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftcrlocal1_FFT

      subroutine fftcrlocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(out) :: y
      double precision, dimension(ny,nsizex), intent(in) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      call starttime_Timer(t0)

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           tempi = x(:,j)
           call realft(tempi,ny,ksign)
           y(:,j) = tempi*scale*2
           !y(:,j) = tempi
      end do
      t_mfft_local1 = t_mfft_local1 + elapsedtime_Timer(t0)

      end subroutine fftcrlocal2_FFT

      subroutine realft(data,n,isign)
      integer isign,n
      real*8 data(n)
      integer i,i1,i2,i3,i4,n2p3
      real*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      end subroutine realft

      subroutine four1(data,nn,isign)
      integer isign,nn
      real*8 data(2*nn)
      integer i,istep,j,m,mmax,n
      real*8 tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end subroutine four1

      subroutine sinft(y,n)
      integer n
      real*8 y(n)
      integer j
      real*8 sum,y1,y2
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0
      do 11 j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
11    continue
      call realft(y,n,+1)
      sum=0.0
      y(1)=0.5*y(1)
      y(2)=0.0
      do 12 j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
12    continue
      return
      end subroutine sinft

    end module FFTclass
