!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! FieldQuantclass: 3D field quantity class in Field module of APPLICATION
!                 layer.
! Version: 1.0-beta
! Author: Ji Qiang, LBNL
! Description: This class defines a 3-D field quantity in the accelerator.
!              The field quantity can be updated at each step. 
! Comments:
!----------------------------------------------------------------
      module FieldQuantclass
        use Timerclass
        use CompDomclass
        use Pgrid2dclass
        use FFTclass
        use Transposeclass
        use PhysConstclass
        use Fldmgerclass
        type FieldQuant
!          private
          !# of mesh points in x and y directions.
          integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal
          !# field quantity array.
          double precision, pointer, dimension(:,:,:) :: FieldQ
        end type FieldQuant
      contains
        !Initialize field class.
        subroutine construct_FieldQuant(this,innx,inny,innz,geom,grid) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz 
        integer :: myid, myidx, myidy, nptot,nproccol,nprocrow
        integer, allocatable, dimension(:,:,:) :: LocalTable

        call getsize_Pgrid2d(grid,nptot,nproccol,nprocrow) 
        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        allocate(LocalTable(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(nproccol.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprocrow.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
  
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal))
        this%FieldQ = 0.0

        deallocate(LocalTable)

        end subroutine construct_FieldQuant
   
        ! set field quantity.
        subroutine set_FieldQuant(this,innx,inny,innz,geom,grid,&
                                  nprx,npry) 
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(inout) :: this
        type (CompDom), intent(in) :: geom
        type (Pgrid2d), intent(in) :: grid
        integer, intent(in) :: innx, inny, innz, nprx, npry 
        integer :: myid, myidx, myidy
        integer, dimension(2,0:nprx-1,0:npry-1)::LocalTable

        call getpost_Pgrid2d(grid,myid,myidy,myidx)

        call getlctabnm_CompDom(geom,LocalTable)

        this%Nx = innx
        this%Ny = inny
        this%Nz = innz
        this%Nxlocal = innx
        if(npry.gt.1) then
          this%Nylocal = LocalTable(2,myidx,myidy)+2
        else
          this%Nylocal = LocalTable(2,myidx,myidy)
        endif
        if(nprx.gt.1) then
          this%Nzlocal = LocalTable(1,myidx,myidy)+2
        else
          this%Nzlocal = LocalTable(1,myidx,myidy)
        endif
        this%Nxlocal = innx
  
        deallocate(this%FieldQ) 
        allocate(this%FieldQ(this%Nxlocal,this%Nylocal,this%Nzlocal)) 
        this%FieldQ = 0.0

        end subroutine set_FieldQuant

        !find the E and B fields in the lab frame 
        !from the potential on the grid of beam frame.
        subroutine gradEB_FieldQuant(innx,inny,innz,temppotent,ptsgeom,&
              grid,Flagbc,gammaz,FlagImage,egxout,egyout,egzout,bgxout,&
              bgyout,bgzout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx, inny, innz, Flagbc,FlagImage
        type (CompDom), intent(in) :: ptsgeom
        double precision,dimension(innx,inny,innz),intent(inout) :: temppotent
        double precision,dimension(innx,inny,innz),intent(inout) :: egxout,&
               egyout,egzout
        double precision,dimension(innx,inny,innz),intent(inout) :: bgxout,&
               bgyout,bgzout
        type (Pgrid2d), intent(in) :: grid
        double precision, intent(in) :: gammaz
        double precision, dimension(3) :: msize
        double precision :: hxi, hyi, hzi
        double precision :: betC
        double precision,dimension(innx,inny,innz) :: egx,&
               egy,egz
        integer :: totnp,nproccol,nprocrow,myid,myidx,myidy
        integer :: i, j, k, yadd, zadd, innp
!        integer :: comm2d,commcol,commrow

        call getsize_Pgrid2d(grid,totnp,nproccol,nprocrow)
        call getpost_Pgrid2d(grid,myid,myidy,myidx)
        if(nproccol.gt.1) then
          yadd = 1
        else
          yadd = 0
        endif
        if(nprocrow.gt.1) then
          zadd = 1
        else
          zadd = 0
        endif

        call getmsize_CompDom(ptsgeom,msize)
        !get the real length for calculation Ex,Ey,Ez
        hxi = 1.0/(msize(1)*Scxlt)
        hyi = 1.0/(msize(2)*Scxlt)
        !gammaz is due to relativistic factor
        hzi = 1.0/(msize(3)*Scxlt*gammaz)

        !call MPI_BARRIER(comm2d,ierr)
        !if(myid.eq.16) then
        !  print*,"before guardexch:"
        !endif

        ! Exchange potential information on guard grid cells.
        if((Flagbc.eq.1)) then ! 3D open 
          if(totnp.ne.1) then
            call guardexch1_Fldmger(temppotent,innx,inny,innz,grid)
          endif
        else
          print*,"no such boundary condition!!!!"
          stop
        endif

! cache optimization here.
        !print*,"yadd: ",yadd,zadd,innx,inny,innz,1.0/hxi,1.0/hyi,1.0/hzi
        !Ex
        !egx = 0.0
        do k = zadd+1, innz-zadd
          do j = yadd+1, inny-yadd
            egx(1,j,k) = hxi*(1.5d0*temppotent(1,j,k)-2.0d0*temppotent(2,j,k)+&
                              0.5d0*temppotent(3,j,k))
            do i = 2, innx-1
              egx(i,j,k) = 0.5*hxi*(temppotent(i-1,j,k)- &
                           temppotent(i+1,j,k))
            enddo
            egx(innx,j,k) = hxi*(-0.5d0*temppotent(innx-2,j,k)+&
                                  2.0d0*temppotent(innx-1,j,k)- &
                                  1.5d0*temppotent(innx,j,k))
          enddo
        enddo

        !Ey
        !egy = 0.0
        if(nproccol.gt.1) then 
          if(myidy.eq.0) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                egy(i,yadd+1,k) = hyi*(1.5d0*temppotent(i,yadd+1,k)- &
                                       2*temppotent(i,yadd+2,k)+ &
                                       0.5d0*temppotent(i,yadd+3,k) )
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+2, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else if(myidy.eq.(nproccol-1)) then
            do k = zadd+1, innz-zadd
              do i = 1, innx
                egy(i,inny-yadd,k) = hyi*(-0.5d0*temppotent(i,inny-yadd-2,k)+ &
                                          2*temppotent(i,inny-yadd-1,k)- &
                                          1.5d0*temppotent(i,inny-yadd,k))
              enddo
            enddo
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd-1
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          else
            do k = zadd+1, innz-zadd
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                               temppotent(i,j+1,k))
                enddo
              enddo
            enddo
          endif
        else
          do k = zadd+1, innz-zadd
            do i = 1, innx
              egy(i,1,k) = hyi*(1.5d0*temppotent(i,1,k)-2*temppotent(i,2,k)+&
                                0.5d0*temppotent(i,3,k) )
            enddo
            do j = 2, inny-1
              do i = 1, innx
                egy(i,j,k) = 0.5*hyi*(temppotent(i,j-1,k)- &
                             temppotent(i,j+1,k))
              enddo
            enddo
            do i = 1, innx
              egy(i,inny,k) = hyi*(-0.5*temppotent(i,inny-2,k)+&
                                    2*temppotent(i,inny-1,k)- &
                                   1.5d0*temppotent(i,inny,k))
            enddo
          enddo
        endif

        !Ez
        !egz = 0.0
        if(nprocrow.gt.1) then 
          if((Flagbc.eq.1)) then ! 3D open
            if(myidx.eq.0) then
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,zadd+1) = hzi*(1.5d0*temppotent(i,j,zadd+1)- &
                                         2*temppotent(i,j,zadd+2)+ &
                                         0.5d0*temppotent(i,j,zadd+3) )
                enddo
              enddo
              do k = zadd+2, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            else if(myidx.eq.(nprocrow-1)) then
              do k = zadd+1, innz-zadd-1
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
              do j = yadd+1, inny-yadd
                do i = 1, innx
                  egz(i,j,innz-zadd) = hzi*(-0.5d0*temppotent(i,j,innz-zadd-2)+&
                                             2*temppotent(i,j,innz-zadd-1)- &
                                             1.5d0*temppotent(i,j,innz-zadd))
                enddo
              enddo
            else
              do k = zadd+1, innz-zadd
                do j = yadd+1, inny-yadd
                  do i = 1, innx
                    egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                 temppotent(i,j,k+1))
                  enddo
                enddo
              enddo
            endif
          else
            print*,"no such boundary condition!!!"
            stop
          endif
        else
          if((Flagbc.eq.1)) then ! 3D open
            do j = yadd+1, inny-yadd
              do i = 1, innx
                !egz(i,j,1) = hzi*(temppotent(i,j,1)-temppotent(i,j,2))
                egz(i,j,1) = hzi*(1.5d0*temppotent(i,j,1)-2*temppotent(i,j,2)+&
                                  0.5d0*temppotent(i,j,3))
              enddo
            enddo
          else
          endif
          do k = 2, innz-1
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,k) = 0.5*hzi*(temppotent(i,j,k-1)- &
                                      temppotent(i,j,k+1))
              enddo
            enddo
          enddo
          if(Flagbc.eq.1) then
            do j = yadd+1, inny-yadd
              do i = 1, innx
                egz(i,j,innz) = hzi*(-0.5d0*temppotent(i,j,innz-2)+&
                                      2*temppotent(i,j,innz-1)- &
                                      1.5d0*temppotent(i,j,innz))
              enddo
            enddo
          else
          endif
        endif

        !Send the E field to the neibhoring guard grid to do the CIC
        !interpolation.
        if(totnp.ne.1) then
          call boundint4_Fldmger(egx,egy,egz,innx,inny,innz,grid)
        endif

        betC = -sqrt(gammaz**2-1.0)/gammaz/Clight
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              egxout(i,j,k) = egxout(i,j,k) + gammaz*egx(i,j,k)
              bgyout(i,j,k) = bgyout(i,j,k) + gammaz*egx(i,j,k)*betC
            enddo
          enddo
        enddo
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              egyout(i,j,k) = egyout(i,j,k) + gammaz*egy(i,j,k)
              bgxout(i,j,k) = bgxout(i,j,k) - gammaz*egy(i,j,k)*betC
            enddo
          enddo
        enddo
        egzout = egzout + egz
        bgzout = 0.0

        !print*,"grad: ",sum(egx),sum(egy),sum(egz),sum(temppotent),&
        !                hxi,hyi,hzi,innx,inny,innz,myidx,myidy

        end subroutine gradEB_FieldQuant

!----------------------------------------------------------------------
! update potential (solving Possion's equation) with 3D isolated 
! boundary conditions.
        subroutine update3Ot_FieldQuant(this,source,fldgeom,grid,nxlc,&
          nylc,nzlc,nprocrow,nproccol,nylcr,nzlcr,gammaz)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nxlc,nylc,nzlc,&
                               nprocrow,nproccol,nylcr,nzlcr
        type (CompDom), intent(in) :: fldgeom
        double precision, dimension(nxlc,nylc,nzlc), intent(in) :: source
        type (Pgrid2d), intent(in) :: grid
        type (FieldQuant), intent(inout) :: this
        double precision, intent(in) :: gammaz
        double precision, dimension(3) :: msize
        double precision :: hx, hy, hz, temp
        integer, dimension(0:nprocrow-1) :: ypzstable,pztable
        integer, dimension(0:nproccol-1) :: xpystable,pytable
        double precision , dimension(nxlc,nylcr,nzlcr) :: rho
        integer :: myid,myidz,myidy,&
                   comm2d,commcol,commrow
        integer :: nxpylc2,nypzlc2
        integer :: nsxy1,nsxy2,nsyz1,nsyz2
        integer :: i,j,k,inxglb,inyglb,inzglb,innx,inny,innz
        integer, dimension(2,0:nprocrow-1,0:nproccol-1)::LocalTable
        integer, dimension(3) :: glmshnm
        integer :: jadd,kadd,ierr
        double precision :: t0

        call getpost_Pgrid2d(grid,myid,myidy,myidz)
        call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)

        call getlctabnm_CompDom(fldgeom,LocalTable)
        
        do i = 0, nprocrow-1
          pztable(i) = LocalTable(1,i,0)
        enddo
        do i = 0, nproccol-1
          pytable(i) = LocalTable(2,0,i)
        enddo

        call getmnum_CompDom(fldgeom,glmshnm)
        inxglb = glmshnm(1) 
        inyglb = glmshnm(2)
        inzglb = glmshnm(3)

        innz = nzlcr
        inny = nylcr
        innx = nxlc

        if(nprocrow.gt.1) then
          kadd = 1
        else
          kadd = 0
        endif
        if(nproccol.gt.1) then
          jadd = 1
        else
          jadd = 0
        endif
        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              rho(i,j,k) = source(i,j+jadd,k+kadd)
            enddo
          enddo
        enddo
        
        ! +1 is from the real to complex fft.
        nsxy1 = (inxglb+1)/nproccol
        nsxy2 = (inxglb+1) - nproccol*nsxy1
        do i = 0, nproccol-1
          if(i.le.(nsxy2-1)) then
            xpystable(i) = nsxy1+1
          else
            xpystable(i) = nsxy1
          endif
        enddo

        nsyz1 = 2*inyglb/nprocrow
        nsyz2 = 2*inyglb - nprocrow*nsyz1
        do i = 0, nprocrow-1
          if(i.le.(nsyz2-1)) then
            ypzstable(i) = nsyz1+1
          else
            ypzstable(i) = nsyz1
          endif
        enddo

        nxpylc2 = xpystable(myidy)
        nypzlc2 = ypzstable(myidz)

        call getmsize_CompDom(fldgeom,msize)
        hx = msize(1)*Scxlt
        hy = msize(2)*Scxlt
        !gammaz is due to the relativistic effect
        hz = msize(3)*gammaz*Scxlt

        ! Open boundary conditions!
        call openBC3D(innx,inny,innz,rho,hx,hy,hz, &
        nxpylc2,nypzlc2,myidz,myidy,nprocrow,nproccol,commrow,commcol,&
        comm2d,pztable,pytable,ypzstable,xpystable,&
        inxglb,inyglb,inzglb)

        do k = 1, innz
          do j = 1, inny
            do i = 1, innx
              this%FieldQ(i,j+jadd,k+kadd) = rho(i,j,k)*hx*hy*hz
            enddo
          enddo
        enddo

        call MPI_BARRIER(comm2d,ierr)
        t_field = t_field + elapsedtime_Timer(t0)

        end subroutine update3Ot_FieldQuant

        ! Solving Poisson's equation with open BCs.
        subroutine openBC3D(innx,inny,innz,rho,hx,hy,hz,&
           nxpylc2,nypzlc2,myidz,myidy,npz,npy,commrow,commcol,comm2d, &
           pztable,pytable,ypzstable,xpystable,inxglb,&
           inyglb,inzglb)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innx,inny,innz,inxglb,inyglb,inzglb
        integer, intent(in) :: nxpylc2,nypzlc2
        integer, intent(in) :: myidz,myidy,npz,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npz-1), intent(in) :: pztable,ypzstable
        integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
        double precision, dimension(innx,inny,innz), intent(inout) :: rho
        double precision, intent(in) :: hx, hy, hz
        double precision :: scalex,scaley,scalez
        integer :: i,j,k,n1,n2,n3,nylc22,nzlc22
        double complex, allocatable, dimension(:,:,:) :: rho2out
        double complex, allocatable, dimension(:,:,:) :: grn
        integer :: ginny,ginnz,ierr

        n1 = 2*inxglb
        n2 = 2*inyglb
        n3 = 2*inzglb

        nylc22 = nxpylc2
        nzlc22 = nypzlc2

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0

        allocate(rho2out(n3,nylc22,nzlc22))

        call fft3d1_FFT(n1,n2,n3,innz,inny,nylc22,nzlc22,&
                1,scalex,scaley,scalez,rho,ypzstable,pztable,&
        xpystable,pytable,npz,commrow,npy,commcol,comm2d,myidz,myidy, &
        rho2out)

        !c compute FFT of the Green function on the grid:
        ! here the +1 is from the unsymmetry of green function
        ! on double-sized grid.
        if(myidz.eq.(npz-1)) then
           ginnz = innz + 1
        else
           ginnz = innz
        endif
        if(myidy.eq.(npy-1)) then
           ginny = inny + 1
        else
           ginny = inny
        endif
        allocate(grn(n3,nylc22,nzlc22))
        call greenf1t(inxglb,inyglb,inzglb,ginnz,ginny,nylc22,nzlc22, &
               hx,hy,hz,myidz,npz,commrow,myidy,npy,commcol,comm2d,&
                  ypzstable,pztable,xpystable,pytable,grn)

        ! multiply transformed charge density and transformed Green 
        ! function:
        do k = 1, nzlc22
          do j = 1, nylc22
            do i = 1, n3
              rho2out(i,j,k) = rho2out(i,j,k)*grn(i,j,k)
            enddo
          enddo
        enddo

        deallocate(grn)

        ! inverse FFT:
        scalex = 1.0/float(n1)
        scaley = 1.0/float(n2)
        scalez = 1.0/float(n3)
        call invfft3d1_FFT(n3,n2,n1,nylc22,nzlc22,inny,innz,&
               -1,scalex,scaley,scalez,rho2out,pztable,ypzstable,&
        pytable,xpystable,npz,commrow,npy,commcol,comm2d,myidz,myidy,&
        rho)

        deallocate(rho2out)

        end subroutine openBC3D

        ! green function for extended array in time domain.
        subroutine greenf1t(nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz,&
                  hx,hy,hz,myidx,npx,commrow,myidy,npy,commcol,comm2d,&
                   xstable,xrtable,ystable,yrtable,grnout)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nx,ny,nz,nsizez,nsizey,nsizexy,nsizeyz
        integer, intent(in) :: myidx,myidy,npx,npy,commrow,commcol,&
                               comm2d
        integer, dimension(0:npx-1),intent(in) :: xstable,xrtable
        integer, dimension(0:npy-1),intent(in) :: ystable,yrtable
        double precision, intent(in) :: hx, hy, hz
        double complex, intent (out), &
                       dimension (2*nz,nsizexy,nsizeyz) :: grnout
        integer, dimension(0:npx-1) :: gxrtable
        integer, dimension(0:npy-1) :: gyrtable
        integer :: i,j,k,kk,ii,jj,iii,jjj,kkk,n1,n2,n3,ns1,ns2
        integer :: nblockx,nblocky,ksign,nxx,ierr
        double precision, dimension (nx+1,nsizey,nsizez) :: grn
        double precision :: scalex,scaley,scalez
        double precision :: t0
        double precision, dimension(2*nx,nsizey) :: tmp1
        double complex, dimension(nx+1,nsizey) :: tmp10
        double complex, dimension(2*ny,nsizexy) :: tmp2
        double complex, dimension(2*nz,nsizexy) :: tmp3
        double complex, allocatable, dimension(:,:,:) :: x1
        double complex, allocatable, dimension(:,:,:) :: x0
        double precision :: recfourpi

        call starttime_Timer(t0)

        recfourpi = 1.0/(8.0*asin(1.0))
        gxrtable = xrtable
        gxrtable(npx-1) = xrtable(npx-1) + 1
        nblockx = 0
        do i = 0, myidx-1
          nblockx = nblockx + gxrtable(i)
        enddo

        gyrtable = yrtable
        gyrtable(npy-1) = yrtable(npy-1) + 1
        nblocky = 0
        do i = 0, myidy-1
          nblocky = nblocky + gyrtable(i)
        enddo

        do k = 1, nsizez
          do j = 1, nsizey
            do i = 1, nx+1
              jj = j + nblocky
              kk = k + nblockx
                kkk = kk - 1
                jjj = jj - 1
                iii = i - 1
              if((iii*iii+jjj*jjj+kkk*kkk) .ne. 0) then
               grn(i,j,k)=recfourpi/sqrt((hx*iii)**2+(hy*jjj)**2+(hz*kkk)**2)
              endif
            enddo
          enddo
        enddo
        if((myidx.eq.0).and.(myidy.eq.0)) then
          if(nsizez.gt.1) then
            grn(1,1,1) = grn(1,1,2)
          else
            grn(1,1,1) = 1.0
          endif
        endif

        scalex = 1.0
        scaley = 1.0
        scalez = 1.0
        n1 = 2*nx
        n2 = 2*ny
        n3 = 2*nz
        ksign = 1

        nxx = n1/2 + 1
        !FFTs along y and z dimensions.
        allocate(x0(n1/2+1,nsizey,nsizez))
        do k = 1, nsizez
          do j = 1, nsizey 
            do i = 1, n1/2+1
              tmp1(i,j) = grn(i,j,k)
            enddo
            do i = n1/2+2, n1
              tmp1(i,j) = grn(n1-i+2,j,k)
            enddo
          enddo

          ! FFTs along z dimensions:
          call fftrclocal_FFT(ksign,scalex,tmp1,n1,nsizey,tmp10)

          do j = 1, nsizey 
            do i = 1, n1/2+1
              x0(i,j,k) = tmp10(i,j)
            enddo
          enddo
        enddo

        allocate(x1(n2/2+1,nsizexy,nsizez))

        ! FFTs along y dimensions:
!        call MPI_BARRIER(commcol,ierr)
        call trans3d_TRANSP(nxx,n2/2+1,nsizexy,nsizey,x0,x1,npy,&
                     ystable,gyrtable,commcol,nsizez)
        deallocate(x0)
        allocate(x0(n2,nsizexy,nsizez))

        do k = 1, nsizez
          do j = 1, nsizexy 
            do i = 1, n2/2+1
              tmp2(i,j) = x1(i,j,k)
            enddo
            do i = n2/2+2,n2
              tmp2(i,j) = x1(n2-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scaley,tmp2,n2,nsizexy)

          do j = 1, nsizexy 
            do i = 1, n2
              x0(i,j,k) = tmp2(i,j) 
            enddo
          enddo
        enddo

        deallocate(x1)
        allocate(x1(n3/2+1,nsizexy,nsizeyz))
!        call MPI_BARRIER(commcol,ierr)
        call trans3d3_TRANSP(n2,nsizexy,nsizez,nsizeyz,x0,x1,npx,&
                      xstable,gxrtable,commrow,myidx,n3/2+1)
        deallocate(x0)

        do k = 1, nsizeyz
          do j = 1, nsizexy
            do i = 1, n3/2+1
              tmp3(i,j) = x1(i,j,k) 
            enddo
            do i = n3/2+2,n3
              tmp3(i,j) = x1(n3-i+2,j,k)
            enddo
          enddo

          call fftlocal_FFT(ksign,scalez,tmp3,n3,nsizexy)

          do j = 1, nsizexy
            do i = 1, n3
              grnout(i,j,k) = tmp3(i,j)
            enddo
          enddo
        enddo

        deallocate(x1)

        t_greenf = t_greenf + elapsedtime_Timer(t0)

        end subroutine greenf1t
!--------------------------------------------------------------------

        subroutine setval_FieldQuant(this,i,j,k,value)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this
        integer, intent(in) :: i, j, k
        double precision, intent(in) :: value

        this%FieldQ(i,j,k) = value

        end subroutine setval_FieldQuant

        double precision function get_FieldQuant(this,i,j,k)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(in) :: i, j, k

        get_FieldQuant = this%FieldQ(i,j,k)

        end function get_FieldQuant

        subroutine getglb_FieldQuant(this,temp)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        type (FieldQuant), intent(out) :: temp
        integer :: i, j, k, lcnz,lcny,lcnx
        double precision :: value

        lcnx = this%Nxlocal
        lcny = this%Nylocal
        lcnz = this%Nzlocal
    
        do k = 1, lcnz
          do j = 1, lcny
            do i = 1, lcnx
              value = get_FieldQuant(this,i,j,k)
              call setval_FieldQuant(temp,i,j,k,value)
            enddo
          enddo 
        enddo

        end subroutine getglb_FieldQuant

        subroutine getlcgrid_FieldQuant(this,nxlc,nylc,nzlc)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(in) :: this
        integer, intent(out) :: nxlc,nylc,nzlc

        nxlc = this%Nxlocal
        nylc = this%Nylocal
        nzlc = this%Nzlocal

        end subroutine getlcgrid_FieldQuant

        subroutine destruct_FieldQuant(this)
        implicit none
        include 'mpif.h'
        type (FieldQuant), intent(out) :: this

        deallocate(this%FieldQ) 

        end subroutine destruct_FieldQuant

      end module FieldQuantclass
