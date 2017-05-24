!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! Depositorclass: Charge deposition class of APPLICATION 
!                 layer.
! Version: 1.0-beta
! Author: Ji Qiang, LBNL
! Description: This class deposit the particles onto computational mesh implementation.
! Comments:
!----------------------------------------------------------------
      module depositorclass
        use Pgrid2dclass
        use CompDomclass
        use BeamBunchclass
        use PhysConstclass
      contains
        !find the charge density of each bunch
        subroutine chgdens_Depositor(pbunch,rho,ptsgeom,pgrid,&
                   gammaz,flagbc,perd,zcent,flagpos)
           type (BeamBunch) :: pbunch
           double precision, dimension(:,:,:) :: rho
           type (CompDom), intent(in) :: ptsgeom
           type (Pgrid2d), intent(in) :: pgrid
           double precision, intent(in) :: gammaz, perd, zcent
           integer, intent(in) :: flagbc,flagpos
           integer, dimension(3) :: lcnum
           integer :: npx,npy,totnp,myid,myidy,myidx,nplc,npglb
           double precision, dimension(3) :: msize
           double precision :: hxi,hyi,hzi,tmp1,tmp2,pchg,pcurr

           nplc = pbunch%Nptlocal
           npglb = pbunch%Npt 
           call getsize_Pgrid2d(pgrid,totnp,npy,npx)
           call getpost_Pgrid2d(pgrid,myid,myidy,myidx)

           call getlcmnum_CompDom(ptsgeom,lcnum)
           innx = lcnum(1)
           if(npy.gt.1) then
             inny = lcnum(2) + 2
           else
             inny = lcnum(2)
           endif
           if(npx.gt.1) then
             innz = lcnum(3) + 2
           else
             innz = lcnum(3)
           endif

           if(flagbc.eq.1) then
             if(flagpos.eq.1) then
               call deposit1p_Depositor(nplc,innx,inny,innz,pbunch%Pts1, &
               rho,ptsgeom,npx,npy,myidx,myidy)
             else
               call deposit1_Depositor(nplc,innx,inny,innz,pbunch%Pts1, &
               rho,ptsgeom,npx,npy,myidx,myidy)
             endif
           endif

           pchg = pbunch%Charge
           pcurr = pbunch%Current
!           tmp1 = pcurr/Scfreq*pchg/abs(pchg)*float(nppos)/npglb
           tmp1 = pcurr/Scfreq*pchg/abs(pchg)/npglb
           call getmsize_CompDom(ptsgeom,msize)
           hxi = 1.0/msize(1)
           hyi = 1.0/msize(2)
           hzi = 1.0/msize(3)
           !gammaz is due to relativistic factor
!           tmp2 = tmp1*hxi*hyi*hzi/(nppos*Scxlt*Scxlt*Scxlt*Epsilon0)/gammaz
           tmp2 = tmp1*hxi*hyi*hzi/(Scxlt*Scxlt*Scxlt*Epsilon0)/gammaz
           !Here, rho is rho/epsilon0
           do k = 1, innz
             do j = 1, inny
               do i = 1, innx
                 rho(i,j,k) = rho(i,j,k)*tmp2
               enddo
             enddo
           enddo
 
           if((Flagbc.eq.1)) then  ! 3D open
             !gather contributions from guard cells.
             call guardsum1_Fldmger(rho,innx,inny,innz,pgrid)
           else
             print*,"no such type of boundary conditions!!!"
           endif
!           print*,"in depo sum(rho): ",sum(rho),tmp1,tmp2,nplc,&
!                  innx,inny,innz,myid

        end subroutine chgdens_Depositor

        ! deposit particles onto grid. (3D open boundary condition)
        subroutine deposit1_Depositor(innp,innx,inny,innz,rays,rho,& 
                                ptsgeom,npx,npy,myidx,myidy)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (6, innp) :: rays
        double precision, intent (out), dimension (innx,inny,innz) :: rho
!        logical, intent (in), dimension (innp) :: msk
        integer :: ix,jx,kx,ix1,jx1,kx1
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        integer :: ngood, i, j, k, kadd, jadd
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        type (CompDom) :: ptsgeom

        call starttime_Timer( t0 )

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        rho=0.
        do i = 1, innp
          ix=(rays(1,i)-xmin)*hxi + 1
          ab=((xmin-rays(1,i))+ix*hx)*hxi
          jx=(rays(3,i)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,i))+(jx-jadd)*hy)*hyi
          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          ix1=ix+1
          if( (ix1.gt.innx).or.(ix.lt.1)) then
            print*,"dep - ix: ",ix,rays(1,i),xmin
            stop
          endif
          jx1=jx+1
          if( (jx1.gt.inny).or.(jx.lt.1)) then
            print*,"dep - jx1: ",jx,rays(3,1),ymin,inny,myidx,myidy
            stop
          endif
          kx1=kx+1
          if( (kx1.gt.innz).or.(kx.lt.1)) then
            print*,"dep - kx: ",kx,rays(5,i),zmin
            stop
          endif
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef
        enddo

        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine deposit1_Depositor

        ! deposit particles onto grid. (3D open boundary condition)
        subroutine deposit1p_Depositor(innp,innx,inny,innz,rays,rho,& 
                                ptsgeom,npx,npy,myidx,myidy)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: innp,innx,inny,innz,npx,npy,myidx,myidy
        double precision, intent (in), dimension (6, innp) :: rays
        double precision, intent (out), dimension (innx,inny,innz) :: rho
!        logical, intent (in), dimension (innp) :: msk
        integer :: ix,jx,kx,ix1,jx1,kx1
        integer, dimension(2,0:npx-1,0:npy-1) :: table
        integer, dimension(0:npx-1) :: xtable
        integer, dimension(0:npy-1) :: ytable
        double precision :: ab,cd,ef
        integer :: ngood, i, j, k, kadd, jadd
        double precision :: t0
        double precision :: hx,hxi,hy,hyi,hz,hzi,xmin,ymin,zmin
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        type (CompDom) :: ptsgeom

        call starttime_Timer( t0 )

        call getmsize_CompDom(ptsgeom,msize)
        hxi = 1.0/msize(1)
        hyi = 1.0/msize(2)
        hzi = 1.0/msize(3)
        hx = msize(1)
        hy = msize(2)
        hz = msize(3)

        call getrange_CompDom(ptsgeom,range)
        xmin = range(1)
        ymin = range(3)
        zmin = range(5)

        call getlctabnm_CompDom(ptsgeom,table)
        xtable(0) = 0
        do i = 1, npx - 1
          xtable(i) = xtable(i-1) + table(1,i-1,0)
        enddo

        ytable(0) = 0
        do i = 1, npy - 1
          ytable(i) = ytable(i-1) + table(2,0,i-1)
        enddo

        if(npx.gt.1) then
          kadd = -xtable(myidx) + 1
        else
          kadd = 0
        endif
        if(npy.gt.1) then
          jadd = -ytable(myidy) + 1
        else
          jadd = 0
        endif

        rho=0.
        do i = 1, innp
          if(rays(5,i).gt.0.0) then !//calculate the space-charge force for particle z > 0
          ix=(rays(1,i)-xmin)*hxi + 1
          ab=((xmin-rays(1,i))+ix*hx)*hxi
          jx=(rays(3,i)-ymin)*hyi + 1 + jadd
          cd=((ymin-rays(3,i))+(jx-jadd)*hy)*hyi
          kx=(rays(5,i)-zmin)*hzi + 1 + kadd
          ef=((zmin-rays(5,i))+(kx-kadd)*hz)*hzi
          ix1=ix+1
          if( (ix1.gt.innx).or.(ix.lt.1)) then
            print*,"dep - ix: ",ix,rays(1,i),xmin
            stop
          endif
          jx1=jx+1
          if( (jx1.gt.inny).or.(jx.lt.1)) then
            print*,"dep - jxp: ",jx,rays(3,1),ymin,inny,myidx,myidy
            stop
          endif
          kx1=kx+1
          if( (kx1.gt.innz).or.(kx.lt.1)) then
            print*,"dep - kx: ",kx,rays(5,i),zmin
            stop
          endif
          ! (i,j,k):
          rho(ix,jx,kx) = rho(ix,jx,kx) + ab*cd*ef
          ! (i,j+1,k):
          rho(ix,jx1,kx) = rho(ix,jx1,kx) + ab*(1.0-cd)*ef
          ! (i,j+1,k+1):
          rho(ix,jx1,kx1) = rho(ix,jx1,kx1)+ab*(1.0-cd)*(1.0-ef)
          ! (i,j,k+1):
          rho(ix,jx,kx1) = rho(ix,jx,kx1)+ab*cd*(1.0-ef)
          ! (i+1,j,k+1):
          rho(ix1,jx,kx1) = rho(ix1,jx,kx1)+(1.0-ab)*cd*(1.0-ef)
          ! (i+1,j+1,k+1):
          rho(ix1,jx1,kx1) = rho(ix1,jx1,kx1)+(1.0-ab)*(1.0-cd)*(1.0-ef)
          ! (i+1,j+1,k):
          rho(ix1,jx1,kx) = rho(ix1,jx1,kx)+(1.0-ab)*(1.0-cd)*ef
          ! (i+1,j,k):
          rho(ix1,jx,kx) = rho(ix1,jx,kx)+(1.0-ab)*cd*ef

          endif
        enddo

        t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

        end subroutine deposit1p_Depositor

      end module depositorclass
