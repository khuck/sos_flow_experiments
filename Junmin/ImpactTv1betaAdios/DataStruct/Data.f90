!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! Dataclass: Field data class in DATA STRUCTURE layer.
! Version: 1.0-beta
! Author: Ji Qiang, LBNL
! Description: This class stores the rf cavity data Ez, Ez', Ez'' on the 
!              axis; Fourier coefficients of Ez on the axis; Ez(r,z),
!              Er(r,z), Htheta(r,z) on the r-z grid plane; and Ex(x,y,z),
!              Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z) on
!              uniform x, y, z grid, and Br(r,z) and Bz(r,z) on the r-z grid.
! Comments: 
!----------------------------------------------------------------
      module Dataclass
        use NumConstclass
!        use mpistub
!        save
!-----------------------------------------------------------------------
! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) directly.
        !number of grid points along x, y, and z direction.
        integer :: NxIntvRfg = 1
        integer :: NyIntvRfg = 1
        integer :: NzIntvRfg = 1
        !range in x, y, and zdirections.
        double precision :: XmaxRfg,XminRfg,YmaxRfg,YminRfg,ZmaxRfg,ZminRfg
        ! discrete Ex(x,y,z), Ey(x,y,z), Ez(x,y,z) and Bx(x,y,z), By(x,y,z), and 
        ! Bz(x,y,z) rf data. Here, the grid is uniform in x, y and z.
        double precision,allocatable,dimension(:,:,:) :: &
               Exgrid,Eygrid,Ezgrid,Bxgrid,Bygrid,Bzgrid
!-----------------------------------------------------------------------
! using the r-z field data (Er,Ez,Htheta) directly.
        !number of grid points along r direction.
        integer :: NrIntvRf = 1
        !number of grid points along z direction.
        integer :: NzIntvRf = 1
        !range in r and z directions.
        double precision :: RmaxRf,RminRf,ZmaxRf,ZminRf
        ! discrete Ez(r,z), Er(r,z) and Htheta(r,z) rf data. Here, the grid
        ! is uniform in r and z.
        double precision,allocatable,dimension(:,:) :: &
               ezdata,erdata,btdata
!-----------------------------------------------------------------------
! using only on-axis field data and its derivities.
        !initial number of grid points on the axis.
        integer, parameter :: Ndataini = 5000
        ! discrete Ez(0,z), Ez'(0,z), Ez''(0,z) rf data.
        double precision,dimension(Ndataini) :: zdat,edat,epdat,eppdat
!----------------------------------------------------------------------
! using the Fourier coefficients
        !number of Fourier expansion coefficients.
        integer, parameter :: NcoefF = 401
        double precision,dimension(NcoefF) :: Fcoef
        !Fcoef(1): constant
        !Fcoef(2k): cosine term
        !Fcoef(2k+1): sine term
!----------------------------------------------------------------------
        ! practical number of grid data on the axis or Fourier coefficients.
        integer :: Ndata
!------------------------------------------------------------------------
!
!       data storage for the t code
        type fielddata
          ! using the x-y-z field data (Ex,Ey,Ez,Bx,By,Bz) 
          !directly for beam line element i.
          !number of grid points along x, y, and z direction.
          integer :: NxIntvRfgt 
          integer :: NyIntvRfgt
          integer :: NzIntvRfgt
          !range in x, y, and zdirections.
          double precision  :: XmaxRfgt,XminRfgt,&
               YmaxRfgt,YminRfgt,ZmaxRfgt,ZminRfgt
          ! discrete Ex(x,y,z,i), Ey(x,y,z,i), Ez(x,y,z,i) and Bx(x,y,z,i), 
          ! By(x,y,z,i), and
          ! Bz(x,y,z,i) rf data. Here, the grid is uniform in x, y and z.
          double precision,pointer,dimension(:,:,:) :: &
               Exgridt,Eygridt,Ezgridt,Bxgridt,Bygridt,Bzgridt
          ! using the r-z field data (Er,Ez,Htheta) directly.
          !number of grid points along r direction.
          integer :: NrIntvRft
          !number of grid points along z direction.
          integer :: NzIntvRft 
          !range in r and z directions.
          double precision :: RmaxRft,RminRft,&
              ZmaxRft,ZminRft
          ! discrete Ezt(r,z,i), Ert(r,z,i) and 
          ! Btheta(r,z,i) rf data. Here, the grid
          ! is uniform in r and z.
          double precision,pointer,dimension(:,:) :: &
               ezdatat,erdatat,btdatat,brdatat,bzdatat
          ! using the analytical function for field data 
          ! (Ex,Ey,Ez,Bx,By,Bz) for beam line element i directly.
          !number of coefficients in Ex.
          integer :: Ncex 
          !number of coefficients in Ey.
          integer :: Ncey 
          !number of coefficients in Ez.
          integer :: Ncez 
          !number of coefficients in Bx.
          integer :: Ncbx
          !number of coefficients in By.
          integer :: Ncby
          !number of coefficients in Bz.
          integer :: Ncbz 
          !coefficients for analytical function description of the rf data.
          double precision,pointer,dimension(:) :: &
              coefex,coefey,coefez,coefbx,coefby,coefbz
          ! using the Fourier coefficients
          !integer, parameter :: NcoefFt = 401
          double precision,dimension(401) :: Fcoeft
          ! practical number of grid data on the axis or Fourier coefficients.
          integer :: Ndatat
        end type fielddata
      contains

!-------------------------------------------------------------------------
! the following are used to store the field data in position z domain, i.e.
! used by Impact-Z code.
        !Initialize the data storage arrays.
        subroutine init_Data()
        implicit none
        include 'mpif.h' 
        integer :: i,j

        NzIntvRf = 1
        NrIntvRf = 1
        allocate(ezdata(NzIntvRf+1,NrIntvRf+1))
        allocate(erdata(NzIntvRf+1,NrIntvRf+1))
        allocate(btdata(NzIntvRf+1,NrIntvRf+1))

        do j = 1, NrIntvRf+1
          do i = 1, NzIntvRf+1
            ezdata(i,j) = 0.0
            erdata(i,j) = 0.0
            btdata(i,j) = 0.0
          enddo
        enddo

        do i = 1, Ndataini
          zdat(i) = 0.0
          edat(i) = 0.0
          epdat(i) = 0.0
          eppdat(i) = 0.0
        enddo

        Ndata = 1
        RminRf = 0.0
        RmaxRf = 1.0
        ZminRf = 0.0
        ZmaxRf = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        NxIntvRfg = 1
        NyIntvRfg = 1
        NzIntvRfg = 1
        allocate(Exgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Eygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Ezgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bxgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bygrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        allocate(Bzgrid(NxIntvRfg+1,NyIntvRfg+1,NzIntvRfg+1))
        Exgrid = 0.0
        Eygrid = 0.0
        Ezgrid = 0.0
        Bxgrid = 0.0
        Bygrid = 0.0
        Bzgrid = 0.0
        XminRfg = 0.0
        XmaxRfg = 1.0
        YminRfg = 0.0
        YmaxRfg = 1.0
        ZminRfg = 0.0
        ZmaxRfg = 1.0

        end subroutine init_Data

        subroutine destruct_Data()
        implicit none
        include 'mpif.h' 

        deallocate(ezdata)
        deallocate(erdata)
        deallocate(btdata)
        deallocate(Exgrid)
        deallocate(Eygrid)
        deallocate(Ezgrid)
        deallocate(Bxgrid)
        deallocate(Bygrid)
        deallocate(Bzgrid)
         
        end subroutine destruct_Data

!-------------------------------------------------------------------------
! the following are used for storing the field data in time t domain 
! used by Impact-T code.
        subroutine initt_Data(this)
        implicit none
        include 'mpif.h' 
        type (fielddata), intent(inout) :: this
        integer :: i,j

!initialization of Er, Ez, Btheta
        this%NzIntvRft = 1
        this%NrIntvRft = 1
        allocate(this%ezdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%erdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%btdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%brdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        allocate(this%bzdatat(this%NzIntvRft+1,this%NrIntvRft+1))
        this%ezdatat = 0.0
        this%erdatat = 0.0
        this%btdatat = 0.0
        this%brdatat = 0.0
        this%bzdatat = 0.0

        this%Ndatat = 1
        this%RminRft = 0.0
        this%RmaxRft = 1.0
        this%ZminRft = 0.0
        this%ZmaxRft = 1.0

!initialization of Ex,Ey,Ez,Bx,By,Bz
        this%NxIntvRfgt = 1
        this%NyIntvRfgt = 1
        this%NzIntvRfgt = 1
        allocate(this%Exgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Eygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Ezgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bxgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bzgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        this%Exgridt = 0.0
        this%Eygridt = 0.0
        this%Ezgridt = 0.0
        this%Bxgridt = 0.0
        this%Bygridt = 0.0
        this%Bzgridt = 0.0
        this%XminRfgt = 0.0
        this%XmaxRfgt = 1.0
        this%YminRfgt = 0.0
        this%YmaxRfgt = 1.0
        this%ZminRfgt = 0.0
        this%ZmaxRfgt = 1.0
! initialization of the coefficients for analytical discription of Ex,Ey,Ez,Bx,By,Bz.
        this%Ncex = 1
        this%Ncey = 1
        this%Ncez = 1
        this%Ncbx = 1
        this%Ncby = 1
        this%Ncbz = 1
        allocate(this%coefex(this%Ncex))
        allocate(this%coefey(this%Ncey))
        allocate(this%coefez(this%Ncez))
        allocate(this%coefbx(this%Ncbx))
        allocate(this%coefby(this%Ncby))
        allocate(this%coefbz(this%Ncbz))
        this%coefex = 0.0
        this%coefey = 0.0
        this%coefez = 0.0
        this%coefbx = 0.0
        this%coefby = 0.0
        this%coefbz = 0.0

        end subroutine initt_Data

        subroutine destructt_Data(this)
        implicit none
        include 'mpif.h' 
        type (fielddata), intent(inout) :: this

        print*,"d1: "
        deallocate(this%ezdatat)
        deallocate(this%erdatat)
        deallocate(this%btdatat)
        deallocate(this%brdatat)
        deallocate(this%bzdatat)
        print*,"d2: "
        deallocate(this%Exgridt)
        print*,"d20: "
        deallocate(this%Eygridt)
        print*,"d21: "
        deallocate(this%Ezgridt)
        print*,"d22: "
        deallocate(this%Bxgridt)
        print*,"d23: "
        deallocate(this%Bygridt)
        print*,"d24: "
        deallocate(this%Bzgridt)
        print*,"d3: "
        deallocate(this%coefex)
        deallocate(this%coefey)
        deallocate(this%coefez)
        deallocate(this%coefbx)
        deallocate(this%coefby)
        deallocate(this%coefbz)
        print*,"d4: "
         
        end subroutine destructt_Data

      end module Dataclass
