!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! AccSimulatorclass: Linear accelerator simulator class in CONTROL layer.
! Version: 1.0-beta
! Author: Ji Qiang
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
        use Pgrid2dclass
        use CompDomclass
        use FieldQuantclass
        use BeamLineElemclass
        use Ptclmgerclass
        use BeamBunchclass
        use Timerclass
        use Inputclass
        use Outputclass
        use Dataclass
        use PhysConstclass
        use NumConstclass
        use Distributionclass
        use Rangerclass
        use Depositorclass
        implicit none
        !# of phase dim., num. total and local particles, int. dist. 
        !and restart switch, error study switch, substep for space-charge
        !switch,# of time step
        integer :: Dim, Flagdist,Rstartflg,Flagerr,&
                            Flagsubstep,ntstep 
        integer, dimension(Nbunchmax) :: Np, Nplocal
        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        !FlagImage: switch flag for image space-charge force calculation: "1" for yes, 
        !otherwise for no. 
        integer :: Nx,Ny,Nz,Nxlocal,Nylocal,Nzlocal,Flagbc,&
                            Nblem,Flagmap,Flagdiag,FlagImage

        !# of processors in column and row direction.
        integer :: npcol, nprow

        !initial # of bunches/bins
        integer :: Nbunch

        !beam current, kin. energy, part. mass, charge, ref. freq., period length, 
        !time step size 
        double precision :: Bcurr,Bkenergy,Bmass,Bcharge,Bfreq,&
                                     Perdlen,dt,xrad,yrad

        !conts. in init. dist.
        integer, parameter :: Ndistparam = 21
        double precision, dimension(Ndistparam) :: distparam

        !2d logical processor array.
        type (Pgrid2d) :: grid2d

        !beam particle object and array.
        type (BeamBunch), dimension(Nbunchmax) :: Ebunch

        !beam charge density and field potential arrays.
        type (FieldQuant) :: Potential

        !geometry object.
        type (CompDom) :: Ageom

        !overlaped external field data array
        type (fielddata), dimension(Maxoverlap) :: fldmp

        !maximum e- emission time
        double precision :: temission
        !number of steps for emission
        integer :: Nemission

        !distance after that to turn off image space-charge
        double precision :: zimage

        !restart time and step
        double precision :: tend 
        integer :: iend,nfileout,ioutend,itszend,isteerend

        !beam line element array.
        type (BPM),target,dimension(Nbpmmax) :: beamln0
        type (DriftTube),target,dimension(Ndriftmax) :: beamln1
        type (BeamLineElem),dimension(Nblemtmax)::Blnelem
        !longitudinal position of each element (min and max).
        double precision, dimension(2,Nblemtmax)::zBlnelem
        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface
      contains
        !set up objects and parameters.
        subroutine init_AccSimulator(time)
        implicit none
        include 'mpif.h'
        integer :: i,test1,test2,j
        integer :: myid,myidx,myidy,ierr,inb,jstp
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,val0,&
        val2, val3,val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24
        double precision :: time
        double precision :: t0,zz
        double precision :: z,phsini
        double precision, dimension(2) :: tmpdr 
        double precision, dimension(5) :: tmpcf 
        double precision, dimension(8) :: tmpbpm 
        double precision, dimension(9) :: tmpquad
        double precision, dimension(10) :: tmpdipole 
        double precision, dimension(11) :: tmprf
        double precision, dimension(12) :: tmpslrf
        double precision, dimension(13) :: tmp13
        double precision, dimension(25) :: tmpdtl
        integer :: iqr,idr,ibpm,iccl,iccdtl,idtl,isc,icf,islrf,isl,idipole,&
          iemfld,iemfldcart,iemfldcyl,iemfldana,ib,imult
        character(20) :: bpname3 = 'hi.bp'


        !start up MPI.
        call init_Input(time)

        ! initialize Timer.
        call construct_Timer(0.0d0)

        call starttime_Timer(t0)

        !Flagmap = 0

!-------------------------------------------------------------------
! get all global input parameters.
        call in_Input(Dim,Np(1),Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,Ndistparam,Bcurr,Bkenergy,Bmass,Bcharge,&
        Bfreq,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
        Flagsubstep,phsini,dt,ntstep,Nbunch,FlagImage,Nemission,&
        temission,zimage)

!-------------------------------------------------------------------
! construct 2D logical processor Cartesian coordinate
        call construct_Pgrid2d(grid2d,MPI_COMM_WORLD,nprow,npcol)
        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        if(myid.eq.0) then
          print*,"!-----------------------------------------------------------"
          print*,"! IMPACT-T Parallel Beam Dynamics Tracking Code: version 1.0beta"
          print*,"! Copyright of The Regents of the University of California"
          print*,"!-----------------------------------------------------------"
          print*,"Bcurr: ",Bcurr,Bkenergy,Bmass,Bcharge,Bfreq
        endif

        !construct Constants class.
        call constructT_PhysConst(dt,Bfreq)

!-------------------------------------------------------------------
! construct computational domain CompDom class and get local geometry 
! information on each processor.
        !if(Rstartflg.eq.1) then
        !  call ingeom_Output(1500,z,inb,jstp,nprow,npcol,Ageom,Nx,Ny,Nz,&
        !                    myidx,myidy)
        !  if(myid.eq.0) print*,"rstart at: ",z,inb,jstp
        !  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !else
          !xrad = 0.1363243029*0.2
          call construct_CompDom(Ageom,distparam,Ndistparam,Flagdist,&
               Nx,Ny,Nz,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
        !endif

!-------------------------------------------------------------------
! initialize Data class.
        do i = 1, Maxoverlap
          call initt_Data(fldmp(i))
        enddo

!-------------------------------------------------------------------
! construct BeamBunch class.
        call construct_BeamBunch(Ebunch(1),Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np(1),phsini)

        !initial value for the output file number
        nfileout = 40
!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val0(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))
        allocate(val9(Nblem),val10(Nblem),val11(Nblem),val12(Nblem))
        allocate(val13(Nblem),val14(Nblem),val15(Nblem),val16(Nblem))
        allocate(val17(Nblem),val18(Nblem),val19(Nblem),val20(Nblem))
        allocate(val21(Nblem),val22(Nblem),val23(Nblem),val24(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val0,val1,&
        val2,val3,&
        val4,val5,val6,val7,val8,val9,val10,val11,val12,val13,val14,&
        val15,val16,val17,val18,val19,val20,val21,val22,val23,val24)

        idr = 0
        ibpm = 0
        zz = 0.0
        !If we allow the starting edge "zz" of one beam line element
        !to be inside the preceding beam, this allows the beam line
        !elment overlaping. 
        do i = 1, Nblem
          zBlnelem(1,i) = val0(i)
          if(bitype(i).lt.0) then
            ibpm = ibpm + 1
            call construct_BPM(beamln0(ibpm),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpbpm(1) = val0(i)
            tmpbpm(2) = val1(i)
            tmpbpm(3) = val2(i)
            tmpbpm(4) = val3(i)
            tmpbpm(5) = val4(i)
            tmpbpm(6) = val5(i)
            tmpbpm(7) = val6(i)
            tmpbpm(8) = val7(i)
            call setparam_BPM(beamln0(ibpm),tmpbpm)
            Blnelem(i) = assign_BeamLineElem(beamln0(ibpm))
            !reset the output file number from BPM type "-3".
            if(bitype(i).eq.(-3)) then
              nfileout = bmpstp(i)
            endif
          else if(bitype(i).eq.0) then
            idr = idr + 1
            call construct_DriftTube(beamln1(idr),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpdr(1) = val0(i)
            tmpdr(2) = val1(i)
            call setparam_DriftTube(beamln1(idr),tmpdr)
            Blnelem(i) = assign_BeamLineElem(beamln1(idr))
          else 
            print*,"element type not available!"
          endif 
          zz = val0(i) + blength(i)
          zBlnelem(2,i) = zz
        enddo
!-------------------------------------------------------------------
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass setting up lattice..."

        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val0)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8,val9)
        deallocate(val10,val11,val12,val13,val14,val15,val16)
        deallocate(val17,val18,val19,val20,val21,val22,val23)
        deallocate(val24)

!-------------------------------------------------------------------
! sample initial particle distribution.
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        tend = 0.0d0
        iend = 0
        ioutend = 0
        itszend = 0
        isteerend = 0
        if(Rstartflg.eq.1) then
!          call inpoint_Output(nfileout+myid,Ebunch(1),tend,iend,ib,nprow,npcol,&
!               Ageom,Nx,Ny,Nz,myidx,myidy,Np(1),ioutend,itszend,isteerend)

           call phase_in(Ebunch(1),bpname3,Nplocal(1),Np(1))
       
          if(myid.eq.0) print*,"restart at: ",tend,iend,ib
        else
          !call sample_Dist(Ebunch(1),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc)
          ib = 1
          call sample_Dist(Ebunch(1),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc,ib,Nbunch)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid.eq.0) print*,"pass generating initial distribution..."

!-------------------------------------------------------------------
! construct FieldQuant class objects.
        call construct_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d)

!        print*,"start to read in the inputs for the other bunches/bins:.."
        !generate the particle distribution for the other bunches/bins
        do ib = 2, Nbunch
          call in_Input(Dim,Np(ib),Nx,Ny,Nz,Flagbc,Flagdist,Rstartflg,&
              Flagmap,distparam,Ndistparam,Bcurr,Bkenergy,Bmass,Bcharge,&
          Bfreq,xrad,yrad,Perdlen,Nblem,npcol,nprow,Flagerr,Flagdiag,&
          Flagsubstep,phsini,dt,ntstep,Nbunch,FlagImage,ib)
          call construct_BeamBunch(Ebunch(ib),Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np(ib),phsini)
          if(Rstartflg.eq.1) then
            call inpoint_Output(ib*nfileout+myid,Ebunch(ib),tend,iend,ib,nprow,npcol,&
                Ageom,Nx,Ny,Nz,myidx,myidy,Np(ib),ioutend,itszend,isteerend)
          else
            !call sample_Dist(Ebunch(ib),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc)
            call sample_Dist(Ebunch(ib),distparam,Ndistparam,Flagdist,Ageom,grid2d,Flagbc,ib,Nbunch)
          endif
        enddo

        !get local particle number and mesh number on each processor.
        do ib = 1, Nbunch
          call getnpt_BeamBunch(Ebunch(ib),Nplocal(ib))
        enddo

        t_init = t_init + elapsedtime_Timer(t0)
        !this is just for finding the driven phase of each cavity
        !tend = phsini

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        implicit none
        include 'mpif.h'
        integer :: i,j,bnseg,bmpstp,bitype
        integer :: myid,myidx,myidy,comm2d,commcol,commrow,npx,npy,&
                   totnp,ierr,nylcr,nzlcr
        integer :: nbal,ibal,istep,ifile
        integer :: ibalend,istepend,nfile
        integer, dimension(3) :: lcgrid
        double precision :: z0,z,tau1,tau2,blength,t0,zend
        double precision, dimension(6) :: lcrange,range,ptrange,sgcenter,grange
        double precision, dimension(3) :: msize
        double precision :: hy,hz,ymin,zmin,piperad,zedge,hzwake
        double precision :: tmp1,tmp2,tmp3,tmp4,rfile
        double precision, allocatable, dimension(:,:,:) :: chgdens,tmppot
        double precision, allocatable, dimension(:,:,:) :: exg,eyg,ezg
        double precision, allocatable, dimension(:,:,:) :: bxg,byg,bzg
        double precision, allocatable, dimension(:,:,:) :: besscoef
        double precision, allocatable, dimension(:,:) :: bessnorm,gml
        integer, allocatable, dimension(:) :: modth,pydisp
        integer :: nmod,k
        !double precision :: sumtest, sumtest2, sumtest3
        double precision, dimension(8) :: drange
        double precision, dimension(3) :: al0,ga0,epson0
        double precision :: realSamplePeriod
        integer :: nsubstep,integerSamplePeriod
        double precision :: zcent,distance,blnLength,dzz
        integer, allocatable, dimension(:,:) :: idrfile
        integer :: ibend,ibstart,isw,ibinit,ibendold,iifile,ii,ibinitold,idbeamln
        double precision :: zmax,t,dtless,gammazavg,curr
        integer :: tmpflag,ib,ibb,ibunch,inib,nplctmp,nptmp,nptottmp
        double precision, allocatable, dimension(:) :: gammaz
        double precision, allocatable, dimension(:,:) :: brange
        double precision :: dGspread
        integer, dimension(Maxoverlap) :: tmpfile
        double precision :: tmpcur,totchrg,r0
        integer :: flagpt2pt,flagpos, flagcathode !//switch for point-to-point calculation and switch for space-charge with z>0
        double precision :: zz,zorgin,zorgin2,vref,gamin,gam
        double precision, dimension(6) :: ptref
        integer :: idbd,idbend,flagbctmp
        integer :: npttmplc,npttmp
        double precision :: ztmp1,ztmp2,deltaz
        integer :: ipt,iptnew
        double precision :: betazini
        !//switch for output the phase space plot: 0 - no output, 1 - output 
        !//6D particle phase space at given time, 2 - output for restart 
        integer :: flagphout,iout
        !//switch for changing time step size: 0 - no change, 1 - change
        !// time step size after tszstart
        integer :: flagtimesz
        double precision :: trstart
        integer, dimension(100) :: nsamp,nfileouttmp
        double precision, dimension(101) :: tszstart,dtnewsize
        double precision, dimension(101) :: tphout,tsteer,xoffset,yoffset,&
                                            pxoffset,pyoffset

        integer :: flagwake,kz,jy,ix,jadd,kadd,kst,jst,nz2,isteer,iizz
        double precision, allocatable, dimension(:,:) :: sendensz,&
                                                         recvdensz
        double precision :: xx,yy,t3dstart,rr,tmger,tmpwk,tstop
        double precision, allocatable, dimension(:,:)  :: tmppts
        integer :: iwk,nwk,flagbtw
        integer  :: nplocal0,np0
        integer :: itsz
        real*8 :: bendlen,zz1,zz2,zwkmin,poscent,zfmin,zfmax,coeftol
        integer :: ldsg,nlsg,nrsg,npsg,msg,ncoefreal,iz

        integer :: flagstep,flagspc,itspc
        double precision, dimension(101) :: tspcstart,vspc
        integer :: my_rank

        double precision :: sec0, sec1

        character(len=8) :: timestepStrfmt = '(I5.5)'
        character(20) :: timestepStr;
        character(20) :: bpname1;
        character(20) :: bpname2 = 'ho.bp'
        character(20) :: bpname3 = 'hox.bp'

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        if(my_rank.eq.0) print *,' Beginning the simulator run!!!!!!!!!!!'

        flagspc = 1
        if(Nemission.gt.0) then
          flagcathode = 1
        else !no cathode model
          flagcathode = 0
          FlagImage = 0
        endif

!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        ibalend = 1 
        istepend = 0
        zend = 0.0
        r0 = 0.0/Scxlt
        flagpt2pt = 0 !"0" not point-2-point, "1" use point-2-point

        call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
        call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)
        call getsize_Pgrid2d(grid2d,totnp,npy,npx)

        call MPI_BARRIER(comm2d,ierr)
        call starttime_Timer(t0)


        nbal = 100
        ibal = ibalend
        istep = istepend
        z = zend

        !assign initial storage for charge density,Ex,Ey,Ez,Bx,By,Bz,image
        !charge potential. 
        allocate(chgdens(1,1,1))
        allocate(exg(1,1,1))
        allocate(eyg(1,1,1))
        allocate(ezg(1,1,1))
        allocate(bxg(1,1,1))
        allocate(byg(1,1,1))
        allocate(bzg(1,1,1))
        allocate(tmppot(1,1,1))

        !no phase space print out
        flagphout = 0
        tphout = 1.0e12 
        !no change of the time step size
        flagtimesz = 0
        tszstart = 1.0e12
        tspcstart = 1.0e12
        dtnewsize = 1.0e12
        trstart = 1.0e12
        t3dstart = 1.0e12
        tsteer = 1.0e12 
        tmger = 1.0e12
        tstop = 1.0e12
       
        iout = 0
        itsz = 0
        itspc = 0
        nwk = 0
        isteer = 0
        !idrfile is used to store the <element type>, <external data file name>,
        !and <id> for the internal data storage of each beamline element  
        allocate(idrfile(3,Nblem))
        idrfile = 1
        do i = 1, Nblem
          call getparam_BeamLineElem(Blnelem(i),blength,bnseg,bmpstp,&
                                     bitype)
          idrfile(1,i) = bitype
          !get external file id for each rf beam line element.
          if(bitype.gt.100) then
            call getparam_BeamLineElem(Blnelem(i),5,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          !get external file for solenoid.
          else if(bitype.eq.3) then
            call getparam_BeamLineElem(Blnelem(i),3,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          else if(bitype.eq.4) then
            call getparam_BeamLineElem(Blnelem(i),4,rfile)
            idrfile(2,i) = int(rfile + 0.1)
          endif
          if(bitype.eq.(-1)) then
            isteer = isteer + 1
            if(isteer.gt.100) then
              print*,"The maximum steering location is 100!!!"
              isteer = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tsteer(isteer))
            call getparam_BeamLineElem(Blnelem(i),4,xoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),5,pxoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),6,yoffset(isteer))
            call getparam_BeamLineElem(Blnelem(i),7,pyoffset(isteer))
            xoffset(isteer) = xoffset(isteer)/Scxlt
            yoffset(isteer) = yoffset(isteer)/Scxlt
          endif
          if(bitype.eq.(-2)) then
            flagphout = 1
            iout = iout + 1
            if(iout.gt.100) then
              print*,"The maximum phase space output is 100!!!"
              iout = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tphout(iout))
            nfileouttmp(iout) = bmpstp 
            nsamp(iout) = bnseg
          endif
          if(bitype.eq.(-3)) then
            flagphout = 2
            call getparam_BeamLineElem(Blnelem(i),3,trstart)
            nfileout = bmpstp 
            nsamp = bnseg
          endif
          if(bitype.eq.(-4)) then
            flagtimesz = 1
            itsz = itsz + 1
            if(itsz.gt.100) then
              print*,"The maximum time step size change is 100!!!"
              itsz = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),3,tszstart(itsz))
            call getparam_BeamLineElem(Blnelem(i),4,dtnewsize(itsz))
          endif
          if(bitype.eq.(-5)) then
            call getparam_BeamLineElem(Blnelem(i),3,t3dstart)
          endif
          if(bitype.eq.(-7)) then
            call getparam_BeamLineElem(Blnelem(i),3,tmger)
          endif
          if(bitype.eq.(-8)) then
            itspc = itspc + 1
            if(itspc.gt.100) then
              print*,"The maximum space-charge change is 100!!!"
              itspc = 100
            endif
            call getparam_BeamLineElem(Blnelem(i),2,vspc(itspc))
            call getparam_BeamLineElem(Blnelem(i),3,tspcstart(itspc))
          endif

          if(bitype.eq.(-99)) then
            call getparam_BeamLineElem(Blnelem(i),3,tstop)
          endif
        enddo

        dtless = 1.0 !dimensionless time step size.
        t = tend
        distance = 0.0
        tmpfile = 0
        !length of the total beamline
        blnLength = zBlnelem(2,Nblem)
        ibinitold = 1
        ibendold = 0
        iifile = 0

        allocate(gammaz(Nbunch))
        allocate(brange(12,Nbunch))
        !count the total current and # of particles and local # of particles for each 
        !bunch or bin
        curr = 0.0
        do ib = 1, Nbunch
          curr = curr + Ebunch(ib)%Current
          Nplocal(ib) = Ebunch(ib)%Nptlocal
          Np(ib) = Ebunch(ib)%Npt
        enddo
        nplocal0 = Nplocal(1)
        np0 = Np(1)
        !output the moments from the average of all bunches at fixed z.
        call diagnostic1avg_Output(t,Ebunch,Nbunch)
        itspc = 0
        iout = ioutend
        itsz = itszend
        isteer = isteerend
        ibunch = Nbunch
        tmpcur = curr
        flagpos = 1
        flagstep = 1
        zorgin = 0.0
        idbd = 1
        !particles behind the cathode will use this beta for emission.
        !betazini = sqrt(1.0-1.0/(1.0+distparam(21)**2))
        betazini = sqrt(1.0-1.0/(1.0+Bkenergy/Bmass)**2)
        !Ebunch(1)%Pts1 = 0.0d0

        !output initial phase space distribution 

        do ib = 1, Nbunch
          !call phaseold_Output(40+ib-1,Ebunch(ib),1)
           !!write (timestepStr,timeStepStrFmt) ib
           !!bpname1='start'//trim(timestepStr)//'.bp'
           !!write(*,*) bpname1, ib, ntstep
           !!call phase_Output(bpname1, 40+ib-1,Ebunch(ib),1, 0, ntstep)
           !!//call phase_Output(bpname1, Nbunch-i,Ebunch(ib),1, 0, ntstep)
        enddo
        
        !!ib=1
        write (timestepStr,timeStepStrFmt) 1
        bpname1='start'//trim(timestepStr)//'.bp'

        !!call phase_Output(bpname1, 40, Ebunch(1), 1, 0, ntstep)
        call phaseBunchOutput(bpname1, Nbunch, Ebunch,1, 0, ntstep)

        dzz = betazini*Clight*dtless*Dt
        zmin = 0.0
        call MPI_BARRIER(comm2d,ierr)
        !print*,"iout: ",iout,tphout(1)
!----------------------------------------------------------------------
! start looping through ntstep timestep.
        !iend is the time step number from last simulation (used in
        !restart function).

        do i = iend+1, ntstep
           !write (timestepStr,timeStepStrFmt) i
           !bpname1='out'//trim(timestepStr)//'.bp'


           !!call phase_Output(bpname1, 40,Ebunch(1),1, i, ntstep)
           call phaseBunchOutput(bpname1, Nbunch, Ebunch,1, i, ntstep)
           
          !steering the beam centroid to the given X, PX,Y, PY values at given time
          !if(t.le.tsteer(isteer+1) .and. (t+dtless*Dt).ge.tsteer(isteer+1)) then
          if(distance.le.tsteer(isteer+1) .and. (distance+dzz).ge.tsteer(isteer+1)) then
            isteer = isteer + 1
            do ib = 1, Nbunch
              !//find the range and center information of each bunch/bin
              call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
              do ipt = 1, Nplocal(ib)
                 Ebunch(ib)%Pts1(1,ipt) = Ebunch(ib)%Pts1(1,ipt) - sgcenter(1) + xoffset(isteer)
                 Ebunch(ib)%Pts1(2,ipt) = Ebunch(ib)%Pts1(2,ipt) - sgcenter(2) + pxoffset(isteer)
                 Ebunch(ib)%Pts1(3,ipt) = Ebunch(ib)%Pts1(3,ipt) - sgcenter(3) + yoffset(isteer)
                 Ebunch(ib)%Pts1(4,ipt) = Ebunch(ib)%Pts1(4,ipt) - sgcenter(4) + pyoffset(isteer)
              enddo
            enddo
          endif

          if(i.le.Nemission) then !first Nemission steps for emission
            dtless = temission/Nemission/Dt
          else
            if(flagstep .eq. 1) then
              dtless = 1.0
            endif
          endif

          !change time step size at tszstart
          !if(t.le.tszstart(itsz+1) .and. (t+dtless*Dt).ge.tszstart(itsz+1)) then
          if(distance.le.tszstart(itsz+1) .and. (distance+dzz).ge.tszstart(itsz+1)) then
            itsz = itsz + 1
            dtless = dtnewsize(itsz)/Dt
            flagstep = 0
          endif

          !change charge space at tspcstart
          !if(t.le.tspcstart(itspc+1) .and. (t+dtless*Dt).ge.tspcstart(itspc+1)) then
          if(distance.le.tspcstart(itspc+1) .and. (distance+dzz).ge.tspcstart(itspc+1)) then
            itspc = itspc + 1
            if(vspc(itspc).gt.0.0) then
              flagspc = 1
            else
              flagspc = -1
            endif
          endif

          !//update particle positions using velocity for half step.
          if(flagcathode.eq.1) then
            do ib = 1, Nbunch   
              call drifthalf_BeamBunch(Ebunch(ib),t,dtless,betazini)
            enddo   
          else
            do ib = 1, Nbunch   
              call drifthalforg_BeamBunch(Ebunch(ib),t,dtless)
            enddo   
          endif
          t = t + 0.5*dtless*Dt

          !ibunch is total number of bunches within the effective computational domain
          !only bunch/bin with zmax>0 is counted as effective bunch/bin 
          if(ibunch.lt.Nbunch) then
            !ibunch = 0
            do ib = 1, Nbunch
              !//find the range and center information of each bunch/bin
              call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
              if(ptrange(6).gt.0.0) then
                !ibunch = ibunch + 1
              endif 
            enddo
          endif

          !only the bunch with zmax>0 is counted as an effective bunch 
          do ib = 1, ibunch
            !//find the range and center of each bunch/bin
            call singlerange(Ebunch(ib)%Pts1,Nplocal(ib),Np(ib),&
                             ptrange,sgcenter)
            !Ebunch(ib)%refptcl(5) = sgcenter(5) + zorgin
            Ebunch(ib)%refptcl(5) = sgcenter(5) 
            gammaz(ib) = sqrt(1.0+sgcenter(6)**2)!//gammaz from <gamma_i betaz_i>
            Ebunch(ib)%refptcl(6) = -gammaz(ib)
            do inib = 1, 6
              brange(inib,ib) = ptrange(inib)
              brange(inib+6,ib) = sgcenter(inib)
            enddo
            !the longitudinal range for space-charge calculation has to be > 0
!            if(flagpos.eq.1) then
              if(ptrange(5).gt.0) then
                brange(5,ib) = ptrange(5)
              else
                brange(5,ib) = 0.0
              endif
!            else
!              brange(5,ib) = ptrange(5)
!            endif
          enddo

          !here, ptrange(5) is the zminlc of the last bunch
          !this is the criterion for all bunches out of the cathode
          if(ptrange(5).gt.0.0  .or. (flagcathode.eq.0) ) flagpos = 0
          !get the global computational domain and center for all effective bunches/bins 
!          print*,"brange: ",brange,Np
          call globalrange(brange,grange,gammazavg,zcent,Np,ibunch) 
!          print*,"zcent: ",zcent

!          print*,"gamz: ",gammaz(1)
          !get the distance of the center of all effective bunches/bins
          distance = zcent*Scxlt
          dzz = sqrt(1.0d0-1.0d0/gammazavg**2)*Clight*dtless*Dt
          !exit if the beam is outside the beamline
          if(distance.gt.blnLength .or. distance.gt.tstop) then
            exit
          endif

          !check the particles outside the computational domain
          do ib = 1, ibunch
            call lost_BeamBunch(Ebunch(ib),xrad,yrad,Perdlen,zcent,&
                                nplctmp,nptmp)
            Nplocal(ib) = nplctmp
            Np(ib) = nptmp
!!            print*,"npt: ",ib,myid,nplctmp,nptmp
          enddo 

          !find the beginning and end beam line elements that the effective
          !bunch particles occupy
          !zmin = grange(5)*Scxlt + zorgin
          !zmax = grange(6)*Scxlt + zorgin
          zmin = grange(5)*Scxlt 
          zmax = grange(6)*Scxlt

          !if(mod(i,5).eq.0) then 
          !   if(myid.eq.0) print*,"zmin: ",zmin,zmax
          !endif

          !using the following way, we can find the id of the min element
          !and the max. element between which the bunch stays. It works even
          !if there are more than one element overlaps each other. However,
          !the beam line element has to be aranged so that zmin(i+1)>=zmin(i). 
          idbeamln = 1
          do ii = 1, Nblem
            if( (zmin.ge.zBlnelem(1,ii)) .and. &
                (zmin.le.zBlnelem(2,ii)) ) then
               idbeamln = ii
               exit !exit the loop from the first element id containing zmin
            endif
          enddo
          ibinit = idbeamln
          !exit the loop from the last element id containing zmax
          do ii = ibinit, Nblem
            if( (zmax.ge.zBlnelem(1,ii)) .and. &
                (zmax.le.zBlnelem(2,ii)) ) then
               idbeamln = ii 
            endif
          enddo
          ibend = idbeamln
          
          !In the following, we will find the id in the internal global data array
          !and assign it to each beam line element.
          !This information is used to get the external field at given position.
          !This should also work for the overlaped beam line element.
!          if(i.eq.1) then
!           ibstart = max(ibendold,ibinit) 
!         else
!           ibstart = max(ibendold,ibinit) + 1
!         endif
          if(ibendold.ge.ibinit) then
            ibstart = ibendold + 1
          else
            ibstart = ibinit
          endif

!          print*,"ibinit: ",ibinit,ibend,ibstart,ibend,ibendold
          idbend = 0
          do ii = ibstart,ibend
            !for element type > 100 or solenoid
            if((idrfile(1,ii).gt.100).or.(idrfile(1,ii).eq.3).or.&
               (idrfile(1,ii).eq.4)) then 
              isw = 1
              !check whether the new element is in the old elements range,
              !which already been read in. 
              do j = ibinitold,ibendold
                if(idrfile(2,ii).ne.idrfile(2,j)) then
                else
                  isw = 0
                  exit
                endif
              enddo
              !check whether the new element file has been read in by the
              !old element.
              do j = 1, mod(iifile-1,Maxoverlap)+1
                if(idrfile(2,ii).ne.tmpfile(j)) then
                else
                  isw = 0
                  idrfile(3,ii) = j
                  exit
                endif
              enddo
              !accumulate the external file 1d in the global data array
!              if(isw.eq.1) iifile = iifile + 1
!              !assign the id in global data array to each beam line element
!              idrfile(3,ii) = mod(iifile-1,Maxoverlap)+1
              if(isw.eq.1) then
                iifile = iifile + 1
                !assign the id in global data array to each beam line element
                idrfile(3,ii) = mod(iifile-1,Maxoverlap)+1
              endif
              tmpfile(idrfile(3,ii)) = idrfile(2,ii)
!              print*,"idrfile: ",ii,idrfile(2,ii),idrfile(3,ii),iifile,Maxiifile
            endif
            if(idrfile(1,ii).eq.4) then
               !this is used to avoid inifinite loop at the exit of the bend
               !where some particles behind the reference particle are still
               !inside the "bend". here, the "bend" is defined as a short section
               !of drift + bend + a short section of drift. 
               if(idbd.ne.ii) then
                 idbend = 1
                 idbd = ii
               else
                 idbend = 0
               endif
            else
               idbend = 0
            endif
          enddo 
          ibinitold = ibinit
          ibendold = ibend

          !calculate space-charge force for curr > 0 
          if((curr.gt.0.0).and.(flagspc.eq.1)) then !solve the Poission equation 

            ! get new boundary from the range of beam particles.
            !the use of grange instead of ptrange is due to multiple bunch
            call update_CompDom(Ageom,grange,grid2d,Flagbc)
            call getlcmnum_CompDom(Ageom,lcgrid)
            Nxlocal = lcgrid(1)
            if(npy.gt.1) then
              Nylocal = lcgrid(2) + 2
            else
              Nylocal = lcgrid(2)
            endif
            if(npx.gt.1) then
              Nzlocal = lcgrid(3) + 2
            else
              Nzlocal = lcgrid(3)
            endif

            call getlcrange_CompDom(Ageom,lcrange)
            !print*,"local range: ",lcrange,myid,myidx,myidy
            !move all effective particles to their local processor 
            if(totnp.gt.1) then
              do ib = 1, ibunch
                nplctmp = Nplocal(ib)
                call ptsmv2_Ptclmger(Ebunch(ib)%Pts1,nplctmp,grid2d,Pdim,&
                                Nplcmax,lcrange)
                Nplocal(ib) = nplctmp
                ! assign new 'Nplocal' local particles on each processor.
                call setnpt_BeamBunch(Ebunch(ib),Nplocal(ib))
              enddo
            endif

            !assign the storage for potential and charge density
            call set_FieldQuant(Potential,Nx,Ny,Nz,Ageom,grid2d,npx,&
                                npy)
            deallocate(chgdens)
            allocate(chgdens(Nxlocal,Nylocal,Nzlocal))

            !//initialize the space-charge fields.
            deallocate(exg)
            deallocate(eyg)
            deallocate(ezg)
            deallocate(bxg)
            deallocate(byg)
            deallocate(bzg)
            allocate(exg(Nxlocal,Nylocal,Nzlocal))
            allocate(eyg(Nxlocal,Nylocal,Nzlocal))
            allocate(ezg(Nxlocal,Nylocal,Nzlocal))
            allocate(bxg(Nxlocal,Nylocal,Nzlocal))
            allocate(byg(Nxlocal,Nylocal,Nzlocal))
            allocate(bzg(Nxlocal,Nylocal,Nzlocal))
            exg = 0.0
            eyg = 0.0
            ezg = 0.0
            bxg = 0.0
            byg = 0.0
            bzg = 0.0

            if(npx.gt.1) then
                nzlcr = Nzlocal-2
                kadd = 1
            else
                nzlcr = Nzlocal
                kadd = 0
            endif
            if(npy.gt.1) then
                nylcr = Nylocal-2
                jadd = 1
            else
                nylcr = Nylocal
                jadd = 0
            endif

            deallocate(tmppot)
            allocate(tmppot(Nxlocal,Nylocal,Nzlocal))
            call getrange_CompDom(Ageom,range)

            !//sum up the space-charge fields from all bunches/bins
            do ib = 1, ibunch
              ! deposit particles onto grid to obtain charge density of each bunch/bin.
              call chgdens_Depositor(Ebunch(ib),chgdens,Ageom,grid2d,&
                                     gammaz(ib),Flagbc,Perdlen,zcent,flagpos)

              !// solve 3D Poisson's equation for each bunch/bin
              if(Flagbc.eq.1) then
                ! solve Poisson's equation using 3D isolated boundary condition.
                call update3Ot_FieldQuant(Potential,chgdens,Ageom,&
                grid2d,Nxlocal,Nylocal,Nzlocal,npx,npy,nylcr,nzlcr,gammaz(ib))
              else
                print*,"no such boundary condition type!!!"
                stop
              endif

              !find the E and B fields in the lab frame from the effective bunch/bin  
              tmpflag = 0

              call gradEB_FieldQuant(Nxlocal,Nylocal,Nzlocal,&
              Potential%FieldQ,Ageom,grid2d,Flagbc,gammaz(ib),tmpflag,&
              exg,eyg,ezg,bxg,byg,bzg)

!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            enddo

            !interpolate the space-charge fields using CIC + external fields to
            !each particle
            !print*,"exg: ",sum(exg),sum(eyg),sum(ezg),sum(bxg),sum(byg),sum(bzg)
            !test wakefield
            !bxg = 0.0
            !byg = 0.0
            !bzg = 0.0
            do ib = 1, ibunch
              call kick2t_BeamBunch(Nplocal(ib),Nxlocal,Nylocal,Nzlocal,&
              Ebunch(ib)%Pts1,exg,eyg,ezg,bxg,byg,bzg,Ageom,npx,npy,myidx,&
              myidy,t,Ebunch(ib)%Charge,Ebunch(ib)%Mass,dtless,Blnelem,&
              zBlnelem,idrfile,Nblem,ibinit,ibend,fldmp,Flagerr)
            enddo
          else !no space-charge fields, only external fields
            do ib = 1, ibunch
              call scatter20t_BeamBunch(Nplocal(ib),Ebunch(ib)%Pts1,t,&
                 Ebunch(ib)%Charge,&
                 Ebunch(ib)%Mass,dtless,Blnelem,zBlnelem,idrfile,Nblem,&
                 ibinit,ibend,fldmp,Flagerr)
            enddo
          endif
          !//update particle positions using new velocity for half step.
          if(flagcathode.eq.1) then
            do ib = 1, Nbunch
              call drifthalf_BeamBunch(Ebunch(ib),t,dtless,betazini)
            enddo
          else
            do ib = 1, Nbunch
              call drifthalforg_BeamBunch(Ebunch(ib),t,dtless)
            enddo
          endif
          t = t + 0.5*dtless*Dt
          !if(mod(i,5).eq.0) then           
          if(mod(i,1).eq.0) then 
             !!call diagnostic1avg_Output(t,Ebunch,Nbunch)
          endif

!          if(t.le.tphout(iout+1) .and. (t+dtless*Dt).ge.tphout(iout+1)) then
          if(distance.le.tphout(iout+1) .and. &
            (distance+dzz).ge.tphout(iout+1)) then
            iout = iout + 1
            !write(*, *), " TEMPORARILY DISABLE PHASE_OUT of ",bpname2
            do ib = 1, Nbunch
              !!call phase_Output(bpname2, nfileouttmp(iout)+ib-1,Ebunch(ib),nsamp(iout))
            enddo
          endif

          if(distance.le.trstart .and. (distance+dzz).ge.trstart) then
            do ib = 1, Nbunch
              call outpoint_Output(ib*nfileout+myid,Ebunch(ib),t,&
                        i,ib,npx,npy,Ageom,iout,itsz,isteer)
            enddo
          endif
         
        enddo  ! end of timestep loop

! final output.
        call MPI_BARRIER(comm2d,ierr)
        !drift back half time step
        do ib = 1, Nbunch   
          call drifthalf_BeamBunch(Ebunch(ib),t,-dtless,betazini)
          do ipt = 1, Nplocal(ib)
            deltaz = blnLength/Scxlt - Ebunch(ib)%Pts1(5,ipt)
            Ebunch(ib)%Pts1(1,ipt) = Ebunch(ib)%Pts1(1,ipt)+&
                     Ebunch(ib)%Pts1(2,ipt)/Ebunch(ib)%Pts1(6,ipt)*deltaz
            Ebunch(ib)%Pts1(3,ipt) = Ebunch(ib)%Pts1(3,ipt)+&
                     Ebunch(ib)%Pts1(4,ipt)/Ebunch(ib)%Pts1(6,ipt)*deltaz
          enddo
        enddo   
        !output six 2-D phase projections.
        !call phase2dold_Output(30,Ebunch,Np)
        !output all particles in 6d phase space at given location blnLength.
        !write(*, *), " TEMPORARILY DISABLE PHASE_OUT of ",bpname3, "NBunch=", Nbunch
        do ib = 1, Nbunch
          i = 50+ib-1 
          !!call phase_Output(bpname3, i,Ebunch(ib),1)
          !!call phaseold_Output(i,Ebunch(ib),1)
        enddo

        t_integ = t_integ + elapsedtime_Timer(t0)
        call showtime_Timer()

        deallocate(chgdens)
        deallocate(tmppot)
        deallocate(exg)
        deallocate(eyg)
        deallocate(ezg)
        deallocate(bxg)
        deallocate(byg)
        deallocate(bzg)
        deallocate(idrfile)
        deallocate(gammaz)
        deallocate(brange)

        !!call phase_in(Ebunch(ib), bpname3);
        end subroutine run_AccSimulator


        subroutine destruct_AccSimulator(time)
        implicit none
        include 'mpif.h'
        double precision :: time
        integer :: ib
 
        do ib = 1, Nbunch
          call destruct_BeamBunch(Ebunch(ib))
        enddo
        call destruct_FieldQuant(Potential)
        call destruct_CompDom(Ageom)

        call end_Output(time)

        end subroutine destruct_AccSimulator

      end module AccSimulatorclass