!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer. 
! Version: 1.0-beta
! Author: Ji Qiang, LBNL
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments:
!----------------------------------------------------------------
      module Outputclass
        use readutil
        use Timerclass
        use BeamBunchclass
        use PhysConstclass
      contains
        ! calculate averaged <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
        ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance from
        ! multiple bunch/bin.
        subroutine diagnostic1avg_Output(z,this,Nbunch)
        implicit none
        include 'mpif.h'
        double precision, intent(in) :: z
        type (BeamBunch), dimension(:), intent(inout) :: this
        integer, intent(in) :: Nbunch
        integer :: innp,nptot
        double precision:: den1,den2,sqsum1,sqsum2,sqsum3,sqsum4,&
                          epsx2,epsy2
        double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
                         xpxfac,ypyfac
        double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
        double precision:: xpxlocal,ypylocal,zpzlocal
        double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
        double precision:: sqsum5local,sqsum6local
        double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
        pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
        z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
        z03,z04,pz03,pz04
        double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
        sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
        double precision:: gam,energy,bet
        integer :: i,my_rank,ierr,j
        double precision:: qmc,xl,xt
        double precision, dimension(6) :: localmax, glmax
        double precision, dimension(29) :: tmplc,tmpgl
        double precision :: t0,lcrmax,glrmax,z0gl,z0avg,testmax,pz0avg
        double precision :: gamlc,gam2lc,gam2avg,tmpgam,gamdel
        integer :: npctmin,npctmax,ib,innpmb,i1,i2
        real*8 :: gamavg,gamgl
        real*8, dimension(3) :: tmp3lc,tmp3gl

        call starttime_Timer(t0)

        qmc = this(1)%Mass/1.0e6
        xl = Scxlt
        xt = Rad2deg

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

        nptot = 0
        innpmb = 0
        z0lc = 0.0
        gamlc = 0.0d0
        pz0lc = 0.0
        tmp3lc = 0.0d0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          innpmb = innpmb + innp
          nptot = nptot + this(ib)%Npt
          do i = 1, innp
            tmp3lc(1) = tmp3lc(1) + this(ib)%Pts1(5,i)
            tmp3lc(2) = tmp3lc(2) + this(ib)%Pts1(6,i)
            tmp3lc(3) = tmp3lc(3) + sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
          enddo
        enddo

        call MPI_ALLREDUCE(tmp3lc,tmp3gl,3,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,MPI_COMM_WORLD,ierr)

        den1 = 1.0d0/nptot
        den2 = den1*den1
        z0avg = tmp3gl(1)/nptot
        pz0avg = tmp3gl(2)/nptot
        gamavg = tmp3gl(3)/nptot
        !print*,"z0avg: ",z0avg,pz0avg,gamavg


        x0lc = 0.0
        px0lc = 0.0
        y0lc = 0.0
        py0lc = 0.0
        pz0lc = 0.0
        sqsum1local = 0.0
        sqsum2local = 0.0
        sqsum3local = 0.0
        sqsum4local = 0.0
        sqsum5local = 0.0
        sqsum6local = 0.0
        xpxlocal = 0.0
        ypylocal = 0.0
        zpzlocal = 0.0
        x0lc3 = 0.0
        x0lc4 = 0.0
        px0lc3 = 0.0
        px0lc4 = 0.0
        y0lc3 = 0.0
        y0lc4 = 0.0
        py0lc3 = 0.0
        py0lc4 = 0.0
        z0lc3 = 0.0
        z0lc4 = 0.0
        pz0lc3 = 0.0
        pz0lc4 = 0.0
        ! for cache optimization.
        if(innp.ne.0) then
          do i = 1, 6
            !localmax(i) = abs(this(1)%Pts1(i,1))
            localmax(i) = 0.0
          enddo
          !lcrmax = this(1)%Pts1(1,1)**2+this(1)%Pts1(3,1)**2
          lcrmax = 0.0
        else
          do i = 1, 6
            localmax(i) = 0.0
          enddo
          lcrmax = 0.0
        endif
!        testmax = 0.0
        gamlc = 0.0
        gam2lc = 0.0
        do ib = 1, Nbunch
          innp = this(ib)%Nptlocal
          do i = 1, innp
            tmpgam = sqrt(1.0+this(ib)%Pts1(2,i)**2+&
               this(ib)%Pts1(4,i)**2+this(ib)%Pts1(6,i)**2)
            gamlc = gamlc + tmpgam
            !gam2lc = gam2lc + tmpgam**2
            gam2lc = gam2lc + (tmpgam-gamavg)**2
            x0lc = x0lc + this(ib)%Pts1(1,i)
            sqsum1local = sqsum1local + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc3 = x0lc3 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)
            x0lc4 = x0lc4 + this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*this(ib)%Pts1(1,i)*&
                    this(ib)%Pts1(1,i)
            xpxlocal = xpxlocal + this(ib)%Pts1(1,i)*this(ib)%Pts1(2,i)
            px0lc = px0lc + this(ib)%Pts1(2,i)
            sqsum2local = sqsum2local + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc3 = px0lc3 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)
            px0lc4 = px0lc4 + this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*this(ib)%Pts1(2,i)*&
                     this(ib)%Pts1(2,i)
            y0lc = y0lc + this(ib)%Pts1(3,i)
            sqsum3local = sqsum3local + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc3 = y0lc3 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)
            y0lc4 = y0lc4 + this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*this(ib)%Pts1(3,i)*&
                    this(ib)%Pts1(3,i)
            ypylocal = ypylocal + this(ib)%Pts1(3,i)*this(ib)%Pts1(4,i)
            py0lc = py0lc + this(ib)%Pts1(4,i)
            sqsum4local = sqsum4local + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc3 = py0lc3 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)
            py0lc4 = py0lc4 + this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*this(ib)%Pts1(4,i)*&
                     this(ib)%Pts1(4,i)
            sqsum5local = sqsum5local + (this(ib)%Pts1(5,i)-z0avg)**2
            z0lc3 = z0lc3 + abs((this(ib)%Pts1(5,i)-z0avg)**3)
            z0lc4 = z0lc4 + (this(ib)%Pts1(5,i)-z0avg)**4
            zpzlocal = zpzlocal + (this(ib)%Pts1(5,i)-z0avg)*(this(ib)%Pts1(6,i)-pz0avg)
            pz0lc = pz0lc + this(ib)%Pts1(6,i)
            sqsum6local = sqsum6local + (this(ib)%Pts1(6,i)-pz0avg)**2 
            pz0lc3 = pz0lc3 + (abs(this(ib)%Pts1(6,i)-pz0avg)**3) 
            pz0lc4 = pz0lc4 + (this(ib)%Pts1(6,i)-pz0avg)**4 
            do j = 1, 4
              if(localmax(j).lt.abs(this(ib)%Pts1(j,i))) then
                 localmax(j) = abs(this(ib)%Pts1(j,i))
              endif
            enddo
            if(localmax(5).lt.abs(this(ib)%Pts1(5,i)-z0avg)) then
               localmax(5) = abs(this(ib)%Pts1(5,i)-z0avg)
               i1 = i
            endif
!            if(testmax .lt. abs(this(ib)%Pts1(5,i)) ) then
!              testmax = this(ib)%Pts1(5,i)
!              i2 = i
!              print*,"testmax: ",i,testmax,this(ib)%Pts1(5,i),&
!                        abs(this(ib)%Pts1(5,i))
!            endif
            if(localmax(6).lt.abs(this(ib)%Pts1(6,i)-pz0avg)) then
                 localmax(6) = abs(this(ib)%Pts1(6,i)-pz0avg)
            endif
            if(lcrmax.lt.(this(ib)%Pts1(1,i)**2+this(ib)%Pts1(3,i)**2)) then
              lcrmax = this(ib)%Pts1(1,i)**2 + this(ib)%Pts1(3,i)**2
            endif
          enddo
        enddo

        !print*,"z0avg: ",z0avg*xl,localmax(5)*xl,&
        !       this(1)%Pts1(5,i1)*xl
        tmplc(1) = x0lc
        tmplc(2) = px0lc
        tmplc(3) = y0lc
        tmplc(4) = py0lc
        tmplc(5) = z0lc
        tmplc(6) = pz0lc
        tmplc(7) = sqsum1local
        tmplc(8) = sqsum2local
        tmplc(9) = sqsum3local
        tmplc(10) = sqsum4local
        tmplc(11) = sqsum5local
        tmplc(12) = sqsum6local
        tmplc(13) = xpxlocal
        tmplc(14) = ypylocal
        tmplc(15) = zpzlocal
        tmplc(16) = x0lc3
        tmplc(17) = x0lc4
        tmplc(18) = px0lc3
        tmplc(19) = px0lc4
        tmplc(20) = y0lc3
        tmplc(21) = y0lc4
        tmplc(22) = py0lc3
        tmplc(23) = py0lc4
        tmplc(24) = z0lc3
        tmplc(25) = z0lc4
        tmplc(26) = pz0lc3
        tmplc(27) = pz0lc4
        tmplc(28) = gamlc
        tmplc(29) = gam2lc
        
        call MPI_REDUCE(tmplc,tmpgl,29,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
                        MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(innpmb,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        if(my_rank.eq.0) then
          x0 = tmpgl(1)*den1
          px0 = tmpgl(2)*den1
          y0 = tmpgl(3)*den1
          py0 = tmpgl(4)*den1
          z0 = tmpgl(5)*den1
          pz0 = tmpgl(6)*den1
          sqx = tmpgl(7)*den1
          sqsum1 = sqx - x0*x0
          sqpx = tmpgl(8)*den1
          sqsum2 = sqpx - px0*px0
          sqy = tmpgl(9)*den1
          sqsum3 = sqy - y0*y0
          sqpy = tmpgl(10)*den1
          sqsum4 = sqpy - py0*py0
          sqz = tmpgl(11)*den1
!          sqsum5 = sqz - z0*z0
          sqsum5 = sqz 
          sqpz = tmpgl(12)*den1
!          sqsum6 = sqpz - pz0*pz0
          sqsum6 = sqpz
          xpx = tmpgl(13)*den1 - x0*px0
          ypy = tmpgl(14)*den1 - y0*py0
!          zpz = tmpgl(15)*den1 - z0*pz0
          zpz = tmpgl(15)*den1 
          cubx = tmpgl(16)*den1
          fthx = tmpgl(17)*den1
          x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
          x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
          cubpx = tmpgl(18)*den1
          fthpx = tmpgl(19)*den1
          px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
          px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
                 3*px0*px0*px0*px0)))
          cuby = tmpgl(20)*den1
          fthy = tmpgl(21)*den1
          y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
          y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
          cubpy = tmpgl(22)*den1
          fthpy = tmpgl(23)*den1
          py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
          py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
                 3*py0*py0*py0*py0)))
          cubz = tmpgl(24)*den1
          fthz = tmpgl(25)*den1
!          z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
          z03 = cubz**(1.0/3.0)
!          z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
          z04 = sqrt(sqrt(fthz))
          cubpz = tmpgl(26)*den1
          fthpz = tmpgl(27)*den1
          pz03 = cubpz**(1.0/3.0)
          pz04 = sqrt(sqrt(fthpz))
          !pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
          !pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
          !       3*pz0*pz0*pz0*pz0)))
          epsx2 = (sqsum1*sqsum2-xpx*xpx)
          epsy2 = (sqsum3*sqsum4-ypy*ypy)
          epsz2 = (sqsum5*sqsum6-zpz*zpz)
          epx = sqrt(max(epsx2,0.0d0))
          epy = sqrt(max(epsy2,0.0d0))
          epz = sqrt(max(epsz2,0.0d0))
          xrms = sqrt(abs(sqsum1))
          pxrms = sqrt(abs(sqsum2))
          yrms = sqrt(abs(sqsum3))
          pyrms = sqrt(abs(sqsum4))
          zrms = sqrt(abs(sqsum5))
          pzrms = sqrt(abs(sqsum6))
          xpxfac = 0.0
          ypyfac = 0.0
          zpzfac = 0.0
          if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
          if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
          if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
          gam = tmpgl(28)*den1
!          gam = sqrt(1.0+px0**2+py0**2+pz0**2)
          energy = qmc*(gam-1.0)
          bet = sqrt(1.0-(1.0/gam)**2)
          gam2avg = tmpgl(29)*den1
          !gamdel = sqrt(abs(gam2avg - gam**2))
          gamdel = sqrt(gam2avg)

          write(*,*) "z0avg=", z0avg, "xl=", xl, "epz=", epz;
          write(*,*) "sqsum5=", sqsum5, "sqsum6=", sqsum6, "zpz=", zpz
          write(*,*) "z03=", z03, "sumcubz=", tmpgl(24), "cubz=", cubz, "fthz=", fthz;
          write(18,100)z,z0avg*xl,gam,energy,bet,sqrt(glrmax)*xl,gamdel

          write(24,102)z,z0avg*xl,x0*xl,xrms*xl,px0,pxrms,-xpx*xl,epx*xl
          write(25,102)z,z0avg*xl,y0*xl,yrms*xl,py0,pyrms,-ypy*xl,epy*xl
          write(26,100)z,z0avg*xl,zrms*xl,pz0,pzrms,-zpz*xl,epz*xl

          write(27,102)z,z0avg*xl,glmax(1)*xl,glmax(2),glmax(3)*xl,&
                   glmax(4),glmax(5)*xl,glmax(6)
          write(28,101)z,z0avg*xl,npctmin,npctmax,nptot
          write(29,102)z,z0avg*xl,x03*xl,px03,y03*xl,py03,z03*xl,&
                   pz03
          write(30,102)z,z0avg*xl,x04*xl,px04,y04*xl,py04,z04*xl,&
                   pz04

          call flush(18)
          call flush(24)
          call flush(25)
          call flush(26)
          call flush(27)
          call flush(28)
          call flush(29)
          call flush(30)
        endif

99      format(6(1x,e16.8))
100      format(7(1x,e18.10))
101     format(1x,e16.8,e16.8,3I10)
102      format(8(1x,e16.8))

        t_diag = t_diag + elapsedtime_Timer(t0)

        end subroutine diagnostic1avg_Output

       !the 6D phase space output has (x(m), px/mc, y(m), py/mc, z(m), pz/mc).
       ! Junmin 
       subroutine phase_Output(bpFileName, bunchID,this,samplePeriod, timestep, totalstep)
         !!use sensei_util
        implicit none
        include 'mpif.h'
        integer, intent(in) :: bunchID
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt,sumnpt,startnpt;
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf
        character(20), intent(in) :: bpFileName;
        character(20) :: senseiConfigFile="./config.xml"

        double precision :: endtime, starttime, duration;

        real*8 ::time = 0.0;
        integer, intent(in) :: timestep, totalstep
        integer :: callSensei = 1

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,sumnpt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
             
        allocate(nptlist(0:np-1))
        nptlist = 0

        !!call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)                        
        call MPI_ALLGATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)                        

        if (my_rank .eq. 0) then
           !write(*,*), "timestep[ = ", timestep, "] sum=", sumnpt, " and nptlist=", nptlist
        endif

        do i=0, np-1           
           if (nptlist(i) .eq. 0) then
              callSensei = 0;
           end if
        end do

        starttime = MPI_WTIME()

        startnpt = 0;
        if (my_rank > 0) then 
           do i=0,my_rank-1
              startnpt = nptlist(i) + startnpt;
           enddo
        endif

        !write (*,*) 'rank=', my_rank, bpFileName, ' start at: ', startnpt, ' numpts = ', this%Nptlocal, this%Pts1(1,1)
        call  utilAdiosOut(bpFileName, MPI_COMM_WORLD, this%Pts1, sumnpt, 6, this%Nptlocal, startnpt, timestep, totalstep, bunchID);                
           
        if (callSensei .gt. -1) then
           !write(*, *) 'rank = ', my_rank, "timestep=", timestep, this%Nptlocal;
           !!call histogramOutInt(MPI_COMM_WORLD, timestep, time, this%Pts1, sumnpt, 6, this%Nptlocal, startnpt, my_rank, np);
        else
           if (my_rank .eq. 1) then
              write(*,*), "skipping timestep=", timestep
           endif
        endif

        call MPI_Barrier (MPI_COMM_WORLD, ierr)

        deallocate(nptlist)


        endtime = MPI_WTIME()
        time = endtime - starttime
        call MPI_REDUCE(time,duration,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        if (my_rank .eq. 0) then
           write(*,*), "phaseout at: ", timestep, " took:", duration
        endif
        call MPI_REDUCE(time,duration,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
        if (my_rank .eq. 0) then
           write(*,*), "phaseout at: ", timestep, " minTime:", duration
        endif
        
        end subroutine phase_Output





       subroutine phaseBunchOutput(bpFileName, Nbunch,this,samplePeriod, timestep, totalstep)
         !!use sensei_util
        implicit none
        include 'mpif.h'
        integer, intent(in) :: Nbunch
        !!type (BeamBunch), intent(in) :: this
        type (BeamBunch), dimension(Nbunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,ib,startnpt, mm;
        integer, dimension(Nbunch)::localnpt, mnpts,sumnpts;

        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf
        character(20), intent(in) :: bpFileName;
        character(20) :: senseiConfigFile="./config.xml"

        double precision :: endtime, starttime, duration;

        real*8 ::time = 0.0;
        integer, intent(in) :: timestep, totalstep
        integer*8 adios_handle;
        integer, dimension(4) :: bunchInput;


        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

        do ib=1, Nbunch
           localnpt(ib)=this(ib)%Nptlocal
        enddo

        call MPI_ALLREDUCE(localnpt,mnpts,  Nbunch,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(localnpt,sumnpts,Nbunch,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0

        starttime = MPI_WTIME()

        bunchInput(1)=timestep;
        bunchInput(2)=totalstep;
        bunchInput(3)=1;
        bunchInput(4)=Nbunch;

        call testMeOpen(bpFileName, MPI_COMM_WORLD, bunchInput, adios_handle)
        do ib=1, Nbunch
           call MPI_ALLGATHER(this(ib)%Nptlocal,1,MPI_INTEGER,nptlist,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)                        

           startnpt = 0;
           if (my_rank > 0) then 
              do i=0,my_rank-1
                 startnpt = nptlist(i) + startnpt;
              enddo
           endif
           
           mm=sumnpts(ib);
           bunchInput(3)=ib;

           !!call utilAdiosBunchOut(bpFileName,MPI_COMM_WORLD,this(ib)%Pts1,mm,6,this(ib)%Nptlocal,startnpt,bunchInput);
           call testMe(bpFileName,MPI_COMM_WORLD,this(ib)%Pts1,mm,6,this(ib)%Nptlocal,startnpt,bunchInput, adios_handle);
           
        enddo
        call MPI_Barrier (MPI_COMM_WORLD, ierr)
        deallocate(nptlist)



        endtime = MPI_WTIME()
        time = endtime - starttime
        call MPI_REDUCE(time,duration,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        if (my_rank .eq. 0) then
           write(*,*), "phaseout at: ", timestep, " took:", duration
        endif
        call MPI_REDUCE(time,duration,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
        if (my_rank .eq. 0) then
           write(*,*), "phaseout at: ", timestep, " minTime:", duration
        endif
        
        end subroutine phaseBunchOutput





        ! Terminate MPI
        subroutine end_Output(time)
        implicit none
        include 'mpif.h'
        !double precision MPI_WTIME
        double precision, intent(inout) :: time
        double precision :: endtime, mtime
        integer :: my_rank,ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        endtime = MPI_WTIME()
        time = endtime - time
        call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
                        MPI_COMM_WORLD,ierr)

        !for measurement of memory
        !call system_stats()

        if(my_rank.eq.0) then
          print*,"time: ",mtime
        endif

!        call MPI_Finalize(ierr)

        end subroutine end_Output

        subroutine inpoint_Output(nfile,this,z,inb,jstp,nprocrow,nproccol,&
                                  geom,nx,ny,nz,myidx,myidy,nptot,iout,itsz,isteer)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(inout) :: this
        double precision, intent(inout) :: z
        type(CompDom), intent(inout) :: geom
        integer, intent(inout) :: inb,jstp
        integer, intent(inout) :: nptot,iout,itsz,isteer
        integer, intent(in) :: nprocrow,nproccol,nx,ny,nz,myidx,myidy
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range
        integer :: nplc,ierr
        character*6 name1
        character*7 name2
        character*8 name3
        character*9 name4
        integer :: i,j,k,l,m,n

        name1 = 'fort.x'
        name2 = 'fort.xx'
        name3 = 'fort.xxx'
        name4 = 'fort.xxxx'

        write (*,*), "=========== inpoint ========"
        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))

!        open(nfile,status="old",form="unformatted")
        if(nfile < 10) then
            name1(6:6) = char(nfile+48)
            open(9,file=name1,status="unknown",form="unformatted")
        else if((nfile.ge.10).and.(nfile.lt.100)) then
            i = nfile/10
            j = nfile - 10*i
            name2(6:6) = char(i+48)
            name2(7:7) = char(j+48)
            open(9,file=name2,status="unknown",form="unformatted")
        else if((nfile.ge.100).and.(nfile.lt.1000)) then
            i = nfile/100
            j = nfile - 100*i
            k = j/10
            l = j - 10*k
            name3(6:6) = char(i+48)
            name3(7:7) = char(k+48)
            name3(8:8) = char(l+48)
            open(9,file=name3,status="unknown",form="unformatted")
        else
            i = nfile/1000
            j = nfile - 1000*i
            k = j/100
            l = j - 100*k
            m = l/10
            n = l - 10*m
            name4(6:6) = char(i+48)
            name4(7:7) = char(k+48)
            name4(8:8) = char(m+48)
            name4(9:9) = char(n+48)
            open(9,file=name4,status="unknown",form="unformatted")
        endif
 
!        print*,"nfile: ",nfile
        read(9)z
        read(9)inb,jstp,iout,itsz,isteer
        read(9)msize(1:3)
        read(9)range(1:6)
!        print*,"zz: ",z,inb,jstp
        read(9)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
        read(9)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        read(9)this%Nptlocal
        read(9)this%refptcl(1:6)
        allocate(this%Pts1(6,this%Nptlocal))
        read(9)this%Pts1(1:6,1:this%Nptlocal)

        close(9)

        geom%Meshsize = msize
        geom%SpatRange = range
        allocate(geom%LcTabrg(4,0:nprocrow-1,0:nproccol-1))
        allocate(geom%LcTabnm(2,0:nprocrow-1,0:nproccol-1))
        geom%lcTabnm = Localnum
        geom%LcTabrg = Localrange

        deallocate(Localnum)
        deallocate(Localrange)

        geom%Meshnum(1) = nx
        geom%Meshnum(2) = ny
        geom%Meshnum(3) = nz

        geom%Mshlocal(3) = geom%LcTabnm(1,myidx,myidy)
        geom%Mshlocal(2) = geom%LcTabnm(2,myidx,myidy)
        geom%Mshlocal(1) = geom%Meshnum(1)

        geom%Sptrnglocal(5) = geom%LcTabrg(1,myidx,myidy)
        geom%Sptrnglocal(6) = geom%LcTabrg(2,myidx,myidy)
        geom%Sptrnglocal(3) = geom%LcTabrg(3,myidx,myidy)
        geom%Sptrnglocal(4) = geom%LcTabrg(4,myidx,myidy)
        geom%Sptrnglocal(1) = geom%SpatRange(1)
        geom%Sptrnglocal(2) = geom%SpatRange(2)

        nplc = this%Nptlocal
        call MPI_ALLREDUCE(nplc,nptot,1,MPI_INTEGER,&
                           MPI_SUM,MPI_COMM_WORLD,ierr)
        this%Npt = nptot

        end subroutine inpoint_Output

        subroutine outpoint_Output(nfile,this,z,inb,jstp,nprocrow,nproccol,&
                                   geom,iout,itsz,isteer)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        double precision, intent(in) :: z
        type(CompDom), intent(in) :: geom
        integer, intent(in) :: inb,jstp,nprocrow,nproccol,&
                               iout,itsz,isteer
        integer, allocatable, dimension(:,:,:) :: Localnum
        double precision, allocatable, dimension(:,:,:) :: Localrange
        double precision, dimension(3) :: msize
        double precision, dimension(6) :: range

        allocate(Localnum(2,0:nprocrow-1,0:nproccol-1))
        call getlctabnm_CompDom(geom,Localnum)

        allocate(Localrange(4,0:nprocrow-1,0:nproccol-1))
        call getlctabrg_CompDom(geom,Localrange)
  
        call getmsize_CompDom(geom,msize)
        call getrange_CompDom(geom,range)

!hjw
! All the following I/O commented out by hjw
!       write(*,*)' Opening file: ',nfile,' in outpoint_Output'
!       open(nfile,status="unknown",form="unformatted")

!       write(nfile)z
!       write(nfile)inb,jstp,iout,itsz,isteer
!       write(nfile)msize(1:3)
!       write(nfile)range(1:6)
!       write(nfile)Localnum(1:2,0:nprocrow-1,0:nproccol-1)
!       write(nfile)Localrange(1:4,0:nprocrow-1,0:nproccol-1)

        deallocate(Localnum)
        deallocate(Localrange)

!       write(nfile)this%Nptlocal
!       write(nfile)this%refptcl(1:6)
!       write(nfile)this%Pts1(1:6,1:this%Nptlocal)

!       close(nfile)

        end subroutine outpoint_Output

subroutine phase_in(this, bpFilename,nptlc,npt)
    implicit none
    include 'mpif.h'

    character(len=20), intent(IN):: bpFileName
    integer my_rank, ierr, np, sumnpt;
    type (BeamBunch), intent(INOUT) :: this
    real*8, dimension(:,:), allocatable :: pts 
    integer :: nptlc,npt,i
    double precision :: endtime, starttime, duration;
    real*8 ::time = 0.0;

    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

    starttime = MPI_WTIME()
        
    call utilAdiosIn(bpFileName, MPI_COMM_WORLD, my_rank, np, pts, nptlc,npt);
    !write(*,*) "rank=", my_rank, pts(1,1), pts(1,2), pts(2,1);

    endtime = MPI_WTIME()
    time = endtime - starttime
    call MPI_REDUCE(time,duration,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    if (my_rank .eq. 0) then
       write(*,*), "phasein: ", bpFilename, " duration:",  duration
    endif


    this%Nptlocal = nptlc
    this%Npt = npt
    allocate(this%Pts1(6,this%Nptlocal))
    do i = 1, nptlc
      this%Pts1(:,i) = pts(:,i)
    enddo

    deallocate(pts);

end subroutine phase_in

!read in the particle data file "particles.bp" using ADIOS
!and put it onto each processor. (by Junmin Gu)
subroutine phaseinadios(this)
    use adios_read_mod
    implicit none
    include 'mpif.h'
    type (BeamBunch), intent(inout) :: this
   
    character(len=50)   :: filename = "particle.bp"

    integer             :: rank, size, i, j, ierr
    integer             :: comm

    integer*8               :: f
    integer                 :: method = ADIOS_READ_METHOD_BP
    integer*8               :: blockSel, boxSel

    ! variables to read in 
    integer                 :: NX, NY, GY

    integer                 :: blocksPP   !! blocks per process
    integer*8, dimension(2) :: boxOffset=0, boxSize=1

    real*8, dimension(:,:), allocatable :: pts
    

    INTEGER, PARAMETER :: N = 6
    INTEGER OTHER_PE,  JERR
    INTEGER SEND, RECV
    INTEGER STATUS(MPI_STATUS_SIZE)
    INTEGER, DIMENSION(N) :: RBUF, SBUF


! =======================================
! read in bp file, otherwise, use default in particle.bp
! =======================================
!    if (command_argument_count() > 0) then
!         call getarg(1, filename);
!    endif


!    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)


! =======================================
! initialize adios read
! =======================================
   call adios_read_init_method (method, comm, "verbose=3", ierr);
   call adios_read_open (f, filename, method, comm, ADIOS_LOCKMODE_NONE, 1.0, ierr);

! =======================================
! use rank in the output filename for each processor
! =======================================

! =======================================
! prepare send buffer (rank 0 reads in bp metadata and figure out workload for 
! each processor
! =======================================
    sbuf(1)=0; ! number of blocks to read
    sbuf(2)=0; ! how many processors will be on one block
    sbuf(3)=0; ! rows to read per block
    sbuf(4)=0;
    sbuf(5)=0;
    sbuf(6)=0;

    rbuf = 0;


! =======================================
! In rank 0, we read out the size of the data (determined by NX/NY/GY), and figures out how much each processor should read out from bp and write to txt
! each procossor will either read out a few blocks (each block contains data from one fortran binary file), or several processors will split out reading on one block.
! the info is passed from rank 0 to other processors. Once received instruction, each processor will start reading out from bp file
! buf(1) is the # of blocks each processor will read. if = 0, will look at buf(2)
! buf(2) is how many processors will work on reading out one block
! buf(3) is the size of rows to read out per processor
! =======================================
if (rank == 0) then
    print *, 'reading from file: ', filename

    call adios_selection_writeblock (blockSel, rank) 

! =======================================
! First get the scalars to calculate the size of the arrays.
! Note that we cannot use adios_get_scalar here because that
!   retrieves the same NX for everyone (from writer rank 0).
! =======================================
    call adios_schedule_read (f, blockSel, "NX", 0, 1, NX, ierr)
    call adios_schedule_read (f, blockSel, "NY", 0, 1, NY, ierr)
    call adios_schedule_read (f, blockSel, "GY", 0, 1, GY, ierr)
    call adios_perform_reads (f, ierr)

!    write (*,'(" rank=",i0," size=",i0," NX=",i0," NY=",i0, " GY=",i0)') rank, size, NX, NY, GY

! =======================================
! share with other processors
! =======================================
    sbuf(4)=NX;
    sbuf(5)=NY;
    sbuf(6)=GY;

! =======================================
! figure out workload
! =======================================
    if (GY/NY >= size) then 
       print *, "======== check blocks ======"
       if (mod(GY/NY, size) /= 0) then
           print *, ' size needs to be a divisor of ', GY/NY, " Exit."
       else 
       	   sbuf(1) = GY/(NY*size);
       endif 
    else 
       print *, GY,"%",size, "=", mod(GY, size)
       if (mod(GY, size) /= 0) then
       	  print *, "size ", size, " needs to divide ", GY, "Exit."
       else 
          if (mod (NY, GY/size) /= 0) then
	      print *, "size ", size, " needs to be multiple of ", GY/NY
	  else 
	      sbuf(1) = 0;
	      sbuf(2) = NY/(GY/size)
	      sbuf(3) = (GY/size)
	  endif 
       endif       
    endif

    if (size > 1) then
       do i=1,size-1
         CALL MPI_SEND(SBUF, N,  MPI_INTEGER, i, 99, MPI_COMM_WORLD, SEND)
       enddo

       IF (SEND /= MPI_SUCCESS) then
       STOP 'BAD SEND on rank 0'
       end if
     endif	

else
! =======================================
! now processor is not of rank 0. 
! receive infor from rank before proceeding
! =======================================
     CALL MPI_RECV(RBUF, N, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, STATUS, RECV)

     IF (RECV /= MPI_SUCCESS) then
        STOP 'BAD receive on rank 1'
     endif
     
     NX=RBUF(4);
     NY=RBUF(5);
     GY=RBUF(6);

     if ((RBUF(1) == 0) .AND. (RBUF(2) /= 0)) then	
     	print *, "process rank/size=", rank, size, RBUF, "  in block: ", rank/RBUF(2)
     else if (RBUF(1) /= 0) then
        print *, "process rank/size=", rank, size, RBUF, "read from blocks"
     else 
        print *, "process rank=", rank, RBUF, "exit."     
     endif 
endif 


! =======================================
! all processors will now process data according 
! to the workload. 
! either reading by blocks or by boundingbox 
! =======================================

    if ((SBUF(1) /= 0) .OR. (RBUF(1) /= 0)) then        
       !===============================================
       ! read out blocks, each block has size: (NX, NY)
       ! a  block was read from a fortran binary file
       !===============================================
       blocksPP = GY/(NY*size);

       print *, 'read block: rank/GY/NY/size/blocksPP=', rank, GY, NY, size, blocksPP, rank*blocksPP, (rank+1)*blocksPP

       if (blocksPP > 0) then
         ! Allocate space for the array to collect block data
         !allocate (pts(NX,NY))
         this%Nptlocal = NY
         print*,"Nptlocal: ",this%Nptlocal,NX
         deallocate(this%Pts1)
         allocate(this%Pts1(NX,this%Nptlocal))

         do j=rank*blocksPP, (rank+1)*blocksPP-1
	    !===============================================
	    ! create the block and read out. block number is calculated according to the way data are written by read.F90 
	    !===============================================
            call adios_selection_writeblock (blockSel, j)
            !call adios_schedule_read (f, blockSel, "particles", 1, 1, pts, ierr)
            call adios_schedule_read (f, blockSel, "particles", 0, 1, this%Pts1, ierr)
            call adios_perform_reads (f, ierr)

            !print*,"this%Pts1:1",this%Pts1(:,1) 
            !print*,"this%Pts1:2",this%Pts1(:,2) 

         end do

         call adios_selection_delete(blockSel);
       else
         print *, "!!Not expected to be here!!", SBUF(1), RBUF(1)
       endif
100    format(9(1x,e16.9))
    else if ((SBUF(2) /= 0) .OR. (RBUF(2) /= 0)) then
       ! ==========================
       ! retrieve part of the block. 
       ! we know in this example that block X (0<=X<=(GY/NY)) has global location in  adios (X*NY:(X*NY+NY-1), 0:9)
       ! BUF(2) is number of processors per block, split even
       ! BUF(3) is rows to read per processor, same across the processors
       ! ==========================  
       
       
       ! ===============
       ! box is identified here according to how bp is constructed in read.bp
       ! ==============
       if (rank > 0) then   
             boxOffset(1) = 0;       boxOffset(2) = (rank/RBUF(2))*NY+RBUF(3)*(mod(rank, RBUF(2)));
             boxSize(1)   = NX;       boxSize(2)   = RBUF(3);
       else  
             boxOffset(1) = 0;       boxOffset(2) = 0;
	     boxSize(1)   = NX;       boxSize(2)   = SBUF(3);
       endif
       
       ! ===================
       ! get data from adios
       ! ===================

       !allocate (pts(NX, boxSize(2)));
       this%Nptlocal = boxSize(2)
       deallocate(this%Pts1)
       allocate(this%Pts1(NX,this%Nptlocal))
       !print*,"Nptlocal: ",this%Nptlocal,NX
       !print*,"boxoff,size: ",boxOffset,boxSize

       call adios_selection_boundingbox(boxsel, 2, boxOffset, boxSize);
       call adios_schedule_read (f, boxSel, "particles", 0, 1, this%Pts1, ierr);
       call adios_perform_reads(f,ierr);

       call adios_selection_delete(boxSel);

    else  
       !
       ! not for this example
       !
    endif

      call adios_read_close (f, ierr)	
      call adios_read_finalize_method (method, ierr);	

      call MPI_Barrier(comm, ierr);

end subroutine phaseinadios

        !the 6D phase space output has (x(m), px/mc, y(m), py/mc, z(m), pz/mc).
        subroutine phaseold_Output(nfile,this,samplePeriod)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nfile
        type (BeamBunch), intent(in) :: this
        integer :: samplePeriod
        integer :: np,my_rank,ierr
        integer status(MPI_STATUS_SIZE)
        integer :: i,j,sixnpt,mnpt
        integer, allocatable, dimension(:) :: nptlist
        double precision, allocatable,dimension(:,:) :: recvbuf

        if (samplePeriod .eq. 0) then
           samplePeriod = 1
        endif

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
        call MPI_ALLREDUCE(this%Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
                        MPI_COMM_WORLD,ierr)

        allocate(nptlist(0:np-1))
        nptlist = 0
        allocate(recvbuf(6,mnpt))
        sixnpt = 6*this%Nptlocal

        call MPI_GATHER(this%Nptlocal,1,MPI_INTEGER,nptlist,1,&
                        MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        nptlist = 6*nptlist

        if(my_rank.eq.0) then
!hjw
! All the following I/O commented out by hjw
!hjw
!         write(*,*)' Opening file: ',nfile,' in phase_Output'
         open(nfile,status='unknown')
          do i = 1, this%Nptlocal,samplePeriod
            write(nfile,100)this%Pts1(1,i),this%Pts1(2,i),this%Pts1(3,i),&
                            this%Pts1(4,i),this%Pts1(5,i),this%Pts1(6,i)
!           write(nfile,100)this%Pts1(1,i)*Scxlt,this%Pts1(2,i),&
!                           this%Pts1(3,i)*Scxlt,&
!                           this%Pts1(4,i),this%Pts1(5,i)*Scxlt,&
!                           this%Pts1(6,i)
          enddo
          do i = 1, np-1
            call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
                          i,1,MPI_COMM_WORLD,status,ierr) 
        
            do j = 1, nptlist(i)/6,samplePeriod
              write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                              recvbuf(4,j),recvbuf(5,j),recvbuf(6,j)
!             write(nfile,100)recvbuf(1,j)*Scxlt,recvbuf(2,j),&
!                             recvbuf(3,j)*Scxlt,&
!                             recvbuf(4,j),recvbuf(5,j)*Scxlt,recvbuf(6,j)
            enddo
          enddo
         close(nfile)
        else
          call MPI_SEND(this%Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
                        MPI_COMM_WORLD,ierr)
        endif

!100     format(6(1x,e17.9))
100     format(6(1x,e20.12))

        deallocate(nptlist)
        deallocate(recvbuf)

        end subroutine phaseold_Output
      end module Outputclass

