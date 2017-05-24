module readutil
contains 

SUBROUTINE utilAdiosIn(bpFileName, comm, rank, numproc, pts, nptlc, npt)
    use adios_read_mod
    character(len=20), intent(IN):: bpFileName  
    integer,  intent(IN):: comm
    integer:: ierr, adios_err, rank, np;
    double precision, dimension(:,:), allocatable, INTENT(OUT):: pts
    integer                 :: method = ADIOS_READ_METHOD_BP
    integer*8               :: f, sel
    integer                 :: GlobalX, GlobalY, start, count, avg;
    integer*8, dimension(2) :: boxOffset=0, boxSize=1
    integer, intent(OUT):: nptlc, npt

    call adios_read_init_method (method, comm, "verbose=3", ierr);
    call adios_read_open (f, bpFilename, method, comm, ADIOS_LOCKMODE_NONE, 1.0, ierr);

    if (f == 0) then       
       write (*,*), "NO such file. ", bpFileName
    else
       !call adios_selection_writeblock (sel, rank)
       call adios_selection_writeblock (sel, 0)
       call adios_schedule_read (f, sel, "GX", 0, 1, GlobalX, ierr)
       call adios_schedule_read (f, sel, "GY", 0, 1, GlobalY, ierr)
       call adios_perform_reads (f, ierr)
       
       !!write (*,*), "phase in: rank=", rank, "total proc= ", numproc, " file=", bpFileName, "Gx/Gy=", GlobalX, GlobalY, f;
       
       start = 0;
       count = GlobalY;
       
       npt = GlobalY
       
       avg = GlobalY/numproc;
       
       if (avg == 0) then 
          allocate (pts(GlobalX, GlobalY));
          sel = 0;
          call adios_schedule_read(f, sel, "particles", 0, 1, pts, ierr);
       else 
          start = rank * avg;
          if (rank < numproc-1) then
             count = avg;
          else 
             count = GlobalY - start;
          endif
          boxOffset(1) = 0; boxOffset(2) = start;
          boxSize(1) = GlobalX;  boxSize(2) = count;
          allocate (pts(GlobalX, count));
          call adios_selection_boundingbox(sel, 2, boxOffset, boxSize);
          call adios_schedule_read(f, sel, "particles", 0, 1, pts, ierr);       
       endif
       call adios_perform_reads(f, ierr);
       
       !write(*,*) "!rank= s", rank, "start/count", start, count,  pts(1,1), pts(1,2), pts(2,1);
       call adios_read_close(f, ierr);
       
       nptlc = count
    endif
    RETURN
 END 
    
SUBROUTINE utilAdiosOut(bpFileName, comm,  pts, total_particles_from_all, num_attr_per_particle, npts, particle_start, timestep, ts)
    character(len=20), intent(IN):: bpFileName
    integer,  intent(IN):: comm
    integer:: ierr, adios_err, rank
    integer*8            :: m_adios_group	
    integer*8            :: adios_handle
    integer*8            :: adios_groupsize, adios_totalsize

    integer, INTENT(IN):: timestep, ts; !ts=totalstep

    !!real*8, dimension(:,:), allocatable, INTENT(IN):: pts	
    real*8, dimension(:,:), INTENT(IN):: pts	
    integer, INTENT(IN):: total_particles_from_all, num_attr_per_particle, npts, particle_start

    integer NX, NY;
    integer :: GlobalX, GlobalY, OffsetX, OffsetY;
    integer estimateMB;

    GlobalX =   num_attr_per_particle;	
    GlobalY = total_particles_from_all;
    OffsetX = 0;
    OffsetY = particle_start;

    NY = npts
    NX = num_attr_per_particle

    call    MPI_comm_rank(comm, rank, ierr);

    !write(*,*) " oh we will write to: [", bpFileName, "]"; 
    write(*,*) "rank=", rank, " on NX/NY= ", NX, NY, NX*NY*8/1000, "GlobalX/Y=", GlobalX, GlobalY, "OffsetX/Y=", OffsetX, OffsetY;
    write(*,*) "rank=", rank, " list 1: ", pts(1:6,1:1)

! ===============================================
! the next four lines prepares adios file writing
! allocates buffer etc. 
! ===============================================
    call adios_init_noxml (comm, adios_err)
    !call adios_allocate_buffer (10, adios_err)
    estimateMB = (NX*NY)/10000;
    estimateMB = esimateMB+5;
    !write(*,*) " oh NX/NY= ", NX, NY, estimateMB;
    !call adios_set_max_buffer_size(estimateMB, adios_err)


! ===============================================
! file is opened for writing
! ===============================================
    if (timestep == 0) then       
       call adios_declare_group (m_adios_group, "restart", "iter", 1, adios_err)
       call adios_select_method (m_adios_group, "FLEXPATH", "", "", adios_err)
       call utilAdiosDefineFile(m_adios_group);

       !call adios_select_method (m_adios_group, "MPI", "", "", adios_err)
       call adios_open (adios_handle, "restart", bpFilename, "w", comm, adios_err)

    else
       !call adios_open (adios_handle, "restart", bpFilename, "a", comm, adios_err)
       call adios_open (adios_handle, "restart", bpFilename, "a", comm, adios_err)

    endif

! ===============================================
! total size is assigned
! ===============================================
    !adios_groupsize = ffp * (4 * 6 + NX*base * 8)
    adios_groupsize = NX * NY * 8 + 200
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)   
    !write (*,*), "oh groupsize=", adios_groupsize, adios_err


    call utilAdiosWriteFile(adios_handle, NX, NY, GlobalX, GlobalY, OffsetX, OffsetY, pts)

!   deallocate(pts);	    

    if (ts == timestep) then
       call adios_finalize (rank, adios_err)
    endif

RETURN
END

SUBROUTINE utilAdiosWriteFile(adios_handle, NX, NY, GlobalX, GlobalY, OffsetX, OffsetY, pts)
      integer*8 , intent(IN)             :: adios_handle
      integer, intent(IN)  ::  NX, NY
      integer, intent(IN)                 :: OffsetX, OffsetY, GlobalX, GlobalY
      !!real*8, dimension(:,:), allocatable, INTENT(IN):: pts
      real*8, dimension(:,:), INTENT(IN):: pts

      integer                 :: adios_err


! ===============================================
! write to adios
! ===============================================

         call adios_write (adios_handle, "NX", NX, adios_err)
         call adios_write (adios_handle, "NY", NY, adios_err)

         call adios_write (adios_handle, "GX", GlobalX, adios_err)
         call adios_write (adios_handle, "GY", GlobalY, adios_err)

         call adios_write (adios_handle, "OX", OffsetX, adios_err)
         call adios_write (adios_handle, "OY", OffsetY, adios_err)
         call adios_write (adios_handle, "particles", pts, adios_err)

    call adios_close (adios_handle, adios_err)

RETURN
END

SUBROUTINE utilAdiosDefineFile(m_adios_group)
      integer*8  , intent(IN)            :: m_adios_group
      integer*8               :: varid
! ===============================================
! Defines variables in adios. Later the variables will be assigned values by each processor
! GX,GY represent gloabl dimention (of all the data)
! NX,NY represent block size (data size from each processor.)
! OX,OY is the offset of block in the global dimention
! e.g. if you have a 2x2 array, and 4 processors, each processor writes 1 element.
! then GX=2, GY=2. NX=1, NY=1 for all processors. (OX,OY) = (0,0) (0,1) (1,0), (1,1) depends on
! which processor. Notice the indices are from 0, instead of 1.
! ===============================================

        ! This example doesn't use varid during writing.
        ! So we simply put 'varid' everywhere.
        ! define a integer
        call adios_define_var (m_adios_group, "NX" ,"", 2 ,"", "", "", varid)
        call adios_define_var (m_adios_group, "NY" ,"", 2 ,"", "", "", varid)

	! define a integer
        call adios_define_var (m_adios_group, "GX" ,"", 2 ,"", "", "", varid)
        call adios_define_var (m_adios_group, "GY" ,"", 2 ,"", "", "", varid)

	! define a integer
        call adios_define_var (m_adios_group, "OX" ,"", 2 ,"", "", "", varid)
        call adios_define_var (m_adios_group, "OY" ,"", 2 ,"", "", "", varid)

	! define a global array
        call adios_define_var (m_adios_group, "particles" ,"", 6 ,"NX,NY", "GX,GY", "OX,OY", varid)

RETURN
END

end module






