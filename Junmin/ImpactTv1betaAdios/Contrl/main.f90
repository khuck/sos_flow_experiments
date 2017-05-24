!----------------------------------------------------------------
! (c) Copyright, 2015 by the Regents of the University of California.
! IMPACT-T: A 3D parallel particle tracking code in time domain
! Version: 1.0-beta
! Author: Ji Qiang
! Description: 
!            IMPACT-T is a fully three-dimensional program to track relativistic 
!	     particles taking into account space charge forces and short-range 
!            longitudinal and transverse wakefields. IMPACT-T is one of the few 
!            codes used in the photoinjector community that has a parallel implementation,
!            making it very useful for high statistics simulations of beam halos and 
!            beam diagnostics. It has a comprehensive set of beamline elements, and 
!            furthermore allows arbitrary overlap of their fields, which gives the 
!            IMPACT-T a capability to model both the standing wave structure and traveling
!            wave structure. It is also unique in its use of space-charge solvers based 
!            on an integrated Green function to efficiently and accurately treat beams 
!            with large aspect ratio, and a shifted Green function to efficiently treat 
!            image charge effects of a cathode. It is also unique in its inclusion of 
!            energy binning in the space-charge calculation to model beams with large 
!            energy spread. IMPACT-T has a flexible data structure that allows particles 
!            to be stored in containers with common characteristics; for photoinjector 
!            simulations the containers represent multiple slices, but in other 
!            applications they could correspond, e.g., to particles of different species. 
!            Together, all these features make IMPACT-T a powerful and versatile tool for 
!            modeling beams in photoinjectors and other systems.
!----------------------------------------------------------------

      program main
      use AccSimulatorclass
      implicit none
      include 'mpif.h'
      double precision :: time
      integer :: ierr

      call MPI_INIT(ierr)

      call construct_AccSimulator(time)
      call run_AccSimulator()
      call destruct_AccSimulator(time)

      call MPI_Finalize(ierr)


      end program main
