!----------------------------------------------------------------
! (c) Copyright, 2004 by the Regents of the University of California.
! BeamLineElemclass: Beam line element base class in Lattice module of 
!                    APPLICATION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL, 1/10/04
! Description: This class defines the base beam line element class for
!              different lattice element class.
! Comments:
!----------------------------------------------------------------
      module BeamLineElemclass
        use DriftTubeclass
        use BPMclass
        type BeamLineElem
!          private
          type (BPM), pointer :: pbpm
          type (DriftTube), pointer :: pdrift
        end type BeamLineElem
        interface assign_BeamLineElem
          module procedure assign_drift,assign_bpm
        end interface
        interface getparam_BeamLineElem
          module procedure getparam1_BeamLineElem, &
                           getparam2_BeamLineElem, &
                           getparam3_BeamLineElem
        end interface
        interface setparam_BeamLineElem
          module procedure setparam1_BeamLineElem, &
                           setparam2_BeamLineElem, &
                           setparam3_BeamLineElem
        end interface
      contains
        function assign_drift(tdrift) result(ppdrift)
        type (BeamLineElem) :: ppdrift
        type (DriftTube), target, intent(in) :: tdrift

        ppdrift%pdrift => tdrift
        nullify(ppdrift%pbpm)

        end function assign_drift

       function assign_bpm(tbpm) result(ppbpm)
        type (BeamLineElem) :: ppbpm
        type (BPM), target, intent(in) :: tbpm
 
        ppbpm%pbpm => tbpm
        nullify(ppbpm%pdrift)
 
        end function assign_bpm

        subroutine getparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,i,blparam)
        endif

        end subroutine getparam1_BeamLineElem
  
        subroutine getparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,blparams)
        endif

        end subroutine getparam2_BeamLineElem

        subroutine getparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,blength,bnseg,bmapstp,&
                                  btype)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,blength,bnseg,bmapstp,btype)
        endif

        end subroutine getparam3_BeamLineElem
       
        subroutine getradius_BeamLineElem(this,piperadius)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: piperadius

        if(associated(this%pdrift)) then
          call getparam_DriftTube(this%pdrift,2,piperadius)
        elseif(associated(this%pbpm)) then
          call getparam_BPM(this%pbpm,2,piperadius)
        endif

        end subroutine getradius_BeamLineElem
  
        subroutine geterr_BeamLineElem(this,xerr,yerr,anglerrx,anglerry,&
                                       anglerrz)
        implicit none 
        type (BeamLineElem), intent(in) :: this
        double precision, intent(out) :: xerr,yerr,anglerrx,anglerry,anglerrz

        if(associated(this%pdrift)) then
          xerr = 0.0
          yerr = 0.0
          anglerrx = 0.0
          anglerry = 0.0
          anglerrz = 0.0
        elseif(associated(this%pbpm)) then
          xerr = 0.0
          yerr = 0.0
          anglerrx = 0.0
          anglerry = 0.0
          anglerrz = 0.0
        endif

        end subroutine geterr_BeamLineElem
  
        subroutine setparam1_BeamLineElem(this,i,blparam)
        implicit none 
        type (BeamLineElem), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: blparam

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,i,blparam)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,i,blparam)
        endif

        end subroutine setparam1_BeamLineElem
  
        subroutine setparam2_BeamLineElem(this,blparams)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, dimension(:), intent(in) :: blparams

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,blparams)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,blparams)
        endif

        end subroutine setparam2_BeamLineElem

        subroutine setparam3_BeamLineElem(this,blength,bnseg,bmapstp,&
                                          btype)
        implicit none
        type (BeamLineElem), intent(inout) :: this
        double precision, intent(in) :: blength
        integer, intent(in) :: bnseg,bmapstp,btype

        if(associated(this%pdrift)) then
          call setparam_DriftTube(this%pdrift,bnseg,bmapstp,&
                                  btype,blength)
        elseif(associated(this%pbpm)) then
          call setparam_BPM(this%pbpm,bnseg,bmapstp,btype,blength)
        endif

        end subroutine setparam3_BeamLineElem
       
        subroutine getfld_BeamLineElem(this,pos,extfld)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld

        if(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        endif

        end subroutine getfld_BeamLineElem

        !get external field with displacement and rotation errors.
        subroutine getflderr_BeamLineElem(this,pos,extfld,dx,dy,anglex,&
                                          angley,anglez)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(in) :: dx,dy,anglex,angley,anglez
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld

        if(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        endif

        end subroutine getflderr_BeamLineElem

        subroutine getaxfldE_BeamLineElem(this,z,ez1,ezp1,ezpp1)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, intent(in) :: z
        double precision, intent(out) :: ez1,ezp1,ezpp1

        if(associated(this%pdrift)) then
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        elseif(associated(this%pbpm)) then
          print*,"no field in BPM!!"
          ez1 = 0.0
          ezp1 = 0.0
          ezpp1 = 0.0
        endif

        end subroutine getaxfldE_BeamLineElem

        subroutine getfldt_BeamLineElem(this,pos,extfld,fldata)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld
        type (fielddata), intent(in) :: fldata

        if(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        endif

        end subroutine getfldt_BeamLineElem

        subroutine getflderrt_BeamLineElem(this,pos,extfld,fldata)
        implicit none
        type (BeamLineElem), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(6), intent(out) :: extfld
        type (fielddata), intent(in) :: fldata

        if(associated(this%pdrift)) then
          call getfld_DriftTube(pos,extfld,this%pdrift)
        elseif(associated(this%pbpm)) then
          !call getfld_BPM(pos,extfld,this%pbpm)
          print*,"no field for BPM!!"
          extfld = 0.0
        endif

        end subroutine getflderrt_BeamLineElem
      end module BeamLineElemclass
