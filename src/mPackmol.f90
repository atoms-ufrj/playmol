!   This file is part of Playmol.
!
!    Playmol is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Playmol is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Playmol. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu at eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module mPackmol

use mGlobal
use mMolecule
use mStruc
use mBox

implicit none

character(sl), parameter :: stdout_name = "/proc/self/fd/1"
character(sl), parameter :: logfile = "packmol.log"

real(rb), parameter, private :: rbtol = 1.e-6_rb

type tPackmol
  integer  :: seed = 1234
  integer  :: nloops = 50
  real(rb) :: retry = 1.0_rb
  real(rb) :: diameter = 2.5_rb
  logical  :: setup
  type(StrucList) :: list = StrucList( "packmol option" )
  contains
    procedure :: run => tPackmol_run
    procedure :: file_names => tPackmol_file_names
    procedure :: total_mass => tPackmol_total_mass
end type tPackmol

interface
  subroutine packmol( inp, stat )
    integer, intent(in)  :: inp
    integer, intent(out) :: stat
  end subroutine packmol 
end interface

contains

  !=================================================================================================

  subroutine tPackmol_file_names( me, inputfile, molfile, mixfile, random_names )
    class(tPackmol), intent(in) :: me
    character(sl),  intent(out) :: inputfile, molfile(:), mixfile
    logical,        intent(in)  :: random_names
    integer :: i
    type(rng) :: random
    inputfile = "packmol.inp"
    do i = 1, size(molfile)
      molfile(i) = "molecule-"//trim(int2str(i))//".xyz"
    end do
    mixfile = "packmol-output.xyz"
    if (random_names) then
      call random % init( me % seed )
      inputfile = "00"//trim(random % letters(6))//"_"//trim(inputfile)
      do i = 1, size(molfile)
        molfile(i) = "00"//trim(random % letters(6))//"_"//trim(molfile(i))
      end do
      mixfile = "00"//trim(random % letters(6))//"_"//trim(mixfile)
    end if
  end subroutine tPackmol_file_names

  !=================================================================================================

  function tPackmol_total_mass( me, molmass, mol ) result( mass )
    class(tPackmol), intent(in) :: me
    real(rb),        intent(in) :: molmass(:)
    type(tMolecule), intent(in) :: mol
    real(rb)                    :: mass

    type(Struc), pointer :: ptr
    integer :: imol, narg, nmol
    character(sl) :: arg(4)

    mass = 0.0_rb
    ptr => me % list % first
    do while (associated(ptr))
      call split( ptr % params, narg, arg )
      imol = mol % index( arg(1) )
      select case (ptr % id(1))
        case ("move","fix")
          mass = mass + molmass(imol)
        case ("copy","pack")
          nmol = str2int(arg(2))
          mass = mass + nmol*molmass(imol)
      end select
      ptr => ptr % next
    end do

  end function tPackmol_total_mass

  !=================================================================================================

  subroutine tPackmol_run( me, mol, lcoord, latom, ldiam, Lbox )
    use iso_fortran_env, only : screen => output_unit
    class(tPackmol), intent(inout) :: me
    type(tMolecule), intent(inout) :: mol
    type(StrucList), intent(inout) :: lcoord, latom, ldiam
    real(rb),        intent(in)    :: Lbox(3)

    integer       :: stat, unit, imol, ifile, iatom, nfiles, file(mol%N), narg
    real(rb)      :: scaling
    logical       :: redirect
    character(sl) :: mixfile, inputfile, ctype, arg(1)

    type(Struc),   pointer     :: ptr, pmol
    integer,       allocatable :: aux(:), molecule(:), natoms(:), first_atom(:)
    real(rb),      allocatable :: diameter(:,:)
    character(sl), allocatable :: molfile(:)

    if (.not.associated(me % list % first)) call error( "no molecules chosen for packing" )

    ! Determine the files associated with each molecule and vice-versa:
    nfiles = 0
    file = 0
    allocate( aux(me % list % count) )
    ptr => me % list % first
    do while (associated(ptr))
      call split( ptr % params, narg, arg )
      imol = mol % index( arg(1) )
      if (all(aux(1:nfiles) /= imol)) then
        nfiles = nfiles + 1
        aux(nfiles) = imol
        file(imol) = nfiles
      end if
      ptr => ptr % next
    end do
    allocate( molecule(nfiles) )
    molecule = aux(1:nfiles)
    deallocate( aux )

    ! Find out the number of atoms of each considered molecule:
    allocate( natoms(nfiles) )
    natoms = 0
    ptr => mol % list % first
    do while (associated(ptr))
      ifile = file(str2int(ptr%params))
      if (ifile > 0) natoms(ifile) = natoms(ifile) + 1
      ptr => ptr % next
    end do

    ! Determine the first atom of each molecule in the coordinate list:
    allocate( first_atom(nfiles) )
    ptr => lcoord % first
    first_atom = 0
    iatom = 0
    do while (associated(ptr) .and. any(first_atom == 0))
      iatom = iatom + 1
      call mol % list % search( ptr % id, pmol )
      ifile = file(str2int(pmol % params))
      if (ifile > 0) then
        if (first_atom(ifile) == 0) first_atom(ifile) = iatom
      end if
      ptr => ptr % next
    end do
    if (any(first_atom == 0)) call error( "all molecules must have predefined coordinates" )

    ! Retrieve file names:
    allocate( molfile(nfiles) )
    call me % file_names( inputfile, molfile, mixfile, random_names = .not.me%setup )

    ! Write down xyz files and determine atom diameters:
    allocate( diameter(nfiles,maxval(natoms)) )
    diameter = me%diameter
    do ifile = 1, nfiles
      if (me%setup) call writeln( "Saving coordinate file", molfile(ifile) )
      open( newunit = unit, file = molfile(ifile), status = "replace" )
      write(unit,'(A)') trim(int2str(natoms(ifile)))
      write(unit,'("# Generated by Playmol on ",A)') trim(now())
      ptr => lcoord % first
      do iatom = 2, first_atom(ifile)
        ptr => ptr % next
      end do
      do iatom = 1, natoms(ifile)
        write(unit,'(A)') trim(ptr % id(1))//" "//trim(ptr % params)
        ctype = latom % parameters( ptr % id(1:1) )
        call ldiam % search( [ctype], pmol )
        if (associated(pmol)) diameter(ifile,iatom) = str2real(pmol % params)
        ptr => ptr % next
      end do
      close(unit)
    end do

    ! Try and redirect the standart output:
    open( newunit = unit, file = stdout_name, status = "old", iostat = stat )
    redirect = stat == 0
    if (redirect) close( unit )

    ! Run playmol with specified parameters:
    scaling = 1.0_rb
    if (me%setup) then

      call writeln( "Saving packmol input file", inputfile )
      call save_input_file

    else if (abs(me%retry - 1.0_rb) < rbtol) then

      call execute_packmol( stat )
      if (stat == 0) then
        call retrieve_coordinates( mixfile )
        call delete_files( [inputfile, mixfile, molfile, logfile] )
        call writeln( "Packmol converged with scaling factor", real2str(scaling) )
      else
        call delete_files( [inputfile, mixfile, trim(mixfile)//"_FORCED", molfile] )
        call error( "Packmol did not converge! See file packmol.log for additional info.")
      end if

    else

      stat = 1
      do while (stat /= 0)
        call execute_packmol( stat )
        if (stat /= 0) scaling = me%retry * scaling
      end do
      call retrieve_coordinates( mixfile )
      call delete_files( [inputfile, mixfile, trim(mixfile)//"_FORCED", molfile, logfile] )
      if (scaling == 1.0_rb) then
        call writeln( "Packmol converged with scaling factor", real2str(scaling) )
      else
        call warning( "Packmol converged with scaling factor", real2str(scaling) )
      end if

    end if

    contains
      !---------------------------------------------------------------------------------------------
      subroutine execute_packmol( stat )
        integer,  intent(out) :: stat
        integer :: inp
        call save_input_file
        open( newunit = inp, file = inputfile, status = "old" )
        call writeln( "Executing Packmol with scaling factor", trim(real2str(scaling))//".", &
                      "Please wait..." )
        if (redirect) open( unit = screen , file = "packmol.log", status = "replace" )
        call packmol( inp, stat )
        if (redirect) then
          close( screen )
          open( unit = screen , file = stdout_name, status = "old" )
        end if
        close(inp)
      end subroutine execute_packmol
      !---------------------------------------------------------------------------------------------
      subroutine save_input_file
        integer :: unit, imol, ifile, iatom, narg
        real(rb) :: D
        character(sl) :: box_limits, arg(4)
        type(Struc), pointer :: ptr
        open( newunit = unit, file = inputfile, status = "replace" )
        write(unit,'("# Generated by Playmol on ",A,/)') trim(now())
        write(unit,'("tolerance ",A)') trim(real2str(scaling*me%diameter))
        write(unit,'("filetype xyz")')
        write(unit,'("seed ",A)') trim(int2str(me%seed))
        write(unit,'("nloop ",A)') trim(int2str(me%nloops))
        write(unit,'("output ",A)') trim(mixfile)
        ptr => me % list % first
        do while (associated(ptr))
          call split( ptr % params, narg, arg )
          imol = mol % index( arg(1) )
          ifile = file(imol)
          write(unit,'(/,"structure ",A)') trim(molfile(ifile))
          select case (ptr % id(1))
          case ("move","fix")
            write(unit,'("  number 1")')
            if (ptr % id(1) == "fix") write(unit,'("  center")')
            write(unit,'("  fixed ",A," 0.0 0.0 0.0")') trim(join(arg(2:4)))
          case ("copy","pack")
            write(unit,'("  number ",A)') trim(arg(2))
            D = scaling*maxval(diameter(ifile,1:natoms(ifile)))
            box_limits = join(real2str([-0.5_rb*(Lbox - D), +0.5_rb*(Lbox - D)]))
            write(unit,'("  inside box ",A)') trim(box_limits)
            if (ptr % id(1) == "copy") then
              write(unit,'("  constrain_rotation x 0.0 0.0")')
              write(unit,'("  constrain_rotation y 0.0 0.0")')
              write(unit,'("  constrain_rotation z 0.0 0.0")')
            end if
          end select
          do iatom = 1, natoms(ifile)
            D = diameter(ifile,iatom)
            if (abs((D - me%diameter)/me%diameter) > rbtol) then
              write(unit,'("  atoms ",A)') trim(int2str(iatom))
              write(unit,'("    radius ",A)') trim(real2str(0.5_rb*scaling*D))
              write(unit,'("  end atoms ")')
            end if
          end do
          write(unit,'("end structure")')
          ptr => ptr % next
        end do
        close(unit)
      end subroutine save_input_file
      !---------------------------------------------------------------------------------------------
      subroutine retrieve_coordinates( mixfile )
        character(*), intent(in) :: mixfile
        integer :: mix, narg, natoms, iatom
        character(sl) :: line, arg(4)
        call lcoord % destroy
        open( newunit = mix, file = mixfile, status = "old" )
        read(mix,'(I12,/)') natoms
        do iatom = 1, natoms
          read(mix,'(A'//csl//')') line
          call split( line, narg, arg )
          call lcoord % add( narg, arg, repeatable = .true., silent = .true. )
        end do
        close( mix, status = "delete" )
      end subroutine retrieve_coordinates
      !---------------------------------------------------------------------------------------------
  end subroutine tPackmol_run

  !=================================================================================================

end module mPackmol
