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

! TO DO:
! 1) unlink command: warn if molecule splitting make extra bonds/angles/dihedrals have atoms in
!    separate molecules
! 2) Allow only one asterisk to be used as wildcard.

module mPlaymol

use mGlobal
use mStruc
use mPackmol
use mMolecule
use mBox
use mCodeFlow
use mFix

implicit none

type tPlaymol
  type(tBox)      :: box
  type(tPackmol)  :: packmol
  type(tCodeFlow) :: code_flow
  type(tMolecule) :: molecules
  type(tFix)      :: atomfix, &
                     typefix

  type(StrucList) :: atom_type_list      = StrucList( "atom type" )
  type(StrucList) :: bond_type_list      = StrucList( "bond type", 2, .true. )
  type(StrucList) :: angle_type_list     = StrucList( "angle type", 3, .true. )
  type(StrucList) :: dihedral_type_list  = StrucList( "dihedral type", 4, .true. )
  type(StrucList) :: improper_type_list  = StrucList( "improper type", 4 )
  type(StrucList) :: mass_list           = StrucList( "mass" )
  type(StrucList) :: atom_masses         = StrucList( "atom mass" )
  type(StrucList) :: diameter_list       = StrucList( "diameter" )
  type(StrucList) :: atom_list           = StrucList( "atom" )
  type(StrucList) :: charge_list         = StrucList( "charge" )
  type(StrucList) :: bond_list           = StrucList( "bond", 2, .true. )
  type(StrucList) :: link_list           = StrucList( "virtual link", 2, .true. )
  type(StrucList) :: angle_list          = StrucList( "angle", 3, .true. )
  type(StrucList) :: dihedral_list       = StrucList( "dihedral", 4, .true. )
  type(StrucList) :: extra_bond_list     = StrucList( "extra bond", 2, .true. )
  type(StrucList) :: extra_angle_list    = StrucList( "extra angle", 3, .true. )
  type(StrucList) :: extra_dihedral_list = StrucList( "extra dihedral", 4, .true. )
  type(StrucList) :: improper_list       = StrucList( "improper", 4 )
  contains
    procedure :: read => tPlaymol_Read
    procedure :: write => tPlaymol_write
    procedure :: write_lammps => tPlaymol_write_lammps
    procedure :: write_lammpstrj => tPlaymol_write_lammpstrj
    procedure :: read_geometry => tPlaymol_read_geometry
    procedure :: write_xyz => tPlaymol_write_xyz
    procedure :: summarize => tPlaymol_summarize
    procedure :: update_structure => tPlaymol_update_structure
    procedure :: search_impropers => tPlaymol_search_impropers
    procedure :: get_types => tPlaymol_get_types
    procedure :: check_coordinates => tPlaymol_check_coordinates
end type tPlaymol

contains

  !=================================================================================================

  recursive subroutine tPlaymol_Read( me, unit, filename )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    character(*),    intent(in)    :: filename
    integer       :: narg
    character(sl) :: arg(50)
    call me % code_flow % next_command( unit, narg, arg )
    do while (narg > 0)
      select case (trim(arg(1)))
        case ("prefix","suffix"); call prefix_suffix_command
        case ("box"); call box_command
        case ("atom_type"); call atom_type_command
        case ("bond_type"); call bond_type_command
        case ("angle_type"); call angle_type_command
        case ("dihedral_type"); call dihedral_type_command
        case ("improper_type"); call improper_type_command
        case ("mass"); call mass_command
        case ("diameter"); call diameter_command
        case ("atom"); call atom_command
        case ("charge"); call charge_command
        case ("bond"); call bond_command
        case ("link"); call link_command
        case ("unlink"); call unlink_command
        case ("extra"); call extra_command
        case ("improper"); call improper_command
        case ("xyz","build"); call build_command
        case ("align"); call align_command
        case ("write"); call write_command
        case ("include"); call include_command
        case ("packmol"); call packmol_command
        case ("reset"); call reset_command
        case ("shell"); call shell_command
        case ("quit")
          if ((narg == 2).and.(arg(2) == "all")) then
            call writeln( "Script", filename, "interrupted by a 'quit all' command" )
            stop
          else if (narg == 1) then
            call writeln( "Script", filename, "interrupted by a 'quit' command" )
            return
          else
            call error( "invalid quit command" )
          end if
        case default; call error( "unknown command", arg(1) )
      end select
      call me % code_flow % next_command( unit, narg, arg )
    end do
    contains
      !---------------------------------------------------------------------------------------------
      subroutine prefix_suffix_command
        if (narg < 3) call error( "invalid", arg(1), "command" )
        call writeln( "Defining new", arg(1), "for ", arg(2), "as", arg(3) )
        select case (arg(2))
          case ("types"); call me % typefix % change( string = arg(3), kind = arg(1) )
          case ("atoms"); call me % atomfix % change( string = arg(3), kind = arg(1) )
          case default; call error( arg(1), "keyword must be 'types' or 'atoms'" )
        end select
        if (has_macros(arg(3))) call error( "invalid", arg(1), arg(3) )
      end subroutine prefix_suffix_command
      !---------------------------------------------------------------------------------------------
      subroutine box_command
        integer :: i
        real(rb) :: scalar, vector(3)
        if (narg < 3) call error( "invalid box command" )
        select case (arg(2))
          case ("density","volume")
            call writeln( "Defining box", arg(2), "as", arg(3) )
            scalar = str2real(arg(3))
            if (narg > 3) then
              if ((narg /= 7).or.(arg(4) /= "aspect")) call error( "invalid box command" )
              call writeln( "Defining box aspect as", join(arg(5:7)) )
              vector = [(str2real(arg(4+i)),i=1,3)]
            else
              vector = 1.0_rb
            end if
          case ("lengths")
            if (narg == 3) then
              call writeln( "Defining box side lengths as", arg(3) )
              vector = str2real(arg(3))
            else if (narg == 5) then
              call writeln( "Defining box side lengths as", join(arg(3:5)) )
              vector = [(str2real(arg(2+i)),i=1,3)]
            else
              call error( "invalid box command" )
            end if
          case ("angles")
            if (narg == 3) then
              call writeln( "Defining box angles as", arg(3) )
              vector = str2real(arg(3))
            else if (narg == 5) then
              call writeln( "Defining box angles as", join(arg(3:5)) )
              vector = [(str2real(arg(2+i)),i=1,3)]
            else
              call error( "invalid box command" )
            end if
          case default
             call error( "invalid box command" )
        end select
        call me % box % define( arg(2), scalar, vector )
      end subroutine box_command
      !---------------------------------------------------------------------------------------------
      subroutine atom_type_command
        call me % typefix % apply( arg(2) )
        call me % atom_type_list % add( narg-1, arg(2:narg) )
        if (has_macros(arg(2))) call error( "invalid atom_type name" )
      end subroutine atom_type_command
      !---------------------------------------------------------------------------------------------
      subroutine bond_type_command
        call me % typefix % apply( arg(2:3) )
        call me % bond_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine bond_type_command
      !---------------------------------------------------------------------------------------------
      subroutine angle_type_command
        call me % typefix % apply( arg(2:4) )
        call me % angle_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine angle_type_command
      !---------------------------------------------------------------------------------------------
      subroutine dihedral_type_command
        call me % typefix % apply( arg(2:5) )
        call me % dihedral_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine dihedral_type_command
      !---------------------------------------------------------------------------------------------
      subroutine improper_type_command
        call me % typefix % apply( arg(2:5) )
        call me % improper_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine improper_type_command
      !---------------------------------------------------------------------------------------------
      subroutine mass_command
        if (narg /= 3) call error( "invalid mass command" )
        call me % typefix % apply( arg(2) )
        call me % mass_list % add( narg-1, arg(2:narg) )
        if (str2real(arg(3)) < 0.0_rb) call error( "invalid mass value" )
      end subroutine mass_command
      !---------------------------------------------------------------------------------------------
      subroutine diameter_command
        if (narg /= 3) call error( "invalid diameter command" )
        call me % typefix % apply( arg(2) )
        call me % diameter_list % add( narg-1, arg(2:narg) )
        if (str2real(arg(3)) < 0.0_rb) call error( "invalid diameter value" )
      end subroutine diameter_command
      !---------------------------------------------------------------------------------------------
      subroutine atom_command
        type(Struc), pointer :: ptr
        if ((narg < 3).or.(narg > 4)) call error( "invalid atom command" )
        call me % atomfix % apply( arg(2) )
        call me % typefix % apply( arg(3) )
        call me % atom_list % add( 2, arg(2:3) )
        if (has_macros(arg(2))) call error( "invalid atom name" )
        if (has_macros(arg(3))) call error( "invalid atom type name" )

        call me % atom_type_list % search( arg(3:3), ptr )
        if (.not.associated(ptr)) call error( "atom type", arg(3), "required, but not found")
        ptr % used = .true.

        call me % mass_list % search( arg(3:3), ptr )
        if (.not.associated(ptr)) call error( "mass of atom type",arg(3), "has not been defined" )
        ptr % used = .true.
        arg(3) = ptr % params
        call me % atom_masses % add( 2, arg(2:3) )

        call me % charge_list % search( arg(2:2), ptr )
        if (associated(ptr)) ptr % used = .true.
        if (narg == 4) then
          if (associated(ptr)) call error( "an applicable charge has already been defined" )
          arg(3) = arg(4)
          call me % charge_list % add( 2, arg(2:3) )
          ptr => me % charge_list % last
        end if

        call me % molecules % add_atom( arg(2) )
      end subroutine atom_command
      !---------------------------------------------------------------------------------------------
      subroutine charge_command
        if (narg /= 3) call error( "invalid charge command" )
        call me % atomfix % apply( arg(2) )
        call me % charge_list % add( narg-1, arg(2:3) )
        if (.not.is_real(arg(3))) call error( "invalid charge value" )
        if (me % atom_list % find( arg(2:2) )) me % charge_list % last % used = .true.
      end subroutine charge_command
      !---------------------------------------------------------------------------------------------
      subroutine add_bond( atom )
        character(sl), intent(inout) :: atom(2)
        call me % bond_list % add( 2, atom, me % atom_list )
        call me % bond_list % handle( atom, me % atom_list, me % bond_type_list, 2 )
        call me % molecules % fuse( atom )
        call me % update_structure()
      end subroutine add_bond
      !---------------------------------------------------------------------------------------------
      subroutine bond_command
        integer :: i
        character(sl) :: central
        if (narg < 3) call error( "invalid bond command" )
        call me % atomfix % apply( arg(2:narg) )
        central = arg(2)
        do i = 3, narg
          arg(2) = central
          arg(3) = arg(i)
          call add_bond( arg(2:3) )
        end do
      end subroutine bond_command
      !---------------------------------------------------------------------------------------------
      subroutine link_command
        if (narg /= 3) call error( "invalid link command" )
        call me % atomfix % apply( arg(2:3) )
        call me % link_list % add( 2, arg(2:3), me % atom_list )
        call me % molecules % fuse( arg(2:3) )
      end subroutine link_command
      !---------------------------------------------------------------------------------------------
      subroutine unlink_command
        if (narg /= 3) call error( "invalid unlink command" )
        call me % atomfix % apply( arg(2:3) )
        call me % link_list % remove( arg(2:3) )
        call me % molecules % split( arg(2:3) )
      end subroutine unlink_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_command
        character(sl) :: choice
        integer :: i, imol, jmol
        if (narg < 2) call error( "invalid extra command" )
        choice = arg(2)
        arg(2:narg-1) = arg(3:narg)
        narg = narg - 1
        select case (choice)
          case ("bond"); call extra_bond_command
          case ("angle"); call extra_angle_command
          case ("dihedral"); call extra_dihedral_command
          case default; call error( "invalid extra command" )
        end select
        imol = str2int(me % molecules % list % parameters( arg(2:2) ))
        do i = 2, narg
          jmol = str2int(me % molecules % list % parameters( arg(i:i) ))
          if (jmol /= imol) call error( "atoms", join(arg(2:narg)), "are not in the same molecule" )
        end do
      end subroutine extra_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_bond_command
        call me % atomfix % apply( arg(2:3) )
        call me % extra_bond_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 3) call error( "invalid extra bond command")
        call me % extra_bond_list % handle( arg(2:3), me%atom_list, me%bond_type_list, 2 )
      end subroutine extra_bond_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_angle_command
        call me % atomfix % apply( arg(2:4) )
        call me % extra_angle_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 4) call error( "invalid extra angle command")
        call me % extra_angle_list % handle( arg(2:4), me%atom_list, me%angle_type_list, 2 )
      end subroutine extra_angle_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_dihedral_command
        call me % atomfix % apply( arg(2:5) )
        call me % extra_dihedral_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 5) call error( "invalid extra dihedral command")
        call me % extra_dihedral_list % handle( arg(2:5), me%atom_list, me%dihedral_type_list, 2 )
      end subroutine extra_dihedral_command
      !---------------------------------------------------------------------------------------------
      subroutine improper_command
        integer :: i, imol, jmol
        if (narg == 2) then
          if (arg(2) /= "search") call error( "invalid improper command" )
          call me % search_impropers()
        else if (narg == 5) then
          call me % atomfix % apply( arg(2:5) )
          call me % improper_list % add( narg-1, arg(2:narg), me % atom_list )
          if (narg /= 5) call error( "invalid improper command")
          call me % improper_list % handle( arg(2:5), me%atom_list, me%improper_type_list, 2 )
          imol = str2int(me % molecules % list % parameters( arg(2:2) ))
          do i = 3, 5
            jmol = str2int(me % molecules % list % parameters( arg(i:i) ))
            if (jmol /= imol) call error( "atoms", join(arg(2:5)), "are not in the same molecule" )
          end do
        else
          call error( "invalid improper command" )
        end if
      end subroutine improper_command
      !---------------------------------------------------------------------------------------------
      subroutine build_command
        logical :: file_exists
        integer :: geo
        if (narg == 1) then
          call writeln("Reading geometric data...")
          call me % read_geometry( unit )
        else
          inquire( file = arg(2), exist = file_exists )
          if (.not.file_exists) call error( "file", arg(2), "not found" )
          open( newunit = geo, file = arg(2), status = "old" )
          call writeln("Reading geometric data from file", trim(arg(2))//"..." )
          call me % read_geometry( geo )
          close( geo )
        end if
      end subroutine build_command
      !---------------------------------------------------------------------------------------------
      subroutine align_command
        integer :: i, imol, natoms(me%molecules%N), N
        real(rb), allocatable :: Mass(:), Coord(:,:)
        integer :: axis(3)
        character :: dir(3) = ["x","y","z"]
        real(rb) :: lb(3), ub(3)
        if (narg /= 4) call error( "invalid align command" )
        imol = me % molecules % index(arg(2))
        do i = 3, 4
          if (len_trim(arg(i)) /= 1) call error( "invalid align command" )
          select case (trim(arg(i)))
            case ("x"); axis(i-2) = 1
            case ("y"); axis(i-2) = 2
            case ("z"); axis(i-2) = 3
            case default; call error( "invalid align command" )
          end select
        end do
        if (axis(1) == axis(2)) call error( "invalid align command" )
        axis(3) = 6 - axis(1) - axis(2)
        call me % check_coordinates( imol )
        call writeln( "Aligning molecule", arg(2), "with axes", arg(3), "and", arg(4) )
        natoms = me % molecules % number_of_atoms()
        N = natoms(imol)
        allocate( Mass(N), Coord(3,N) )
        Mass = 1.0_rb
        call me % molecules % coordinates( imol, N, Coord, me % molecules % xyz, 1 )
        lb = minval(Coord,dim=2)
        ub = maxval(Coord,dim=2)
        do i = 1, 3
          call writeln( "Bounds in direction", dir(i)//":", join(real2str([lb(i),ub(i)])) )
        end do
        call writeln( "Bounding box lengths: ", join(real2str(ub-lb)) )
        call me % molecules % align( N, Mass, Coord, axis )
        call me % molecules % coordinates( imol, N, Coord, me % molecules % xyz, 2 )
        lb = minval(Coord,dim=2)
        ub = maxval(Coord,dim=2)
        do i = 1, 3
          call writeln( "--> New bounds in direction", dir(i)//":", join(real2str([lb(i),ub(i)])) )
        end do
        call writeln( "--> New bounding box lengths: ", join(real2str(ub-lb)) )
        deallocate( Mass, Coord )
      end subroutine align_command
      !---------------------------------------------------------------------------------------------
      subroutine write_command
        integer :: unit
        if ((narg < 2).or.(narg > 3)) call error( "invalid write command" )
        if (.not.any(arg(2) == [character(9)::"playmol","lammps","summary","xyz","lammpstrj"]) ) then
          call error( "invalid write command" )
        end if
        if (narg == 3) then
          open( newunit = unit, file = arg(3), status = "replace" )
          call writeln( "Writing data to file", arg(3), "in", arg(2), "format..." )
        else
          unit = stdout
        end if
        select case (arg(2))
          case ("playmol"); call me % write( unit )
          case ("lammps");  call me % write_lammps( unit )
          case ("summary"); call me % summarize( unit )
          case ("xyz"); call me % write_xyz( unit )
          case ("lammpstrj"); call me % write_lammpstrj( unit )
        end select
        if (unit /= stdout) close(unit)
      end subroutine write_command
      !---------------------------------------------------------------------------------------------
      subroutine include_command
        logical :: file_exists
        integer :: input
        if (narg < 2) call error( "invalid include command" )
        inquire( file = arg(2), exist = file_exists )
        if (.not.file_exists) call error( "file", arg(2), "does not exist" )
        open( newunit = input, file = arg(2), status = "old" )
        call me % Read( input, trim(arg(2)) )
        close(input)
      end subroutine include_command
      !---------------------------------------------------------------------------------------------
      subroutine reset_lists( lists )
        integer, intent(in) :: lists(:)
        ! Lists directly defined by the user:
        if (any(lists ==  1)) call me % atom_type_list % destroy
        if (any(lists ==  2)) call me % mass_list % destroy
        if (any(lists ==  3)) call me % diameter_list % destroy
        if (any(lists ==  4)) call me % bond_type_list % destroy
        if (any(lists ==  5)) call me % angle_type_list % destroy
        if (any(lists ==  6)) call me % dihedral_type_list % destroy
        if (any(lists ==  7)) call me % improper_type_list % destroy
        if (any(lists ==  8)) call me % atom_list % destroy
        if (any(lists ==  9)) call me % charge_list % destroy
        if (any(lists == 10)) call me % bond_list % destroy
        if (any(lists == 11)) call me % angle_list % destroy
        if (any(lists == 12)) call me % dihedral_list % destroy
        if (any(lists == 13)) call me % extra_bond_list % destroy
        if (any(lists == 14)) call me % extra_angle_list % destroy
        if (any(lists == 15)) call me % extra_dihedral_list % destroy
        if (any(lists == 16)) call me % improper_list % destroy
        if (any(lists == 17)) call me % link_list % destroy
        if (any(lists == 18)) call me % molecules % xyz % destroy
        if (any(lists == 19)) call me % packmol % list % destroy
        ! Lists automatically defined by Playmol:
        if (any(lists == 20)) call me % molecules % bonds % destroy
        if (any(lists == 21)) call me % molecules % list % destroy
        if (any(lists == 22)) call me % atom_masses % destroy
      end subroutine reset_lists
      !---------------------------------------------------------------------------------------------
      subroutine reset_command
        integer :: i
        type(Struc), pointer :: ptr
        if (narg /= 2) call error( "invalid reset command" )
        select case (arg(2))
          case ("link")
            ptr => me % link_list % first
            do while (associated(ptr))
              call me % molecules % split( ptr%id )
              ptr => ptr % next
            end do
          case ("bond")
            ptr => me % molecules % list % first
            do i = 1, me % molecules % list % count
              ptr % params = int2str(i)
              ptr => ptr % next
            end do
        end select
        select case (arg(2))
          case ("all");            call reset_lists( [(i,i= 1,22)] )
          case ("atom");           call reset_lists( [(i,i= 8,22)] )
          case ("charge");         call reset_lists( [9] )
          case ("bond");           call reset_lists( [(i,i=10,17),20] )
          case ("angle");          call reset_lists( [11,14] )
          case ("dihedral");       call reset_lists( [12,15] )
          case ("extra_bond");     call reset_lists( [13] )
          case ("extra_angle");    call reset_lists( [14] )
          case ("extra_dihedral"); call reset_lists( [15] )
          case ("improper");       call reset_lists( [16] )
          case ("link");           call reset_lists( [17] )
          case ("xyz");            call reset_lists( [18] )
          case ("packmol");        call reset_lists( [19] )
          case default;            call error( "invalid reset command" )
        end select
      end subroutine reset_command
      !---------------------------------------------------------------------------------------------
      subroutine shell_command
        integer :: stat
        call writeln( "Executing shell command: ", join(arg(2:narg)) )
        call execute_command_line( join(arg(2:narg)), exitstat = stat )
        if (stat /= 0) call error( "unsuccessful shell command" )
      end subroutine shell_command
      !---------------------------------------------------------------------------------------------
      subroutine packmol_command
        integer :: iarg, imol, nopts, i, n
        real(rb) :: pos(3), mass(me%molecules%N)
        character(sl) :: action
        type(Struc), pointer :: ptr
        action = ""
        if (narg < 3) call error( "invalid packmol command" )
        iarg = 2
        do while (iarg < narg)
          select case (arg(iarg))

            case ("seed")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol seed to", arg(iarg+1) )
              me % packmol % seed = str2int(arg(iarg+1))
              iarg = iarg + 2

            case ("tolerance","diameter")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol default diameter to", arg(iarg+1) )
              me % packmol % diameter = str2real(arg(iarg+1))
              if (me % packmol % diameter < 0.0_rb) call error( "invalid parameter value" )
              iarg = iarg + 2

            case ("nloops")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol nloops parameter to", arg(iarg+1) )
              me % packmol % nloops = str2int(arg(iarg+1))
              if (me % packmol % nloops <= 0) call error( "invalid parameter value" )
              iarg = iarg + 2

            case ("retry")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol retry parameter to", arg(iarg+1) )
              me % packmol % retry = str2real(arg(iarg+1))
              if ((me % packmol % retry <= 0.0_rb) .or. &
                  (me % packmol % retry > 1.0_rb)) call error( "invalid parameter value" )
              iarg = iarg + 2

            case ("fix","move","copy","pack")
              select case (arg(iarg))
                case ("move","fix"); nopts = 3
                case ("copy","pack"); nopts = 1
                case default; call error( "invalid packmol command" )
              end select
              if (narg < iarg+nopts+1) call error( "invalid packmol command" )
              call me % packmol % list % add( nopts+2, arg(iarg:iarg+1+nopts), repeatable = .true. )
              imol = me % molecules % index( arg(iarg+1) )
              select case (arg(iarg+1))
                case ("move","fix")
                  pos = [(str2real(arg(iarg+1+i)),i=1,3)]
                case ("copy","pack")
                  n = str2int(arg(iarg+2))
                  if (n <= 0) call error( "invalid parameter value" )
              end select
              iarg = iarg + 2 + nopts

            case ("action")
              if ((action /= "").or.(narg < iarg+1)) call error( "invalid packmol command" )
              action = arg(iarg+1)
              if ((action /= "setup").and.(action /= "execute")) then
                call error( "invalid packmol command: unknown action", action )
              end if
              if (.not. me % box % exists()) then
                call error( "packmol action keyword requires previous box definition" )
              end if
              me % packmol % setup = action == "setup"
              iarg = iarg + 2

            case default
              call error( "invalid packmol command: unknown keyword", arg(iarg) )

          end select
        end do
        if (action /= "") then
          call writeln( "Packmol invoked with action <"//trim(action)//">" )
          ptr => me % packmol % list % first
          do while (associated(ptr))
            call split( ptr%params, narg, arg(1:1) )
            imol = me % molecules % index( arg(1) )
            call me % check_coordinates( imol )
            ptr => ptr % next
          end do
          mass = me % molecules % per_molecule( me % atom_masses )
          call me % box % compute( me % packmol % total_mass( mass, me % molecules ) )
          call writeln( "Box lengths are ", join(real2str(me % box % length)))
          call me % packmol % run( me % molecules, me % molecules % xyz, &
                                   me % atom_list, me % diameter_list,   &
                                   me % box % length )
        end if
      end subroutine packmol_command
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_Read

  !=================================================================================================

  subroutine tPlaymol_read_geometry( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: N, i, narg, natoms
    character(sl) :: arg(7)
    integer,       allocatable :: ndata(:)
    character(sl), allocatable :: data(:,:)
    call me % code_flow % next_command( unit, narg, arg )
    if (narg == 0) call error( "expected geometric data not provided" )
    N = str2int( arg(1) )
    call writeln( "Number =", int2str(N) )
    allocate( data(N,7), ndata(N) )
    do i = 1, N
      call me % code_flow % next_command( unit, narg, arg )
      call writeln( int2str(i), ":", join(arg(1:narg)) )
      if (narg == 0) call error( "expected geometric data not completed" )
      select case (narg)
        case (3); natoms = 2 ! Bond
        case (4); natoms = 1 ! Coordinates
        case (5); natoms = 3 ! Bond and angle
        case (7); natoms = 4 ! Bond, angle, and dihedral
        case default; natoms = 0
      end select
      call me % atomfix % apply( arg(1:natoms) )
      ndata(i) = narg
      data(i,1:narg) = arg(1:narg)
    end do
    call me % molecules % set_geometry( data, ndata )
  end subroutine tPlaymol_read_geometry

  !=================================================================================================

  subroutine tPlaymol_write_xyz( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    type(Struc), pointer :: ptr
    write(unit,'(A)') trim(int2str(me % molecules % xyz % count))
    write(unit,'("# Generated by Playmol on ",A)') trim(now())
    ptr => me % molecules % xyz % first
    do while (associated(ptr))
      write(unit,'(A)') trim(ptr % id(1))//" "//trim(ptr % params)
      ptr => ptr % next
    end do
  end subroutine tPlaymol_write_xyz

  !=================================================================================================

  subroutine tPlaymol_get_types( me, atom, atom_type )
    class(tPlaymol),  intent(in)  :: me
    character(sl), intent(in)  :: atom(:)
    character(sl), intent(out) :: atom_type(:)
    integer :: i, narg
    type(Struc), pointer :: ptr
    do i = 1, size(atom)
      call me % atom_list % search( atom(i:i), ptr )
      if (.not.associated(ptr)) call error( "unknown atom", atom(i) )
      call split( ptr % params, narg, atom_type(i:i) )
    end do
  end subroutine tPlaymol_get_types

  !=================================================================================================

  subroutine tPlaymol_update_structure( me )
    class(tPlaymol), intent(inout) :: me
    type(Struc), pointer :: b1, b2, bnew
    integer :: i, j
    character(sl) :: angle(3), atom(4), new1, new2
    bnew => me % bond_list % last
    new1 = bnew%id(1)
    new2 = bnew%id(2)
    b1 => me % bond_list % first
    do while (associated(b1 % next))
      if (b1%id(1) == new1) then
        i = 2; j = 1
      else if (b1%id(2) == new1) then
        i = 1; j = 1
      else if (b1%id(1) == new2) then
        i = 2; j = 2
      else if (b1%id(2) == new2) then
        i = 1; j = 2
      else
        i = 0
      end if
      if (i /= 0) then
        angle(1) = b1%id(i)
        angle(2) = bnew%id(j)
        angle(3) = bnew%id(3-j)
        atom(1:3) = angle
        call me % angle_list % handle( atom(1:3), me % atom_list, me % angle_type_list, 1 )
        b2 => me % bond_list % first
        do while (.not.associated(b2,b1))
          if (b2%id(1) == angle(1)) then
            atom = [b2%id(2),angle]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          else if (b2%id(2) == angle(1)) then
            atom = [b2%id(1),angle]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          end if
          b2 => b2 % next
        end do
        b2 => b1 % next
        do while (associated(b2 % next))
          if (b2%id(1) == angle(3)) then
            atom = [angle,b2%id(2)]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          else if (b2%id(2) == angle(3)) then
            atom = [angle,b2%id(1)]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          end if
          if (b2%id(1) == angle(1)) then
            atom = [b2%id(2),angle]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          else if (b2%id(2) == angle(1)) then
            atom = [b2%id(1),angle]
            call me % dihedral_list % handle( atom, me % atom_list, me % dihedral_type_list, 1 )
          end if
          b2 => b2 % next
        end do
      end if
      b1 => b1 % next
    end do
  end subroutine tPlaymol_update_structure

  !=================================================================================================

  subroutine tPlaymol_search_impropers( me )
    class(tPlaymol), intent(inout) :: me
    integer :: nimp
    type(Struc), pointer :: b1, b2, b3
    character(sl) :: atom(3,2), central
    call writeln( "Searching for impropers..." )
    nimp = me % improper_list % count
    b1 => me % bond_list % first
    if (.not.associated(b1)) return
    if (.not.associated(b1 % next)) return
    do while (associated(b1 % next % next))
      atom(1,:) = b1%id
      b2 => b1 % next
      do while (associated(b2 % next))
        atom(2,:) = b2%id
        if (any(atom(2,:) == atom(1,1))) then
          if (atom(2,2) == atom(1,1)) call str_swap( atom(2,1), atom(2,2) )
        else if (any(atom(2,:) == atom(1,2))) then
          call str_swap( atom(1,1), atom(1,2) )
          if (atom(2,2) == atom(1,1)) call str_swap( atom(2,1), atom(2,2) )
        end if
        if (atom(1,1) == atom(2,1)) then
          central = atom(1,1)
          b3 => b2 % next
          do while (associated(b3))
            atom(3,:) = b3%id
            if (any(atom(3,:) == central)) then
              if (atom(3,2) == central) call str_swap( atom(3,1), atom(3,2) )
              call check_improper( atom(1,2), atom(2,2), central, atom(3,2) )
              call check_improper( atom(1,2), atom(3,2), central, atom(2,2) )
              call check_improper( atom(2,2), atom(1,2), central, atom(3,2) )
              call check_improper( atom(2,2), atom(3,2), central, atom(1,2) )
              call check_improper( atom(3,2), atom(1,2), central, atom(2,2) )
              call check_improper( atom(3,2), atom(2,2), central, atom(1,2) )
            end if
            b3 => b3 % next
          end do
        end if
        b2 => b2 % next
      end do
      b1 => b1 % next
    end do
    if (me % improper_list % count == nimp) call writeln( "--> no impropers found." )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine check_improper( a1, a2, a3, a4 )
        character(sl), intent(in) :: a1, a2, a3, a4
        character(sl) :: id(4)
        id = [a1, a2, a3, a4]
        call me % improper_list % handle( id, me % atom_list, me % improper_type_list, 3 )
      end subroutine check_improper
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_search_impropers

  !=================================================================================================

  subroutine tPlaymol_write( me, unit )
    class(tPlaymol), intent(in) :: me
    integer,      intent(in) :: unit
    integer :: N
    character(sl) :: CN
    type(Struc), pointer :: current
    call me % atom_type_list % print( unit, used_only = .true. )
    call me % mass_list % print( unit, used_only = .true. )
    call me % diameter_list % print( unit, used_only = .true. )
    call me % bond_type_list % print( unit, used_only = .true. )
    call me % angle_type_list % print( unit, used_only = .true. )
    call me % dihedral_type_list % print( unit, used_only = .true. )
    call me % improper_type_list % print( unit, used_only = .true. )
    call me % atom_list % print( unit )
    call me % charge_list % print( unit, used_only = .true. )
    call me % bond_list % print( unit )
    call me % angle_list % print( unit, comment = .true. )
    call me % dihedral_list % print( unit, comment = .true. )
    call me % extra_bond_list % print( unit )
    call me % extra_angle_list % print( unit )
    call me % extra_dihedral_list % print( unit )
    call me % improper_list % print( unit )
    call me % link_list % print( unit )
    N = me % molecules % xyz % count
    if (N > 0) then
      write(CN,*) N
      write(unit,'("build"/,A,/,"# atom x y z")') trim(adjustl(CN))
      current => me % molecules % xyz % first
      do while (associated(current))
        write(unit,'(A,X,A)') trim(current % id(1)), trim(current % params)
        current => current % next
      end do
    end if
  end subroutine tPlaymol_write

  !=================================================================================================

  subroutine tPlaymol_summarize( me, output )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: output
    integer :: i, unit, imol, molcount(me%molecules%N), natoms(me%molecules%N)
    integer, allocatable :: output_unit(:)
    real(rb) :: mass(me%molecules%N), charge(me%molecules%N)
    type(Struc), pointer :: ptr
    character(sl) :: atoms
    natoms = me % molecules % number_of_atoms()
    molcount = me % molecules % count()
    charge = me % molecules % per_molecule( me % charge_list )
    mass = me % molecules % per_molecule( me % atom_masses )
    if ((output == stdout).and.(logunit /= 0)) then
      allocate( output_unit(2) )
      output_unit(1) = stdout
      output_unit(2) = logunit
    else
      allocate( output_unit(1) )
      output_unit(1) = output
    end if
    do i = 1, size(output_unit)
      unit = output_unit(i)
      write(unit,'(A,/,"SUMMARY",/,A,/,"Specified:")') repeat("-",80), repeat("-",80)
      call flush_data( unit, "atom type", me % atom_type_list % count )
      call flush_data( unit, "bond type", me % bond_type_list % count )
      call flush_data( unit, "angle type", me % angle_type_list % count )
      call flush_data( unit, "dihedral type", me % dihedral_type_list % count )
      call flush_data( unit, "improper type", me % improper_type_list % count )
      call flush_data( unit, "atom", me % atom_list % count )
      call flush_data( unit, "bond", me % bond_list % count )
      call flush_data( unit, "virtual link", me % link_list % count )
      call flush_data( unit, "extra bond", me % extra_bond_list % count )
      call flush_data( unit, "extra angle", me % extra_angle_list % count )
      call flush_data( unit, "extra dihedral", me % extra_dihedral_list % count )
      call flush_data( unit, "improper", me % improper_list % count )
      write(unit,'(/,A)') "Detected:"
      call flush_data( unit, "angle", me % angle_list % count )
      call flush_data( unit, "dihedral", me % dihedral_list % count )
      call flush_data( unit, "molecule", me%molecules%N )
      write(unit,'(/,A)') "Effectively used:"
      call flush_data( unit, "atom type", me % atom_type_list % count_used() )
      call flush_data( unit, "bond type", me % bond_type_list % count_used() )
      call flush_data( unit, "angle type", me % angle_type_list % count_used() )
      call flush_data( unit, "dihedral type", me % dihedral_type_list % count_used() )
      call flush_data( unit, "improper type", me % improper_type_list % count_used() )
      do imol = 1, me%molecules%N
        write(unit,'(/,"Molecule[",A,"]:")') trim(int2str(imol))
        write(unit,'("- Amount: ",A)') trim(int2str(molcount(imol)))
        write(unit,'("- Mass: ",A)') trim(real2str(mass(imol)))
        write(unit,'("- Charge: ",A)') trim(real2str(charge(imol)))
        write(unit,'("- Number of atoms: ",A)') trim(int2str(natoms(imol)))
        atoms = "- Atoms:"
        ptr => me % molecules % list % first
        do while (associated(ptr))
          if (str2int(ptr % params) == imol) then
            if (len_trim(atoms) + len_trim(ptr % id(1)) > 79) then
              write(unit,'(A)') trim(atoms)
              atoms = ""
            end if
            atoms = trim(atoms)//" "//ptr%id(1)
          end if
          ptr => ptr % next
        end do
        write(unit,'(A)') trim(atoms)
      end do
      write(unit,'(/,"System:")')
      write(unit,'("- Molecules: ",A)') trim(int2str(sum(molcount)))
      write(unit,'("- Total mass: ",A)') trim(real2str(sum(molcount*mass)))
      if (me % box % exists()) then
        call me % box % compute( sum(molcount*mass) )
        write(unit,'("- Box lengths: ",A)') trim(join(real2str(me % box % length)))
        if (me % box % def_type == 4) then
          write(unit,'("- Box angles: ",A)') trim(join(real2str(me % box % angle)))
        end if
        write(unit,'("- Box density: ",A)') trim(real2str(me % box % density))
      end if
      write(unit,'("- Residual charge: ",A)') trim(real2str(sum(molcount*charge)))
      write(unit,'(A)') repeat("-",80)
    end do
    contains
      subroutine flush_data( unit, title, amount )
        integer,      intent(in) :: unit
        character(*), intent(in) :: title
        integer,      intent(in) :: amount
        if (amount == 1) then
          write(unit,'(I5,X,A,".")') amount, title
        else if (amount > 1) then
          write(unit,'(I5,X,A,"s.")') amount, title
        end if
      end subroutine flush_data
  end subroutine tPlaymol_summarize

  !=================================================================================================

  subroutine tPlaymol_write_lammps( me, unit )
    class(tPlaymol),  intent(inout) :: me
    integer,          intent(in)    :: unit
    integer :: i, j
    type StrucHolder
      integer :: multiplicity = 0
      integer :: mol
      character(sl) :: charge
      character(sl), allocatable :: atoms(:), atom_types(:)
      integer, allocatable :: atom_in_mol(:), itype(:)
    end type StrucHolder
    type(StrucHolder) :: atom(me % atom_list % count), bond(me % bond_list % count),    &
                         ang(me % angle_list % count), dih(me % dihedral_list % count), &
                         imp(me % improper_list % count)
    type TypeHolder
      integer :: index
      character(sl), allocatable :: types(:)
      character(sl) :: params
      character(sl) :: mass
    end type TypeHolder
    type(TypeHolder), allocatable :: atom_types(:), bond_types(:), &
                                     ang_types(:),  dih_types(:),  &
                                     imp_types(:)
    type tCounter
      integer :: atoms, bonds, angs, dihs, imps, mols
    end type tCounter
    type(tCounter) :: n(me%molecules%N), total(me%molecules%N)
    type(Struc), pointer :: ptr
    integer, allocatable :: mol_index(:)
    ! Molecules:
    n%mols = me % molecules % count()
    ! Atoms:
    call analyze( atom, n%atoms, total%atoms, atom_types, &
                  me%atom_list, me%atom_type_list, n%mols )
    ! Bonds:
    call attach( me%bond_list, me%extra_bond_list )
    call analyze( bond, n%bonds, total%bonds, bond_types,  &
                  me%bond_list, me%bond_type_list, n%mols, &
                  atom )
    call detach( me%bond_list, me%extra_bond_list )
    ! Angles:
    call attach( me%angle_list, me%extra_angle_list )
    call analyze( ang, n%angs, total%angs, ang_types,        &
                  me%angle_list, me%angle_type_list, n%mols, &
                  atom )
    call detach( me%angle_list, me%extra_angle_list )
    ! Dihedrals:
    call attach( me%dihedral_list, me%extra_dihedral_list )
    call analyze( dih, n%dihs, total%dihs, dih_types,              &
                  me%dihedral_list, me%dihedral_type_list, n%mols, &
                  atom )
    call detach( me%dihedral_list, me%extra_dihedral_list )
    ! Impropers:
    call analyze( imp, n%imps, total%imps, imp_types,              &
                  me%improper_list, me%improper_type_list, n%mols, &
                  atom )
    ! Molecule indices:
    allocate( mol_index(sum(n%mols)) )
    ptr => me % molecules % xyz % first
    do i = 1, size(mol_index)
      mol_index(i) = str2int(me % molecules % list % parameters( ptr % id(1:1) ))
      do j = 1, n(mol_index(i))%atoms
        ptr => ptr % next
      end do
    end do
    ! Write LAMMPS data file:
    write(unit,'("LAMMPS data file",/,"# Generated by Playmol on ",A,/)') trim(now())
    call write_count( size(atom_types), "atom types"     )
    call write_count( size(bond_types), "bond types"     )
    call write_count( size(ang_types),  "angle types"    )
    call write_count( size(dih_types),  "dihedral types" )
    call write_count( size(imp_types),  "improper types" )
    write(unit,'()')
    call write_count( sum(n%mols * total%atoms), "atoms"     )
    call write_count( sum(n%mols * total%bonds), "bonds"     )
    call write_count( sum(n%mols * total%angs),  "angles"    )
    call write_count( sum(n%mols * total%dihs),  "dihedrals" )
    call write_count( sum(n%mols * total%imps),  "impropers" )
    if (me % box % exists()) call write_box_limits
    call write_masses
    call write_type( "Pair Coeffs",     atom_types, me % atom_type_list     )
    call write_type( "Bond Coeffs",     bond_types, me % bond_type_list     )
    call write_type( "Angle Coeffs",    ang_types,  me % angle_type_list    )
    call write_type( "Dihedral Coeffs", dih_types,  me % dihedral_type_list )
    call write_type( "Improper Coeffs", imp_types,  me % improper_type_list )
    call write_atoms( n%atoms )
    call write_structure( "Bonds",     bond, n%bonds, mol_index, n%atoms )
    call write_structure( "Angles",    ang,  n%angs,  mol_index, n%atoms )
    call write_structure( "Dihedrals", dih,  n%dihs,  mol_index, n%atoms )
    call write_structure( "Impropers", imp,  n%imps,  mol_index, n%atoms )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine analyze( structure, permol, total, type_map, list, typelist, nmols, atom )
        type(StrucList),   intent(in)               :: list, typelist
        integer,           intent(in)               :: nmols(me%molecules%N)
        type(StrucHolder), intent(out)              :: structure(list%count)
        type(StrucHolder), intent(in),  optional    :: atom(:)
        integer,           intent(out)              :: permol(me%molecules%N)
        integer,           intent(out)              :: total(me%molecules%N)
        type(TypeHolder),  intent(out), allocatable :: type_map(:)
        integer :: i, j, k, m, n, itype, imol, imax
        type(Struc), pointer :: ptr
        integer, allocatable :: aux(:)
        if (list%count > 0) then
          ! Identify the atoms and types in each structure and
          ! count the number of structures in each molecule:
          permol = 0
          ptr => list%first
          m = list%number
          do i = 1, list%count
            associate (s => structure(i))
              allocate(s%atoms(m), s%atom_types(m), s%atom_in_mol(m), s%itype(0) )
              S%atoms = ptr%id
              call me % get_types( ptr%id, s%atom_types )
              imol = str2int(me % molecules % list % parameters( ptr%id(1:1) ))
              s%mol = imol
              permol(imol) = permol(imol) + 1
              if (present(atom)) then
                do j = 1, m
                  k = 1
                  do while (atom(k)%atoms(1) /= s%atoms(j))
                    k = k + 1
                  end do
                  s%atom_in_mol(j) = atom(k)%atom_in_mol(1)
                end do
              end if
              ptr => ptr%next
            end associate
          end do
          ! Sort structures according to the molecules they belong to:
          structure = structure(sorted( structure%mol ))
          ! If structure = atom, determine atom_in_mol and charge:
          if (.not.present(atom)) then
            imol = 0
            do i = 1, list%count
              associate(s => structure(i))
                if (s%mol > imol) then
                  j = 0
                  imol = s%mol
                end if
                j = j + 1
                s%atom_in_mol = j
                s%charge = me % charge_list % parameters( s%atoms, default = "0.0" )
              end associate
            end do
          end if
          ! Determine how many times the same structure appears in a molecule:
          total = 0
          ptr => typelist % first
          do itype = 1, typelist % count
            if (ptr % used) then
              do i = 1, list%count
                associate (s => structure(i))
                  if (ptr % match_id( s%atom_types, typelist%two_way )) then
                    n = s%multiplicity + 1
                    allocate( aux(n) )
                    aux = [s%itype, itype]
                    call move_alloc( aux, s%itype )
                    s%multiplicity = n
                    total(s%mol) = total(s%mol) + 1
                  end if
                end associate
              end do
            end if
            ptr => ptr % next
          end do
          ! Create a map for the indices of actually used types:
          imax = maxval([(maxval(structure(i)%itype),i=1,list%count)])
          allocate( aux(imax) )
          aux = 0
          forall (i=1:list%count, nmols(structure(i)%mol) > 0) aux(structure(i)%itype) = 1
          allocate( type_map(count(aux == 1)) )
          type_map%index = pack([(i,i=1,imax)],aux == 1)
          ! Retrieve the parameters of each used type:
          do i = 1, size(type_map)
            allocate( type_map(i) % types(m) )
            ptr => typelist % point_to( type_map(i)%index )
            type_map(i) % types = ptr % id
            type_map(i) % params = ptr % params
            if (.not.present(atom)) type_map(i) % mass = me % mass_list % parameters( ptr % id )
          end do
          ! Apply an inverse map to the indices:
          aux(type_map%index) = [(i,i=1,size(type_map))]
          forall (i=1:list%count) structure(i)%itype = aux(structure(i)%itype)
        else
          permol = 0
          total = 0
          allocate( type_map(0) )
        end if
      end subroutine analyze
      !---------------------------------------------------------------------------------------------
      subroutine attach( list, extra_list )
        type(StrucList), intent(inout) :: list, extra_list
        if (associated(list % last)) then
          list % last % next => extra_list % first
        else
          list % first => extra_list % first
        end if
        list%count = list%count + extra_list%count
      end subroutine attach
      !---------------------------------------------------------------------------------------------
      subroutine detach( list, extra_list )
        type(StrucList), intent(inout) :: list, extra_list
        if (associated(list % last)) then
          list % last % next => null()
        else
          list % first => null()
        end if
        list%count = list%count - extra_list%count
      end subroutine detach
      !---------------------------------------------------------------------------------------------
      subroutine write_count( n, name )
        integer,      intent(in) :: n
        character(*), intent(in) :: name
        if (n > 0) then
          write(unit,'(A,X,A)') trim(int2str(n)), name
          call writeln( int2str(n), name )
        end if
      end subroutine write_count
      !---------------------------------------------------------------------------------------------
      subroutine write_box_limits
        real(rb), parameter :: Deg2Rad = 0.01745329252_rb, tol = 1.0e-10_rb
        real(rb) :: lim(2) = [-0.5_rb, 0.5_rb]
        real(rb) :: L(3), cosine(3), lx, ly, lz, xy, xz, yz
        call me % box % compute( sum(n%mols*me % molecules % per_molecule( me % atom_masses )) )
        L = me%box%length
        write(unit,'()')
        if (me % box % def_type /= 4) then ! Orthogonal box
          write(unit,'(A," xlo xhi")') trim(join(real2str(lim*L(1))))
          write(unit,'(A," ylo yhi")') trim(join(real2str(lim*L(2))))
          write(unit,'(A," zlo zhi")') trim(join(real2str(lim*L(3))))
        else ! Triclinic box
          cosine = cos(Deg2Rad*me%box%angle)
          lx = L(1)
          xy = L(2)*cosine(3); if (abs(xy) < tol) xy = 0.0_rb
          xz = L(3)*cosine(2); if (abs(xz) < tol) xz = 0.0_rb
          ly = sqrt(L(2)**2 - xy**2)
          yz = (L(2)*L(3)*cosine(1) - xy*xz)/ly; if (abs(yz) < tol) yz = 0.0_rb
          lz = sqrt(L(3)**2 - xz**2 - yz**2)
          write(unit,'(A," xlo xhi")') trim(join(real2str(lim*lx)))
          write(unit,'(A," ylo yhi")') trim(join(real2str(lim*ly)))
          write(unit,'(A," zlo zhi")') trim(join(real2str(lim*lz)))
          write(unit,'(A," xy xz yz")') trim(join(real2str([xy,xz,yz])))
        end if
      end subroutine write_box_limits
      !---------------------------------------------------------------------------------------------
      subroutine write_type( title, types, list )
        character(*),      intent(in) :: title
        type(TypeHolder),  intent(in) :: types(:)
        type(StrucList),   intent(in) :: list
        integer :: i
        if (size(types) > 0) write(unit,'(/,A,/)') title
        do i = 1, size(types)
          write(unit,'(A)') trim(join([int2str(i), types(i)%params, "#", types(i)%types]))
        end do
      end subroutine write_type
      !---------------------------------------------------------------------------------------------
      subroutine write_masses
        integer :: i
        write(unit,'(/,"Masses",/)')
        do i = 1, size(atom_types)
          write(unit,'(A)') trim(join([int2str(i), atom_types(i)%mass, "#", atom_types(i)%types]))
        end do
      end subroutine write_masses
      !---------------------------------------------------------------------------------------------
      subroutine write_atoms( natoms )
        integer, intent(in) :: natoms(:)
        type(Struc), pointer :: patom
        integer :: i, j, k, m, katom, imol, kmol, prev
        character(sl) :: cstruc
        integer, allocatable :: iatom(:)
        character(sl), allocatable :: catom(:), xyz(:)
        if (any(n%mols*n%atoms > 0)) then
          write(unit,'(/,"Atoms",/)')
          m = maxval(n%atoms,n%mols > 0)
          allocate( iatom(m), catom(m), xyz(m) )
          patom => me % molecules % xyz % first
          katom = 0
          do kmol = 1, size(mol_index)
            imol = mol_index(kmol)
            prev = sum(natoms(1:imol-1))
            m = natoms(imol)
            forall (j=1:m) catom(j) = atom(prev+j)%atoms(1)
            do j = 1, m
              iatom(j:j) = pack([(k,k=1,m)],catom(1:m) == patom%id(1))
              xyz(j) = patom%params
              patom => patom % next
            end do
            xyz(1:m) = xyz(sorted(iatom(1:m)))
            do j = 1, m
              katom = katom + 1
              i = prev + j
              cstruc = join(int2str([katom, kmol, atom(i)%itype]))
              write(unit,'(A)') trim(join([cstruc, atom(i)%charge, xyz(j), "#", catom(j)]))
            end do
          end do
        end if
      end subroutine write_atoms
      !---------------------------------------------------------------------------------------------
      subroutine write_structure( title, structure, permol, mol_index, natoms )
        character(*),      intent(in) :: title
        type(StrucHolder), intent(in) :: structure(:)
        integer,           intent(in) :: permol(me%molecules%N), mol_index(:), natoms(:)
        integer :: i, j, k, kmol, imol, sprev, aprev
        character(sl) :: cstruc
        if (any(n%mols*permol > 0)) then
          write(unit,'(/,A,/)') title
          k = 0
          aprev = 0
          do kmol = 1, size(mol_index)
            imol = mol_index(kmol)
            sprev = sum(permol(1:imol-1))
            do i = 1, permol(imol)
              associate (s => structure(sprev+i))
                do j = 1, s%multiplicity
                  k = k + 1
                  cstruc = join(int2str([k, s%itype(j), aprev + s%atom_in_mol]))
                  write(unit,'(A)') trim(join([cstruc, "#", s%atoms]))
                end do
              end associate
            end do
            aprev = aprev + natoms(imol)
          end do
        end if
      end subroutine write_structure
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_write_lammps

  !=================================================================================================

  subroutine tPlaymol_write_lammps_old( me, unit )
    class(tPlaymol),  intent(inout) :: me
    integer,          intent(in)    :: unit
    integer :: nb, na, nd, ni, natoms(me%molecules%N)
    logical :: bond_exists, angle_exists, dihedral_exists
    natoms = me % molecules % number_of_atoms()
    call link_extra( bond_exists, me % bond_list, me % extra_bond_list )
    call link_extra( angle_exists, me % angle_list, me % extra_angle_list )
    call link_extra( dihedral_exists, me % dihedral_list, me % extra_dihedral_list )
    write(unit,'("LAMMPS data file",/,"# Generated by Playmol on ",A,/)') trim(now())
    call write_count( me % atom_type_list % count_used(), "atom types" )
    call write_count( me % bond_type_list % count_used(), "bond types" )
    call write_count( me % angle_type_list % count_used(), "angle types" )
    call write_count( me % dihedral_type_list % count_used(), "dihedral types" )
    call write_count( me % improper_type_list % count_used(), "improper types" )
    write(unit,'()')
    call write_count( me % molecules % xyz % count, "atoms" )
    call handle_struc( "bonds", me % bond_list, me % bond_type_list, nb )
    call write_count( nb, "bonds" )
    call handle_struc( "angles", me % angle_list, me % angle_type_list, na )
    call write_count( na, "angles" )
    call handle_struc( "dihedrals", me % dihedral_list, me % dihedral_type_list, nd )
    call write_count( nd, "dihedrals" )
    call handle_struc( "impropers", me % improper_list, me % improper_type_list, ni )
    call write_count( ni, "impropers" )
    if (me % box % exists()) call write_box_limits
    call write_masses
    call write_type( "Pair Coeffs", me % atom_type_list )
    call write_type( "Bond Coeffs", me % bond_type_list )
    call write_type( "Angle Coeffs", me % angle_type_list )
    call write_type( "Dihedral Coeffs", me % dihedral_type_list )
    call write_type( "Improper Coeffs", me % improper_type_list )
    call write_atoms
    if (nb > 0) call handle_struc( "Bonds", me % bond_list, me % bond_type_list )
    if (na > 0) call handle_struc( "Angles", me % angle_list, me % angle_type_list )
    if (nd > 0) call handle_struc( "Dihedrals", me % dihedral_list, me % dihedral_type_list )
    if (ni > 0) call handle_struc( "Impropers", me % improper_list, me % improper_type_list )
    call unlink_extra( bond_exists, me % bond_list )
    call unlink_extra( angle_exists, me % angle_list )
    call unlink_extra( dihedral_exists, me % dihedral_list )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine link_extra( exists, list, extra_list )
        logical,         intent(inout) :: exists
        type(StrucList), intent(inout) :: list, extra_list
        exists = associated(list % first)
        if (exists) then
          list % last % next => extra_list % first
        else
          list % first => extra_list % first
        end if
      end subroutine link_extra
      !---------------------------------------------------------------------------------------------
      subroutine unlink_extra( exists, list )
        logical,         intent(in)    :: exists
        type(StrucList), intent(inout) :: list
        if (exists) then
          list % last % next => null()
        else
          list % first => null()
        end if
      end subroutine unlink_extra
      !---------------------------------------------------------------------------------------------
      subroutine write_count( n, name )
        integer,      intent(in) :: n
        character(*), intent(in) :: name
        if (n > 0) then
          write(unit,'(A,X,A)') trim(int2str(n)), name
          call writeln( int2str(n), name )
        end if
      end subroutine write_count
      !---------------------------------------------------------------------------------------------
      subroutine write_box_limits
        integer :: i, molcount(me%molecules%N)
        real(rb) :: mass(me%molecules%N)
        real(rb), parameter :: Deg2Rad = 0.01745329251994329577_rb
        character :: dir(3) = ["x","y","z"]
        character(sl) :: limits
        real(rb) :: a, b, c, alpha, beta, gamma, lx, ly, lz, xy, xz, yz, L(3)
        molcount = me % molecules % count()
        mass = me % molecules % per_molecule( me % atom_masses )
        call me % box % compute( sum(molcount*mass) )
        write(unit,'()')
        if (me % box % def_type /= 4) then
          do i = 1, 3
            limits = join(real2str( me%box%length(i)*[-0.5_rb,+0.5_rb] ))
            write(unit,'(A," ",A,"lo ",A,"hi")') trim(limits), dir(i), dir(i)
          end do
        else ! Triclinic box
          a = me%box%length(1)
          b = me%box%length(2)
          c = me%box%length(3)
          alpha = Deg2Rad*me%box%angle(1)
          beta = Deg2Rad*me%box%angle(2)
          gamma = Deg2Rad*me%box%angle(3)
          lx = a
          xy = b*cos(gamma)
          xz = c*cos(beta)
          ly = sqrt(b**2 - xy**2)
          yz = (b*c*cos(alpha) - xy*xz)/ly
          lz = sqrt(c**2 - xz**2 - yz**2)
          L = [lx,ly,lz]
          do i = 1, 3
            limits = join(real2str( L(i)*[-0.5_rb,+0.5_rb] ))
            write(unit,'(A," ",A,"lo ",A,"hi")') trim(limits), dir(i), dir(i)
          end do
          L = [xy,xz,yz]
          where (abs(L) < 1.e-10_rb) L = 0.0_rb
          write(unit,'(A," xy xz yz")') trim(join([real2str(L)]))
        end if
      end subroutine write_box_limits
      !---------------------------------------------------------------------------------------------
      subroutine write_type( name, list )
        character(*),     intent(in) :: name
        class(StrucList), intent(in) :: list
        type(Struc), pointer :: current
        integer :: i
        if (list % count_used() > 0) then
        current => list % first
          write(unit,'(/,A,/)') name
          i = 0
          do while (associated(current))
            if (current % used) then
              i = i + 1
              write(unit,'(A,X,A," # ",A)') trim(int2str(i)), trim(current % params), &
                                            trim(join(current % id))
            end if
            current => current % next
          end do
        end if
      end subroutine write_type
      !---------------------------------------------------------------------------------------------
      subroutine write_masses
        type(Struc), pointer :: current
        character(sl) :: mass
        integer :: i
        if (associated(me % mass_list % first)) then
          current => me % atom_type_list % first
          i = 0
          write(unit,'(/,"Masses",/)')
          do while (associated(current))
            if (current % used) then
              i = i + 1
              mass = me % mass_list % parameters( current % id )
              if (mass == "") call error( "undefined mass for atom type", current%id(1) )
              write(unit,'(A,X,A," # ",A)') trim(int2str(i)), trim(mass), trim(join(current % id))
            end if
            current => current % next
          end do
        end if
      end subroutine write_masses
      !---------------------------------------------------------------------------------------------
      subroutine write_atoms
        type(Struc), pointer :: current, atom_type
        integer :: iatom, itype, i, imol, jmol, narg
        character(sl) :: arg(1), charge
        logical :: found
        write(unit,'(/,"Atoms",/)')
        current => me % molecules % xyz % first
        iatom = 0
        jmol = 0
        do while (associated(current))
          jmol = jmol + 1
          imol = str2int(me % molecules % list % parameters( current % id ) )
          do i = 1, natoms(imol)
            iatom = iatom + 1
            call split( me % atom_list % parameters( current % id ), narg, arg )
            atom_type => me % atom_type_list % first
            itype = 0
            found = .false.
            do while (associated(atom_type).and.(.not.found))
              if (atom_type % used) then
                itype = itype + 1
                found = atom_type % match_id( arg, two_way = .false. )
              end if
              atom_type => atom_type % next
            end do
            charge = me % charge_list % parameters( current % id )
            if (charge == "") charge = "0.0"
            write(unit,'(3(A,X),"# ",A)') trim(join(int2str([iatom,jmol,itype]))), trim(charge), &
                                          trim(current%params), trim(current % id(1))
            current => current % next
          end do
        end do
      end subroutine write_atoms
      !---------------------------------------------------------------------------------------------
      subroutine handle_struc( name, list, type_list, count )
        character(*),    intent(in)            :: name
        type(StrucList), intent(in)            :: list, type_list
        integer,         intent(out), optional :: count
        integer :: i, imol, jmol, istruc, itype, last_atom
        type(Struc), pointer :: coord, struct, struc_type
        character(sl) :: line, types(list%number), atoms(list%number)
        character(sl), allocatable :: atom_id(:)
        istruc = 0
        coord => me % molecules % xyz % first
        if (.not.present(count)) then
          call writeln( "Writing down structures:", name )
          write(unit,'(/,A,/)') name
        end if
        allocate(atom_id(maxval(natoms)))
        last_atom = 0
        do while (associated(coord))
          imol = str2int(me % molecules % list % parameters( coord%id ) )
          do i = 1, natoms(imol)
            atom_id(i) = coord % id(1)
            coord => coord % next
          end do
          struct => list % first
          do while (associated(struct))
            jmol = str2int(me % molecules % list % parameters(struct%id(1:1)))
            if (jmol == imol) then
              call me % get_types( struct%id, types )
              itype = 1
              struc_type => type_list % first
              do while (associated(struc_type))
                if (struc_type % used) then
                  if (struc_type % match_id( types, type_list % two_way)) then
                    istruc = istruc + 1
                    if (.not.present(count)) then
                      if (type_list % two_way) then
                        if (struc_type % match_id( types, .false.)) then
                          atoms = struct%id
                        else
                          atoms = struct%id(list%number:1:-1)
                        end if
                      else
                        atoms = struct%id
                      end if
                      line = join(int2str([istruc, itype, last_atom + str_find(atoms, atom_id)]))
                      write(unit,'(A," # ",A)') trim(line), trim(join(atoms))
                    end if
                  end if
                  itype = itype + 1
                end if
                struc_type => struc_type % next
              end do
            end if
            struct => struct % next
          end do
          last_atom = last_atom + natoms(imol)
        end do
        if (present(count)) count = istruc
      end subroutine handle_struc
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_write_lammps_old

  !=================================================================================================

  subroutine tPlaymol_write_lammpstrj( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    integer :: molcount(me%molecules%N), iatom, itype, i, imol, jmol, narg, natoms(me%molecules%N)
    real(rb) :: mass(me%molecules%N)
    character(sl) :: limits, arg(1)
    type(Struc), pointer :: current, atom_type
    logical :: found
    if (.not.me % box % exists()) call error( "simulation box has not been defined" )
    write(unit,'("ITEM: TIMESTEP",/,"0")')
    write(unit,'("ITEM: NUMBER OF ATOMS")')
    write(unit,'(A)') trim(int2str(me % molecules % xyz % count))
    write(unit,'("ITEM: BOX BOUNDS pp pp pp")')
    molcount = me % molecules % count()
    mass = me % molecules % per_molecule( me % atom_masses )
    call me % box % compute( sum(molcount*mass) )
    do i = 1, 3
      limits = join(real2str( me%box%length(i)*[-0.5_rb,+0.5_rb] ))
      write(unit,'(A)') trim(limits)
    end do
!    write(unit,'("ITEM: ATOMS id mol type x y z ix iy iz ")')
    write(unit,'("ITEM: ATOMS id mol type x y z")')
    natoms = me % molecules % number_of_atoms()
    current => me % molecules % xyz % first
    iatom = 0
    jmol = 0
    do while (associated(current))
      jmol = jmol + 1
      imol = str2int(me % molecules % list % parameters( current % id ) )
      do i = 1, natoms(imol)
        iatom = iatom + 1
        call split( me % atom_list % parameters( current % id ), narg, arg )
        atom_type => me % atom_type_list % first
        itype = 0
        found = .false.
        do while (associated(atom_type).and.(.not.found))
          if (atom_type % used) then
            itype = itype + 1
            found = atom_type % match_id( arg, two_way = .false. )
          end if
          if (.not.found) atom_type => atom_type % next
        end do
!        write(unit,'(3(A,X),"0 0 0")') trim(join(int2str([iatom,jmol,itype]))), trim(current%params)
        write(unit,'(3(A,X))') trim(join([int2str([iatom,jmol]),atom_type%id(1)])), trim(current%params)
        current => current % next
      end do
    end do
  end subroutine tPlaymol_write_lammpstrj

  !=================================================================================================

  subroutine tPlaymol_check_coordinates( me, imol )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: imol
    type(Struc), pointer :: coord, atom, first, current
    logical :: found
    integer :: natoms
    coord => me % molecules % xyz % first
    found = .false.
    do while (associated(coord).and.(.not.found))
      found = str2int(me % molecules % list % parameters( coord%id ) ) == imol
      coord => coord % next
    end do
    if (.not.found) then
      call warning( "No predefined coordinates for molecule", int2str(imol) )
      atom => me % molecules % list % first
      natoms = 0
      do while (associated(atom))
        if (str2int(atom%params) == imol) then
          if (natoms == 0) then
            allocate( first )
            current => first
          else
            allocate( current % next )
            current => current % next
          end if
          natoms = natoms + 1
          current % params = atom % id(1)
        end if
        atom => atom % next
      end do
      select case (natoms)
        case (1); call build_monoatomic_molecule
        case (2); call build_diatomic_molecule
        case (3); call build_triatomic_molecule
        case default; call error( "cannot guess coordinates" )
      end select
    end if
    contains
      !-------------------------------------------------------------------------
      subroutine build_monoatomic_molecule
        character(sl) :: arg(4)
        call warning( "  Ok: molecule", int2str(imol), "is monoatomic" )
        arg = [first%params, "0.0", "0.0", "0.0"]
        call me % molecules % xyz % add( 4, arg )
      end subroutine build_monoatomic_molecule
      !-------------------------------------------------------------------------
      subroutine build_diatomic_molecule
        integer :: narg
        character(sl) :: atoms(2), types(2), arg(10), R0
        call warning( "  molecule", int2str(imol), "is diatomic" )
        atoms = [first % params, first % next % params]
        call me % get_types( atoms, types )
        call split( me % bond_type_list % parameters( types ), narg, arg )
        if (narg /= 2) call error( "cannot guess coordinates" )
        R0 = arg(2)
        if (.not.is_real(R0)) call error( "cannot guess coordinates" )
        call warning( "  considering harmonic potential (if wrong, please provide coordinates)" )
        call writeln( "R0 = ", R0 )
        arg(1:4) = [atoms(1), "0.0", "0.0", "0.0"]
        call me % molecules % xyz % add( 4, arg )
        arg(1:4) = [atoms(2), R0, "0.0", "0.0"]
        call me % molecules % xyz % add( 4, arg )
      end subroutine build_diatomic_molecule
      !-------------------------------------------------------------------------
      subroutine build_triatomic_molecule
        integer :: narg
        character(sl) :: atoms(3), types(3), arg(10), R0A, R0B, Theta0
        real(rb) :: X3, Y3
        call warning( "  molecule", int2str(imol), "is triatomic" )
        atoms = [first % params, first % next % params, first % next % next % params ]
        call me % bond_list % search( atoms(1:2), current )
        if (associated(current)) then
          call me % bond_list % search( atoms(2:3), current )
          if (.not.associated(current)) atoms = atoms([2,1,3])
        else
          atoms = atoms([1,3,2])
        end if
        call me % get_types( atoms, types )
        call split( me % bond_type_list % parameters( types(1:2) ), narg, arg )
        if (narg /= 2) call error( "cannot guess coordinates" )
        R0A = arg(2)
        call split( me % bond_type_list % parameters( types(2:3) ), narg, arg )
        if (narg /= 2) call error( "cannot guess coordinates" )
        R0B = arg(2)
        call split( me % angle_type_list % parameters( types ), narg, arg )
        if (narg /= 2) call error( "cannot guess coordinates" )
        Theta0 = arg(2)
        if (is_real(R0A).and.is_real(R0B).and.is_real(Theta0)) then
        call warning( "  considering harmonic potentials (if wrong, please provide coordinates)" )
          call writeln( "R0(",atoms(1),"-",atoms(2),") = ", R0A )
          call writeln( "R0(",atoms(2),"-",atoms(3),") = ", R0B )
          call writeln( "Theta0(",atoms(1),"-",atoms(2),"-",atoms(3),") = ", Theta0, "degrees" )
        else
          call error( "cannot guess coordinates" )
        end if
        arg(1:4) = [atoms(1), "0.0", "0.0", "0.0"]
        call me % molecules % xyz % add( 4, arg )
        arg(1:4) = [atoms(2), R0A, "0.0", "0.0"]
        call me % molecules % xyz % add( 4, arg )
        X3 = str2real(R0A) - str2real(R0B)*cos(0.0174532925_rb*str2real(Theta0))
        Y3 = str2real(R0B)*sin(0.0174532925_rb*str2real(Theta0))
        arg(1:4) = [atoms(3), real2str(X3), real2str(Y3), "0.0"]
        call me % molecules % xyz % add( 4, arg )
      end subroutine build_triatomic_molecule
      !-------------------------------------------------------------------------
  end subroutine tPlaymol_check_coordinates

  !=================================================================================================

end module mPlaymol
