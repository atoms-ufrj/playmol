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
! 1) Accept double negation in mathematical expressions so that -$i will work when $i is negative.
! 2) Allow only one asterisk to be used as wildcard.
! 3) unlink command: warn if molecule splitting make extra bonds/angles/dihedrals have atoms in
!    separate molecules and stop if rigid bodies have atoms in separate molecules

module mPlaymol

use mGlobal
use mStruc
use mPackmol
use mMolecule
use mBox
use mCodeFlow
use mFix

implicit none

type, private :: tVelocity
  logical  :: active = .false.
  integer  :: seed
  real(rb) :: kT
end type tVelocity

type tPlaymol
  type(tBox)      :: box
  type(tPackmol)  :: packmol
  type(tCodeFlow) :: code_flow
  type(tMolecule) :: molecules
  type(tFix)      :: atomfix, typefix
  type(tVelocity) :: velocity

  type(StrucList) :: atom_type_list      = StrucList( "atom type" )
  type(StrucList) :: bond_type_list      = StrucList( "bond type", 2 )
  type(StrucList) :: angle_type_list     = StrucList( "angle type", 3 )
  type(StrucList) :: dihedral_type_list  = StrucList( "dihedral type", 4 )
  type(StrucList) :: improper_type_list  = StrucList( "improper type", 4 )
  type(StrucList) :: mass_list           = StrucList( "mass" )
  type(StrucList) :: atom_masses         = StrucList( "atom mass" )
  type(StrucList) :: diameter_list       = StrucList( "diameter" )
  type(StrucList) :: atom_list           = StrucList( "atom" )
  type(StrucList) :: charge_list         = StrucList( "charge" )
  type(StrucList) :: bond_list           = StrucList( "bond", 2 )
  type(StrucList) :: link_list           = StrucList( "virtual link", 2 )
  type(StrucList) :: angle_list          = StrucList( "angle", 3 )
  type(StrucList) :: dihedral_list       = StrucList( "dihedral", 4, .false., .true. )
  type(StrucList) :: improper_list       = StrucList( "improper", 4, .false. )
  type(StrucList) :: body_list           = StrucList( "rigid body" )
  type(StrucList) :: atom_bodies         = StrucList( "atom body" )
  type(StrucList) :: mixing_rule_list    = StrucList( "mixing rule", 2 )

  contains
    procedure :: read => tPlaymol_Read
    procedure :: write => tPlaymol_write
    procedure :: write_internals => tPlaymol_write_internals
    procedure :: analyze_struct => tPlaymol_analyze_struct
    procedure :: bodies_in_molecules => tPlaymol_bodies_in_molecules
    procedure :: write_lammps => tPlaymol_write_lammps
    procedure :: write_lammpstrj => tPlaymol_write_lammpstrj
    procedure :: write_emdee => tPlaymol_write_emdee
    procedure :: read_geometry => tPlaymol_read_geometry
    procedure :: write_xyz => tPlaymol_write_xyz
    procedure :: summarize => tPlaymol_write_summary
    procedure :: update_structure => tPlaymol_update_structure
    procedure :: search_impropers => tPlaymol_search_impropers
    procedure :: get_types => tPlaymol_get_types
    procedure :: check_coordinates => tPlaymol_check_coordinates
end type tPlaymol

type StrucHolder
  integer :: multiplicity = 0
  integer :: mol, body
  character(sl) :: charge
  character(sl), allocatable :: atoms(:), atom_types(:)
  integer, allocatable :: atom_in_mol(:), itype(:), ibody(:), wildcards(:)
end type StrucHolder

type TypeHolder
  integer :: index
  character(sl) :: types
  character(sl) :: model
  character(sl) :: params
  character(sl) :: mass
end type TypeHolder

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
        case ("velocity"); call velocity_command
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
        case ("body"); call body_command
        case ("mixing_rule"); call mixing_rule_command
        case ("link"); call link_command
        case ("unlink"); call unlink_command
        case ("extra"); call extra_command
        case ("improper"); call improper_command
        case ("xyz","build"); call build_command
        case ("align"); call align_command
        case ("write"); call write_command
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
        ptr % usable = .true.

        call me % mass_list % search( arg(3:3), ptr )
        if (.not.associated(ptr)) call error( "mass of atom type",arg(3), "has not been defined" )
        ptr % usable = .true.
        arg(3) = ptr % params
        call me % atom_masses % add( 2, arg(2:3), silent = .true. )

        call me % charge_list % search( arg(2:2), ptr )
        if (associated(ptr)) ptr % usable = .true.
        if (narg == 4) then
          if (associated(ptr)) call error( "an applicable charge has already been defined" )
          arg(3) = arg(4)
          call me % charge_list % add( 2, arg(2:3) )
          ptr => me % charge_list % last
        end if

        call me % molecules % add_atom( arg(2) )
        call me % atom_bodies % add( 2, [arg(2),"0"], silent = .true. )

      end subroutine atom_command
      !---------------------------------------------------------------------------------------------
      subroutine charge_command
        if (narg /= 3) call error( "invalid charge command" )
        call me % atomfix % apply( arg(2) )
        call me % charge_list % add( narg-1, arg(2:3) )
        if (.not.is_real(arg(3))) call error( "invalid charge value" )
        if (me % atom_list % find( arg(2:2) )) me % charge_list % last % usable = .true.
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
      subroutine body_command
        integer :: nbodies, i, j
        character(sl) :: imol, jmol
        type(Struc), pointer :: ptr
        if (narg < 3) call error( "invalid body command" )
        nbodies = me % body_list % count + 1
        arg(1) = int2str(nbodies)
        call me % atomfix % apply( arg(2:narg) )
        call me % body_list % add( narg, arg )
        if (any(has_macros(arg(2:narg)))) call error( "atom names cannot contain macros" )
        do i = 2, narg-1
          do j = i+1, narg
            if (arg(i) == arg(j)) call error( "repeated atom name", arg(i) )
          end do
        end do
        do i = 2, narg
          call me % atom_bodies % search( arg(i:i), ptr )
          if (.not.associated(ptr)) call error( "there is no atom identified as", arg(i) )
          if (ptr%params == "0") then
            ptr%params = int2str(nbodies)
          else
            call error( "atom", arg(i), "has already been added to body", ptr%params )
          end if
        end do
        imol = me % molecules % list % parameters( arg(2:2) )
        do i = 3, narg
          jmol = me % molecules % list % parameters( arg(i:i) )
          if (jmol /= imol) call error( "not all atoms are in the same molecule" )
        end do
      end subroutine body_command
      !---------------------------------------------------------------------------------------------
      subroutine mixing_rule_command
        call me % typefix % apply( arg(2:3) )
        call me % mixing_rule_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine mixing_rule_command
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
        if (choice /= "bond") then
          imol = str2int(me % molecules % list % parameters( arg(2:2) ))
          do i = 2, narg
            jmol = str2int(me % molecules % list % parameters( arg(i:i) ))
            if (jmol /= imol) call error("atoms", join(arg(2:narg)), "are not in the same molecule")
          end do
        end if
      end subroutine extra_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_bond_command
        call me % atomfix % apply( arg(2:3) )
        call me % bond_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 3) call error( "invalid extra bond command")
        call me % bond_list % handle( arg(2:3), me%atom_list, me%bond_type_list, 2 )
        call me % molecules % fuse( arg(2:3) )
      end subroutine extra_bond_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_angle_command
        call me % atomfix % apply( arg(2:4) )
        call me % angle_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 4) call error( "invalid extra angle command")
        call me % angle_list % handle( arg(2:4), me%atom_list, me%angle_type_list, 2 )
      end subroutine extra_angle_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_dihedral_command
        call me % atomfix % apply( arg(2:5) )
        call me % dihedral_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 5) call error( "invalid extra dihedral command")
        call me % dihedral_list % handle( arg(2:5), me%atom_list, me%dihedral_type_list, 2 )
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
          if (narg /= 5) call error( "invalid improper command" )
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
        character(sl), allocatable :: atom(:)
        integer :: axis(3)
        character :: dir(3) = ["x","y","z"]
        real(rb) :: lb(3), ub(3)
        if (narg /= 4) call error( "invalid align command" )
        call me % atomfix % apply( arg(2) )
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
        allocate( Mass(N), Coord(3,N), atom(N) )
        Mass = 1.0
        call me % molecules % coordinates( imol, N, Coord, atom, 1 )
        lb = minval(Coord,dim=2)
        ub = maxval(Coord,dim=2)
        do i = 1, 3
          call writeln( "Bounds in direction", dir(i)//":", join(real2str([lb(i),ub(i)])) )
        end do
        call writeln( "Bounding box lengths: ", join(real2str(ub-lb)) )
        call me % molecules % align( N, Mass, Coord, axis )
        call me % molecules % coordinates( imol, N, Coord, atom, 2 )
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
        integer :: unit, nspec
        character(sl) :: formats(8) = ["playmol   ", "lammps    ", "lmp/models", "summary   ", &
                                       "xyz       ", "lammpstrj ", "emdee     ", &
                                       "internals "  ]
        if (narg < 2) call error( "invalid write command" )
        if (.not.any(formats == arg(2))) call error( "invalid format", arg(2), "in write command" )
        select case (arg(2))
          case ("internals"); nspec = 3
          case default; nspec = 0
        end select
        if (narg == nspec+2) then
          call writeln( "Writing data in", arg(2), "format..." )
          unit = stdout
        else if (narg == nspec+3) then
          call writeln( "Writing data to file", arg(narg), "in", arg(2), "format..." )
          open( newunit = unit, file = arg(narg), status = "replace" )
        else
          call error( "invalid write command" )
        end if
        select case (arg(2))
          case ("playmol"); call me % write( unit )
          case ("lammps"); call me % write_lammps( unit, models = .false. )
          case ("lmp/models"); call me % write_lammps( unit, models = .true. )
          case ("summary"); call me % summarize( unit )
          case ("xyz"); call me % write_xyz( unit )
          case ("lammpstrj"); call me % write_lammpstrj( unit )
          case ("emdee"); call me % write_emdee( unit )
          case ("internals"); call me % write_internals( unit, arg(3:2+nspec) )
        end select
        if (unit /= stdout) close(unit)
      end subroutine write_command
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
        if (any(lists == 13)) call me % improper_list % destroy
        if (any(lists == 14)) call me % link_list % destroy
        if (any(lists == 15)) call me % molecules % xyz % destroy
        if (any(lists == 16)) call me % packmol % list % destroy
        ! Lists automatically defined by Playmol:
        if (any(lists == 17)) call me % molecules % bonds % destroy
        if (any(lists == 18)) call me % molecules % list % destroy
        if (any(lists == 19)) call me % atom_masses % destroy
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
          case ("all");      call reset_lists( [(i,i=1,19)] )
          case ("atom");     call reset_lists( [(i,i=8,19)] )
          case ("charge");   call reset_lists( [9] )
          case ("bond");     call reset_lists( [(i,i=10,14),17] )
          case ("angle");    call reset_lists( [11] )
          case ("dihedral"); call reset_lists( [12] )
          case ("improper"); call reset_lists( [13] )
          case ("link");     call reset_lists( [14] )
          case ("xyz");      call reset_lists( [15] )
          case ("packmol");  call reset_lists( [16] )
          case default;      call error( "invalid reset command" )
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
              call me % atomfix % apply( arg(iarg+1) )
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
      subroutine velocity_command
        call writeln( "Setting velocities with parameters seed =", arg(2), "and kT = ", arg(3) )
        if (narg /= 3) call error( "invalid velocity command" )
        me % velocity % active = .true.
        me % velocity % seed = str2int(arg(2))
        me % velocity % kT = str2real(arg(3))
        if (me % velocity % kT < 0.0_rb) call error( "invalid kT value" )
      end subroutine velocity_command
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_Read

  !=================================================================================================

  subroutine tPlaymol_read_geometry( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: N, i, narg
    character(sl) :: arg(7)
    integer,       allocatable :: ndata(:)
    character(sl), allocatable :: data(:,:)
    call me % code_flow % next_command( unit, narg, arg )
    if (narg == 0) call error( "expected geometric data not provided" )
    N = str2int( arg(1) )
    call writeln( "Number of data to be read:", int2str(N) )
    allocate( data(N,7), ndata(N) )
    do i = 1, N
      call me % code_flow % next_command( unit, narg, arg )
      call writeln( int2str(i), ":", join(arg(1:narg)) )
      if (narg == 0) call error( "expected geometric data not completed" )
      select case (narg)
        case (3)
          call me % atomfix % apply( arg(1:2) ) ! Bond
        case (4)
          call me % atomfix % apply( arg(1) )   ! Coordinates
        case (5)
          call me % atomfix % apply( arg(1:2) ) ! Bond and angle
          call me % atomfix % apply( arg(4) )
        case (7)
          call me % atomfix % apply( arg(1:2) ) ! Bond, angle, and dihedral
          call me % atomfix % apply( arg(4) )
          call me % atomfix % apply( arg(6) )
      end select
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
    character(sl) :: atom(4), atom_type(3)
    logical :: found
    call writeln( "Searching for impropers..." )
    nimp = me % improper_list % count
    b1 => me % bond_list % first
    if (.not.associated(b1)) return
    if (.not.associated(b1 % next)) return
    do while (associated(b1 % next % next))
      b2 => b1 % next
      do while (associated(b2 % next))
        if (any(b2%id == b1%id(1))) then
          atom(1:2) = b1%id
        else if (any(b2%id == b1%id(2))) then
          atom(1:2) = b1%id(2:1:-1)
        else
          atom(1) = ""
        end if
        if (atom(1) /= "") then
          atom(3) = merge(b2%id(1), b2%id(2), b2%id(1) /= atom(1))
          b3 => b2 % next
          do while (associated(b3))
            if (any(b3%id == atom(1))) then ! Branch found
              atom(4) = merge(b3%id(1), b3%id(2), b3%id(1) /= atom(1))
              call me % get_types( atom(2:4), atom_type )
              call sort( atom(2:4), atom_type )
              call check_improper( atom )
              if (.not.found) call check_improper( atom([1,2,4,3]) )
              if (.not.found) call check_improper( atom([1,3,2,4]) )
              if (.not.found) call check_improper( atom([1,3,4,2]) )
              if (.not.found) call check_improper( atom([1,4,2,3]) )
              if (.not.found) call check_improper( atom([1,4,3,2]) )
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
      subroutine check_improper( id )
        character(sl), intent(in) :: id(4)
        call me % improper_list % handle( id, me%atom_list, me%improper_type_list, 3, found )
      end subroutine check_improper
      !---------------------------------------------------------------------------------------------
      subroutine sort( a, b )
        character(sl), intent(inout) :: a(3), b(3)
        if (.not.alphabetical( b(1), b(2) )) then
          call str_swap( a(1), a(2) )
          call str_swap( b(1), b(2) )
        end if
        if (.not.alphabetical( b(2), b(3) )) then
          call str_swap( a(2), a(3) )
          call str_swap( b(2), b(3) )
        end if
        if (.not.alphabetical( b(1), b(2) )) then
          call str_swap( a(1), a(2) )
          call str_swap( b(1), b(2) )
        end if
      end subroutine sort
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

  subroutine tPlaymol_analyze_struct( me, structure, permol, total, type_map, list, typelist, &
                                      nmols, models, atom )
    class(tPlaymol),   intent(inout)            :: me
    type(StrucList),   intent(in)               :: list, typelist
    integer,           intent(in)               :: nmols(me%molecules%N)
    type(StrucHolder), intent(out)              :: structure(list%count)
    logical,           intent(in)               :: models
    type(StrucHolder), intent(in),  optional    :: atom(:)
    integer,           intent(out)              :: permol(me%molecules%N)
    integer,           intent(out)              :: total(me%molecules%N)
    type(TypeHolder),  intent(out), allocatable :: type_map(:)
    integer :: i, j, k, m, itype, imol, imax, narg
    logical :: match
    character(sl) :: arg(10)
    type(Struc), pointer :: ptr
    integer, allocatable :: aux(:)
    logical, allocatable :: valid(:)
    if (list%count > 0) then
      ! Identify the atoms and types in each structure and
      ! count the number of structures in each molecule:
      permol = 0
      ptr => list%first
      m = list%number
      do i = 1, list%count
        associate (s => structure(i))
          allocate( s%atoms(m), s%ibody(m), s%atom_types(m), s%atom_in_mol(m) )
          allocate( s%itype(0), s%wildcards(0) )
          s%atoms = ptr%id
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
              s%ibody(j) = atom(k)%body
            end do
          end if
          ptr => ptr%next
        end associate
      end do

      ! Sort structures according to the molecules they belong to:
      structure = structure(sorted( structure%mol ))

      ! If structure = atom, determine atom_in_mol, body and charge:
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
            s%body = str2int(me % atom_bodies % parameters( s%atoms ))
            s%charge = me % charge_list % parameters( s%atoms, default = "0.0" )
          end associate
        end do
      end if

      ! Determine how many times the same structure appears in a molecule:
      total = 0
      ptr => typelist % first
      do itype = 1, typelist % count
        if (ptr % usable) then
          do i = 1, list%count
            associate (s => structure(i))
              if (list%bothways) then
                match = ptr % match_id( s%atom_types )
              else
                match = all(match_str( ptr%id, s%atom_types ))
              end if
              if (match) then
                aux = [s%itype, itype]
                call move_alloc( aux, s%itype )
                aux = [s%wildcards, count(has_macros(ptr%id))]
                call move_alloc( aux, s%wildcards)
                s%multiplicity = s%multiplicity + 1
                total(s%mol) = total(s%mol) + 1
              end if
            end associate
          end do
        end if
        ptr => ptr % next
      end do

      ! Resolve ambiguities:
      do i = 1, list%count
        associate (s => structure(i))
          if (s%multiplicity > 1) then
            if (any(s%wildcards /= s%wildcards(1))) then
              valid = s%wildcards == minval(s%wildcards)
              call warning( "ambiguity resolved for", list%name, join(s%atoms), &
                            "( atom types", join(s%atom_types), ")" )
              aux = pack(s%itype,valid)
              ptr => typelist % first
              do itype = 1, typelist % count
                if (any(s%itype == itype)) then
                  if (any(aux == itype)) then
                    call writeln( "-->", typelist%name, join(ptr%id), "assigned" )
                  else
                    call writeln( "-->", typelist%name, join(ptr%id), "ignored" )
                  end if
                end if
                ptr => ptr % next
              end do
              s%itype = aux
              s%multiplicity = size(aux)
              total(s%mol) = total(s%mol) - size(valid) + count(valid)
            end if
          end if
        end associate
      end do

      ! Create a map for the indices of actually used types:
      imax = maxval([(maxval(structure(i)%itype),i=1,list%count)])
      aux = [(0,i=1,imax)]
      forall (i=1:list%count, nmols(structure(i)%mol) > 0) aux(structure(i)%itype) = 1
      allocate( type_map(count(aux == 1)) )
      type_map%index = pack([(i,i=1,imax)],aux == 1)

      ! Retrieve the parameters of each used type:
      do i = 1, size(type_map)
        ptr => typelist % point_to( type_map(i)%index )
        type_map(i) % types = join(ptr % id)
        if (models) then
          call split( ptr % params, narg, arg )
          type_map(i) % model = arg(1)
          type_map(i) % params = join(arg(2:narg))
        else
          type_map(i) % model = ""
          type_map(i) % params = ptr % params
        end if
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
  end subroutine tPlaymol_analyze_struct

  !=================================================================================================

  function tPlaymol_bodies_in_molecules( me ) result( body_in_mol )
    class(tPlaymol), intent(in) :: me
    integer                     :: body_in_mol(me%body_list%count)
    integer :: nbodies(me%molecules%N), ibody, imol
    character(sl) :: atom(1)
    type(Struc), pointer :: ptr
    nbodies = 0
    ptr => me % body_list % first
    do ibody = 1, me % body_list % count
      call split( ptr%params, imol, atom )
      imol = str2int(me % molecules % list % parameters( atom(1:1) ))
      nbodies(imol) = nbodies(imol) + 1
      body_in_mol(ibody) = nbodies(imol)
      ptr => ptr % next
    end do
  end function tPlaymol_bodies_in_molecules

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
!        if (narg /= 2) call error( "cannot guess coordinates" )
        R0 = arg(narg)
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
!        if (narg /= 2) call error( "cannot guess coordinates" )
        R0A = arg(narg)
        call split( me % bond_type_list % parameters( types(2:3) ), narg, arg )
!        if (narg /= 2) call error( "cannot guess coordinates" )
        R0B = arg(narg)
        call split( me % angle_type_list % parameters( types ), narg, arg )
!        if (narg /= 2) call error( "cannot guess coordinates" )
        Theta0 = arg(narg)
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
        X3 = str2real(R0A) - str2real(R0B)*cosine(str2real(Theta0))
        Y3 = str2real(R0B)*sine(str2real(Theta0))
        arg(1:4) = [atoms(3), real2str(X3), real2str(Y3), "0.0"]
        call me % molecules % xyz % add( 4, arg )
      end subroutine build_triatomic_molecule
      !-------------------------------------------------------------------------
  end subroutine tPlaymol_check_coordinates

  !=================================================================================================

  function velocities( N, seed, kT, mass ) result( V )
    integer,  intent(in) :: N, seed
    real(rb), intent(in) :: kT, mass(N)
    real(rb)             :: V(3,N)
    integer :: i
    real(rb) :: Vcm(3), stdev, Mtotal
    type(rng) :: random
    call random % init( seed )
    Mtotal = 0.0_rb
    Vcm = 0.0_rb
    do i = 1, N
      stdev = sqrt(kT/mass(i))
      V(:,i) = stdev*[random % normal(), random % normal(), random % normal()]
      Vcm = Vcm + mass(i)*V(:,i)
      Mtotal = Mtotal + mass(i)
    end do
    Vcm = Vcm/Mtotal
    forall (i=1:N) V(:,i) = V(:,i) - Vcm
    V = sqrt(kT*(3*N-3)/sum([(mass(i)*V(:,i)**2,i=1,N)]))*V
  end function velocities

  !=================================================================================================
  include "write_lammps.f90"
  !=================================================================================================
  include "write_lammpstrj.f90"
  !=================================================================================================
  include "write_emdee.f90"
  !=================================================================================================
  include "write_summary.f90"
  !=================================================================================================
  include "write_internals.f90"
  !=================================================================================================

end module mPlaymol
