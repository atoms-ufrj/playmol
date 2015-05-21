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

! TO DO: describe the command 'ionic_bond' in the manual.

module mPlaymol

use mGlobal
use mStruc
use mPackmol
use mBox

implicit none

type tPlaymol
  type(tBox) :: box
  integer  :: nmol = 0
  type(StrucList) :: atom_type_list = StrucList( name = "atom type", number = 1 )
  type(StrucList) :: bond_type_list = StrucList( name = "bond type", number = 2 )
  type(StrucList) :: angle_type_list = StrucList( name = "angle type", number = 3 )
  type(StrucList) :: dihedral_type_list = StrucList( name = "dihedral type", number = 4 )
  type(StrucList) :: improper_type_list = StrucList( name = "improper type", number = 4, & 
                                                     two_way = .false. )
  type(StrucList) :: mass_list = StrucList( name = "mass", number = 1 )
  type(StrucList) :: atom_list = StrucList( name = "atom", number = 1 )
  type(StrucList) :: charge_list = StrucList( name = "charge", number = 1 )
  type(StrucList) :: bond_list = StrucList( name = "bond", number = 2 )
  type(StrucList) :: ionic_bond_list = StrucList( name = "ionic bond", number = 2 )
  type(StrucList) :: angle_list = StrucList( name = "angle", number = 3 )
  type(StrucList) :: dihedral_list = StrucList( name = "dihedral", number = 4 )
  type(StrucList) :: extra_dihedral_list = StrucList( name = "extra dihedral", number = 4 )
  type(StrucList) :: improper_list = StrucList( name = "improper", number = 4 )
  type(StrucList) :: molecule_list = StrucList( name = "molecule", number = 1 )
  type(StrucList) :: coordinate_list = StrucList( name = "position of atom", number = 1 )
  type(StrucList) :: packmol_list = StrucList( name = "packmol", number = 1 )
  contains
    procedure :: read => tPlaymol_Read
    procedure :: write => tPlaymol_write
    procedure :: write_lammps => tPlaymol_write_lammps
    procedure :: write_lammpstrj => tPlaymol_write_lammpstrj
    procedure :: read_xyz => tPlaymol_read_xyz
    procedure :: write_xyz => tPlaymol_write_xyz
    procedure :: summarize => tPlaymol_summarize
    procedure :: fuse_molecules => tPlaymol_fuse_molecules
    procedure :: count_molecules => tPlaymol_count_molecules
    procedure :: update_structure => tPlaymol_update_structure
    procedure :: search_impropers => tPlaymol_search_impropers
    procedure :: get_types => tPlaymol_get_types
    procedure :: atoms_in_molecules => tPlaymol_atoms_in_molecules
    procedure :: check_coordinates => tPlaymol_check_coordinates
    procedure :: molecule_coordinates => tPlaymol_molecule_coordinates
end type tPlaymol

contains

  !=================================================================================================

  recursive subroutine tPlaymol_Read( me, unit, filename )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    character(*),    intent(in)    :: filename
    integer       :: narg
    character(sl) :: arg(50)
    call next_command( unit, narg, arg )
    do while (narg > 0)
      select case (trim(arg(1)))
        case ("prefix"); call prefix_command
        case ("suffix"); call suffix_command
        case ("box"); call box_command
        case ("atom_type"); call atom_type_command
        case ("bond_type"); call bond_type_command
        case ("angle_type"); call angle_type_command
        case ("dihedral_type"); call dihedral_type_command
        case ("improper_type"); call improper_type_command
        case ("mass"); call mass_command
        case ("atom"); call atom_command
        case ("charge"); call charge_command
        case ("bond"); call bond_command
        case ("ionic_bond"); call ionic_bond_command
        case ("extra_dihedral"); call extra_dihedral_command
        case ("improper"); call improper_command
        case ("xyz"); call xyz_command
        case ("align"); call align_command
        case ("write"); call write_command
        case ("include"); call include_command
        case ("packmol"); call packmol_command
        case ("reset"); call reset_command
        case ("shell"); call shell_command
        case ("quit")
          if ((narg == 2).and.(arg(2) == "all")) then
            call writeln( "Script", filename, "interrupted by 'quit all' command" )
            stop
          else if (narg == 1) then
            call writeln( "Script", filename, "interrupted by 'quit' command" )
            return
          else
            call error( "invalid quit command" )
          end if
        case default; call error( "unknown command", arg(1) )
      end select
      call next_command( unit, narg, arg )
    end do
    contains
      !---------------------------------------------------------------------------------------------
      subroutine prefix_command
        if (narg < 3) call error( "invalid prefix command" )
        call writeln( "Defining new prefix for ", arg(2), ":", arg(3) )
        if (arg(3) == "none") arg(3) = ""
        select case (arg(2))
          case ("types"); me % atom_type_list % prefix = arg(3)
          case ("atoms"); me % atom_list % prefix = arg(3)
          case default; call error( "prefix keyword must be 'types' or 'atoms'" )
        end select
        if (has_macros(arg(3))) call error( "invalid prefix", arg(3) )
      end subroutine prefix_command
      !---------------------------------------------------------------------------------------------
      subroutine suffix_command
        if (narg < 3) call error( "invalid suffix command" )
        call writeln( "Defining new suffix for ", arg(2), ":", arg(3) )
        if (arg(3) == "none") arg(3) = ""
        select case (arg(2))
          case ("types"); me % atom_type_list % suffix = arg(3)
          case ("atoms"); me % atom_list % suffix = arg(3)
          case default; call error( "suffix keyword must be 'types' or 'atoms'" )
        end select
        if (has_macros(arg(3))) call error( "invalid suffix", arg(3) )
      end subroutine suffix_command
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
          case default
             call error( "invalid box command" )
        end select
        call me % box % define( arg(2), scalar, vector )
      end subroutine box_command
      !---------------------------------------------------------------------------------------------
      subroutine atom_type_command
        call me % atom_type_list % add( narg-1, arg(2:narg) )
        if (has_macros(arg(2))) call error( "invalid atom_type name" )
      end subroutine atom_type_command
      !---------------------------------------------------------------------------------------------
      subroutine bond_type_command
        call me % bond_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine bond_type_command
      !---------------------------------------------------------------------------------------------
      subroutine angle_type_command
        call me % angle_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine angle_type_command
      !---------------------------------------------------------------------------------------------
      subroutine dihedral_type_command
        call me % dihedral_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine dihedral_type_command
      !---------------------------------------------------------------------------------------------
      subroutine improper_type_command
        call me % improper_type_list % add( narg-1, arg(2:narg), me % atom_type_list, .true. )
      end subroutine improper_type_command
      !---------------------------------------------------------------------------------------------
      subroutine mass_command
        call me % mass_list % add( narg-1, arg(2:narg) )
      end subroutine mass_command
      !---------------------------------------------------------------------------------------------
      subroutine atom_command
        type(Struc), pointer :: type_ptr
        if ((narg < 3).or.(narg > 4)) call error( "invalid atom command" )
        arg(3) = trim(me % atom_type_list % prefix) // trim(arg(3)) // me % atom_type_list % suffix
        call me % atom_list % add( 2, arg(2:3) )
        if (has_macros(arg(2))) call error( "invalid atom name" )
        if (has_macros(arg(3))) call error( "invalid atom type name" )
        call me % atom_type_list % search( arg(3:3), type_ptr )
        if (associated(type_ptr)) then
          type_ptr % used = .true.
        else
          call error( "atom type", arg(3), "required, but not found")
        end if
        if (.not.me % mass_list % find(arg(3:3))) call error( "atom type",arg(3), "has no mass" )
        me%nmol = me%nmol + 1
        arg(3) = int2str( me%nmol )
        call me % molecule_list % add( 2, arg(2:3) )
        if (narg == 4) then
          arg(3) = arg(4)
          call me % charge_list % add( 2, arg(2:3) )
        end if
      end subroutine atom_command
      !---------------------------------------------------------------------------------------------
      subroutine charge_command
        call me % charge_list % add( narg-1, arg(2:narg) )
      end subroutine charge_command
      !---------------------------------------------------------------------------------------------
      subroutine add_bond( atom )
        character(sl), intent(inout) :: atom(2)
        call me % bond_list % add( 2, atom, me % atom_list )
        call me % bond_list % handle( atom, me % atom_list, me % bond_type_list, 2 )
        if (atom(1) == atom(2)) call error( "atom", atom(1), "cannot bind to itself" )
        call me % fuse_molecules( atom )
        call me % update_structure()
      end subroutine add_bond
      !---------------------------------------------------------------------------------------------
      subroutine bond_command
        integer :: i
        character(sl) :: central
        if (narg < 3) call error( "invalid bond command" )
        central = arg(2)
        do i = 3, narg
          arg(2) = central
          arg(3) = arg(i)
          call add_bond( arg(2:3) )
        end do
      end subroutine bond_command
      !---------------------------------------------------------------------------------------------
      subroutine add_ionic_bond( atom )
        character(sl), intent(inout) :: atom(2)
        call me % ionic_bond_list % add( 2, atom, me % atom_list )
        if (atom(1) == atom(2)) call error( "atom", atom(1), "cannot bind to itself" )
        call me % fuse_molecules( atom )
      end subroutine add_ionic_bond
      !---------------------------------------------------------------------------------------------
      subroutine ionic_bond_command
        integer :: i
        character(sl) :: central
        if (narg < 3) call error( "invalid ionic_bond command" )
        central = arg(2)
        do i = 3, narg
          arg(2) = central
          arg(3) = arg(i)
          call add_ionic_bond( arg(2:3) )
        end do
      end subroutine ionic_bond_command
      !---------------------------------------------------------------------------------------------
      subroutine extra_dihedral_command
        integer :: i, imol, jmol
        call me % extra_dihedral_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 5) call error( "invalid extra_dihedral command")
        call me % extra_dihedral_list % handle( arg(2:5), me%atom_list, me%dihedral_type_list, 2 )
        imol = str2int(me % molecule_list % parameters( arg(2:2) ))
        do i = 3, 5
          jmol = str2int(me % molecule_list % parameters( arg(i:i) ))
          if (jmol /= imol) call error( "atoms", join(arg(2:5)), "are not in the same molecule" )
        end do
      end subroutine extra_dihedral_command
      !---------------------------------------------------------------------------------------------
      subroutine improper_command
        integer :: i, imol, jmol
        if (narg == 2) then
          if (arg(2) /= "search") call error( "invalid improper command" )
          call me % search_impropers()
        else if (narg == 5) then
          call me % improper_list % add( narg-1, arg(2:narg), me % atom_list )
          if (narg /= 5) call error( "invalid improper command")
          call me % improper_list % handle( arg(2:5), me%atom_list, me%improper_type_list, 2 )
          imol = str2int(me % molecule_list % parameters( arg(2:2) ))
          do i = 3, 5
            jmol = str2int(me % molecule_list % parameters( arg(i:i) ))
            if (jmol /= imol) call error( "atoms", join(arg(2:5)), "are not in the same molecule" )
          end do
        else
          call error( "invalid improper command" )
        end if
      end subroutine improper_command
      !---------------------------------------------------------------------------------------------
      subroutine xyz_command
        logical :: file_exists
        integer :: xyz
        if (narg == 1) then
          call writeln("Reading xyz data...")
          call me % read_xyz( unit )
        else
          inquire( file = arg(2), exist = file_exists )
          if (.not.file_exists) call error( "file", arg(2), "not found" )
          open( newunit = xyz, file = arg(2), status = "old" )
          call writeln("Reading xyz data from file", trim(arg(2))//"..." )
          call me % read_xyz( xyz )
          close( xyz )
        end if
      end subroutine xyz_command
      !---------------------------------------------------------------------------------------------
      subroutine align_command
        use mAlign
        integer :: i, imol, natoms(me%nmol), N
        real(rb), allocatable :: Mass(:), Coord(:,:)
        integer :: axis(3)
        character :: dir(3) = ["x","y","z"]
        real(rb) :: lb(3), ub(3)
        if (narg /= 4) call error( "invalid align command" )
        imol = str2int(arg(2))
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
        if (imol < 1) call error( "wrong molecule index" )
        if (imol > me%nmol) call error( "last molecule is", int2str(me%nmol) )
        call me % check_coordinates( imol )
        call writeln( "Aligning molecule", arg(2), "with axes", arg(3), "and", arg(4) )
        natoms = me % atoms_in_molecules()
        N = natoms(imol)
        allocate( Mass(N), Coord(3,N) )
        Mass = 1.0_rb
        call me % molecule_coordinates( imol, N, Coord, 1 )
        call align_molecule( N, Mass, Coord, axis )
        call me % molecule_coordinates( imol, N, Coord, 2 )
        lb = minval(Coord,dim=2)
        ub = maxval(Coord,dim=2)
        do i = 1, 3
          call writeln( "New bounds in direction", dir(i)//":", join(real2str([lb(i),ub(i)])) )
        end do
        do i = 1, 3
          call writeln( "New length in direction", dir(i)//":", real2str(ub(i)-lb(i)) )
        end do
        deallocate( Mass, Coord )
      end subroutine align_command
      !---------------------------------------------------------------------------------------------
      subroutine write_command
        integer :: unit
        if ((narg < 2).or.(narg > 4)) call error( "invalid write command" )
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
        if (any(lists ==  1)) call me % atom_type_list % destroy
        if (any(lists ==  2)) call me % bond_type_list % destroy
        if (any(lists ==  3)) call me % angle_type_list % destroy
        if (any(lists ==  4)) call me % dihedral_type_list % destroy
        if (any(lists ==  5)) call me % improper_type_list % destroy
        if (any(lists ==  6)) call me % mass_list % destroy
        if (any(lists ==  7)) call me % atom_list % destroy
        if (any(lists ==  8)) call me % charge_list % destroy
        if (any(lists ==  9)) call me % bond_list % destroy
        if (any(lists == 10)) call me % angle_list % destroy
        if (any(lists == 11)) call me % dihedral_list % destroy
        if (any(lists == 12)) call me % improper_list % destroy
        if (any(lists == 13)) call me % molecule_list % destroy
        if (any(lists == 14)) call me % coordinate_list % destroy
        if (any(lists == 15)) call me % packmol_list % destroy
      end subroutine reset_lists
      !---------------------------------------------------------------------------------------------
      subroutine reset_command
        integer :: i
        if (narg /= 2) call error( "invalid reset command" )
        select case (arg(2))
          case ("all"); call reset_lists( [(i,i=1,15)] )
          case ("atoms"); call reset_lists( [(i,i=7,15)] )
          case ("bonds"); call reset_lists( [9,10,11,13,14,15] )
          case ("bond_types"); call me % bond_type_list % destroy
          case ("angle_types"); call me % angle_type_list % destroy
          case ("dihedral_types"); call me % dihedral_type_list % destroy
          case ("improper_types"); call me % improper_type_list % destroy
          case ("impropers"); call me % improper_list % destroy
          case ("xyz"); call me % coordinate_list % destroy
          case ("packmol"); call me % packmol_list % destroy
          case default; call error( "invalid reset command" )
        end select
      end subroutine reset_command
      !---------------------------------------------------------------------------------------------
      subroutine shell_command
        integer :: stat
        call writeln( "Executing shell command: ", join(arg(2:)) )
        call execute_command_line( join(arg(2:)), exitstat = stat )
        if (stat /= 0) call error( "unsuccessful shell command" )
      end subroutine shell_command
      !---------------------------------------------------------------------------------------------
      subroutine packmol_command
        integer :: iarg, imol, nopts, i, n
        real(rb) :: pos(3), mass(me%nmol)
        character(sl) :: action
        action = ""
        if (narg < 3) call error( "invalid packmol command" )
        iarg = 2
        do while (iarg < narg)
          select case (arg(iarg))

            case ("seed")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol seed to", arg(iarg+1) )
              seed = str2int(arg(iarg+1))
              iarg = iarg + 2

            case ("tolerance")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol tolerance to", arg(iarg+1) )
              tolerance = str2real(arg(iarg+1))
              iarg = iarg + 2

            case ("nloops")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol nloops parameter to", arg(iarg+1) )
              nloops = str2int(arg(iarg+1))
              if (nloops <= 0) call error( "invalid parameter value" )
              iarg = iarg + 2

            case ("retry")
              if (narg < iarg+1) call error( "invalid packmol command" )
              call writeln( "Setting packmol retry parameter to", arg(iarg+1) )
              change = str2real(arg(iarg+1))
              if ((change <= 0.0_rb).or.(change > 1.0_rb)) call error( "invalid parameter value" )
              iarg = iarg + 2

            case ("fix","move","copy","pack")
              if (narg < iarg+2) call error( "invalid packmol command" )
              select case (arg(iarg))
                case ("move","fix"); nopts = 3
                case ("copy","pack"); nopts = 1
                case default; call error( "invalid packmol command" )
              end select
              call me % packmol_list % add( nopts+2, arg(iarg:iarg+1+nopts), repeatable = .true. )
              imol = str2int(arg(iarg+1))
              if ((imol < 1).or.(imol > me%nmol)) call error( "last molecule is", int2str(me%nmol) )
              call me % check_coordinates( imol )
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
              iarg = iarg + 2

            case default
              call error( "invalid packmol command: unknown keyword", arg(iarg) )

          end select
        end do
        if (action /= "") then
          call me % count_molecules( mass = mass )
          call me % box % compute( packmol_total_mass( me % packmol_list, mass ) )
          call writeln( "Box lengths are ", join(real2str(me % box % length)))
          call run_packmol( me % packmol_list, me % molecule_list, me % coordinate_list, &
                            me % nmol, me % atoms_in_molecules(), me % box % length, &
                            seed, tolerance, action )
        end if
      end subroutine packmol_command
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_Read

  !=================================================================================================

  subroutine tPlaymol_read_xyz( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: N, i, narg, imol, iatom
    character(sl) :: arg(5), line, catom
    integer :: natoms(me % nmol)
    logical :: new_molecule
    type(Struc), pointer :: atom
    character(sl), allocatable :: prev(:)
    natoms = me % atoms_in_molecules()
    allocate( prev(maxval(natoms)) )
    call next_command( unit, narg, arg )
    if (narg > 0) then
      call writeln( "Number of coordinates: ", arg(1) )
      N = str2int( arg(1) )
      read(unit,'(A'//csl//')') line
      new_molecule = .true.
      do i = 1, N
        read(unit,'(A'//csl//')') line
        call split( line, narg, arg )
        if (narg /= 4) call error( "invalid xyz file format" )
        catom = trim(me % atom_list % prefix)//trim(arg(1))//trim(me % atom_list % suffix)
        call me % molecule_list % search( [catom], atom )
        if (associated(atom)) then
          if (new_molecule) then
            imol = str2int(atom%params)
            iatom = 1
            call writeln( "Reading", int2str(natoms(imol)), &
                          "atom coordinates of molecule", trim(int2str(imol))//":" )
          else if (trim(atom%params) /= trim(int2str(imol))) then
            call error( "atom", catom, "does not belong to molecule", int2str(imol) )
          else if (any(str_find([catom],prev(1:iatom)) > 0)) then
            call error( "coordinates of atom", catom, "have already been defined" )
          else
            iatom = iatom + 1
          end if
        end if
        call me % coordinate_list % add( narg, arg, me % atom_list, repeatable = .true. )
        if (.not.all(is_real(arg(2:4)))) call error( "invalid coordinate" )
        prev(iatom) = arg(1)
        new_molecule = iatom == natoms(imol)
      end do
      if (.not.new_molecule) then
        call error( "coordinates of molecule", int2str(imol), "are incomplete" )
      end if
    end if
  end subroutine tPlaymol_read_xyz

  !=================================================================================================

  subroutine tPlaymol_write_xyz( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    type(Struc), pointer :: ptr
    write(unit,'(A)') trim(int2str(me % coordinate_list % count()))
    write(unit,'("# Generated by Playmol on ",A)') trim(now())
    ptr => me % coordinate_list % first
    do while (associated(ptr))
      write(unit,'(A)') trim(ptr % id(1))//" "//trim(ptr % params)
      ptr => ptr % next
    end do
  end subroutine tPlaymol_write_xyz

  !=================================================================================================

  subroutine tPlaymol_fuse_molecules( me, atom )
    class(tPlaymol),    intent(inout) :: me
    character(sl),   intent(in)    :: atom(2)
    integer :: i, mol(2), imin, imax
    type(Struc), pointer :: ptr
    do i = 1, 2
      call me % molecule_list % search( atom(i:i), ptr )
      mol(i) = str2int( ptr % params )
    end do
    imin = minval(mol)
    imax = maxval(mol)
    if (imin == imax) then
      call writeln( "Closing cycle(s) in molecule", int2str(imin) )
    else
      call writeln( "Fusing molecule", int2str(imax), "to molecule", int2str(imin) )
      call rename_molecule( trim(int2str(imax)), trim(int2str(imin)) )
      call writeln( "Lowering indices of molecules", int2str(imax+1), "to", int2str(me%nmol) )
      do i = imax+1, me%nmol
        call rename_molecule( trim(int2str(i)), trim(int2str(i-1)) )
      end do
      me%nmol = me%nmol - 1
    end if	
    contains
      !---------------------------------------------------------------------------------------------
      subroutine rename_molecule( old, new )
        character(*), intent(in) :: old, new
        ptr => me % molecule_list % first
        do while (associated(ptr))
          if (ptr % params == old) ptr % params = new
          ptr => ptr % next
        end do
      end subroutine rename_molecule
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_fuse_molecules

  !=================================================================================================

  subroutine tPlaymol_count_molecules( me, number, mass, charge )
    class(tPlaymol),  intent(in)            :: me
    integer,       intent(out), optional :: number(me % nmol)
    real(rb),      intent(out), optional :: mass(me % nmol), charge(me % nmol)
    integer  :: imol, iatom, natoms(me % nmol)
    character(sl) :: atom(1), atom_type(1)
    type(Struc), pointer :: ptr
    natoms = me % atoms_in_molecules()
    if (present(number)) then
      number = 0
      ptr => me % coordinate_list % first
      do while (associated(ptr))
        imol = str2int(me % molecule_list % parameters( ptr%id ) )
        number(imol) = number(imol) + 1
        do iatom = 1, natoms(imol)
          ptr => ptr % next
        end do
      end do
    end if
    if (present(mass)) then
      mass = 0.0_rb
      ptr => me % molecule_list % first
      do while (associated(ptr))
        atom = ptr % id
        imol = str2int( ptr % params )
        call me % get_types( atom, atom_type )
        mass(imol) = mass(imol) + str2real( me % mass_list % parameters( atom_type ) )
        ptr => ptr % next
      end do
    end if
    if (present(charge)) then
      charge = 0.0_rb
      if (associated(me % charge_list % first)) then
        ptr => me % molecule_list % first
        do while (associated(ptr))
          atom = ptr % id
          imol = str2int( ptr % params )
          call me % get_types( atom, atom_type )
          charge(imol) = charge(imol) + str2real( me % charge_list % parameters( atom ) )
          ptr => ptr % next
        end do
      end if
    end if
  end subroutine tPlaymol_count_molecules

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
    type(Struc), pointer :: b1, b2, b3
    character(sl) :: atom(3,2), central
    call writeln( "Searching for impropers..." )
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

  function tPlaymol_atoms_in_molecules( me ) result( natoms )
    class(tPlaymol), intent(in) :: me
    integer                  :: natoms(me%nmol)
    type(Struc), pointer :: atom
    integer :: i
    natoms = 0
    atom => me % molecule_list % first
    do while (associated(atom))
      i = str2int(atom%params)
      natoms(i) = natoms(i) + 1
      atom => atom % next
    end do
  end function tPlaymol_atoms_in_molecules

  !=================================================================================================

  subroutine tPlaymol_write( me, unit )
    class(tPlaymol), intent(in) :: me
    integer,      intent(in) :: unit
    integer :: N
    character(sl) :: CN
    type(Struc), pointer :: current
    call me % atom_type_list % print( unit )
    call me % mass_list % print( unit )
    call me % bond_type_list % print( unit )
    call me % angle_type_list % print( unit )
    call me % dihedral_type_list % print( unit )
    call me % atom_list % print( unit )
    call me % charge_list % print( unit )
    call me % bond_list % print( unit )
    call me % angle_list % print( unit, comment = .true. )
    call me % dihedral_list % print( unit, comment = .true. )
    call me % extra_dihedral_list % print( unit )
    call me % improper_list % print( unit )
    N = me % coordinate_list % count()
    if (N > 0) then
      write(CN,*) N
      write(unit,'("xyz"/,A,/,"# atom x y z")') trim(adjustl(CN))
      current => me % coordinate_list % first
      do while (associated(current))
        write(unit,'(A,X,A)') trim(current % id(1)), trim(current % params)
        current => current % next
      end do
    end if
  end subroutine tPlaymol_write

  !=================================================================================================

  subroutine tPlaymol_summarize( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer :: imol, molcount(me % nmol), natoms(me % nmol)
    real(rb) :: mass(me % nmol), charge(me % nmol)
    type(Struc), pointer :: ptr
    character(sl) :: atoms
    natoms = me % atoms_in_molecules()
    write(unit,'(A)') repeat("-",80)
    write(unit,'(A)') "SUMMARY"
    write(unit,'(A)') repeat("-",80)
    write(unit,'(A)') "Specified:"
    call flush( "atom type", me % atom_type_list % count() )
    call flush( "bond type", me % bond_type_list % count() )
    call flush( "angle type", me % angle_type_list % count() )
    call flush( "dihedral type", me % dihedral_type_list % count() )
    call flush( "improper type", me % improper_type_list % count() )
    call flush( "atom", me % atom_list % count() )
    call flush( "bond", me % bond_list % count() )
    call flush( "ionic bond", me % ionic_bond_list % count() )
    call flush( "extra dihedral", me % extra_dihedral_list % count() )
    call flush( "improper", me % improper_list % count() )
    write(unit,'(/,A)') "Detected:"
    call flush( "angle", me % angle_list % count() )
    call flush( "dihedral", me % dihedral_list % count() )
    call flush( "molecule", me % nmol )
    write(unit,'(/,A)') "Effectively used:"
    call flush( "atom type", me % atom_type_list % count_used() )
    call flush( "bond type", me % bond_type_list % count_used() )
    call flush( "angle type", me % angle_type_list % count_used() )
    call flush( "dihedral type", me % dihedral_type_list % count_used() )
    call flush( "improper type", me % improper_type_list % count_used() )
    call me % count_molecules( molcount, mass, charge )
    do imol = 1, me % nmol
      write(unit,'(/,"Molecule[",A,"]:")') trim(int2str(imol))
      write(unit,'("- Amount: ",A)') trim(int2str(molcount(imol)))
      write(unit,'("- Mass: ",A)') trim(real2str(mass(imol)))
      write(unit,'("- Charge: ",A)') trim(real2str(charge(imol)))
      write(unit,'("- Number of atoms: ",A)') trim(int2str(natoms(imol)))
      atoms = "- Atoms:"
      ptr => me % molecule_list % first
      do while (associated(ptr))
        if (str2int(ptr % params) == imol) then
          if (len_trim(atoms) + len_trim(ptr % id(1)) > 79) then
            write(unit,'(A)') trim(atoms)
            atoms = ""
          end if
          atoms = trim(atoms)//" "//ptr % id(1)
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
      write(unit,'("- Box density: ",A)') trim(real2str(me % box % density))
    end if
    write(unit,'("- Residual charge: ",A)') trim(real2str(sum(molcount*charge)))
    write(unit,'(A)') repeat("-",80)
    contains
      subroutine flush( title, amount )
        character(*), intent(in) :: title
        integer,      intent(in) :: amount
        if (amount == 1) then
          write(unit,'(I5,X,A,".")') amount, title
        else if (amount > 1) then
          write(unit,'(I5,X,A,"s.")') amount, title
        end if
      end subroutine flush
  end subroutine tPlaymol_summarize

  !=================================================================================================

  subroutine tPlaymol_write_lammps( me, unit )
    class(tPlaymol),  intent(inout) :: me
    integer,          intent(in)    :: unit
    integer :: nb, na, nd, ni, natoms(me % nmol)
    natoms = me % atoms_in_molecules()
    if (associated(me % dihedral_list % first)) then
      me % dihedral_list % last % next => me % extra_dihedral_list % first
    end if
    write(unit,'("LAMMPS data file",/,"# Generated by Playmol on ",A,/)') trim(now())
    call write_count( me % atom_type_list % count_used(), "atom types" )
    call write_count( me % bond_type_list % count_used(), "bond types" )
    call write_count( me % angle_type_list % count_used(), "angle types" )
    call write_count( me % dihedral_type_list % count_used(), "dihedral types" )
    call write_count( me % improper_type_list % count_used(), "improper types" )
    write(unit,'()')
    call write_count( me % coordinate_list % count(), "atoms" )
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
    if (associated(me % dihedral_list % first)) then
      me % dihedral_list % last % next => null()
    end if
    contains
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
        integer :: i, molcount(me % nmol)
        real(rb) :: mass(me % nmol)
        character :: dir(3) = ["x","y","z"]
        character(sl) :: limits
        call me % count_molecules( molcount, mass )
        call me % box % compute( sum(molcount*mass) )
        write(unit,'()')
        do i = 1, 3
          limits = join(real2str( me%box%length(i)*[-0.5_rb,+0.5_rb] ))
          write(unit,'(A," ",A,"lo ",A,"hi")') trim(limits), dir(i), dir(i)
        end do
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
        current => me % coordinate_list % first
        iatom = 0
        jmol = 0
        do while (associated(current))
          jmol = jmol + 1
          imol = str2int(me % molecule_list % parameters( current % id ) )
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
        coord => me % coordinate_list % first
        if (.not.present(count)) then
          call writeln( "Writing down structures:", name )
          write(unit,'(/,A,/)') name
        end if
        allocate(atom_id(maxval(natoms)))
        last_atom = 0
        do while (associated(coord))
          imol = str2int(me % molecule_list % parameters( coord%id ) )
          do i = 1, natoms(imol)
            atom_id(i) = coord % id(1)
            coord => coord % next
          end do
          struct => list % first
          do while (associated(struct))
            jmol = str2int(me % molecule_list % parameters(struct%id(1:1)))
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
  end subroutine tPlaymol_write_lammps

  !=================================================================================================

  subroutine tPlaymol_write_lammpstrj( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    integer :: molcount(me % nmol), iatom, itype, i, imol, jmol, narg, natoms(me % nmol)
    real(rb) :: mass(me % nmol)
    character(sl) :: limits, arg(1)
    type(Struc), pointer :: current, atom_type
    logical :: found
    if (.not.me % box % exists()) call error( "simulation box has not been defined" )
    write(unit,'("ITEM: TIMESTEP",/,"0")')
    write(unit,'("ITEM: NUMBER OF ATOMS")')
    write(unit,'(A)') trim(int2str(me % coordinate_list % count()))
    write(unit,'("ITEM: BOX BOUNDS pp pp pp")')
    call me % count_molecules( molcount, mass )
    call me % box % compute( sum(molcount*mass) )
    do i = 1, 3
      limits = join(real2str( me%box%length(i)*[-0.5_rb,+0.5_rb] ))
      write(unit,'(A)') trim(limits)
    end do
    write(unit,'("ITEM: ATOMS id mol type x y z ix iy iz ")')
    natoms = me % atoms_in_molecules()
    current => me % coordinate_list % first
    iatom = 0
    jmol = 0
    do while (associated(current))
      jmol = jmol + 1
      imol = str2int(me % molecule_list % parameters( current % id ) )
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
        write(unit,'(3(A,X),"0 0 0")') trim(join(int2str([iatom,jmol,itype]))), trim(current%params)
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
    coord => me % coordinate_list % first
    found = .false.
    do while (associated(coord).and.(.not.found))
      found = str2int(me % molecule_list % parameters( coord%id ) ) == imol
      coord => coord % next
    end do
    if (.not.found) then
      call warning( "No predefined coordinates for molecule", int2str(imol) )
      atom => me % molecule_list % first
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
        call writeln( "Molecule", int2str(imol), "is monoatomic:" )
        arg = [first%params, "0.0", "0.0", "0.0"]
        call me % coordinate_list % add( 4, arg )
      end subroutine build_monoatomic_molecule
      !-------------------------------------------------------------------------
      subroutine build_diatomic_molecule
        integer :: narg
        character(sl) :: atoms(2), types(2), arg(10), R0
        call writeln( "Molecule", int2str(imol), "is diatomic:" )
        atoms = [first % params, first % next % params]
        call me % get_types( atoms, types )
        call split( me % bond_type_list % parameters( types ), narg, arg )
        if (narg /= 2) call error( "cannot guess coordinates" )
        R0 = arg(2)
        if (.not.is_real(R0)) call error( "cannot guess coordinates" )
        call writeln( "Considering harmonic potential with R0 = ", R0 )
        call warning( "if this is wrong, please provide coordinates manually via xyz" )
        arg(1:4) = [atoms(1), "0.0", "0.0", "0.0"]
        call me % coordinate_list % add( 4, arg )
        arg(1:4) = [atoms(2), R0, "0.0", "0.0"]
        call me % coordinate_list % add( 4, arg )
      end subroutine build_diatomic_molecule
      !-------------------------------------------------------------------------
      subroutine build_triatomic_molecule
        integer :: narg
        character(sl) :: atoms(3), types(3), arg(10), R0A, R0B, Theta0
        real(rb) :: X3, Y3
        call writeln( "Molecule", int2str(imol), "is triatomic:" )
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
          call writeln( "Considering harmonic potentials:" )
          call writeln( "R0(",atoms(1),"-",atoms(2),") = ", R0A )
          call writeln( "R0(",atoms(2),"-",atoms(3),") = ", R0B )
          call writeln( "Theta0(",atoms(1),"-",atoms(2),"-",atoms(3),") = ", Theta0, "degrees" )
          call warning( "if this is wrong, please provide coordinates manually via xyz" )
        else
          call error( "cannot guess coordinates" )
        end if
        arg(1:4) = [atoms(1), "0.0", "0.0", "0.0"]
        call me % coordinate_list % add( 4, arg )
        arg(1:4) = [atoms(2), R0A, "0.0", "0.0"]
        call me % coordinate_list % add( 4, arg )
        X3 = str2real(R0A) - str2real(R0B)*cos(0.0174532925_rb*str2real(Theta0))
        Y3 = str2real(R0B)*sin(0.0174532925_rb*str2real(Theta0))
        arg(1:4) = [atoms(3), real2str(X3), real2str(Y3), "0.0"]
        call me % coordinate_list % add( 4, arg )
      end subroutine build_triatomic_molecule
      !-------------------------------------------------------------------------
  end subroutine tPlaymol_check_coordinates
  !=================================================================================================
  subroutine tPlaymol_molecule_coordinates( me, imol, N, Coord, option )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: imol, N, option
    real(rb),        intent(inout) :: Coord(3,N)
    type(Struc), pointer :: current
    logical :: found
    integer :: i, narg
    character(sl) :: arg(3)
    current => me % coordinate_list % first
    found = .false.
    do while (associated(current).and.(.not.found))
      found = str2int(me % molecule_list % parameters( current%id )) == imol
      if (.not.found) current => current % next
    end do
    do i = 1, N
      if (option == 1) then ! Retrieve coordinates and masses:
        call split( current % params, narg, arg )
        Coord(:,i) = [str2real(arg(1)), str2real(arg(2)), str2real(arg(3))]
      else ! Set coordinates:
        current % params = join(real2str(Coord(:,i)))
      end if
      current => current % next
    end do
  end subroutine tPlaymol_molecule_coordinates
  !=================================================================================================
end module mPlaymol
