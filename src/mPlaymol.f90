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
  type(StrucList) :: dihedral_type_list = StrucList( name = "dihedral type", number = 4, two_way = .false. )
  type(StrucList) :: improper_type_list = StrucList( name = "improper type", number = 4, two_way = .false. )
  type(StrucList) :: mass_list = StrucList( name = "mass", number = 1 )
  type(StrucList) :: pair_list = StrucList( name = "pair", number = 2 )
  type(StrucList) :: atom_list = StrucList( name = "atom", number = 1 )
  type(StrucList) :: charge_list = StrucList( name = "charge", number = 1 )
  type(StrucList) :: bond_list = StrucList( name = "bond", number = 2 )
  type(StrucList) :: angle_list = StrucList( name = "angle", number = 3 )
  type(StrucList) :: dihedral_list = StrucList( name = "dihedral", number = 4 )
  type(StrucList) :: improper_list = StrucList( name = "improper", number = 4 )
  type(StrucList) :: molecule_list = StrucList( name = "molecule", number = 1 )
  type(StrucList) :: coordinate_list = StrucList( name = "position of atom", number = 1 )
  type(StrucList) :: packmol_list = StrucList( name = "packmol molecule", number = 1 )
  contains
    procedure :: read => tPlaymol_Read
    procedure :: write => tPlaymol_write
    procedure :: write_lammps => tPlaymol_write_lammps
    procedure :: read_xyz => tPlaymol_read_xyz
    procedure :: write_xyz => tPlaymol_write_xyz
    procedure :: summarize => tPlaymol_summarize
    procedure :: fuse_molecules => tPlaymol_fuse_molecules
    procedure :: count_molecules => tPlaymol_count_molecules
    procedure :: update_structure => tPlaymol_update_structure
    procedure :: get_types => tPlaymol_get_types
    procedure :: check_types => tPlaymol_check_types
    procedure :: atoms_in_molecules => tPlaymol_atoms_in_molecules
end type tPlaymol

contains

  !=================================================================================================

  recursive subroutine tPlaymol_Read( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: narg
    character(sl) :: arg(50)
    call next_command( unit, narg, arg )
    do while (narg > 0)
      select case (trim(arg(1)))
        case ("prefix"); call prefix_command
        case ("box"); call box_command
        case ("atom_type"); call atom_type_command
        case ("bond_type"); call bond_type_command
        case ("angle_type"); call angle_type_command
        case ("dihedral_type"); call dihedral_type_command
        case ("improper_type"); call improper_type_command
        case ("mass"); call mass_command
!        case ("pair"); call pair_command
        case ("atom"); call atom_command
        case ("charge"); call charge_command
        case ("bond"); call bond_command
        case ("improper"); call improper_command
        case ("xyz"); call xyz_command
        case ("write"); call write_command
        case ("include"); call include_command
        case ("packmol"); call packmol_command
        case ("quit"); stop "interrupted by quit command"
        case default; call error( "unknown command", arg(1) )
      end select
      call next_command( unit, narg, arg )
    end do
    contains
      !---------------------------------------------------------------------------------------------
      subroutine prefix_command
        if (narg < 3) call error( "invalid prefix command" )
        if (arg(3) == "none") arg(3) = ""
        select case (arg(2))
          case ("types"); me % atom_type_list % prefix = arg(3)
          case ("atoms"); me % atom_list % prefix = arg(3)
          case default; call error( "prefix keyword must be 'types' or 'atoms'" )
        end select
        if (has_macros(arg(3))) call error( "invalid prefix", arg(3) )
      end subroutine prefix_command
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
      subroutine pair_command
        call me % pair_list % add( narg-1, arg(2:narg), me % atom_type_list )
      end subroutine pair_command
      !---------------------------------------------------------------------------------------------
      subroutine atom_command
        if (narg >= 3) arg(3) = trim(me % atom_type_list % prefix) // arg(3)
        call me % atom_list % add( narg-1, arg(2:narg) )
        if (has_macros(arg(2))) call error( "invalid atom name" )
        call me % check_types( arg(2:2), me % atom_type_list )
        if (.not.me % mass_list % find(arg(3:3))) call error( "atom type",arg(3), "has no mass" )
        me%nmol = me%nmol + 1
        arg(3) = int2str( me%nmol )
        call me % molecule_list % add( 2, arg(2:3) )
      end subroutine atom_command
      !---------------------------------------------------------------------------------------------
      subroutine charge_command
        call me % charge_list % add( narg-1, arg(2:narg), me % atom_list )
      end subroutine charge_command
      !---------------------------------------------------------------------------------------------
      subroutine bond_command
        call me % bond_list % add( narg-1, arg(2:narg), me % atom_list )
        if (arg(2) == arg(3)) call error( "atom", arg(2), "cannot bind to itself" )
        call me % check_types( arg(2:3), me % bond_type_list )
        call me % fuse_molecules( arg(2:3) )
        call me % update_structure()
      end subroutine bond_command
      !---------------------------------------------------------------------------------------------
      subroutine improper_command
        integer :: i, imol, jmol
        call me % improper_list % add( narg-1, arg(2:narg), me % atom_list )
        if (narg /= 5) call error( "invalid improper command")
        call me % check_types( arg(2:5), me % improper_type_list )
        imol = str2int(me % molecule_list % parameters( arg(2:2) ))
        do i = 3, 5
          jmol = str2int(me % molecule_list % parameters( arg(i:i) ))
          if (jmol /= imol) call error( "atoms", join(arg(2:5)), "are not in the same molecule" )
        end do
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
      subroutine write_command
        integer :: unit
        if ((narg < 2).or.(narg > 4)) call error( "invalid write command" )
        if (.not.any(arg(2) == [character(7)::"playmol","lammps","summary","xyz"]) ) then
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
        call me % Read( input )
        close(input)
      end subroutine include_command
      !---------------------------------------------------------------------------------------------
      subroutine packmol_command
        integer :: seed, iarg, imol, last, nopts, i, n, molcount(me%nmol)
        real(rb) :: tol, pos(3), mass(me%nmol)
        character(sl) :: action
        if (.not. me % box % exists()) then
          call error( "packmol run command requires previous box definition" )
        end if
        if (narg < 8) call error( "invalid packmol command" )
        action = arg(narg)
        if (all(["setup  ","execute","persist"] /= action)) then
          call error( "invalid packmol command - last argument must be an action" )
        end if
        call writeln( "Configuring packmol with tolerance", arg(2), "and seed", arg(3) )
        tol = str2real(arg(2))
        seed = str2int(arg(3))
        if (arg(4) /= "molecule") call error( "invalid packmol command" )
        iarg = 4
        last = narg-1
        molcount = 0
        do while (arg(iarg) == "molecule")
          if (iarg+2 > last) call error( "invalid packmol command" )
          select case (arg(iarg+2))
            case ("move","fix"); nopts = 3
            case ("copy","pack"); nopts = 1
            case default; call error( "invalid packmol command" )
          end select
          call me % packmol_list % add( nopts+2, arg(iarg+1:iarg+2+nopts), repeatable = .true. )
          imol = str2int(arg(iarg+1))
          if ((imol < 1).or.(imol > me%nmol)) call error( "last molecule is", int2str(me%nmol) )
          select case (arg(iarg+2))
            case ("move","fix")
              pos = [(str2real(arg(iarg+2+i)),i=1,3)]
              molcount(imol) = molcount(imol) + 1
            case ("copy","pack")
              n = str2int(arg(iarg+3))
              molcount(imol) = molcount(imol) + n
          end select
          iarg = iarg + 3 + nopts
        end do
        if (iarg /= narg) call error( "invalid packmol command" )
        call me % count_molecules( mass = mass )
        call me % box % compute( sum(mass*molcount) )
        call run_packmol( me % packmol_list, me % molecule_list, me % coordinate_list, &
                          me % nmol, me % atoms_in_molecules(), me % box % length, &
                          seed, tol, action )
        me % nmol = sum(molcount)
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
        catom = trim(me % atom_list % prefix)//arg(1)
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
    write(unit,'("# Generated by playmol")')
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
    if (present(mass).or.present(charge)) then
      ptr => me % molecule_list % first
      do while (associated(ptr))
        atom = ptr % id
        imol = str2int( ptr % params )
        call me % get_types( atom, atom_type )
        if (present(mass)) &
          mass(imol) = mass(imol) + str2real( me % mass_list % parameters( atom_type ) )
        if (present(charge)) &
          charge(imol) = charge(imol) + str2real( me % charge_list % parameters( atom ) )
        ptr => ptr % next
      end do
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
      call split( ptr % params, narg, atom_type(i:i) )
    end do
  end subroutine tPlaymol_get_types

  !=================================================================================================

  subroutine tPlaymol_check_types( me, atom, list )
    class(tPlaymol),    intent(in) :: me
    character(sl),   intent(in) :: atom(:)
    type(StrucList), intent(in) :: list
    character(sl) :: atom_type(list % number)
    call me % get_types( atom, atom_type)
    if (.not.list % find( atom_type )) then
      call error( list%name, join(atom_type), "required, but not found" )
    end if
  end subroutine tPlaymol_check_types

  !=================================================================================================

  subroutine tPlaymol_update_structure( me )
    class(tPlaymol), intent(inout) :: me
    type(Struc), pointer :: b1, b2
    integer :: i, j
    logical :: match(2,2)
    character(sl) :: atoms(4)
    b1 => me % bond_list % first
    do while (associated(b1 % next))
      b2 => b1 % next
      do while (associated(b2 % next))
        do i = 1, 2
          do j = 1, 2
            if (me % bond_list % last % match_id( [b1%id(i), b2%id(j)], two_way = .true. )) then
              atoms = [ b1%id(3-i), b1%id(i), b2%id(j), b2%id(3-j) ]
              call add( me % dihedral_list, me % dihedral_type_list )
            end if
          end do
        end do
        b2 => b2 % next
      end do
      forall (i=1:2, j=1:2) match(i,j) = b1%id(i) == b2%id(j)
      if (any(match)) then
        call find_true_element( i, j )
        atoms = [ b1%id(3-i), b1%id(i), b2%id(3-j), "" ]
        call add( me % angle_list, me % angle_type_list )
      else
        do i = 1, 2
          do j = 1, 2
            if (me % bond_list % find( [b1%id(i), b2%id(j)] )) then
              atoms = [ b1%id(3-i), b1%id(i), b2%id(j), b2%id(3-j) ]
              call add( me % dihedral_list, me % dihedral_type_list )
            end if
          end do
        end do
      end if
      b1 => b1 % next
    end do
    contains
      !---------------------------------------------------------------------------------------------
      subroutine find_true_element( i, j )
        integer, intent(out) :: i, j
        do j = 1, 2
          do i = 1, 2
            if (match(i,j)) return
          end do
        end do
      end subroutine
      !---------------------------------------------------------------------------------------------
      subroutine add( list, type_list )
        type(StrucList), intent(inout) :: list
        type(StrucList), intent(in)    :: type_list
        character(sl) :: types(type_list%number)
        type(Struc), pointer :: ptr
        integer :: n
        n = list%number
        call me % get_types( atoms(1:n), types )
        ptr => type_list % member( types )
        if (associated(ptr)) then
          if (.not.all(match_str( types, ptr%id ))) atoms(1:n) = atoms(n:1:-1)
        end if
        call list % add( list%number, atoms )
        if (.not.associated(ptr)) then
          call warning( "undefined", list%name, "type for atoms", join(atoms) )
        end if
      end subroutine add
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_update_structure

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
    call me % pair_list % print( unit )
    call me % bond_type_list % print( unit )
    call me % angle_type_list % print( unit )
    call me % dihedral_type_list % print( unit )
    call me % atom_list % print( unit )
    call me % charge_list % print( unit )
    call me % bond_list % print( unit )
    call me % angle_list % print( unit, comment = .true. )
    call me % dihedral_list % print( unit, comment = .true. )
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
    integer :: imol, molcount(me % nmol)
    real(rb) :: mass(me % nmol), charge(me % nmol)
    type(Struc), pointer :: ptr
    character(sl) :: atoms
    write(unit,'(A)') repeat("-",80)
    write(unit,'(A)') "SUMMARY"
    write(unit,'(A)') repeat("-",80)
    write(unit,'(A)') "Specified:"
    write(unit,'(I5," atom type(s).")') me % atom_type_list % count()
    write(unit,'(I5," bond type(s).")') me % bond_type_list % count()
    write(unit,'(I5," angle type(s).")') me % angle_type_list % count()
    write(unit,'(I5," dihedral type(s).")') me % dihedral_type_list % count()
    write(unit,'(I5," improper type(s).")') me % improper_type_list % count()
    write(unit,'(I5," atom(s).")') me % atom_list % count()
    write(unit,'(I5," bonds(s).")') me % bond_list % count()
    write(unit,'(I5," improper(s).")') me % improper_list % count()
    write(unit,'(/,A)') "Detected:"
    write(unit,'(I5," angle(s).")') me % angle_list % count()
    write(unit,'(I5," dihedral(s).")') me % dihedral_list % count()
    write(unit,'(I5," molecule(s).")') me % nmol
    call me % count_molecules( molcount, mass, charge )
    do imol = 1, me % nmol
      write(unit,'(/,"Molecule[",A,"]:")') trim(int2str(imol))
      write(unit,'("- Amount: ",A)') trim(int2str(molcount(imol)))
      write(unit,'("- Mass: ",A)') trim(real2str(mass(imol)))
      write(unit,'("- Charge: ",A)') trim(real2str(charge(imol)))
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
    write(unit,'("- Residual charge: ",A)') trim(real2str(sum(molcount*charge)))
    write(unit,'(A)') repeat("-",80)
  end subroutine tPlaymol_summarize

  !=================================================================================================

  subroutine tPlaymol_write_lammps( me, unit )
    class(tPlaymol),  intent(inout) :: me
    integer,       intent(in)    :: unit
    integer :: nb, na, nd, ni, natoms(me % nmol)
    natoms = me % atoms_in_molecules()
    write(unit,'("LAMMPS data file",/,"# Generated by playmol",/)')
    call write_count( me % atom_type_list % count(), "atom types" )
    call write_count( me % bond_type_list % count(), "bond types" )
    call write_count( me % angle_type_list % count(), "angle types" )
    call write_count( me % dihedral_type_list % count(), "dihedral types" )
    call write_count( me % improper_type_list % count(), "improper types" )
    write(unit,'()')
    call write_count( me % coordinate_list % count(), "atoms" )
    call handle_struc( "bonds", me % bond_list, me % bond_type_list, nb )
    call handle_struc( "angles", me % angle_list, me % angle_type_list, na )
    call handle_struc( "dihedrals", me % dihedral_list, me % dihedral_type_list, nd )
    call handle_struc( "impropers", me % improper_list, me % improper_type_list, ni )
    call write_count( nb, "bonds" )
    call write_count( na, "angles" )
    call write_count( nd, "dihedrals" )
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
    contains
      !---------------------------------------------------------------------------------------------
      subroutine write_count( n, name )
        integer,      intent(in) :: n
        character(*), intent(in) :: name
        if (n > 0) write(unit,'(A,X,A)') trim(int2str(n)), name
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
        current => list % first
        if (associated(current)) then
          write(unit,'(/,A,/)') name
          i = 0
          do while (associated(current))
            i = i + 1
            write(unit,'(A,X,A," # ",A)') trim(int2str(i)), trim(current % params), &
                                          trim(join(current % id))
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
            i = i + 1
            mass = me % mass_list % parameters( current % id )
            if (mass == "") call error( "undefined mass for atom type", current%id(1) )
            write(unit,'(A,X,A," # ",A)') trim(int2str(i)), trim(mass), trim(join(current % id))
            current => current % next
          end do
        end if
      end subroutine write_masses
      !---------------------------------------------------------------------------------------------
      subroutine write_atoms
        type(Struc), pointer :: current
        integer :: iatom, itype, i, imol, jmol, narg
        character(sl) :: arg(1), charge
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
            itype = me % atom_type_list % index( arg )
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
        character(sl) :: line, types(list%number)
        character(sl), allocatable :: atom_id(:)
        istruc = 0
        coord => me % coordinate_list % first
        if (.not.present(count)) write(unit,'(/,A,/)') name
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
                if (struc_type % match_id(types,two_way=.true.)) then
                  istruc = istruc + 1
                  line = join(int2str([istruc, itype, last_atom + str_find(struct%id, atom_id)]))
                  if (.not.present(count)) then
                    write(unit,'(A," # ",A)') trim(line), trim(join(struct%id))
                  end if
                end if
                itype = itype + 1
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

end module mPlaymol
