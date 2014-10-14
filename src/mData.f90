module mData

use mGlobal
use mStruc

implicit none

type tData
  integer  :: nmol = 0
  real(rb) :: density = 0.0_rb, aspect(3) = 1.0_rb
  type(StrucList) :: atom_type_list = StrucList( name = "atom_type", number = 1 )
  type(StrucList) :: bond_type_list = StrucList( name = "bond_type", number = 2 )
  type(StrucList) :: angle_type_list = StrucList( name = "angle_type", number = 3 )
  type(StrucList) :: dihedral_type_list = StrucList( name = "dihedral_type", number = 4, two_way = .false. )
  type(StrucList) :: improper_type_list = StrucList( name = "improper_type", number = 4, two_way = .false. )
  type(StrucList) :: mass_list = StrucList( name = "mass", number = 1 )
  type(StrucList) :: pair_list = StrucList( name = "pair", number = 2 )
  type(StrucList) :: atom_list = StrucList( name = "atom", number = 1 )
  type(StrucList) :: charge_list = StrucList( name = "charge", number = 1 )
  type(StrucList) :: bond_list = StrucList( name = "bond", number = 2 )
  type(StrucList) :: angle_list = StrucList( name = "angle", number = 3 )
  type(StrucList) :: dihedral_list = StrucList( name = "dihedral", number = 4 )
  type(StrucList) :: improper_list = StrucList( name = "improper", number = 4 )
  type(StrucList) :: molecule_list = StrucList( name = "molecule", number = 1 )
  type(StrucList) :: coordinate_list = StrucList( name = "coordinate", number = 1 )
  contains
    procedure :: read => tData_Read
    procedure :: write => tData_write
    procedure :: write_lammps => tData_write_lammps
    procedure :: read_xyz => tData_read_xyz
    procedure :: summarize => tData_summarize
    procedure :: fuse_molecules => tData_fuse_molecules
    procedure :: count_molecules => tData_count_molecules
    procedure :: update_structure => tData_update_structure
    procedure :: get_types => tData_get_types
    procedure :: check_types => tData_check_types
    procedure :: atoms_in_molecules => tData_atoms_in_molecules
end type tData

contains

  !=================================================================================================

  recursive subroutine tData_Read( me, unit )
    class(tData), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: narg
    character(sl) :: arg(10)
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
        case ("pair"); call pair_command
        case ("atom"); call atom_command
        case ("charge"); call charge_command
        case ("bond"); call bond_command
        case ("improper"); call improper_command
        case ("xyz"); call xyz_command
        case ("write"); call write_command
        case ("include"); call include_command
        case ("print_summary"); call me % summarize
        case ("quit"); stop
        case default; call error( "unrecognizable command", arg(1) )
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
        select case (narg)
          case (2)
            me % density = str2real(arg(2))
          case (5)
            me % density = str2real(arg(2))
            me % aspect = [(str2real(arg(2+i)),i=1,3)]
          case default
            call error( "invalid box command" )
        end select
        call writeln( "Defining box density =", arg(2) )
        if (narg == 5) call writeln( "Defining box aspect = ", join(arg(3:5)) )
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
        call me % mass_list % add( narg-1, arg(2:narg), me % atom_type_list )
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
          call me % read_xyz( unit )
        else
          inquire( file = arg(2), exist = file_exists )
          if (.not.file_exists) call error( "file", arg(2), "does not exist" )
          open( newunit = xyz, file = arg(2), status = "old" )
          call me % read_xyz( xyz )
          close( xyz )
        end if
      end subroutine xyz_command
      !---------------------------------------------------------------------------------------------
      subroutine write_command
        integer :: unit
        if (narg < 3) call error( "invalid write command")
        if ((arg(2) /= "playmol").and.(arg(2) /= "lammps")) call error( "invalid write command" )
        call writeln( "Writing data to file", arg(3), "in", arg(2), "format..." )
        open( newunit = unit, file = arg(3), status = "replace" )
        select case (arg(2))
          case ("playmol"); call me % write( unit )
          case ("lammps");  call me % write_lammps( unit )
        end select
        close(unit)
      end subroutine write_command
      !---------------------------------------------------------------------------------------------
      subroutine include_command
        logical :: file_exists
        integer :: input
        if (narg == 1) call error( "invalid include command" )
        inquire( file = arg(2), exist = file_exists )
        if (.not.file_exists) call error( "file", arg(2), "does not exist" )
        open( newunit = input, file = arg(2), status = "old" )
        call me % Read( input )
        close(input)
      end subroutine include_command
      !---------------------------------------------------------------------------------------------
  end subroutine tData_Read

  !=================================================================================================

  subroutine tData_read_xyz( me, unit )
    class(tData), intent(inout) :: me
    integer,      intent(in)    :: unit
    integer       :: N, i, narg, imol, iatom
    character(sl) :: arg(4), line, catom
    integer :: natoms(me % nmol)
    logical :: new_molecule
    type(Struc), pointer :: atom
    character(sl), allocatable :: prev(:)
    natoms = me % atoms_in_molecules()
    allocate( prev(maxval(natoms)) )
    call writeln("Reading xyz data...")
    call next_command( unit, narg, arg )
    if (narg > 0) then
      call writeln( "Number of coordinates: ", arg(1) )
      N = str2int( arg(1) )
      read(unit,'(A'//csl//')') line
      new_molecule = .true.
      do i = 1, N
        read(unit,'(A'//csl//')') line
        call split( line, narg, arg )
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
        prev(iatom) = arg(1)
        new_molecule = iatom == natoms(imol)
      end do
      if (.not.new_molecule) then
        call error( "coordinates of molecule", int2str(imol), "are incomplete" )
      end if
    end if
  end subroutine tData_read_xyz

  !=================================================================================================

  subroutine tData_summarize( me )
    class(tData), intent(inout) :: me
    integer :: imol, molcount(me % nmol)
    real(rb) :: mass(me % nmol), charge(me % nmol)
    call writeln( repeat("-",80) )
    call writeln( "SUMMARY" )
    call writeln( repeat("-",80) )
    call writeln( "Specified:")
    call writeln( "*", int2str(me % atom_type_list % count()), "atom type(s)")
    call writeln( "*", int2str(me % bond_type_list % count()), "bond type(s)" )
    call writeln( "*", int2str(me % angle_type_list % count()), "angle type(s)" )
    call writeln( "*", int2str(me % dihedral_type_list % count()), "dihedral type(s)" )
    call writeln( "*", int2str(me % improper_type_list % count()), "improper type(s)" )
    call writeln( "*", int2str(me % atom_list % count()), "atom(s)" )
    call writeln( "*", int2str(me % bond_list % count()),"bond(s)" )
    call writeln( "*", int2str(me % improper_list % count()), "improper(s)" )
    call writeln( "Detected:")
    call writeln( "*", int2str(me % angle_list % count()), "angle(s)" )
    call writeln( "*", int2str(me % dihedral_list % count()), "dihedral(s)" )
    call writeln( "*", int2str(me % nmol), "molecule(s)" )
    call me % count_molecules( molcount, mass, charge )
    if (me % nmol > 1) then
      call writeln
      do imol = 1, me % nmol
        call writeln( "Molecule[", int2str(imol),"]:" )
        call writeln( "- Amount:", int2str(molcount(imol)) )
        call writeln( "- Mass:", real2str(mass(imol)) )
        call writeln( "- Charge:", real2str(charge(imol)) )
      end do
    end if
    call writeln
    call writeln( "All molecules:" )
    call writeln( "- Amount:", int2str(sum(molcount)) )
    call writeln( "- Mass:", real2str(sum(molcount*mass)) )
    call writeln( "- Charge:", real2str(sum(molcount*charge)) )
    call writeln( repeat("-",80) )
  end subroutine tData_summarize

  !=================================================================================================

  subroutine tData_fuse_molecules( me, atom )
    class(tData),    intent(inout) :: me
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
  end subroutine tData_fuse_molecules

  !=================================================================================================

  subroutine tData_count_molecules( me, number, mass, charge )
    class(tData),  intent(in)            :: me
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
  end subroutine tData_count_molecules

  !=================================================================================================

  subroutine tData_get_types( me, atom, atom_type )
    class(tData),  intent(in)  :: me
    character(sl), intent(in)  :: atom(:)
    character(sl), intent(out) :: atom_type(:)
    integer :: i, narg
    type(Struc), pointer :: ptr
    do i = 1, size(atom)
      call me % atom_list % search( atom(i:i), ptr )
      call split( ptr % params, narg, atom_type(i:i) )
    end do
  end subroutine tData_get_types

  !=================================================================================================

  subroutine tData_check_types( me, atom, list )
    class(tData),    intent(in) :: me
    character(sl),   intent(in) :: atom(:)
    type(StrucList), intent(in) :: list
    character(sl) :: atom_type(list % number)
    call me % get_types( atom, atom_type)
    if (.not.list % find( atom_type )) then
      call error( list%name, join(atom_type), "required, but not found" )
    end if
  end subroutine tData_check_types

  !=================================================================================================

  subroutine tData_update_structure( me )
    class(tData), intent(inout) :: me
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
  end subroutine tData_update_structure

  !=================================================================================================

  function tData_atoms_in_molecules( me ) result( natoms )
    class(tData), intent(in) :: me
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
  end function tData_atoms_in_molecules

  !=================================================================================================

  subroutine tData_write( me, unit )
    class(tData), intent(in) :: me
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
  end subroutine tData_write

  !=================================================================================================

  subroutine tData_write_lammps( me, unit )
    class(tData),  intent(in) :: me
    integer,       intent(in) :: unit
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
    if (me % density > 0.0_rb) call write_box_limits
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
        integer  :: molcount(me % nmol)
        real(rb) :: mass(me % nmol), V, L
        call me % count_molecules( molcount, mass )
        V = sum(molcount*mass) / me % density
        L = (V/product(me % aspect))**(1.0_rb/3.0_rb)
        write(unit,'(/,A," xlo xhi")') "0 "//trim(real2str(L * me % aspect(1)))
        write(unit,'(  A," ylo yhi")') "0 "//trim(real2str(L * me % aspect(2)))
        write(unit,'(  A," zlo zhi")') "0 "//trim(real2str(L * me % aspect(3)))
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
  end subroutine tData_write_lammps

  !=================================================================================================

end module mData
