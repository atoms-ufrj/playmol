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

module mMolecule

use mGlobal
use mStruc

implicit none

type tMolecule
  integer :: N = 0
  type(StrucList) :: list = StrucList( "molecule" )
  type(StrucList) :: bonds = StrucList( "bond", 2, .true. )
  type(StrucList) :: xyz = StrucList( "coordinate" )
  contains
    procedure :: add_atom => tMolecule_add_atom
    procedure :: index => tMolecule_index
    procedure :: fuse => tMolecule_fuse
    procedure :: split => tMolecule_split
    procedure :: number_of_atoms => tMolecule_number_of_atoms
    procedure :: count => tMolecule_count
    procedure :: per_molecule => tMolecule_per_molecule
    procedure :: set_geometry => tMolecule_set_geometry
    procedure :: coordinates => tMolecule_coordinates
    procedure :: align => tMolecule_align
end type tMolecule

contains

  !=================================================================================================

  subroutine tMolecule_add_atom( me, atom )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(in)    :: atom
    character(sl) :: arg(2)
    me%N = me%N + 1
    arg = [atom,int2str(me%N)]
    call me % list % add( 2, arg )
  end subroutine tMolecule_add_atom

  !=================================================================================================

  function tMolecule_index( me, string ) result( imol )
    class(tMolecule), intent(in) :: me
    character(sl),    intent(in) :: string
    integer                      :: imol
    integer       :: last
    character(sl) :: atom
    type(Struc), pointer :: ptr
    last = len_trim(string)
    if ((last > 5).and.(string(1:4) == "mol(").and.(string(last:last) == ")")) then
      atom = string(5:last-1)
      if (has_macros(atom)) call error( "invalid atom name", atom )
      call me % list % search( [atom], ptr )
      if (.not.associated(ptr)) call error( "atom", atom, "does not exist" )
      imol = str2int( ptr % params )
    else if (is_int(string)) then
      imol = str2int(string)
      if (imol < 1) call error( "molecule index cannot be zero or negative" )
      if (imol > me%N) call error( "molecule index must lie between 1 and ", int2str(me%N) )
    else
      call error( "molecule definition as", string, "is invalid" )
    end if
  end function tMolecule_index

  !=================================================================================================

  subroutine tMolecule_fuse( me, atoms )
    class(tMolecule), intent(inout)        :: me
    character(sl),    intent(inout)        :: atoms(2)
    integer :: i, mol(2), imin, imax
    type(Struc), pointer :: ptr
    if (atoms(1) == atoms(2)) call error( "atom", atoms(1), "cannot bind to itself" )
    do i = 1, 2
      call me % list % search( atoms(i:i), ptr )
      mol(i) = str2int( ptr % params )
    end do
    imin = minval(mol)
    imax = maxval(mol)
    if (imin == imax) then
      call writeln( "Closing cycle(s) in molecule", int2str(imin) )
    else
      call writeln( "Fusing molecule", int2str(imax), "to molecule", int2str(imin) )
      call rename_molecule( trim(int2str(imax)), trim(int2str(imin)) )
      call writeln( "Lowering indices of molecules", int2str(imax+1), "to", int2str(me%N) )
      do i = imax+1, me%N
        call rename_molecule( trim(int2str(i)), trim(int2str(i-1)) )
      end do
      me%N = me%N - 1
    end if
    call me % bonds % add( 2, atoms, silent = .true. )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine rename_molecule( old, new )
        character(*), intent(in) :: old, new
        ptr => me % list % first
        do while (associated(ptr))
          if (ptr % params == old) ptr % params = new
          ptr => ptr % next
        end do
      end subroutine rename_molecule
      !---------------------------------------------------------------------------------------------
  end subroutine tMolecule_fuse

  !=================================================================================================

  subroutine tMolecule_split( me, arg )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(in)    :: arg(2)
    integer :: i, j, k, imin, imax, natoms
    character(sl) :: cmol, newmol
    type(Struc), pointer :: ptr, pcoords
    integer, allocatable :: mol(:), sort(:)
    character(sl), allocatable :: atom(:), id(:), params(:)
    ! Remove the bond between these atoms (if any) from the list:
    call me % bonds % remove( arg, silent = .true. )
    ! Determine the molecule two which these atoms belong:
    call me % list % search( arg(1:1), ptr )
    cmol = me % list % parameters( arg(1:1) )
    ! Determine the number of atoms of this molecule:
    natoms = 0
    ptr => me % list % first
    do while (associated(ptr))
      if (ptr%params == cmol) natoms = natoms + 1
      ptr => ptr % next
    end do
    ! Save the atoms as if they were monoatomic molecules:
    allocate( mol(natoms), atom(natoms) )
    ptr => me % list % first
    do i = 1, natoms
      mol(i) = i
      do while (ptr%params /= cmol)
        ptr => ptr % next
      end do
      atom(i) = ptr%id(1)
      ptr => ptr % next
    end do
    ! Search for non-deleted bonds and reunite the atoms as molecules:
    do i = 1, natoms-1
      do j = i+1, natoms
        if (me % bonds % find([atom(i),atom(j)])) then
          if (mol(i) /= mol(j)) then
            imin = min(mol(i),mol(j))
            imax = max(mol(i),mol(j))
            forall (k=1:natoms,mol(k) == imax) mol(k) = imin
            forall (k=1:natoms,mol(k) > imax) mol(k) = mol(k) - 1
          end if
        end if
      end do
    end do
    ! If the removed bond used to define two separate regions, there will now be two molecules:
    if (any(mol == 2)) then
      me%N = me%N + 1
      newmol = int2str(me%N)
      call writeln( "Splitting molecule", cmol, "into molecules", cmol, "and", newmol )
      ptr => me % list % first
      do while (associated(ptr))
        if (any((atom == ptr%id(1)).and.(mol == 2))) ptr%params = newmol
        ptr => ptr % next
      end do
      ! Sort the arrays in accordance with the order in mol:
      call writeln( "Reorganizing coordinates...", advance = .false. )
      allocate( sort(natoms) )
      sort = sorted( mol )
      atom = atom(sort)
      allocate( id(natoms), params(natoms) )
      ptr => me % xyz % first
      do while (associated(ptr))
        if (any(atom == ptr%id(1))) then
          pcoords => ptr
          ! Get the current order of atoms:
          do i = 1, natoms
            id(i) = ptr % id(1)
            params(i) = ptr % params
            ptr => ptr % next
          end do
          ! Reset the atoms in the new order:
          do i = 1, natoms
            pcoords % id(1) = atom(i)
            j = 1
            do while (id(j) /= atom(i))
              j = j + 1
            end do
            pcoords % params = params(j)
            pcoords => pcoords % next
          end do
        else
          ptr => ptr % next
        end if
      end do
      call writeln( " done." )
    else
      call writeln( "Opening cycle(s) in molecule", cmol )
    end if
  end subroutine tMolecule_split

  !=================================================================================================

  function tMolecule_number_of_atoms( me ) result( natoms )
    class(tMolecule), intent(in) :: me
    integer                      :: natoms(me%N)
    type(Struc), pointer :: atom
    integer :: i
    natoms = 0
    atom => me % list % first
    do while (associated(atom))
      i = str2int(atom%params)
      natoms(i) = natoms(i) + 1
      atom => atom % next
    end do
  end function tMolecule_number_of_atoms

  !=================================================================================================

  function tMolecule_count( me ) result( nmols )
    class(tMolecule), intent(in) :: me
    integer                      :: nmols(me%N)
    integer  :: imol, iatom, natoms(me%N)
    type(Struc), pointer :: ptr
    natoms = me % number_of_atoms()
    nmols = 0
    ptr => me % xyz % first
    do while (associated(ptr))
      imol = str2int(me % list % parameters( ptr%id ) )
      nmols(imol) = nmols(imol) + 1
      do iatom = 1, natoms(imol)
        ptr => ptr % next
      end do
    end do
  end function tMolecule_count

  !=================================================================================================

  function tMolecule_per_molecule( me, val_list ) result( total )
    class(tMolecule), intent(in) :: me
    type(StrucList),  intent(in) :: val_list
    real(rb)                     :: total(me%N)
    integer :: imol
    type(Struc), pointer :: ptr, pval
    total = 0.0_rb
    pval => val_list % first
    do while (associated(pval))
      ptr => me % list % first
      do while (associated(ptr))
        if (ptr % match_id(pval % id,.false.)) then
          imol = str2int( ptr % params )
          total(imol) = total(imol) + str2real( pval % params )
        end if
        ptr => ptr % next
      end do
      pval => pval % next
    end do
  end function tMolecule_per_molecule

  !=================================================================================================

  subroutine tMolecule_set_geometry( me, data, ndata )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(in)    :: data(:,:)
    integer,          intent(in)    :: ndata(size(data,1))
    integer       :: N, i, j, narg, imol, iatom, ind(3)
    character(sl) :: arg(size(data,2)), catom
    integer :: natoms(me%N)
    logical :: new_molecule
    type(Struc), pointer :: atom
    character(sl), allocatable :: prev(:)
    character(sl), allocatable :: name(:)
    real(rb), allocatable :: R(:,:)
    real(rb) :: L, theta, phi, R1(3), R2(3), R3(3), x(3), y(3), z(3)
    natoms = me % number_of_atoms()
    allocate( prev(maxval(natoms)) )
    N = size(ndata)
    call writeln( "Number of provided geometric data: ", int2str(N) )
    allocate( name(N), R(3,N) )
    new_molecule = .true.
    do i = 1, N
      narg = ndata(i)
      arg = data(i,:)
      if ((narg < 3).or.(narg > 7).or.(narg == 6)) call error( "invalid geometric info format" )
      catom = arg(1)
      call me % list % search( [catom], atom )
      if (.not.associated(atom)) call error( "invalid atom", catom )
      if (new_molecule) then
        imol = str2int(atom%params)
        iatom = 1
        call writeln( "Processing", int2str(natoms(imol)), &
                      "geometric data for molecule", trim(int2str(imol))//":" )
      else if (trim(atom%params) /= trim(int2str(imol))) then
        call error( "atom", catom, "does not belong to molecule", int2str(imol) )
      else if (any(str_find([catom],prev(1:iatom)) > 0)) then
        call error( "repeated geometric info for atom", catom )
      else
        iatom = iatom + 1
      end if
      call writeln( "Data provided for atom", catom, ":", join(arg(2:narg)) )
      name(i) = catom
      select case (narg)
        case (3) ! Bond
          call check_atoms( arg(2:2), ind(1:1) )
          L = str2real(arg(3))
          x = real([1,0,0],rb)
          R(:,i) = R(:,ind(1)) + L*x
        case (4) ! Coordinates
          do j = 1, 3
            R(j,i) = str2real(arg(j+1))
          end do
        case (5) ! Bond and angle
          call check_atoms( arg(2:3), ind(1:2) )
          L = str2real(arg(4))
          theta = str2real(arg(5))
          R1 = R(:,ind(1))
          R2 = R(:,ind(2))
          x = (R1 - R2)/norm(R1 - R2)
          y = real([0,1,0],rb)
          if (abs(x(2)-1.0_rb) < 0.01_rb) y = real([1,0,0],rb)
          y = y - scalar(y,x)*x
          y = y / norm(y)
          R(:,i) = R1 + L*(cosine(180-theta)*x + sine(180-theta)*y)
        case (7) ! Bond, angle, and dihedral
          call check_atoms( arg(2:4), ind(1:3) )
          L = str2real(arg(5))
          theta = str2real(arg(6))
          phi = str2real(arg(7))
          R1 = R(:,ind(1))
          R2 = R(:,ind(2))
          R3 = R(:,ind(3))
          x = (R1 - R2)/norm(R1 - R2)
          y = R3 - R2 - scalar(R3 - R2,x)*x
          y = y / norm(y)
          z = cross(x,y)
          R(:,i) = R1 + L*(cosine(180-theta)*x + sine(180-theta)*(cosine(phi)*y + sine(phi)*z))
        case default
          call error( "invalid geometric info" )
      end select
      arg(2:4) = real2str(R(:,i))
      call me % xyz % add( 4, arg(1:4), me % list, repeatable = .true. )
      prev(iatom) = catom
      new_molecule = iatom == natoms(imol)
    end do
    if (.not.new_molecule) then
      call error( "geometric info for molecule", int2str(imol), "is incomplete" )
    end if
    contains
      subroutine check_atoms( atom, ind )
        character(sl), intent(in) :: atom(:)
        integer,       intent(inout) :: ind(size(atom))
        integer :: j, k
        do j = 1, size(atom)
          if (.not. me % list % find( [atom(j)])) call error( "invalid atom", atom(j) )
          if (any(name(1:i-1) == atom(j))) then
            k = 1
            do while (name(k) /= atom(j))
              k = k + 1
            end do
            ind(j) = k
          else
            call error( "no coordinates defined for atom", atom(j) )
          end if
        end do
      end subroutine
  end subroutine tMolecule_set_geometry

  !=================================================================================================

  subroutine tMolecule_coordinates( me, imol, N, Coord, atom, option )
    class(tMolecule), intent(inout) :: me
    integer,          intent(in)    :: imol, N, option
    real(rb),         intent(inout) :: Coord(3,N)
    character(sl),    intent(inout) :: atom(N)
    type(Struc), pointer :: current
    logical :: found
    integer :: i, narg
    character(sl) :: arg(3)
    current => me % xyz % first
    found = .false.
    do while (associated(current).and.(.not.found))
      found = str2int(me % list % parameters( current%id )) == imol
      if (.not.found) current => current % next
    end do
    if (.not.found) call error( "no coordinates for molecule", int2str(imol) )
    do i = 1, N
      if (option == 1) then ! Retrieve coordinates:
        atom(i) = current % id(1)
        call split( current % params, narg, arg )
        Coord(:,i) = [str2real(arg(1)), str2real(arg(2)), str2real(arg(3))]
      else ! Set coordinates:
        current % id(1) = atom(i)
        current % params = join(real2str(Coord(:,i)))
      end if
      current => current % next
    end do
  end subroutine tMolecule_coordinates

  !=================================================================================================

  subroutine tMolecule_align( me, N, Mass, Coord, axis )
    class(tMolecule), intent(inout) :: me
    integer,  intent(in)    :: N, axis(3)
    real(rb), intent(in)    :: Mass(N)
    real(rb), intent(inout) :: Coord(3,N)
    integer :: i
    real(rb) :: MolMass, Rcm(3), inertia(3,3), delta(3,N), MoI(3), A(3,3)
    if (N == 1) then
      Coord = 0.0_rb
    else
      ! Compute molecular mass and the center-of-mass position:
      MolMass = sum(Mass)
      forall (i=1:3) Rcm(i) = sum(Mass*Coord(i,:))/MolMass
      ! Compute inertia tensor:
      inertia = 0.0_rb
      do i = 1, N
        delta(:,i) = Coord(:,i) - Rcm
        call add_inertia( inertia, Mass(i), delta(:,i) )
      end do
      ! Diagonalize the inertia tensor and compute rotation matrix:
      MoI(axis) = eigenvalues( inertia )
      A = eigenvectors( inertia, MoI )
      A = transpose(A(:,sort_vector(MoI)))
      ! Recalculate positions in the body-fixed frame:
      forall (i=1:N) Coord(:,i) = matmul( A, delta(:,i) )
    end if
    contains
    !-----------------------------------------------------------------------------------------------
    subroutine add_inertia( inertia, mass, delta )
      real(rb), intent(inout) :: inertia(3,3)
      real(rb), intent(in)    :: mass, delta(3)
      inertia(1,1) = inertia(1,1) + mass*(delta(2)**2 + delta(3)**2)
      inertia(2,2) = inertia(2,2) + mass*(delta(1)**2 + delta(3)**2)
      inertia(3,3) = inertia(3,3) + mass*(delta(1)**2 + delta(2)**2)
      inertia(1,2) = inertia(1,2) - mass*delta(1)*delta(2)
      inertia(1,3) = inertia(1,3) - mass*delta(1)*delta(3)
      inertia(2,3) = inertia(2,3) - mass*delta(2)*delta(3)
      inertia(2,1) = inertia(1,2)
      inertia(3,1) = inertia(1,3)
      inertia(3,2) = inertia(2,3)
    end subroutine add_inertia
    !-----------------------------------------------------------------------------------------------
    function eigenvalues( Matrix ) result( Lambda )
      real(rb), intent(in) :: Matrix(3,3)
      real(rb)             :: Lambda(3)
      real(rb), parameter :: SQRT3 = 1.7320508075688773_rb
      real(rb) :: A(3,3), M, C1, C0
      real(rb) :: DE, DD, EE, FF
      real(rb) :: P, SQRTP, Q, C, S, PHI
      A = Matrix
      DE = A(1,2)*A(2,3)
      DD = A(1,2)**2
      EE = A(2,3)**2
      FF = A(1,3)**2
      M  = A(1,1) + A(2,2) + A(3,3)
      C1 = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
      C0 = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0_rb * A(1,3)*DE
      P = M**2 - 3.0_rb * C1
      Q = M*(P - 1.5_rb*C1) - 13.5_rb*C0
      SQRTP = sqrt(abs(P))
      PHI = 27.0_rb*(0.25_rb * C1**2 * (P - C1) + C0 * (Q + 6.75_rb*C0))
      PHI = atan2(sqrt(abs(PHI)),Q)/3.0_rb
      C = SQRTP*cos(PHI)
      S = SQRTP*sin(PHI)/SQRT3
      Lambda(2) = (M - C)/3.0_rb
      Lambda(3) = Lambda(2) + S
      Lambda(1) = Lambda(2) + C
      Lambda(2) = Lambda(2) - S
    end function eigenvalues
    !-----------------------------------------------------------------------------------------------
    function eigenvectors( Matrix, W ) result( Q )
      real(rb), intent(in) :: Matrix(3,3), W(3)
      real(rb)             :: Q(3,3)
      real(rb), parameter :: EPS = 2.2204460492503131e-16_rb
      real(rb) :: A(3,3), NORM, N1, N2, N1TMP, N2TMP
      real(rb) :: THRESH, ERROR, WMAX, F, T
      integer :: I, J
      logical :: SUCCESS
      A = Matrix
      WMAX   = max(abs(W(1)), abs(W(2)), abs(W(3)))
      THRESH = (8.0_rb * EPS * WMAX)**2
      N1TMP = A(1,2)**2 + A(1,3)**2
      N2TMP = A(1,2)**2 + A(2,3)**2
      Q(1,1) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
      Q(1,2) = Q(1,1)
      Q(2,1) = A(1,3) * A(1,2) - A(2,3) * A(1,1)
      Q(2,2) = Q(2,1)
      Q(3,2) = A(1,2)**2
      A(1,1) = A(1,1) - W(1)
      A(2,2) = A(2,2) - W(1)
      Q(1,1) = Q(1,2) + A(1,3) * W(1)
      Q(2,1) = Q(2,2) + A(2,3) * W(1)
      Q(3,1) = A(1,1) * A(2,2) - Q(3,2)
      NORM = Q(1,1)**2 + Q(2,1)**2 + Q(3,1)**2
      N1 = N1TMP + A(1,1)**2
      N2 = N2TMP + A(2,2)**2
      ERROR = N1 * N2
      if (N1 <= THRESH) then
        Q(:,1) = [1.0_rb,0.0_rb,0.0_rb]
      else if (N2 <= THRESH) then
        Q(:,1) = [0.0_rb,1.0_rb,0.0_rb]
      else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
        T = abs(A(1,2))
        F = -A(1,1) / A(1,2)
        if (abs(A(2,2)) > T) then
          T = abs(A(2,2))
          F = -A(1,2) / A(2,2)
        end if
        if (abs(A(2,3)) < T) then
          F = -A(1,3) / A(2,3)
        end if
        NORM = 1.0_rb / sqrt(1.0_rb + F**2)
        Q(:,1) = [NORM,F*NORM,0.0_rb]
      else
        NORM = sqrt(1.0_rb / NORM)
        Q(:,1) = Q(:,1) * NORM
      end if
      T = W(1) - W(2)
      if (abs(T) < 8.0_rb * EPS * WMAX) then
        A(1,1) = A(1,1) + T
        A(2,2) = A(2,2) + T
        Q(1,2) = Q(1,2) + A(1,3) * W(2)
        Q(2,2) = Q(2,2) + A(2,3) * W(2)
        Q(3,2) = A(1,1) * A(2,2) - Q(3,2)
        NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
        N1 = N1TMP + A(1,1)**2
        N2 = N2TMP + A(2,2)**2
        ERROR = N1 * N2
        if (N1 <= THRESH) then
          Q(:,2) = [1.0_rb,0.0_rb,0.0_rb]
        else if (N2 <= THRESH) then
          Q(:,2) = [0.0_rb,1.0_rb,0.0_rb]
        else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
          T = abs(A(1,2))
          F = -A(1,1) / A(1,2)
          if (abs(A(2,2)) > T) then
            T = abs(A(2,2))
            F = -A(1,2) / A(2,2)
          end if
          if (abs(A(2,3)) > T) then
            F = -A(1,3) / A(2,3)
          end if
          NORM = 1.0_rb / sqrt(1.0_rb + F**2)
          Q(:,2) = [NORM,F*NORM,0.0_rb]
        else
          NORM = sqrt(1.0_rb / NORM)
          Q(:,2) = Q(:,2) * NORM
        end if
      else
        A(2,1) = A(1,2)
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)
        A(1,1) = A(1,1) + W(1)
        A(2,2) = A(2,2) + W(1)
        do I = 1, 3
          A(I,I) = A(I,I) - W(2)
          N1 = A(1,I)**2 + A(2,I)**2 + A(3,I)**2
          if (N1 > THRESH) then
            Q(1,2) = Q(2,1) * A(3,I) - Q(3,1) * A(2,I)
            Q(2,2) = Q(3,1) * A(1,I) - Q(1,1) * A(3,I)
            Q(3,2) = Q(1,1) * A(2,I) - Q(2,1) * A(1,I)
            NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
            SUCCESS = NORM <= (256.0_rb * EPS)**2 * N1
            if (.not.SUCCESS) then
              NORM = sqrt(1.0_rb / NORM)
              Q(:, 2) = Q(:, 2) * NORM
              exit
            end if
          end if
        end do
        if (SUCCESS) then
          do J = 1, 3
            if (Q(J,1) <= 0.0_rb) then
              I = 1 + mod(J,3)
              NORM = 1.0_rb / sqrt(Q(J,1)**2 + Q(I,1)**2)
              Q(J,2) = Q(I,1) * NORM
              Q(I,2) = -Q(J,1) * NORM
              Q(I,2) = 0.0_rb
              exit
            end if
          end do
        end if
      end if
      Q(:,3) = cross_product(Q(:,1),Q(:,2))
    end function eigenvectors
    !-----------------------------------------------------------------------------------------------
    function determinant( A ) result( det )
      real(rb), intent(in) :: A(3,3)
      real(rb)             :: det
      det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + &
            A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
    end function determinant
    !-----------------------------------------------------------------------------------------------
    function sort_vector( a ) result( b )
      real(rb), intent(in) :: a(3)
      integer              :: b(3)
      integer :: imin, imax
      imin = minloc(a,dim=1)
      imax = maxloc(a,dim=1)
      if (imin == imax) then
        b = [1,2,3]
      else
        b = [imin,6-imin-imax,imax]
      end if
    end function sort_vector
    !-----------------------------------------------------------------------------------------------
    pure function cross_product(a, b) result( c )
      real(rb), intent(in) :: a(3), b(3)
      real(rb)             :: c(3)
      c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
    end function cross_product
    !-----------------------------------------------------------------------------------------------
    function inv3x3( A ) result( AINV )
      real(rb), intent(in) :: A(3,3)
      real(rb)             :: AINV(3,3)
      real(rb), parameter :: EPS = 1.0E-10_rb
      real(rb) :: DET, COFACTOR(3,3)
      DET = determinant( A )
      if (abs(DET) <= EPS) call error( "could not invert 3x3 matrix" )
      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
      AINV = transpose(COFACTOR) / DET
    end function inv3x3
    !-----------------------------------------------------------------------------------------------
  end subroutine  tMolecule_align
  !=================================================================================================
end module mMolecule
