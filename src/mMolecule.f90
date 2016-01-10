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
  type(StrucList) :: list = StrucList( name = "molecule", number = 1 )
  type(StrucList) :: xyz = StrucList( name = "coordinate", number = 1 )
  contains
    procedure :: index => tMolecule_index
    procedure :: coordinates => tMolecule_coordinates
    procedure :: align => tMolecule_align
end type tMolecule

contains

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

  subroutine tMolecule_coordinates( me, imol, N, Coord, lcoord, option )
    class(tMolecule), intent(inout) :: me
    integer,          intent(in)    :: imol, N, option
    real(rb),         intent(inout) :: Coord(3,N)
    type(StrucList),  intent(inout) :: lcoord
    type(Struc), pointer :: current
    logical :: found
    integer :: i, narg
    character(sl) :: arg(3)
    current => lcoord % first
    found = .false.
    do while (associated(current).and.(.not.found))
      found = str2int(me % list % parameters( current%id )) == imol
      if (.not.found) current => current % next
    end do
    do i = 1, N
      if (option == 1) then ! Retrieve coordinates:
        call split( current % params, narg, arg )
        Coord(:,i) = [str2real(arg(1)), str2real(arg(2)), str2real(arg(3))]
      else ! Set coordinates:
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
