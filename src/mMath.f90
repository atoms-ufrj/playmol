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

module mMath

use mGlobal

implicit none

contains

  !=================================================================================================
  ! Numerical diagonalization of 3x3 matrcies
  ! Copyright (C) 2006  Joachim Kopp
  !=================================================================================================

  pure subroutine diagonalization( matrix, q, w )
    real(rb), intent(in)  :: matrix(3,3)
    real(rb), intent(out) :: q(3,3), w(3)

    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
    ! matrix A using Cardano's method for the eigenvalues and an analytical
    ! method based on vector cross products for the eigenvectors.
    ! Only the diagonal and upper triangular parts of A need to contain
    ! meaningful values. However, all of A may be used as temporary storage
    ! and may hence be destroyed.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   matrix: The symmetric input matrix
    !   q: Storage buffer for eigenvectors
    !   w: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------

    real(rb), parameter :: zero = 0.0_rb, one = 1.0_rb, two = 2.0_rb, three = 3.0_rb
    real(rb), parameter :: half = one/two, third = one/three
    real(rb), parameter :: eps = epsilon(one)
    real(rb), parameter :: sqrt3 = sqrt(3.0_rb)

    integer  :: i, j
    real(rb) :: a(3,3), norm, n1, n2, w1, w2, w3, thresh, t, wmax8eps
    real(rb) :: m, c1, c0, de, dd, ee, ff, p, sqrtp, r, c, s, phi
    logical  :: success

    a = matrix

    ! Calculate the eigenvalues of a symmetric 3x3 matrix a using Cardano's
    ! analytical algorithm. Only the diagonal and upper triangular parts of A are
    ! accessed. The access is read-only.
    de = a(1,2)*a(2,3)
    dd = a(1,2)**2
    ee = a(2,3)**2
    ff = a(1,3)**2
    m  = a(1,1) + a(2,2) + a(3,3)
    c1 = (a(1,1)*a(2,2) + a(1,1)*a(3,3) + a(2,2)*a(3,3)) - (dd + ee + ff)
    c0 = 27.0_rb*(a(3,3)*dd + a(1,1)*ee + a(2,2)*ff - a(1,1)*a(2,2)*a(3,3) - two*a(1,3)*de)

    p = m*m - 3.0_rb*c1
    r = m*(p - 1.5_rb*c1) - half*c0
    sqrtp = sqrt(abs(p))
    phi = third*atan2(sqrt(abs(6.75_rb*c1*c1*(p - c1) + c0*(r + 0.25_rb*c0))),r)

    c = sqrtp*cos(phi)
    s = (one/sqrt3)*sqrtp*sin(phi)

    p = third*(m - c)
    w1 = p + c
    w2 = p - s
    w3 = p + s

    ! Sort eigenvalues:
    if (abs(w1) > abs(w3)) call swap( w1, w3 )
    if (abs(w1) > abs(w2)) call swap( w1, w2 )
    if (abs(w2) > abs(w3)) call swap( w2, w3 )
    w = [w1, w2, w3]

    wmax8eps = 8.0_rb*eps*abs(w1)
    thresh = wmax8eps**2

    ! Prepare calculation of eigenvectors
    n1 = a(1,2)**2 + a(1,3)**2
    n2 = a(1,2)**2 + a(2,3)**2
    q(1,1) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    q(1,2) = q(1,1)
    q(2,1) = a(1,3)*a(1,2) - a(2,3)*a(1,1)
    q(2,2) = q(2,1)
    q(3,2) = a(1,2)**2

    ! Calculate first eigenvector by the formula v(1) = (A - lambda(1)).e1 x (A - lambda(1)).e2
    a(1,1) = a(1,1) - w1
    a(2,2) = a(2,2) - w1
    q(:,1) = [q(1,2) + a(1,3)*w1, q(2,2) + a(2,3)*w1, a(1,1)*a(2,2) - q(3,2)]
    call compute_eigenvector( q(:,1), a, n1, n2 )

    ! Prepare calculation of second eigenvector     
    t = w1 - w2

    ! Is this eigenvalue degenerate?
    if (abs(t) > wmax8eps) then

      ! For non-degenerate eigenvalue, calculate second eigenvector by the formula
      !         v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
      a(1,1) = a(1,1) + t
      a(2,2) = a(2,2) + t
      q(:,2) = [q(1,2) + a(1,3)*w2, q(2,2) + a(2,3)*w2, a(1,1)*a(2,2) - q(3,2)]
      call compute_eigenvector( q(:,2), a, n1, n2 )

    else

      ! For degenerate eigenvalue, calculate second eigenvector according to
      !         v[1] = v(1) x (A - lambda[1]).e[i]

      ! This would really get too complicated if we could not assume all of A to
      !       contain meaningful values.
      a(2,1) = a(1,2)
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      a(1,1) = a(1,1) + w1
      a(2,2) = a(2,2) + w1
      i = 0
      success = .false.
      do while ((i < 3).and.(.not.success))
        i = i + 1
        a(i,i) = a(i,i) - w2
        n1 = sum(a(:,i)**2)
        success = n1 > thresh
        if (success) then
          q(:,2) = cross_product( q(:,1), a(:,i) )
          norm = sum(q(:,2)**2)
          success = norm > (256.0_rb*eps)**2*n1
          if (success) q(:,2) = q(:,2)*sqrt(one/norm)
        end if
      end do

      ! This means that any vector orthogonal to v(1) is an EV.
      if (.not.success) then
        i = 1
        do while (q(i,1) == zero)
          i = i + 1
        end do
        j = 1 + mod(i,3)
        norm = one/sqrt(q(i,1)**2 + q(j,1)**2)
        q(i,2) =  q(j,1)*norm
        q(j,2) = -q(i,1)*norm
        q(1+mod(i+1,3),2) = zero
      end if
    end if

    ! Calculate third eigenvector according to v[2] = v(1) x v[1]
    q(:,3) = cross_product( q(:,1), q(:,2) )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine compute_eigenvector( q, a, n1tmp, n2tmp )
        real(rb), intent(inout) :: q(3)
        real(rb), intent(in)    :: a(3,3), n1tmp, n2tmp

        real(rb) :: norm, n1, n2, error, t, f

        norm = sum(q**2)
        n1 = n1tmp + a(1,1)**2
        n2 = n2tmp + a(2,2)**2
        error = n1*n2

        ! If the first column is zero, then (1, 0, 0) is an eigenvector
        if (n1 <= thresh) then
          q = [one, zero, zero]

        ! If the second column is zero, then (0, 1, 0) is an eigenvector
        else if (n2 <= thresh) then
          q = [zero, one, zero]

        ! If angle between A(*,1) and A(*,2) is too small, don't use
        !  cross product, but calculate v ~ (1, -A0/A1, 0)
        else if (norm < (64.0_rb*eps)**2*error) then
          t = abs(a(1,2))
          f = -a(1,1)/a(1,2)
          if (abs(a(2,2)) > t) then
            t = abs(a(2,2))
            f = -a(1,2)/a(2,2)
          end if
          if (abs(a(2,3)) > t) f = -a(1,3)/a(2,3)
          norm = one/sqrt(one + f**2)
          q = [norm, f*norm, zero]

        ! This is the standard branch
        else
          q = q*sqrt(one/norm)
        end if

      end subroutine compute_eigenvector
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine swap( a, b )
        real(rb), intent(inout) :: a, b
        real(rb) :: c
        c = a; a = b; b = c
      end subroutine swap
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine diagonalization

  !=================================================================================================

  pure function cross_product(a, b) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross_product

  !=================================================================================================

  function gaussian_elimination( a, b ) result( c )
    real(rb), intent(in) :: a(:,:), b(:)
    real(rb)             :: c(size(b))
    c = solve_wbs(ge_wpp(a,b))
    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function solve_wbs(u) result(x) ! solve with backward substitution
        real(rb), intent(in) :: u(:,:)
        real(rb)             :: x(size(u,1))
        integer :: i, n
        n = size(u,1)
        forall (i=n:1:-1) x(i) = ( u(i,n+1) - sum(u(i,i+1:n)*x(i+1:n)) ) / u(i,i)
      end function
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function ge_wpp(a, b) result(u) ! gaussian eliminate with partial pivoting
        real(rb), intent(in) :: a(:,:), b(:)
        real(rb)             :: u(size(b),size(b)+1)
        integer  :: i, j, n, p
        real(rb) :: upi
        n = size(a,1)
        u(1:n,1:n) = a
        u(:,n+1) = b
        do j=1,n
          p = maxloc(abs(u(j:n,j)),1) + j-1 ! maxloc returns indices between (1,n-j+1)
          if (p /= j) u([p,j],j) = u([j,p],j)
          u(j+1:,j) = u(j+1:,j)/u(j,j)
          do i=j+1,n+1
            upi = u(p,i)
            if (p /= j) u([p,j],i) = u([j,p],i)
            u(j+1:n,i) = u(j+1:n,i) - upi*u(j+1:n,j)
          end do
        end do
      end function
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function gaussian_elimination

  !=================================================================================================

  logical function colinear(a, b, c)
    real(rb), intent(in) :: a(:), b(size(a)), c(size(a))
    real(rb), parameter :: tol = 1.0E-4
    real(rb) :: u(size(a)), v(size(a))
    u = b - a
    v = c - a
    colinear = abs(sum(u*v)) < tol*sqrt(sum(u*u)*sum(v*v))
  end function colinear

  !=================================================================================================

  function unit_vector( a, b ) result( u )
    real(rb), intent(in) :: a(:)
    real(rb), intent(in), optional :: b(size(a))
    real(rb) :: u(size(a))
    if (present(b)) then
      u = b - a
      u = u/sqrt(sum(u*u))
    else
      u = a/sqrt(sum(a*a))
    end if
  end function unit_vector

  !=================================================================================================

end module mMath
