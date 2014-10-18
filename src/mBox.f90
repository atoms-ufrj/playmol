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

module mBox

use mGlobal

implicit none

type tBox
  integer  :: def_type = 0    ! 0 - none
                              ! 1 - density and aspect
                              ! 2 - volume and aspect
                              ! 3 - lengths
  real(rb) :: volume
  real(rb) :: density
  real(rb) :: length(3)
  real(rb) :: aspect(3) = 1.0_rb
  contains
    procedure :: define => tBox_define
    procedure :: exists => tBox_exists
    procedure :: compute => tBox_compute
end type tBox

contains

  !=================================================================================================

  subroutine tBox_define( box, prop, scalar, vector )
    class(tBox),   intent(inout)        :: box
    character(sl), intent(in)           :: prop
    real(rb),      intent(in)           :: scalar
    real(rb),      intent(in), optional :: vector(3)
    if (any(vector <= 0.0_rb)) call error( "invalid value(s)" )
    select case (prop)
      case ("density")
        if (scalar <= 0.0_rb) call error( "invalid value" )
        box % def_type = 1
        box % density = scalar
        box % aspect = vector
      case ("volume")
        if (scalar <= 0.0_rb) call error( "invalid value" )
        box % def_type = 2
        box % density = scalar
        box % aspect = vector
      case ("size")
        box % def_type = 3
        box % length = vector
    end select
  end subroutine tBox_define

  !=================================================================================================

  function tBox_exists( box ) result( exists )
    class(tBox), intent(in) :: box
    logical                 :: exists
    exists = box % def_type /= 0
  end function tBox_exists

  !=================================================================================================

  subroutine tBox_compute( box, mass )
    class(tBox), intent(inout) :: box
    real(rb),    intent(in)    :: mass
    real(rb) :: L
    if (.not.box % exists()) call error( "box has not been specified" )
    select case (box % def_type)
      case (1)
        box % volume = mass / box % density
        L = (box % volume / product(box % aspect))**(1.0_rb/3.0_rb)
        box % length = L * box % aspect
      case (2)
        box % density = mass / box % volume
        L = (box % volume / product(box % aspect))**(1.0_rb/3.0_rb)
        box % length = L * box % aspect
      case (3)
        box % volume = product(box % length)
        box % density = mass / box % volume
        box % aspect = box % length / box % length(1)
    end select
  end subroutine tBox_compute

  !=================================================================================================

end module mBox
