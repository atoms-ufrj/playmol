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

module mFix

use mGlobal

implicit none

type tFix
  character(sl) :: prefix = ""
  character(sl) :: suffix = ""
  contains
    procedure :: change => tFix_change
    procedure :: apply => tFix_apply
end type tFix

contains

  !=================================================================================================

  subroutine tFix_change( me, string, kind )
    class(tFix),   intent(out) :: me
    character(sl), intent(in)  :: string, kind
    character(sl) :: empty
    empty = ""
    select case (kind)
      case ("prefix"); me % prefix = merge(string,empty,string /= "none")
      case ("suffix"); me % suffix = merge(string,empty,string /= "none")
      case default; call error( "kind must be 'prefix' or 'suffix'" )
    end select
  end subroutine tFix_change

  !=================================================================================================

  elemental subroutine tFix_apply( me, string )
    class(tFix),   intent(in)    :: me
    character(sl), intent(inout) :: string
    string = trim(me % prefix) // trim(string) // trim(me % suffix)
  end subroutine tFix_apply

  !=================================================================================================

end module mFix
