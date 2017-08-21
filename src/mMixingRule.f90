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

module mMixingRule

use mGlobal
use mString

implicit none

contains

  !=================================================================================================

  character(sl) function apply_rule( a, b, rule ) result( mix )
    character(sl), intent(in) :: a, b, rule
    if ((rule == "arithmetic").or.(rule == "geometric")) then
      if ((.not.is_real(a)) .or. (.not.is_real(b))) then
        call error( "cannot apply mixing rule ", rule, "to parameters", &
                    '"'//trim(a)//'" and "'//trim(b)//'"' )
      end if
    end if
    select case (rule)
      case ("arithmetic")
        mix = real2str(0.5_rb*(str2real(a) + str2real(b)))
      case ("geometric")
        mix = real2str(sqrt(str2real(a)*str2real(b)))
      case default
        mix = rule
      end select
  end function apply_rule

  !=================================================================================================

end module mMixingRule
