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
    select case (rule)
      case ("arithmetic")
        
        mix = real2str(0.5_rb*(numeric(a) + numeric(b)))
      case ("geometric")
        mix = real2str(sqrt(numeric(a)*numeric(b)))
      case default
        mix = rule
      end select

      contains
        function numeric(str) result(val)
          character(sl), intent(in) :: str
          real(rb)                  :: val
          if (str == "") then
            val = 0.0
          else if (is_real(str)) then
            val = str2real(str)
          else
            call error( "cannot apply mixing rule ", rule, "with parameters", a, "and", b )
          end if
        end function numeric 
  end function apply_rule

  !=================================================================================================

end module mMixingRule
