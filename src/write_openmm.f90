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

  subroutine tPlaymol_write_openmm( me, unit )
    class(tPlaymol),  intent(inout) :: me
    integer,          intent(in)    :: unit

    integer :: i
    type(Struc), pointer :: current
    character(sl) :: mass, element

    write(unit,'("<ForceField>")')

    write(unit,'("<AtomTypes>")')
    current => me % atom_type_list % first
    do while (associated(current))
      if (current % usable) then
        i = str2int(me % element_list % parameters( current%id ))
        element = me % elements(i)
        mass = me % mass_list % parameters( current%id )
        if (mass == "") mass = real2str(me % masses(i))
        write(unit,'("<Type name=""",A,""" class=""",A,""" element=""",A,""" mass=""",A,"""/>")') &
          trim(current%id(1)), trim(current%id(1)), trim(element), trim(mass)
      end if
      current => current % next
    end do
    write(unit,'("</AtomTypes>")')

    write(unit,'("</ForceField>")')

  end subroutine tPlaymol_write_openmm
