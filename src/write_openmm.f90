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

  subroutine tPlaymol_write_openmm( me, unit, allow_ua )
    class(tPlaymol),  intent(inout)        :: me
    integer,          intent(in)           :: unit
    logical,          intent(in), optional :: allow_ua

    integer :: i
    logical :: united_atom
    character(sl) :: name, class, element, mass
    type(Struc), pointer :: current

    united_atom = present(allow_ua)
    if (united_atom) united_atom = allow_ua

    write(unit,'("<ForceField>")')

    write(unit,'(2X,"<AtomTypes>")')
    current => me % atom_type_list % first
    do while (associated(current))
      if (current % usable) then
        name = " name = """//trim(current%id(1))//""""
        class = " class = """//trim(current%id(1))//""""
        call me % element_and_mass( current%id(1), element, mass )
        element = " element = """//trim(element)//""""
        mass = " mass = """//trim(mass)//""""
        write(unit,'(4X,"<Type",A,"/>")') trim(name)//trim(class)//trim(element)//trim(mass)
      end if
      current => current % next
    end do
    write(unit,'(2X,"</AtomTypes>")')

    write(unit,'("</ForceField>")')

  end subroutine tPlaymol_write_openmm
