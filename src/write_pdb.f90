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

  subroutine tPlaymol_write_pdb( me, unit )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    integer :: iatom, i, j, imol, jmol, narg, n
    integer :: natoms(me%molecules%N)
    character(2)  :: element
    character(sl) :: limits, atom_type(1), xyz(3)
    character(sl), allocatable :: atom_id(:)
    integer,       allocatable :: atom_index(:)
    type(Struc), pointer :: current, atom_element, bond

    if (me % box % exists()) then
      write(unit,'("CRYST1",3F9.3,3F7.2)') me%box%length, me%box%angle
    else
      call warning( "created PDB file does not contain box info" )
    end if

    natoms = me % molecules % number_of_atoms()
    iatom = 0
    jmol = 0
    current => me % molecules % xyz % first
    do while (associated(current))
      jmol = jmol + 1
      imol = str2int(me % molecules % list % parameters( current % id ) )
      do i = 1, natoms(imol)
        iatom = iatom + 1
        call split( me % atom_list % parameters( current % id ), narg, atom_type )
        call split( current%params, narg, xyz )
        element = me % element_list % parameters( atom_type )
        write(unit,'("HETATM",I5,X,A4,X,I3,2X,I4,4X,3F8.3,"  1.0   0.0 ",A2)') &
          iatom, adjustl(atom_type(1)(1:4)), jmol, jmol, (str2real(xyz(j)),j=1,3), element
        current => current % next
      end do
    end do

    iatom = 0
    allocate( atom_id(maxval(natoms)), atom_index(maxval(natoms)) )
    current => me % molecules % xyz % first
    do while (associated(current))
      imol = str2int(me % molecules % list % parameters( current % id ) )
      n = natoms(imol)
      do i = 1, n
        iatom = iatom + 1
        atom_id(i) = current % id(1)
        atom_index(i) = iatom
        current => current % next
      end do
      bond => me % bond_list % first
      do while (associated(bond))
        do i = 1, n
          if (atom_id(i) == bond % id(1)) then
            do j = 1, n
              if (atom_id(j) == bond % id(2)) exit
            end do
            write(unit,'("CONECT",2I5)') atom_index(i), atom_index(j)   
            exit
          end if
        end do
        bond => bond % next
      end do
    end do
    write(unit,'("END")')
  end subroutine tPlaymol_write_pdb
