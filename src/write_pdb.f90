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
    integer :: iatom, i, j, imol, jmol, narg, n, indx
    integer :: natoms(me%molecules%N)
    real(rb) :: HL(3)
    character(sl) :: atom, element, residue, xyz(3), atom_name
    character(sl), allocatable :: atom_id(:)
    integer,       allocatable :: atom_index(:)
    logical,       allocatable :: water(:)
    type(Struc), pointer :: current, bond, ptr

    ! Identify water molecules:
    water = me % is_water()

    ! Write out box info:
    write(unit,'("REMARK Generated by Playmol on ",A)') trim(now())
    if (me % box % exists()) then
      write(unit,'("CRYST1",3F9.3,3F7.2," P 1",11X,"1 ")') me%box%length, me%box%angle
      HL = 0.5_rb*me%box%length
    else
      call warning( "created PDB file does not contain box info" )
      HL = 0.0_rb
    end if

    ! Write out atomic info:
    iatom = 0
    jmol = 0
    current => me % molecules % xyz % first
    natoms = me % molecules % number_of_atoms()
    do while (associated(current))
      jmol = jmol + 1
      imol = str2int(me % molecules % list % parameters( current % id ) )
      if (water(imol)) then
        residue = "HOH"
      else
        residue = letterCode( imol - count(water(1:imol-1)) )
      end if
      do i = 1, natoms(imol)
        iatom = iatom + 1
        call me % atom_list % search( current%id, ptr, indx )
        ptr => me % raw_atom_list % point_to( indx )
        atom_name = ptr%params
        ptr => me % atom_elements % point_to( indx )
        element = ptr%params
        if (element == "UA") element = ""
        call split( current%params, narg, xyz )
        write(unit,'("HETATM",I5,X,A4,X,A3,2X,I4,4X,3F8.3,"  1.0   0.0 ",A2)') &
          iatom, & !Atom serial number
          adjustl(atom_name(1:4)), & ! Atom name
          residue, & ! Residue name
          jmol, & ! Residue sequence number
          (str2real(xyz(j)) + HL(j),j=1,3), & ! Coordinates
          element ! Element symbol
        current => current % next
      end do
    end do

    ! Write out connectivity info:
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
