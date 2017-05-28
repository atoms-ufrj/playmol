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

  subroutine tPlaymol_write_summary( me, output )
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: output
    integer :: i, unit, imol, molcount(me%molecules%N), natoms(me%molecules%N)
    integer, allocatable :: output_unit(:)
    real(rb) :: mass(me%molecules%N), charge(me%molecules%N)
    type(Struc), pointer :: ptr
    character(sl) :: atoms
    natoms = me % molecules % number_of_atoms()
    molcount = me % molecules % count()
    charge = me % molecules % per_molecule( me % charge_list )
    mass = me % molecules % per_molecule( me % atom_masses )
    if ((output == stdout).and.(logunit /= 0)) then
      allocate( output_unit(2) )
      output_unit(1) = stdout
      output_unit(2) = logunit
    else
      allocate( output_unit(1) )
      output_unit(1) = output
    end if
    do i = 1, size(output_unit)
      unit = output_unit(i)
      write(unit,'(A,/,"SUMMARY",/,A,/,"Specified:")') repeat("-",80), repeat("-",80)
      call flush_data( unit, "atom type", me % atom_type_list % count )
      call flush_data( unit, "bond type", me % bond_type_list % count )
      call flush_data( unit, "angle type", me % angle_type_list % count )
      call flush_data( unit, "dihedral type", me % dihedral_type_list % count )
      call flush_data( unit, "improper type", me % improper_type_list % count )
      call flush_data( unit, "atom", me % atom_list % count )
      call flush_data( unit, "bond", me % bond_list % count )
      call flush_data( unit, "virtual link", me % link_list % count )
      call flush_data( unit, "improper", me % improper_list % count )
      call flush_data( unit, "body", me % body_list % count )
      write(unit,'(/,A)') "Detected:"
      call flush_data( unit, "angle", me % angle_list % count )
      call flush_data( unit, "dihedral", me % dihedral_list % count )
      call flush_data( unit, "molecule", me%molecules%N )
      write(unit,'(/,A)') "Effectively used:"
      call flush_data( unit, "atom type", me % atom_type_list % count_used() )
      call flush_data( unit, "bond type", me % bond_type_list % count_used() )
      call flush_data( unit, "angle type", me % angle_type_list % count_used() )
      call flush_data( unit, "dihedral type", me % dihedral_type_list % count_used() )
      call flush_data( unit, "improper type", me % improper_type_list % count_used() )
      do imol = 1, me%molecules%N
        write(unit,'(/,"Molecule[",A,"]:")') trim(int2str(imol))
        write(unit,'("- Amount: ",A)') trim(int2str(molcount(imol)))
        write(unit,'("- Mass: ",A)') trim(real2str(mass(imol)))
        write(unit,'("- Charge: ",A)') trim(real2str(charge(imol)))
        write(unit,'("- Number of atoms: ",A)') trim(int2str(natoms(imol)))
        atoms = "- Atoms:"
        ptr => me % molecules % list % first
        do while (associated(ptr))
          if (str2int(ptr % params) == imol) then
            if (len_trim(atoms) + len_trim(ptr % id(1)) > 79) then
              write(unit,'(A)') trim(atoms)
              atoms = ""
            end if
            atoms = trim(atoms)//" "//ptr%id(1)
          end if
          ptr => ptr % next
        end do
        write(unit,'(A)') trim(atoms)
      end do
      write(unit,'(/,"System:")')
      write(unit,'("- Molecules: ",A)') trim(int2str(sum(molcount)))
      write(unit,'("- Total mass: ",A)') trim(real2str(sum(molcount*mass)))
      if (me % box % exists()) then
        call me % box % compute( sum(molcount*mass) )
        write(unit,'("- Box lengths: ",A)') trim(join(real2str(me % box % length)))
        if (me % box % def_type == 4) then
          write(unit,'("- Box angles: ",A)') trim(join(real2str(me % box % angle)))
        end if
        write(unit,'("- Box density: ",A)') trim(real2str(me % box % density))
      end if
      write(unit,'("- Residual charge: ",A)') trim(real2str(sum(molcount*charge)))
      write(unit,'(A)') repeat("-",80)
    end do
    contains
      subroutine flush_data( unit, title, amount )
        integer,      intent(in) :: unit
        character(*), intent(in) :: title
        integer,      intent(in) :: amount
        if (amount == 1) then
          write(unit,'(I5,X,A,".")') amount, title
        else if (amount > 1) then
          write(unit,'(I5,X,A,"s.")') amount, title
        end if
      end subroutine flush_data
  end subroutine tPlaymol_write_summary

