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

  subroutine tPlaymol_write_openmm( me, unit, keywords )
    class(tPlaymol), intent(inout)        :: me
    integer,         intent(in)           :: unit
    character(*),    intent(in)           :: keywords

    integer :: imol, nmols(me % molecules % N), natoms(me % molecules % N)
    real(rb) :: length, energy
    logical :: guess, water(me % molecules % N)

    call process( keywords )

    write(unit,'("<ForceField>")')
    call write_atom_types()

    write(unit,'("  <Residues>")')
    water = me % is_water()
    nmols = me % molecules % count()
    natoms = me % molecules % number_of_atoms()
    do imol = 1, me % molecules % N
      if (nmols(imol) > 0) call write_residue( imol, natoms(imol) )
    end do
    write(unit,'("  </Residues>")')

    call write_bond_types()
    write(unit,'("</ForceField>")')

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine process( keywords )
        character(*), intent(in) :: keywords
        integer :: i, narg
        character(sl) :: keyword, value, arg(40)
        call split( keywords, narg, arg)
        if (mod(narg, 2) == 1) call error( "invalid write openmm command" )
        do i = 1, narg/2
          keyword = arg(2*i-1)
          value = arg(2*i)
          select case (keyword)
            case ("length")
              length = str2real(value)
            case ("energy")
              energy = str2real(value)
            case ("elements")
              if (.not.any(value == ["yes", "no "])) call error( "invalid write openmm command" )
              guess = (value == "yes")
            case default
              call error( "invalid write openmm command" )
          end select
        end do
      end subroutine process
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine write_atom_types()
        character(sl) :: atom_type, name, class, element, mass
        type(Struc), pointer :: current

        write(unit,'(2X,"<AtomTypes>")')
        current => me % atom_type_list % first
        do while (associated(current))
          if (current % usable) then
            atom_type = current%id(1)
            name = item("name", atom_type)
            class = item("class", atom_type)
            call me % element_and_mass( atom_type, element, mass )
            if ((element == "EP").or.((element == "UA").and.(.not.guess))) then
              element = ""
            else if (guess) then
              element = item("element", element_guess( mass ))
            else
              element = item("element", element)
            end if
            mass = item("mass", mass)
            write(unit,'(4X,"<Type ",A,"/>")') trim(join([name, class, element, mass]))
          end if
          current => current % next
        end do
        write(unit,'(2X,"</AtomTypes>")')
      end subroutine write_atom_types
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine write_bond_types()
        integer :: narg
        real(rb) :: K, r0
        character(sl) :: arg(20)
        type(Struc), pointer :: current

        write(unit,'(2X,"<HarmonicBondForce>")')
        current => me % bond_type_list % first
        do while (associated(current))
          if (current % usable) then
            call split( current%params, narg, arg )
            if ((narg == 2).and.all(is_real(arg(1:2)))) then
              K = 2.0_rb*str2real(arg(1)) * energy/length**2
              r0 = str2real(arg(2)) * length
            else if (arg(1) == "harmonic") then
              K = 2.0_rb*str2real(arg(2)) * energy/length**2
              r0 = str2real(arg(3)) * length
            else if (arg(1) == "zero") then
              current => current % next
              cycle
            else
              call error( "bond model", arg(1), "not implemented" )
            end if
            call write_items(4, "Bond", [character(sl) :: "type1", "type2", "length", "k"], &
                                        [current%id, real2str([r0, K])])
          end if
          current => current % next
        end do
        write(unit,'(2X,"</HarmonicBondForce>")')
      end subroutine write_bond_types
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine write_residue( imol, natoms )
        integer, intent(in) :: imol, natoms

        integer :: i, indx
        character(sl) :: string, imol_c
        character(sl) :: atom_type(natoms), charge(natoms), element(natoms), xyz(3)
        integer, allocatable :: pos(:)
        character(sl), allocatable :: atom(:), raw_atom(:)
        type(Struc), pointer :: current, ptr

        imol_c = int2str(imol)
        allocate( atom(natoms), raw_atom(natoms) )

        if (water(imol)) then
          string = "HOH"
        else
          string = letterCode( imol - count(water(1:imol-1)) )
        end if
        write(unit,'(4X,"<Residue ",A,">")') trim(item("name", string))

        i = 0
        current => me % molecules % list % first
        do while (associated(current))
          if (current % params == imol_c) then
            i = i + 1
            atom(i) = current % id(1)
            charge(i) = me % charge_list % parameters( current%id, default = "0" )
            call me % atom_list % search( current%id, ptr, indx )
            atom_type(i) = ptr%params
            call me % element_and_mass( ptr%params, element(i), string )
            ptr => me % raw_atom_list % point_to( indx )
            raw_atom(i) = ptr%params
          end if
          current => current % next
        end do

        ! Atoms:
        do i = 1, natoms
          call write_items(6, "Atom", [character(sl) :: "name", "type", "charge"], &
                                      [raw_atom(i), atom_type(i), charge(i)])
        end do

        ! Virtual sites:
        do i = 1, natoms
          if (element(i) == "EP") call virtual_site(i, atom, raw_atom)
        end do

        ! Bonds:
        atom = pack(atom, element /= "EP")
        raw_atom = pack(raw_atom, element /= "EP")
        current => me % bond_list % first
        do while (associated(current))
          pos = pack([(i,i=1,size(atom))], (atom == current%id(1)).or.(atom == current%id(2)))
          if (size(pos) == 2) then
            call write_items(6, "Bond", ["atomName1", "atomName2"], raw_atom(pos))
          end if
          current => current % next
        end do

        write(unit,'(4X,"</Residue>")')
      end subroutine write_residue
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine virtual_site( i, atom, raw_atom )
        integer, intent(in) :: i
        character(sl), intent(in) :: atom(:), raw_atom(:)

        real(rb), parameter :: tol = 1.0E-4_rb

        integer :: j, k, n
        character(sl) :: string
        integer :: partner(3)
        real(rb) :: a(3,3), b(3), w(3), axis(3,3)
        character(sl) :: xyz(3), average(3)
        integer, allocatable :: pos(:)
        character(sl), allocatable :: properties(:), values(:)
        type(Struc), pointer :: current, ptr

        string = me % molecules % xyz % parameters( [atom(i)] )
        call split(string, k, xyz)
        b = [(str2real(xyz(k)),k=1,3)]
        n = 0
        current => me % link_list % first
        do while (associated(current).and.(n < 4))
          pos = pack([2,1], current%id == atom(i))
          if (size(pos) > 0) then
            n = n + 1
            string = me % molecules % xyz % parameters( [current%id(pos(1))] )
            call split(string, k, xyz)
            a(:,n) = [(str2real(xyz(k)),k=1,3)]
            partner(n:n) = pack([(k,k=1,size(atom))], atom == current%id(pos(1)))
          end if
          current => current % next
        end do
        if (n == 2) then
          ! Test for colinearity:
          if (colinear(a(:,1), a(:,2), b)) then
            w(1:2) = gaussian_elimination( a(1:2,1:2), b(1:2) )
            average(1:2) = [character(sl) :: "1", "2"]
          else
            call error( "VirtualSite type average2 requires colinearity")
          end if
        else if (n == 3) then
          axis(:,1) = unit_vector(a(:,1), a(:,2))
          axis(:,2) = unit_vector(a(:,1), a(:,3))
          axis(:,3) = cross_product(axis(:,1), axis(:,2))
          w = gaussian_elimination( axis, b - a(:,1) )
          if (abs(w(3)) < tol) then
            w = gaussian_elimination( a, b )
            average = [character(sl) :: "1", "2", "3"]
          else
            average = [character(sl) :: "12", "13", "Cross"]
          end if
        else
          call error( "Extra particle must be linked to 2 or 3 atoms" )
        end if
        properties = [character(sl) :: "type", "siteName", &
                      ("atomName"//int2str(k), k=1, n), &
                      ("weight"//average(k), k=1, n)]
        values = [character(sl) :: "average"//average(n), raw_atom(i), &
                  raw_atom(partner), (float2str(w(k)),k=1,3)]
        call write_items(6, "VirtualSite", properties, values)
      end subroutine virtual_site
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine write_items( ident, title, property, value )
        integer,      intent(in) :: ident
        character(*), intent(in) :: title, property(:), value(:)
        integer :: i
        write(unit,'("'//repeat(" ",ident)//'","<",A,X,A,"/>")') trim(title), &
          trim(join([(item(property(i), value(i)), i=1, size(property))]))
      end subroutine write_items
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental character(sl) function item( property, value )
        character(*), intent(in) :: property, value
        item = trim(property)//"="""//trim(value)//""""
      end function item
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      function element_guess( mass ) result( element )
        character(sl), intent(in) :: mass
        character(sl)             :: element
        element = me%elements(minloc(abs(me%masses - str2real(mass)), dim = 1))
      end function element_guess
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine tPlaymol_write_openmm
