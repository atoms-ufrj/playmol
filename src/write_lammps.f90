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

  subroutine tPlaymol_write_lammps( me, unit, models )
    class(tPlaymol),  intent(inout) :: me
    integer,          intent(in)    :: unit
    logical,          intent(in)    :: models
    integer :: i, j
    type(StrucHolder) :: atom(me % atom_list % count), bond(me % bond_list % count),    &
                         ang(me % angle_list % count), dih(me % dihedral_list % count), &
                         imp(me % improper_list % count)
    type(TypeHolder), allocatable :: atom_types(:), bond_types(:), &
                                     ang_types(:),  dih_types(:),  &
                                     imp_types(:)
    type tCounter
      integer :: atoms, bonds, angs, dihs, imps, mols
    end type tCounter
    type(tCounter) :: n(me%molecules%N), total(me%molecules%N)
    type(Struc), pointer :: ptr
    integer, allocatable :: mol_index(:)
    ! Molecules:
    n%mols = me % molecules % count()
    ! Atoms:
    call me % analyze_struct( atom, n%atoms, total%atoms, atom_types, &
                              me%atom_list, me%atom_type_list, n%mols, models )
    ! Bonds:
    call me % analyze_struct( bond, n%bonds, total%bonds, bond_types, &
                              me%bond_list, me%bond_type_list, n%mols, models, atom )
    ! Angles:
    call me % analyze_struct( ang, n%angs, total%angs, ang_types, &
                              me%angle_list, me%angle_type_list, n%mols, models, atom )
    ! Dihedrals:
    call me % analyze_struct( dih, n%dihs, total%dihs, dih_types, &
                              me%dihedral_list, me%dihedral_type_list, n%mols, models, atom )
    ! Impropers:
    call me % analyze_struct( imp, n%imps, total%imps, imp_types, &
                              me%improper_list, me%improper_type_list, n%mols, models, atom )
    ! Molecule indices:
    allocate( mol_index(sum(n%mols)) )
    ptr => me % molecules % xyz % first
    do i = 1, size(mol_index)
      mol_index(i) = str2int(me % molecules % list % parameters( ptr % id(1:1) ))
      do j = 1, n(mol_index(i))%atoms
        ptr => ptr % next
      end do
    end do
    ! Write LAMMPS data file:
    write(unit,'("LAMMPS data file",/,"# Generated by Playmol on ",A,/)') trim(now())
    call write_count( size(atom_types), "atom types"     )
    call write_count( size(bond_types), "bond types"     )
    call write_count( size(ang_types),  "angle types"    )
    call write_count( size(dih_types),  "dihedral types" )
    call write_count( size(imp_types),  "improper types" )
    write(unit,'()')
    call write_count( sum(n%mols * total%atoms), "atoms"     )
    call write_count( sum(n%mols * total%bonds), "bonds"     )
    call write_count( sum(n%mols * total%angs),  "angles"    )
    call write_count( sum(n%mols * total%dihs),  "dihedrals" )
    call write_count( sum(n%mols * total%imps),  "impropers" )
    if (me % box % exists()) call write_box_limits
    call write_masses
    if (me % mixing_rule_list % count == 0) then
      call write_type( "Pair Coeffs", atom_types )
    else
      call write_ij_pair_coeffs( atom_types, me % atom_type_list )
    end if
    call write_type( "Bond Coeffs",     bond_types )
    call write_type( "Angle Coeffs",    ang_types  )
    call write_type( "Dihedral Coeffs", dih_types  )
    call write_type( "Improper Coeffs", imp_types  )
    call write_atoms( n%atoms )
    call write_structure( "Bonds",     bond, n%bonds, mol_index, n%atoms )
    call write_structure( "Angles",    ang,  n%angs,  mol_index, n%atoms )
    call write_structure( "Dihedrals", dih,  n%dihs,  mol_index, n%atoms )
    call write_structure( "Impropers", imp,  n%imps,  mol_index, n%atoms )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine write_count( n, name )
        integer,      intent(in) :: n
        character(*), intent(in) :: name
        if (n > 0) then
          write(unit,'(A,X,A)') trim(int2str(n)), name
          call writeln( int2str(n), name )
        end if
      end subroutine write_count
      !---------------------------------------------------------------------------------------------
      subroutine write_box_limits
        real(rb), parameter :: tol = 1.0e-10_rb
        real(rb) :: lim(2) = [-0.5_rb, 0.5_rb]
        real(rb) :: L(3), cos_theta(3), lx, ly, lz, xy, xz, yz
        call me % box % compute( sum(n%mols*me % molecules % per_molecule( me % atom_masses )) )
        L = me%box%length
        write(unit,'()')
        if (me % box % def_type /= 4) then ! Orthogonal box
          write(unit,'(A," xlo xhi")') trim(join(real2str(lim*L(1))))
          write(unit,'(A," ylo yhi")') trim(join(real2str(lim*L(2))))
          write(unit,'(A," zlo zhi")') trim(join(real2str(lim*L(3))))
        else ! Triclinic box
          cos_theta = cosine(me%box%angle)
          lx = L(1)
          xy = L(2)*cos_theta(3); if (abs(xy) < tol) xy = 0.0_rb
          xz = L(3)*cos_theta(2); if (abs(xz) < tol) xz = 0.0_rb
          ly = sqrt(L(2)**2 - xy**2)
          yz = (L(2)*L(3)*cos_theta(1) - xy*xz)/ly; if (abs(yz) < tol) yz = 0.0_rb
          lz = sqrt(L(3)**2 - xz**2 - yz**2)
          write(unit,'(A,X,"xlo xhi")') trim(join(real2str(lim*lx)))
          write(unit,'(A,X,"ylo yhi")') trim(join(real2str(lim*ly)))
          write(unit,'(A,X,"zlo zhi")') trim(join(real2str(lim*lz)))
          write(unit,'(A,X,"xy xz yz")') trim(join(real2str([xy,xz,yz])))
        end if
      end subroutine write_box_limits
      !---------------------------------------------------------------------------------------------
      subroutine write_ij_pair_coeffs( types, list )
        use mMixingRule
        type(TypeHolder),  intent(inout) :: types(:)
        type(StrucList),   intent(in)    :: list
        integer :: i, j, k, m, ntypes, npairs, narg, found, first
        character(sl) :: rule, arg(20), itype(20), jtype(20)
        character(sl), allocatable :: pair(:), model(:)
        type(Struc), pointer :: ptr
        if (size(types) > 0) then
          ntypes = size(types)
          npairs = ntypes*(ntypes - 1)/2 + ntypes
          allocate( pair(npairs), model(npairs) )
          k = 0
          found = 0
          first = merge(2,1,models)
          do i = 1, ntypes
            call split( list % parameters( [types(i)%types] ), narg, itype )
            itype(narg+1:) = ""
            k = k + 1
            found = found + 1
            if (models) model(k) = itype(1)
            pair(k) = join(itype(first:narg))
            do j = i+1, ntypes
              call split( list % parameters( [types(j)%types] ), narg, jtype )
              jtype(narg+1:) = ""
              k = k + 1
              call me % mixing_rule_list % search_exact( [types(i)%types,types(j)%types], ptr )
              if (.not.associated(ptr)) then
                call me % mixing_rule_list % search( [types(i)%types,types(j)%types], ptr )
              end if
              if (associated(ptr)) then
                rule = ptr % params
                found = found + 1
                pair(k) = ""
                call split( rule, narg, arg )
                if (models) model(k) = arg(1)
                do m = first, narg
                  pair(k) = trim(pair(k))//" "//trim(apply_rule( itype(m), jtype(m), arg(m) ))
                end do
              end if
            end do
          end do
          if (found == 0) then
            call write_type( "Pair Coeffs", types )
          else if (found == npairs) then
            if (models) then
              if (any(model(2:) /= model(1))) then
                write(unit,'(/,"PairIJ Coeffs # hybrid",/)')
              else
                write(unit,'(/,"PairIJ Coeffs # ",A,/)') trim(model(1))
                model = ""
              end if
            else
              write(unit,'(/,"PairIJ Coeffs",/)')
              model = ""
            end if
            k = 0
            do i = 1, ntypes
              do j = i, ntypes
                k = k + 1
                write(unit,'(A," # ",A)') trim(join([int2str(i), int2str(j), model(k), pair(k)])), &
                                          trim(join([types(i)%types,types(j)%types]))
              end do
            end do
          else
            call error("not all required mixing rules were defined")
          end if
        end if
      end subroutine write_ij_pair_coeffs
      !---------------------------------------------------------------------------------------------
      subroutine write_type( title, types )
        character(*),      intent(in)    :: title
        type(TypeHolder),  intent(inout) :: types(:)
        integer :: i
        if (size(types) > 0) then
          if (models) then
            if (any(types(2:)%model /= types(1)%model)) then
              write(unit,'(/,A," # hybrid",/)') title
            else
              write(unit,'(/,A," # ",A,/)') title, trim(types(1)%model)
              types%model = ""
            end if
          else
            write(unit,'(/,A,/)') title
            types%model = ""
          end if
          do i = 1, size(types)
            write(unit,'(A," # ",A)') trim(join([int2str(i), types(i)%model, types(i)%params])), &
                                      trim(types(i)%types)
          end do
        end if
      end subroutine write_type
      !---------------------------------------------------------------------------------------------
      subroutine write_masses
        integer :: i
        write(unit,'(/,"Masses",/)')
        do i = 1, size(atom_types)
          write(unit,'(A)') trim(join([int2str(i), atom_types(i)%mass, "#", atom_types(i)%types]))
        end do
      end subroutine write_masses
      !---------------------------------------------------------------------------------------------
      subroutine write_atoms( natoms )
        integer, intent(in) :: natoms(:)
        type(Struc), pointer :: patom
        integer :: i, j, k, m, Na, katom, imol, kmol, prev, ibody
        character(sl) :: cstruc
        integer,  allocatable :: iatom(:), body(:), atom_body(:)
        real(rb), allocatable :: V(:,:), atom_mass(:), type_mass(:)
        character(sl), allocatable :: catom(:), xyz(:)
        Na = sum(n%mols*total%atoms)
        if (Na > 0) then
          write(unit,'(/,"Atoms",/)')
          m = maxval(n%atoms,n%mols > 0)
          allocate( iatom(m), catom(Na), xyz(m), atom_mass(Na), body(Na) )
          type_mass = me % atom_masses % convert_to_real()
          patom => me % molecules % xyz % first
          katom = 0
          do kmol = 1, size(mol_index)
            imol = mol_index(kmol)
            prev = sum(natoms(1:imol-1))
            m = natoms(imol)
            forall (j=1:m) catom(katom+j) = atom(prev+j)%atoms(1)
            do j = 1, m
              iatom(j:j) = pack([(k,k=1,m)],catom(1:m) == patom%id(1))
              xyz(j) = patom%params
              patom => patom % next
            end do
            xyz(1:m) = xyz(sorted(iatom(1:m)))
            do j = 1, m
              katom = katom + 1
              i = prev + j
              cstruc = join(int2str([katom, kmol, atom(i)%itype]))
              write(unit,'(A)') trim(join([cstruc, atom(i)%charge, xyz(j), "#", catom(katom)]))
              atom_mass(katom) = type_mass(atom(i)%itype(1))
              body(katom) = atom(i)%body
            end do
          end do
          if (any(body /= 0)) then
            write(unit,'(/,"BodyTags",/)')
            katom = 0
            ibody = 0
            do kmol = 1, size(mol_index)
              m = natoms(mol_index(kmol))
              allocate( atom_body(m) )
              do i = 1, m
                k = body(katom+i)
                if ((k == 0).or.all(body(katom+1:katom+i-1) /= k)) then
                  ibody = ibody + 1
                  atom_body(i) = ibody
                else
                  j = 1
                  do while (body(katom+j) /= k)
                    j = j + 1
                  end do
                  atom_body(i) = atom_body(j)
                end if
              end do
              do i = 1, m
                katom = katom + 1
                write(unit,'(A)') trim(join([int2str([katom,atom_body(i)]),"#",catom(katom)]))
              end do
              deallocate( atom_body )
            end do
          end if
          if (me%velocity%active) then
            write(unit,'(/,"Velocities",/)')
            V = velocities( sum(n%mols*total%atoms), me%velocity%seed, me%velocity%kT, atom_mass )
            k = 0
            do kmol = 1, size(mol_index)
              do j = 1, natoms(mol_index(kmol))
                k = k + 1
                write(unit,'(A)') trim(join([int2str(k),real2str(V(:,k)), "#", catom(k)]))
              end do
            end do
          end if
        end if
      end subroutine write_atoms
      !---------------------------------------------------------------------------------------------
      subroutine write_structure( title, structure, permol, mol_index, natoms )
        character(*),      intent(in) :: title
        type(StrucHolder), intent(in) :: structure(:)
        integer,           intent(in) :: permol(me%molecules%N), mol_index(:), natoms(:)
        integer :: i, j, k, kmol, imol, sprev, aprev
        character(sl) :: cstruc
        if (any(structure%multiplicity > 0)) then
          write(unit,'(/,A,/)') title
          k = 0
          aprev = 0
          do kmol = 1, size(mol_index)
            imol = mol_index(kmol)
            sprev = sum(permol(1:imol-1))
            do i = 1, permol(imol)
              associate (s => structure(sprev+i))
                do j = 1, s%multiplicity
                  k = k + 1
                  cstruc = join(int2str([k, s%itype(j), aprev + s%atom_in_mol]))
                  write(unit,'(A)') trim(join([cstruc, "#", s%atoms]))
                end do
              end associate
            end do
            aprev = aprev + natoms(imol)
          end do
        end if
      end subroutine write_structure
      !---------------------------------------------------------------------------------------------
  end subroutine tPlaymol_write_lammps

