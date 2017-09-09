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

module mStruc

use mGlobal
use mString

implicit none

type Struc
  character(sl), allocatable :: id(:)
  character(sl)              :: params
  type(Struc), pointer       :: next => null()
  logical                    :: usable = .false.
  contains
    procedure :: init => Struc_init
    procedure :: match_id => Struc_match_id
end type Struc

type StrucList
  character(sl)        :: name = ""
  integer              :: number = 1
  logical              :: reversible = .true.
  integer              :: count = 0
  type(Struc), pointer :: first => null()
  type(Struc), pointer :: last  => null()
  contains
    procedure :: add => StrucList_add
    procedure :: handle => StrucList_handle
    procedure :: search => StrucList_search
    procedure :: parameters => StrucList_params_from_id
    procedure :: point_to => StrucList_point_to
    procedure :: index => StrucList_index
    procedure :: find => StrucList_find
    procedure :: count_used => StrucList_count_used
    procedure :: print => StrucList_print
    procedure :: remove => StrucList_remove
    procedure :: convert_to_real => StrucList_convert_to_real
    procedure :: destroy => StrucList_destroy
end type StrucList

contains

  !=================================================================================================

  subroutine Struc_init( me, narg, arg, number )
    class(Struc), intent(out) :: me
    integer,      intent(in)  :: narg
    character(*), intent(in)  :: arg(:)
    integer,      intent(in)  :: number
    allocate( me % id(number) )
    me % id = arg(1:number)
    me % params = join( arg(number+1:narg) )
  end subroutine Struc_init

  !=================================================================================================

  function Struc_match_id( me, id ) result( match )
    class(Struc),  intent(in) :: me
    character(sl), intent(in) :: id(:)
    logical                   :: match
    match = all(match_str( me % id, id ))
    if (size(id) > 1) match = match .or. all(match_str( me % id, id(size(id):1:-1) ))
  end function Struc_match_id

  !=================================================================================================

  subroutine StrucList_add( me, narg, arg, allowedIDs, repeatable, silent )
    class(StrucList),           intent(inout) :: me
    integer,                    intent(in)    :: narg
    character(*),               intent(in)    :: arg(:)
    class(StrucList), optional, intent(in)    :: allowedIDs
    logical,          optional, intent(in)    :: repeatable, silent
    integer :: i, n
    logical :: unique, print
    type(Struc), pointer :: ptr
    if (present(silent)) then
      print = .not.silent
    else
      print = .true.
    end if
    n = me % number
    if (narg == n) then
      if (print) call writeln( "Adding", me%name, join(arg(1:n)) )
    else if (narg > n) then
      if (print) then
        call writeln( "Adding", me%name, join(arg(1:n)), "with parameters", join(arg(n+1:narg)) )
      end if
    else
      call error( "invalid definition for", me%name, join(arg(1:narg)) )
    end if
    if (present(allowedIDs)) then
      do i = 1, n
        if (.not. allowedIDs % find(arg(i:i))) then
          call error( "undefined", allowedIDs%name, arg(i), "in", me%name, "definition" )
        end if
      end do
    end if
    if (present(repeatable)) then
      unique = .not.repeatable
    else
      unique = .true.
    end if
    if (unique) then
      call me % search( arg(1:n), ptr )
      if (associated(ptr)) then
        if (me%reversible) then
          call error( "conflicts with", me%name, join( ptr%id ) )
        else if (all(match_str( arg(1:n), ptr%id ))) then
          call error( "conflicts with", me%name, join( ptr%id ) )
        end if
      end if
    end if
    if (associated(me % last)) then
      allocate( me % last % next )
      me % last => me % last % next
    else
      allocate( me % last )
      me % first => me % last
    end if
    call me % last % init( narg, arg, me % number )
    me % count = me % count + 1
  end subroutine StrucList_add

  !=================================================================================================

  subroutine StrucList_handle( me, id, id_list, type_list, option )
    class(StrucList),           intent(inout) :: me
    character(*),               intent(inout) :: id(:)
    class(StrucList), optional, intent(in)    :: id_list, type_list
    integer,                    intent(in)    :: option
    ! Search for types, mark found types as usable, and:
    ! option = 1: add, WARN if no types were found
    ! option = 2: STOP if no types were found
    ! option = 3: add only if types were found
    integer :: i, n
    type(Struc), pointer :: current
    character(sl) :: types(me%number)
    logical :: direct, reverse, found_direct, found_reverse, found_two_way, found, direct_first

    do i = 1, me%number
      call id_list % search( id(i:i), current )
      if (.not.associated(current)) call error( "unknown", id_list % name, id(i) )
      call split( current % params, n, types(i:i) )
    end do

    current => type_list % first
    found = .false.
    do while (associated(current).and.(.not.found))
      direct = all(match_str( current%id, types ))
      reverse = all(match_str( current%id(me%number:1:-1), types ))
      found = direct .neqv. reverse
      current => current % next
    end do
    direct_first = direct .or. (.not.found)

    current => type_list % first
    found_direct = .false.
    found_reverse = .false.
    do while (associated(current))
      direct = all(match_str( current%id, types ))
      reverse = all(match_str( current%id(me%number:1:-1), types ))
      if (direct_first) then
        if (direct) then
          found_direct = .true.
          current % usable = .true.
          found_two_way = reverse
        else if (reverse) then
          found_reverse = .true.
          current % usable = .true.
        end if
      else
        if (reverse) then
          found_reverse = .true.
          current % usable = .true.
          found_two_way = direct
        else if (direct) then
          found_direct = .true.
          current % usable = .true.
        end if
      end if
      current => current % next
    end do
    found = found_direct .or. found_reverse
    select case (option)
      case (1) ! add, WARN if no types were found
        if (me%reversible) then
          call me % add( me%number, id )
        else
          if (found_direct.or.(.not.found)) call me % add( me%number, id )
          if (found_reverse) call me % add( me%number, id(me%number:1:-1) )
          if (found_direct .and. found_reverse .and. found_two_way) then
            call warning( "direction conflict for", me%name, join(id), "( types", join(types), ")" )
            current => type_list % first
            do while (associated(current))
              if (current%usable) then
                direct = all(match_str( current%id, types ))
                reverse = all(match_str( current%id(me%number:1:-1), types ))
                if (direct .and. reverse) then
                  call writeln( "-->", type_list%name, join(current%id), "( two-way )" )
                else if (direct) then
                  call writeln( "-->", type_list%name, join(current%id), "( direct  )" )
                else if (reverse) then
                  call writeln( "-->", type_list%name, join(current%id), "( reverse )" )
                end if
              end if
              current => current % next
            end do
          end if
        end if
        if (.not.found) then
          call warning("undefined",type_list%name,"for",id_list%name,join(id),"(",join(types),")")
        end if
      case (2) ! STOP if no types were found
        if (.not.found) then
          call error("undefined",type_list%name,"for",id_list%name,join(id),"(",join(types),")")
        end if
      case (3) ! add only if types were found
        if (found) then
          if (me%reversible) then
            call me % add( me%number, id )
          else
            if (found_direct) call me % add( me%number, id )
            if (found_reverse) call me % add( me%number, id(me%number:1:-1) )
          end if
        end if
    end select
  end subroutine StrucList_handle

  !=================================================================================================

  subroutine StrucList_search( me, id, ptr, index )
    class(StrucList),     intent(in)            :: me
    character(*),         intent(in)            :: id(:)
    type(Struc), pointer, intent(out), optional :: ptr
    integer,              intent(out), optional :: index
    logical :: found
    integer :: i
    type(Struc), pointer :: current
    found = .false.
    current => me % first
    i = 0
    do while (associated(current).and.(.not.found))
      i = i + 1
      found = current % match_id( id )
      if (.not.found) current => current % next
    end do
    if (found) then
      if (present(index)) index = i
      if (present(ptr)) ptr => current
    else
      if (present(index)) index = 0
      if (present(ptr)) ptr => null()
    end if
  end subroutine StrucList_search

  !=================================================================================================

  function StrucList_params_from_id( me, id, default ) result( params )
    class(StrucList), intent(in)           :: me
    character(*),     intent(in)           :: id(:)
    character(*),     intent(in), optional :: default
    character(sl)                          :: params
    type(Struc), pointer :: ptr
    call me % search( id, ptr )
    if (associated(ptr)) then
      params = ptr % params
    else
      if (present(default)) then
        params = default
      else
        params = ""
      end if
    end if
  end function StrucList_params_from_id

  !=================================================================================================

  function StrucList_point_to( me, index ) result( ptr )
    class(StrucList), intent(in) :: me
    integer,          intent(in) :: index
    type(Struc), pointer         :: ptr
    integer :: i
    if ((index < 1).or.(index > me%count)) then
      ptr => null()
    else
      ptr => me % first
      do i = 2, index
        ptr => ptr % next
      end do
    end if
  end function StrucList_point_to

  !=================================================================================================

  function StrucList_find( me, id ) result( found )
    class(StrucList), intent(in) :: me
    character(*),     intent(in) :: id(:)
    logical                      :: found
    type(Struc), pointer :: ptr
    call me % search( id, ptr )
    found = associated( ptr )
  end function StrucList_find

  !=================================================================================================

  function StrucList_index( me, id ) result( index )
    class(StrucList), intent(in) :: me
    character(*),     intent(in) :: id(:)
    integer                      :: index
    call me % search( id, index = index )
  end function StrucList_index

  !=================================================================================================

  function StrucList_count( me, valids_only ) result( N )
    class(StrucList), intent(in)  :: me
    integer                       :: N
    logical, intent(in), optional :: valids_only
    type(Struc), pointer :: current
    integer :: Nv
    current => me % first
    N = 0
    Nv = 0
    do while (associated(current))
      N = N + 1
      if (current % usable) Nv = Nv + 1
      current => current % next
    end do
    if (present(valids_only)) N = Nv
  end function StrucList_count

  !=================================================================================================

  function StrucList_count_used( me ) result( N )
    class(StrucList), intent(in) :: me
    integer                      :: N
    type(Struc), pointer :: current
    current => me % first
    N = 0
    do while (associated(current))
      if (current % usable) N = N + 1
      current => current % next
    end do
  end function StrucList_count_used

  !=================================================================================================

  subroutine StrucList_print( me, unit, comment, used_only )
    class(StrucList), intent(in)           :: me
    integer,          intent(in)           :: unit
    logical,          intent(in), optional :: comment, used_only
    type(Struc), pointer :: current
    character(sl) :: name
    logical :: test_used
    name = me % name
    if (present(comment)) then
      if (comment) name = "# "//trim(me % name)
    end if
    test_used = present(used_only)
    if (test_used) test_used = used_only
    current => me % first
    do while (associated(current))
      if ((.not.test_used).or.(test_used.and.current%usable)) then
        write(unit,'(A,X,A,X,A)') trim(name), trim(join(current%id)), trim(current%params)
      end if
      current => current % next
    end do
  end subroutine StrucList_print

  !=================================================================================================

  subroutine StrucList_remove( me, arg, silent )
    class(StrucList), intent(inout)        :: me
    character(sl),    intent(in)           :: arg(me%number)
    logical,          intent(in), optional :: silent
    integer :: i, index
    type(Struc), pointer :: current, aux
    if (present(silent)) then
      if (.not.silent) call writeln( "Removing", me%name, join(arg) )
    else
      call writeln( "Removing", me%name, join(arg) )
    end if
    call me % search( arg, index = index )
    if (index == 0) call error( me%name, join(arg), "does not exist" )
    current => me % first
    if (index == 1) then
      me % first => me % first % next
      deallocate( current )
    else
      do i = 1, index-2
        current => current % next
      end do
      aux => current % next
      current % next => aux % next
      deallocate( aux )
    end if
    me % count = me % count - 1
  end subroutine StrucList_remove

  !=================================================================================================

  function StrucList_convert_to_real( me ) result( value )
    class(StrucList), intent(in) :: me
    real(rb)                     :: value(me%count)
    integer :: i
    type(Struc), pointer :: current
    current => me % first
    i = 0
    do while (associated(current))
      i = i + 1
      value(i) = str2real( current % params )
      current => current % next
    end do
  end function StrucList_convert_to_real

  !=================================================================================================

  subroutine StrucList_destroy( me )
    class(StrucList), intent(inout) :: me
    type(Struc), pointer :: current, aux
    current => me % first
    if (associated(current)) call writeln( "Deleting ", me%name, "list..." )
    do while (associated(current))
      aux => current
      current => current % next
      call writeln( "Deleting ", me%name, join(aux % id) )
      deallocate(aux)
    end do
    me % first => null()
    me % last => null()
    me % count = 0
  end subroutine StrucList_destroy

  !=================================================================================================

end module mStruc
