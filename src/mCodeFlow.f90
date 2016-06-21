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

module mCodeFlow

use mGlobal
use mStruc
use mParser

implicit none

type, private :: tCommand
  character(sl) :: content
  type(tCommand), pointer :: next => null()
end type

type tCodeFlow
  type(tCommand), pointer :: commands => null()
  type(StrucList) :: variable_list = StrucList( name = "variable", number = 1 )
  contains
    procedure :: next_command => tCodeFlow_next_command
end type

contains

  !=================================================================================================

  recursive subroutine tCodeFlow_next_command( me, unit, narg, arg )
    class(tCodeFlow), intent(inout) :: me
    integer,          intent(in)    :: unit
    integer,          intent(out)   :: narg
    character(sl),    intent(inout) :: arg(:)
    character(sl) :: command
    command = next_command_line( unit )
    call replace_variables
    call evaluate_expressions
    call split( command, narg, arg )
    if (narg > 0) then
      select case (arg(1))
        case ("for")
          call evaluate_for_command
          call me % next_command( unit, narg, arg )
        case ("define")
          call evaluate_define_command
          call me % next_command( unit, narg, arg )
        case ("if")
          call evaluate_if_command
          call me % next_command( unit, narg, arg )
        case ("include")
          call evaluate_include_command
          call me % next_command( unit, narg, arg )
      end select
    end if
    contains
      !---------------------------------------------------------------------------------------------
      subroutine evaluate_define_command
        type(Struc), pointer :: ptr
        if ((narg < 4).or.(arg(3) /= "as")) call error( "invalid define command" )
        if (.not.is_variable(trim(arg(2)))) call error( "invalid variable name", arg(2) )
        call me % variable_list % search( arg(2:2), ptr )
        if (associated(ptr)) then
          call writeln( "Redefining variable", join(arg(2:narg)) )
          ptr % params = join(arg(4:narg))
        else
          arg(3) = arg(2)
          call me % variable_list % add( narg-2, arg(3:narg) )
        end if
      end subroutine evaluate_define_command
      !---------------------------------------------------------------------------------------------
      subroutine evaluate_for_command
        integer :: ifirst, ilast, step, open_for
        character(sl) :: variable, first, last, now, after, commandline
        logical :: integers, execute
        type(tCommand), pointer :: a_first, a_curr, b_first, b_curr
        if (narg < 3) call error( "invalid for...next definition" )
        variable = arg(2)
        if (.not.is_variable(trim(variable))) call error( "invalid variable name", variable )
        if (arg(3) == "in") then
          execute = narg > 3
          if (execute) then
            now = "define "//trim(variable)//" as "//trim(arg(4))
            after = join([arg(1:3),arg(5:narg)])
          end if
        else if (arg(3) == "from") then
          if (narg /= 6) call error( "invalid for...next definition" )
          first = arg(4)
          select case (arg(5))
            case ("to"); step = 1
            case ("downto"); step = -1
            case default; call error( "invalid for...next definition" )
          end select
          last = arg(6)
          integers = is_int(first) .and. is_int(last)
          if (integers) then
            ifirst = str2int(first)
            ilast = str2int(last)
            arg(4) = int2str(ifirst + step)
          else if ((len_trim(first) == 1).and.(len_trim(last) == 1)) then
            ifirst = ichar(first(1:1))
            ilast = ichar(last(1:1))
            arg(4) = char(ifirst + step)
          else
            call error( "invalid for...next definition" )
          end if
          execute = step*ifirst <= step*ilast
          if (execute) then
            now = "define "//trim(variable)//" as "//trim(first)
            after = join(arg(1:6))
          end if
        else
          call error( "invalid for...next definition" )
        end if
        if (execute) then
          allocate( a_first, b_first )
          a_first % content = now
          b_first % content = after
          a_curr => a_first
          b_curr => b_first
        end if
        open_for = 1
        do while (open_for > 0)
          commandline = next_command_line( unit )
          if (commandline == "") then
            open_for = -1
          else if (commandline == "next") then
            open_for = open_for - 1
          else if (commandline(1:3) == "for") then
            open_for = open_for + 1
          end if
          if ((open_for > 0) .and. execute) then
            allocate( a_curr % next, b_curr % next )
            a_curr => a_curr % next
            b_curr => b_curr % next
            a_curr % content = commandline
            b_curr % content = commandline
          end if
        end do
        if (open_for == -1) call error( "unfinished for...next loop" )
        if (execute) then
          allocate( b_curr % next )
          b_curr => b_curr % next
          b_curr % content = "next"
          b_curr % next => me % commands
          a_curr % next => b_first
          me % commands => a_first
        end if
      end subroutine evaluate_for_command
      !---------------------------------------------------------------------------------------------
      subroutine evaluate_if_command
        logical :: execute
        type(tCommand), pointer :: first, curr
        integer :: open_if
        character(sl) :: commandline
        logical :: else_found
        if (narg < 3) call error( "invalid if/then/else definition" )
        if (arg(3) /= "then") call error( "invalid if/then/else definition" )
        if (.not.any(arg(2)==["1","0"])) call  error( "invalid if/then/else clause" )
        execute = arg(2) == "1"
        first => null()
        open_if = 1
        else_found = .false.
        do while (.true.)
          commandline = next_command_line( unit )
          if (commandline == "") call error( "unfinished if/then/else construct" )
          if (open_if == 1) then
            if (commandline == "else") then
              if (else_found) call error( "invalid if/then/else definition" )
              execute = .not.execute
              else_found = .true.
              cycle
            else if (commandline == "endif") then
              exit
            end if
          end if
          if (execute) then
            if (associated(first)) then
              allocate( curr % next )
              curr => curr % next
            else
              allocate( first )
              curr => first
            end if
            curr % content = commandline
          end if
          if (commandline(1:2) == "if") then
            open_if = open_if + 1
          else if (commandline == "endif") then
            open_if = open_if - 1
          end if
        end do
        if (associated(first)) then
          curr % next => me % commands
          me % commands => first
        end if
      end subroutine evaluate_if_command
      !---------------------------------------------------------------------------------------------
      subroutine evaluate_include_command
        logical :: file_exists
        integer :: input
        character(sl) :: commandline
        type(tCommand), pointer :: first, curr
        if (narg < 2) call error( "invalid include command" )
        inquire( file = arg(2), exist = file_exists )
        if (.not.file_exists) call error( "file", arg(2), "does not exist" )
        call writeln( "Including file", arg(2) )
        open( newunit = input, file = arg(2), status = "old" )
        call read_command( input, commandline )
        if (commandline /= "") then
          allocate( first )
          curr => first
          curr % content = commandline
          call read_command( input, commandline )
          do while (commandline /= "")
            allocate( curr % next )
            curr => curr % next
            curr % content = commandline
            call read_command( input, commandline )
          end do
          curr % next => me % commands
          me % commands => first
        end if
        close(input)
      end subroutine evaluate_include_command
      !---------------------------------------------------------------------------------------------
      subroutine replace_variables
        integer :: N, first, last
        character(sl) :: vname
        type(Struc), pointer :: ptr
        first = index(trim(command),"$",back=.true.)
        do while (first > 0)
          N = len_trim(command)
          if (first == N) call error( "invalid variable" )
          if (command(first+1:first+1) == "{") then
            last = first + index(command(first+1:N),"}")
            if (last <= first+2) call error( "invalid variable" )
            if (.not.is_variable(command(first+2:last-1))) &
              call error( "invalid variable", command(first+2:last-1) )
            vname = command(first+2:last-1)
          else
            last = first
            do while (is_variable(command(first+1:last+1)) .and. (last < N))
              last = last + 1
            end do
            if (.not.is_variable(command(first+1:last))) &
              call error( "invalid variable", command(first+1:last) )
            vname = command(first+1:last)
          end if
          call me % variable_list % search( [vname], ptr )
          if (associated(ptr)) then
            command = command(1:first-1)//trim(ptr%params)//command(last+1:N)
          else
            call error( "undefined variable", vname )
          end if
          first = index(trim(command),"$",back=.true.)
        end do
      end subroutine replace_variables
      !---------------------------------------------------------------------------------------------
      subroutine evaluate_expressions
        integer :: N, first, last
        real(rb) :: value
        type(tParser) :: Comp
        first = index(trim(command),"{",back=.true.)
        do while (first > 0)
          N = len_trim(command)
          last = first + index(command(first+1:N),"}")
          if (last == first) call error( "unfinished math expression" )
          if (last == first+1) then
            command = command(1:first-1)//command(last+1:N)
          else
            call Comp % parse( command(first+1:last-1) )
            value = Comp % evaluate()
            if (Comp % Is_Integer) then
              command = command(1:first-1)//trim(int2str(nint(value)))//command(last+1:N)
            else
              command = command(1:first-1)//trim(real2str(value))//command(last+1:N)
            end if
          end if
          first = index(trim(command),"{",back=.true.)
        end do
      end subroutine evaluate_expressions
      !---------------------------------------------------------------------------------------------
      function next_command_line( unit ) result( content )
        integer, intent(in)     :: unit
        character(sl)           :: content
        type(tCommand), pointer :: aux
        if (associated(me % commands)) then
          content = me % commands % content
          aux => me % commands
          me % commands => me % commands % next
          deallocate( aux )
        else
          call read_command( unit, content )
        end if
      end function next_command_line
      !---------------------------------------------------------------------------------------------
  end subroutine tCodeFlow_next_command

  !=================================================================================================

end module mCodeFlow
