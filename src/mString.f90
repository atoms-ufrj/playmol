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

module mString

use mGlobal

implicit none

character(2), parameter, private :: delimiters = achar(32)//achar(9)
character,    parameter, private :: comment_mark = "#"

interface join
  module procedure :: join_with_space
  module procedure :: join_with_sep
end interface

contains

  !=================================================================================================

  elemental function is_variable( str ) result( yes )
    character(*), intent(in) :: str
    logical                  :: yes
    character :: C
    integer :: i, N
    i = 0
    N = len(str)
    yes = N > 0
    do while (yes.and.(i < N))
      i = i + 1
      C = str(i:i)
      yes = ((C >= "A").and.(C <= "Z")) .or. &
            ((C >= "a").and.(C <= "z")) .or. &
            ((C >= "0").and.(C <= "9")) .or. &
            (C == "_")
    end do
    if (yes) yes = (str(1:1) >= "A").and.(str(1:1) /= "_")
  end function is_variable

  !=================================================================================================

  elemental function has_macros( str ) result( yes )
    character(sl), intent(in) :: str
    logical                   :: yes
    yes = scan(trim(str),"*?") > 0
  end function has_macros

  !=================================================================================================

  elemental function replace_macro( a ) result( b )
   character(*), intent(in) :: a
   character(len(a))        :: b
   integer :: i
   b = ""
   do i = 1, len_trim(a)
     if (a(i:i) == "?") then
       b = trim(b)//"x"
     else if (a(i:i) /= "*") then
       b = trim(b)//a(i:i)
     end if
   end do
  end function replace_macro

  !=================================================================================================

  function replace( string, a, b ) result( new )
    character(*), intent(in)    :: string
    character,    intent(in)    :: a, b
    character(len_trim(string)) :: new
    integer :: k
    new = trim(string)
    forall (k=1:len_trim(new),new(k:k)==a) new(k:k) = b
  end function replace

  !=================================================================================================

  elemental function match_str( a, b ) result( match )
    character(*), intent(in) :: a, b
    logical                  :: match
    logical :: has_macros_a
    has_macros_a = has_macros(a)
    if (has_macros_a.and.has_macros(b)) then
      match = match_wild( a, replace_macro( b ) ).or.match_wild( b, replace_macro( a ) )
    else if (has_macros_a) then
      match = match_wild( a, b )
    else
      match = match_wild( b, a )
    end if
  end function match_str

  !=================================================================================================

  function str_find( names, list ) result( index )
    character(*), intent(in) :: names(:), list(:)
    integer                  :: index(size(names))
    integer :: i, j, n
    logical :: found
    n = size(list)
    do i = 1, size(names)
      j = 0
      found = .false.
      do while (.not.found .and. (j < n))
        j = j + 1
        found = match_str(names(i),list(j))
      end do
      if (found) then
        index(i) = j
      else
        index(i) = 0
      end if
    end do
  end function str_find

  !=================================================================================================

  subroutine clean( str )
    character(*), intent(inout) :: str
    integer :: m, n
    m = scan(trim(str),comment_mark)
    if (m > 0) str = str(1:m-1)
    m = verify(str,delimiters)
    if (m == 0) then
      str = ""
    else
      n = verify(str,delimiters,back=.true.)
      str = str(m:n)
    end if
  end subroutine clean

  !=================================================================================================

  subroutine split( str, narg, arg )
    character(*), intent(in)  :: str
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    logical :: letter, word
    integer :: i, wlen
    narg = 0
    wlen = 0
    word = .false.
    do i = 1, len_trim(str)
      letter = scan(str(i:i),delimiters) == 0
      if (word) then
        if (letter) then
          wlen = wlen + 1
          arg(narg)(wlen:wlen) = str(i:i)
        else
          arg(narg) = arg(narg)(1:wlen)
          if (narg == size(arg)) return
          word = .false.
        end if
      else
        if (letter) then
          narg = narg + 1
          wlen = 1
          arg(narg)(wlen:wlen) = str(i:i)
          word = .true.
        end if
      end if
    end do
    if (word) arg(narg) = arg(narg)(1:wlen)
  end subroutine split

  !=================================================================================================

  function str2cchar( str ) result( array )
    use iso_c_binding, only : c_char
    character(*), intent(in) :: str
    character(kind=c_char)   :: array(len(str))
    integer :: i
    forall(i=1:len(str)) array(i) = str(i:i)
  end function str2cchar

  !=================================================================================================

  function join_with_space( arg ) result( str )
    character(*), intent(in) :: arg(:)
    character(sl)            :: str
    integer :: i, narg
    narg = size(arg)
    if (narg == 0) then
      str = ""
    else
      str = arg(1)
      do i = 2, narg
        str = trim(str)//" "//arg(i)
      end do
    end if
  end function join_with_space

  !=================================================================================================

  function join_with_sep( arg, sep ) result( str )
    character(*), intent(in) :: arg(:)
    character,    intent(in) :: sep
    character(sl)            :: str
    integer :: i, narg
    narg = size(arg)
    if (narg == 0) then
      str = ""
    else
      str = arg(1)
      do i = 2, narg
        str = trim(str)//sep//arg(i)
      end do
    end if
  end function join_with_sep

  !=================================================================================================

  subroutine str_swap( a, b )
    character(sl), intent(inout) :: a, b
    character(sl) :: aux
    aux = a
    a = b
    b = aux
  end subroutine str_swap

  !=================================================================================================

  subroutine read_command( unit, command )
    integer,      intent(in)  :: unit
    character(*), intent(out) :: command
    command = ""
    call add_string( command )
    command = adjustl(command)
    contains
      recursive subroutine add_string( command )
        character(*), intent(inout) :: command
        integer       :: ioerr, last
        character(sl) :: line
        read(unit,'(A'//csl//')',iostat=ioerr) line
        call clean( line )
        do while ((ioerr == 0).and.(line == ""))
          read(unit,'(A'//csl//')',iostat=ioerr) line
          call clean( line )
        end do
        if (ioerr == 0) then
          last = len_trim(line)
          if (line(last:last) == "&") then
            command = trim(command)//" "//line(1:last-1)
            call add_string( command )
          else
            command = trim(command)//" "//line(1:last)
          end if
        end if
      end subroutine add_string
  end subroutine read_command

  !=================================================================================================

  function str2int( str ) result( i )
    character(*), intent(in) :: str
    integer                  :: i
    integer :: ioerr
    read(str,*,iostat=ioerr) i
    if (ioerr /= 0) call error( str, "is not a valid integer number" )
  end function str2int

  !=================================================================================================

  function str2real( str ) result( r )
    character(*), intent(in) :: str
    real(rb)                 :: r
    integer :: ioerr
    read(str,*,iostat=ioerr) r
    if (ioerr /= 0) call error( str, "is not a valid real number" )
  end function str2real

  !=================================================================================================

  elemental function int2str( i ) result( str )
    integer, intent(in) :: i
    character(sl)       :: str
    write(str,*) i
    str = adjustl(str)
  end function int2str

  !=================================================================================================

  elemental function real2str( a ) result( str )
    real(rb), intent(in) :: a
    character(sl)        :: str
    real(4) :: b
    b = a
    if (abs(b) < epsilon(b)) then
      write(str,*) 0.0_4
    else
      write(str,*) b
    end if
    str = adjustl(str)
  end function real2str

  !=================================================================================================

  elemental function float2str( a ) result( str )
    real(rb), intent(in) :: a
    character(sl)        :: str
    write(str,*) real(a,4)
    str = adjustl(str)
  end function float2str

  !=================================================================================================

  elemental function is_int( arg ) result( ok )
    character(*), intent(in) :: arg
    logical                  :: ok
    integer :: ioerr, i
    read(arg,*,iostat=ioerr) i
    ok = ioerr == 0
  end function is_int

  !=================================================================================================

  elemental function is_real( arg ) result( ok )
    character(sl), intent(in) :: arg
    logical                   :: ok
    integer  :: ioerr
    real(rb) :: r
    read(arg,*,iostat=ioerr) r
    ok = ioerr == 0
  end function is_real

  !=================================================================================================

  elemental logical function match_wild(pattern, string)

    ! Compare given string for match to pattern which may contain wildcard characters:
    ! "?" matching any one character, and
    ! "*" matching any zero or more characters.
    ! Both strings may have trailing spaces which are ignored.
    ! Authors: clive page, userid: cgp  domain: le.ac.uk, 2003 (original code)
    !          rolf sander, 2005 (bug fixes and pattern preprocessing)
    ! Minor bug fixed by clive page, 2005 nov 29, bad comment fixed 2005 dec 2.
    ! Serious bug fixed by  robert h mcclanahan, 2011 april 11th
    !
    ! * Reproduced under the terms of the GNU GPL version 3

    implicit none

    character(len=*), intent(in) :: pattern ! pattern may contain * and ?
    character(len=*), intent(in) :: string  ! string to be compared
    integer :: lenp, lenp2, lens, n, p2, p, s
    integer :: n_question, n_asterisk
    logical :: found

    character(len=len(pattern)) :: pattern2
    lens = len_trim(string)
    lenp = len_trim(pattern)

    ! if the pattern is empty, always return true
    if (lenp == 0) then
      match_wild = .true.
      return
    endif

    ! the pattern must be preprocessed. all consecutive occurrences of
    ! one or more question marks ('?') and asterisks ('*') are sorted and
    ! compressed. the result is stored in pattern2.

    pattern2(:)=''
    p  = 1 ! current position in pattern
    p2 = 1 ! current position in pattern2
    do
      if ((pattern(p:p) == '?').or.(pattern(p:p) == '*')) then
        ! a special character was found in the pattern
        n_question = 0
        n_asterisk = 0
        do while (p <= lenp)
          ! count the consecutive question marks and asterisks
          if ((pattern(p:p) /= '?').and.(pattern(p:p) /= '*')) exit
          if (pattern(p:p) == '?') n_question = n_question + 1
          if (pattern(p:p) == '*') n_asterisk = n_asterisk + 1
          p = p + 1
        enddo
        if (n_question>0) then ! first, all the question marks
          pattern2(p2:p2+n_question-1) = repeat('?',n_question)
          p2 = p2 + n_question
        endif
        if (n_asterisk>0) then ! next, the asterisk (only one!)
          pattern2(p2:p2) = '*'
          p2 = p2 + 1
        endif
      else
      ! just a normal character
        pattern2(p2:p2) = pattern(p:p)
        p2 = p2 + 1
        p = p + 1
      endif
      if (p > lenp) exit
    enddo
    !!   lenp2 = p2 - 1
    lenp2 = len_trim(pattern2)

    ! the modified wildcard in pattern2 is compared to the string:
    p2 = 1
    s = 1
    match_wild = .false.
    do
      if (pattern2(p2:p2) == '?') then
        ! accept any char in string
        p2 = p2 + 1
        s = s + 1
      elseif (pattern2(p2:p2) == "*") then
        p2 = p2 + 1
        if (p2 > lenp2) then
          ! anything goes in rest of string
          match_wild = .true.
          exit ! .true.
        else
          ! search string for char at p2
          n = index(string(s:), pattern2(p2:p2))
          if (n == 0) exit  ! .false.
          s = n + s - 1
        endif
      elseif (pattern2(p2:p2) == string(s:s)) then
        ! single char match
        p2 = p2 + 1
        s = s + 1
     else
        ! non-match
        !       exit ! .false.
        ! previous line buggy because failure to match one character in the pattern
        ! does not mean that a match won't be found later. back up through pattern string
        ! until first wildcard character is found and start over with the exact character
        ! match. if the end of the string is reached, then return .false.
        !      04/11/2011 robert mcclanahan    robert.mcclanahan   <<at>>   aecc.com
        found = .false.
        do while (p2 > 0 .and. .not. found)
          p2 = p2 - 1
          if (p2 == 0) exit  !  .false.
          if (pattern(p2:p2) == '*' .or. pattern(p2:p2) == '?') found = .true.
        end do
        s = s + 1
      endif

      if (p2 > lenp2 .and. s > lens) then
        ! end of both pattern2 and string
        match_wild = .true.
        exit ! .true.
      endif

      if (s > lens .and. p2 == lenp) then
        if(pattern2(p2:p2) == "*") then
          ! "*" at end of pattern2 represents an empty string
          match_wild = .true.
          exit
        endif
      endif

      if (p2 > lenp2 .or. s > lens) then
        ! end of either pattern2 or string
        exit ! .false.
      endif
    enddo

  end function match_wild

  !=================================================================================================

  function now() result( string )
    character(sl) :: string
    character(3) :: month(12) = ['Jan','Feb','Mar','Apr','May','Jun', &
                                 'Jul','Aug','Sep','Oct','Nov','Dec']
    integer :: value(8)
    character(2) :: day, hour, min, sec
    character(4) :: year
    call date_and_time( values = value )
    write(day,'(I2)') value(3)
    write(year,'(I4)') value(1)
    write(hour,'(I2)') value(5)
    write(min,'(I2)') value(6); if (min(1:1) == " ") min(1:1) = "0"
    write(sec,'(I2)') value(7); if (sec(1:1) == " ") sec(1:1) = "0"
    string = trim(adjustl(day))//"-"//month(value(2))//"-"//year//" at "//hour//":"//min//":"//sec
  end function now

  !=================================================================================================

  function alphabetical( p, q ) result(lexicalLess)
    character (len = *), intent (in) :: p, q
    logical :: lexicalLess
    integer :: kq, k
    kq = 1
    do k = 1, max(len_trim(p), len_trim(q))
      if (UpperCase(p(k:k)) == UpperCase(q(k:k)) ) then
        cycle
      else
        kq = k
        exit
      end if
    end do
    lexicalLess = UpperCase(p(kq:kq)) < UpperCase(q(kq:kq))
  end function alphabetical

  !=================================================================================================

  function UpperCase(letter) result(L)
    character (len = *), intent (in) :: letter
    character (len = 1) :: L
    character (len = 26), parameter :: Lower = "abcdefghijklmnopqrstuvwxyz", &
       Upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    integer :: k
    k = index(Lower, letter)
    if (k > 0) then
      L = Upper(k:k)
    else
      L = letter
    end if
  end function UpperCase

  !=================================================================================================

end module mString
