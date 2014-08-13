module mString

use mGlobal

implicit none

character(2), parameter, private :: delimiters = achar(32)//achar(9)
character,    parameter, private :: comment_mark = "#"

contains

  !=================================================================================================

  elemental function has_macros( str ) result( yes )
    character(sl), intent(in) :: str
    logical                   :: yes
    yes = scan(trim(str),"*?") > 0
  end function has_macros

  !=================================================================================================

  elemental function match_str( a, b ) result( match )
    character(*), intent(in) :: a, b
    logical                  :: match
    logical :: has_macros_a
    has_macros_a = has_macros(a)
    if (has_macros_a.and.has_macros(b)) then
      match = trim(a) == trim(b)
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

  function join( arg ) result( str )
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
  end function join

  !=================================================================================================

  subroutine str_swap( a, b )
    character(sl), intent(inout) :: a, b
    character(sl) :: aux
    aux = a
    a = b
    b = aux
  end subroutine str_swap

  !=================================================================================================

  subroutine next_command( unit, narg, arg )
    integer,      intent(in)  :: unit
    integer,      intent(out) :: narg
    character(*), intent(out) :: arg(:)
    integer       :: ioerr
    character(sl) :: line
    read(unit,'(A'//csl//')',iostat=ioerr) line
    call clean( line )
    do while ((ioerr == 0).and.(line == ""))
      read(unit,'(A'//csl//')',iostat=ioerr) line
      call clean( line )
    end do
    if (ioerr == 0) then
      call split( line, narg, arg )
    else
      narg = 0
    end if
  end subroutine next_command

  !=================================================================================================

  function str2int( str ) result( i )
    character(*), intent(in) :: str
    integer                  :: i
    integer :: ioerr
    read(str,*,iostat=ioerr) i
    if (ioerr /= 0) call error( "bad integer" )
  end function str2int

  !=================================================================================================

  function str2real( str ) result( r )
    character(*), intent(in) :: str
    real(rb)                 :: r
    integer :: ioerr
    read(str,*,iostat=ioerr) r
    if (ioerr /= 0) call error( "bad real number" )
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
    write(str,*) a
    str = adjustl(str)
  end function real2str

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

end module mString
