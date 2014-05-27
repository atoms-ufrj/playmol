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
    integer :: n, i
    i = 1
    match = .true.
    n = min(len_trim(a),len_trim(b))
    do while (match.and.(i <= n))
      if ((a(i:i) == "*").or.(b(i:i) == "*")) return
      match = (a(i:i) == b(i:i)).or.(a(i:i) == "?").or.(b(i:i) == "?")
      i = i + 1
    end do
    if (match) match = len_trim(a) == len_trim(b)
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

end module mString
