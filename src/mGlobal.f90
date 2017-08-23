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

module mGlobal

implicit none

integer,      parameter :: rb = 8      !< Default number of bytes for real numbers
integer,      parameter :: sl = 256    !< Default character string length
character(3), parameter :: csl = "256" !< String with default character string length

integer :: stdout = 6                  !< Standard output unit
integer :: logunit = 0                 !< Output unit for logging

!> A simple random number generator:
type rng
  integer, private :: jsr
  logical, private :: generate = .true.
  real(8), private :: saved
  contains
    procedure :: init => rng_init
    procedure :: i32 => rng_i32
    procedure :: uni => i32rng_uniform
    procedure :: normal => i32rng_normal
    procedure :: letters => rng_letters
end type rng

type tWarning
  character(sl) :: msg
  type(tWarning), pointer :: next => null()
end type tWarning

type(tWarning), pointer :: first_warning => null(), &
                           last_warning => null()

contains

  !=================================================================================================

  subroutine init_log( file )
    character(*), intent(in) :: file
    open( newunit = logunit, file = file, status = "replace" )
  end subroutine init_log

  !=================================================================================================

  subroutine stop_log
    close( logunit )
    logunit = 0
  end subroutine stop_log

  !=================================================================================================

  subroutine end_line()
    write(stdout,'()')
    if (logunit /= 0) write(logunit,'()')
  end subroutine end_line

  !=================================================================================================

  subroutine writeln( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, advance )
    character(*), intent(in), optional :: msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    logical,      intent(in), optional :: advance
    character(sl) :: emsg
    logical :: wait
    emsg = trim(msg)
    if (present(msg1)) emsg = trim(emsg)//" "//trim(msg1)
    if (present(msg2)) emsg = trim(emsg)//" "//trim(msg2)
    if (present(msg3)) emsg = trim(emsg)//" "//trim(msg3)
    if (present(msg4)) emsg = trim(emsg)//" "//trim(msg4)
    if (present(msg5)) emsg = trim(emsg)//" "//trim(msg5)
    if (present(msg6)) emsg = trim(emsg)//" "//trim(msg6)
    if (present(msg7)) emsg = trim(emsg)//" "//trim(msg7)
    if (present(msg8)) emsg = trim(emsg)//" "//trim(msg8)
    if (present(msg9)) emsg = trim(emsg)//" "//trim(msg9)
    wait = present(advance)
    if (wait) wait = .not.advance
    if (wait) then
      write(stdout,'(A)',advance="no") trim(emsg)
      if (logunit /= 0) write(logunit,'(A)',advance="no") trim(emsg)
    else
      write(stdout,'(A)') trim(emsg)
      if (logunit /= 0) write(logunit,'(A)') trim(emsg)
    end if
  end subroutine writeln

  !=================================================================================================

  subroutine warning( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    character(sl) :: wmsg
    if (.not.associated(first_warning)) then
      allocate( first_warning )
      last_warning => first_warning
    else
      allocate( last_warning % next )
      last_warning => last_warning % next
    end if
    wmsg = trim(msg)
    if (present(msg1)) wmsg = trim(wmsg)//" "//trim(msg1)
    if (present(msg2)) wmsg = trim(wmsg)//" "//trim(msg2)
    if (present(msg3)) wmsg = trim(wmsg)//" "//trim(msg3)
    if (present(msg4)) wmsg = trim(wmsg)//" "//trim(msg4)
    if (present(msg5)) wmsg = trim(wmsg)//" "//trim(msg5)
    if (present(msg6)) wmsg = trim(wmsg)//" "//trim(msg6)
    if (present(msg7)) wmsg = trim(wmsg)//" "//trim(msg7)
    if (present(msg8)) wmsg = trim(wmsg)//" "//trim(msg8)
    if (present(msg9)) wmsg = trim(wmsg)//" "//trim(msg9)
    write(stdout,'(A)') achar(27)//"[1;96mWARNING: "//trim(wmsg)//achar(27)//"[0m"
    if (logunit /= 0) write(logunit,'(A)') "WARNING: "//trim(wmsg)
    last_warning % msg = wmsg
  end subroutine warning

  !=================================================================================================

  subroutine reprint_warnings
    type(tWarning), pointer :: ptr
    ptr => first_warning
    if (associated(ptr)) then
      call writeln()
      write(stdout,'(A)') achar(27)//"[1;96m*** SUMMARY OF WARNINGS ***"//achar(27)//"[0m"
      if (logunit /= 0) write(logunit,'("********** SUMMARY OF WARNINGS **********")')
    end if
    do while (associated(ptr))
      call writeln( ">", ptr % msg )
      ptr => ptr % next
    end do
  end subroutine reprint_warnings

  !=================================================================================================

  subroutine error( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    character(sl) :: emsg
    emsg = trim(msg)
    if (present(msg1)) emsg = trim(emsg)//" "//trim(msg1)
    if (present(msg2)) emsg = trim(emsg)//" "//trim(msg2)
    if (present(msg3)) emsg = trim(emsg)//" "//trim(msg3)
    if (present(msg4)) emsg = trim(emsg)//" "//trim(msg4)
    if (present(msg5)) emsg = trim(emsg)//" "//trim(msg5)
    if (present(msg6)) emsg = trim(emsg)//" "//trim(msg6)
    if (present(msg7)) emsg = trim(emsg)//" "//trim(msg7)
    if (present(msg8)) emsg = trim(emsg)//" "//trim(msg8)
    if (present(msg9)) emsg = trim(emsg)//" "//trim(msg9)
    write(stdout,'(A)') achar(27)//"[1;91mERROR"//achar(27)//"[0m: "//trim(emsg)
    if (logunit /= 0) write(logunit,'(A)') "ERROR: "//trim(emsg)
    stop
  end subroutine error

  !=================================================================================================

  subroutine delete_files( filename )
    character(sl),   intent(in) :: filename(:)
    integer :: ifile, unit, ierr
    do ifile = 1, size(filename)
      open( newunit = unit, file = filename(ifile), status = "old", iostat = ierr )
      if (ierr == 0) close(unit, status = "delete")
    end do
  end subroutine delete_files

  !=================================================================================================

  subroutine rng_init( a, seed )
    class(rng), intent(inout) :: a
    integer,    intent(in)    :: seed
    a%jsr = seed
  end subroutine rng_init

  !=================================================================================================

  function rng_i32( a ) result( i32 )
    class(rng), intent(inout) :: a
    integer                   :: i32
    integer :: jz
    jz = a%jsr
    a%jsr = ieor(a%jsr,ishft(a%jsr, 13))
    a%jsr = ieor(a%jsr,ishft(a%jsr,-17))
    a%jsr = ieor(a%jsr,ishft(a%jsr,  5))
    i32 = jz + a%jsr
  end function rng_i32

  !=================================================================================================

  function i32rng_uniform( a ) result( uni )
    class(rng), intent(inout) :: a
    real(8)                   :: uni
    uni = 0.5_8 + 0.2328306e-9_8 * a%i32()
  end function i32rng_uniform

  !=================================================================================================

  function i32rng_normal( a ) result( normal )
    class(rng), intent(inout) :: a
    real(8)                   :: normal
    integer  :: i1, i2
    real(rb) :: a1, a2
    if (a%generate) then
      i1 = -huge(1)
      do while (i1 == -huge(1))
        i1 = a%i32()
        i2 = a%i32()
      end do
      a1 = sqrt(-2.0_8 * log(0.5_8 + 0.2328306e-9_8*i1))
      a2 = 6.283185307179586477_8*(0.5_8 + 0.2328306e-9_8*i2)
      a%saved = a1*cos(a2)
      normal  = a1*sin(a2)
    else
      normal = a%saved
    end if
    a%generate = .not. a%generate
  end function i32rng_normal

  !=================================================================================================

  function rng_letters( a, n ) result( word )
    class(rng), intent(inout) :: a
    integer,    intent(in)    :: n
    character(sl)             :: word
    integer :: i
    word = ""
    do i = 1, n
      word = trim(word)//achar(97+mod(abs(a % i32()),26))
    end do
  end function rng_letters

  !=================================================================================================

  real(rb) function scalar( u, v )
    real(rb), intent(in) :: u(:), v(:)
    scalar = sum(u*v)
  end function scalar

  !=================================================================================================

  real(rb) function norm( v )
    real(rb), intent(in) :: v(:)
    norm = sqrt(scalar(v,v))
  end function norm

  !=================================================================================================

  function cross( a, b ) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross

  !=================================================================================================

  elemental real(rb) function cosine( theta )
    real(rb), intent(in) :: theta
    cosine = cos(0.01745329251994329577_rb*theta)
  end function cosine

  !=================================================================================================

  elemental real(rb) function sine( theta )
    real(rb), intent(in) :: theta
    sine = sin(0.01745329251994329577_rb*theta)
  end function sine

  !=================================================================================================

  function sorted( x ) result( ind )
    integer, intent(in) :: x(:)
    integer             :: ind(size(x))
    integer :: a(size(x)), N, i, j
    logical :: swapped
    N = size(x)
    a = x
    ind = [(i,i=1,N)]
    do j = N-1, 1, -1
      swapped = .false.
      do i = 1, j
        if (a(i) > a(i+1)) then ! increasing order
!        if (a(i) < a(i+1)) then ! decreasing order
          call swap( a(i), a(i+1) )
          call swap( ind(i), ind(i+1) )
          swapped = .true.
        end if
      end do
      if (.not. swapped) exit
    end do
    contains
      subroutine swap( a, b )
        integer, intent(inout) :: a, b
        integer :: temp
        temp = a
        a = b
        b = temp
      end subroutine swap
  end function sorted

  !=================================================================================================

end module mGlobal
