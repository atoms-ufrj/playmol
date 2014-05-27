!> A module containing global variables and subroutines
!! @author Charlles R. A. Abreu (abreu at eq.ufrj.br)
!! @date 04/10/2014 11:39:15 AM
module mGlobal

integer,      parameter :: rb = 8      !< Default number of bytes for real numbers
integer,      parameter :: sl = 256    !< Default character string length
character(3), parameter :: csl = "256" !< String with default character string length

integer :: stdout = 6                  !< Standard output unit
integer :: logunit = 0                 !< Output unit for logging

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

  subroutine write_msg( prefix, msg )
    character(*), intent(in) :: prefix, msg
    write(stdout,'(A,A)',advance='no') prefix, msg
    if (logunit /= 0) write(logunit,'(A,A)',advance='no') prefix, msg
  end subroutine write_msg

  !=================================================================================================

  subroutine end_line()
    write(stdout,'()')
    if (logunit /= 0) write(logunit,'()')
  end subroutine end_line

  !=================================================================================================

  subroutine writeln( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, advance )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    logical,      intent(in), optional :: advance
    call write_msg( "", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    if (present(advance)) then
      if (advance) call end_line
    else
      call end_line
    end if
  end subroutine writeln

  !=================================================================================================

  subroutine warning( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    call write_msg( "WARNING: ", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    call end_line
  end subroutine warning

  !=================================================================================================

  subroutine error( msg, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9 )
    character(*), intent(in)           :: msg
    character(*), intent(in), optional :: msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9
    call write_msg( "ERROR: ", trim(msg) )
    if (present(msg1)) call write_msg( " ", trim(msg1) )
    if (present(msg2)) call write_msg( " ", trim(msg2) )
    if (present(msg3)) call write_msg( " ", trim(msg3) )
    if (present(msg4)) call write_msg( " ", trim(msg4) )
    if (present(msg5)) call write_msg( " ", trim(msg5) )
    if (present(msg6)) call write_msg( " ", trim(msg6) )
    if (present(msg7)) call write_msg( " ", trim(msg7) )
    if (present(msg8)) call write_msg( " ", trim(msg8) )
    if (present(msg9)) call write_msg( " ", trim(msg9) )
    call end_line
    stop
  end subroutine error

  !=================================================================================================

end module mGlobal
