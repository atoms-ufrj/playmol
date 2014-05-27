program playmol
use mData
use mGlobal
implicit none
type(tData)   :: System
integer       :: i, inp
character(sl) :: infile
if (iargc() == 0) call error( "Usage: playmol <file-1> <file-2> ..." )
call init_log( file = "playmol.log" )
do i = 1, iargc()
  call getarg( i, infile )
  open( newunit = inp, file = infile, status = "old" )
  write(*,'("Reading file ",A,"...")') trim(infile)
  call System % Read( inp )
  close(inp)
end do
call stop_log
end program playmol
