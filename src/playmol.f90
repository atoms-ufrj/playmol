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
!            Federal University of Rio de Janeiro

program playmol
use mPlaymol
use mGlobal
implicit none
integer :: i, first, inp
character(sl) :: infile
type(tPlaymol) :: System
call writeln( "Playmol (Version: 28 Oct 2019)" )
if (iargc() == 0) call error( "Usage: playmol [-i] <file-1> <file-2> ..." )
call init_log( file = "playmol.log" )
call getarg( 1, infile )
if (infile == "-i") then
  System % detect_angles_and_dihedrals = .false.
  call getarg( 2, infile )
  first = 2
else
  first = 1
end if
if (infile == "-") then
  call System % Read( 5, trim(infile) )
else
  do i = first, iargc()
    call getarg( i, infile )
    open( newunit = inp, file = infile, status = "old" )
    write(*,'("Reading file ",A,"...")') trim(infile)
    call System % Read( inp, trim(infile) )
    close(inp)
  end do
end if
call reprint_warnings
call stop_log
end program playmol
