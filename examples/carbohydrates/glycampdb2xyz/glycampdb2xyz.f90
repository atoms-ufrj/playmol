program glycampdb2xyz

use mGlobal
use mString

implicit none

integer :: narg
character(sl) :: arg(11), reaction, cskip

integer :: nitems, residue, skip
logical :: found
type Item
  integer :: index
  character(sl) :: atom, xyz
  type(Item), pointer :: next => null()
end type Item
type(Item), pointer :: first => null(), current, other

! Read items of the pdb file:
call next_command( 5, narg, arg )
nitems = 0
residue = 1
do while (narg > 0)
  if (arg(1) == "HETATM") then
    if (associated(first)) then
      allocate( current % next )
      current => current % next
    else
      allocate( first )
      current => first
    end if
    current % atom = arg(3)
    current % xyz = trim(join(arg(6:8)))
    other => first
    found = .false.
    do while (associated(other % next).and.(.not.found))
      found = (other % atom == current % atom).and.(other % index == residue)
      other => other % next
    end do
    if (found) residue = residue + 1
    current % index = residue
    nitems = nitems + 1
  end if
  call next_command( 5, narg, arg )
end do

if (.not.associated(first)) call error( "the pdb is empty" )

narg = command_argument_count()
if (narg > 2) call error( "wrong number of command-line arguments" )
if (narg > 0) then
  call get_command_argument( 1, reaction )
  if (narg == 2) then
    call get_command_argument( 2, cskip )
    skip = str2int(cskip)
  else
    skip = 0
  end if
  call perform_reaction( first, nitems, trim(reaction), skip )
end if

! Write xyz file:
write(*,'(A,/)') trim(int2str(nitems))
current => first
do while (associated(current))
  call writeln( trim(current % atom)//achar(96 + current % index), current % xyz )
  current => current % next
end do

contains
  !-------------------------------------------------------------------------------------------------
  subroutine perform_reaction( first, nitems, option, skip )
    type(Item), pointer, intent(inout) :: first
    integer,             intent(inout) :: nitems
    character(*),        intent(in)    :: option
    integer,             intent(in)    :: skip

    type(Item), pointer :: current, previous
    integer       :: i, j, residue, ncopies
    real(rb)      :: theta, L, R1(3), R2(3), R3(3), R4(3), D(3), x(3), y(3), z(3)
    character(10) :: coord(3)
    character(3)  :: atom(4)
    character(3), allocatable :: erase(:)
    real(rb),     allocatable :: phi(:)

    select case (option)
      case ("deacetylation")
        allocate( erase(7) )
        erase = [ "H2N", "C2N", "O2N", "CME", "H3M", "H2M", "H1M" ]
        atom = [ "C1", "C2", "N2", "HN" ]
        ncopies = 3
        allocate( phi(ncopies) )
        phi = real([60,180,300],rb)
        theta = 109.2
        L = 1.01
      case ("oxidation")
        allocate( erase(4) )
        erase = [ "H61", "H62", "O6 ", "H6O" ]
        atom = [ "O5", "C5", "C6", "O6" ]
        ncopies = 2
        allocate( phi(ncopies) )
        phi = real([0,180],rb)
        theta = 115.0
        L = 1.25
      case default
        call error( "unknown reation - ", option )
    end select

    do while (associated(first).and.(any(erase == first % atom)))
      current => first
      first => first % next
      deallocate( current )
      nitems = nitems - 1
    end do

    if (associated(first)) then
      previous => first
      current => first % next
      do while (associated(current))
        residue = current % index
        if (mod(residue,skip+1) == 0) then
          if (any(erase == current % atom)) then
            previous % next => current % next
            deallocate( current )
            nitems = nitems - 1
          else
            previous => current
            if (current % atom == atom(1)) then
              R1 = get_vector(current % xyz)
            else if (current % atom == atom(2)) then
              R2 = get_vector(current % xyz)
            else if (current % atom == atom(3)) then
              R3 = get_vector(current % xyz)
              x = (R3 - R2)/norm(R3 - R2)
              D = (R1 - R2) - scalar(R1 - R2,x)*x
              y = D/norm(D)
              z = cross(x,y)
              do i = 1, ncopies
                R4 = R3 + L*(cosine(180-theta)*x + sine(180-theta)*(cosine(phi(i))*y + sine(phi(i))*z))
                current => current % next
                allocate( previous % next )
                previous % next % next => current
                current => previous % next
                current % atom = trim(atom(4))//achar(64+i)
                current % index = residue
                do j = 1, 3
                  write(coord(j),'(F10.4)') R4(j)
                  coord(j) = adjustl(coord(j))
                end do
                current % xyz = join(coord)
                nitems = nitems + 1
                previous => current
              end do 
            end if
          end if
        else
          previous => current
        end if
        current => previous % next
      end do
    end if
    deallocate( erase, phi )
  end subroutine perform_reaction
  !-------------------------------------------------------------------------------------------------
  function get_vector( a ) result( b )
    character(*), intent(in) :: a
    real(rb)                 :: b(3)
    character(sl) :: ca(3)
    integer :: i
    call split( a, i, ca )
    do i = 1, 3
      b(i) = str2real(ca(i))
    end do
  end function get_vector
  !-------------------------------------------------------------------------------------------------
  real(rb) function scalar( u, v )
    real(rb), intent(in) :: u(:), v(:)
    scalar = sum(u*v)
  end function scalar
  !-------------------------------------------------------------------------------------------------
  real(rb) function norm( v )
    real(rb), intent(in) :: v(:)
    norm = sqrt(scalar(v,v))
  end function norm
  !-------------------------------------------------------------------------------------------------
  function cross( a, b ) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross
  !-------------------------------------------------------------------------------------------------
  elemental real(rb) function cosine( theta )
    real(rb), intent(in) :: theta
    cosine = cos(0.01745329251994329577_rb*theta)
  end function cosine
  !-------------------------------------------------------------------------------------------------
  elemental real(rb) function sine( theta )
    real(rb), intent(in) :: theta
    sine = sin(0.01745329251994329577_rb*theta)
  end function sine
  !-------------------------------------------------------------------------------------------------
end program glycampdb2xyz
