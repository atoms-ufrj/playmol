module mPackmol

use mGlobal
use mStruc

implicit none

contains

  !=================================================================================================

  subroutine Run_Packmol( mol_list, coord_list, seed, tol, nmol, mass, density, aspect )
    type(StrucList), intent(inout) :: mol_list, coord_list
    integer,         intent(in)    :: seed, nmol(:)
    real(rb),        intent(in)    :: tol, mass(:), density, aspect(3)
  end subroutine Run_Packmol

  !=================================================================================================


end module mPackmol
