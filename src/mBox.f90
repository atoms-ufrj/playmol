module mBox

use mGlobal

implicit none

type tBox
  logical :: has_volume = .false.
  logical :: has_density = .false.
  real(rb) :: volume
  real(rb) :: density
  real(rb) :: length(3)
  real(rb) :: aspect(3) = 1.0_rb
  contains
    procedure :: exists => tBox_exists
    procedure :: define_density => tBox_define_density
    procedure :: define_volume => tBox_define_volume
    procedure :: compute_lengths => tBox_compute_lengths
end type tBox

contains

  !=================================================================================================

  function tBox_exists( box ) result( exists )
    class(tBox), intent(in) :: box
    logical                 :: exists
    exists = box % has_density .neqv. box % has_volume
  end function tBox_exists

  !=================================================================================================

  subroutine tBox_define_density( box, density )
    class(tBox), intent(inout) :: box
    real(rb),    intent(in)    :: density
    if (density <= 0.0_rb) call error( "invalid density value" )
    box % density = density
    box % has_density = .true.
    box % has_volume = .false.
  end subroutine tBox_define_density

  !=================================================================================================

  subroutine tBox_define_volume( box, volume )
    class(tBox), intent(inout) :: box
    real(rb),    intent(in)    :: volume
    if (volume <= 0.0_rb) call error( "invalid volume value" )
    box % volume = volume
    box % has_volume = .true.
    box % has_density = .false.
  end subroutine tBox_define_volume

  !=================================================================================================

  subroutine tBox_compute_lengths( box, mass )
    class(tBox), intent(inout) :: box
    real(rb),    intent(in)    :: mass
    if (.not.box % exists()) call error( "incorrect box specification" )
    if (box % has_density) then
      box % volume = mass / box % density
    else
      box % density = mass / box % volume
    end if
    box % length = (box % volume / product(box % aspect))**(1.0_rb/3.0_rb) * box % aspect
  end subroutine tBox_compute_lengths

  !=================================================================================================

end module mBox
