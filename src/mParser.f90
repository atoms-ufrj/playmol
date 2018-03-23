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
!
!   The present module is a modification of Roland Schmehl's "Fortran 90
!   function parser v1.1" (please see the Copyright statement at the end
!   of this file). The main modification are:
!
!   - Object-oriented version: encapsulation of the data structure in a
!     new class "tParser"
!
!   - Implementation of integer arithmetics, inclusion of new operator %
!     and inclusion of new functions int, nint, ceil, and floor.
!
!   - Change in the precedence order of binary operators, so that two or
!     more operators can have same level of precedence
!
!   - Implementation of boolean arithmetics, including binary relational
!     operators <, <=, >, >=, ==, and !=, binary logical operators & and
!     |, as well as unary logical function not().

module mParser

use mGlobal

implicit none

integer, parameter, private :: is = 1 ! Data type of bytecode

real(rb), parameter, private :: rbtol = 1.0e-8_rb

integer(is), parameter, private :: cImmed   = 1_is,  &
                                   cNeg     = 2_is,  &
                                   cAdd     = 3_is,  &
                                   cSub     = 4_is,  &
                                   cMul     = 5_is,  &
                                   cDiv     = 6_is,  &
                                   cMod     = 7_is,  &
                                   cPow     = 8_is,  &
                                   cLT      = 9_is,  &
                                   cLE      = 10_is, &
                                   cGT      = 11_is, &
                                   cGE      = 12_is, &
                                   cEQ      = 13_is, &
                                   cNE      = 14_is, &
                                   cAnd     = 15_is, &
                                   cOr      = 16_is, &
                                   cNot     = 17_is, &
                                   cAbs     = 18_is, &
                                   cExp     = 19_is, &
                                   cLog10   = 20_is, &
                                   cLog     = 21_is, &
                                   cSqrt    = 22_is, &
                                   cSinh    = 23_is, &
                                   cCosh    = 24_is, &
                                   cTanh    = 25_is, &
                                   cSin     = 26_is, &
                                   cCos     = 27_is, &
                                   cTan     = 28_is, &
                                   cAsin    = 29_is, &
                                   cAcos    = 30_is, &
                                   cAtan    = 31_is, &
                                   cInt     = 32_is, &
                                   cNint    = 33_is, &
                                   cCeil    = 34_is, &
                                   cFloor   = 35_is, &
                                   VarBegin = 36_is

character, parameter, private :: Ops(cAdd:cOr) =  ["+", &
                                                   "-", &
                                                   "*", &
                                                   "/", &
                                                   "%", &
                                                   "^", &
                                                   "<", &
                                             char(243), &
                                                   ">", &
                                             char(242), &
                                             char(240), &
                                             char(216), &
                                                   "&", &
                                                   "|"  ]

integer,   parameter, private :: Prec(cAdd:cOr) = [ 4, & ! +
                                                    4, & ! -
                                                    5, & ! *
                                                    5, & ! /
                                                    5, & ! %
                                                    6, & ! ^
                                                    3, & ! <
                                                    3, & ! <=
                                                    3, & ! >
                                                    3, & ! >=
                                                    3, & ! ==
                                                    3, & ! !=
                                                    2, & ! &
                                                    1  ] ! |

character(5), parameter, private :: Funcs(cNot:cFloor) = ["not  ", &
                                                          "abs  ", &
                                                          "exp  ", &
                                                          "log  ", &
                                                          "ln   ", &
                                                          "sqrt ", &
                                                          "sinh ", &
                                                          "cosh ", &
                                                          "tanh ", &
                                                          "sin  ", &
                                                          "cos  ", &
                                                          "tan  ", &
                                                          "asin ", &
                                                          "acos ", &
                                                          "atan ", &
                                                          "int  ", &
                                                          "nint ", &
                                                          "ceil ", &
                                                          "floor"  ]

type tParser
  private
  logical, public :: Is_Integer
  integer :: NRealVar, NIntVar
  integer :: ByteCodesize
  integer :: Immedsize
  integer :: Stacksize, StackPtr
  integer :: EvalErrType
  integer(is), allocatable :: ByteCode(:)
  real(rb),    allocatable :: Immed(:)
  logical,     allocatable :: IntImmed(:)
  real(rb),    allocatable :: Stack(:)
  logical,     allocatable :: IntStack(:)
  contains
    procedure :: parse => tParser_parse
    procedure :: evaluate => tParser_evaluate
    procedure :: error_msg => tParser_error_msg
    procedure :: Compile => tParser_Compile
    procedure :: AddCompiledByte => tParser_AddCompiledByte
    procedure :: MathItemIndex => tParser_MathItemIndex
    procedure :: CompileSubstr => tParser_CompileSubstr
    procedure :: destroy => tParser_Destroy
end type tParser

! Module subroutines:
private :: tParser_parse,           &
           CheckSyntax,             &
           RemoveSpaces,            &
           Replace,                 &
           tParser_Compile,         &
           tParser_AddCompiledByte, &
           tParser_CompileSubstr,   &
           LowCase

! Module functions:
private :: tParser_evaluate,      &
           tParser_error_msg,     &
           OperatorIndex,         &
           MathFunctionIndex,     &
           VariableIndex,         &
           tParser_MathItemIndex, &
           CompletelyEnclosed,    &
           IsBinaryOp,            &
           realNum,               &
           IsInt

contains

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Parse function string FuncStr and compile it into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_parse( Comp, FuncStr, RealVar, IntVar )
    class(tParser), intent(inout)        :: Comp
    character(*),   intent(in)           :: FuncStr    ! function string
    character(*),   intent(in), optional :: RealVar(:) ! Real variable names
    character(*),   intent(in), optional :: IntVar(:)  ! Integer variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    character(len(FuncStr)) :: Func                    ! function string, local use
    character(sl), allocatable :: Var(:)
    integer :: Nr, Ni, ipos(len_trim(FuncStr))
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    call Comp % destroy
    if (present(RealVar)) then; Nr = size(RealVar); else; Nr = 0; end if
    if (present(IntVar )) then; Ni = size(IntVar);  else; Ni = 0; end if
    allocate( Var(Nr+Ni) )
    if (present(RealVar)) Var(1:Nr) = RealVar
    if (present(IntVar))  Var(Nr+1:Nr+Ni) = IntVar
    Comp%NRealVar = Nr
    Comp%NIntVar  = Ni
    Func = FuncStr                                  ! Local copy of function string
    call Replace( ["<=",">=","==","!="], Ops([cLE,cGE,cEQ,cNE]), Func)
    call RemoveSpaces( Func, ipos )                 ! Condense function string
    call CheckSyntax( Func, FuncStr, Var, ipos )
    call Comp % Compile( Func, Var )                ! Compile into bytecode
  end subroutine tParser_parse

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Evaluate bytecode of function for the values passed in arrays realval and intval
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_evaluate( Comp, realval, intval, ierr ) result ( res )
    class(tParser),     intent(inout) :: Comp
    real(rb), optional, intent(in)    :: realval(:)  ! Real variable values
    integer,  optional, intent(in)    :: intval(:)   ! Integer variable values
    integer,  optional, intent(out)   :: ierr        ! Error code
    real(rb)                          :: res         ! result
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: IP, & ! Instruction pointer
               DP, & ! Data pointer
               SP    ! Stack pointer
    integer :: BC, Nr, Ni, ivar
    character(4) :: C
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (present(realval)) then; Nr = size(realval); else; Nr = 0; end if
    if (present(intval))  then; Ni = size(intval);  else; Ni = 0; end if
    if (Nr /= Comp%NRealVar) then
      write(C,'(I4)') Comp%NRealVar
      call error( "values of", C, "real variables are required" )
    end if
    if (Ni /= Comp%NIntVar) then
      write(C,'(I4)') Comp%NIntVar
      call error( "values of", C, "integer variables are required" )
    end if
    DP = 1
    SP = 0
    Comp%EvalErrType = 0
    do IP = 1, Comp%ByteCodesize

      BC = Comp%ByteCode(IP)
      select case (BC)
        case (cExp:cAtan);  Comp%IntStack(SP) = .false.
        case (cInt:cFloor); Comp%IntStack(SP) = .true.
        case (cLT:cOr);     Comp%IntStack(SP-1) = .true.
      end select
      select case (BC)

        case (cImmed)
          SP = SP + 1
          Comp%Stack(SP) = Comp%Immed(DP)
          Comp%IntStack(SP) = Comp%IntImmed(DP)
          DP = DP + 1

        case (cNeg)
          Comp%Stack(SP) = -Comp%Stack(SP)

        case (cAdd)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) + Comp%Stack(SP)
          SP = SP - 1

        case (cSub)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) - Comp%Stack(SP)
          SP = SP - 1

        case (cMul)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) * Comp%Stack(SP)
          SP = SP - 1

        case (cDiv)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          if (is_zero(SP)) then
            Comp%EvalErrType = 1
          else if (Comp%IntStack(SP-1)) then
            Comp%Stack(SP-1) = nint(Comp%Stack(SP-1))/nint(Comp%Stack(SP))
          else
            Comp%Stack(SP-1) = Comp%Stack(SP-1)/Comp%Stack(SP)
          end if
          SP = SP - 1

        case (cMod)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          if (is_zero(SP)) then
            Comp%EvalErrType = 1
          else if (Comp%IntStack(SP-1)) then
            Comp%Stack(SP-1) = mod(nint(Comp%Stack(SP-1)),nint(Comp%Stack(SP)))
          else
            Comp%Stack(SP-1) = mod(Comp%Stack(SP-1),Comp%Stack(SP))
          end if
          SP = SP - 1

        case (cPow)
          Comp%IntStack(SP-1) = Comp%IntStack(SP-1) .and. Comp%IntStack(SP)
          if (Comp%IntStack(SP-1)) then
            Comp%Stack(SP-1) = nint(Comp%Stack(SP-1))**nint(Comp%Stack(SP))
          else if (Comp%IntStack(SP)) then
            Comp%Stack(SP-1) = Comp%Stack(SP-1)**nint(Comp%Stack(SP))
          else
            Comp%Stack(SP-1) = Comp%Stack(SP-1)**Comp%Stack(SP)
          end if
          SP = SP - 1

        case (cLT)
          Comp%Stack(SP-1) = merge(one,zero,Comp%Stack(SP-1) < Comp%Stack(SP))
          SP = SP - 1

        case (cLE)
          Comp%Stack(SP-1) = merge(one,zero,Comp%Stack(SP-1) <= Comp%Stack(SP))
          SP = SP - 1

        case (cGT)
          Comp%Stack(SP-1) = merge(one,zero,Comp%Stack(SP-1) > Comp%Stack(SP))
          SP = SP - 1

        case (cGE)
          Comp%Stack(SP-1) = merge(one,zero,Comp%Stack(SP-1) >= Comp%Stack(SP))
          SP = SP - 1

        case (cEQ)
          Comp%Stack(SP-1) = merge(one,zero,abs(Comp%Stack(SP-1) - Comp%Stack(SP)) <= rbtol)
          SP = SP - 1

        case (cNE)
          Comp%Stack(SP-1) = merge(one,zero,abs(Comp%Stack(SP-1) - Comp%Stack(SP)) > rbtol)
          SP = SP - 1

        case (cAnd)
          if (Bool(SP-1).and.Bool(SP)) then
            Comp%Stack(SP-1) = Comp%Stack(SP-1)*Comp%Stack(SP)
          else
            Comp%EvalErrType = 5
          end if
          SP = SP - 1

        case (cOr)
          if (Bool(SP-1).and.Bool(SP)) then
            Comp%Stack(SP-1) = min( one, Comp%Stack(SP-1)+Comp%Stack(SP) )
          else
            Comp%EvalErrType = 5
          end if
          SP = SP - 1

        case (cNot)
          if (Bool(SP)) then
            Comp%Stack(SP) = one - Comp%Stack(SP)
          else
            Comp%EvalErrType = 5
          end if

        case (cAbs)
          Comp%Stack(SP) = abs(Comp%Stack(SP))

        case (cExp)
          Comp%Stack(SP) = exp(Comp%Stack(SP))

        case (cLog10)
          if (Comp%Stack(SP) <= zero) then
            Comp%EvalErrType = 3
          else
            Comp%Stack(SP) = log10(Comp%Stack(SP))
          end if

        case (cLog)
          if (Comp%Stack(SP) <= zero) then
            Comp%EvalErrType = 3
          else
            Comp%Stack(SP) = log(Comp%Stack(SP))
          end if

        case (cSqrt)
          if (Comp%Stack(SP) < zero) then
            Comp%EvalErrType = 3
          else
          Comp%Stack(SP) = sqrt(Comp%Stack(SP))
          end if

        case (cSinh)
          Comp%Stack(SP) = sinh(Comp%Stack(SP))

        case (cCosh)
          Comp%Stack(SP) = cosh(Comp%Stack(SP))

        case (cTanh)
          Comp%Stack(SP) = tanh(Comp%Stack(SP))

        case (cSin)
          Comp%Stack(SP) = sin(Comp%Stack(SP))

        case (cCos)
          Comp%Stack(SP) = cos(Comp%Stack(SP))

        case (cTan)
          Comp%Stack(SP) = tan(Comp%Stack(SP))

        case (cAsin)
          if (abs(Comp%Stack(SP)) > one) then
            Comp%EvalErrType = 4
          else
            Comp%Stack(SP) = asin(Comp%Stack(SP))
          end if

        case (cAcos)
          if (abs(Comp%Stack(SP)) > one) then
            Comp%EvalErrType = 4
          else
            Comp%Stack(SP) = acos(Comp%Stack(SP))
          end if

        case (cAtan)
          Comp%Stack(SP) = atan(Comp%Stack(SP))

        case (cInt)
          Comp%Stack(SP) = aint(Comp%Stack(SP))

        case (cNint)
          Comp%Stack(SP) = anint(Comp%Stack(SP))

        case (cCeil)
          Comp%Stack(SP) = ceiling(Comp%Stack(SP))

        case (cFloor)
          Comp%Stack(SP) = floor(Comp%Stack(SP))

        case default
          SP = SP + 1
          ivar = BC - VarBegin + 1
          Comp%IntStack(SP) = ivar > Nr
          if (Comp%IntStack(SP)) then
            Comp%Stack(SP) = intval(ivar-Nr)
          else
            Comp%Stack(SP) = realval(ivar)
          end if
        end select

        if (Comp%EvalErrType > 0) then
          if (present(ierr)) then
            res = zero
            ierr = Comp%EvalErrType
            return
          else
            call error( Comp % error_msg() )
          end if
        end if
    end do
    res = Comp%Stack(1)
    Comp%Is_Integer = Comp%IntStack(1)
    contains
      !- -------- --------- --------- --------- --------- --------- --------- -------
      logical function bool( i )
        integer, intent(in) :: i
        bool = Comp%IntStack(i)
        if (bool) bool = (nint(Comp%Stack(i)) == 1) .or. (nint(Comp%Stack(i)) == 0)
      end function bool
      !- -------- --------- --------- --------- --------- --------- --------- -------
      logical function is_zero( i )
        integer, intent(in) :: i
        if (Comp%IntStack(i)) then
          is_zero = nint(Comp%Stack(i)) == 0
        else
          is_zero = Comp%Stack(i) == zero
        end if
      end function is_zero
      !- -------- --------- --------- --------- --------- --------- --------- -------
  end function tParser_evaluate

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! destroy tParser internal structure
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_destroy( Comp )
    class(tParser), intent(inout) :: Comp
    if (allocated(Comp%ByteCode)) then
      deallocate( Comp%ByteCode, Comp%Immed, Comp%IntImmed, Comp%Stack, Comp%IntStack )
    end if
  end subroutine tParser_destroy

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return error message
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_error_msg( Comp ) result( msg )
    class(tParser), intent(in) :: Comp
    character(sl)              :: msg
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    character(*), parameter :: m(5) = [ "division by zero                ", &
                                        "negative SQRT argument          ", &
                                        "non-positive LOG argument       ", &
                                        "illegal argument of ASIN or ACOS", &
                                        "illegal logical relation        "  ]
    if ((Comp%EvalErrType < 1) .or. (Comp%EvalErrType > size(m))) then
       msg = ""
    else
       msg = m(Comp%EvalErrType)
    end if
  end function tParser_error_msg

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check syntax of function string,  returns 0 if syntax is ok
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine CheckSyntax( Func, FuncStr, Var, ipos )
    character(*), intent(in) :: Func         ! function string without spaces
    character(*), intent(in) :: FuncStr      ! Original function string
    character(*), intent(in) :: Var(:)       ! Array with variable names
    integer,      intent(in) :: ipos(:)
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: n
    character   :: c
    real(rb)    :: r
    logical     :: err
    integer     :: ParCnt, & ! Parenthesis counter
                   j, ib, in, lFunc
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    j = 1
    ParCnt = 0
    lFunc = len_trim(Func)
    step: do
      if (j > lFunc) call ErrMsg(j, FuncStr)
      c = Func(j:j)
      ! Check for valid operand (must appear)
      if ((c == "-") .or. (c == "+")) then           ! Check for leading - or +
        j = j + 1
        if (j > lFunc) call ErrMsg(j, FuncStr, "Missing operand")
        c = Func(j:j)
        if (any(c == Ops)) call ErrMsg(j, FuncStr, "Multiple operators")
      end if
      n = MathFunctionIndex(Func(j:))
      if (n > 0) then                                ! Check for math function
        j = j + len_trim(Funcs(n))
        if (j > lFunc) call ErrMsg(j, FuncStr, "Missing function argument")
        c = Func(j:j)
        if (c /= "(") call ErrMsg(j, FuncStr, "Missing opening parenthesis")
      end if
      if (c == "(") then                             ! Check for opening parenthesis
        ParCnt = ParCnt + 1
        j = j + 1
        cycle step
      end if
      if (scan(c,"0123456789.") > 0) then            ! Check for number
        r = realNum(Func(j:),ib,in,err)
        if (err) call ErrMsg(j, FuncStr, "Invalid number format: "//Func(j+ib-1:j+in-2))
        j = j + in - 1
        if (j > lFunc) exit
        c = Func(j:j)
      else                                           ! Check for variable
        n = VariableIndex (Func(j:),Var,ib,in)
        if (n == 0) call ErrMsg(j, FuncStr, "Invalid element: "//Func(j+ib-1:j+in-2))
        j = j + in - 1
        if (j > lFunc) exit
        c = Func(j:j)
      end if
      do while (c == ")")                            ! Check for closing parenthesis
        ParCnt = ParCnt - 1
        if (ParCnt < 0) call ErrMsg(j, FuncStr, "Mismatched parenthesis")
        if (Func(j-1:j-1) == "(") call ErrMsg(j-1, FuncStr, "Empty parentheses")
        j = j + 1
        if (j > lFunc) exit
        c = Func(j:j)
      end do
      ! Now, we have a legal operand: A legal operator or end of string must follow
      if (j > lFunc) exit
      if (any(c == Ops)) then                        ! Check for multiple operators
        if (j+1 > lFunc) call ErrMsg(j, FuncStr)
        if (any(Func(j+1:j+1) == Ops)) call ErrMsg(j+1, FuncStr, "Multiple operators")
      else                                           ! Check for next operand
        call ErrMsg(j, FuncStr, "Missing operator")
      end if
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      j = j + 1
    end do step
    if (ParCnt > 0) call ErrMsg(j, FuncStr, "Missing )")

    contains
      !- -------- --------- --------- --------- --------- --------- --------- -------
      ! Print error message and terminate program
      subroutine ErrMsg( j, FuncStr, Msg )
        integer,                intent(in) :: j
        character(*),           intent(in) :: FuncStr       ! Original function string
        character(*), optional, intent(in) :: Msg
        call writeln( "Error: wrong syntax" )
        if (present(Msg)) then
           call writeln( ">", Msg )
        end if
        call writeln( ">", FuncStr )
        call writeln( ">", repeat(" ",ipos(j)-1)//"^" )
        stop
      end subroutine ErrMsg

  end subroutine CheckSyntax

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return operator index
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function OperatorIndex( c ) result (n)
    character, intent(in) :: c
    integer(is)           :: n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer (is) :: j
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    do j = cAdd, cOr
      if (c == Ops(j)) then
        n = j
        exit
      end if
    end do
  end function OperatorIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return index of math function beginnig at 1st position of string str
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function MathFunctionIndex( str ) result( n )
    character(*), intent(in) :: str
    integer(is)              :: n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: j
    integer     :: k
    character(len(Funcs)) :: fun
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0_is
    j = cNot - 1_is
    do while (j < VarBegin)                             ! Check all math functions
       j = j + 1_is
       k = min(len_trim(Funcs(j)), len(str))   
       call LowCase (str(1:k), fun)
       if (fun == Funcs(j)) then                        ! Compare lower case letters
          n = j                                         ! Found a matching function
          exit
       end if
    end do
  end function MathFunctionIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return index of variable at begin of string str (returns 0 if no variable found)
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function VariableIndex( str, Var, ibegin, inext ) result (n)
    character(*),      intent(in)  :: str       ! String
    character(*),      intent(in)  :: Var(:)    ! Array with variable names
    integer, optional, intent(out) :: ibegin, & ! Start position of variable name
                                      inext     ! Position of character after name
    integer(is)                    :: n         ! index of variable
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j, ib, in, lstr
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0_is
    lstr = len_trim(str)
    if (lstr > 0) then
       do ib = 1, lstr                          ! Search for first character in str
          if (str(ib:ib) /= " ") exit           ! When lstr>0 at least 1 char in str
       end do                        
       do in = ib, lstr                         ! Search for name terminators
          if (scan(str(in:in),"+-*/^) ") > 0) exit
       end do
       do j = 1, size(Var)
          if (str(ib:in-1) == Var(j)) then                     
             n = int(j,is)                      ! Variable name found
             exit
          end if
       end do
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
  end function VariableIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Remove Spaces from string, remember positions of characters in old string
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine RemoveSpaces( str, ipos )
    character(*), intent(inout) :: str
    integer,      intent(inout) :: ipos(:)
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: k, lstr
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    lstr = len_trim(str)
    ipos = [(k,k=1,lstr)]
    k = 1
    do while (str(k:lstr) /= " ")                             
       if (str(k:k) == " ") then
          str(k:lstr)  = str(k+1:lstr)//" "    ! Move 1 character to left
          ipos(k:lstr) = [ipos(k+1:lstr), 0]   ! Move 1 element to left
          k = k-1
       end if
       k = k+1
    end do
  end subroutine RemoveSpaces

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Replace ALL appearances of character set ca in string str by character set cb
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine Replace( ca, cb, str )
    character(2), intent(in)    :: ca(:)
    character,    intent(in)    :: cb(size(ca))           ! size(cb) must be size(ca)
    character(*), intent(inout) :: str
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer  :: i, j
    do i = 1, size(ca)
      do j = 1, len_trim(str)-1
        if (str(j:j+1) == ca(i)) str(j:j+1) = cb(i)//" "
      end do
    end do
  end subroutine Replace

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Compile i-th function string F into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_Compile( Comp, F, Var )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F            ! function string
    character(*),   intent(in)    :: Var(:)       ! Array with variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: istat
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    call Comp % destroy
    Comp%ByteCodesize = 0
    Comp%Immedsize    = 0
    Comp%Stacksize    = 0
    Comp%StackPtr     = 0
    call Comp % CompileSubstr(F,1,len_trim(F),Var) ! Compile string to determine size
    allocate( Comp%ByteCode(Comp%ByteCodesize), &
              Comp%Immed(Comp%Immedsize),       &
              Comp%IntImmed(Comp%Immedsize),       &
              Comp%Stack(Comp%Stacksize),       &
              Comp%IntStack(Comp%Stacksize),    &
              STAT = istat                      )
    if (istat /= 0) then
      call error( "Memmory allocation for byte code failed" )
    else
      Comp%ByteCodesize = 0
      Comp%Immedsize    = 0
      Comp%Stacksize    = 0
      Comp%StackPtr     = 0
      call Comp % CompileSubstr(F,1,len_trim(F),Var) ! Compile string into bytecode
    end if
  end subroutine tParser_Compile

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Add compiled byte to bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_AddCompiledByte( Comp, b )
    class(tParser), intent(inout) :: Comp
    integer(is),    intent(in)    :: b                   ! Value of byte to be added
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    Comp%ByteCodesize = Comp%ByteCodesize + 1
    if (allocated(Comp%ByteCode)) Comp%ByteCode(Comp%ByteCodesize) = b
  end subroutine tParser_AddCompiledByte

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return math item index, if item is real number, enter it into Comp-structure
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_MathItemIndex( Comp, F, Var ) result( n )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F         ! function substring
    character(*),   intent(in)    :: Var(:)    ! Array with variable names
    integer(is)                   :: n         ! Byte value of math item
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    if (scan(F(1:1),"0123456789.") > 0) then            ! Check for begin of a number
      Comp%Immedsize = Comp%Immedsize + 1
      if (allocated(Comp%Immed)) then
        Comp%Immed(Comp%Immedsize) = RealNum(F)
        Comp%IntImmed(Comp%Immedsize) = IsInt(F)
      end if
      n = cImmed
    else                                                ! Check for a variable
      n = VariableIndex(F, Var)
      if (n > 0) n = VarBegin + n - 1_is
    end if
  end function tParser_MathItemIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function CompletelyEnclosed( F, b, e ) result (res)
    character(*), intent(in) :: F                   ! function substring
    integer,      intent(in) :: b,e                 ! First and last pos. of substring
    logical                  :: res
    !---- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j, k
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    res = .false.
    if (F(b:b) == "(" .and. F(e:e) == ")") then
      k = 0
      do j = b+1, e-1
        if (F(j:j) == "(") then
          k = k+1
        else if (F(j:j) == ")") then
          k = k-1
        end if
        if (k < 0) exit
      end do
      if (k == 0) res = .true.                      ! All opened parenthesis closed
    end if
  end function CompletelyEnclosed

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Compile function string F into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  recursive subroutine tParser_CompileSubstr( Comp, F, b, e, Var )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F         ! function substring
    integer,        intent(in)    :: b,e       ! Begin and end position substring
    character(*),   intent(in)    :: Var(:)    ! Array with variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: n, indx
    integer     :: b2, j, k, iprec
    character(*), parameter :: calpha = "abcdefghijklmnopqrstuvwxyz"// &
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    ! Check for special cases of substring
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (F(b:b) == "+") then                              ! case 1: F(b:e) = "+..."
!     write(*,*)"1. F(b:e) = "+...""
      call Comp % CompileSubstr( F, b+1, e, Var )
      return
    else if (CompletelyEnclosed( F, b, e )) then         ! case 2: F(b:e) = "(...)"
!     write(*,*)"2. F(b:e) = "(...)""
      call Comp % CompileSubstr( F, b+1, e-1, Var )
      return
    else if (scan(F(b:b),calpha) > 0) then        
      n = MathFunctionIndex( F(b:e) )
      if (n > 0) then
        b2 = b + index(F(b:e),"(") - 1
        if (CompletelyEnclosed( F, b2, e )) then         ! case 3: F(b:e) = "fcn(...)"
!         write(*,*)"3. F(b:e) = "fcn(...)""
          call Comp % CompileSubstr( F, b2+1, e-1, Var )
          call Comp % AddCompiledByte( n )
          return
        end if
      end if
    else if (F(b:b) == "-") then
      if (CompletelyEnclosed(F, b+1, e)) then           ! case 4: F(b:e) = "-(...)"
!       write(*,*)"4. F(b:e) = "-(...)""
        call Comp % CompileSubstr( F, b+2, e-1, Var )
        call Comp % AddCompiledByte( cNeg )
        return
      else if (scan(F(b+1:b+1),calpha) > 0) then
        n = MathFunctionIndex(F(b+1:e))
        if (n > 0) then
          b2 = b + index(F(b+1:e),"(")
          if (CompletelyEnclosed(F, b2, e)) then      ! case 5: F(b:e) = "-fcn(...)"
!           write(*,*)"5. F(b:e) = "-fcn(...)""
            call Comp % CompileSubstr( F, b2+1, e-1, Var )
            call Comp % AddCompiledByte( n )
            call Comp % AddCompiledByte( cNeg )
            return
          end if
        end if
      end if
    end if
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    do iprec = 1, maxval(Prec)
      k = 0
      do j = e, b, -1
        if (F(j:j) == ")") then
           k = k + 1
        else if (F(j:j) == "(") then
           k = k - 1
        end if
        if (k == 0) then
          indx = OperatorIndex(F(j:j))
          if (indx /= 0) then
            if ((Prec(indx) == iprec) .and. IsBinaryOp(j, F)) then
              if ((indx >= cMul) .and. (F(b:b) == "-")) then ! case 6: F(b:e) = "-...Op..." with Op > -
!               write(*,*)"6. F(b:e) = "-...Op..." with Op > -"
                call Comp % CompileSubstr( F, b+1, e, Var )
                call Comp % AddCompiledByte( cNeg )
                return
              else                                                        ! case 7: F(b:e) = "...BinOp..."
!               write(*,*)"7. Binary operator",F(j:j)
                call Comp % CompileSubstr( F, b, j-1, Var )
                call Comp % CompileSubstr( F, j+1, e, Var )
                call Comp % AddCompiledByte( indx )
                Comp%StackPtr = Comp%StackPtr - 1
                return
              end if
            end if
          end if
        end if
      end do
    end do
    ! Check for remaining items, i.e. variables or explicit numbers
    b2 = b
    if (F(b:b) == "-") b2 = b2 + 1
    n = Comp % MathItemIndex( F(b2:e), Var )
!   write(*,*)"8. AddCompiledByte ",n
    call Comp % AddCompiledByte( n )
    Comp%StackPtr = Comp%StackPtr + 1
    if (Comp%StackPtr > Comp%Stacksize) Comp%Stacksize = Comp%Stacksize + 1
    if (b2 > b) call Comp % AddCompiledByte( cNeg )
  end subroutine tParser_CompileSubstr

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check if operator F(j:j) in string F is binary operator
  ! Special cases already covered elsewhere:              (that is corrected in v1.1)
  ! - operator character F(j:j) is first character of string (j=1)
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function IsBinaryOp( j, F ) result( res )
    integer,      intent(in) :: j                       ! Position of Operator
    character(*), intent(in) :: F                       ! String
    logical                  :: res                     ! result
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: k
    logical :: Dflag,Pflag
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    res=.true.
    if ((F(j:j) == "+") .or. (F(j:j) == "-")) then       ! Plus or minus sign:
      if (j == 1) then                                   ! - leading unary operator ?
        res = .false.
      else if (scan(F(j-1:j-1),"+-*/%^(") > 0) then      ! - other unary operator ?
        res = .false.
      else if (scan(F(j+1:j+1),"0123456789") > 0 .and. & ! - in exponent of real number ?
               scan(F(j-1:j-1),"eEdD")       > 0) then
        Dflag = .false.
        Pflag = .false.
        k = j - 1
        do while (k > 1)                                 !   step to the left in mantissa 
          k = k - 1
          if (scan(F(k:k),"0123456789") > 0) then
            Dflag=.true.
          else if (F(k:k) == ".") then
            if (Pflag) then
              exit                                       !   * exit: 2nd appearance of "."
            else
              Pflag=.true.                               !   * mark 1st appearance of "."
            end if
          else
            exit                                         !   * all other characters
          end if
        end do
        if (Dflag .and. (k == 1 .or. scan(F(k:k),"+-*/%^(") > 0)) res = .false.
      end if
    end if
  end function IsBinaryOp

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function realNum( str, ibegin, inext, error ) result (res)
    character(*),      intent(in)  :: str           ! String
    integer, optional, intent(out) :: ibegin,     & ! Start position of real number
                                      inext         ! 1st character after real number
    logical, optional, intent(out) :: error         ! Error flag
    real(rb)                       :: res           ! real number
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer  :: ib,in,istat
    logical  :: Bflag,               & ! .T. at begin of number in str
                InMan,               & ! .T. in mantissa of number
                Pflag,               & ! .T. after 1st "." encountered
                Eflag,               & ! .T. at exponent identifier "eEdD"
                InExp,               & ! .T. in exponent of number
                DInMan,              & ! .T. if at least 1 digit in mant.
                DInExp,              & ! .T. if at least 1 digit in exp.
                err                    ! Local error flag
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    Bflag = .true.
    InMan = .false.
    Pflag = .false.
    Eflag = .false.
    InExp = .false.
    DInMan = .false.
    DInExp = .false.
    ib = 1
    in = 1
    do while (in <= len_trim(str))
      select case (str(in:in))
      case (" ")                               ! Only leading blanks permitted
        ib = ib + 1
        if (InMan .or. Eflag .or. InExp) exit
      case ("+","-")                           ! Permitted only
        if (Bflag) then           
          InMan = .true.
          Bflag = .false.                      ! - at beginning of mantissa
        else if (Eflag) then               
          InExp = .true.
          Eflag = .false.                      ! - at beginning of exponent
        else
          exit                                 ! - otherwise stop
        end if
      case ("0":"9")                           ! Mark
        if (Bflag) then           
          InMan = .true.
          Bflag = .false.                      ! - beginning of mantissa
        else if (Eflag) then               
           InExp = .true.
           Eflag = .false.                     ! - beginning of exponent
        end if
        if (InMan) DInMan=.true.               ! Mantissa contains digit
        if (InExp) DInExp=.true.               ! Exponent contains digit
      case (".")
        if (Bflag) then
          Pflag = .true.                       ! - mark 1st appearance of "."
          InMan = .true.
          Bflag = .false.                      !   mark beginning of mantissa
        else if (InMan .and..not.Pflag) then
          Pflag = .true.                       ! - mark 1st appearance of "."
        else
           exit                                ! - otherwise stop
        end if
      case ("e","E","d","D")                   ! Permitted only
        if (InMan) then
           Eflag = .true.
           InMan = .false.                     ! - following mantissa
        else
           exit                                ! - otherwise stop
        end if
      case default
        exit                                   ! stop at all other characters
      end select
      in = in + 1
    end do
    err = (ib > in-1) .or. (.not.DInMan) .or. ((Eflag.or.InExp).and..not.DInExp)
    if (err) then
      res = zero
    else
      read(str(ib:in-1),*,iostat=istat) res
      err = istat /= 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(error))  error  = err
  end function realNum

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check integer number format: [blanks][+|-][nnn][blanks]
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function IsInt( str ) result( res )
    character(*),      intent(in)  :: str           ! String
    logical                        :: res           ! True if str is an integer
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    character(len_trim(str)) :: C
    integer :: k, n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    C = adjustl(str)
    n = len_trim(C)
    if (n == 0) then
      res = .false.
    else if (n == 1) then
      res = scan(C(1:1),"0123456789") /= 0
    else
      res = scan(C(1:1),"+-0123456789") /= 0
      k = 1
      do while (res .and. (k < n))
        k = k + 1
        res = scan(C(k:k),"+-0123456789") /= 0
      end do
    end if
  end function IsInt

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Transform upper case letters in str1 into lower case letters, result is str2
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine LowCase( str1, str2 )
    character(*), intent(in)  :: str1
    character(*), intent(out) :: str2
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j, k
    character(*), parameter :: lc = "abcdefghijklmnopqrstuvwxyz"
    character(*), parameter :: uc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    str2 = str1
    do j = 1, len_trim(str1)
      k = index(uc,str1(j:j))
      if (k > 0) str2(j:j) = lc(k:k)
    end do
  end subroutine LowCase

end module mParser

!------- -------- --------- --------- --------- --------- --------- --------- -------
! Fortran 90 function parser v1.1
!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! * Redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
!        
! * Redistributions of source code must retain the above copyright notice, this list
! of conditions and the following disclaimer.
!        
! * Redistributions in binary form must reproduce the above copyright notice, this
! list of conditions and the following disclaimer in the documentation and/or other
! materials provided with the distribution.
!        
! * Neither the name of the copyright holder nor the names of its contributors may
! be used to endorse or promote products derived from this software without
! specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
!
! The source code is available from http://fparser.sourceforge.net
!
! Roland Schmehl <roland.schmehl at alumni.uni-karlsruhe.de>
!------- -------- --------- --------- --------- --------- --------- --------- -------
